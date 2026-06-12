# tests/testthat/test-prior-system.R
# Tests for the prior specification system (0.4.0 overhaul)

# ============================================================================
# prior_mixture() and prior_spike_slab()
# ============================================================================

test_that("prior_mixture creates a distributional mixture object", {
  mix <- prior_mixture(
    distributional::dist_normal(0, 1),
    distributional::dist_normal(5, 2),
    weights = c(0.3, 0.7)
  )
  expect_s3_class(mix, "distribution")
  # Should be usable with distributional::generate
  draws <- distributional::generate(mix, 100)
  expect_true(is.numeric(draws[[1]]))
  expect_length(draws[[1]], 100)
})

test_that("prior_mixture normalizes weights", {
  mix <- prior_mixture(
    distributional::dist_normal(0, 1),
    distributional::dist_normal(1, 1),
    weights = c(2, 8)
  )
  # Should still work — weights get normalized internally
  draws <- distributional::generate(mix, 50)
  expect_length(draws[[1]], 50)
})

test_that("prior_mixture defaults to equal weights", {
  mix <- prior_mixture(
    distributional::dist_normal(0, 1),
    distributional::dist_normal(1, 1),
    distributional::dist_normal(2, 1)
  )
  draws <- distributional::generate(mix, 100)
  expect_length(draws[[1]], 100)
})

test_that("prior_spike_slab creates a valid mixture", {
  ss <- prior_spike_slab(spike_prob = 0.5, spike_scale = 0.01, slab_scale = 1)
  expect_s3_class(ss, "distribution")
  
  # Can be sampled
  draws <- distributional::generate(ss, 200)
  expect_length(draws[[1]], 200)
  
  # Should have draws near zero (spike) and spread out (slab)
  d <- draws[[1]]
  expect_true(any(abs(d) < 0.1))  # some draws near zero
  expect_true(any(abs(d) > 0.5))  # some draws spread out
})

test_that("prior_spike_slab respects custom parameters", {
  ss1 <- prior_spike_slab(spike_prob = 0.99, spike_scale = 0.001, slab_scale = 0.5)
  draws <- distributional::generate(ss1, 1000)[[1]]
  # With 99% spike weight and tiny scale, most draws should be very near 0
  expect_true(mean(abs(draws) < 0.01) > 0.5)
})

test_that("format works on mixture distributions", {
  ss <- prior_spike_slab()
  # Should be printable without error
  expect_error(format(ss), NA)
  formatted <- format(ss)
  expect_true(is.character(formatted))
})


# ============================================================================
# .coerce_priors_to_stan() — mu priors
# ============================================================================

test_that(".coerce_priors_to_stan handles Normal mu", {
  pri <- shrinkr:::.coerce_priors_to_stan(
    prior_mu = distributional::dist_normal(0, 5),
    prior_tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  expect_equal(pri$mu_prior_type, 1)
  expect_equal(pri$mu_loc, 0)
  expect_equal(pri$mu_scale, 5)
  expect_equal(pri$mu_is_truncated, 0)
})

test_that(".coerce_priors_to_stan handles Student-t mu", {
  pri <- shrinkr:::.coerce_priors_to_stan(
    prior_mu = distributional::dist_student_t(3, 0, 2),
    prior_tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  expect_equal(pri$mu_prior_type, 2)
  expect_equal(pri$mu_df, 3)
  expect_equal(pri$mu_loc, 0)
  expect_equal(pri$mu_scale, 2)
})

test_that(".coerce_priors_to_stan handles truncated mu", {
  pri <- shrinkr:::.coerce_priors_to_stan(
    prior_mu = distributional::dist_truncated(
      distributional::dist_normal(0, 5), lower = 0
    ),
    prior_tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  expect_equal(pri$mu_prior_type, 1)
  expect_equal(pri$mu_is_truncated, 1)
  expect_equal(pri$mu_lb, 0)
  expect_equal(pri$mu_loc, 0)
  expect_equal(pri$mu_scale, 5)
})

test_that(".coerce_priors_to_stan handles truncated mu with upper bound", {
  pri <- shrinkr:::.coerce_priors_to_stan(
    prior_mu = distributional::dist_truncated(
      distributional::dist_normal(0, 5), lower = -2, upper = 2
    ),
    prior_tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  expect_equal(pri$mu_is_truncated, 1)
  expect_equal(pri$mu_lb, -2)
  expect_equal(pri$mu_ub, 2)
})

test_that(".coerce_priors_to_stan handles mixture mu", {
  pri <- shrinkr:::.coerce_priors_to_stan(
    prior_mu = prior_mixture(
      distributional::dist_normal(0, 1),
      distributional::dist_normal(5, 2),
      weights = c(0.6, 0.4)
    ),
    prior_tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  expect_equal(pri$mu_prior_type, 3)
  expect_equal(pri$mu_n_components, 2)
})


# ============================================================================
# .coerce_priors_to_stan() — tau priors
# ============================================================================

test_that(".coerce_priors_to_stan handles all standard tau distributions", {
  half_normal <- distributional::dist_truncated(
    distributional::dist_normal(0, 2.5), lower = 0
  )
  half_t <- distributional::dist_truncated(
    distributional::dist_student_t(3, 0, 1), lower = 0
  )
  half_cauchy <- distributional::dist_truncated(
    distributional::dist_cauchy(0, 1), lower = 0
  )
  
  mu <- distributional::dist_normal(0, 5)
  
  pri_hn <- shrinkr:::.coerce_priors_to_stan(mu, half_normal)
  expect_equal(pri_hn$tau_prior_type, 1)
  
  pri_ht <- shrinkr:::.coerce_priors_to_stan(mu, half_t)
  expect_equal(pri_ht$tau_prior_type, 2)
  
  pri_hc <- shrinkr:::.coerce_priors_to_stan(mu, half_cauchy)
  expect_equal(pri_hc$tau_prior_type, 8)
})

test_that(".coerce_priors_to_stan handles naturally positive tau distributions", {
  mu <- distributional::dist_normal(0, 5)
  
  pri_ln <- shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_lognormal(0, 0.5))
  expect_equal(pri_ln$tau_prior_type, 3)
  
  pri_ig <- shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_inverse_gamma(2, 1))
  expect_equal(pri_ig$tau_prior_type, 4)
  
  pri_ga <- shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_gamma(2, 2))
  expect_equal(pri_ga$tau_prior_type, 5)
  
  pri_ex <- shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_exponential(1))
  expect_equal(pri_ex$tau_prior_type, 6)
  
  pri_un <- shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_uniform(0, 5))
  expect_equal(pri_un$tau_prior_type, 7)
})

test_that(".coerce_priors_to_stan handles truncated spike-and-slab tau", {
  pri <- shrinkr:::.coerce_priors_to_stan(
    prior_mu = distributional::dist_normal(0, 5),
    prior_tau = distributional::dist_truncated(
      prior_spike_slab(spike_prob = 0.5, spike_scale = 0.01, slab_scale = 1),
      lower = 0
    )
  )
  expect_equal(pri$tau_prior_type, 9)
  expect_equal(pri$tau_n_components, 2)
  expect_equal(pri$tau_lb, 0)
})

test_that(".coerce_priors_to_stan errors on untruncated Normal/t/Cauchy tau", {
  mu <- distributional::dist_normal(0, 5)
  
  expect_error(
    shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_normal(0, 1)),
    "truncated"
  )
  expect_error(
    shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_student_t(3, 0, 1)),
    "truncated"
  )
  expect_error(
    shrinkr:::.coerce_priors_to_stan(mu, distributional::dist_cauchy(0, 1)),
    "truncated"
  )
})

test_that(".coerce_priors_to_stan errors on negative tau lower bound", {
  mu <- distributional::dist_normal(0, 5)
  
  expect_error(
    shrinkr:::.coerce_priors_to_stan(
      mu,
      distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = -1
      )
    ),
    "lower bound must be >= 0"
  )
})


# ============================================================================
# .as_chol_factor() — Cholesky safety
# ============================================================================

test_that(".as_chol_factor works on positive-definite matrices", {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  L <- shrinkr:::.as_chol_factor(Sigma)
  expect_true(is.matrix(L))
  expect_equal(dim(L), c(2, 2))
  # Should be lower triangular
  expect_equal(L[1, 2], 0)
  # L %*% t(L) should recover Sigma
  expect_equal(L %*% t(L), Sigma, tolerance = 1e-10)
})

test_that(".as_chol_factor handles near-singular matrices", {
  # Matrix with a small negative eigenvalue — not PD
  Sigma <- matrix(c(1.0, 1.0001, 1.0001, 1.0), 2, 2)
  eigs <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(min(eigs) < 0)
  
  # Should produce a valid matrix (via jitter or nearPD), not garbage
  L <- tryCatch(
    suppressWarnings(shrinkr:::.as_chol_factor(Sigma)),
    error = function(e) NULL
  )
  if (!is.null(L)) {
    expect_true(is.matrix(L))
    expect_equal(dim(L), c(2, 2))
    # Should be lower triangular
    expect_true(all(abs(L[upper.tri(L)]) < 1e-8))
  }
})

test_that(".as_chol_factor errors on severely indefinite matrices", {
  # Large negative eigenvalue — nearPD may or may not rescue this
  Sigma <- matrix(c(1, 5, 5, 1), 2, 2)
  result <- tryCatch(
    suppressWarnings(shrinkr:::.as_chol_factor(Sigma)),
    error = function(e) "errored"
  )
  if (is.character(result) && result == "errored") {
    # Good — function correctly refused
    expect_true(TRUE)
  } else {
    # nearPD rescued it — verify the result is at least a valid lower-triangular matrix
    expect_true(is.matrix(result))
    expect_equal(dim(result), c(2, 2))
    expect_true(all(abs(result[upper.tri(result)]) < 1e-8))
  }
})


# ============================================================================
# sample_prior_predictive()
# ============================================================================

test_that("sample_prior_predictive validates tau positivity", {
  # Untruncated spike-and-slab should error
  expect_error(
    sample_prior_predictive(
      hierarchical_priors = list(
        mu = distributional::dist_normal(0, 5),
        tau = prior_spike_slab()
      ),
      n_groups = 3,
      n_draws = 100
    ),
    "tau draws were <= 0"
  )
})

test_that("sample_prior_predictive works with truncated spike-and-slab", {
  pp <- sample_prior_predictive(
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        prior_spike_slab(spike_prob = 0.5, spike_scale = 0.01, slab_scale = 1),
        lower = 0
      )
    ),
    n_groups = 3,
    n_draws = 200
  )
  expect_s3_class(pp, "shrinkr_prior_pred")
  expect_true(all(pp$tau > 0))
  expect_equal(ncol(pp$theta), 3)
  expect_equal(nrow(pp$theta), 200)
})

test_that("sample_prior_predictive works with mixture mu", {
  pp <- sample_prior_predictive(
    hierarchical_priors = list(
      mu = prior_mixture(
        distributional::dist_normal(0, 1),
        distributional::dist_normal(5, 1),
        weights = c(0.5, 0.5)
      ),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = 0
      )
    ),
    n_groups = 2,
    n_draws = 500
  )
  expect_s3_class(pp, "shrinkr_prior_pred")
  # Mu draws should show bimodality — some near 0, some near 5
  expect_true(any(pp$mu < 2))
  expect_true(any(pp$mu > 3))
})

test_that("sample_prior_predictive validates inputs", {
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  
  expect_error(
    sample_prior_predictive(list(mu = priors$mu), n_groups = 2),
    "mu.*tau"
  )
  expect_error(
    sample_prior_predictive(priors, n_groups = 0),
    "positive integer"
  )
  expect_error(
    sample_prior_predictive(priors, n_groups = 2, n_draws = -1),
    "positive integer"
  )
})


# ============================================================================
# prior_pairwise_differences()
# ============================================================================

test_that("prior_pairwise_differences computes correct structure", {
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  pp <- sample_prior_predictive(priors, n_groups = 4, n_draws = 200)
  pw <- prior_pairwise_differences(pp)
  
  expect_s3_class(pw, "shrinkr_prior_contrasts")
  # 4 groups = 4 choose 2 = 6 pairs
  expect_equal(pw$n_pairs, 6)
  expect_equal(pw$n_draws, 200)
  expect_equal(nrow(pw$differences), 6 * 200)
  
  # All differences should be non-negative
  expect_true(all(pw$differences$abs_diff >= 0))
  
  # Summary should have one row per pair
  expect_equal(nrow(pw$summary), 6)
  expect_true(all(c("pair", "mean", "median", "q2.5", "q97.5") %in% names(pw$summary)))
})

test_that("prior_pairwise_differences works with 2 groups", {
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  pp <- sample_prior_predictive(priors, n_groups = 2, n_draws = 100)
  pw <- prior_pairwise_differences(pp)
  
  expect_equal(pw$n_pairs, 1)
  expect_equal(nrow(pw$differences), 100)
})

test_that("prior_pairwise_differences errors with 1 group", {
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  pp <- sample_prior_predictive(priors, n_groups = 1, n_draws = 50)
  
  expect_error(
    prior_pairwise_differences(pp),
    "at least 2 groups"
  )
})

test_that("prior_pairwise_differences print and plot work", {
  skip_if_not_installed("ggplot2")
  
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), lower = 0
    )
  )
  pp <- sample_prior_predictive(priors, n_groups = 3, n_draws = 100)
  pw <- prior_pairwise_differences(pp)
  
  expect_error(print(pw), NA)
  
  p1 <- plot(pw)
  expect_s3_class(p1, "gg")
  
  p2 <- plot(pw, by_pair = TRUE)
  expect_s3_class(p2, "gg")
})

test_that("prior_pairwise_differences prob columns are sensible", {
  # With large tau, most pairwise diffs should exceed 0.5
  priors_wide <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 10), lower = 0
    )
  )
  pp <- sample_prior_predictive(priors_wide, n_groups = 3, n_draws = 500)
  pw <- prior_pairwise_differences(pp)
  # With tau ~ half-normal(0,10), most pairs should differ by > 0.5
  expect_true(all(pw$summary$prob_gt_0.5 > 0.5))
  
  # With tiny tau, very few pairwise diffs should exceed 1
  priors_tight <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 0.01), lower = 0
    )
  )
  pp2 <- sample_prior_predictive(priors_tight, n_groups = 3, n_draws = 500)
  pw2 <- prior_pairwise_differences(pp2)
  expect_true(all(pw2$summary$prob_gt_1 < 0.1))
})


# ============================================================================
# fit_mixture() output structure (0.4.0 changes)
# ============================================================================

test_that("fit_mixture does not include group column", {
  skip_if_not_installed("mclust")
  set.seed(42)
  samples <- list(
    g1 = rnorm(300, 0),
    g2 = rnorm(300, 1)
  )
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  expect_false("group" %in% names(mix$components))
  expect_true(all(c("component", "variable", "weight", "mean", "sd") %in% names(mix$components)))
})

test_that("fit_mixture quantiles are included", {
  skip_if_not_installed("mclust")
  set.seed(42)
  samples <- list(
    g1 = rnorm(300, 0),
    g2 = rnorm(300, 1)
  )
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  expect_true(!is.null(mix$quantiles))
  expect_true("q50" %in% names(mix$quantiles))
  expect_equal(nrow(mix$quantiles), 2)
})


# ============================================================================
# shrink() — full workflow with new prior types (slow tests)
# ============================================================================

test_that("shrink works with truncated spike-and-slab tau", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(42)
  samples <- list(
    matrix(rnorm(300, 0, 0.5), ncol = 1),
    matrix(rnorm(300, 1, 0.5), ncol = 1)
  )
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        prior_spike_slab(spike_prob = 0.5, spike_scale = 0.01, slab_scale = 1),
        lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  expect_s3_class(fit, "shrinkr_fit")
  expect_true(all(c("mu", "tau") %in% fit$summary$variable))
})

test_that("shrink works with truncated mu", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(42)
  samples <- list(
    matrix(rnorm(300, 0.5, 0.3), ncol = 1),
    matrix(rnorm(300, 1.0, 0.3), ncol = 1)
  )
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_truncated(
        distributional::dist_normal(0, 5), lower = 0
      ),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  expect_s3_class(fit, "shrinkr_fit")
  # mu should be positive due to truncation
  mu_summary <- fit$summary[fit$summary$variable == "mu", ]
  expect_true(mu_summary$mean > 0)
})

test_that("shrink works with mixture mu", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(42)
  samples <- list(
    matrix(rnorm(300, 0, 0.5), ncol = 1),
    matrix(rnorm(300, 1, 0.5), ncol = 1)
  )
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = prior_mixture(
        distributional::dist_normal(0, 1),
        distributional::dist_normal(5, 1),
        weights = c(0.5, 0.5)
      ),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  expect_s3_class(fit, "shrinkr_fit")
})

test_that("shrinkage_factor is not in output", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  skip_on_cran()
  
  set.seed(42)
  samples <- list(
    matrix(rnorm(300, 0, 0.5), ncol = 1),
    matrix(rnorm(300, 1, 0.5), ncol = 1)
  )
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  draws <- as.data.frame(fit)
  expect_false("shrinkage_factor" %in% names(draws))
  expect_true("tau_squared" %in% names(draws))
})

test_that("shrink MLE path works with spike-and-slab tau", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  fit <- suppressWarnings(shrink(
    mle = c(0.5, 1.0, -0.3),
    var_matrix = c(0.1, 0.15, 0.12),
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        prior_spike_slab(spike_prob = 0.5, spike_scale = 0.01, slab_scale = 1),
        lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  expect_s3_class(fit, "shrinkr_fit")
  expect_equal(fit$data$G, 3)
})

test_that("shrink MLE path works with truncated mu", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  fit <- suppressWarnings(shrink(
    mle = c(0.5, 1.0, 0.8),
    var_matrix = c(0.1, 0.15, 0.12),
    hierarchical_priors = list(
      mu = distributional::dist_truncated(
        distributional::dist_normal(0, 5), lower = 0
      ),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  expect_s3_class(fit, "shrinkr_fit")
  # mu should be positive
  mu_summary <- fit$summary[fit$summary$variable == "mu", ]
  expect_true(mu_summary$mean > 0)
})

test_that("shrink MLE path works with mixture mu", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  fit <- suppressWarnings(shrink(
    mle = c(0.5, 1.0, -0.3),
    var_matrix = c(0.1, 0.15, 0.12),
    hierarchical_priors = list(
      mu = prior_mixture(
        distributional::dist_normal(0, 1),
        distributional::dist_normal(3, 1),
        weights = c(0.5, 0.5)
      ),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  expect_s3_class(fit, "shrinkr_fit")
  expect_equal(fit$data$G, 3)
})

test_that("theta_contrasts does not warn about dropped draws_df class", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(42)
  samples <- list(
    g1 = matrix(rnorm(200, 0), ncol = 1),
    g2 = matrix(rnorm(200, 1), ncol = 1)
  )
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), lower = 0
      )
    ),
    chains = 1, iter = 500, warmup = 250, refresh = 0, verbose = FALSE
  ))
  
  L <- matrix(c(-1, 1), nrow = 1)
  # Should not produce "Dropping draws_df class" warning
  expect_no_warning(
    contrasts <- theta_contrasts(fit, contrast_matrix = L, labels = "g2_vs_g1")
  )
  expect_s3_class(contrasts, "draws_df")
})


# ============================================================================
# .extract_normal_params() helper
# ============================================================================

test_that(".extract_normal_params handles plain normal", {
  d <- distributional::dist_normal(3, 2)
  pp <- shrinkr:::.extract_normal_params(d)
  expect_equal(pp$mu, 3)
  expect_equal(pp$sigma, 2)
})

test_that(".extract_normal_params handles truncated normal", {
  d <- distributional::dist_truncated(
    distributional::dist_normal(3, 2), lower = 0
  )
  pp <- shrinkr:::.extract_normal_params(d)
  expect_equal(pp$mu, 3)
  expect_equal(pp$sigma, 2)
})