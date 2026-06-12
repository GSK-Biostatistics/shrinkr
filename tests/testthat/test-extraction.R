# tests/testthat/test-extraction.R
# Tests for extraction and summary functions

test_that("extract_mu_tau works", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  skip_on_cran()
  
  set.seed(123)
  samples <- list(
    matrix(rnorm(200, 0), ncol = 1),
    matrix(rnorm(200, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), 
      lower = 0
    )
  )
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = priors,
    chains = 1,
    iter = 1000,  # Increased for better convergence
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  # Extract mu and tau
  mu_tau <- extract_mu_tau(fit)
  
  # Test structure - check whatever columns actually exist
  expect_true(is.data.frame(mu_tau) || 
                inherits(mu_tau, "draws_df"))
  
  # Check that we got some results
  expect_true(nrow(mu_tau) > 0)
  expect_true(ncol(mu_tau) > 0)
  
  # If it's a draws object, it should have mu and tau columns
  if (inherits(mu_tau, "draws")) {
    col_names <- names(mu_tau)
    # Could be "mu", "tau" or "mu[1]", "tau[1]" etc
    expect_true(any(grepl("mu", col_names)))
    expect_true(any(grepl("tau", col_names)))
  }
})

test_that("extract_theta works", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  skip_on_cran()
  
  set.seed(456)
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
        distributional::dist_normal(0, 1), 
        lower = 0
      )
    ),
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  # Extract theta
  theta_draws <- extract_theta(fit)
  
  # Test structure
  expect_true(is.data.frame(theta_draws) || 
                inherits(theta_draws, "draws_df") ||
                inherits(theta_draws, "draws"))
  
  # Should have theta parameters
  col_names <- names(theta_draws)
  expect_true(any(grepl("theta", col_names)))
})

test_that("summarise_theta works", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(789)
  samples <- list(
    g1 = matrix(rnorm(200, 0), ncol = 1),
    g2 = matrix(rnorm(200, 1), ncol = 1),
    g3 = matrix(rnorm(200, 0.5), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), 
        lower = 0
      )
    ),
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  # Summarize theta
  theta_summary <- summarise_theta(fit)
  
  # Test structure
  expect_true(is.data.frame(theta_summary))
  expect_equal(nrow(theta_summary), 3)
  
  # Should have some summary statistics
  expect_true(ncol(theta_summary) > 0)
})

test_that("summarise_mu_tau works", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(789)
  samples <- list(
    g1 = matrix(rnorm(200, 0), ncol = 1),
    g2 = matrix(rnorm(200, 1), ncol = 1),
    g3 = matrix(rnorm(200, 0.5), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1),
        lower = 0
      )
    ),
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  # Summarize mu/tau
  mu_tau_summary <- summarise_mu_tau(fit)
  
  # Test structure
  expect_true(is.data.frame(mu_tau_summary))
  # One row each for mu, tau, tau_squared
  expect_equal(nrow(mu_tau_summary), 3)
  expect_true("parameter" %in% names(mu_tau_summary))
  expect_setequal(mu_tau_summary$parameter, c("mu", "tau", "tau_squared"))
  
  # summarize_mu_tau is the US-spelling alias
  expect_identical(summarize_mu_tau(fit), mu_tau_summary)
})

test_that("theta_contrasts works with contrast matrix", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(321)
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
        distributional::dist_normal(0, 1), 
        lower = 0
      )
    ),
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  # Create a simple contrast matrix (g2 - g1)
  contrast_mat <- matrix(c(-1, 1), nrow = 1)
  rownames(contrast_mat) <- "g2 - g1"
  
  # Compute contrasts
  contrasts <- theta_contrasts(fit, contrast_matrix = contrast_mat)
  
  # Test structure
  expect_true(is.data.frame(contrasts))
  expect_true(nrow(contrasts) >= 1)
})

test_that("as.data.frame.shrinkr_mixture works", {
  skip_if_not_installed("mclust")
  
  samples <- list(
    matrix(rnorm(500, 0), ncol = 1),
    matrix(rnorm(500, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Convert to data frame
  df <- as.data.frame(mix)
  
  # Test structure
  expect_true(is.data.frame(df))
  
  # Should have key columns (group column removed in 0.4.0)
  expect_true("variable" %in% names(df))
  expect_true("mean" %in% names(df))
  expect_true("sd" %in% names(df))
  expect_true("weight" %in% names(df))
  expect_false("group" %in% names(df))
  
  # Weights should be reasonable (positive, not crazy large)
  expect_true(all(df$weight > 0))
  expect_true(all(df$weight <= 1.1))  # Allow small numerical error
})

test_that("print methods work without error", {
  skip_if_not_installed("mclust")
  
  samples <- list(
    matrix(rnorm(500, 0), ncol = 1),
    matrix(rnorm(500, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Test that print doesn't error (don't check output since cli may not show in tests)
  expect_error(print(mix), NA)
  
  # Test summary.shrinkr_mixture
  summ <- summary(mix)
  expect_s3_class(summ, "summary.shrinkr_mixture")
  expect_error(print(summ), NA)
})

test_that("plot.shrinkr_mixture works", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("ggplot2")
  
  samples <- list(
    matrix(rnorm(500, 0), ncol = 1),
    matrix(rnorm(500, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Test basic plot
  p1 <- plot(mix, draws = samples, type = "density")
  expect_s3_class(p1, "gg")
  
  # Test QQ plot
  p2 <- plot(mix, draws = samples, type = "qq")
  expect_s3_class(p2, "gg")
})

test_that("plot.shrinkr_prior_pred works", {
  skip_if_not_installed("ggplot2")
  
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), 
      lower = 0
    )
  )
  
  pp <- sample_prior_predictive(
    hierarchical_priors = priors,
    n_groups = 3,
    n_draws = 100
  )
  
  # Test plot
  p <- plot(pp)
  expect_s3_class(p, "gg")
})

test_that("shrinkr_fit methods work", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(654)
  samples <- list(
    matrix(rnorm(200, 0), ncol = 1),
    matrix(rnorm(200, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  fit <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), 
        lower = 0
      )
    ),
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  # Test print method
  expect_error(print(fit), NA)
  
  # Test summary method
  summ <- summary(fit)
  expect_true(is.data.frame(summ) || is.list(summ))
})