# tests/testthat/test-shrinkr.R
# Basic workflow and core functionality tests

test_that("fit_mixture works with list of matrices", {
  skip_if_not_installed("mclust")
  
  set.seed(123)
  samples <- list(
    matrix(rnorm(500, 0, 0.5), ncol = 1),
    matrix(rnorm(500, 1, 0.5), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
  
  # Test structure
  expect_s3_class(mix, "shrinkr_mixture")
  expect_true(is.list(mix))
  expect_true("K" %in% names(mix))
  expect_true("vars" %in% names(mix))
  expect_true("weights" %in% names(mix))
  
  # Test values
  expect_equal(length(mix$vars), 2)
  expect_true(mix$K >= 1)
  expect_true(all(mix$weights > 0))
  expect_equal(sum(mix$weights), 1, tolerance = 1e-6)
})

test_that("fit_mixture works with data frames", {
  skip_if_not_installed("mclust")
  
  set.seed(456)
  samples_df <- data.frame(
    group1 = rnorm(500, 0),
    group2 = rnorm(500, 1)
  )
  
  mix <- fit_mixture(samples_df, K_max = 2, verbose = FALSE)
  
  expect_s3_class(mix, "shrinkr_mixture")
  expect_equal(length(mix$vars), 2)
})

test_that("fit_mixture validates inputs", {
  skip_if_not_installed("mclust")
  
  # Too few variables
  expect_error(
    fit_mixture(list(matrix(rnorm(100), ncol = 1)), K_max = 1),
    "at least two variables"
  )
  
  # Invalid K_max
  expect_error(
    fit_mixture(list(rnorm(100), rnorm(100)), K_max = 0),
    "must be a single integer"
  )
})

test_that("shrink works with mixture input", {
  skip_if_not_installed("rstan")
  skip_on_cran()  # Slow Stan fitting
  
  set.seed(789)
  samples <- list(
    matrix(rnorm(200, 0), ncol = 1),
    matrix(rnorm(200, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Use distributional package for priors
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
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  # Test structure
  expect_s3_class(fit, "shrinkr_fit")
  expect_true("fit" %in% names(fit))
  expect_true("data" %in% names(fit))
  expect_true("summary" %in% names(fit))
  
  # Test that we got results
  expect_true(is.data.frame(fit$summary))
  expect_true(nrow(fit$summary) > 0)
})

test_that("shrink works with MLE input", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), 
      lower = 0
    )
  )
  
  fit <- suppressWarnings(shrink(
    mle = c(0.5, 1.0, 0.3),
    var_matrix = c(0.1, 0.15, 0.12),
    hierarchical_priors = priors,
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  expect_s3_class(fit, "shrinkr_fit")
  expect_equal(fit$data$G, 3)
})

test_that("shrink validates inputs", {
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), 
      lower = 0
    )
  )
  
  # Missing both mixture and mle - match pattern not exact text (cli formatting)
  expect_error(
    shrink(hierarchical_priors = priors),
    "provide either.*mixture.*mle"
  )
  
  # Missing priors
  samples <- list(matrix(rnorm(100), ncol = 1), matrix(rnorm(100), ncol = 1))
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  expect_error(
    shrink(mixture = mix, hierarchical_priors = list()),
    "must be a list with elements"
  )
  
  # Non-distributional priors - match pattern
  expect_error(
    shrink(mixture = mix, hierarchical_priors = list(mu = 0, tau = 1)),
    "distributional objects"
  )
})

test_that("mixture priors work", {
  # Test spike and slab (now returns distributional::dist_mixture)
  ss <- prior_spike_slab(spike_prob = 0.5)
  expect_s3_class(ss, "distribution")
  # Should be sampleable
  draws <- distributional::generate(ss, 50)
  expect_length(draws[[1]], 50)
  
  # Test custom mixture
  mix <- prior_mixture(
    distributional::dist_normal(0, 0.5),
    distributional::dist_normal(0, 1),
    weights = c(0.3, 0.7)
  )
  expect_s3_class(mix, "distribution")
  draws2 <- distributional::generate(mix, 50)
  expect_length(draws2[[1]], 50)
})

test_that("sample_prior_predictive works", {
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_truncated(
      distributional::dist_normal(0, 1), 
      lower = 0
    )
  )
  
  pp <- sample_prior_predictive(
    hierarchical_priors = priors,
    n_groups = 4,
    n_draws = 100
  )
  
  # Test structure
  expect_s3_class(pp, "shrinkr_prior_pred")
  expect_true("mu" %in% names(pp))
  expect_true("tau" %in% names(pp))
  expect_true("theta" %in% names(pp))
  
  # Test dimensions
  expect_length(pp$mu, 100)
  expect_length(pp$tau, 100)
  expect_equal(dim(pp$theta), c(100, 4))
  
  # Test values
  expect_true(all(is.finite(pp$mu)))
  expect_true(all(pp$tau > 0))
  expect_true(all(is.finite(pp$theta)))
})

test_that("shrink with centered parameterization", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(135)
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
  
  # Centered parameterization
  fit_c <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = priors,
    centered = TRUE,
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  
  expect_s3_class(fit_c, "shrinkr_fit")
  expect_equal(fit_c$data$centered, 1)
})