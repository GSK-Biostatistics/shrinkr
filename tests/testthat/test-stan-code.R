# tests/testthat/test-stan-code.R
# Tests for Stan code generation (shrinkr's unique contribution)

# NOTE: stan_code() is an internal generic used by shrink() during prior
# translation; it is not exported. We don't test it directly — instead we
# verify that the priors it supports (normal, student-t, lognormal,
# inverse-gamma, uniform, gamma, exponential) all work end-to-end through
# shrink() on the full pipeline.

test_that("shrink accepts various prior distributions", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(123)
  samples <- list(
    matrix(rnorm(200, 0), ncol = 1),
    matrix(rnorm(200, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Test that various priors work through the full workflow
  prior_list <- list(
    normal = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        distributional::dist_normal(0, 1), 
        lower = 0
      )
    ),
    student_t = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_truncated(
        distributional::dist_student_t(3, 0, 1), 
        lower = 0
      )
    ),
    lognormal = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_lognormal(0, 0.5)
    )
  )
  
  # Test that each prior type works (suppress MCMC warnings since we're testing functionality not convergence)
  for (prior_name in names(prior_list)) {
    expect_error({
      suppressWarnings({
        fit <- shrink(
          mixture = mix,
          hierarchical_priors = prior_list[[prior_name]],
          chains = 1,
          iter = 500,
          warmup = 250,
          refresh = 0,
          verbose = FALSE
        )
      })
    }, NA)  # NA means no error expected
  }
})

test_that("shrink handles inverse-gamma prior", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(456)
  samples <- list(
    matrix(rnorm(200, 0), ncol = 1),
    matrix(rnorm(200, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Inverse-gamma prior
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_inverse_gamma(2, 1)
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
  
  expect_s3_class(fit, "shrinkr_fit")
  expect_true(is.list(fit$data))
})

test_that("shrink handles uniform prior", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(789)
  samples <- list(
    matrix(rnorm(200, 0), ncol = 1),
    matrix(rnorm(200, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Uniform prior
  priors <- list(
    mu = distributional::dist_normal(0, 5),
    tau = distributional::dist_uniform(0, 5)
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
  
  expect_s3_class(fit, "shrinkr_fit")
})

test_that("shrink handles gamma and exponential priors", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  
  set.seed(321)
  samples <- list(
    matrix(rnorm(200, 0), ncol = 1),
    matrix(rnorm(200, 1), ncol = 1)
  )
  
  mix <- fit_mixture(samples, K_max = 1, verbose = FALSE)
  
  # Gamma prior
  fit_gamma <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_gamma(2, 2)
    ),
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  expect_s3_class(fit_gamma, "shrinkr_fit")
  
  # Exponential prior
  fit_exp <- suppressWarnings(shrink(
    mixture = mix,
    hierarchical_priors = list(
      mu = distributional::dist_normal(0, 5),
      tau = distributional::dist_exponential(1)
    ),
    chains = 1,
    iter = 1000,
    warmup = 500,
    refresh = 0,
    verbose = FALSE
  ))
  expect_s3_class(fit_exp, "shrinkr_fit")
})