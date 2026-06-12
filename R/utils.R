#' Internal null coalescing operator
#'
#' Returns `y` if `x` is `NULL`, otherwise returns `x`.
#' Used internally to simplify default argument handling.
#'
#' @name grapes-or-or-grapes
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Coerce posterior draws to a data frame with stable names
#'
#' @description
#' Accepts the same input shapes as [fit_mixture()] and returns a
#' `data.frame` with column names that match the names created during
#' fitting (so `plot_mixture()` and `fit_mixture()` stay perfectly aligned).
#'
#' Supported inputs:
#' * `posterior::draws_*` -> coerced via `posterior::as_draws_df()`
#' * `data.frame` / `matrix` -> returned (or converted) as a data frame
#' * `list` of numeric vectors/matrices:
#'   - vector -> one column named by the list name (or `group{j}` if unnamed)
#'   - 1-col matrix -> one column named by the list name
#'   - multi-col matrix -> columns named `{list_name}_{1..p}`
#'
#' @param samples Posterior samples.
#' @return A `data.frame` with stable, deterministic column names.
#' @keywords internal
#' @noRd
.coerce_draws_df <- function(x) {
  # Case 1: Already a data frame
  if (is.data.frame(x)) {
    return(x)
  }
  
  # Case 2: List of vectors/matrices (common from extract_draws)
  if (is.list(x)) {
    all_numeric <- all(vapply(x, function(el) {
      is.numeric(el) || (is.matrix(el) && is.numeric(el))
    }, logical(1)))
    
    if (!all_numeric) {
      stop("List elements must be numeric vectors or matrices")
    }
    
    vectors <- lapply(x, function(el) {
      if (is.matrix(el)) {
        if (ncol(el) == 1) {
          as.vector(el)
        } else {
          stop("Each list element should be univariate (vector or single-column matrix)")
        }
      } else {
        as.vector(el)
      }
    })
    
    lengths <- vapply(vectors, length, integer(1))
    if (length(unique(lengths)) != 1) {
      stop("All posterior samples must have the same length. Found lengths: ",
           paste(unique(lengths), collapse = ", "))
    }
    
    df <- as.data.frame(vectors)
    
    if (!is.null(names(x)) && length(names(x)) == length(vectors)) {
      colnames(df) <- names(x)
    } else {
      colnames(df) <- paste0("var", seq_along(vectors))
    }
    
    return(df)
  }
  
  # Case 3: Matrix
  if (is.matrix(x)) {
    df <- as.data.frame(x)
    if (is.null(colnames(df))) {
      colnames(df) <- paste0("var", seq_len(ncol(df)))
    }
    return(df)
  }
  
  # Case 4: Single vector
  if (is.vector(x)) {
    stop("Single vector provided. fit_mixture() requires at least 2 variables for shrinkage.")
  }
  
  stop("Unsupported input type. Provide a data.frame, matrix, or named list of vectors/matrices.")
}

#' Internal: stop if fewer than 2 groups
#' @keywords internal
.stop_if_lt_two_groups <- function(groups) {
  G <- length(unique(groups))
  if (G < 2) {
    stop("[shrinkr] Shrinkage requires at least 2 groups; found ", G, ".")
  }
}

#' Prepare Stan data from mixture
#' @keywords internal
.prep_bhm_data_from_mixture <- function(mixture, pri, centered) {
  if (!is.null(mixture$weights) && !is.null(mixture$covs) && !is.null(mixture$means)) {
    w <- as.numeric(mixture$weights)
    if (any(w < 0)) stop("Mixture weights must be non-negative.")
    if (!isTRUE(all.equal(sum(w), 1, tol = 1e-8))) w <- w / sum(w)
    K <- length(w)
    
    if (is.list(mixture$means)) {
      m_mat <- do.call(rbind, mixture$means)
    } else {
      m_mat <- as.matrix(mixture$means)
      if (nrow(m_mat) != K) m_mat <- t(m_mat)
      if (nrow(m_mat) != K) stop("Couldn't coerce mixture means to K x G.")
    }
    G <- ncol(m_mat)
    L <- lapply(seq_len(K), function(k) .as_chol_factor(mixture$covs[[k]]))
    
  } else if (!is.null(mixture$components) && !is.null(mixture$covs)) {
    comp <- mixture$components
    if (!all(c("component", "variable", "mean", "weight") %in% names(comp))) {
      stop("Mixture components table lacks required columns: component/variable/mean/weight.")
    }
    K    <- length(unique(comp$component))
    vars <- if ("vars" %in% names(mixture) && !is.null(mixture$vars)) mixture$vars else unique(comp$variable)
    G    <- length(vars)
    
    w <- vapply(seq_len(K), function(k) {
      wk <- unique(comp$weight[comp$component == k])[1]
      as.numeric(wk)
    }, numeric(1))
    if (!isTRUE(all.equal(sum(w), 1, tol = 1e-8))) w <- w / sum(w)
    
    m_mat <- matrix(NA_real_, nrow = K, ncol = G)
    for (k in seq_len(K)) {
      mk <- comp[comp$component == k, ]
      mk <- mk[match(vars, mk$variable), ]
      m_mat[k, ] <- mk$mean
    }
    
    if (!is.list(mixture$covs) || length(mixture$covs) != K) {
      stop("mixture$covs must be a list of length K with each element a GxG covariance matrix.")
    }
    L <- lapply(mixture$covs, .as_chol_factor)
    
  } else {
    stop("Unsupported mixture structure. Provide either ($weights,$means,$covs) or ($components,$covs).")
  }
  
  m_arr <- lapply(seq_len(K), function(k) as.numeric(m_mat[k, ]))
  
  list(
    G = as.integer(ncol(m_mat)),
    K = as.integer(K),
    w = array(as.numeric(w), dim = length(w)),
    m = m_arr,
    L = L,
    
    mu_prior_type = as.integer(pri$mu_prior_type),
    mu_loc = pri$mu_loc,
    mu_scale = pri$mu_scale,
    mu_df = pri$mu_df,
    mu_n_components = as.integer(pri$mu_n_components),
    mu_mix_weights = array(pri$mu_mix_weights, dim = pri$mu_n_components),
    mu_mix_locs = array(pri$mu_mix_locs, dim = pri$mu_n_components),
    mu_mix_scales = array(pri$mu_mix_scales, dim = pri$mu_n_components),
    mu_lb = pri$mu_lb,
    mu_ub = pri$mu_ub,
    mu_is_truncated = as.integer(pri$mu_is_truncated),
    
    tau_prior_type = as.integer(pri$tau_prior_type),
    tau_params = pri$tau_params,
    tau_n_components = as.integer(pri$tau_n_components),
    tau_mix_weights = array(pri$tau_mix_weights, dim = pri$tau_n_components),
    tau_mix_locs = array(pri$tau_mix_locs, dim = pri$tau_n_components),
    tau_mix_scales = array(pri$tau_mix_scales, dim = pri$tau_n_components),
    tau_lb = pri$tau_lb,
    tau_ub = pri$tau_ub,
    
    custom_prior_type = as.integer(pri$custom_prior_type),
    custom_params = pri$custom_params,
    
    centered = as.integer(centered)
  )
}

#' Prepare Stan data from mle and cov
#' @keywords internal
.prep_bhm_data_from_mle <- function(mle, var_matrix, pri, centered) {
  if (!is.numeric(mle) || is.matrix(mle)) stop("`mle` must be a numeric vector.")
  G <- length(mle)
  
  if (is.vector(var_matrix) && !is.matrix(var_matrix)) {
    if (length(var_matrix) != G) stop("Variance vector length must match length(mle).")
    if (any(var_matrix <= 0))   stop("All variances must be positive.")
    Sigma <- diag(as.numeric(var_matrix))
  } else if (is.matrix(var_matrix)) {
    if (!all(dim(var_matrix) == c(G, G))) stop("`var_matrix` must be G x G.")
    if (!isTRUE(all.equal(var_matrix, t(var_matrix), tolerance = 1e-8))) {
      stop("`var_matrix` must be symmetric (covariance).")
    }
    eig <- eigen(var_matrix, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) <= 0) stop("`var_matrix` must be positive-definite.")
    Sigma <- var_matrix
  } else {
    stop("`var_matrix` must be a variance vector or covariance matrix.")
  }
  
  L1 <- .as_chol_factor(Sigma)
  
  list(
    G = as.integer(G),
    K = 1L,
    w = array(1.0, dim = 1L),
    m = list(as.numeric(mle)),
    L = list(L1),
    
    mu_prior_type = as.integer(pri$mu_prior_type),
    mu_loc = pri$mu_loc,
    mu_scale = pri$mu_scale,
    mu_df = pri$mu_df,
    mu_n_components = as.integer(pri$mu_n_components),
    mu_mix_weights = array(pri$mu_mix_weights, dim = pri$mu_n_components),
    mu_mix_locs = array(pri$mu_mix_locs, dim = pri$mu_n_components),
    mu_mix_scales = array(pri$mu_mix_scales, dim = pri$mu_n_components),
    mu_lb = pri$mu_lb,
    mu_ub = pri$mu_ub,
    mu_is_truncated = as.integer(pri$mu_is_truncated),
    
    tau_prior_type = as.integer(pri$tau_prior_type),
    tau_params = pri$tau_params,
    tau_n_components = as.integer(pri$tau_n_components),
    tau_mix_weights = array(pri$tau_mix_weights, dim = pri$tau_n_components),
    tau_mix_locs = array(pri$tau_mix_locs, dim = pri$tau_n_components),
    tau_mix_scales = array(pri$tau_mix_scales, dim = pri$tau_n_components),
    tau_lb = pri$tau_lb,
    tau_ub = pri$tau_ub,
    
    custom_prior_type = as.integer(pri$custom_prior_type),
    custom_params = pri$custom_params,
    
    centered = as.integer(centered)
  )
}

#' Internal: convert covariance to Cholesky factor for Stan
#'
#' Uses a jitter-and-retry strategy when the input matrix is not quite
#' positive-definite. Falls back to [Matrix::nearPD()] if available, and
#' errors rather than passing a non-lower-triangular factor to Stan.
#'
#' @keywords internal
.as_chol_factor <- function(Sigma) {
  # Try direct Cholesky first
  L <- tryCatch(
    t(chol(Sigma)),
    error = function(e) NULL
  )
  if (!is.null(L)) return(as.matrix(L))
  
  # Jitter-and-retry strategy
  jitters <- c(1e-10, 1e-8, 1e-6, 1e-4)
  n <- nrow(Sigma)
  for (jit in jitters) {
    L <- tryCatch(
      t(chol(Sigma + diag(jit, n))),
      error = function(e) NULL
    )
    if (!is.null(L)) {
      warning("Added jitter (", jit, ") to covariance matrix for Cholesky decomposition.",
              call. = FALSE)
      return(as.matrix(L))
    }
  }
  
  # Try nearPD if Matrix package is available
  if (requireNamespace("Matrix", quietly = TRUE)) {
    Sigma_pd <- tryCatch({
      as.matrix(Matrix::nearPD(Sigma, corr = FALSE, keepDiag = TRUE)$mat)
    }, error = function(e) NULL)
    if (!is.null(Sigma_pd)) {
      L <- tryCatch(
        t(chol(Sigma_pd)),
        error = function(e) NULL
      )
      if (!is.null(L)) {
        warning("Used nearPD() to make covariance positive-definite for Cholesky.",
                call. = FALSE)
        return(as.matrix(L))
      }
    }
  }
  
  stop("Could not compute Cholesky factor of covariance matrix. ",
       "The matrix may be severely ill-conditioned or indefinite.",
       call. = FALSE)
}

#' Extract mean and sd from a mixture component (possibly truncated normal)
#' @keywords internal
#' @noRd
.extract_normal_params <- function(comp) {
  cj <- if (inherits(comp, "distribution")) unclass(comp)[[1]] else comp
  # Unwrap truncation if present
  if (inherits(cj, "dist_truncated") && !is.null(cj$dist)) {
    cj <- cj$dist
  }
  if (!is.null(cj$mu) && !is.null(cj$sigma)) {
    return(list(mu = cj$mu, sigma = cj$sigma))
  }
  # Fallback via parameters()
  pj <- tryCatch(distributional::parameters(comp), error = function(e) NULL)
  if (!is.null(pj) && !is.null(pj$mu) && !is.null(pj$sigma)) {
    return(list(mu = pj$mu, sigma = pj$sigma))
  }
  list(mu = 0, sigma = 1)
}

#' Convert distributional priors to Stan format
#'
#' Handles Normal, Student-t, Cauchy, Gamma, Exponential, Lognormal,
#' Inverse-Gamma, Uniform, truncated distributions, and mixture priors
#' (including spike-and-slab) for both mu and tau.
#'
#' @keywords internal
.coerce_priors_to_stan <- function(prior_mu, prior_tau) {
  
  out <- list(
    # mu parameters
    mu_prior_type = 1,
    mu_loc = NA_real_,
    mu_scale = NA_real_,
    mu_df = NA_real_,
    mu_n_components = 1,
    mu_mix_weights = 1,
    mu_mix_locs = 0,
    mu_mix_scales = 10,
    mu_lb = -Inf,
    mu_ub = Inf,
    mu_is_truncated = 0,
    
    # tau parameters
    tau_prior_type = 1,
    tau_params = rep(NA_real_, 6),
    tau_lb = NA_real_,
    tau_ub = NA_real_,
    tau_n_components = 1,
    tau_mix_weights = 1,
    tau_mix_locs = 0,
    tau_mix_scales = 2.5,
    
    # custom prior parameters
    custom_prior_type = 0,
    custom_params = rep(0, 10)
  )
  
  # =========================================================================
  # Process mu prior
  # =========================================================================
  if (!inherits(prior_mu, "distribution")) {
    stop("mu prior must be a distributional object")
  }
  
  mu_raw <- unclass(prior_mu)[[1]]
  
  # Check if mu is truncated — now supported
  mu_is_truncated <- inherits(mu_raw, "dist_truncated")
  
  if (mu_is_truncated) {
    out$mu_is_truncated <- 1
    if (!is.null(mu_raw$lower) && is.finite(mu_raw$lower)) {
      out$mu_lb <- mu_raw$lower
    }
    if (!is.null(mu_raw$upper) && is.finite(mu_raw$upper)) {
      out$mu_ub <- mu_raw$upper
    }
    mu_raw <- mu_raw$dist
  }
  
  # Check for mixture prior on mu
  # distributional::dist_mixture uses $dist (list) and $w (weights)
  # Our old custom class used $components and $weights
  mu_is_mix <- inherits(mu_raw, "dist_mixture") ||
    (!is.null(mu_raw$dist) && !is.null(mu_raw$w) && is.list(mu_raw$dist)) ||
    (!is.null(mu_raw$components) && !is.null(mu_raw$weights))
  
  if (mu_is_mix) {
    out$mu_prior_type <- 3  # Mixture
    # Handle both distributional and legacy field names
    comps <- mu_raw$dist %||% mu_raw$components
    wts <- mu_raw$w %||% mu_raw$weights
    n_comp <- length(comps)
    out$mu_n_components <- n_comp
    
    mix_locs <- numeric(n_comp)
    mix_scales <- numeric(n_comp)
    for (j in seq_len(n_comp)) {
      pp <- .extract_normal_params(comps[[j]])
      mix_locs[j] <- pp$mu
      mix_scales[j] <- pp$sigma
    }
    out$mu_mix_weights <- wts
    out$mu_mix_locs <- mix_locs
    out$mu_mix_scales <- mix_scales
    out$mu_loc <- sum(wts * mix_locs)
    out$mu_scale <- max(mix_scales)
    out$mu_df <- 1
    
  } else {
    # Parse standard mu distribution
    mu_parsed <- FALSE
    
    if (inherits(mu_raw, "dist_student_t")) {
      out$mu_prior_type <- 2
      out$mu_df <- mu_raw$df
      out$mu_loc <- mu_raw$mu
      out$mu_scale <- mu_raw$sigma
      mu_parsed <- TRUE
    } else if (inherits(mu_raw, "dist_normal")) {
      out$mu_prior_type <- 1
      out$mu_loc <- mu_raw$mu
      out$mu_scale <- mu_raw$sigma
      mu_parsed <- TRUE
    }
    
    if (!mu_parsed) {
      stop("mu prior must be Normal, Student-t, or a mixture. Got: ", format(prior_mu)[1])
    }
  }
  
  # =========================================================================
  # Process tau prior
  # =========================================================================
  tau_raw <- unclass(prior_tau)[[1]]
  
  is_truncated <- inherits(tau_raw, "dist_truncated")
  
  if (is_truncated) {
    if (!is.null(tau_raw$lower)) out$tau_lb <- tau_raw$lower
    if (!is.null(tau_raw$upper)) out$tau_ub <- tau_raw$upper
    tau_raw <- tau_raw$dist
  }
  
  # Check for mixture prior on tau
  # distributional::dist_mixture uses $dist (list) and $w (weights)
  # Legacy custom class used $components and $weights
  tau_is_mix <- inherits(tau_raw, "dist_mixture") ||
    (!is.null(tau_raw$dist) && !is.null(tau_raw$w) && is.list(tau_raw$dist)) ||
    (!is.null(tau_raw$components) && !is.null(tau_raw$weights))
  
  if (tau_is_mix) {
    out$tau_prior_type <- 9  # Mixture
    comps <- tau_raw$dist %||% tau_raw$components
    wts <- tau_raw$w %||% tau_raw$weights
    n_comp <- length(comps)
    out$tau_n_components <- n_comp
    
    mix_locs <- numeric(n_comp)
    mix_scales <- numeric(n_comp)
    for (j in seq_len(n_comp)) {
      pp <- .extract_normal_params(comps[[j]])
      mix_locs[j] <- pp$mu
      mix_scales[j] <- pp$sigma
    }
    out$tau_mix_weights <- wts
    out$tau_mix_locs <- mix_locs
    out$tau_mix_scales <- mix_scales
    out$tau_params <- rep(0, 6)
    if (is.na(out$tau_lb)) out$tau_lb <- 0
    if (is.na(out$tau_ub)) out$tau_ub <- Inf
    out$tau_params[6] <- out$tau_lb
    tau_parsed <- TRUE
    
  } else {
    tau_parsed <- FALSE
  }
  
  if (!tau_parsed) {
    if (inherits(tau_raw, "dist_student_t") &&
        !is.null(tau_raw$mu) && !is.null(tau_raw$sigma) && !is.null(tau_raw$df)) {
      out$tau_prior_type <- 2
      out$tau_params[1] <- tau_raw$mu
      out$tau_params[2] <- tau_raw$sigma
      out$tau_params[3] <- tau_raw$df
      tau_parsed <- TRUE
    } else if (inherits(tau_raw, "dist_normal") &&
               !is.null(tau_raw$mu) && !is.null(tau_raw$sigma)) {
      out$tau_prior_type <- 1
      out$tau_params[1] <- tau_raw$mu
      out$tau_params[2] <- tau_raw$sigma
      tau_parsed <- TRUE
    } else if (inherits(tau_raw, "dist_cauchy") &&
               (!is.null(tau_raw$location) || !is.null(tau_raw$x0)) &&
               (!is.null(tau_raw$scale) || !is.null(tau_raw$gamma))) {
      out$tau_prior_type <- 8
      out$tau_params[1] <- tau_raw$location %||% tau_raw$x0 %||% 0
      out$tau_params[2] <- tau_raw$scale %||% tau_raw$gamma %||% 1
      tau_parsed <- TRUE
    } else if (inherits(tau_raw, "dist_exponential") && !is.null(tau_raw$rate)) {
      out$tau_prior_type <- 6
      out$tau_params[1] <- tau_raw$rate
      tau_parsed <- TRUE
    } else if (inherits(tau_raw, "dist_lognormal") &&
               !is.null(tau_raw$mu) && !is.null(tau_raw$sigma)) {
      out$tau_prior_type <- 3
      out$tau_params[4] <- tau_raw$mu
      out$tau_params[5] <- tau_raw$sigma
      tau_parsed <- TRUE
    } else if (inherits(tau_raw, "dist_inverse_gamma")) {
      # distributional::dist_inverse_gamma(shape, rate) stores the parameters
      # in fields named `s` (shape) and `r` (rate). Stan's
      # inv_gamma_lpdf(x | alpha, beta) uses (shape, rate) under the same
      # density (beta^alpha/Gamma(alpha)) * x^(-alpha-1) * exp(-beta/x), so
      # the values pass through unchanged.
      params <- tryCatch(
        distributional::parameters(tau_raw),
        error = function(e) NULL
      )
      if (!is.null(params) && !is.null(params$s) && !is.null(params$r)) {
        out$tau_prior_type <- 4
        out$tau_params[1] <- params$s
        out$tau_params[2] <- params$r
        tau_parsed <- TRUE
      }
    } else if (inherits(tau_raw, "dist_gamma") &&
               !is.null(tau_raw$shape) && !is.null(tau_raw$rate)) {
      out$tau_prior_type <- 5
      out$tau_params[1] <- tau_raw$shape
      out$tau_params[2] <- tau_raw$rate
      tau_parsed <- TRUE
    } else if (inherits(tau_raw, "dist_uniform") &&
               !is.null(tau_raw$l) && !is.null(tau_raw$u)) {
      # distributional::dist_uniform(min, max) stores bounds as `l` and `u`
      out$tau_prior_type <- 7
      out$tau_lb <- tau_raw$l
      out$tau_ub <- tau_raw$u
      tau_parsed <- TRUE
    }
  }
  
  if (!tau_parsed) {
    stop("Unsupported tau prior type: ", format(prior_tau)[1])
  }
  
  # Validate tau prior is appropriate (must be positive)
  distributions_needing_truncation <- c(1, 2, 8)
  
  if (out$tau_prior_type %in% distributions_needing_truncation && !is_truncated) {
    stop("tau prior must be truncated with lower bound >= 0. ",
         "Untruncated Normal, Student-t, and Cauchy priors can produce negative values, ",
         "which is not valid for a scale parameter. ",
         "Use dist_truncated(your_prior, lower = 0) to fix this.")
  }
  
  if (is_truncated && !is.na(out$tau_lb) && out$tau_lb < 0) {
    stop("tau prior lower bound must be >= 0. Got lower bound: ", out$tau_lb)
  }
  
  if (is.na(out$tau_lb)) out$tau_lb <- 0
  if (is.na(out$tau_ub)) out$tau_ub <- Inf
  
  if (is.na(out$mu_loc)) out$mu_loc <- 0
  if (is.na(out$mu_scale)) out$mu_scale <- 1
  if (is.na(out$mu_df)) out$mu_df <- 1
  
  out$tau_params[is.na(out$tau_params)] <- 0
  out$tau_params[6] <- out$tau_lb
  
  return(out)
}