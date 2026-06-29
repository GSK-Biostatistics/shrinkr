#' Fit Gaussian mixture models to posterior samples
#'
#' @description
#' Fits a multivariate Gaussian mixture model (GMM) jointly across all supplied
#' variables using \pkg{mclust}. The function is intended for approximating
#' posterior draws from Bayesian models and produces a tidy component table
#' suitable for visualization and shrinkage modeling.
#'

#' @param samples Posterior samples in one of the following formats:
#'   \itemize{
#'     \item **Data frame or matrix** (recommended): Each column represents one 
#'           group/variable to shrink, with rows as posterior draws. Column names 
#'           are used as variable labels. Example: output from 
#'           `posterior::as_draws_df()` or a matrix where columns = groups.
#'     \item **Named list of vectors**: Each list element contains posterior 
#'           samples for one group (as a numeric vector). Names are used as 
#'           variable labels. All vectors must have the same length.
#'     \item **Named list of matrices**: Each list element is a matrix of 
#'           posterior samples for one group. For univariate parameters, these 
#'           should be single-column matrices (n × 1). All matrices must have 
#'           the same number of rows (draws).
#'   }
#'   At least **two** groups/variables are required for hierarchical shrinkage. 
#'   Non-numeric columns in data frames are automatically dropped. Rows with 
#'   missing values are removed before fitting.
#' @param K_max Integer \eqn{(\ge 1)}. Maximum number of mixture components to consider
#'   during model selection by BIC. Internally capped at `n - 1` for stability.
#' @param verbose Logical. If `TRUE`, progress and diagnostic messages are printed.
#' @param model_names Optional character vector of \pkg{mclust} covariance model
#'   codes to consider (e.g., `"EII"`, `"VVV"`, etc.). If `NULL` (default), all
#'   models appropriate for the data dimension are considered by \pkg{mclust}.
#'   This is a convenience wrapper for \code{modelNames} in \code{mclust::Mclust()}.
#' @param ... Additional arguments forwarded to \code{mclust::Mclust()}, such as
#'   \code{prior}, \code{initialization}, \code{control = mclust::emControl()},
#'   \code{warn}, or \code{verbose}. Use these to fine-tune EM control, priors,
#'   or initialization.
#'
#' @returns A list with class "shrinkr_mixture" containing:
#' \itemize{
#'   \item `components` — data frame with columns: `group`, `component`, `variable`,
#'         `weight`, `mean`, and `sd` (marginal standard deviations).
#'   \item `K` — number of mixture components selected.
#'   \item `vars` — character vector of variable names.
#'   \item `weights` — vector of component weights (mixing proportions).
#'   \item `covs` — list of component covariance matrices (p × p each).
#'   \item `model_name` — selected mclust covariance structure (e.g., "VVV").
#'   \item `bic` — Bayesian Information Criterion for the fitted model.
#'   \item `n_samples` — number of samples used in fitting.
#'   \item `n_vars` — number of variables.
#'   \item `mclust_fit` — the complete mclust model object (for advanced use).
#'   \item `diagnostics` — list with sample size details, removed rows, and
#'         quality warnings.
#' }
#'
#' @details
#' - Requires at least **two** variables; shrinkage across a single variable is not meaningful.
#' - Rows with any missing values are removed before fitting.
#' - Component weights are normalized to sum to one.
#' - The `sd` column reports marginal standard deviations (square roots of
#'   diagonal entries) from each component covariance matrix.
#'
#' @section Covariance structures in \pkg{mclust}:
#' \pkg{mclust} parameterizes component covariances via eigen-decomposition and
#' offers a set of model families controlling volume (V), shape (S), and
#' orientation (O). Common codes include (non-exhaustive):
#'
#' \describe{
#'   \item{\strong{Spherical}}{
#'     \code{"EII"}: equal volume, spherical \cr
#'     \code{"VII"}: variable volume, spherical
#'   }
#'   \item{\strong{Diagonal}}{
#'     \code{"EEI"}: equal volume & shape (axis-aligned) \cr
#'     \code{"VEI"}: variable volume, equal shape \cr
#'     \code{"EVI"}: equal volume, variable shape \cr
#'     \code{"VVI"}: variable volume & shape
#'   }
#'   \item{\strong{Ellipsoidal (full covariance)}}{
#'     \code{"EEE"}: equal volume, shape, orientation \cr
#'     \code{"EEV"}: equal volume & shape, variable orientation \cr
#'     \code{"VEV"}: variable volume, equal shape, variable orientation \cr
#'     \code{"VVV"}: variable volume, shape, and orientation (most flexible)
#'   }
#' }
#'
#' If \code{model_names = NULL} (default), \code{mclust::Mclust()} selects among
#' the models appropriate for the data dimension via BIC. In practice, this lets
#' the data decide between parsimonious structures (e.g., \code{"EII"}, \code{"VVI"})
#' and fully flexible ones (e.g., \code{"VVV"}). You can restrict or expand the
#' search space by supplying \code{model_names}.
#'
#' @seealso
#' \code{\link{plot.shrinkr_mixture}} for visualizing marginal fits,
#' \code{\link{as.data.frame.shrinkr_mixture}} for extracting component data
#'
#' @examples
#' set.seed(1)
#' samples <- list(
#'   group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1),
#'   group3 = matrix(rnorm(100, 1.0, 0.5), ncol = 1)
#' )
#' mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
#' summary(mix)
#'
#' @import mclust
#' @importFrom stats cov complete.cases var
#' @export
fit_mixture <- function(samples,
                        K_max = 5L,
                        verbose = FALSE,
                        model_names = NULL,
                        ...) {
  
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  # --- Input checks ----------------------------------------------------------
  if (length(K_max) != 1L || !is.finite(K_max) ||
      K_max < 1L || K_max != as.integer(K_max)) {
    cli::cli_abort("{.arg K_max} must be a single integer >= 1")
  }
  
  if (verbose) cli::cli_alert_info("Coercing inputs...")
  X <- .coerce_draws_df(samples)
  stopifnot(is.data.frame(X))
  
  # Require >= 2 variables
  if (ncol(X) < 2L) {
    cli::cli_abort("{.fn fit_mixture} requires at least two variables; shrinkage is undefined with one")
  }
  
  # Keep only numeric columns
  num_cols <- vapply(X, is.numeric, logical(1))
  if (!all(num_cols)) {
    X <- X[, num_cols, drop = FALSE]
    if (verbose) cli::cli_alert_info("Dropped non-numeric columns")
  }
  if (ncol(X) < 2L) {
    cli::cli_abort("After dropping non-numeric columns, fewer than two variables remain")
  }
  
  # Drop rows with any NA
  keep <- stats::complete.cases(X)
  n_removed <- sum(!keep)
  if (n_removed) X <- X[keep, , drop = FALSE]
  n <- nrow(X)
  if (n < 2L) cli::cli_abort("Not enough complete rows after removing NA (need >= 2)")
  
  vars <- colnames(X)
  p    <- ncol(X)
  
  if (!requireNamespace("mclust", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg mclust} is required. Install with: {.code install.packages('mclust')}")
  }
  
  # Cap number of components by n - 1
  G_max_eff <- max(1L, min(as.integer(K_max), n - 1L))
  if (G_max_eff < K_max && verbose) {
    cli::cli_alert_info("Capping {.arg K_max} from {K_max} to {G_max_eff} (n = {n})")
  }
  G_range <- 1:G_max_eff
  
  # Heuristic warning for ill-conditioning
  if (n < p && verbose) {
    cli::cli_alert_warning("n ({n}) < p ({p}); covariances may be ill-conditioned")
  }
  
  # --- Fit multivariate GMM ---------------------------------------------------
  if (verbose) cli::cli_alert_info("Fitting multivariate Gaussian mixture...")
  # map convenience argument to Mclust's `modelNames`
  mnames <- model_names
  fit <- mclust::Mclust(
    data = as.matrix(X),
    G = G_range,
    modelNames = mnames,
    ...
  )
  
  K   <- fit$G
  pis <- as.numeric(fit$parameters$pro); pis <- pis / sum(pis)
  mus <- as.matrix(fit$parameters$mean)
  
  # Normalize sigma to p x p x K array (handle K = 1)
  sig <- fit$parameters$variance$sigma
  if (is.list(sig)) sig <- simplify2array(sig)
  if (length(dim(sig)) == 2L) sig <- array(sig, dim = c(p, p, 1L))
  
  # Post-fit singularity check
  singular_flags <- logical(dim(sig)[3])
  for (k in seq_len(dim(sig)[3])) {
    ev <- tryCatch(
      eigen(sig[, , k], symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) rep(NA_real_, p)
    )
    singular_flags[k] <- any(!is.finite(ev)) || min(ev, na.rm = TRUE) < 1e-8
  }
  if (any(singular_flags) && verbose) {
    cli::cli_alert_warning("Near-singular component covariance(s): {paste(which(singular_flags), collapse = ', ')}")
  }
  
  cov_list <- lapply(seq_len(dim(sig)[3]), function(k) sig[, , k, drop = FALSE][, , 1])
  
  # --- Build long component table --------------------------------------------
  components_df <- do.call(rbind, lapply(seq_len(K), function(k) {
    sd_k <- sqrt(diag(cov_list[[k]]))
    data.frame(
      component = k,
      variable  = vars,
      weight    = pis[k],
      mean      = mus[, k],
      sd        = sd_k,
      row.names = NULL, check.names = FALSE
    )
  }))
  
  # --- compute quantiles for credible interval  --------------------------------------------
  q_probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  quantiles_df <- as.data.frame(t(apply(X, 2, quantile, probs = q_probs)))
  colnames(quantiles_df) <- paste0("q", q_probs * 100)
  quantiles_df$variable <- vars
  quantiles_df <- quantiles_df[, c("variable", paste0("q", q_probs * 100))]
  
  
  # --- Build clean output structure ------------------------------------------
  result <- list(
    # Main outputs (user-facing)
    components = components_df,
    K = as.integer(K),
    vars = vars,
    weights = pis,
    covs = cov_list,
    quantiles = quantiles_df,
    
    # Model info
    model_name = fit$modelName,
    bic = fit$bic,
    loglik = fit$loglik,
    
    # Sample info
    n_samples = n,
    n_vars = p,
    
    # Advanced/technical
    mclust_fit = fit,
    diagnostics = list(
      n_removed_na = n_removed,
      K_max_requested = as.integer(K_max),
      K_max_effective = G_max_eff,
      singular_components = which(singular_flags),
      model_searched = if (is.null(mnames)) "all" else mnames
    )
  )
  
  class(result) <- c("shrinkr_mixture", "list")
  result
}
