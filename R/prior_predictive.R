#' Sample from prior predictive distribution
#'
#' @description
#' Generates samples from the prior predictive distribution for the hierarchical
#' shrinkage model. Useful for prior elicitation and sensitivity analysis.
#'
#' The generative process is:
#' \enumerate{
#'   \item Sample mu from p(mu)
#'   \item Sample tau from p(tau)
#'   \item Sample theta_i ~ N(mu, tau) for each group i
#' }
#'
#' @param hierarchical_priors Named list with `mu` and `tau` distributional objects
#'   from the `distributional` package.
#' @param n_groups Integer; number of subgroups (G).
#' @param n_draws Integer; number of prior predictive samples to draw. Default 1000.
#' @param group_names Optional character vector of length `n_groups` to label groups.
#'
#' @return A list with class "shrinkr_prior_pred" containing:
#'   \item{mu}{Vector of mu draws}
#'   \item{tau}{Vector of tau draws}
#'   \item{theta}{Matrix of theta draws (n_draws x n_groups)}
#'   \item{implied_range}{Vector of ranges (max - min) of theta across groups for each draw}
#'   \item{implied_sd}{Vector of standard deviations of theta across groups for each draw}
#'   \item{group_names}{Group labels}
#'   \item{n_draws}{Number of draws}
#'   \item{n_groups}{Number of groups}
#'   \item{priors}{The hierarchical_priors specification used}
#'
#' @seealso
#' \code{\link{shrink}} for fitting the hierarchical model,
#' \code{\link{plot.shrinkr_prior_pred}} for visualizing prior predictive samples
#'
#' @export
#' 
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' median(prior_pred$implied_range)
#' head(as.data.frame(prior_pred))
sample_prior_predictive <- function(hierarchical_priors,
                                    n_groups,
                                    n_draws = 1000,
                                    group_names = NULL
) {
  
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  # Validate inputs
  if (!is.list(hierarchical_priors) || 
      !all(c("mu", "tau") %in% names(hierarchical_priors))) {
    cli::cli_abort("{.arg hierarchical_priors} must be a list with elements {.field mu} and {.field tau}")
  }
  
  if (!inherits(hierarchical_priors$mu, "distribution") || 
      !inherits(hierarchical_priors$tau, "distribution")) {
    cli::cli_abort("{.field hierarchical_priors$mu} and {.field hierarchical_priors$tau} must be distributional objects")
  }
  
  if (!is.numeric(n_groups) || length(n_groups) != 1 || n_groups < 1) {
    cli::cli_abort("{.arg n_groups} must be a positive integer")
  }
  n_groups <- as.integer(n_groups)
  
  if (!is.numeric(n_draws) || length(n_draws) != 1 || n_draws < 1) {
    cli::cli_abort("{.arg n_draws} must be a positive integer")
  }
  n_draws <- as.integer(n_draws)
  
  
  # Determine group names
  if (is.null(group_names)) {
    group_names <- paste0("group", seq_len(n_groups))
  } else if (length(group_names) != n_groups) {
    cli::cli_abort("{.arg group_names} must have length n_groups = {n_groups}")
  }
  
  # Sample mu from prior
  mu_draws <- distributional::generate(hierarchical_priors$mu, n_draws)
  # Extract vector from list if needed
  if (is.list(mu_draws)) {
    mu_draws <- mu_draws[[1]]
  }
  if (is.matrix(mu_draws)) {
    mu_draws <- as.vector(mu_draws)
  }
  
  # Sample tau from prior
  tau_draws <- distributional::generate(hierarchical_priors$tau, n_draws)
  # Extract vector from list if needed
  if (is.list(tau_draws)) {
    tau_draws <- tau_draws[[1]]
  }
  if (is.matrix(tau_draws)) {
    tau_draws <- as.vector(tau_draws)
  }
  
  # Validate tau draws
  if (!is.numeric(tau_draws) || any(!is.finite(tau_draws))) {
    cli::cli_abort("tau prior produced invalid samples (non-numeric or non-finite)")
  }
  
  # tau is a scale parameter — must be strictly positive
  if (any(tau_draws <= 0)) {
    n_neg <- sum(tau_draws <= 0)
    cli::cli_abort(
      c("x" = "{n_neg}/{n_draws} tau draws were <= 0.",
        "i" = "tau is a scale parameter and must be positive.",
        "i" = "Truncate your prior with {.code dist_truncated(..., lower = 0)}.",
        "i" = "For spike-and-slab: {.code dist_truncated(prior_spike_slab(...), lower = 0)}")
    )
  }
  
  # Sample theta ~ N(mu, tau) for each group
  theta_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_groups)
  for (i in seq_len(n_draws)) {
    theta_matrix[i, ] <- rnorm(n_groups, mean = mu_draws[i], sd = tau_draws[i])
  }
  colnames(theta_matrix) <- group_names
  
  # Compute implied summaries across groups for each draw
  implied_range <- apply(theta_matrix, 1, function(row) diff(range(row)))
  implied_sd <- apply(theta_matrix, 1, sd)
  
  # Create result object
  result <- list(
    mu = mu_draws,
    tau = tau_draws,
    theta = theta_matrix,
    implied_range = implied_range,
    implied_sd = implied_sd,
    group_names = group_names,
    n_draws = n_draws,
    n_groups = n_groups,
    priors = hierarchical_priors
  )
  
  class(result) <- c("shrinkr_prior_pred", "list")
  result
}

#' Convert prior predictive samples to data frame
#'
#' @description
#' Converts prior predictive samples to a tidy long-format data frame suitable
#' for analysis and visualization with tidyverse tools.
#'
#' @param x A `shrinkr_prior_pred` object from \code{\link{sample_prior_predictive}}.
#' @param row.names Ignored (for S3 consistency).
#' @param optional Ignored (for S3 consistency).
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame (or tibble if tibble package is available) with columns:
#'   \item{.draw}{Draw number (1 to n_draws)}
#'   \item{group}{Group name}
#'   \item{theta}{Sampled group-level effect}
#'   \item{mu}{Sampled global mean for this draw}
#'   \item{tau}{Sampled heterogeneity parameter for this draw}
#'
#' @seealso \code{\link{sample_prior_predictive}} for generating prior predictive samples
#'
#' @export
#'
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' df <- as.data.frame(prior_pred)
#' head(df)
as.data.frame.shrinkr_prior_pred <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot(inherits(x, "shrinkr_prior_pred"))
  
  # Expand theta to long format
  n_draws <- x$n_draws
  n_groups <- x$n_groups
  
  result <- data.frame(
    .draw = rep(seq_len(n_draws), each = n_groups),
    group = rep(x$group_names, times = n_draws),
    theta = as.vector(t(x$theta)),
    mu = rep(x$mu, each = n_groups),
    tau = rep(x$tau, each = n_groups),
    stringsAsFactors = FALSE
  )
  
  # Convert to tibble if available
  if (requireNamespace("tibble", quietly = TRUE)) {
    result <- tibble::as_tibble(result)
  }
  
  result
}


#' Plot prior predictive samples
#'
#' @description
#' Visualizes the prior predictive distribution for hyperparameters (mu and tau)
#' and subgroup effects (theta). Requires the ggplot2 package.
#'
#' @param x A `shrinkr_prior_pred` object from \code{\link{sample_prior_predictive}}.
#' @param type Character; type of plot. Options:
#'   \itemize{
#'     \item `"both"` - Both hyperparameter and theta plots (default)
#'     \item `"hyperparameters"` - Density plots for mu and tau only
#'     \item `"theta"` - Violin plots for theta by group only
#'   }
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object, or a list of two ggplot2 objects if `type = "both"`
#'   and patchwork package is not available. If patchwork is available and
#'   `type = "both"`, returns a combined plot.
#'
#' @seealso
#' \code{\link{sample_prior_predictive}} for generating samples,
#' \code{\link{as.data.frame.shrinkr_prior_pred}} for extracting data
#'
#' @export
#' @method plot shrinkr_prior_pred
#'
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_student_t(3, 0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' plot(prior_pred, type = "hyperparameters")
plot.shrinkr_prior_pred <- function(x, type = c("both", "hyperparameters", "theta"), ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting. Install with: install.packages('ggplot2')")
  }
  
  stopifnot(inherits(x, "shrinkr_prior_pred"))
  type <- match.arg(type)
  
  # Hyperparameter plot
  if (type %in% c("both", "hyperparameters")) {
    hyper_df <- data.frame(
      parameter = rep(c("mu", "tau"), each = x$n_draws),
      value = c(x$mu, x$tau)
    )
    
    p_hyper <- ggplot2::ggplot(hyper_df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "skyblue", alpha = 0.6) +
      ggplot2::facet_wrap(~parameter, scales = "free") +
      ggplot2::labs(
        title = "Prior predictive: Hyperparameters",
        x = "Value",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
  }
  
  # Theta plot
  if (type %in% c("both", "theta")) {
    theta_df <- as.data.frame(x)
    
    p_theta <- ggplot2::ggplot(theta_df, ggplot2::aes(x = theta, y = group)) +
      ggplot2::geom_violin(fill = "lightcoral", alpha = 0.6) +
      ggplot2::labs(
        title = "Prior predictive: Subgroup effects (theta)",
        x = "Effect size",
        y = "Group"
      ) +
      ggplot2::theme_minimal()
  }
  
  # Return appropriate plot(s)
  if (type == "both") {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(p_hyper / p_theta)
    } else {
      message("Install 'patchwork' for combined plots. Returning list instead.")
      return(list(hyperparameters = p_hyper, theta = p_theta))
    }
  } else if (type == "hyperparameters") {
    return(p_hyper)
  } else {
    return(p_theta)
  }
}


#' Print method for prior predictive samples
#'
#' @description
#' Displays summary information about prior predictive samples in a readable format.
#'
#' @param x A `shrinkr_prior_pred` object from \code{\link{sample_prior_predictive}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @seealso
#' \code{\link{sample_prior_predictive}} for generating samples,
#' \code{\link{summary.shrinkr_prior_pred}} for detailed summaries
#'
#' @export
#'
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' print(prior_pred)
print.shrinkr_prior_pred <- function(x, ...) {
  stopifnot(inherits(x, "shrinkr_prior_pred"))
  
  cat("== Prior Predictive Samples =========================\n\n")
  cat("Draws:  ", x$n_draws, "\n")
  cat("Groups: ", x$n_groups, "\n\n")
  
  cat("Prior specifications:\n")
  cat("  mu:  ", format(x$priors$mu), "\n")
  cat("  tau: ", format(x$priors$tau), "\n\n")
  
  cat("Hyperparameter summaries:\n")
  cat("  mu:  mean =", round(mean(x$mu), 3), 
      ", sd =", round(sd(x$mu), 3),
      ", range = [", round(min(x$mu), 3), ",", round(max(x$mu), 3), "]\n")
  cat("  tau: mean =", round(mean(x$tau), 3), 
      ", sd =", round(sd(x$tau), 3),
      ", range = [", round(min(x$tau), 3), ",", round(max(x$tau), 3), "]\n\n")
  
  cat("Theta (by group):\n")
  theta_means <- colMeans(x$theta)
  theta_sds <- apply(x$theta, 2, sd)
  for (i in seq_len(min(5, x$n_groups))) {
    cat("  ", x$group_names[i], ": mean =", round(theta_means[i], 3),
        ", sd =", round(theta_sds[i], 3), "\n")
  }
  if (x$n_groups > 5) {
    cat("  ... and", x$n_groups - 5, "more groups\n")
  }
  
  cat("\n-----------------------------------------------------\n")
  cat("Use plot() to visualize\n")
  cat("Use as.data.frame() for tidy format\n")
  cat("Use summary() for detailed statistics\n")
  
  invisible(x)
}


#' Summary statistics for prior predictive samples
#'
#' @description
#' Computes comprehensive summary statistics for hyperparameters (mu and tau)
#' and theta parameters from prior predictive samples.
#'
#' @param object A `shrinkr_prior_pred` object from \code{\link{sample_prior_predictive}}.
#' @param probs Numeric vector of quantiles to compute. Default is 
#'   c(0.025, 0.5, 0.975) for 95% credible intervals.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with class "summary.shrinkr_prior_pred" containing:
#'   \item{hyperparameters}{Data frame with summary statistics for mu and tau}
#'   \item{theta}{Data frame with summary statistics for each group's theta}
#'   
#'   Each data frame includes columns for mean, sd, and the requested quantiles.
#'
#' @seealso
#' \code{\link{sample_prior_predictive}} for generating samples,
#' \code{\link{print.shrinkr_prior_pred}} for quick overview
#' @method summary shrinkr_prior_pred
#' @export
#'
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' summ <- summary(prior_pred)
#' summ$theta
summary.shrinkr_prior_pred <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  stopifnot(inherits(object, "shrinkr_prior_pred"))
  
  # Validate probs
  if (!is.numeric(probs) || any(probs < 0) || any(probs > 1)) {
    stop("`probs` must be numeric values between 0 and 1")
  }
  
  # Summarize hyperparameters
  hyper_summary <- data.frame(
    parameter = c("mu", "tau"),
    mean = c(mean(object$mu), mean(object$tau)),
    sd = c(sd(object$mu), sd(object$tau)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Add quantiles
  for (p in probs) {
    col_name <- paste0("q", formatC(p * 100, format = "f", digits = 1))
    hyper_summary[[col_name]] <- c(
      quantile(object$mu, probs = p),
      quantile(object$tau, probs = p)
    )
  }
  
  # Summarize theta
  theta_summary <- data.frame(
    group = object$group_names,
    mean = colMeans(object$theta),
    sd = apply(object$theta, 2, sd),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Add quantiles for theta
  for (p in probs) {
    col_name <- paste0("q", formatC(p * 100, format = "f", digits = 1))
    theta_summary[[col_name]] <- apply(object$theta, 2, quantile, probs = p)
  }
  
  # Convert to tibble if available
  if (requireNamespace("tibble", quietly = TRUE)) {
    hyper_summary <- tibble::as_tibble(hyper_summary)
    theta_summary <- tibble::as_tibble(theta_summary)
  }
  
  result <- list(
    hyperparameters = hyper_summary,
    theta = theta_summary,
    n_draws = object$n_draws,
    n_groups = object$n_groups
  )
  
  class(result) <- "summary.shrinkr_prior_pred"
  result
}


#' Print summary of prior predictive samples
#'
#' @description
#' Displays formatted summary statistics for prior predictive samples.
#'
#' @param x A summary object from \code{\link{summary.shrinkr_prior_pred}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @seealso \code{\link{summary.shrinkr_prior_pred}} for creating summaries
#'
#' @export
#'
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' print(summary(prior_pred))
print.summary.shrinkr_prior_pred <- function(x, digits = 3, ...) {
  stopifnot(inherits(x, "summary.shrinkr_prior_pred"))
  
  cat("== Prior Predictive Summary =========================\n\n")
  cat("Based on", x$n_draws, "draws for", x$n_groups, "groups\n\n")
  
  cat("Hyperparameters:\n")
  print(x$hyperparameters, digits = digits, row.names = FALSE)
  
  cat("\nSubgroup effects (theta):\n")
  print(x$theta, digits = digits, row.names = FALSE)
  
  cat("\n-----------------------------------------------------\n")
  
  invisible(x)
}


#' Compute prior predictive pairwise differences |theta_i - theta_j|
#'
#' @description
#' Computes the prior-implied distribution of absolute pairwise differences
#' between subgroup effects. This is useful for calibrating priors: if your
#' prior implies that subgroup differences of 5 units are common but clinical
#' relevance starts at 0.5, your prior may be too diffuse.
#'
#' This implements the recommendation from the SDSIH Vignettes Library:
#' inspect the prior distribution of \eqn{|\theta_i - \theta_j|} when
#' choosing priors for Bayesian hierarchical models.
#'
#' @param prior_pred A `shrinkr_prior_pred` object from
#'   \code{\link{sample_prior_predictive}}.
#'
#' @return A list with class `"shrinkr_prior_contrasts"` containing:
#'   \item{differences}{Data frame with columns `pair`, `abs_diff`, and `.draw`}
#'   \item{summary}{Data frame with per-pair summary statistics}
#'   \item{overall_summary}{Named numeric vector of quantiles across all pairs}
#'   \item{n_pairs}{Number of unique pairs}
#'   \item{n_draws}{Number of prior predictive draws}
#'   \item{group_names}{Group labels used}
#'
#' @seealso
#' \code{\link{sample_prior_predictive}} for generating prior predictive samples,
#' \code{\link{plot.shrinkr_prior_contrasts}} for visualizing the result
#'
#' @export
#'
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' pw <- prior_pairwise_differences(prior_pred)
#' print(pw)
#' pw$overall_summary
prior_pairwise_differences <- function(prior_pred) {
  stopifnot(inherits(prior_pred, "shrinkr_prior_pred"))
  
  theta <- prior_pred$theta
  gnames <- prior_pred$group_names
  n_groups <- prior_pred$n_groups
  n_draws <- prior_pred$n_draws
  
  if (n_groups < 2) {
    stop("Need at least 2 groups for pairwise differences.")
  }
  
  # Generate all unique pairs
  pairs <- utils::combn(n_groups, 2)
  n_pairs <- ncol(pairs)
  
  # Compute |theta_i - theta_j| for each pair and draw
  diff_list <- vector("list", n_pairs)
  for (p in seq_len(n_pairs)) {
    i <- pairs[1, p]
    j <- pairs[2, p]
    pair_label <- paste0(gnames[i], " vs ", gnames[j])
    abs_diffs <- abs(theta[, i] - theta[, j])
    diff_list[[p]] <- data.frame(
      pair = pair_label,
      abs_diff = abs_diffs,
      .draw = seq_len(n_draws),
      stringsAsFactors = FALSE
    )
  }
  
  differences <- do.call(rbind, diff_list)
  
  # Summary per pair
  pair_labels <- unique(differences$pair)
  summary_df <- do.call(rbind, lapply(pair_labels, function(pl) {
    d <- differences$abs_diff[differences$pair == pl]
    data.frame(
      pair = pl,
      mean = mean(d),
      sd = sd(d),
      median = stats::median(d),
      q2.5 = stats::quantile(d, 0.025),
      q25 = stats::quantile(d, 0.25),
      q75 = stats::quantile(d, 0.75),
      q97.5 = stats::quantile(d, 0.975),
      prob_gt_0.5 = mean(d > 0.5),
      prob_gt_1 = mean(d > 1),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
  
  # Overall summary across all pairs
  all_diffs <- differences$abs_diff
  overall <- stats::quantile(all_diffs, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(overall) <- paste0("q", c(2.5, 25, 50, 75, 97.5))
  
  if (requireNamespace("tibble", quietly = TRUE)) {
    differences <- tibble::as_tibble(differences)
    summary_df <- tibble::as_tibble(summary_df)
  }
  
  result <- list(
    differences = differences,
    summary = summary_df,
    overall_summary = overall,
    n_pairs = n_pairs,
    n_draws = n_draws,
    group_names = gnames
  )
  
  class(result) <- c("shrinkr_prior_contrasts", "list")
  result
}


#' Print method for prior pairwise contrasts
#'
#' @param x A `shrinkr_prior_contrasts` object.
#' @param digits Number of digits to display. Default 3.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.shrinkr_prior_contrasts <- function(x, digits = 3, ...) {
  stopifnot(inherits(x, "shrinkr_prior_contrasts"))
  
  cat("== Prior Predictive: Pairwise |theta_i - theta_j| ==\n\n")
  cat("Groups: ", length(x$group_names), "\n")
  cat("Pairs:  ", x$n_pairs, "\n")
  cat("Draws:  ", x$n_draws, "\n\n")
  
  cat("Overall quantiles of |theta_i - theta_j|:\n")
  cat("  ", paste(names(x$overall_summary), "=",
                  round(x$overall_summary, digits), collapse = ", "), "\n\n")
  
  cat("Per-pair summary:\n")
  print(x$summary[, c("pair", "median", "q2.5", "q97.5", "prob_gt_0.5", "prob_gt_1")],
        digits = digits, row.names = FALSE)
  
  cat("\n-----------------------------------------------------\n")
  cat("Use plot() to visualize\n")
  
  invisible(x)
}


#' Plot prior predictive pairwise differences
#'
#' @description
#' Creates a density plot of \eqn{|\theta_i - \theta_j|}
#' from prior predictive samples. Useful for calibrating hierarchical priors.
#' Styling matches [plot.shrinkr_prior_pred()] for visual consistency.
#'
#' @param x A `shrinkr_prior_contrasts` object from
#'   \code{\link{prior_pairwise_differences}}.
#' @param by_pair Logical; if `TRUE`, facet by pair. If `FALSE` (default),
#'   pool all pairwise differences into a single plot.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object.
#'
#' @seealso \code{\link{prior_pairwise_differences}}
#'
#' @export
#' @method plot shrinkr_prior_contrasts
#'
#' @examples
#' priors <- list(
#'   mu = distributional::dist_normal(0, 5),
#'   tau = distributional::dist_truncated(distributional::dist_normal(0, 1), lower = 0)
#' )
#' prior_pred <- sample_prior_predictive(priors, n_groups = 3, n_draws = 50)
#' pw <- prior_pairwise_differences(prior_pred)
#' plot(pw)
plot.shrinkr_prior_contrasts <- function(x,
                                         by_pair = FALSE,
                                         ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting. Install with: install.packages('ggplot2')")
  }
  stopifnot(inherits(x, "shrinkr_prior_contrasts"))
  
  df <- x$differences
  
  if (by_pair && x$n_pairs > 1) {
    # By-pair: violin-style density per pair, matching theta plot (lightcoral)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = abs_diff, y = pair)) +
      ggplot2::geom_violin(fill = "lightcoral", alpha = 0.6) +
      ggplot2::labs(
        title = expression("Prior predictive: Pairwise " * "|" * theta[i] - theta[j] * "|"),
        x = expression("|" * theta[i] - theta[j] * "|"),
        y = "Pair"
      ) +
      ggplot2::theme_minimal()
  } else {
    # Pooled: single density, matching hyperparameter plot (skyblue)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = abs_diff)) +
      ggplot2::geom_density(fill = "skyblue", alpha = 0.6) +
      ggplot2::labs(
        title = expression("Prior predictive: Pairwise " * "|" * theta[i] - theta[j] * "|"),
        x = expression("|" * theta[i] - theta[j] * "|"),
        y = "Density"
      ) +
      ggplot2::theme_minimal()
  }
  
  p
}
