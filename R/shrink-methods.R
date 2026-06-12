# shrink-methods.R

#' Print method for shrinkr_fit
#'
#' @description
#' Displays a compact summary of the fitted model including dimensions,
#' hyperparameter estimates, and diagnostics.
#'
#' @param x A `shrinkr_fit` object.
#' @param digits Number of digits to display. Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @seealso
#' [shrink()] for fitting models,
#' [summarise_theta()] for detailed theta estimates
#'
#' @export
print.shrinkr_fit <- function(x, digits = 3, ...) {
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(x, "shrinkr_fit"))
  
  cli::cli_rule("Shrinkr Fit")
  cli::cli_text("")
  
  cli::cli_h3("Model dimensions")
  cli::cli_ul(c(
    "Groups (G): {x$data$G}",
    "Mixture components: {x$data$K}",
    "Parameterization: {if (isTRUE(x$data$centered)) 'centered' else 'non-centered'}"
  ))
  cli::cli_text("")
  
  # Hyperparameters summary
  cli::cli_h3("Hyperparameters (mu, tau)")
  mt_summary <- tryCatch({
    draws <- extract_mu_tau(x)
    posterior::summarise_draws(
      draws,
      "mean", 
      "sd",
      ~posterior::quantile2(., probs = c(0.025, 0.5, 0.975)),
      "rhat"
    )
  }, error = function(e) NULL)
  
  if (!is.null(mt_summary)) {
    print(mt_summary, digits = digits)
  } else {
    cli::cli_alert_warning("Could not extract summaries")
  }
  
  # Diagnostics
  cli::cli_text("")
  cli::cli_h3("Diagnostics")
  cli::cli_ul(c(
    "Divergences: {x$diagnostics$n_divergent}",
    "Max tree depth: {x$diagnostics$max_treedepth}",
    "Max leapfrog: {x$diagnostics$n_leapfrog}"
  ))
  
  cli::cli_text("")
  cli::cli_rule()
  cli::cli_alert_info("Use {.fn summary} for detailed theta estimates")
  cli::cli_alert_info("Use {.fn as.data.frame} for full posterior draws")
  cli::cli_alert_info("Use {.fn plot} to visualize shrinkage")
  
  invisible(x)
}


#' Summary method for shrinkr_fit
#'
#' @description
#' Comprehensive posterior summary including hyperparameters and all
#' subgroup effects with convergence diagnostics.
#'
#' @param object A `shrinkr_fit` object.
#' @param probs Numeric vector of quantiles to compute. Default is 
#'   `c(0.025, 0.5, 0.975)`.
#' @param group_names Optional character vector to label groups.
#' @param digits Number of digits to display. Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list with `mu_tau` and `theta` summary tables.
#'   Prints formatted output.
#'
#' @seealso
#' [shrink()] for fitting models,
#' [summarise_theta()] for theta-only summaries
#' @method summary shrinkr_fit
#' @export
summary.shrinkr_fit <- function(object, 
                                probs = c(0.025, 0.5, 0.975),
                                group_names = NULL,
                                digits = 3,
                                ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg posterior} required. Install with: {.code install.packages('posterior')}")
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(object, "shrinkr_fit"))
  
  cli::cli_rule("Shrinkr Model Summary")
  cli::cli_text("")
  
  # Hyperparameters
  cli::cli_h3("Hyperparameters (mu, tau)")
  mt_summary <- summarise_mu_tau(object, probs = probs)
  print(mt_summary, digits = digits)
  
  cli::cli_text("")
  
  # Theta parameters
  cli::cli_h3("Subgroup effects (theta)")
  theta_summary <- summarise_theta(object, probs = probs, group_names = group_names)
  print(theta_summary, digits = digits)
  
  cli::cli_text("")
  cli::cli_rule()
  
  invisible(list(mu_tau = mt_summary, theta = theta_summary))
}


#' Plot shrinkage fit
#'
#' @description
#' Visualizes the hierarchical shrinkage model fit. Creates either:
#' \itemize{
#'   \item `"shrinkage"` — Shows pre/post-shrunk estimates with arrows
#'   \item `"diagnostics"` — Multi-panel view with hyperparameters and shrinkage factor
#' }
#'
#' The shrinkage plot displays:
#' \itemize{
#'   \item **Pre-shrunk estimates** (hollow circles) — from stage 1 mixture or MLEs
#'   \item **Post-shrunk estimates** (filled circles) — posterior means of theta
#'   \item **Global mean** (dashed line) — posterior mean of mu
#'   \item **Credible intervals** (optional) — for shrunken estimates
#'   \item **Arrows** (optional) — showing direction and magnitude of shrinkage
#' }
#'
#' @param x A `shrinkr_fit` object from [shrink()].
#' @param type Character; type of plot. Options:
#'   \itemize{
#'     \item `"shrinkage"` - Basic shrinkage visualization (default)
#'     \item `"diagnostics"` - Multi-panel with hyperparameters and shrinkage factor
#'   }
#' @param group_names Optional character vector of length G to label groups.
#'   If `NULL`, uses names from `x$data$vars` or defaults to "group1", etc.
#' @param show_arrows Logical; draw arrows from pre-shrunk to post-shrunk estimates?
#'   Default `FALSE`. Only applies when `type = "shrinkage"`.
#' @param show_intervals Logical; show credible intervals for both pre-shrunk and
#'   post-shrunk estimates? Default `TRUE`. Only applies when `type = "shrinkage"`.
#' @param interval_prob Numeric; probability mass for credible intervals.
#'   Default 0.95 for 95% intervals. Only applies when `type = "shrinkage"`.
#' @param point_size Numeric; size of points. Default 3.
#' @param arrow_alpha Numeric; transparency of arrows (0-1). Default 0.6.
#'   Only applies when `show_arrows = TRUE`.
#' @param dodge_width Numeric; horizontal spacing between pre-shrunk and 
#'   post-shrunk estimates in the side-by-side display. Default 0.3.
#'   Larger values increase separation between estimate types.
#' @param title Character; plot title. If `NULL`, uses default title.
#' @param subtitle Character; plot subtitle. If `NULL`, auto-generates from 
#'   global mean and tau.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object (for `type = "shrinkage"`), or a patchwork object/list
#'   (for `type = "diagnostics"`).
#'
#' @seealso 
#' [shrink()] for fitting models,
#' [extract_mu_tau()] for hyperparameter draws
#'
#' @export
#' @method plot shrinkr_fit
#'
#' @examples
#' \dontrun{
#' library(distributional)
#' 
#' # Fit model
#' priors <- list(
#'   mu = dist_normal(0, 5),
#'   tau = dist_truncated(dist_normal(0, 2.5), lower = 0)
#' )
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Basic shrinkage plot with side-by-side estimates
#' plot(fit)
#' plot(fit, type = "shrinkage")
#' 
#' # Full diagnostics
#' plot(fit, type = "diagnostics")
#' 
#' # Customized shrinkage plot
#' plot(
#'   fit,
#'   group_names = c("Control", "Low Dose", "Med Dose", "High Dose"),
#'   show_arrows = TRUE,
#'   interval_prob = 0.9,
#'   dodge_width = 0.4
#' )
#' 
#' # Minimal version
#' plot(fit, show_arrows = FALSE, show_intervals = FALSE)
#' }
plot.shrinkr_fit <- function(x,
                             type = c("shrinkage", "diagnostics"),
                             group_names = NULL,
                             show_arrows = FALSE,
                             show_intervals = TRUE,
                             interval_prob = 0.95,
                             point_size = 3,
                             arrow_alpha = 0.6,
                             dodge_width = 0.3,
                             title = NULL,
                             subtitle = NULL,
                             ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} required. Install with: {.code install.packages('ggplot2')}")
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(x, "shrinkr_fit"))
  type <- match.arg(type)
  
  if (type == "diagnostics") {
    return(.plot_shrink_diagnostics(x, group_names))
  }
  
  # Get dimensions
  G <- x$data$G
  
  # Determine group names
  if (is.null(group_names)) {
    if (!is.null(x$data$vars) && length(x$data$vars) == G) {
      group_names <- x$data$vars
    } else {
      group_names <- paste0("Group ", seq_len(G))
    }
  } else if (length(group_names) != G) {
    cli::cli_abort("{.arg group_names} must have length G = {G}")
  }
  
  # Get pre-shrunk estimates and intervals from stored quantiles
  if (is.null(x$data$quantiles)) {
    cli::cli_abort("Pre-shrunk quantiles not found in fitted object")
  }
  
  quants <- x$data$quantiles
  if (nrow(quants) != G) {
    cli::cli_abort("Quantiles dimension mismatch: expected {G} rows, got {nrow(quants)}")
  }
  
  # Extract pre-shrunk estimates (use median if available, otherwise weighted mean from mixture)
  if ("q50" %in% names(quants)) {
    pre_mean <- quants$q50
  } else {
    # Fallback to computing from mixture components
    if (x$data$K == 1) {
      pre_mean <- as.numeric(x$data$m[[1]])
    } else {
      w <- as.numeric(x$data$w)
      m_mat <- do.call(rbind, x$data$m)
      pre_mean <- colSums(w * m_mat)
    }
  }
  
  # Get pre-shrunk credible intervals based on requested probability
  alpha <- 1 - interval_prob
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2
  
  # Map requested quantiles to available columns
  lower_col <- paste0("q", lower_q * 100)
  upper_col <- paste0("q", upper_q * 100)
  
  # If exact quantile not available, use closest
  available_q <- grep("^q[0-9.]+$", names(quants), value = TRUE)
  if (!lower_col %in% available_q) {
    warning("Requested lower quantile not available, using q2.5")
    lower_col <- "q2.5"
  }
  if (!upper_col %in% available_q) {
    warning("Requested upper quantile not available, using q97.5")
    upper_col <- "q97.5"
  }
  
  pre_lower <- quants[[lower_col]]
  pre_upper <- quants[[upper_col]]
  
  # Get post-shrunk estimates and intervals
  theta_summary <- x$summary[grepl("^theta\\[", x$summary$variable), ]
  theta_indices <- as.integer(gsub("theta\\[(\\d+)\\]", "\\1", theta_summary$variable))
  theta_summary <- theta_summary[theta_indices <= G, ]
  
  if (nrow(theta_summary) != G) {
    cli::cli_abort("Expected {G} theta parameters, found {nrow(theta_summary)}")
  }
  
  post_mean <- theta_summary$mean
  post_lower <- theta_summary$q2.5
  post_upper <- theta_summary$q97.5
  
  # Get global mean
  mu_mean <- x$summary[x$summary$variable == "mu", "mean"]
  
  # Create long-format data for side-by-side plotting
  plot_data <- data.frame(
    group = rep(factor(group_names, levels = group_names), 2),
    group_idx = rep(seq_len(G), 2),
    estimate_type = factor(
      rep(c("Pre-shrunk", "Post-shrunk"), each = G),
      levels = c("Pre-shrunk", "Post-shrunk")
    ),
    mean = c(pre_mean, post_mean),
    lower = c(pre_lower, post_lower),
    upper = c(pre_upper, post_upper)
  )
  
  # Calculate x positions for side-by-side display
  plot_data$x_pos <- plot_data$group_idx + 
    ifelse(plot_data$estimate_type == "Pre-shrunk", -dodge_width/2, dodge_width/2)
  
  # Set title and subtitle
  if (is.null(title)) {
    title <- "Hierarchical Shrinkage"
  }
  
  if (is.null(subtitle)) {
    tau_mean <- x$summary[x$summary$variable == "tau", "mean"]
    subtitle <- sprintf("Global mean: %.2f, tau: %.2f", mu_mean, tau_mean)
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_pos, y = mean, 
                                               fill = estimate_type, 
                                               color = estimate_type)) +
    ggplot2::geom_hline(yintercept = mu_mean, 
                        linetype = "dashed", 
                        color = "gray40",
                        linewidth = 0.8) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  # Add intervals if requested
  if (show_intervals) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.15,
      position = ggplot2::position_identity(),
      alpha = 0.7
    )
  }
  
  # Add points
  p <- p +
    ggplot2::geom_point(
      size = point_size,
      shape = 21,
      color = "black"
    )
  
  # Add arrows if requested
  if (show_arrows) {
    arrow_data <- data.frame(
      group_idx = seq_len(G),
      x = seq_len(G) - dodge_width/2,
      xend = seq_len(G) + dodge_width/2,
      y = pre_mean,
      yend = post_mean
    )
    
    p <- p + ggplot2::geom_segment(
      data = arrow_data,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed"),
      alpha = arrow_alpha,
      color = "gray50",
      linewidth = 0.5,
      inherit.aes = FALSE
    )
  }
  
  # Labels and legend
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq_len(G),
      labels = group_names
    ) +
    ggplot2::scale_fill_manual(
      values = c("Pre-shrunk" = "white", "Post-shrunk" = "steelblue"),
      name = ""
    ) +
    ggplot2::scale_color_manual(
      values = c("Pre-shrunk" = "gray30", "Post-shrunk" = "steelblue"),
      name = ""
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "",
      y = "Estimate"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4)),
      color = "none"
    )
  
  p
}

#' Internal diagnostics plot
#' @keywords internal
#' @noRd
.plot_shrink_diagnostics <- function(x, group_names = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for diagnostics plots")
  }
  
  # Try to load patchwork for multi-panel layout
  has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
  
  # Extract posterior draws
  posterior <- rstan::extract(x$fit, pars = c("mu", "tau", "theta"))
  
  G <- x$data$G
  
  # Determine group names
  if (is.null(group_names)) {
    if (!is.null(x$data$vars) && length(x$data$vars) == G) {
      group_names <- x$data$vars
    } else {
      group_names <- paste0("Group ", seq_len(G))
    }
  }
  
  # 1. Hyperparameters plot (mu and tau combined)
  hyper_df <- data.frame(
    parameter = rep(c("mu", "tau"), each = nrow(posterior$mu)),
    value = c(as.numeric(posterior$mu), as.numeric(posterior$tau))
  )
  
  p_hyper <- ggplot2::ggplot(hyper_df, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = "skyblue", alpha = 0.6) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::labs(
      title = "Posterior: Hyperparameters",
      x = "Value",
      y = "Density"
    ) +
    ggplot2::theme_minimal()
  
  # 2. Theta posteriors (group effects) - only first G
  theta_mat <- posterior$theta[, 1:G, drop = FALSE]
  theta_long <- data.frame(
    group = rep(group_names, each = nrow(theta_mat)),
    theta = as.numeric(theta_mat),
    group_idx = rep(seq_len(G), each = nrow(theta_mat))
  )
  theta_long$group <- factor(theta_long$group, levels = group_names)
  
  p_theta <- ggplot2::ggplot(theta_long, ggplot2::aes(x = theta, y = group)) +
    ggplot2::geom_violin(fill = "lightcoral", alpha = 0.6) +
    ggplot2::geom_vline(xintercept = mean(posterior$mu), 
                        linetype = "dashed", 
                        color = "gray40") +
    ggplot2::labs(
      title = "Posterior: Subgroup effects (theta)",
      x = "Effect size",
      y = ""
    ) +
    ggplot2::theme_minimal()
  
  # Combine plots
  if (has_patchwork) {
    combined <- p_hyper / p_theta
    return(combined)
  } else {
    # Return list if patchwork not available
    message("Install 'patchwork' for combined plot layout: install.packages('patchwork')")
    return(list(
      hyperparameters = p_hyper,
      theta = p_theta
    ))
  }
} 

#' Extract mu and tau parameters
#'
#' @description
#' Extracts posterior draws for the hyperparameters mu (global mean) and 
#' tau (heterogeneity standard deviation) from a fitted shrinkage model.
#'
#' @param x A `shrinkr_fit` object from [shrink()].
#' @param ... Additional arguments (currently unused).
#'
#' @return A `posterior::draws_df` with columns:
#'   \describe{
#'     \item{`.chain`}{Chain index}
#'     \item{`.iteration`}{Iteration within chain}
#'     \item{`.draw`}{Overall draw index}
#'     \item{`mu`}{Global mean parameter}
#'     \item{`tau`}{Heterogeneity parameter}
#'     \item{`tau_squared`}{Variance (tau^2)}
#'   }
#'
#' @seealso
#' [shrink()] for fitting models,
#' [summarise_mu_tau()] for summary statistics,
#' [as_draws_df.shrinkr_fit()] for all parameters
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Extract hyperparameter draws
#' mu_tau <- extract_mu_tau(fit)
#' 
#' # Summarize
#' summarise_mu_tau(fit)
#' 
#' # Visualize
#' library(bayesplot)
#' mcmc_pairs(mu_tau, pars = c("mu", "tau"))
#' }
extract_mu_tau <- function(x, ...) {
  # Extract all draws directly from Stan fit (not via S3 method)
  draws <- posterior::as_draws_df(x$fit, ...)
  
  # Subset to only mu, tau, and tau_squared
  posterior::subset_draws(draws, variable = c("mu", "tau", "tau_squared"))
}


#' Extract theta (group-level effect) parameters
#'
#' @description
#' Extracts posterior draws for the group-level effects (theta parameters) 
#' from a fitted shrinkage model. This is the hierarchically shrunk version
#' of the subgroup effects.
#'
#' @param x A `shrinkr_fit` object from [shrink()].
#' @param ... Additional arguments passed to [as_draws_df.shrinkr_fit()].
#'
#' @return A `posterior::draws_df` with columns:
#'   \describe{
#'     \item{`.chain`}{Chain index}
#'     \item{`.iteration`}{Iteration within chain}
#'     \item{`.draw`}{Overall draw index}
#'     \item{`theta[1]`, `theta[2]`, ...}{Group-level effects}
#'   }
#'
#' @seealso
#' [shrink()] for fitting models,
#' [extract_mu_tau()] for hyperparameters,
#' [summarise_theta()] for summary statistics,
#' [theta_contrasts()] for pairwise comparisons
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Extract theta draws
#' theta_draws <- extract_theta(fit)
#' 
#' # Summarize
#' library(posterior)
#' summarise_draws(theta_draws)
#' 
#' # Visualize
#' library(bayesplot)
#' mcmc_intervals(theta_draws)
#' 
#' # Compare to summaries
#' theta_summary <- summarise_theta(fit)
#' print(theta_summary)
#' }
extract_theta <- function(x, ...) {
  stopifnot(inherits(x, "shrinkr_fit"))
  
  # Get number of groups
  G <- x$data$G
  
  # Construct exact parameter names for theta (not theta_c)
  theta_vars <- paste0("theta[", seq_len(G), "]")
  
  # Extract only these specific variables
  as_draws_df.shrinkr_fit(x, variables = theta_vars, ...)
}

#' Convert shrinkr_fit to draws_df
#'
#' @description
#' Extracts posterior draws in tidy format using the `posterior` package.
#' By default returns user-facing parameters (mu, tau, theta, etc.) and 
#' excludes internal parameterization details. Set `include_internals = TRUE`
#' to access all parameters including theta_c and z.
#'
#' @param x A `shrinkr_fit` object from [shrink()].
#' @param variables Character vector of parameter names to extract. Options include:
#'   \itemize{
#'     \item `"mu"` - Global mean
#'     \item `"tau"` - Heterogeneity SD
#'     \item `"tau_squared"` - Heterogeneity variance
#'     \item `"theta"` or `"theta[i]"` - Subgroup effects
#'   }
#'   If `NULL` (default), returns all user-facing parameters.
#' @param include_internals Logical; if `TRUE`, includes internal Stan parameters
#'   (`theta_c`, `z`) used for parameterization. Default `FALSE`. Only applies
#'   when `variables = NULL`.
#' @param ... Additional arguments passed to `posterior::as_draws_df()`.
#'
#' @return A `posterior::draws_df` with columns for chain, iteration, draw,
#'   and requested parameters.
#'
#' @seealso
#' [shrink()] for fitting models,
#' [extract_mu_tau()] for hyperparameters only
#'
#' @importFrom posterior as_draws_df
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # User-facing parameters only (default)
#' all_draws <- as_draws_df(fit)
#' variables(all_draws)  # mu, tau, theta[1], ..., tau_squared
#' 
#' # Include internal parameters for diagnostics
#' all_draws_internal <- as_draws_df(fit, include_internals = TRUE)
#' variables(all_draws_internal)  # includes theta_c, z
#' 
#' # Just theta parameters
#' theta_draws <- as_draws_df(fit, variables = "theta")
#' 
#' # Specific thetas
#' theta12_draws <- as_draws_df(fit, variables = c("theta[1]", "theta[2]"))
#' 
#' # Work with draws
#' library(posterior)
#' summarise_draws(all_draws)
#' }
as_draws_df.shrinkr_fit <- function(x, 
                                    variables = NULL,
                                    include_internals = FALSE,
                                    ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg posterior} required. Install with: {.code install.packages('posterior')}")
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(x, "shrinkr_fit"))
  
  # Extract from Stan fit
  draws <- posterior::as_draws_df(x$fit, variables = variables, ...)
  
  # Filter out internal parameters unless requested
  if (!include_internals && is.null(variables)) {
    # Remove theta_c and z (internal parameterization variables)
    internal_params <- c("theta_c", "z")
    all_vars <- posterior::variables(draws)
    
    # Find variables to exclude (theta_c[...] and z[...])
    exclude_pattern <- paste0("^(", paste(internal_params, collapse = "|"), ")\\[")
    vars_to_exclude <- grep(exclude_pattern, all_vars, value = TRUE)
    
    if (length(vars_to_exclude) > 0) {
      # Keep only non-internal variables
      vars_to_keep <- setdiff(all_vars, vars_to_exclude)
      draws <- posterior::subset_draws(draws, variable = vars_to_keep)
    }
  }
  
  draws
}


#' Convert shrinkr_fit to data.frame
#'
#' @description
#' Extracts posterior draws as a regular data frame. This is a convenience
#' wrapper around `as_draws_df()` that returns a plain data.frame.
#'
#' @param x A `shrinkr_fit` object from [shrink()].
#' @param row.names NULL or character vector giving row names.
#' @param optional Logical; if `TRUE`, setting row names and converting
#'   column names is optional.
#' @param variables Character vector of parameter names to extract.
#'   If `NULL`, returns all user-facing parameters (excludes internals).
#' @param include_internals Logical; if `TRUE`, includes internal Stan parameters.
#'   Default `FALSE`.
#' @param ... Additional arguments passed to `as_draws_df()`.
#'
#' @return A data.frame with columns for chain, iteration, draw, and
#'   requested parameters.
#'
#' @seealso
#' [as_draws_df.shrinkr_fit()] for posterior package format,
#' [extract_mu_tau()] for hyperparameters only
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Extract as data frame
#' draws_df <- as.data.frame(fit)
#' head(draws_df)
#' 
#' # Just mu and tau
#' mu_tau_df <- as.data.frame(fit, variables = c("mu", "tau"))
#' }
as.data.frame.shrinkr_fit <- function(x, 
                                      row.names = NULL, 
                                      optional = FALSE, 
                                      variables = NULL,
                                      include_internals = FALSE,
                                      ...) {
  as.data.frame(
    as_draws_df.shrinkr_fit(x, variables = variables, include_internals = include_internals, ...),
    row.names = row.names,
    optional = optional
  )
}


#' Summarize theta parameters by group
#'
#' @description
#' Computes posterior summaries for subgroup effects (theta parameters).
#' Returns a data frame with one row per group containing posterior means,
#' standard deviations, quantiles, and convergence diagnostics.
#' 
#' This is a focused alternative to `summary(fit)`, which returns summaries
#' for all parameters including mu and tau.
#'
#' @param fit A `shrinkr_fit` object from [shrink()].
#' @param probs Numeric vector of quantiles to compute. Default is 
#'   `c(0.025, 0.5, 0.975)` for 95% credible intervals.
#' @param group_names Optional character vector of length G to label groups.
#'   If `NULL`, uses names from `fit$data$vars` or defaults to "group1", etc.
#' @param measures Optional character vector or list of summary measures to compute. 
#'   If `NULL`, uses mean, sd, and convergence diagnostics.
#'
#' @return A data frame (tibble if available) with one row per group and columns:
#'   \describe{
#'     \item{`group`}{Group identifier}
#'     \item{`mean`}{Posterior mean}
#'     \item{`sd`}{Posterior standard deviation}
#'     \item{`q2.5`, `q50`, `q97.5`}{Quantiles (or custom quantiles from `probs`)}
#'     \item{`rhat`}{R-hat convergence diagnostic}
#'     \item{`ess_bulk`}{Effective sample size (bulk)}
#'     \item{`ess_tail`}{Effective sample size (tail)}
#'   }
#'
#' @seealso
#' [shrink()] for fitting models,
#' [summarise_mu_tau()] for hyperparameter summaries,
#' [theta_contrasts()] for computing contrasts
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Basic summary
#' theta_summary <- summarise_theta(fit)
#' print(theta_summary)
#' 
#' # Custom quantiles
#' theta_summary <- summarise_theta(fit, probs = c(0.05, 0.5, 0.95))
#' 
#' # With custom group names
#' theta_summary <- summarise_theta(
#'   fit, 
#'   group_names = c("Control", "Treatment A", "Treatment B")
#' )
#' 
#' # Custom measures
#' theta_summary <- summarise_theta(fit, measures = c("mean", "median", "mad"))
#' }
summarise_theta <- function(fit, 
                            probs = c(0.025, 0.5, 0.975),
                            group_names = NULL,
                            measures = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg posterior} required. Install with: {.code install.packages('posterior')}")
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(fit, "shrinkr_fit"))
  
  G <- fit$data$G
  
  # Determine group names
  if (is.null(group_names)) {
    if (!is.null(fit$data$vars) && length(fit$data$vars) == G) {
      group_names <- fit$data$vars
    } else {
      group_names <- paste0("group", seq_len(G))
    }
  } else if (length(group_names) != G) {
    cli::cli_abort("{.arg group_names} must have length G = {G}")
  }
  
  # Extract theta draws using extract_theta (automatically gets only first G)
  draws_df <- extract_theta(fit)
  
  # Create explicit variable names for the G theta parameters
  theta_vars <- paste0("theta[", seq_len(G), "]")
  
  # Subset to only these variables (in case there are extras)
  draws_df <- posterior::subset_draws(draws_df, variable = theta_vars)
  
  # Compute summaries
  if (is.null(measures)) {
    summary_tbl <- posterior::summarise_draws(
      draws_df,
      "mean",
      "sd",
      ~posterior::quantile2(., probs = probs)
    )
    
    # Add convergence diagnostics if available (only if multiple chains)
    if (posterior::nchains(draws_df) > 1) {
      tryCatch({
        conv_stats <- posterior::summarise_draws(
          draws_df,
          "rhat",
          "ess_bulk",
          "ess_tail"
        )
        summary_tbl <- merge(summary_tbl, conv_stats[, c("variable", "rhat", "ess_bulk", "ess_tail")], 
                             by = "variable", all.x = TRUE)
      }, error = function(e) {
        # Convergence stats failed
      })
    }
  } else {
    summary_tbl <- posterior::summarise_draws(
      draws_df,
      measures,
      ~posterior::quantile2(., probs = probs)
    )
  }
  
  # Verify we have exactly G rows
  if (nrow(summary_tbl) != G) {
    cli::cli_warn("Expected {G} theta parameters but got {nrow(summary_tbl)} rows. Truncating to first {G}.")
    summary_tbl <- summary_tbl[seq_len(G), ]
  }
  
  # Add group names
  summary_tbl$group <- group_names
  
  # Reorder columns to put group first, then remove variable
  cols_to_keep <- setdiff(names(summary_tbl), "variable")
  summary_tbl <- summary_tbl[, c("group", setdiff(cols_to_keep, "group"))]
  
  # Convert to tibble if available
  if (requireNamespace("tibble", quietly = TRUE)) {
    summary_tbl <- tibble::as_tibble(summary_tbl)
  }
  
  summary_tbl
}


#' @rdname summarise_theta
#' @export
summarize_theta <- summarise_theta


#' Summarize mu and tau hyperparameters
#'
#' @description
#' Computes posterior summaries for the hierarchical hyperparameters
#' (`mu`, `tau`, and `tau_squared`). Returns a data frame with one row per
#' parameter containing posterior means, standard deviations, quantiles, and
#' convergence diagnostics.
#'
#' This is a focused alternative to `summary(fit)`, which returns summaries
#' for all parameters including theta.
#'
#' @param fit A `shrinkr_fit` object from [shrink()].
#' @param probs Numeric vector of quantiles to compute. Default is
#'   `c(0.025, 0.5, 0.975)` for 95% credible intervals.
#' @param measures Optional character vector or list of summary measures to compute.
#'   If `NULL`, uses mean, sd, and convergence diagnostics.
#'
#' @return A data frame (tibble if available) with one row per parameter and columns:
#'   \describe{
#'     \item{`parameter`}{Parameter name (`mu`, `tau`, or `tau_squared`)}
#'     \item{`mean`}{Posterior mean}
#'     \item{`sd`}{Posterior standard deviation}
#'     \item{`q2.5`, `q50`, `q97.5`}{Quantiles (or custom quantiles from `probs`)}
#'     \item{`rhat`}{R-hat convergence diagnostic}
#'     \item{`ess_bulk`}{Effective sample size (bulk)}
#'     \item{`ess_tail`}{Effective sample size (tail)}
#'   }
#'
#' @seealso
#' [shrink()] for fitting models,
#' [extract_mu_tau()] for raw hyperparameter draws,
#' [summarise_theta()] for group-level summaries
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#'
#' # Basic summary
#' mu_tau_summary <- summarise_mu_tau(fit)
#' print(mu_tau_summary)
#'
#' # Custom quantiles
#' mu_tau_summary <- summarise_mu_tau(fit, probs = c(0.05, 0.5, 0.95))
#'
#' # Custom measures
#' mu_tau_summary <- summarise_mu_tau(fit, measures = c("mean", "median", "mad"))
#' }
summarise_mu_tau <- function(fit,
                             probs = c(0.025, 0.5, 0.975),
                             measures = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg posterior} required. Install with: {.code install.packages('posterior')}")
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(fit, "shrinkr_fit"))
  
  # Extract mu/tau draws using extract_mu_tau
  draws_df <- extract_mu_tau(fit)
  
  mu_tau_vars <- c("mu", "tau", "tau_squared")
  
  # Subset to only these variables (in case there are extras)
  draws_df <- posterior::subset_draws(draws_df, variable = mu_tau_vars)
  
  # Compute summaries
  if (is.null(measures)) {
    summary_tbl <- posterior::summarise_draws(
      draws_df,
      "mean",
      "sd",
      ~posterior::quantile2(., probs = probs)
    )
    
    # Add convergence diagnostics if available (only if multiple chains)
    if (posterior::nchains(draws_df) > 1) {
      tryCatch({
        conv_stats <- posterior::summarise_draws(
          draws_df,
          "rhat",
          "ess_bulk",
          "ess_tail"
        )
        summary_tbl <- merge(summary_tbl, conv_stats[, c("variable", "rhat", "ess_bulk", "ess_tail")],
                             by = "variable", all.x = TRUE)
      }, error = function(e) {
        # Convergence stats failed
      })
    }
  } else {
    summary_tbl <- posterior::summarise_draws(
      draws_df,
      measures,
      ~posterior::quantile2(., probs = probs)
    )
  }
  
  # Rename `variable` column to `parameter` and put it first
  names(summary_tbl)[names(summary_tbl) == "variable"] <- "parameter"
  cols_to_keep <- c("parameter", setdiff(names(summary_tbl), "parameter"))
  summary_tbl <- summary_tbl[, cols_to_keep]
  
  # Preserve merge order (merge can re-sort alphabetically): restore mu, tau, tau_squared order
  summary_tbl <- summary_tbl[match(mu_tau_vars, summary_tbl$parameter), ]
  rownames(summary_tbl) <- NULL
  
  # Convert to tibble if available
  if (requireNamespace("tibble", quietly = TRUE)) {
    summary_tbl <- tibble::as_tibble(summary_tbl)
  }
  
  summary_tbl
}


#' @rdname summarise_mu_tau
#' @export
summarize_mu_tau <- summarise_mu_tau


#' Linear combinations of theta
#'
#' @description
#' Computes posterior draws for linear combinations of subgroup effects.
#' Useful for pairwise contrasts (e.g., treatment vs control), weighted 
#' averages, or any custom linear estimand involving theta parameters.
#'
#' @param fit A `shrinkr_fit` object from [shrink()].
#' @param contrast_matrix A numeric matrix L with `ncol(L) = G` (number of groups)
#'   and `nrow(L) = M` (number of contrasts). Each row defines one linear combination:
#'   \deqn{contrast_i = L_{i1}\theta_1 + L_{i2}\theta_2 + \ldots + L_{iG}\theta_G}
#' @param labels Optional character vector of length M to name the contrasts.
#'   If `NULL`, uses "contrast1", "contrast2", etc.
#'
#' @return A `posterior::draws_df` with columns `.chain`, `.iteration`, `.draw`,
#'   and one column per contrast.
#'
#' @seealso
#' [shrink()] for fitting models,
#' [summarise_theta()] for basic theta summaries
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Pairwise contrast: group 2 vs group 1
#' L <- matrix(c(-1, 1, 0, 0), nrow = 1)
#' contrast <- theta_contrasts(fit, L, labels = "Trt2_vs_Trt1")
#' 
#' # Multiple contrasts
#' L <- rbind(
#'   c(-1, 1, 0, 0),   # Group 2 vs 1
#'   c(-1, 0, 1, 0),   # Group 3 vs 1
#'   c(0, -1, 1, 0)    # Group 3 vs 2
#' )
#' contrasts <- theta_contrasts(
#'   fit, L,
#'   labels = c("2vs1", "3vs1", "3vs2")
#' )
#' 
#' # Summarize contrasts
#' library(posterior)
#' summarise_draws(contrasts)
#' 
#' # Visualize
#' library(bayesplot)
#' mcmc_areas(contrasts, regex_pars = ".*")
#' 
#' # Probability that contrast > 0
#' mean(contrasts$`2vs1` > 0)
#' }
theta_contrasts <- function(fit, contrast_matrix, labels = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg posterior} required. Install with: {.code install.packages('posterior')}")
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(fit, "shrinkr_fit"))
  
  if (!is.matrix(contrast_matrix)) {
    cli::cli_abort("{.arg contrast_matrix} must be a numeric matrix")
  }
  
  G <- fit$data$G
  M <- nrow(contrast_matrix)
  
  if (ncol(contrast_matrix) != G) {
    cli::cli_abort("{.code ncol(contrast_matrix)} must equal number of groups G = {G}")
  }
  
  # Create labels
  if (is.null(labels)) {
    labels <- paste0("contrast", seq_len(M))
  } else if (length(labels) != M) {
    cli::cli_abort("{.arg labels} must have length {.code nrow(contrast_matrix)} = {M}")
  }
  
  # Extract theta draws using extract_theta
  draws_df <- extract_theta(fit)
  
  # Get theta columns in correct order
  theta_vars <- paste0("theta[", seq_len(G), "]")
  
  # Verify all variables exist
  if (!all(theta_vars %in% posterior::variables(draws_df))) {
    cli::cli_abort("Could not find all theta parameters in draws")
  }
  
  # Extract as plain matrix for linear algebra (avoids draws_df metadata warnings)
  theta_matrix <- posterior::as_draws_matrix(draws_df)
  theta_matrix <- theta_matrix[, theta_vars, drop = FALSE]
  
  # Verify dimensions
  if (ncol(theta_matrix) != G) {
    cli::cli_abort("Expected {G} theta columns, found {ncol(theta_matrix)}")
  }
  
  # Compute contrasts: theta_matrix is n_draws x G, contrast_matrix is M x G
  # Result should be n_draws x M
  contrasts <- unclass(theta_matrix) %*% t(contrast_matrix)
  colnames(contrasts) <- labels
  
  # Build result as plain data frame, then convert to draws_df
  meta <- as.data.frame(draws_df)[, c(".chain", ".iteration", ".draw"), drop = FALSE]
  result <- cbind(meta, as.data.frame(contrasts))
  
  posterior::as_draws_df(result)
}