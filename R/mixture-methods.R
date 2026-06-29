# mixture-methods.R

#' Convert mixture fit to data frame
#'
#' @description
#' Converts a fitted mixture model to a tidy long-format data frame.
#' This is essentially an accessor for the `components` element.
#'
#' @param x A `shrinkr_mixture` object from \code{\link{fit_mixture}}.
#' @param row.names Ignored (for S3 consistency).
#' @param optional Ignored (for S3 consistency).
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame (or tibble if available) containing the component
#'   specifications with columns:
#'   \item{component}{Component number (1 to K)}
#'   \item{variable}{Variable name}
#'   \item{weight}{Component weight (mixing proportion)}
#'   \item{mean}{Component mean for this variable}
#'   \item{sd}{Component marginal standard deviation for this variable}
#'
#' @seealso \code{\link{fit_mixture}} for fitting mixture models
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' samples <- list(
#'   group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
#' )
#' mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
#' df <- as.data.frame(mix)
#' head(df)
as.data.frame.shrinkr_mixture <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot(inherits(x, "shrinkr_mixture"))
  
  result <- x$components
  
  # Convert to tibble if available
  if (requireNamespace("tibble", quietly = TRUE)) {
    result <- tibble::as_tibble(result)
  }
  
  result
}


#' Print method for mixture fits
#'
#' @description
#' Displays summary information about a fitted mixture model in a readable format.
#'
#' @param x A `shrinkr_mixture` object from \code{\link{fit_mixture}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @seealso
#' \code{\link{fit_mixture}} for fitting mixture models,
#' \code{\link{summary.shrinkr_mixture}} for detailed summaries
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' samples <- list(
#'   group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
#' )
#' mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
#' print(mix)
print.shrinkr_mixture <- function(x, ...) {
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  stopifnot(inherits(x, "shrinkr_mixture"))
  
  cli::cli_rule("Gaussian Mixture Model Fit")
  cli::cli_text("")
  
  cli::cli_text("Components: {x$K}")
  cli::cli_text("Variables:  {x$n_vars} ({paste(head(x$vars, 3), collapse = ', ')}{if (x$n_vars > 3) ', ...' else ''})")
  cli::cli_text("Samples:    {x$n_samples}")
  
  if (x$diagnostics$n_removed_na > 0) {
    cli::cli_text("  ({x$diagnostics$n_removed_na} rows removed due to NA)")
  }
  
  cli::cli_text("")
  cli::cli_text("Model: {x$model_name}")
  cli::cli_text("BIC:   {round(x$bic, 2)}")
  cli::cli_text("")
  
  cli::cli_h3("Component weights")
  for (k in seq_len(x$K)) {
    cli::cli_text("  k = {k}: {round(x$weights[k], 4)}")
  }
  
  if (length(x$diagnostics$singular_components) > 0) {
    cli::cli_text("")
    cli::cli_alert_warning("Near-singular covariance in component(s): {paste(x$diagnostics$singular_components, collapse = ', ')}")
  }
  
  cli::cli_text("")
  cli::cli_rule()
  cli::cli_alert_info("Use {.fn plot} to visualize marginal fits")
  cli::cli_alert_info("Use {.fn as.data.frame} for component specifications")
  cli::cli_alert_info("Use {.fn summary} for detailed statistics")
  
  invisible(x)
}

#' Summary statistics for mixture fits
#'
#' @description
#' Computes comprehensive summary statistics for fitted mixture components,
#' including component-wise and variable-wise summaries.
#'
#' @param object A `shrinkr_mixture` object from \code{\link{fit_mixture}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with class "summary.shrinkr_mixture" containing:
#'   \item{components}{Data frame with component weights and sizes}
#'   \item{by_variable}{Data frame with per-variable mixture summaries}
#'   \item{model_info}{List with model selection details}
#'   \item{diagnostics}{Fit diagnostics}
#'
#' @seealso
#' \code{\link{fit_mixture}} for fitting mixture models,
#' \code{\link{print.shrinkr_mixture}} for quick overview
#' @method summary shrinkr_mixture
#' @export
#'
#' @examples
#' set.seed(1)
#' samples <- list(
#'   group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
#' )
#' mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
#' summ <- summary(mix)
#' summ$components
summary.shrinkr_mixture <- function(object, ...) {
  stopifnot(inherits(object, "shrinkr_mixture"))
  
  # Component-level summary
  comp_summary <- data.frame(
    component = seq_len(object$K),
    weight = object$weights,
    stringsAsFactors = FALSE
  )
  
  # Variable-level summary (weighted means and SDs across components)
  var_summary <- do.call(rbind, lapply(object$vars, function(v) {
    comp_v <- object$components[object$components$variable == v, ]
    weighted_mean <- sum(comp_v$weight * comp_v$mean)
    weighted_var <- sum(comp_v$weight * (comp_v$sd^2 + comp_v$mean^2)) - weighted_mean^2
    weighted_sd <- sqrt(pmax(0, weighted_var))
    
    data.frame(
      variable = v,
      weighted_mean = weighted_mean,
      weighted_sd = weighted_sd,
      range_mean = diff(range(comp_v$mean)),
      range_sd = diff(range(comp_v$sd)),
      stringsAsFactors = FALSE
    )
  }))
  
  # Model info
  model_info <- list(
    model_name = object$model_name,
    loglik = object$loglik,
    bic = object$bic,
    K = object$K,
    n_parameters = object$mclust_fit$df
  )
  
  # Convert to tibble if available
  if (requireNamespace("tibble", quietly = TRUE)) {
    comp_summary <- tibble::as_tibble(comp_summary)
    var_summary <- tibble::as_tibble(var_summary)
  }
  
  result <- list(
    components = comp_summary,
    by_variable = var_summary,
    model_info = model_info,
    diagnostics = object$diagnostics
  )
  
  class(result) <- "summary.shrinkr_mixture"
  result
}


#' Print summary of mixture fit
#'
#' @description
#' Displays formatted summary statistics for a fitted mixture model.
#'
#' @param x A summary object from \code{\link{summary.shrinkr_mixture}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @seealso \code{\link{summary.shrinkr_mixture}} for creating summaries
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' samples <- list(
#'   group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
#' )
#' mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
#' summ <- summary(mix)
#' print(summ)
print.summary.shrinkr_mixture <- function(x, digits = 3, ...) {
  stopifnot(inherits(x, "summary.shrinkr_mixture"))
  
  cat("== Mixture Model Summary ============================\n\n")
  cat("Model:", x$model_info$model_name, "\n")
  cat("Components:", x$model_info$K, "\n")
  cat("Log-likelihood:", round(x$model_info$loglik, digits), "\n")
  cat("BIC:", round(x$model_info$bic, digits), "\n")
  cat("Parameters:", x$model_info$n_parameters, "\n\n")
  
  cat("Component weights:\n")
  print(x$components, digits = digits, row.names = FALSE)
  
  cat("\nVariable-wise summaries (weighted across components):\n")
  print(x$by_variable, digits = digits, row.names = FALSE)
  
  if (x$diagnostics$n_removed_na > 0) {
    cat("\nNote:", x$diagnostics$n_removed_na, "rows removed due to missing values\n")
  }
  
  if (length(x$diagnostics$singular_components) > 0) {
    cat("\nWarning: Near-singular covariance in component(s):",
        paste(x$diagnostics$singular_components, collapse = ", "), "\n")
  }
  
  cat("\n-----------------------------------------------------\n")
  
  invisible(x)
}


#' Plot fitted marginal densities or QQ plots for mixture models
#'
#' @description
#' Overlays fitted **marginal** mixture densities from a `shrinkr_mixture`
#' (returned by [fit_mixture()]) on top of the **observed samples** for
#' selected variables, OR creates QQ plots comparing empirical vs fitted quantiles.
#' The function uses the *same coercion logic* as [fit_mixture()] (via an internal
#' `.coerce_draws_df()` helper), ensuring that variable names line up even when
#' users pass a list of matrices.
#'
#' **Important:** For multivariate joint fits, this produces **marginal** 
#' overlays (one panel per variable when faceting). Each marginal density is 
#' computed by summing the weighted component densities for that variable.
#'
#' @param x A `shrinkr_mixture` object from [fit_mixture()].
#' @param draws Optional samples to show as histogram/KDE or for QQ plot. Accepts any
#'   input shape supported by [fit_mixture()]. When `NULL`, only fitted curves are drawn
#'   (QQ plot requires draws).
#' @param variables Character vector of variables to plot. Defaults to all
#'   variables in `x$components$variable`. Variable names must match the names
#'   created by the fitter (and by `.coerce_draws_df()`).
#' @param type One of `c("density","qq")`. Default `"density"` shows density overlay;
#'   `"qq"` creates quantile-quantile plots comparing empirical vs fitted quantiles.
#' @param overlay One of `c("hist","kde","both","none")`. Default `"hist"`.
#'   Only applies when `type = "density"`.
#' @param bins Integer number of bins for the histogram (default `50`).
#' @param kde_bw Bandwidth for `stats::density()`; `NULL` uses the default.
#'   Ignored unless `overlay` is `"kde"` or `"both"`.
#' @param show_components Logical; if `TRUE` (default) overlays per-component
#'   curves using component weights, means, and **marginal** SDs from
#'   `x$components`. Only applies when `type = "density"`.
#' @param facet Logical; if `TRUE` (default) facet by variable when plotting
#'   more than one variable.
#' @param n_points Integer; number of x grid points for evaluating densities
#'   (default `501`). For QQ plots, this controls the number of quantiles to compare.
#' @param verbose Logical; print brief matching diagnostics.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' ## Density plots
#' The total marginal density for each variable \eqn{j} is computed as
#' \deqn{f_j(x)=\sum_{k=1}^{K} w_k \,\phi\!\left(x \mid \mu_{jk},\, \sigma_{jk}\right),}
#' using per-component marginal SDs (`sd`) already stored in `x$components`.
#'
#' The plotting range per variable is taken from the sample range if available
#' (with 5% padding), otherwise from `mean +/- 4*sd` across that variable's
#' components--avoiding non-finite `seq()` errors when samples are absent.
#'
#' ## QQ plots
#' When `type = "qq"`, the function creates quantile-quantile plots by:
#' 1. Computing empirical quantiles from the observed data
#' 2. Computing theoretical quantiles from the fitted mixture CDF via numerical inversion
#' 3. Plotting empirical vs theoretical quantiles with a 45-degree reference line
#'
#' Points falling on the reference line indicate good agreement between the fitted
#' mixture and the data. Systematic deviations suggest model misfit.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' set.seed(1)
#' samples <- list(
#'   group1 = matrix(rnorm(100, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(100, 0.5, 0.5), ncol = 1)
#' )
#' mix <- fit_mixture(samples, K_max = 2, verbose = FALSE)
#' plot(mix, draws = samples, type = "density", variables = c("group1", "group2"))
#'
#' @seealso [fit_mixture()] for fitting mixture models
#'
#' @importFrom posterior as_draws_df
#' @importFrom stats dnorm density pnorm quantile uniroot
#' @export
#' @method plot shrinkr_mixture
plot.shrinkr_mixture <- function(x,
                                 draws = NULL,
                                 variables = NULL,
                                 type = c("density", "qq"),
                                 overlay = c("hist","kde","both","none"),
                                 bins = 50,
                                 kde_bw = NULL,
                                 show_components = TRUE,
                                 facet = TRUE,
                                 n_points = 501,
                                 verbose = FALSE,
                                 ...) {
  type <- match.arg(type)
  overlay <- match.arg(overlay)
  
  # QQ plots require draws
  if (type == "qq" && is.null(draws)) {
    stop("QQ plots require `draws` to be provided.")
  }
  
  # Check is now redundant (dispatch ensures correct class) but keep for clarity
  stopifnot(inherits(x, "shrinkr_mixture"))
  
  comps <- x$components
  needed <- c("component","variable","weight","mean","sd")
  if (!all(needed %in% names(comps))) {
    stop("`x$components` is missing required columns: ",
         paste(setdiff(needed, names(comps)), collapse = ", "))
  }
  
  # ----------------- choose variables -----------------
  all_vars <- unique(comps$variable)
  if (is.null(variables)) variables <- all_vars
  variables <- intersect(variables, all_vars)
  if (!length(variables))
    stop("No requested `variables` found in x$components$variable.")
  comps <- comps[comps$variable %in% variables, , drop = FALSE]
  
  # ----------------- coerce draws (mirror fitter) -----
  df_draws <- NULL
  if (!is.null(draws)) {
    df_draws <- .coerce_draws_df(draws)
    if (isTRUE(verbose)) {
      missing <- setdiff(variables, names(df_draws))
      if (length(missing)) message("No sample column(s) for: ", paste(missing, collapse = ", "))
    }
  }
  
  # ----------------- dispatch to appropriate plot type -----
  if (type == "qq") {
    return(.plot_mixture_qq(comps, df_draws, variables, facet, n_points, verbose))
  } else {
    return(.plot_mixture_density(comps, df_draws, variables, overlay, bins, 
                                 kde_bw, show_components, facet, n_points, verbose))
  }
}

# Helper: Compute mixture CDF for a single variable
.mixture_cdf <- function(x, comps_v) {
  sum(vapply(seq_len(nrow(comps_v)), function(i) {
    comps_v$weight[i] * stats::pnorm(x, mean = comps_v$mean[i], sd = comps_v$sd[i])
  }, numeric(1)))
}

# Helper: Invert mixture CDF to get quantile
.mixture_quantile <- function(p, comps_v, tol = 1e-6) {
  # Find x such that CDF(x) = p
  # Start with a reasonable range based on component parameters
  mu_range <- range(comps_v$mean)
  sd_max <- max(comps_v$sd)
  lower <- mu_range[1] - 5 * sd_max
  upper <- mu_range[2] + 5 * sd_max
  
  tryCatch({
    stats::uniroot(
      function(x) .mixture_cdf(x, comps_v) - p,
      lower = lower,
      upper = upper,
      tol = tol
    )$root
  }, error = function(e) {
    # If uniroot fails, return weighted mean as fallback
    sum(comps_v$weight * comps_v$mean)
  })
}

# Helper: Create QQ plot
.plot_mixture_qq <- function(comps, df_draws, variables, facet, n_points, verbose) {
  
  build_qq <- function(vname) {
    comps_v <- comps[comps$variable == vname, , drop = FALSE]
    comps_v <- comps_v[order(comps_v$component), , drop = FALSE]
    
    # Get empirical samples
    if (is.null(df_draws) || !vname %in% names(df_draws)) {
      if (verbose) message("Skipping QQ plot for '", vname, "': no sample data")
      return(NULL)
    }
    
    xv <- suppressWarnings(as.numeric(df_draws[[vname]]))
    xv <- xv[is.finite(xv)]
    if (length(xv) == 0) {
      if (verbose) message("Skipping QQ plot for '", vname, "': no valid samples")
      return(NULL)
    }
    
    # Compute quantiles
    probs <- seq(0.01, 0.99, length.out = min(n_points, length(xv)))
    empirical_q <- stats::quantile(xv, probs = probs, names = FALSE)
    
    # Compute theoretical quantiles from mixture
    theoretical_q <- vapply(probs, function(p) {
      .mixture_quantile(p, comps_v)
    }, numeric(1))
    
    tibble::tibble(
      empirical = empirical_q,
      theoretical = theoretical_q,
      variable = vname
    )
  }
  
  qq_data <- purrr::map_dfr(variables, build_qq)
  
  if (is.null(qq_data) || nrow(qq_data) == 0) {
    stop("No QQ plot data could be generated. Check that samples are available.")
  }
  
  if (isTRUE(verbose)) {
    message("Creating QQ plots for: ", paste(variables, collapse = ", "))
  }
  
  # Create plot
  p <- ggplot2::ggplot(qq_data, ggplot2::aes(x = theoretical, y = empirical)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, 
                         color = "gray40", linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_point(color = "steelblue", size = 1.5, alpha = 0.6) +
    ggplot2::labs(
      x = "Theoretical Quantiles (Mixture)",
      y = "Empirical Quantiles (Data)",
      title = "QQ Plot: Mixture vs Empirical",
      subtitle = "Points should fall on diagonal line for good fit"
    ) +
    ggplot2::theme_minimal()
  
  if (facet && length(variables) > 1L) {
    p <- p + ggplot2::facet_wrap(dplyr::vars(variable), scales = "free")
  } else {
    # Only use coord_equal when not faceting with free scales
    p <- p + ggplot2::coord_equal()
  }
  
  p
}

# Helper: Create density plot (original functionality)
.plot_mixture_density <- function(comps, df_draws, variables, overlay, bins,
                                  kde_bw, show_components, facet, n_points, verbose) {
  
  build_one <- function(vname) {
    comps_v <- comps[comps$variable == vname, , drop = FALSE]
    comps_v <- comps_v[order(comps_v$component), , drop = FALSE]
    
    # Range: from samples if present; else from mixture params
    if (!is.null(df_draws) && vname %in% names(df_draws)) {
      xv <- suppressWarnings(as.numeric(df_draws[[vname]]))
      xv <- xv[is.finite(xv)]
    } else {
      xv <- numeric(0)
    }
    if (length(xv)) {
      r <- range(xv); pad <- diff(r) * 0.05; if (!is.finite(pad) || pad == 0) pad <- 1
      rng <- c(r[1] - pad, r[2] + pad)
    } else {
      mu <- comps_v$mean; sd <- comps_v$sd
      lo <- min(mu - 4 * sd, na.rm = TRUE)
      hi <- max(mu + 4 * sd, na.rm = TRUE)
      if (!is.finite(lo) || !is.finite(hi) || lo >= hi)
        stop("Invalid plotting range for '", vname, "'.")
      rng <- c(lo, hi)
    }
    grid <- seq(rng[1], rng[2], length.out = n_points)
    
    # Total marginal density from per-component marginals
    dens <- rowSums(vapply(seq_len(nrow(comps_v)), function(i) {
      comps_v$weight[i] * stats::dnorm(grid, mean = comps_v$mean[i], sd = comps_v$sd[i])
    }, numeric(length(grid))))
    dens_df <- tibble::tibble(x = grid, density = dens, variable = vname)
    
    comp_df <- NULL
    if (isTRUE(show_components)) {
      comp_df <- purrr::map_dfr(seq_len(nrow(comps_v)), function(i) {
        tibble::tibble(
          x = grid,
          comp_density = comps_v$weight[i] *
            stats::dnorm(grid, mean = comps_v$mean[i], sd = comps_v$sd[i]),
          component = paste0("k=", comps_v$component[i]),
          variable  = vname
        )
      })
    }
    
    # Sample overlays
    draws_hist <- NULL
    draws_kde  <- NULL
    if (length(xv)) {
      if (overlay %in% c("hist","both")) {
        draws_hist <- tibble::tibble(x = xv, variable = vname)
      }
      if (overlay %in% c("kde","both")) {
        kd <- stats::density(xv, bw = kde_bw, n = max(512, n_points),
                             from = rng[1], to = rng[2])
        draws_kde <- tibble::tibble(x = kd$x, density = kd$y, variable = vname)
      }
    }
    
    list(dens = dens_df, comp = comp_df, hist = draws_hist, kde = draws_kde)
  }
  
  built   <- purrr::map(variables, build_one)
  dens_df <- dplyr::bind_rows(purrr::map(built, "dens"))
  comp_df <- dplyr::bind_rows(purrr::compact(purrr::map(built, "comp")))
  hist_df <- dplyr::bind_rows(purrr::compact(purrr::map(built, "hist")))
  kde_df  <- dplyr::bind_rows(purrr::compact(purrr::map(built, "kde")))
  
  if (isTRUE(verbose)) message("Plotting variables: ", paste(variables, collapse = ", "))
  
  # ----------------- ggplot (matching prior_pred styling) ------
  p <- ggplot2::ggplot()
  
  # Histogram in lightcoral (matching theta violin from prior_pred)
  if (nrow(hist_df)) {
    p <- p + ggplot2::geom_histogram(
      data = hist_df,
      ggplot2::aes(x = x, y = ggplot2::after_stat(density)),
      bins = bins, 
      fill = "lightcoral", 
      alpha = 0.25
    )
  }
  
  # KDE in coral (darker shade for visibility)
  if (nrow(kde_df)) {
    p <- p + ggplot2::geom_line(
      data = kde_df,
      ggplot2::aes(x = x, y = density),
      color = "coral3",
      linewidth = 0.7, 
      linetype = "dashed"
    )
  }
  
  # Main marginal density in steelblue (matching skyblue family)
  p <- p + ggplot2::geom_line(
    data = dens_df,
    ggplot2::aes(x = x, y = density),
    color = "steelblue",
    linewidth = 1
  )
  
  # Component densities in skyblue (matching hyperparameter plots)
  if (isTRUE(show_components) && nrow(comp_df)) {
    p <- p + ggplot2::geom_line(
      data = comp_df,
      ggplot2::aes(x = x, y = comp_density, group = component),
      color = "skyblue",
      linewidth = 0.4, 
      alpha = 0.6
    )
  }
  
  if (facet && length(variables) > 1L) {
    p <- p + ggplot2::facet_wrap(dplyr::vars(variable), scales = "free")
  } else if (length(variables) == 1L) {
    p <- p + ggplot2::labs(title = paste("Marginal mixture fit:", variables))
  }
  
  p +
    ggplot2::labs(
      x = "Value", 
      y = "Density",
      title = if (length(variables) > 1L) "Marginal mixture fits" else NULL,
      subtitle = switch(overlay,
                        hist = "Observed samples (histogram) + fitted marginal density",
                        kde  = "Observed samples (KDE, dashed) + fitted marginal density",
                        both = "Observed samples (histogram & KDE) + fitted marginal density",
                        none = "Fitted marginal density")
    ) +
    ggplot2::theme_minimal()
}
