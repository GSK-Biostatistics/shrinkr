#' Bayesian Hierarchical Shrinkage Model
#'
#' @description
#' Applies hierarchical shrinkage to group-specific estimates using a two-stage
#' Bayesian approach. Takes either a Gaussian mixture approximation of Stage 1 
#' posteriors or point estimates with variance, and applies a Normal hierarchical
#' model with flexible hyperpriors.
#'
#' @param mixture A `shrinkr_mixture` object from [fit_mixture()]. Contains the 
#'   Gaussian mixture approximation of Stage 1 posteriors. Either `mixture` or 
#'   both `mle` and `var_matrix` must be provided.
#' @param mle Numeric vector of group point estimates. Used when `mixture` is NULL.
#' @param var_matrix Numeric vector of variances (length `G`) or covariance 
#'   matrix (`G × G`). Required when `mle` is provided.
#' @param hierarchical_priors Named list with `mu` and `tau` priors as 
#'   `distributional` objects. Defaults to weakly informative priors:
#'   - `mu`: Global mean, `dist_normal(0, 5)`
#'   - `tau`: Between-group SD, `dist_truncated(dist_normal(0, 2.5), lower = 0)`
#'   
#'   Supported distributions for `mu`: Normal, Student-t, mixture priors,
#'   and truncated versions of these (e.g.,
#'   `dist_truncated(dist_normal(0, 5), lower = 0)`).
#'   
#'   Supported distributions for `tau`: Normal (truncated), Student-t (truncated),
#'   Cauchy (truncated), Lognormal, Gamma, Inverse-Gamma, Exponential, Uniform,
#'   and mixture priors (including spike-and-slab via [prior_spike_slab()]).
#' @param centered Logical; use centered (`TRUE`) or non-centered (`FALSE`, default)
#'   parameterization. Non-centered is more efficient when heterogeneity is small.
#' @param verbose Logical; print progress messages (default `TRUE`).
#' @param ... Additional arguments passed to [rstan::sampling()]:
#'   - `chains`: Number of chains (default 4)
#'   - `iter`: Iterations per chain (default 2000)
#'   - `warmup`: Warmup iterations (default iter/2)
#'   - `cores`: Cores for parallel sampling
#'   - `seed`: Random seed
#'   - `control`: List of sampler controls (e.g., `list(adapt_delta = 0.95)`)
#'
#' @return A `shrinkr_fit` object (list) containing:
#' \item{fit}{Stan model object}
#' \item{data}{Data list used for fitting}
#' \item{summary}{Parameter summaries (mean, sd, quantiles, Rhat, ESS)}
#' \item{diagnostics}{Sampler diagnostics (divergences, treedepth)}
#' \item{priors}{Prior specifications used}
#'
#' @details
#' ## Model Specification
#' 
#' **Hierarchical model (Stage 2):**
#' \deqn{\theta_g \mid \mu, \tau \sim \text{Normal}(\mu, \tau^2), \quad g = 1, \ldots, G}
#' \deqn{\mu \sim \pi(\mu)}
#' \deqn{\tau \sim \pi(\tau)}
#' 
#' **Stage 1 likelihood (Gaussian mixture approximation):**
#' \deqn{\theta_g \mid D_g \sim q_g(\theta_g) \approx \sum_{k=1}^K w_k \, \text{MVN}(\mu_k, \Sigma_k)}
#' 
#' **Full posterior:**
#' \deqn{\pi(\theta, \mu, \tau \mid D) \propto 
#'   \left[\prod_{g=1}^G q_g(\theta_g)\right] 
#'   \left[\prod_{g=1}^G \text{Normal}(\theta_g \mid \mu, \tau^2)\right] 
#'   \pi(\mu) \pi(\tau)}
#' 
#' where \eqn{q_g(\theta_g)} approximates the Stage 1 posterior for group \eqn{g}.
#' 
#' ## What's Fixed vs. Flexible
#' 
#' **Fixed:**
#' - Hierarchical distribution: \eqn{\theta_g \mid \mu, \tau \sim \text{Normal}(\mu, \tau^2)}
#' 
#' **Flexible:**
#' - Hyperpriors \eqn{\pi(\mu)} and \eqn{\pi(\tau)}: Normal, Student-t, Cauchy, 
#'   Lognormal, Gamma, Inverse-Gamma, Exponential, Uniform, mixtures (including
#'   spike-and-slab), and truncated versions
#' - Stage 1 posteriors: Can be non-Normal (handled by mixture approximation)
#' 
#' ## Critical Requirements
#' 
#' 1. **Stage 1 must use flat/uninformative priors** on \eqn{\theta_g}
#'    - Ensures two-stage = one-stage hierarchical model
#'    - Stan: Don't specify prior (defaults to flat)
#'    - JAGS/NIMBLE: Use very wide priors
#' 
#' 2. **Verify mixture quality:** `plot(mixture, draws = samples)`
#'    - Check density overlays and QQ plots
#'    - Poor approximation → biased shrinkage
#' 
#' 3. **Check prior implications:** `sample_prior_predictive(hierarchical_priors)`
#'    - Understand what priors imply before fitting
#'    - Avoid prior-data conflicts
#' 
#' 4. **Minimum 2 groups required** for heterogeneity estimation
#' 
#' ## Common Prior Choices for τ
#' 
#' \itemize{
#'   \item Half-Normal: `dist_truncated(dist_normal(0, s), lower = 0)` - Weakly informative
#'   \item Half-t: `dist_truncated(dist_student_t(df, 0, s), lower = 0)` - Heavier tails
#'   \item Half-Cauchy: `dist_truncated(dist_cauchy(0, s), lower = 0)` - Very diffuse
#'   \item Uniform: `dist_uniform(0, U)` - Bounded heterogeneity
#'   \item Inverse-Gamma: `dist_inverse_gamma(a, b)` - Traditional choice
#' }
#' 
#' See `vignette("getting_started")` for complete workflow,
#' `vignette("brms_integration")` for real examples, and package README for
#' mathematical justification.
#'
#' @examples
#' \dontrun{
#' library(distributional)
#' 
#' # Simulate Stage 1 posteriors
#' samples <- list(
#'   group1 = matrix(rnorm(2000, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(2000, 0.5, 0.5), ncol = 1),
#'   group3 = matrix(rnorm(2000, 1.0, 0.5), ncol = 1)
#' )
#' 
#' # Fit mixture approximation
#' mix <- fit_mixture(samples, K_max = 3)
#' plot(mix, draws = samples)  # Check quality
#' 
#' # Apply shrinkage with default priors
#' fit <- shrink(mixture = mix)
#' print(fit)
#' plot(fit)
#' 
#' # Custom priors (half-t for heavier tails)
#' priors <- list(
#'   mu = dist_normal(0, 10),
#'   tau = dist_truncated(dist_student_t(3, 0, 2.5), lower = 0)
#' )
#' fit2 <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Using point estimates only
#' fit3 <- shrink(
#'   mle = c(0.5, 1.2, -0.3),
#'   var_matrix = c(0.1, 0.15, 0.12),
#'   hierarchical_priors = priors
#' )
#' 
#' # Extract results
#' mu_tau <- extract_mu_tau(fit)
#' mu_tau_summary <- summarize_mu_tau(fit)
#' theta_summary <- summarize_theta(fit)
#' }
#'
#' @seealso 
#' **Workflow functions:**
#' [fit_mixture()], [sample_prior_predictive()]
#' 
#' **Extract results:**
#' [extract_mu_tau()], [extract_theta()], [summarize_mu_tau()], [summarize_theta()], [theta_contrasts()]
#' 
#' **Visualization:**
#' [plot.shrinkr_fit()], [plot.shrinkr_mixture()]
#' 
#' **Vignettes:**
#' - `vignette("getting_started")` - Complete workflow with Stan example
#' - `vignette("brms_integration")` - Survival analysis example
#'
#' @export
shrink <- function(mixture = NULL,
                   mle = NULL,
                   var_matrix = NULL,
                   hierarchical_priors = list(
                     mu = distributional::dist_normal(0, 5),
                     tau = distributional::dist_truncated(
                       distributional::dist_normal(0, 2.5), 
                       lower = 0
                     )
                   ),
                   centered = FALSE,
                   verbose = TRUE,
                   ...) {
  
  # Check dependencies
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' required. Install with: install.packages('cli')")
  }
  
  if (!requireNamespace("rstan", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg rstan} is required. Install with: {.code install.packages('rstan')}")
  }
  
  if (!requireNamespace("distributional", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg distributional} is required for prior specification")
  }
  
  # Validate hierarchical_priors structure
  if (!is.list(hierarchical_priors) || !all(c("mu", "tau") %in% names(hierarchical_priors))) {
    cli::cli_abort(
      c("x" = "{.arg hierarchical_priors} must be a list with elements {.field mu} and {.field tau}",
        "i" = "Example: {.code list(mu = dist_normal(0, 5), tau = dist_truncated(dist_normal(0, 2.5), lower = 0))}")
    )
  }
  
  # Validate prior types
  if (!inherits(hierarchical_priors$mu, "distribution") || 
      !inherits(hierarchical_priors$tau, "distribution")) {
    cli::cli_abort(
      c("x" = "{.field hierarchical_priors$mu} and {.field hierarchical_priors$tau} must be distributional objects",
        "i" = "Use functions from the {.pkg distributional} package like {.fn dist_normal}, {.fn dist_student_t}, etc.")
    )
  }
  
  # Convert priors to Stan format
  pri <- .coerce_priors_to_stan(
    prior_mu  = hierarchical_priors$mu,
    prior_tau = hierarchical_priors$tau
  )
  
  # Prepare data based on input type
  if (!is.null(mixture)) {
    if (!inherits(mixture, "shrinkr_mixture")) {
      cli::cli_abort("{.arg mixture} must be a {.cls shrinkr_mixture} object from {.fn fit_mixture}")
    }
    data_list <- .prep_bhm_data_from_mixture(mixture = mixture, pri = pri, centered = centered)
    G <- data_list$G
    vars <- if (!is.null(mixture$vars)) mixture$vars else NULL
    
    # Get quantiles from mixture object
    plot_quantiles <- if (!is.null(mixture$quantiles)) {
      mixture$quantiles
    } else {
      # Fallback: compute from mixture components (shouldn't happen if fit_mixture is updated)
      NULL
    }
    
  } else {
    if (is.null(mle) || is.null(var_matrix)) {
      cli::cli_abort(
        c("x" = "Must provide either {.arg mixture} OR both {.arg mle} and {.arg var_matrix}",
          "i" = "Use {.fn fit_mixture} to create a mixture, or provide point estimates with variance")
      )
    }
    data_list <- .prep_bhm_data_from_mle(mle = mle, var_matrix = var_matrix, pri = pri, centered = centered)
    G <- data_list$G
    
    # Preserve subgroup names when mle is a named vector
    vars <- names(mle)
    if (is.null(vars) || length(vars) != G || any(is.na(vars)) || any(!nzchar(vars))) {
      vars <- paste0("group", seq_len(G))
    }
    data_list$vars <- vars
    
    # Compute quantiles from MLE + variance using normal approximation
    se <- if (is.matrix(var_matrix)) sqrt(diag(var_matrix)) else sqrt(var_matrix)
    
    q_probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    quantiles_list <- lapply(q_probs, function(p) {
      qnorm(p, mean = mle, sd = se)
    })
    
    plot_quantiles <- as.data.frame(do.call(cbind, quantiles_list))
    colnames(plot_quantiles) <- paste0("q", q_probs * 100)
    plot_quantiles$variable <- vars
    plot_quantiles <- plot_quantiles[, c("variable", paste0("q", q_probs * 100))]
  }
  
  # Check for minimum groups
  if (G < 2L) {
    cli::cli_abort(
      c("x" = "Shrinkage requires at least 2 groups",
        "i" = "Found {G} group{?s}. Hierarchical models need multiple groups to estimate between-group variance.")
    )
  }
  
  # Load or compile Stan model
  if (verbose) cli::cli_alert_info("Preparing Stan model...")
  
  # Check if precompiled model exists
  if (exists("stanmodels", where = asNamespace("shrinkr"))) {
    mod <- stanmodels$stage2_shrinkage
  } else {
    # Fall back to compiling from source
    stan_file <- system.file("stan/stage2_shrinkage.stan", package = "shrinkr")
    if (!nzchar(stan_file) || !file.exists(stan_file)) {
      stan_file <- "inst/stan/stage2_shrinkage.stan"  # Development fallback
    }
    if (verbose) cli::cli_alert_info("Compiling Stan model (this may take a minute)...")
    mod <- rstan::stan_model(stan_file)
  }
  
  # Fit the model
  if (verbose) cli::cli_alert_info("Running MCMC sampler...")
  
  fit <- rstan::sampling(
    object  = mod,
    data    = data_list,
    refresh = if (verbose) getOption("shrinkr.refresh", 100L) else 0,
    ...
  )
  
  # Extract parameter summaries
  pars <- c("mu", "tau", "theta", "tau_squared")
  s <- tryCatch(
    rstan::summary(fit, pars = pars)$summary,
    error = function(e) rstan::summary(fit)$summary
  )
  
  # Create tidy summary dataframe
  summ <- data.frame(
    variable = rownames(s),
    mean     = s[, "mean"],
    sd       = s[, "sd"],
    q2.5     = if ("2.5%"  %in% colnames(s)) s[, "2.5%"]  else NA_real_,
    q25      = if ("25%"   %in% colnames(s)) s[, "25%"]   else NA_real_,
    q50      = if ("50%"   %in% colnames(s)) s[, "50%"]   else NA_real_,
    q75      = if ("75%"   %in% colnames(s)) s[, "75%"]   else NA_real_,
    q97.5    = if ("97.5%" %in% colnames(s)) s[, "97.5%"] else NA_real_,
    n_eff    = if ("n_eff" %in% colnames(s)) s[, "n_eff"] else NA_real_,
    Rhat     = if ("Rhat"  %in% colnames(s)) s[, "Rhat"]  else NA_real_,
    row.names   = NULL,
    check.names = FALSE
  )
  
  # Extract diagnostics
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  sp_all <- tryCatch(do.call(rbind, sp), error = function(e) NULL)
  
  diagnostics <- list(
    n_divergent   = if (!is.null(sp_all) && "divergent__" %in% colnames(sp_all))
      sum(sp_all[, "divergent__"]) else NA_integer_,
    max_treedepth = if (!is.null(sp_all) && "treedepth__" %in% colnames(sp_all))
      max(sp_all[, "treedepth__"]) else NA_integer_,
    n_leapfrog    = if (!is.null(sp_all) && "n_leapfrog__" %in% colnames(sp_all))
      max(sp_all[, "n_leapfrog__"]) else NA_integer_
  )
  
  # Add diagnostic warnings if needed
  if (!is.na(diagnostics$n_divergent) && diagnostics$n_divergent > 0) {
    cli::cli_alert_warning("Found {diagnostics$n_divergent} divergent transition{?s}. Consider increasing adapt_delta.")
  }
  
  # Construct and return shrinkr_fit object
  result <- list(
    fit         = fit,
    data        = c(data_list, list(vars = vars, quantiles = plot_quantiles)),
    summary     = summ,
    diagnostics = diagnostics,
    priors      = hierarchical_priors
  ) 
  
  class(result) <- c("shrinkr_fit", "list")
  
  if (verbose) {
    cli::cli_alert_success("Shrinkage complete!")
  }
  
  result
}
