# prior_system.R
# Mixture-prior helper functions.
#
# Note on prior dispatch architecture: shrinkr's Stan model
# (inst/stan/stage2_shrinkage.stan) handles all prior translation internally
# via integer prior codes (`tau_prior_type`, `mu_prior_type`) and parameter
# arrays. The R side translates user-supplied distributional objects into
# those integer codes and parameters in `R/utils.R`. Earlier versions of
# this file contained an `stan_code()` S3 generic with methods for emitting
# Stan sampling-statement strings; that machinery was never used by the
# runtime workflow and was removed in v0.4.3.

#' @importFrom distributional parameters
NULL

# ============================================================================
# Mixture priors
# ============================================================================

#' Create a mixture prior
#'
#' Creates a mixture distribution using [distributional::dist_mixture()].
#' All standard distributional operations (sampling, density, quantiles,
#' formatting) work automatically.
#'
#' @param ... Component distributions (from distributional package)
#' @param weights Mixture weights (normalized automatically if not summing to 1)
#' @return A `distributional` mixture distribution object
#' @export
#' @examples
#' mix <- prior_mixture(
#'   distributional::dist_normal(0, 0.1),
#'   distributional::dist_normal(0, 1),
#'   weights = c(0.7, 0.3)
#' )
prior_mixture <- function(..., weights = NULL) {
  components <- list(...)
  n_comp <- length(components)
  
  if (is.null(weights)) {
    weights <- rep(1 / n_comp, n_comp)
  } else {
    weights <- weights / sum(weights)
  }
  
  # Use distributional's native dist_mixture
  # It expects individual distribution objects as ..., plus a weights argument
  do.call(distributional::dist_mixture, c(components, list(weights = weights)))
}

#' Spike-and-slab prior for testing homogeneity
#'
#' Creates a mixture of two Normal distributions. Since tau is a scale parameter,
#' this must be wrapped in [distributional::dist_truncated()] with `lower = 0`
#' before passing to [shrink()]:
#'
#' ```
#' tau_prior <- dist_truncated(prior_spike_slab(), lower = 0)
#' ```
#'
#' @param spike_location Location of the spike (default 0)
#' @param spike_scale Scale of the spike component (default 0.01)
#' @param slab_scale Scale of the slab component (default 1)
#' @param spike_prob Probability of the spike component (default 0.5)
#' @return A spike-and-slab mixture distribution
#' @export
#' @examples
#' tau_prior <- distributional::dist_truncated(
#'   prior_spike_slab(spike_prob = 0.5, spike_scale = 0.01, slab_scale = 1),
#'   lower = 0
#' )
prior_spike_slab <- function(spike_location = 0, spike_scale = 0.01,
                             slab_scale = 1, spike_prob = 0.5) {
  prior_mixture(
    distributional::dist_normal(spike_location, spike_scale),
    distributional::dist_normal(spike_location, slab_scale),
    weights = c(spike_prob, 1 - spike_prob)
  )
}
