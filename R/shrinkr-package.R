#' shrinkr: Modular Bayesian Hierarchical Shrinkage Models
#'
#' @description
#' The `shrinkr` package provides a flexible framework for two-stage Bayesian 
#' hierarchical modeling. It enables post-hoc shrinkage of subgroup-specific 
#' posterior estimates from any Bayesian model, with support for diverse prior 
#' specifications and diagnostic tools.
#'
#' @section Key Features:
#' 
#' **Two-Stage Workflow:**
#' \itemize{
#'   \item Stage 1: Fit any Bayesian model without shrinkage
#'   \item Stage 2: Apply hierarchical shrinkage with flexible priors
#' }
#' 
#' **Flexible Priors:**
#' \itemize{
#'   \item Standard families (Normal, Student-t, Cauchy, Lognormal)
#'   \item Heavy-tailed (Inverse-Gamma, Half-Cauchy, Half-t)
#'   \item Bounded (Uniform)
#'   \item Mixture priors (spike-and-slab)
#'   \item Truncated distributions
#' }
#' 
#' **Input Methods:**
#' \itemize{
#'   \item Full posterior samples (via mixture approximation)
#'   \item Point estimates + variance/covariance
#' }
#'
#' @section Main Functions:
#' 
#' **Core Workflow:**
#' \itemize{
#'   \item \code{\link{fit_mixture}()}: Approximate Stage 1 posteriors with Gaussian mixture
#'   \item \code{\link{shrink}()}: Main user interface for hierarchical shrinkage
#' }
#' 
#' **Prior Specification:**
#' \itemize{
#'   \item \code{\link{prior_spike_slab}()}: Create spike-and-slab mixture prior
#'   \item \code{\link{prior_mixture}()}: Create custom mixture prior
#'   \item \code{\link{sample_prior_predictive}()}: Generate prior predictive samples for checking
#' }
#' 
#' **Extraction & Visualization:**
#' \itemize{
#'   \item \code{\link{extract_mu_tau}()}: Extract hyperparameter draws
#'   \item \code{\link{extract_theta}()}: Extract group-level draws
#'   \item \code{\link{summarise_mu_tau}()}: Summarize hyperparameters
#'   \item \code{\link{summarise_theta}()}: Summarize group-level estimates
#'   \item \code{\link{theta_contrasts}()}: Compute linear combinations of theta
#'   \item \code{plot()}: Visualize shrinkage effect and mixture approximation quality
#' }
#'
#' @section Getting Started:
#' 
#' See \code{vignette("getting_started", package = "shrinkr")} for a basic workflow,
#' or \code{vignette("brms_integration", package = "shrinkr")} for a survival analysis example.
#'
#' @section Use Cases:
#' 
#' \itemize{
#'   \item **Meta-analysis:** Shrink study-specific effects
#'   \item **Clinical trials:** Borrow information across subgroups or historical controls
#'   \item **Genomics:** Regularize gene-specific effects
#'   \item **Simulation studies:** Compare shrinkage methods systematically
#' }
#'
#' @section Package Options:
#' 
#' \itemize{
#'   \item \code{shrinkr.refresh}: Controls Stan sampling progress output (default: 100)
#' }
#'
#' @docType package
#' @name shrinkr-package
#' @aliases shrinkr
#' 
#' @references
#' Maronge, J. M. (2026). shrinkr: Modular Bayesian Hierarchical Shrinkage Models. 
#' R package version 0.4.3.
#'
#' @seealso
#' \itemize{
#'   \item Stan: \url{https://mc-stan.org/}
#'   \item distributional package: \url{https://pkg.mitchelloharawild.com/distributional/}
#' }
#'
## usethis namespace: start
#' @import Rcpp
#' @import methods  
#' @import rstan
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib shrinkr, .registration = TRUE
## usethis namespace: end
#' @examples
#' \dontrun{
#' library(shrinkr)
#' library(distributional)
#' 
#' # Generate Stage 1 posterior samples
#' set.seed(123)
#' samples <- list(
#'   group1 = matrix(rnorm(2000, 0.0, 0.5), ncol = 1),
#'   group2 = matrix(rnorm(2000, 0.5, 0.5), ncol = 1),
#'   group3 = matrix(rnorm(2000, 1.0, 0.5), ncol = 1)
#' )
#' 
#' # Fit mixture approximation
#' mix <- fit_mixture(samples, K_max = 3)
#' 
#' # Apply hierarchical shrinkage
#' priors <- list(
#'   mu = dist_normal(0, 5),
#'   tau = dist_truncated(dist_student_t(3, 0, 1), lower = 0)
#' )
#' 
#' fit <- shrink(mixture = mix, hierarchical_priors = priors)
#' 
#' # Examine results
#' summary(fit)
#' plot(fit)
#' extract_mu_tau(fit)
#' }
"_PACKAGE"