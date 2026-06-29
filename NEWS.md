# shrinkr 0.4.5

## Documentation

* **Improved CRAN compatibility.** Updated package metadata and documentation
  following CRAN review comments, including expanding Markov chain Monte Carlo
  (MCMC) on first use in the `DESCRIPTION` file.

* **Modernized examples.** Reworked documentation examples so lightweight
  examples are executable and self-contained. Replaced unnecessary
  `\dontrun{}` blocks with executable examples or `\donttest{}` where
  appropriate.

* **Protected computationally intensive examples.** Kept examples that require
  fitting Bayesian hierarchical models with Stan inside `\dontrun{}` because
  they are not suitable for routine CRAN example checks.

* **Improved method documentation examples.** Added small example objects for
  selected `shrinkr_fit` and mixture-method documentation so examples no longer
  depend on objects created in other help pages.

## Vignettes

* **Restored graphical parameters.** Updated the `federated_learning` vignette
  to save and restore graphical settings after calls to `par()`, avoiding
  persistent changes to the user's plotting environment.

---

# shrinkr 0.4.4

## New Features

* **New vignette: `vignette("map_prior_with_beastt")`.** Demonstrates how to
  build a robust meta-analytic-predictive (MAP) prior by pairing shrinkr with
  the `beastt` package. shrinkr runs the hierarchical meta-analysis across
  historical control arms and returns the MAP as a `distributional` object,
  which `beastt` then robustifies (`robustify_norm()`) and combines with the
  internal control arm (`calc_post_norm()`) to form the control posterior and
  report an effective sample size. The worked example uses a continuous outcome
  with known SD, and the prior is calibrated with
  `sample_prior_predictive()` / `prior_pairwise_differences()`.

## Documentation

* Added the new MAP-prior vignette to the website articles index and the README
  Documentation section.
* Minor documentation fixes ahead of CRAN submission: package version
  references synchronized to 0.4.4, and small cross-reference and formatting
  tidy-ups in roxygen docs and man pages.

---

# shrinkr 0.4.3
 
## Bug Fixes

* **Fixed inverse-gamma prior translation in `R/utils.R`.** The branch
  handling `dist_inverse_gamma` priors on `tau` was reading from
  `params$shape` and `params$scale`, but `distributional::parameters()`
  returns those fields as `s` and `r`. The check `!is.null(params$shape) &&
  !is.null(params$scale)` was always FALSE, so the branch silently failed
  and downstream code relied on a format-string regex fallback to rescue
  the parsing. End-to-end posteriors were correct (the fallback happened to
  parse `format()` output identically to what the explicit branch should
  have produced), but the code path was fragile. The branch now reads
  `params$s` / `params$r` from `parameters(tau_raw)` (the unwrapped
  distribution, not the possibly-truncated wrapper) and dispatches
  correctly without needing the fallback.
* **Fixed uniform prior translation in `R/utils.R`.** The `dist_uniform`
  branch was reading from `tau_raw$min` and `tau_raw$max`, but
  distributional stores the bounds as `l` and `u`. Same pattern as the IG
  bug: branch silently failed, format-string fallback rescued it. The
  branch now reads `tau_raw$l` / `tau_raw$u` and dispatches correctly.

## Internal Changes

* **Removed unused `stan_code()` generic and methods** from `R/prior_system.R`.
  The generic and its 9 distributional-family methods were vestiges of an
  earlier, abandoned design where Stan code was assembled at runtime from
  prior strings. shrinkr's actual workflow handles all prior dispatch
  internally inside the precompiled Stan model
  (`inst/stan/stage2_shrinkage.stan`) via integer prior codes and parameter
  arrays. The dead code path was never reached at runtime; removing it
  drops ~120 lines of source, 9 `S3method` entries from `NAMESPACE`, and
  10 `.Rd` files from `man/`. `prior_mixture()` and `prior_spike_slab()`
  remain unchanged.
* **Removed the format-string regex fallback in `R/utils.R`.** Roughly 85
  lines across the `mu` and `tau` translation paths used `format(prior)`
  string parsing as a last-resort dispatch when the explicit class-based
  branches missed. With the inverse-gamma fix above, every supported prior
  now dispatches via class-based checks, making the fallback redundant.
  Removing it leaves a cleaner, more transparent failure mode: an
  unsupported prior raises a clear error rather than silently passing
  through a regex parser.

---

# shrinkr 0.4.2

## License

* **shrinkr is now released under GPL (>= 3).** Previous versions did not carry
  a formal open-source license. The `LICENSE` field in `DESCRIPTION` now reads
  `GPL (>= 3)` and the full license text is shipped in `LICENSE`. This
  formalizes the terms under which the package may be used, modified, and
  redistributed.

## Breaking Changes

* **Renamed `compute_theta_contrasts()` to `theta_contrasts()`** for consistency
  with the other `theta`-family functions (`extract_theta()`, `summarise_theta()`).
  The old name has been removed; update calls accordingly.

## New Features

* **`summarise_mu_tau()` / `summarize_mu_tau()`**: Summarize hyperparameters
  (`mu`, `tau`, `tau_squared`) in the same format as `summarise_theta()`. Returns
  posterior mean, sd, quantiles, and convergence diagnostics (`rhat`, `ess_bulk`,
  `ess_tail`) when multiple chains are available. Closes a naming-symmetry gap:
  previously `extract_mu_tau()` existed but had no summary counterpart.

## Internal Changes

* **Migrated Stan array syntax to the `array[N] type name` form** in
  `inst/stan/stage2_shrinkage.stan`. The function signature for
  `builtin_tau_prior_lpdf()` and the `tau_params` / `custom_params` data
  declarations previously used the pre-2.26 syntax (`real[] params`,
  `real tau_params[6]`, `real custom_params[10]`), which emitted stanc3
  deprecation warnings and will become errors in a future Stan release. No
  user-visible change in behavior; R-side data passing is unaffected.

---

# shrinkr 0.4.1

## Bug Fixes

* **Fixed Stan model sampling speed regression**: Restored `std_normal()` prior on
  the inactive parameterization vector (`z` when centered, `theta_c` when
  non-centered). In 0.4.0, the inactive vector had no sampling statement, giving
  it an implicit improper uniform prior. This caused Stan's NUTS
  sampler to waste effort exploring unbounded space, significantly slowing
  sampling. The fix constrains the inactive vector without affecting the model.

## Documentation

* Minor formatting improvements to README.

---

# shrinkr 0.4.0

## Breaking Changes

* **Removed `parse_priors()`**: The legacy formula-based prior interface has been
  removed. Use `distributional` objects (e.g., `dist_normal()`, `dist_student_t()`)
  directly with `shrink()` via the `hierarchical_priors` argument.

* **Removed `shrinkage_factor`** from Stan generated quantities. This quantity
  (`tau^2 / (tau^2 + 1)`) was misleading because it assumed standardized Stage 1
  variance. Use `tau_squared` directly for heterogeneity assessment.

* **Removed `group` column** from `fit_mixture()` output. The column was always
  `1L` and served no purpose. The `variable` column identifies subgroups.

* **`prior_mixture()` and `prior_spike_slab()` now return native `distributional`
  objects** via `distributional::dist_mixture()`. This replaces the previous custom
  `dist_mixture` class. All `distributional` operations (sampling, density, format,
  quantiles) work automatically. Code that accessed `$components` or `$weights`
  fields directly will need updating.

* **`prior_spike_slab()` must be truncated for tau**: Since tau is a scale parameter,
  spike-and-slab priors must be wrapped with `dist_truncated(..., lower = 0)`.
  `sample_prior_predictive()` now errors (rather than warns) if tau draws are negative.
  Example: `tau = dist_truncated(prior_spike_slab(), lower = 0)`.

## New Features

* **Mixture prior support now works end-to-end**: `prior_mixture()` and
  `prior_spike_slab()` now correctly flow through to Stan for both `mu`
  and `tau` priors. Previously, these were exported but not wired into
  the `.coerce_priors_to_stan()` pipeline.

* **Truncated mu priors**: `mu` now supports truncation via
  `dist_truncated(dist_normal(0, 5), lower = 0)`. This enables constraining the
  global mean to a scientifically meaningful range.

* **`prior_pairwise_differences()`**: New function to compute the prior-implied
  distribution of |theta_i - theta_j| from prior predictive samples. Includes
  `print()` and `plot()` S3 methods. The pooled view uses a skyblue density
  (matching the hyperparameter panel) and the by-pair view uses lightcoral violins
  (matching the theta panel). Useful for calibrating priors following the
  recommendation to inspect pairwise difference distributions before fitting.

* **Re-exported `as_draws_df()`** from the `posterior` package so users can call
  `as_draws_df(fit)` without explicitly loading `posterior`.

## Bug Fixes

* **Safe Cholesky factorization**: `.as_chol_factor()` now uses a jitter-and-retry
  strategy (1e-10 through 1e-4) with `Matrix::nearPD()` fallback instead of passing
  a non-lower-triangular eigen decomposition to Stan. The old fallback could silently
  produce invalid `cholesky_factor_cov` input.

* **Stan model no longer exposes unused parameters**: The inactive parameterization
  (`z` when centered, `theta_c` when non-centered) is constrained with a
  `std_normal()` prior so it doesn't appear as a meaningful parameter in output.

* **Fixed all vignette references**: Standardized to underscore naming
  (`getting_started`, `brms_integration`, `tidy_bayesian_workflow`,
  `federated_learning`) throughout R docs, man pages, vignettes, and README.
  Removed references to nonexistent vignettes (`mathematical-foundation`,
  `two_stage_demo`, `tidybayes-integration`).

* **Mixture component extraction handles truncated normals**: The internal
  `.extract_normal_params()` helper correctly unwraps `dist_truncated(dist_normal())`
  to extract the underlying mean and sd for Stan mixture priors. This is needed
  for `dist_truncated(prior_spike_slab(), lower = 0)` to work.

## Testing

* **New `test-prior-system.R`** with 41 tests covering:
  - `prior_mixture()` and `prior_spike_slab()` object creation and sampling
  - `.coerce_priors_to_stan()` for all mu types (Normal, Student-t, truncated, mixture)
  - `.coerce_priors_to_stan()` for all tau types (8 built-in + mixture + error cases)
  - `.as_chol_factor()` safety (positive-definite, near-singular, bad matrices)
  - `sample_prior_predictive()` validation and tau positivity enforcement
  - `prior_pairwise_differences()` structure, edge cases, and probability calibration
  - `fit_mixture()` output structure (no `group` column, quantiles present)
  - End-to-end `shrink()` with spike-and-slab tau, truncated mu, mixture mu
  - Verification that `shrinkage_factor` is absent from output

* **Updated existing tests** for 0.4.0 API changes:
  - `test-shrinkr.R`: Updated mixture prior tests for `distributional::dist_mixture`
  - `test-extraction.R`: Updated `as.data.frame.shrinkr_mixture` to expect no `group` column

## Documentation

* **README substantially improved**:
  - Defined theta on first mention in Quick Start
  - Clarified what `K_max` controls and how mclust selects K via BIC
  - Explained `implied_range` in prior predictive section
  - Added "Why a Mixture of Normals?" section with universal approximation
    justification and K_max sensitivity guidance
  - Noted that standard rstan/bayesplot diagnostic tools work on `fit$fit`
  - Added examples for truncated mu, spike-and-slab (with truncation), and
    pairwise contrasts
  - Highlighted use cases where normal approximation is insufficient (small
    subgroups, rare events)
  - Fixed vignette names in Documentation section

* **Getting started vignette**: Added section on `prior_pairwise_differences()`
  with explanation of clinical calibration for pairwise |theta_i - theta_j|.

* Package version references synchronized to 0.4.0 throughout.

---

# shrinkr 0.3.2

* **Vignettes hosted on website only**: Small package changes to host vignettes on website and not include in build.
  - This will improve CRAN compatibility later


# shrinkr 0.3.1

## New Features

* **Federated Learning Vignette**: Added comprehensive vignette demonstrating shrinkr's 
  application to federated learning scenarios (`vignette("federated_learning")`). 
  The vignette covers:
  - Multi-hospital mortality prediction example with privacy-preserving analysis
  - Two federated learning paths: full posteriors vs. summary statistics
  - Critical guidance on Bayesian CLT assumptions for summary statistics
  - Normality diagnostics with concrete counter-examples when CLT fails
  - Privacy comparison showing data transfer requirements (200 KB vs. 16 bytes)
  - Best practices for data governance, quality control, and transparent reporting
  - Advanced scenarios including heterogeneous models and iterative updates

## Improvements

* **Plot Layout Enhancement**: Improved `plot.shrinkr_fit()` default visualization:
  - Optimized spacing and arrangement of plot elements for better readability
  - Enhanced visual distinction between Stage 1 (independent) and Stage 2 (shrunken) estimates
  - Clearer presentation of global mean and credible intervals
  - More intuitive display of shrinkage effects across groups

## Documentation

* Added detailed warnings about Bayesian CLT assumptions when using summary statistics
  instead of full posterior samples
* Emphasized mixture approximation (Path A) as the safer default for federated learning
* Provided explicit normality checking procedures with QQ plots and diagnostic statistics
* Clarified when each federated learning path is appropriate vs. when it can lead to bias

---

# shrinkr 0.3.0

## First Internal Release

This release focuses on core functionality with a clean, stable API for internal use at GSK.

## Core Features

* **Two-stage Bayesian hierarchical shrinkage framework**: Separate Stage 1 model fitting from Stage 2 hierarchical shrinkage
* **Flexible prior specification** via the `distributional` package for standard prior families:
  - Normal, Student-t, Cauchy, Lognormal
  - Heavy-tailed: Inverse-Gamma, Half-Cauchy, Half-t
  - Bounded: Uniform
  - Truncated distributions
  - Mixture priors: spike-and-slab

## Main Functions

* `fit_mixture()`: Approximate Stage 1 posteriors with Gaussian mixture models
* `shrink()`: Apply hierarchical shrinkage with custom priors
* `prior_spike_slab()`: Create spike-and-slab mixture priors for testing homogeneity
* `prior_mixture()`: Create custom mixture priors with multiple components
* `sample_prior_predictive()`: Generate prior predictive samples for prior elicitation

## Prior Specification & Checking

* Specify priors using `distributional` package syntax
* Test prior implications with `sample_prior_predictive()`
* Visualize prior predictive distributions
* Check implied ranges and spread of effects before fitting

## Extraction & Analysis

* `extract_mu_tau()`: Extract global parameter estimates (mean and heterogeneity)
* `summarise_theta()` / `summarize_theta()`: Summarize group-level estimates with credible intervals
* `compute_theta_contrasts()`: Compute pairwise group comparisons
* `plot()` methods: Visualize shrinkage effects and mixture approximation quality

## Input Methods

* **Full posterior samples** (recommended): Provide mixture approximation that captures posterior uncertainty
* **Point estimates**: Use MLE + variance vector or covariance matrix when full posteriors unavailable

## Integration

* Works with **tidybayes** and **posterior** packages for tidy workflows
* Compatible with outputs from **Stan**, **NIMBLE**, **JAGS**, and **brms**
* Supports both univariate and multivariate parameters

## Documentation

* Comprehensive vignettes:
  - Getting Started: Basic workflow and examples
  - Integration with tidyBayes ecosystem
  - Working with brms models
* Complete function documentation with examples

## Package Quality

* Comprehensive test suite
* All functions documented with roxygen2
* Proper NAMESPACE management
* Precompiled Stan models for efficiency

## Technical Details

* Built on **rstan** for Markov chain Monte Carlo (MCMC) sampling
* Uses **mclust** for mixture model fitting
* Integrates with **posterior** and **distributional** packages
* Supports centered and non-centered parameterizations
* Vectorized Stan code for computational efficiency

---

# shrinkr 0.2.0 (Development)

* Prototype version with extended prior system
* Experimental features for prior selection

---

# shrinkr 0.1.0 (Development)

* Initial development version
* Core two-stage framework established