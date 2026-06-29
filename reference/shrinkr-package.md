# shrinkr: Modular Bayesian Hierarchical Shrinkage Models

The `shrinkr` package provides a flexible framework for two-stage
Bayesian hierarchical modeling. It enables post-hoc shrinkage of
subgroup-specific posterior estimates from any Bayesian model, with
support for diverse prior specifications and diagnostic tools.

## Key Features

**Two-Stage Workflow:**

- Stage 1: Fit any Bayesian model without shrinkage

- Stage 2: Apply hierarchical shrinkage with flexible priors

**Flexible Priors:**

- Standard families (Normal, Student-t, Cauchy, Lognormal)

- Heavy-tailed (Inverse-Gamma, Half-Cauchy, Half-t)

- Bounded (Uniform)

- Mixture priors (spike-and-slab)

- Truncated distributions

**Input Methods:**

- Full posterior samples (via mixture approximation)

- Point estimates + variance/covariance

## Main Functions

**Core Workflow:**

- [`fit_mixture()`](fit_mixture.md): Approximate Stage 1 posteriors with
  Gaussian mixture

- [`shrink()`](shrink.md): Main user interface for hierarchical
  shrinkage

**Prior Specification:**

- [`prior_spike_slab()`](prior_spike_slab.md): Create spike-and-slab
  mixture prior

- [`prior_mixture()`](prior_mixture.md): Create custom mixture prior

- [`sample_prior_predictive()`](sample_prior_predictive.md): Generate
  prior predictive samples for checking

**Extraction & Visualization:**

- [`extract_mu_tau()`](extract_mu_tau.md): Extract hyperparameter draws

- [`extract_theta()`](extract_theta.md): Extract group-level draws

- [`summarise_mu_tau()`](summarise_mu_tau.md): Summarize hyperparameters

- [`summarise_theta()`](summarise_theta.md): Summarize group-level
  estimates

- [`theta_contrasts()`](theta_contrasts.md): Compute linear combinations
  of theta

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): Visualize
  shrinkage effect and mixture approximation quality

## Getting Started

See
[`vignette("getting_started", package = "shrinkr")`](../articles/getting_started.md)
for a basic workflow, or
[`vignette("brms_integration", package = "shrinkr")`](../articles/brms_integration.md)
for a survival analysis example.

## Use Cases

- **Meta-analysis:** Shrink study-specific effects

- **Clinical trials:** Borrow information across subgroups or historical
  controls

- **Genomics:** Regularize gene-specific effects

- **Simulation studies:** Compare shrinkage methods systematically

## Package Options

- `shrinkr.refresh`: Controls Stan sampling progress output (default:
  100)

## References

Maronge, J. M. (2026). shrinkr: Modular Bayesian Hierarchical Shrinkage
Models. R package version 0.4.3.

## See also

- Stan: <https://mc-stan.org/>

- distributional package:
  <https://pkg.mitchelloharawild.com/distributional/>

## Author

**Maintainer**: Jacob M. Maronge <jacob.m.maronge@gsk.com>
([ORCID](https://orcid.org/0000-0003-3606-0841))

Authors:

- Jacob M. Maronge <jacob.m.maronge@gsk.com>
  ([ORCID](https://orcid.org/0000-0003-3606-0841))

Other contributors:

- GlaxoSmithKline Research & Development Limited \[copyright holder,
  funder\]

- Trustees of Columbia University (R/stanmodels.R, configure,
  configure.win) \[copyright holder\]

## Examples

``` r
if (FALSE) { # \dontrun{
# This example fits a Stan model, so it is not run during package checks.
library(shrinkr)
priors <- list(
  mu = distributional::dist_normal(0, 5),
  tau = distributional::dist_truncated(distributional::dist_student_t(3, 0, 1), lower = 0)
)
fit <- shrink(
  mle = c(0.0, 0.5, 1.0),
  var_matrix = c(0.25, 0.25, 0.25),
  hierarchical_priors = priors,
  iter = 1000, chains = 2, seed = 1
)
summary(fit)
} # }
```
