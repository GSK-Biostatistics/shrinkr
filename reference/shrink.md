# Bayesian Hierarchical Shrinkage Model

Applies hierarchical shrinkage to group-specific estimates using a
two-stage Bayesian approach. Takes either a Gaussian mixture
approximation of Stage 1 posteriors or point estimates with variance,
and applies a Normal hierarchical model with flexible hyperpriors.

## Usage

``` r
shrink(
  mixture = NULL,
  mle = NULL,
  var_matrix = NULL,
  hierarchical_priors = list(mu = distributional::dist_normal(0, 5), tau =
    distributional::dist_truncated(distributional::dist_normal(0, 2.5), lower = 0)),
  centered = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- mixture:

  A `shrinkr_mixture` object from [`fit_mixture()`](fit_mixture.md).
  Contains the Gaussian mixture approximation of Stage 1 posteriors.
  Either `mixture` or both `mle` and `var_matrix` must be provided.

- mle:

  Numeric vector of group point estimates. Used when `mixture` is NULL.

- var_matrix:

  Numeric vector of variances (length `G`) or covariance matrix
  (`G × G`). Required when `mle` is provided.

- hierarchical_priors:

  Named list with `mu` and `tau` priors as `distributional` objects.
  Defaults to weakly informative priors:

  - `mu`: Global mean, `dist_normal(0, 5)`

  - `tau`: Between-group SD,
    `dist_truncated(dist_normal(0, 2.5), lower = 0)`

  Supported distributions for `mu`: Normal, Student-t, mixture priors,
  and truncated versions of these (e.g.,
  `dist_truncated(dist_normal(0, 5), lower = 0)`).

  Supported distributions for `tau`: Normal (truncated), Student-t
  (truncated), Cauchy (truncated), Lognormal, Gamma, Inverse-Gamma,
  Exponential, Uniform, and mixture priors (including spike-and-slab via
  [`prior_spike_slab()`](prior_spike_slab.md)).

- centered:

  Logical; use centered (`TRUE`) or non-centered (`FALSE`, default)
  parameterization. Non-centered is more efficient when heterogeneity is
  small.

- verbose:

  Logical; print progress messages (default `TRUE`).

- ...:

  Additional arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html):

  - `chains`: Number of chains (default 4)

  - `iter`: Iterations per chain (default 2000)

  - `warmup`: Warmup iterations (default iter/2)

  - `cores`: Cores for parallel sampling

  - `seed`: Random seed

  - `control`: List of sampler controls (e.g.,
    `list(adapt_delta = 0.95)`)

## Value

A `shrinkr_fit` object (list) containing:

- fit:

  Stan model object

- data:

  Data list used for fitting

- summary:

  Parameter summaries (mean, sd, quantiles, Rhat, ESS)

- diagnostics:

  Sampler diagnostics (divergences, treedepth)

- priors:

  Prior specifications used

## Details

### Model Specification

**Hierarchical model (Stage 2):** \$\$\theta_g \mid \mu, \tau \sim
\text{Normal}(\mu, \tau^2), \quad g = 1, \ldots, G\$\$ \$\$\mu \sim
\pi(\mu)\$\$ \$\$\tau \sim \pi(\tau)\$\$

**Stage 1 likelihood (Gaussian mixture approximation):** \$\$\theta_g
\mid D_g \sim q_g(\theta_g) \approx \sum\_{k=1}^K w_k \\
\text{MVN}(\mu_k, \Sigma_k)\$\$

**Full posterior:** \$\$\pi(\theta, \mu, \tau \mid D) \propto
\left\[\prod\_{g=1}^G q_g(\theta_g)\right\] \left\[\prod\_{g=1}^G
\text{Normal}(\theta_g \mid \mu, \tau^2)\right\] \pi(\mu) \pi(\tau)\$\$

where \\q_g(\theta_g)\\ approximates the Stage 1 posterior for group
\\g\\.

### What's Fixed vs. Flexible

**Fixed:**

- Hierarchical distribution: \\\theta_g \mid \mu, \tau \sim
  \text{Normal}(\mu, \tau^2)\\

**Flexible:**

- Hyperpriors \\\pi(\mu)\\ and \\\pi(\tau)\\: Normal, Student-t, Cauchy,
  Lognormal, Gamma, Inverse-Gamma, Exponential, Uniform, mixtures
  (including spike-and-slab), and truncated versions

- Stage 1 posteriors: Can be non-Normal (handled by mixture
  approximation)

### Critical Requirements

1.  **Stage 1 must use flat/uninformative priors** on \\\theta_g\\

    - Ensures two-stage = one-stage hierarchical model

    - Stan: Don't specify prior (defaults to flat)

    - JAGS/NIMBLE: Use very wide priors

2.  **Verify mixture quality:** `plot(mixture, draws = samples)`

    - Check density overlays and QQ plots

    - Poor approximation → biased shrinkage

3.  **Check prior implications:**
    `sample_prior_predictive(hierarchical_priors)`

    - Understand what priors imply before fitting

    - Avoid prior-data conflicts

4.  **Minimum 2 groups required** for heterogeneity estimation

### Common Prior Choices for \\\tau\\

- Half-Normal: `dist_truncated(dist_normal(0, s), lower = 0)` - Weakly
  informative

- Half-t: `dist_truncated(dist_student_t(df, 0, s), lower = 0)` -
  Heavier tails

- Half-Cauchy: `dist_truncated(dist_cauchy(0, s), lower = 0)` - Very
  diffuse

- Uniform: `dist_uniform(0, U)` - Bounded heterogeneity

- Inverse-Gamma: `dist_inverse_gamma(a, b)` - Traditional choice

See [`vignette("getting_started")`](../articles/getting_started.md) for
complete workflow,
[`vignette("brms_integration")`](../articles/brms_integration.md) for
real examples, and package README for mathematical justification.

## See also

**Workflow functions:** [`fit_mixture()`](fit_mixture.md),
[`sample_prior_predictive()`](sample_prior_predictive.md)

**Extract results:** [`extract_mu_tau()`](extract_mu_tau.md),
[`extract_theta()`](extract_theta.md),
[`summarize_mu_tau()`](summarise_mu_tau.md),
[`summarize_theta()`](summarise_theta.md),
[`theta_contrasts()`](theta_contrasts.md)

**Visualization:** [`plot.shrinkr_fit()`](plot.shrinkr_fit.md),
[`plot.shrinkr_mixture()`](plot.shrinkr_mixture.md)

**Vignettes:**

- [`vignette("getting_started")`](../articles/getting_started.md) -
  Complete workflow with Stan example

- [`vignette("brms_integration")`](../articles/brms_integration.md) -
  Survival analysis example

## Examples

``` r
if (FALSE) { # \dontrun{
library(distributional)

# Simulate Stage 1 posteriors
samples <- list(
  group1 = matrix(rnorm(2000, 0.0, 0.5), ncol = 1),
  group2 = matrix(rnorm(2000, 0.5, 0.5), ncol = 1),
  group3 = matrix(rnorm(2000, 1.0, 0.5), ncol = 1)
)

# Fit mixture approximation
mix <- fit_mixture(samples, K_max = 3)
plot(mix, draws = samples)  # Check quality

# Apply shrinkage with default priors
fit <- shrink(mixture = mix)
print(fit)
plot(fit)

# Custom priors (half-t for heavier tails)
priors <- list(
  mu = dist_normal(0, 10),
  tau = dist_truncated(dist_student_t(3, 0, 2.5), lower = 0)
)
fit2 <- shrink(mixture = mix, hierarchical_priors = priors)

# Using point estimates only
fit3 <- shrink(
  mle = c(0.5, 1.2, -0.3),
  var_matrix = c(0.1, 0.15, 0.12),
  hierarchical_priors = priors
)

# Extract results
mu_tau <- extract_mu_tau(fit)
mu_tau_summary <- summarize_mu_tau(fit)
theta_summary <- summarize_theta(fit)
} # }
```
