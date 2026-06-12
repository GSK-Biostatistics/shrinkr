# Fit Gaussian mixture models to posterior samples

Fits a multivariate Gaussian mixture model (GMM) jointly across all
supplied variables using mclust. The function is intended for
approximating posterior draws from Bayesian models and produces a tidy
component table suitable for visualization and shrinkage modeling.

## Usage

``` r
fit_mixture(samples, K_max = 5L, verbose = FALSE, model_names = NULL, ...)
```

## Arguments

- samples:

  Posterior samples in one of the following formats:

  - **Data frame or matrix** (recommended): Each column represents one
    group/variable to shrink, with rows as posterior draws. Column names
    are used as variable labels. Example: output from
    [`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
    or a matrix where columns = groups.

  - **Named list of vectors**: Each list element contains posterior
    samples for one group (as a numeric vector). Names are used as
    variable labels. All vectors must have the same length.

  - **Named list of matrices**: Each list element is a matrix of
    posterior samples for one group. For univariate parameters, these
    should be single-column matrices (n × 1). All matrices must have the
    same number of rows (draws).

  At least **two** groups/variables are required for hierarchical
  shrinkage. Non-numeric columns in data frames are automatically
  dropped. Rows with missing values are removed before fitting.

- K_max:

  Integer (≥ 1). Maximum number of mixture components to consider during
  model selection by BIC. Internally capped at `n - 1` for stability.

- verbose:

  Logical. If `TRUE`, progress and diagnostic messages are printed.

- model_names:

  Optional character vector of mclust covariance model codes to consider
  (e.g., `"EII"`, `"VVV"`, etc.). If `NULL` (default), all models
  appropriate for the data dimension are considered by mclust. This is a
  convenience wrapper for `modelNames` in
  [`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html).

- ...:

  Additional arguments forwarded to
  [`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html),
  such as `prior`, `initialization`, `control = mclust::emControl()`,
  `warn`, or `verbose`. Use these to fine-tune EM control, priors, or
  initialization.

## Value

A list with class "shrinkr_mixture" containing:

- `components` — data frame with columns: `group`, `component`,
  `variable`, `weight`, `mean`, and `sd` (marginal standard deviations).

- `K` — number of mixture components selected.

- `vars` — character vector of variable names.

- `weights` — vector of component weights (mixing proportions).

- `covs` — list of component covariance matrices (p × p each).

- `model_name` — selected mclust covariance structure (e.g., "VVV").

- `bic` — Bayesian Information Criterion for the fitted model.

- `n_samples` — number of samples used in fitting.

- `n_vars` — number of variables.

- `mclust_fit` — the complete mclust model object (for advanced use).

- `diagnostics` — list with sample size details, removed rows, and
  quality warnings.

## Details

- Requires at least **two** variables; shrinkage across a single
  variable is not meaningful.

- Rows with any missing values are removed before fitting.

- Component weights are normalized to sum to one.

- The `sd` column reports marginal standard deviations (square roots of
  diagonal entries) from each component covariance matrix.

## Covariance structures in mclust

mclust parameterizes component covariances via eigen-decomposition and
offers a set of model families controlling volume (V), shape (S), and
orientation (O). Common codes include (non-exhaustive):

- **Spherical**:

  `"EII"`: equal volume, spherical  
  `"VII"`: variable volume, spherical

- **Diagonal**:

  `"EEI"`: equal volume & shape (axis-aligned)  
  `"VEI"`: variable volume, equal shape  
  `"EVI"`: equal volume, variable shape  
  `"VVI"`: variable volume & shape

- **Ellipsoidal (full covariance)**:

  `"EEE"`: equal volume, shape, orientation  
  `"EEV"`: equal volume & shape, variable orientation  
  `"VEV"`: variable volume, equal shape, variable orientation  
  `"VVV"`: variable volume, shape, and orientation (most flexible)

If `model_names = NULL` (default),
[`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html)
selects among the models appropriate for the data dimension via BIC. In
practice, this lets the data decide between parsimonious structures
(e.g., `"EII"`, `"VVI"`) and fully flexible ones (e.g., `"VVV"`). You
can restrict or expand the search space by supplying `model_names`.

## See also

[`plot.shrinkr_mixture`](plot.shrinkr_mixture.md) for visualizing
marginal fits,
[`as.data.frame.shrinkr_mixture`](as.data.frame.shrinkr_mixture.md) for
extracting component data

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with three groups
samples <- list(
  group1 = matrix(rnorm(2000, 0.0, 0.5), ncol = 1),
  group2 = matrix(rnorm(2000, 0.5, 0.5), ncol = 1),
  group3 = matrix(rnorm(2000, 1.0, 0.5), ncol = 1)
)

# Let mclust choose K and covariance model via BIC:
mix <- fit_mixture(samples, K_max = 5, verbose = TRUE)

# View the fit
mix
summary(mix)
plot(mix, draws = samples)

# Restrict to diagonal models only:
mix_diag <- fit_mixture(samples, K_max = 5, model_names = c("EEI","VVI"))

# Pass EM controls and initialization through ...
mix_tuned <- fit_mixture(
  samples, K_max = 5,
  control = mclust::emControl(eps = 1e-6, itmax = 5e3),
  initialization = list(hcPairs = TRUE)
)
} # }
```
