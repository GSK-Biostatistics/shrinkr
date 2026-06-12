# Print method for mixture fits

Displays summary information about a fitted mixture model in a readable
format.

## Usage

``` r
# S3 method for class 'shrinkr_mixture'
print(x, ...)
```

## Arguments

- x:

  A `shrinkr_mixture` object from [`fit_mixture`](fit_mixture.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## See also

[`fit_mixture`](fit_mixture.md) for fitting mixture models,
[`summary.shrinkr_mixture`](summary.shrinkr_mixture.md) for detailed
summaries

## Examples

``` r
if (FALSE) { # \dontrun{
samples <- list(
  matrix(rnorm(2000, 0.0, 0.5), ncol = 1),
  matrix(rnorm(2000, 0.5, 0.5), ncol = 1),
  matrix(rnorm(2000, 1.0, 0.5), ncol = 1)
)

mix <- fit_mixture(samples, K_max = 5)

# Print summary
print(mix)

# Or just type the object name
mix
} # }
```
