# Print summary of mixture fit

Displays formatted summary statistics for a fitted mixture model.

## Usage

``` r
# S3 method for class 'summary.shrinkr_mixture'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A summary object from
  [`summary.shrinkr_mixture`](summary.shrinkr_mixture.md).

- digits:

  Number of decimal digits to display. Default is 3.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## See also

[`summary.shrinkr_mixture`](summary.shrinkr_mixture.md) for creating
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
summ <- summary(mix)

# Print summary
print(summ)

# With more decimal places
print(summ, digits = 4)
} # }
```
