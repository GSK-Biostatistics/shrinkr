# Print method for shrinkr_fit

Displays a compact summary of the fitted model including dimensions,
hyperparameter estimates, and diagnostics.

## Usage

``` r
# S3 method for class 'shrinkr_fit'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `shrinkr_fit` object.

- digits:

  Number of digits to display. Default is 3.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

## See also

[`shrink()`](shrink.md) for fitting models,
[`summarise_theta()`](summarise_theta.md) for detailed theta estimates
