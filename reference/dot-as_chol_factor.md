# Internal: convert covariance to Cholesky factor for Stan

Uses a jitter-and-retry strategy when the input matrix is not quite
positive-definite. Falls back to
[`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html) if
available, and errors rather than passing a non-lower-triangular factor
to Stan.

## Usage

``` r
.as_chol_factor(Sigma)
```
