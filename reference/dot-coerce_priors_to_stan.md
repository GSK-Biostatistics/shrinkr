# Convert distributional priors to Stan format

Handles Normal, Student-t, Cauchy, Gamma, Exponential, Lognormal,
Inverse-Gamma, Uniform, truncated distributions, and mixture priors
(including spike-and-slab) for both mu and tau.

## Usage

``` r
.coerce_priors_to_stan(prior_mu, prior_tau)
```
