# Working with shrinkr in the Tidy Bayesian Ecosystem

``` r

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  warning = FALSE,
  message = FALSE
)
```

## Overview

**shrinkr** is designed to work seamlessly with modern R workflows. This
vignette shows practical examples of using shrinkr with:

- **tidyverse**: dplyr, ggplot2, tidyr for data manipulation and
  visualization
- **posterior**: Working with MCMC draws, computing summaries and
  diagnostics
- **bayesplot**: MCMC diagnostic plots (trace plots, pairs plots, etc.)
- **tidybayes**: Tidy manipulation of Bayesian posteriors
- **ggdist**: Modern distribution visualizations

``` r

library(shrinkr)
library(distributional)
library(MASS)

# Bayesian ecosystem
library(posterior)
library(bayesplot)
library(tidybayes)
library(ggdist)

# Tidyverse
library(dplyr)
library(tidyr)
library(ggplot2)

theme_set(theme_minimal(base_size = 12))
```

## Example: Multi-Region Clinical Trial

Imagine a clinical trial run across 5 regions testing a new treatment.
We have Stage 1 posterior samples from region-specific analyses.

### Simulate Stage 1 Results

``` r

set.seed(1104)

# True effects (unknown in practice)
true_effects <- c(0.45, 0.60, 0.38, -0.10, 0.65)
region_names <- c("North", "South", "East", "West", "Central")

# Simulate posterior samples from Stage 1
samples_list <- lapply(1:5, function(i) {
  matrix(rnorm(2000, true_effects[i], 0.20), ncol = 1)
})
names(samples_list) <- region_names
```

### Fit shrinkr Model

``` r

# Fit mixture approximation
mix <- fit_mixture(samples_list, K_max = 3, verbose = FALSE)

# Specify hierarchical priors
priors <- list(
  mu = dist_normal(0, 5),
  tau = dist_truncated(dist_student_t(3, 0, 1), lower = 0)
)

# Run hierarchical shrinkage
fit <- shrink(
  mixture = mix,
  hierarchical_priors = priors,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 1,
  seed = 2024,
  refresh = 0
)
#> 
#> SAMPLING FOR MODEL 'stage2_shrinkage' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.14 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  100 / 2000 [  5%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  300 / 2000 [ 15%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  500 / 2000 [ 25%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  700 / 2000 [ 35%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration:  900 / 2000 [ 45%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1100 / 2000 [ 55%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1300 / 2000 [ 65%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1500 / 2000 [ 75%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1700 / 2000 [ 85%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1900 / 2000 [ 95%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.054 seconds (Warm-up)
#> Chain 1:                0.048 seconds (Sampling)
#> Chain 1:                0.102 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'stage2_shrinkage' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 5e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  100 / 2000 [  5%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  300 / 2000 [ 15%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  500 / 2000 [ 25%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  700 / 2000 [ 35%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration:  900 / 2000 [ 45%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1100 / 2000 [ 55%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1300 / 2000 [ 65%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1500 / 2000 [ 75%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1700 / 2000 [ 85%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1900 / 2000 [ 95%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.057 seconds (Warm-up)
#> Chain 2:                0.038 seconds (Sampling)
#> Chain 2:                0.095 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'stage2_shrinkage' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 5e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  100 / 2000 [  5%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  300 / 2000 [ 15%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  500 / 2000 [ 25%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  700 / 2000 [ 35%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration:  900 / 2000 [ 45%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1100 / 2000 [ 55%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1300 / 2000 [ 65%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1500 / 2000 [ 75%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1700 / 2000 [ 85%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 1900 / 2000 [ 95%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.052 seconds (Warm-up)
#> Chain 3:                0.046 seconds (Sampling)
#> Chain 3:                0.098 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'stage2_shrinkage' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 5e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  100 / 2000 [  5%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  300 / 2000 [ 15%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  500 / 2000 [ 25%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  700 / 2000 [ 35%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration:  900 / 2000 [ 45%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1100 / 2000 [ 55%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1300 / 2000 [ 65%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1500 / 2000 [ 75%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1700 / 2000 [ 85%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 1900 / 2000 [ 95%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.053 seconds (Warm-up)
#> Chain 4:                0.047 seconds (Sampling)
#> Chain 4:                0.1 seconds (Total)
#> Chain 4:
```

## Working with posterior Package

The **posterior** package provides the foundation for working with MCMC
draws.

### Extract Draws

``` r

# Extract all parameters as draws_df
draws <- as_draws_df(fit)

# See what's available
variables(draws)
#> [1] "mu"          "tau"         "theta[1]"    "theta[2]"    "theta[3]"   
#> [6] "theta[4]"    "theta[5]"    "tau_squared" "lp__"

# Extract specific parameters
mu_tau_draws <- extract_mu_tau(fit)
theta_draws <- extract_theta(fit)
```

### Basic Summaries

``` r

# Quick summary of all parameters
summarize_draws(draws)
#> # A tibble: 9 × 10
#>   variable     mean  median    sd    mad       q5    q95  rhat ess_bulk ess_tail
#>   <chr>       <dbl>   <dbl> <dbl>  <dbl>    <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 mu         0.395   0.391  0.190 0.160   9.20e-2  0.710 1.00     1071.    1044.
#> 2 tau        0.315   0.267  0.222 0.173   4.25e-2  0.747 1.00      806.    1000.
#> 3 theta[1]   0.432   0.430  0.165 0.156   1.66e-1  0.717 1.000    4168.    2974.
#> 4 theta[2]   0.513   0.502  0.170 0.172   2.51e-1  0.805 1.00     3332.    3620.
#> 5 theta[3]   0.376   0.375  0.169 0.161   9.25e-2  0.653 1.00     5308.    3523.
#> 6 theta[4]   0.0997  0.105  0.210 0.221  -2.60e-1  0.427 1.00     1788.    2266.
#> 7 theta[5]   0.541   0.533  0.178 0.179   2.69e-1  0.850 1.00     2889.    3636.
#> 8 tau_squa…  0.148   0.0711 0.236 0.0811  1.81e-3  0.559 1.00      806.    1000.
#> 9 lp__      -6.47   -6.09   3.11  3.14   -1.22e+1 -2.13  1.00      936.    1611.

# Focus on theta parameters
summarize_draws(theta_draws, mean, sd, median, mad, ~quantile(.x, c(0.025, 0.975)))
#> # A tibble: 19 × 7
#>    variable       mean    sd  median    mad     `2.5%` `97.5%`
#>    <chr>         <dbl> <dbl>   <dbl>  <dbl>      <dbl>   <dbl>
#>  1 mu           0.395  0.190  0.391  0.160    0.0186     0.788
#>  2 tau          0.315  0.222  0.267  0.173    0.0199     0.892
#>  3 theta_c[1]  -0.0126 1.00  -0.0199 0.996   -1.97       1.94 
#>  4 theta_c[2]  -0.0244 1.01  -0.0259 1.00    -2.01       1.93 
#>  5 theta_c[3]   0.0116 1.02   0.0182 1.00    -2.03       1.96 
#>  6 theta_c[4]  -0.0210 0.984 -0.0255 0.989   -1.96       1.89 
#>  7 theta_c[5]   0.0175 0.999  0.0294 1.02    -1.95       1.89 
#>  8 z[1]         0.147  0.735  0.123  0.699   -1.28       1.63 
#>  9 z[2]         0.411  0.750  0.397  0.703   -1.07       1.92 
#> 10 z[3]        -0.0603 0.741 -0.0496 0.727   -1.52       1.42 
#> 11 z[4]        -0.985  0.776 -0.993  0.717   -2.48       0.528
#> 12 z[5]         0.518  0.721  0.502  0.699   -0.850      1.96 
#> 13 theta[1]     0.432  0.165  0.430  0.156    0.116      0.766
#> 14 theta[2]     0.513  0.170  0.502  0.172    0.200      0.864
#> 15 theta[3]     0.376  0.169  0.375  0.161    0.0323     0.710
#> 16 theta[4]     0.0997 0.210  0.105  0.221   -0.320      0.477
#> 17 theta[5]     0.541  0.178  0.533  0.179    0.224      0.912
#> 18 tau_squared  0.148  0.236  0.0711 0.0811   0.000396   0.796
#> 19 lp__        -6.47   3.11  -6.09   3.14   -13.6       -1.54

# Convergence diagnostics
summarize_draws(draws, default_convergence_measures())
#> # A tibble: 9 × 4
#>   variable     rhat ess_bulk ess_tail
#>   <chr>       <dbl>    <dbl>    <dbl>
#> 1 mu          1.00     1071.    1044.
#> 2 tau         1.00      806.    1000.
#> 3 theta[1]    1.000    4168.    2974.
#> 4 theta[2]    1.00     3332.    3620.
#> 5 theta[3]    1.00     5308.    3523.
#> 6 theta[4]    1.00     1788.    2266.
#> 7 theta[5]    1.00     2889.    3636.
#> 8 tau_squared 1.00      806.    1000.
#> 9 lp__        1.00      936.    1611.

# Custom summaries
summarise_draws(
  theta_draws,
  mean,
  sd,
  prob_positive = ~mean(.x > 0),
  prob_large = ~mean(.x > 0.5)
)
#> # A tibble: 19 × 5
#>    variable       mean    sd prob_positive prob_large
#>    <chr>         <dbl> <dbl>         <dbl>      <dbl>
#>  1 mu           0.395  0.190        0.978      0.255 
#>  2 tau          0.315  0.222        1          0.162 
#>  3 theta_c[1]  -0.0126 1.00         0.492      0.303 
#>  4 theta_c[2]  -0.0244 1.01         0.490      0.302 
#>  5 theta_c[3]   0.0116 1.02         0.507      0.318 
#>  6 theta_c[4]  -0.0210 0.984        0.492      0.301 
#>  7 theta_c[5]   0.0175 0.999        0.512      0.321 
#>  8 z[1]         0.147  0.735        0.574      0.303 
#>  9 z[2]         0.411  0.750        0.72       0.444 
#> 10 z[3]        -0.0603 0.741        0.472      0.218 
#> 11 z[4]        -0.985  0.776        0.0905     0.0265
#> 12 z[5]         0.518  0.721        0.771      0.501 
#> 13 theta[1]     0.432  0.165        0.995      0.324 
#> 14 theta[2]     0.513  0.170        0.999      0.504 
#> 15 theta[3]     0.376  0.169        0.982      0.224 
#> 16 theta[4]     0.0997 0.210        0.680      0.016 
#> 17 theta[5]     0.541  0.178        1.000      0.572 
#> 18 tau_squared  0.148  0.236        1          0.058 
#> 19 lp__        -6.47   3.11         0.0025     0.0005
```

### Check Convergence

``` r

# Check Rhat for all parameters
all_rhats <- summarise_draws(draws, "rhat")
max(all_rhats$rhat, na.rm = TRUE)
#> [1] 1.003493

# Check effective sample size
summarise_draws(draws, "ess_bulk", "ess_tail") %>%
  filter(ess_bulk < 400 | ess_tail < 400)
#> # A tibble: 0 × 3
#> # ℹ 3 variables: variable <chr>, ess_bulk <dbl>, ess_tail <dbl>

# Detailed diagnostics for specific parameters
summarise_draws(
  subset_draws(draws, variable = c("mu", "tau")),
  default_convergence_measures()
)
#> # A tibble: 2 × 4
#>   variable  rhat ess_bulk ess_tail
#>   <chr>    <dbl>    <dbl>    <dbl>
#> 1 mu        1.00    1071.    1044.
#> 2 tau       1.00     806.    1000.
```

## Diagnostic Plots with bayesplot

**bayesplot** provides essential MCMC diagnostic visualizations.

### Trace Plots

Check for mixing and stationarity:

``` r

# Check hyperparameters
mcmc_trace(draws, pars = c("mu", "tau", "tau_squared"))
```

![](tidy_bayesian_workflow_files/figure-html/trace_plots-1.png)

``` r


# Check first few thetas
mcmc_trace(draws, regex_pars = "theta\\[[1-3]\\]")
```

![](tidy_bayesian_workflow_files/figure-html/trace_plots-2.png)

``` r


# All thetas at once (if not too many)
mcmc_trace(draws, regex_pars = "theta")
```

![](tidy_bayesian_workflow_files/figure-html/trace_plots-3.png)

### Density Plots

Compare chains and check for multimodality:

``` r

# Overlay densities from different chains
mcmc_dens_overlay(draws, pars = c("mu", "tau"))
```

![](tidy_bayesian_workflow_files/figure-html/density_plots-1.png)

``` r


# Individual densities
mcmc_dens(draws, pars = c("mu", "tau", "tau_squared"))
```

![](tidy_bayesian_workflow_files/figure-html/density_plots-2.png)

``` r


# Compare all thetas
mcmc_dens_overlay(draws, regex_pars = "theta")
```

![](tidy_bayesian_workflow_files/figure-html/density_plots-3.png)

### Interval Plots

Visualize posterior uncertainties:

``` r

# All thetas with 50% and 95% intervals
mcmc_intervals(draws, regex_pars = "theta", prob = 0.5, prob_outer = 0.95)
```

![](tidy_bayesian_workflow_files/figure-html/interval_plots-1.png)

``` r


# With point estimates
mcmc_intervals_data(draws, regex_pars = "theta") %>%
  ggplot(aes(y = parameter)) +
  geom_pointrange(aes(x = m, xmin = ll, xmax = hh)) +
  geom_point(aes(x = m), size = 3) +
  labs(title = "Posterior Intervals for Regional Effects", x = "Effect Size", y = NULL)
```

![](tidy_bayesian_workflow_files/figure-html/interval_plots-2.png)

### Area Plots

Density plots with shaded intervals:

``` r

# Hyperparameters
mcmc_areas(draws, pars = c("mu", "tau"), prob = 0.95, prob_outer = 0.99)
```

![](tidy_bayesian_workflow_files/figure-html/area_plots-1.png)

``` r


# All thetas
mcmc_areas(draws, regex_pars = "theta", prob = 0.8)
```

![](tidy_bayesian_workflow_files/figure-html/area_plots-2.png)

## Tidy Analysis with tidybayes

**tidybayes** makes it easy to manipulate and visualize posteriors using
tidy principles.

### Spread and Gather Draws

``` r

# Gather theta parameters into long format
theta_tidy <- draws %>%
  gather_draws(theta[region]) %>%
  mutate(region = region_names[region])

head(theta_tidy)
#> # A tibble: 6 × 6
#> # Groups:   region, .variable [1]
#>   region .chain .iteration .draw .variable .value
#>   <chr>   <int>      <int> <int> <chr>      <dbl>
#> 1 North       1          1     1 theta      0.432
#> 2 North       1          2     2 theta      0.439
#> 3 North       1          3     3 theta      0.390
#> 4 North       1          4     4 theta      0.443
#> 5 North       1          5     5 theta      0.246
#> 6 North       1          6     6 theta      0.623

# Spread into wide format
theta_wide <- draws %>%
  spread_draws(theta[region]) %>%
  mutate(region = region_names[region])

head(theta_wide)
#> # A tibble: 6 × 5
#> # Groups:   region [1]
#>   region theta .chain .iteration .draw
#>   <chr>  <dbl>  <int>      <int> <int>
#> 1 North  0.432      1          1     1
#> 2 North  0.439      1          2     2
#> 3 North  0.390      1          3     3
#> 4 North  0.443      1          4     4
#> 5 North  0.246      1          5     5
#> 6 North  0.623      1          6     6
```

### Point and Interval Summaries

``` r

# Median and 95% quantile intervals
theta_tidy %>%
  group_by(region) %>%
  median_qi(.value, .width = 0.95)
#> # A tibble: 5 × 7
#>   region  .value  .lower .upper .width .point .interval
#>   <chr>    <dbl>   <dbl>  <dbl>  <dbl> <chr>  <chr>    
#> 1 Central  0.533  0.224   0.912   0.95 median qi       
#> 2 East     0.375  0.0323  0.710   0.95 median qi       
#> 3 North    0.430  0.116   0.766   0.95 median qi       
#> 4 South    0.502  0.200   0.864   0.95 median qi       
#> 5 West     0.105 -0.320   0.477   0.95 median qi

# Multiple interval widths
theta_tidy %>%
  group_by(region) %>%
  median_qi(.value, .width = c(0.5, 0.8, 0.95))
#> # A tibble: 15 × 7
#>    region  .value  .lower .upper .width .point .interval
#>    <chr>    <dbl>   <dbl>  <dbl>  <dbl> <chr>  <chr>    
#>  1 Central  0.533  0.415   0.656   0.5  median qi       
#>  2 East     0.375  0.268   0.486   0.5  median qi       
#>  3 North    0.430  0.326   0.536   0.5  median qi       
#>  4 South    0.502  0.395   0.629   0.5  median qi       
#>  5 West     0.105 -0.0442  0.254   0.5  median qi       
#>  6 Central  0.533  0.324   0.778   0.8  median qi       
#>  7 East     0.375  0.166   0.592   0.8  median qi       
#>  8 North    0.430  0.230   0.645   0.8  median qi       
#>  9 South    0.502  0.301   0.734   0.8  median qi       
#> 10 West     0.105 -0.175   0.372   0.8  median qi       
#> 11 Central  0.533  0.224   0.912   0.95 median qi       
#> 12 East     0.375  0.0323  0.710   0.95 median qi       
#> 13 North    0.430  0.116   0.766   0.95 median qi       
#> 14 South    0.502  0.200   0.864   0.95 median qi       
#> 15 West     0.105 -0.320   0.477   0.95 median qi

# Mean and HDI (highest density interval)
theta_tidy %>%
  group_by(region) %>%
  mean_hdi(.value, .width = 0.95)
#> # A tibble: 5 × 7
#>   region  .value  .lower .upper .width .point .interval
#>   <chr>    <dbl>   <dbl>  <dbl>  <dbl> <chr>  <chr>    
#> 1 Central 0.541   0.210   0.899   0.95 mean   hdi      
#> 2 East    0.376   0.0340  0.714   0.95 mean   hdi      
#> 3 North   0.432   0.108   0.756   0.95 mean   hdi      
#> 4 South   0.513   0.195   0.851   0.95 mean   hdi      
#> 5 West    0.0997 -0.297   0.502   0.95 mean   hdi
```

### Custom Summaries with dplyr

``` r

# Probability of positive effect
theta_tidy %>%
  group_by(region) %>%
  summarise(
    mean_effect = mean(.value),
    sd_effect = sd(.value),
    prob_positive = mean(.value > 0),
    prob_clinically_meaningful = mean(.value > 0.3),
    .groups = "drop"
  ) %>%
  arrange(desc(prob_positive))
#> # A tibble: 5 × 5
#>   region  mean_effect sd_effect prob_positive prob_clinically_meaningful
#>   <chr>         <dbl>     <dbl>         <dbl>                      <dbl>
#> 1 Central      0.541      0.178         1.000                      0.923
#> 2 South        0.513      0.170         0.999                      0.901
#> 3 North        0.432      0.165         0.995                      0.802
#> 4 East         0.376      0.169         0.982                      0.676
#> 5 West         0.0997     0.210         0.680                      0.184
```

### Computing Contrasts

``` r

# Method 1: Using shrinkr's built-in function
L <- rbind(
  "South - North" = c(-1, 1, 0, 0, 0),
  "Central - North" = c(-1, 0, 0, 0, 1),
  "South - West" = c(0, 1, 0, -1, 0)
)
contrasts <- theta_contrasts(fit, L, labels = rownames(L))
summarise_draws(contrasts)
#> # A tibble: 3 × 10
#>   variable         mean median    sd   mad      q5   q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 South - North  0.0804 0.0630 0.218 0.195 -0.257  0.462  1.00    4920.    3318.
#> 2 Central - Nor… 0.109  0.0928 0.222 0.200 -0.241  0.489  1.00    3635.    3398.
#> 3 South - West   0.413  0.404  0.286 0.315 -0.0105 0.904  1.00    1602.    1993.

# Method 2: Using tidybayes compare_levels
theta_wide %>%
  compare_levels(theta, by = region, comparison = "pairwise") %>%
  group_by(region) %>%
  median_qi(theta) %>%
  arrange(desc(theta))
#> # A tibble: 10 × 7
#>    region            theta .lower .upper .width .point .interval
#>    <chr>             <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
#>  1 South - East     0.112  -0.275 0.613    0.95 median qi       
#>  2 South - North    0.0630 -0.332 0.543    0.95 median qi       
#>  3 North - East     0.0419 -0.397 0.519    0.95 median qi       
#>  4 South - Central -0.0168 -0.483 0.403    0.95 median qi       
#>  5 North - Central -0.0928 -0.577 0.317    0.95 median qi       
#>  6 East - Central  -0.138  -0.668 0.245    0.95 median qi       
#>  7 West - East     -0.253  -0.816 0.132    0.95 median qi       
#>  8 West - North    -0.317  -0.892 0.0860   0.95 median qi       
#>  9 West - South    -0.404  -0.989 0.0471   0.95 median qi       
#> 10 West - Central  -0.430  -1.06  0.0286   0.95 median qi
```

## Modern Visualizations with ggdist

**ggdist** provides publication-ready distribution visualizations.

### Halfeye Plots

Eye + interval visualization:

``` r

theta_tidy %>%
  ggplot(aes(y = region, x = .value)) +
  stat_halfeye(
    .width = c(0.66, 0.95),
    fill = "steelblue"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Regional Treatment Effects",
    subtitle = "Posterior distributions with median and 66%/95% intervals",
    x = "Treatment Effect",
    y = NULL
  )
```

![](tidy_bayesian_workflow_files/figure-html/halfeye-1.png)

### Slab + Interval

Density with separate interval layer:

``` r

theta_tidy %>%
  ggplot(aes(y = region, x = .value)) +
  stat_slab(aes(fill_ramp = after_stat(level)), fill = "steelblue", alpha = 0.8) +
  stat_pointinterval(.width = c(0.66, 0.95), position = position_nudge(y = -0.15)) +
  scale_fill_ramp_discrete(range = c(1, 0.2), guide = "none") +
  labs(
    title = "Posterior Densities with Quantile Intervals",
    x = "Treatment Effect",
    y = NULL
  )
```

![](tidy_bayesian_workflow_files/figure-html/slab_interval-1.png)

### Quantile Dotplots

Each dot = quantile of the distribution:

``` r

theta_tidy %>%
  ggplot(aes(y = region, x = .value)) +
  stat_dots(quantiles = 100) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Quantile Dotplots",
    subtitle = "Each dot represents 1% of the posterior",
    x = "Treatment Effect",
    y = NULL
  )
```

![](tidy_bayesian_workflow_files/figure-html/dotplot-1.png)

### Gradient Intervals

Continuous representation of uncertainty:

``` r

theta_tidy %>%
  ggplot(aes(y = region, x = .value)) +
  stat_gradientinterval(.width = ppoints(50)) +
  scale_color_brewer(palette = "Blues", guide = "none") +
  labs(
    title = "Gradient Interval Representation",
    x = "Treatment Effect",
    y = NULL
  )
```

![](tidy_bayesian_workflow_files/figure-html/gradient-1.png)

## Comparing Pre- and Post-Shrinkage

### Extract Both Estimates

``` r

# Get pre-shrunk estimates from mixture
pre_shrunk <- summarise_theta(fit) %>%
  mutate(type = "Pre-shrunk")

# Get post-shrunk estimates
post_shrunk <- summarise_theta(fit) %>%
  mutate(type = "Post-shrunk")

# Or use shrinkr's built-in plot
plot(fit, group_names = region_names)
```

![](tidy_bayesian_workflow_files/figure-html/compare_shrinkage-1.png)

### Custom Comparison Plot

``` r

# Get the hierarchical mean (mu)
mu_draws <- draws %>% spread_draws(mu)
mu_mean <- mean(mu_draws$mu)

# Combine with Stage 1 samples
stage1_draws <- lapply(seq_along(samples_list), function(i) {
  data.frame(
    region = region_names[i],
    .value = samples_list[[i]][,1],
    type = "Stage 1"
  )
}) %>% bind_rows()

stage2_draws <- theta_tidy %>%
  mutate(type = "Stage 2 (Shrunk)")

# Plot side by side
bind_rows(stage1_draws, stage2_draws) %>%
  ggplot(aes(y = region, x = .value, fill = type)) +
  stat_halfeye(alpha = 0.7, position = position_dodge(width = 0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  geom_vline(xintercept = mu_mean, linetype = "solid", alpha = 0.8, 
             color = "darkred", linewidth = 1) +
  annotate("text", x = mu_mean, y = 0.5, 
           label = sprintf("Global mean (μ) = %.2f", mu_mean),
           hjust = -0.1, color = "darkred", size = 3.5) +
  scale_fill_manual(values = c("Stage 1" = "gray70", "Stage 2 (Shrunk)" = "steelblue")) +
  labs(
    title = "Stage 1 vs Stage 2: Effect of Hierarchical Shrinkage",
    subtitle = "Stage 2 estimates are pulled toward the global mean",
    x = "Treatment Effect",
    y = NULL,
    fill = NULL
  ) +
  theme(legend.position = "bottom")
```

![](tidy_bayesian_workflow_files/figure-html/custom_comparison-1.png)

## Complete Workflow Example

Here’s a typical analysis workflow using tidy principles:

``` r

# 1. Extract and prepare data
analysis_data <- draws %>%
  spread_draws(mu, tau, theta[i]) %>%
  mutate(region = region_names[i])

# 2. Compute summaries
summary_table <- analysis_data %>%
  group_by(region) %>%
  summarise(
    mean = mean(theta),
    median = median(theta),
    sd = sd(theta),
    q025 = quantile(theta, 0.025),
    q975 = quantile(theta, 0.975),
    prob_positive = mean(theta > 0),
    prob_clinically_important = mean(theta > 0.3),
    .groups = "drop"
  ) %>%
  arrange(desc(median))

print(summary_table)
#> # A tibble: 5 × 8
#>   region    mean median    sd    q025  q975 prob_positive prob_clinically_impo…¹
#>   <chr>    <dbl>  <dbl> <dbl>   <dbl> <dbl>         <dbl>                  <dbl>
#> 1 Central 0.541   0.533 0.178  0.224  0.912         1.000                  0.923
#> 2 South   0.513   0.502 0.170  0.200  0.864         0.999                  0.901
#> 3 North   0.432   0.430 0.165  0.116  0.766         0.995                  0.802
#> 4 East    0.376   0.375 0.169  0.0323 0.710         0.982                  0.676
#> 5 West    0.0997  0.105 0.210 -0.320  0.477         0.680                  0.184
#> # ℹ abbreviated name: ¹​prob_clinically_important

# 3. Create advanced figure
library(patchwork)

p1 <- analysis_data %>%
  ggplot(aes(y = reorder(region, theta), x = theta)) +
  stat_halfeye(.width = c(0.66, 0.95), fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "A. Regional Treatment Effects",
    x = "Effect Size",
    y = NULL
  )

p2 <- analysis_data %>%
  dplyr::ungroup() %>%
  dplyr::select(mu, tau, .draw) %>%
  dplyr::distinct() %>%
  tidyr::pivot_longer(cols = c(mu, tau), names_to = "name", values_to = "value") %>%
  ggplot(aes(x = value, fill = name)) +
  stat_halfeye(alpha = 0.7) +
  facet_wrap(~name, scales = "free", labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "B. Hyperparameters",
    x = "Value",
    y = "Density"
  ) +
  theme(legend.position = "none")

p3 <- analysis_data %>%
  dplyr::ungroup() %>%
  dplyr::select(.draw, region, theta) %>%
  compare_levels(theta, by = region) %>%
  ggplot(aes(y = region, x = theta)) +
  stat_halfeye(fill = "coral", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "C. Pairwise Regional Comparisons",
    x = "Difference in Effect Size",
    y = NULL
  )

p4 <- analysis_data %>%
  dplyr::ungroup() %>%
  dplyr::select(.draw, mu, tau) %>%
  dplyr::distinct() %>%
  ggplot(aes(x = mu, y = tau)) +
  geom_hex(bins = 30) +
  stat_ellipse(level = 0.95, color = "red", linewidth = 1) +
  scale_fill_viridis_c() +
  labs(
    title = "D. Hyperparameter Correlation",
    x = expression(mu~"(global mean)"),
    y = expression(tau~"(heterogeneity)")
  )

(p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Complete Bayesian Shrinkage Analysis",
    subtitle = sprintf(
      "Global effect: %.2f [%.2f, %.2f] | Heterogeneity (tau): %.2f",
      median(analysis_data$mu),
      quantile(analysis_data$mu, 0.025),
      quantile(analysis_data$mu, 0.975),
      median(analysis_data$tau)
    )
  )
```

![](tidy_bayesian_workflow_files/figure-html/complete_workflow-1.png)

## Advanced: Custom Analyses

### Probability Statements

``` r

# Which region is best?
analysis_data %>%
  group_by(.draw) %>%
  slice_max(theta, n = 1) %>%
  ungroup() %>%
  count(region) %>%
  mutate(probability = n / sum(n)) %>%
  arrange(desc(probability))
#> # A tibble: 5 × 3
#>   region      n probability
#>   <chr>   <int>       <dbl>
#> 1 Central  1571      0.393 
#> 2 South    1321      0.330 
#> 3 North     668      0.167 
#> 4 East      397      0.0992
#> 5 West       43      0.0107

# Alternative: probability each region is best
analysis_data %>%
  group_by(.draw) %>%
  mutate(rank = rank(-theta)) %>%
  ungroup() %>%
  group_by(region) %>%
  summarise(
    prob_best = mean(rank == 1),
    prob_top2 = mean(rank <= 2),
    mean_rank = mean(rank),
    .groups = "drop"
  ) %>%
  arrange(mean_rank)
#> # A tibble: 5 × 4
#>   region  prob_best prob_top2 mean_rank
#>   <chr>       <dbl>     <dbl>     <dbl>
#> 1 Central    0.393     0.69        2.06
#> 2 South      0.330     0.618       2.24
#> 3 North      0.167     0.395       2.81
#> 4 East       0.0992    0.260       3.23
#> 5 West       0.0107    0.0368      4.65

# Pairwise comparisons: Probability that South > North
# Create wide format for comparisons
theta_wide_for_contrasts <- analysis_data %>%
  ungroup() %>%
  dplyr::select(.draw, region, theta) %>%
  tidyr::pivot_wider(names_from = region, values_from = theta)

theta_wide_for_contrasts %>%
  summarise(
    prob_south_beats_north = mean(South > North),
    prob_south_beats_north_by_02 = mean((South - North) > 0.2),
    prob_central_beats_all = mean(
      Central > North & Central > South & 
      Central > East & Central > West
    )
  )
#> # A tibble: 1 × 3
#>   prob_south_beats_north prob_south_beats_north_by_02 prob_central_beats_all
#>                    <dbl>                        <dbl>                  <dbl>
#> 1                  0.636                        0.270                  0.393
```

### Tail Probabilities

``` r

# Classify effects into categories
theta_tidy %>%
  group_by(region) %>%
  summarise(
    prob_harm = mean(.value < -0.1),
    prob_null = mean(abs(.value) < 0.1),
    prob_small_benefit = mean(.value > 0.1 & .value < 0.3),
    prob_large_benefit = mean(.value > 0.3),
    .groups = "drop"
  ) %>%
  arrange(desc(prob_large_benefit))
#> # A tibble: 5 × 5
#>   region  prob_harm prob_null prob_small_benefit prob_large_benefit
#>   <chr>       <dbl>     <dbl>              <dbl>              <dbl>
#> 1 Central   0         0.003               0.0742              0.923
#> 2 South     0         0.00475             0.094               0.901
#> 3 North     0.001     0.0192              0.178               0.802
#> 4 East      0.00275   0.0508              0.271               0.676
#> 5 West      0.179     0.314               0.324               0.184

# Visualize classification
theta_tidy %>%
  mutate(
    category = case_when(
      .value < -0.1 ~ "Harm",
      abs(.value) < 0.1 ~ "Null",
      .value > 0.1 & .value < 0.3 ~ "Small Benefit",
      .value > 0.3 ~ "Large Benefit"
    )
  ) %>%
  count(region, category) %>%
  group_by(region) %>%
  mutate(probability = n / sum(n)) %>%
  ggplot(aes(x = probability, y = region, fill = category)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(
      "Harm" = "red",
      "Null" = "gray",
      "Small Benefit" = "lightblue",
      "Large Benefit" = "darkblue"
    )
  ) +
  labs(
    title = "Classification of Treatment Effects",
    x = "Probability",
    y = NULL,
    fill = "Effect Category"
  ) +
  theme(legend.position = "bottom")
```

![](tidy_bayesian_workflow_files/figure-html/tail_probs-1.png)

### Ranking Analysis

``` r

# Compute ranks for each draw
rank_data <- analysis_data %>%
  group_by(.draw) %>%
  mutate(rank = rank(-theta)) %>%
  ungroup()

# Summary statistics
rank_summary <- rank_data %>%
  group_by(region) %>%
  summarise(
    mean_rank = mean(rank),
    median_rank = median(rank),
    prob_rank1 = mean(rank == 1),
    prob_rank2 = mean(rank == 2),
    prob_top3 = mean(rank <= 3),
    .groups = "drop"
  ) %>%
  arrange(mean_rank)

print(rank_summary)
#> # A tibble: 5 × 6
#>   region  mean_rank median_rank prob_rank1 prob_rank2 prob_top3
#>   <chr>       <dbl>       <dbl>      <dbl>      <dbl>     <dbl>
#> 1 Central      2.06           2     0.393       0.297    0.876 
#> 2 South        2.24           2     0.330       0.288    0.841 
#> 3 North        2.81           3     0.167       0.228    0.687 
#> 4 East         3.23           3     0.0992      0.161    0.513 
#> 5 West         4.65           5     0.0107      0.026    0.0828

# Visualize ranking distribution
rank_data %>%
  ggplot(aes(x = rank, y = reorder(region, -theta))) +
  stat_dots(quantiles = 100) +
  scale_x_continuous(breaks = 1:5) +
  labs(
    title = "Ranking Distribution",
    subtitle = "Each dot represents 1% of posterior draws",
    x = "Rank (1 = best, 5 = worst)",
    y = NULL
  )
```

![](tidy_bayesian_workflow_files/figure-html/ranking-1.png)

``` r


# Alternative: bar chart of ranking probabilities
rank_data %>%
  count(region, rank) %>%
  group_by(region) %>%
  mutate(probability = n / sum(n)) %>%
  ggplot(aes(x = rank, y = probability, fill = region)) +
  geom_col() +
  facet_wrap(~region, ncol = 1) +
  scale_x_continuous(breaks = 1:5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Probability of Each Rank by Region",
    x = "Rank (1 = best)",
    y = "Probability"
  ) +
  theme(legend.position = "none")
```

![](tidy_bayesian_workflow_files/figure-html/ranking-2.png)

## Further Reading

- **posterior** package: <https://mc-stan.org/posterior/>
- **bayesplot** package: <https://mc-stan.org/bayesplot/>
- **tidybayes** package: <http://mjskay.github.io/tidybayes/>
- **ggdist** package: <https://mjskay.github.io/ggdist/>
