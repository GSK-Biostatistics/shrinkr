# Spike-and-slab prior for testing homogeneity

Creates a mixture of two Normal distributions. Since tau is a scale
parameter, this must be wrapped in
[`distributional::dist_truncated()`](https://pkg.mitchelloharawild.com/distributional/reference/dist_truncated.html)
with `lower = 0` before passing to [`shrink()`](shrink.md):

## Usage

``` r
prior_spike_slab(
  spike_location = 0,
  spike_scale = 0.01,
  slab_scale = 1,
  spike_prob = 0.5
)
```

## Arguments

- spike_location:

  Location of the spike (default 0)

- spike_scale:

  Scale of the spike component (default 0.01)

- slab_scale:

  Scale of the slab component (default 1)

- spike_prob:

  Probability of the spike component (default 0.5)

## Value

A spike-and-slab mixture distribution

## Details

    tau_prior <- dist_truncated(prior_spike_slab(), lower = 0)

## Examples

``` r
tau_prior <- distributional::dist_truncated(
  prior_spike_slab(spike_prob = 0.5, spike_scale = 0.01, slab_scale = 1),
  lower = 0
)
```
