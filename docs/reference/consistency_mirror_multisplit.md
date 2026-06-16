# Compute Mirrors across Splits

The mirror of `consistency_mirror` works on effects derived from a
single split. In practice, we will want to have mirror statistics across
a multiple splits. This is a small wrapper of that function that
computes mirrors for effects available along a list.

## Usage

``` r
consistency_mirror_multisplit(effects)
```

## Arguments

- effects:

  A list of arrays containing estimated partial dependence effects. The
  list indexes different splits. Within each list element, the expected
  dimensions are n_taxa x time_lag x random_split_index.

## Examples

``` r
effects <- list()
effects[[1]] <- array(rnorm(1000), dim = c(250, 2, 2))
effects[[2]] <- array(rnorm(1000), dim = c(250, 2, 2))
m <- consistency_mirror_multisplit(effects)
```
