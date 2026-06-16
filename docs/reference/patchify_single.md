# Sliding Windows for a ts_inter_single object

This creates sliding windows of intervention and community composition
features for a single subject. It returns a list giving the "patchified"
data, which can then be organized into matrices for prediction.

## Usage

``` r
patchify_single(ts_inter, p = 2, q = 3)
```

## Arguments

- ts_inter:

  An object of class `ts_inter_single` over which to generate sliding
  windows.

- p:

  The number of time lags to. use in the sliding window for the
  microbiome features.

- q:

  The number of time lags to use in the sliding window for
  interventions.

## Examples

``` r
data(sim_ts)
patches <- patchify_single(sim_ts[[1]])
head(patches$x[[1]])
#>      S1_T1 S1_T2
#> tax1     5     7
#> tax2    23    43
#> tax3    17    40
#> tax4    30    12
#> tax5     3     2
#> tax6   115    20
head(patches$w[[1]])
#>    S1_T1 S1_T2 S1_T3
#> P1     0     0     0
head(patches$y[[1]])
#>      S1_T4
#> tax1     3
#> tax2     3
#> tax3    25
#> tax4     7
#> tax5    13
#> tax6     1
```
