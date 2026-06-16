# Randomly subset a `ts_inter` object

There are times we want to randomly sample patches from across
timepoints and subjects. This is a helper function that does this random
sampling adapted to `ts_inter` objects, subsetting taxa and intervention
timepoints simultaneously.

## Usage

``` r
sample_ts(ts, n, patch_len = 5, intervention_len = 0)
```

## Arguments

- ts:

  A `ts_inter` object whose random patches we want to sample.

- n:

  The number of randomly sampled time series to return.

- patch_len:

  The length of the randomly sampled patches

- intervention_len:

  We will pad the @time slot by an extra intervention_len set of
  timepoints. Defaults to 0.

## Examples

``` r
data(sim_ts)
sample_ts(sim_ts, 10, 4, 1)
#> # A ts_inter object | 100 taxa | 10 subjects | 4.00 ± 0.00 timepoints:
#> 
#> S21:
#>      S21_T22 S21_T23 S21_T24 S21_T25  
#> tax1      10       8       7       6 …
#> tax2      13      17       1      19 …
#> tax3      13      14       7       2 …
#> tax4       6      14      10      12 …
#>            ⋮       ⋮       ⋮       ⋮ ⋱
#> 
#> S28:
#>      S28_T4 S28_T5 S28_T6 S28_T7  
#> tax1      2      1      3     12 …
#> tax2     10      6     14      0 …
#> tax3     26      2     23     16 …
#> tax4     28      8      4      4 …
#>           ⋮      ⋮      ⋮      ⋮ ⋱
#> 
#> S18:
#>      S18_T18 S18_T19 S18_T20 S18_T21  
#> tax1       9       8       4       1 …
#> tax2       8       6      14      44 …
#> tax3       2      12      16       2 …
#> tax4       5      11       7      41 …
#>            ⋮       ⋮       ⋮       ⋮ ⋱
#> 
#> and 7 more subjects.
```
