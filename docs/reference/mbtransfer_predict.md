# Prediction method for `mbtransfer` objects

Prediction method for `mbtransfer` objects

## Usage

``` r
mbtransfer_predict(object, newdata)
```

## Arguments

- object:

  An object of class `mbtransfer_model`, as generated using the
  `mbtransfer` function. This includes trained boosting models for every
  taxon, stored within the `@parameters` slot.

- newdata:

  A new `ts_inter` object over which to perform prediction. This method
  will make predictions for every timepoint that appears in the
  `@interventions` slot but not the `@values`.

## Examples

``` r
data(sim_ts)
fit <- mbtransfer(sim_ts)
ts_subset <- subset_values(sim_ts, 1:25)
predict(fit, ts_subset)
#> # A ts_inter object | 100 taxa | 50 subjects | 27.00 ± 0.00 timepoints:
#> 
#> S1:
#>      S1_T1 S1_T2 S1_T3 S1_T4  
#> tax1     5     7     5     3 …
#> tax2    23    43    19     3 …
#> tax3    17    40    42    25 …
#> tax4    30    12    15     7 …
#>          ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> S2:
#>      S2_T1 S2_T2 S2_T3 S2_T4  
#> tax1     6     4    11     8 …
#> tax2    53    12    15     1 …
#> tax3     6     8    23     5 …
#> tax4    24    11    12    15 …
#>          ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> S3:
#>      S3_T1 S3_T2 S3_T3 S3_T4  
#> tax1     5     3     6     5 …
#> tax2     9    13    26     3 …
#> tax3    60    13     4     6 …
#> tax4    22    17     6     8 …
#>          ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> and 47 more subjects.
```
