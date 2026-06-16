# Subset values of a ts object

This is a helper to filter the number of timepoints across all `values`
slots in a `ts_inter` class.

## Usage

``` r
subset_values(ts, ix)
```

## Examples

``` r
data(sim_ts)
subset_values(sim_ts, 1:2)
#> # A ts_inter object | 100 taxa | 50 subjects | 2.00 ± 0.00 timepoints:
#> 
#> S1:
#>      S1_T1 S1_T2  
#> tax1     5     7 …
#> tax2    23    43 …
#>          ⋮     ⋮ ⋱
#> 
#> S2:
#>      S2_T1 S2_T2  
#> tax1     6     4 …
#> tax2    53    12 …
#>          ⋮     ⋮ ⋱
#> 
#> S3:
#>      S3_T1 S3_T2  
#> tax1     5     3 …
#> tax2     9    13 …
#>          ⋮     ⋮ ⋱
#> 
#> and 47 more subjects.
subset_values(sim_ts, 4:8)
#> # A ts_inter object | 100 taxa | 50 subjects | 5.00 ± 0.00 timepoints:
#> 
#> S1:
#>      S1_T4 S1_T5 S1_T6 S1_T7  
#> tax1     3     6     2     7 …
#> tax2     3     7    14     2 …
#> tax3    25     8    13    29 …
#> tax4     7    23    13    13 …
#>          ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> S2:
#>      S2_T4 S2_T5 S2_T6 S2_T7  
#> tax1     8     4     6     8 …
#> tax2     1     1     5    10 …
#> tax3     5    13    12     9 …
#> tax4    15    14    28     8 …
#>          ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> S3:
#>      S3_T4 S3_T5 S3_T6 S3_T7  
#> tax1     5    11    10     6 …
#> tax2     3     3     7     5 …
#> tax3     6    17     2     9 …
#> tax4     8    18     3    16 …
#>          ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> and 47 more subjects.
```
