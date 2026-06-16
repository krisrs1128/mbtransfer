# An S4 class representing a collection of time series under environmental shifts

While `ts_inter_single` represents shifts for a single subject,
`ts_inter` includes shifts across a collection of subjects. It includes
an additional slot, `subject_data` to store non-time-varying subject
metadata.

## Slots

- `series`:

  A list of `ts_inter_single` objects, each representing interventions
  and compositional responses for a subject.

- `subject_data`:

  An optional data.frame storing static host-level metadata associated
  with each series.

## Examples

``` r
data(sim_ts)
sim_ts
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
