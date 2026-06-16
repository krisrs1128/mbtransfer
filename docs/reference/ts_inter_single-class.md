# An S4 class to represent a single subject's series

Taxonomic abundances, binary/continuous interventions, and underlying
timepoints for a single subject are stored in this unified data
structure. A `ts_inter` object is list of these `ts_inter_single`
classes (together wtih static subject data).

## Slots

- `values`:

  A matrix whose rows are taxa and whose columns are samples.

- `time`:

  A vector of the timepoints associated with the samples in `values`.

- `interventions`:

  A matrix whose rows are perturbations, columns are samples, and values
  are either binary interventions or continuous input series,
  representing the value of the exogenous influence.

## Examples

``` r
data(sim_ts)
sim_ts[[1]]
#> # A ts_inter_single object | 1 taxa | 27 timepoints:
#> taxa:
#>      S1_T1 S1_T2 S1_T3 S1_T4  
#> tax1     5     7     5     3 …
#> tax2    23    43    19     3 …
#> tax3    17    40    42    25 …
#> tax4    30    12    15     7 …
#>          ⋮     ⋮     ⋮     ⋮ ⋱
#> interventions:
#>    S1_T1 S1_T2 S1_T3 S1_T4  
#> P1     0     0     0     0 …
```
