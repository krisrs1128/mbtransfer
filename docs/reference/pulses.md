# Hypothetical Pulse Interventions

Create counterfactual perturbations with specific 0/1 intervention and
total length patterns.

## Usage

``` r
pulses(p_states, lags = 1, L = 3, w_star = c(0, 1))
```

## Arguments

- lags:

  The time length of any interventions that we generate.

- L:

  The total length returned, including both intervention and
  non-intervention timepoints

- w_star:

  The unique values to include in and out of the intervention. Defaults
  to 1/0, respectively.

- p_stats:

  A list specifying the row names and whether they should include
  interventions. Any name set to TRUE will include intervention effects,
  while those with FALSE will not.

## Examples

``` r
library(mbtransfer)
pulses(c("P1" = TRUE), 1, 4)
#> [[1]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    0    0    0
#> 
#> [[2]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    0    0    1
#> 
pulses(c("P1" = TRUE), 1:3, 4)
#> [[1]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    0    0    0
#> 
#> [[2]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    0    0    0
#> 
#> [[3]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    0    0    0
#> 
#> [[4]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    0    0    1
#> 
#> [[5]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    0    1    0
#> 
#> [[6]]
#>    Tn_1 Tn_2 Tn_3 Tn_4
#> P1    0    1    0    0
#> 
```
