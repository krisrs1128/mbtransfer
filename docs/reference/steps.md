# Hypothetical Step Interventions

Hypothetical Step Interventions

## Usage

``` r
steps(p_states, starts = 1, lengths = 1:3, L = 3, w_star = c(0, 1))
```

## Examples

``` r
steps(c("P1" = TRUE), 1:3, 2:3, 4)
#> [[1]]
#>    t1 t2 t3 t4
#> P1  0  0  0  0
#> 
#> [[2]]
#>    t1 t2 t3 t4
#> P1  1  1  0  0
#> 
#> [[3]]
#>    t1 t2 t3 t4
#> P1  1  1  1  0
#> 
#> [[4]]
#>    t1 t2 t3 t4
#> P1  0  1  1  0
#> 
#> [[5]]
#>    t1 t2 t3 t4
#> P1  0  1  1  1
#> 
#> [[6]]
#>    t1 t2 t3 t4
#> P1  0  0  1  1
#> 
```
