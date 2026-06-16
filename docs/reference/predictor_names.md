# Create uniform column names for patchified df

Create uniform column names for patchified df

## Usage

``` r
predictor_names(x_dim, w_dim)
```

## Examples

``` r
data(sim_ts)
predictor_names(c(2, 3), c(2, 2))
#>  [1] "taxon1_lag3"        "taxon2_lag3"        "taxon1_lag2"       
#>  [4] "taxon2_lag2"        "taxon1_lag1"        "taxon2_lag1"       
#>  [7] "intervention1_lag1" "intervention2_lag1" "intervention1_lag0"
#> [10] "intervention2_lag0"
predictor_names(c(5, 2), c(1, 1))
#>  [1] "taxon1_lag2"        "taxon2_lag2"        "taxon3_lag2"       
#>  [4] "taxon4_lag2"        "taxon5_lag2"        "taxon1_lag1"       
#>  [7] "taxon2_lag1"        "taxon3_lag1"        "taxon4_lag1"       
#> [10] "taxon5_lag1"        "intervention1_lag0"
```
