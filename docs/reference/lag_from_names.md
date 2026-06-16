# Detect the time lags based on column names

This is a helper function for deciding on P and Q in transfer function
models using just the names of an example set of covariates.

## Usage

``` r
lag_from_names(names, group = "taxon")
```

## Examples

``` r
lag_from_names(c("taxon1_lag1", "taxon1_lag2", "taxon1_lag3"))
#> [1] 3
```
