# Detect the time lags over which a model was trained

This is a helper for determining the P, Q associated with a fitted
mbtransfer model.

## Usage

``` r
time_lags(fit)
```

## Examples

``` r
data(sim_ts)
fit <- mbtransfer(sim_ts)
time_lags(fit@parameters[[1]])
#> [1] 1 1
```
