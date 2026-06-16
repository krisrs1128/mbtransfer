# Prediction method for mbtransfer models

Prediction method for mbtransfer models

## Usage

``` r
# S4 method for class 'mbtransfer_model'
predict(object, newdata)
```

## Examples

``` r
data(sim_ts)
fit <- mbtransfer(sim_ts, 1, 1)
sim_sub <- subset_values(sim_ts, 1:26) # remove last timepoint
y_hat <- predict(fit, sim_sub)
plot(values(sim_ts[[1]])[, 27], values(y_hat[[1]])[, 27], xlab = "Truth", ylab = "Predicted")
```
