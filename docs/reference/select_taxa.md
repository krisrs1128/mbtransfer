# Significant Taxa using Mirror Statistics

This selects taxa through a data splitting procedure. It retrains an
mbtransfer model across random splits of the data. Partial dependence
profiles for pairs of splits are compared with one another – if they
agree for a given feature, then that feature is considered more likely
to have a true effect on the response. This function only supports
inference of the effects of interventions on taxa responses, but the
same principle could apply to estimate significant relationships between
taxa.

## Usage

``` r
select_taxa(ts, w0, w1, tr_fun, qvalue = 0.2, ...)
```

## Arguments

- ts:

  An object of class `ts_inter` containing the time-varying microbiome
  community, environmental interventions, and static host features data.
  The columns for each element of the `values` matrix are expected to be
  consecutive timepoints in an changing community.
  [`ts_from_dfs()`](ts_from_dfs.md) is a helper to create these objects
  from `data.frames` with the relevant information.

- w0:

  One of the counterfactuals with which to compute partial dependence
  profiles. See `steps` or `pulses` for helpers in generating these
  counterfactuals. The procedure concatenates these counterfactuals to
  the end of the series and computes the difference in the forecasts.

- w1:

  One of the counterfactuals with which to compute partial dependence
  profiles. See `steps` or `pulses` for helpers in generating these
  counterfactuals. The procedure concatenates these counterfactuals to
  the end of the series and computes the difference in the forecasts.

- tr_fun:

  A function that can be used to train over random splits. In the
  examples in this package, we use an anonymous function that fills our
  chosen hyperparameters in `mbtransfer`. For example \\x) mbtransfer(x,
  P = 1, Q = 1) will fit a lag-1 transfer function model on all the
  random splits.

- qvalue:

  The target False Discovery Rate. Defaults to 0.2.

## Examples

``` r
data(sim_ts)
w0 <- cbind(sim_ts[[1]]@interventions, matrix(0, nrow = 1, ncol = 3))
w1 <- cbind(sim_ts[[1]]@interventions, matrix(1, nrow = 1, ncol = 3))
S <- select_taxa(sim_ts, w0, w1, \(x) mbtransfer(x, 1, 1, nrounds = 10), n_splits = 1)
#> Training models for split 1/1
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
```
