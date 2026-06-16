# Transfer Function Model for `ts_inter` objects

This is the main prediction function in the `mbtransfer` package. Given
an object of class `ts_inter` (see [`ts_from_dfs()`](ts_from_dfs.md)),
this will fit a collection of linear gradient boosting-based transfer
function models. The output is an object of class `mbtransfer_model`.
Each component boosting model is contained in the `@parameters` slot,
which is a list whose j^th element is the model for the j^th taxon (row)
within each `ts`'s values matrix.

## Usage

``` r
mbtransfer(
  ts_inter,
  P = 1,
  Q = 1,
  nrounds = 500,
  early_stopping_rounds = 5,
  verbose = 0,
  lambda = 0.1,
  alpha = 0.01,
  eta = 0.1,
  interactions = "search",
  nthreads = 0,
  ...
)
```

## Arguments

- ts_inter:

  An object of class `ts_inter` containing the time-varying microbiome
  community, environmental interventions, and static host features data.
  The columns for each element of the `values` matrix are expected to be
  consecutive timepoints in an changing community.
  [`ts_from_dfs()`](ts_from_dfs.md) is a helper to create these objects
  from `data.frames` with the relevant information.

- P:

  The number of historical community composition timepoints to consider
  when making predictions.

- Q:

  The number of historical intervention timepoints to consider when
  making predictions.

- nrounds:

  The maximum number of rounds for each taxon's gradient boosting model.
  Smaller values will lead to faster training, but at the risk of poorer
  fits. Defaults to 500.

- early_stopping_rounds:

  If the loss function does not improve after this many rounds, then the
  model is assumed to have converged and training is stopped. Defaults
  to 5.

- verbose:

  Should information about each gradient boosting model's performance be
  printed? Allowable values are 2 (all information), 1 (some
  information), and 0 (no information, default).

- lambda:

  The l2-regularization value in the linear gradient boosting model.
  Defaults to 1e-1.

- alpha:

  The l1-regularization value in the linear gradient boosting model.
  Defaults to 1e-2. This value generally leads to less sparse fits,
  which creates useful variation in potential downstream mirror
  statistics calculations.

- eta:

  The learning rate. Defaults to 0.1. This is slower than the default in
  xgboost (0.3) but has been found to improve stability when needing to
  train on taxa with a wide range of abundances.

## Examples

``` r
data(sim_ts)
fit <- mbtransfer(sim_ts)
fit@parameters[[1]]
#> ##### xgb.Booster
#> call:
#>   xgb.train(params = params, data = dmat, nrounds = nrounds, evals = list(train = dmat), 
#>     verbose = verbose, early_stopping_rounds = early_stopping_rounds)
#> # of features: 195 
#> # of rounds:  500 
#> xgb.attributes:
#>    best_iteration, best_score 
#> callbacks:
#>    early_stop, evaluation_log 
#> evaluation_log:
#>   iter train_rmse
#>  <num>      <num>
#>      1   4.702258
#>      2   4.658111
#>    ---        ---
#>    499   4.290243
#>    500   4.290227
```
