# Generate Counterfactual versions of a ts_inter object

mbtransfer makes predictions starting from ts_inter objects, so in order
to simulate counterfactuals we need to provide alternative versions of
those objects as input. This function truncates an existing ts object at
start_ix and creates new intervention series according to the values of
w0 and w1 (starting from the truncation point).

## Usage

``` r
counterfactual_ts(ts, w0, w1, start_ix = NULL)
```

## Arguments

- ts:

  The starting ts_inter object from which to build the '
  counterfactuals. By default, this is the last timepoint in the current
  series.

- w0:

  The first version of the intervention series to consider.

- w1:

  The second version of the intervention series to consider.

- start_ix:

  The truncation position for the original ts. Defaults to no
  truncation, which appends new interventions to the end of the existing
  series.

## Value

Two ts_inter objects with interventions corresponding to w0 and w1.

## Examples

``` r
data(sim_ts)
ws <- steps(c("P1" = TRUE), 2, 3, 4)
ts_star <- counterfactual_ts(sim_ts, ws[[1]], ws[[2]], start_ix = 10)
```
