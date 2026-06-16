# Predictions for a single timepoint and subject

This make predictions for all taxa for a single timestep ahead in one
subject. It loops over the trained boosting models for each taxon
predicts a single value for each.

## Usage

``` r
model_predict_step(ts_inter, fit, lags, subject = NULL, interactions = NULL)
```

## Arguments

- ts_inter:

  A new `ts_inter_single` object over which to perform prediction. This
  method will make predictions for every timepoint that appears in the
  `@interventions` slot but not the `@values`. This is assumed to be a
  single subject from a larger `ts_inter` object.

- fit:

  An object of class `mbtransfer_model`, as generated using the
  `mbtransfer` function. This includes trained boosting models for every
  taxon, stored within the `@parameters` slot.

- lags:

  A vector specifying `P` and `Q` in the trained mbtransfer model.

- subject:

  A static data frame of subject-level variables. This will be
  concatenated to time-varying intervention and taxonomic covariates
  when making predictions. This is analogous to the training process.
