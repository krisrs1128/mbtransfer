# Prediction for a Single Subject

This loops over all timepoints for a single subject and makes
predictions by comparing the number of timepoints in @interventions and
in @values. The gap will be filled in one step at a time using
`mbtransfer_predict_step()`.

## Usage

``` r
model_predict_single(fit, ts_inter, lags, subject = NULL, interactions = NULL)
```

## Arguments

- fit:

  An object of class `mbtransfer_model`, as generated using the
  `mbtransfer` function. This includes trained boosting models for every
  taxon, stored within the `@parameters` slot.

- ts_inter:

  A new `ts_inter_single` object over which to perform prediction. This
  method will make predictions for every timepoint that appears in the
  `@interventions` slot but not the `@values`. This is assumed to be a
  single subject from a larger `ts_inter` object.

- lags:

  A vector specifying `P` and `Q` in the trained mbtransfer model.

- subject:

  A static data frame of subject-level variables. This will be
  concatenated to time-varying intervention and taxonomic covariates
  when making predictions. This is analogous to the training process.
