# Add Interactions to Predictor Matrix

Add Interactions to Predictor Matrix

## Usage

``` r
append_interactions(x, interactions = NULL)
```

## Arguments

- x:

  The original matrix containing the source of the interaction terms.

- interactions:

  An n_interactions x 2 matrix containing interactions to add to the
  original matrix x.

## Value

The matrix x with new columns associated wtih each interaction.
