# Generate a Predictor Matrix

`patchify_df` is used to generate training data. At test time, we don't
need all the patches, but we do need to construct regressors for the
immediate past. This will be matched to the model's original training
data and used to generate predictions.

## Usage

``` r
predictors(ts_inter, lags, subject, interactions = NULL)
```

## Arguments

- ts_inter:

  An object of class `ts_inter_single` over which to generate sliding
  windows.

- lags:

  A vector specifying `P` and `Q` in the trained mbtransfer model.

- subject:

  A static data frame of subject-level variables. This will be
  concatenated to time-varying intervention and taxonomic covariates
  when making predictions. This is analogous to the training process.

## Examples

``` r
data(sim_ts)
ts <- subset_values(sim_ts, 1:5)
predictors(ts[[1]], c(2, 2), NULL)
#>      taxon1_lag2 taxon2_lag2 taxon3_lag2 taxon4_lag2 taxon5_lag2 taxon6_lag2
#> [1,]           3           3          25           7          13           1
#>      taxon7_lag2 taxon8_lag2 taxon9_lag2 taxon10_lag2 taxon11_lag2 taxon12_lag2
#> [1,]           6          19          14           15           21           22
#>      taxon13_lag2 taxon14_lag2 taxon15_lag2 taxon16_lag2 taxon17_lag2
#> [1,]           21            5           18            8            6
#>      taxon18_lag2 taxon19_lag2 taxon20_lag2 taxon21_lag2 taxon22_lag2
#> [1,]            4           10            0            4           13
#>      taxon23_lag2 taxon24_lag2 taxon25_lag2 taxon26_lag2 taxon27_lag2
#> [1,]           21            1            7            7            6
#>      taxon28_lag2 taxon29_lag2 taxon30_lag2 taxon31_lag2 taxon32_lag2
#> [1,]            5            0            5            4            8
#>      taxon33_lag2 taxon34_lag2 taxon35_lag2 taxon36_lag2 taxon37_lag2
#> [1,]           10           10            1           11            4
#>      taxon38_lag2 taxon39_lag2 taxon40_lag2 taxon41_lag2 taxon42_lag2
#> [1,]            7            1            7            8            7
#>      taxon43_lag2 taxon44_lag2 taxon45_lag2 taxon46_lag2 taxon47_lag2
#> [1,]            2            9            7            6            7
#>      taxon48_lag2 taxon49_lag2 taxon50_lag2 taxon51_lag2 taxon52_lag2
#> [1,]           24           17           19           11            3
#>      taxon53_lag2 taxon54_lag2 taxon55_lag2 taxon56_lag2 taxon57_lag2
#> [1,]           10           32            8            5           17
#>      taxon58_lag2 taxon59_lag2 taxon60_lag2 taxon61_lag2 taxon62_lag2
#> [1,]            8           14            9           17           12
#>      taxon63_lag2 taxon64_lag2 taxon65_lag2 taxon66_lag2 taxon67_lag2
#> [1,]            2           13            9           11           21
#>      taxon68_lag2 taxon69_lag2 taxon70_lag2 taxon71_lag2 taxon72_lag2
#> [1,]            9            3           10            4            5
#>      taxon73_lag2 taxon74_lag2 taxon75_lag2 taxon76_lag2 taxon77_lag2
#> [1,]           16            2            3            4           16
#>      taxon78_lag2 taxon79_lag2 taxon80_lag2 taxon81_lag2 taxon82_lag2
#> [1,]            3            8            0           12            5
#>      taxon83_lag2 taxon84_lag2 taxon85_lag2 taxon86_lag2 taxon87_lag2
#> [1,]           16           12            2            6           11
#>      taxon88_lag2 taxon89_lag2 taxon90_lag2 taxon91_lag2 taxon92_lag2
#> [1,]            5            3            1           13            8
#>      taxon93_lag2 taxon94_lag2 taxon95_lag2 taxon96_lag2 taxon97_lag2
#> [1,]           31            6           31           20           16
#>      taxon98_lag2 taxon99_lag2 taxon100_lag2 taxon1_lag1 taxon2_lag1
#> [1,]            0           10             1           6           7
#>      taxon3_lag1 taxon4_lag1 taxon5_lag1 taxon6_lag1 taxon7_lag1 taxon8_lag1
#> [1,]           8          23           7          15          27          13
#>      taxon9_lag1 taxon10_lag1 taxon11_lag1 taxon12_lag1 taxon13_lag1
#> [1,]          15           33           11           14           12
#>      taxon14_lag1 taxon15_lag1 taxon16_lag1 taxon17_lag1 taxon18_lag1
#> [1,]            7           13            1           11            7
#>      taxon19_lag1 taxon20_lag1 taxon21_lag1 taxon22_lag1 taxon23_lag1
#> [1,]            9            5            8            1            3
#>      taxon24_lag1 taxon25_lag1 taxon26_lag1 taxon27_lag1 taxon28_lag1
#> [1,]            4            4            4           22           12
#>      taxon29_lag1 taxon30_lag1 taxon31_lag1 taxon32_lag1 taxon33_lag1
#> [1,]            0            6           14            0            2
#>      taxon34_lag1 taxon35_lag1 taxon36_lag1 taxon37_lag1 taxon38_lag1
#> [1,]           13           14           17           14            4
#>      taxon39_lag1 taxon40_lag1 taxon41_lag1 taxon42_lag1 taxon43_lag1
#> [1,]            2            7           10           16           10
#>      taxon44_lag1 taxon45_lag1 taxon46_lag1 taxon47_lag1 taxon48_lag1
#> [1,]           32           21            6            5            8
#>      taxon49_lag1 taxon50_lag1 taxon51_lag1 taxon52_lag1 taxon53_lag1
#> [1,]           12           19            4           19            5
#>      taxon54_lag1 taxon55_lag1 taxon56_lag1 taxon57_lag1 taxon58_lag1
#> [1,]           26           16           17            8            3
#>      taxon59_lag1 taxon60_lag1 taxon61_lag1 taxon62_lag1 taxon63_lag1
#> [1,]           17           26            9            7            4
#>      taxon64_lag1 taxon65_lag1 taxon66_lag1 taxon67_lag1 taxon68_lag1
#> [1,]           25           11            6            4           11
#>      taxon69_lag1 taxon70_lag1 taxon71_lag1 taxon72_lag1 taxon73_lag1
#> [1,]            2            7           31           13           18
#>      taxon74_lag1 taxon75_lag1 taxon76_lag1 taxon77_lag1 taxon78_lag1
#> [1,]           16            5           21            2           22
#>      taxon79_lag1 taxon80_lag1 taxon81_lag1 taxon82_lag1 taxon83_lag1
#> [1,]           14           11            9           11            9
#>      taxon84_lag1 taxon85_lag1 taxon86_lag1 taxon87_lag1 taxon88_lag1
#> [1,]           10            0           12            6           12
#>      taxon89_lag1 taxon90_lag1 taxon91_lag1 taxon92_lag1 taxon93_lag1
#> [1,]            8           17            4           25           14
#>      taxon94_lag1 taxon95_lag1 taxon96_lag1 taxon97_lag1 taxon98_lag1
#> [1,]            8            9            4           15           12
#>      taxon99_lag1 taxon100_lag1 intervention1_lag1 intervention1_lag0
#> [1,]            7            16                  0                  0
```
