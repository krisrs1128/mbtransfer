# Joined ts dfs

Joined ts dfs

## Usage

``` r
pivot_ts(ts)
```

## Examples

``` r
data(sim_ts)
pivoted <- pivot_ts(sim_ts)
#> Joining with `by = join_by(sample)`
#> Joining with `by = join_by(sample)`
#> Joining with `by = join_by(subject)`
head(pivoted)
#> # A tibble: 6 × 7
#>   taxon sample value subject  time    P1    V1
#>   <chr> <glue> <dbl> <chr>   <dbl> <dbl> <dbl>
#> 1 tax1  S1_T1      5 S1        -11     0 0.330
#> 2 tax1  S1_T2      7 S1        -10     0 0.330
#> 3 tax1  S1_T3      5 S1         -9     0 0.330
#> 4 tax1  S1_T4      3 S1         -8     0 0.330
#> 5 tax1  S1_T5      6 S1         -7     0 0.330
#> 6 tax1  S1_T6      2 S1         -6     0 0.330
```
