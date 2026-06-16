# Convert `patchify_single` output into a `data.frame`

This generates a data.frame of sliding window predictors and responses
by wrapping `patchify_single`. It is in a format that can directly be
used by different regression methods.

## Usage

``` r
patchify_single_df(ts_inter, p, q)
```

## Arguments

- ts_inter:

  An object of class `ts_inter_single` over which to generate sliding
  windows.

- p:

  The number of time lags to. use in the sliding window for the
  microbiome features.

- q:

  The number of time lags to use in the sliding window for
  interventions.

## Examples

``` r
data(sim_ts)
patches <- patchify_single_df(sim_ts[[1]], 2, 2)
head(patches$x)
#>      taxon1_lag2 taxon2_lag2 taxon3_lag2 taxon4_lag2 taxon5_lag2 taxon6_lag2
#> [1,]           5          23          17          30           3         115
#> [2,]           7          43          40          12           2          20
#> [3,]           5          19          42          15           3          23
#> [4,]           3           3          25           7          13           1
#> [5,]           6           7           8          23           7          15
#> [6,]           2          14          13          13           4           2
#>      taxon7_lag2 taxon8_lag2 taxon9_lag2 taxon10_lag2 taxon11_lag2 taxon12_lag2
#> [1,]           8           2          14           42           21            7
#> [2,]          42           3           8            9            5           53
#> [3,]          22          12           4           13            8           47
#> [4,]           6          19          14           15           21           22
#> [5,]          27          13          15           33           11           14
#> [6,]           9          14          22           19           13           14
#>      taxon13_lag2 taxon14_lag2 taxon15_lag2 taxon16_lag2 taxon17_lag2
#> [1,]           12            2           10           12            4
#> [2,]            8            4           18            6            7
#> [3,]            2           23           18            2            7
#> [4,]           21            5           18            8            6
#> [5,]           12            7           13            1           11
#> [6,]            4            6            7            3           11
#>      taxon18_lag2 taxon19_lag2 taxon20_lag2 taxon21_lag2 taxon22_lag2
#> [1,]           20            5            1            3            0
#> [2,]            5            2           27            3            2
#> [3,]            9            9           21            5           48
#> [4,]            4           10            0            4           13
#> [5,]            7            9            5            8            1
#> [6,]            5           24            7            1            6
#>      taxon23_lag2 taxon24_lag2 taxon25_lag2 taxon26_lag2 taxon27_lag2
#> [1,]           22            1            5           35           10
#> [2,]            1            7            8            5            2
#> [3,]            3           29            2           12            8
#> [4,]           21            1            7            7            6
#> [5,]            3            4            4            4           22
#> [6,]           12            6           21            8            5
#>      taxon28_lag2 taxon29_lag2 taxon30_lag2 taxon31_lag2 taxon32_lag2
#> [1,]            7           29            3           14            5
#> [2,]           13            4            5            1            0
#> [3,]            6            4           14            7            1
#> [4,]            5            0            5            4            8
#> [5,]           12            0            6           14            0
#> [6,]           10           10            6            2            2
#>      taxon33_lag2 taxon34_lag2 taxon35_lag2 taxon36_lag2 taxon37_lag2
#> [1,]            6            3            7           10           48
#> [2,]            5            2           13            1           14
#> [3,]            2            8            5            9            8
#> [4,]           10           10            1           11            4
#> [5,]            2           13           14           17           14
#> [6,]            7           15            4            5           21
#>      taxon38_lag2 taxon39_lag2 taxon40_lag2 taxon41_lag2 taxon42_lag2
#> [1,]            6            8           10            6           19
#> [2,]            9            4           11            8           56
#> [3,]            7           13            4            4            5
#> [4,]            7            1            7            8            7
#> [5,]            4            2            7           10           16
#> [6,]            4            7            6            2            8
#>      taxon43_lag2 taxon44_lag2 taxon45_lag2 taxon46_lag2 taxon47_lag2
#> [1,]           18           13           44           23            2
#> [2,]           30           22           13           79            1
#> [3,]           26           36           10           11            2
#> [4,]            2            9            7            6            7
#> [5,]           10           32           21            6            5
#> [6,]            3            0            8           30            9
#>      taxon48_lag2 taxon49_lag2 taxon50_lag2 taxon51_lag2 taxon52_lag2
#> [1,]            4            8           17            0            4
#> [2,]            7            9            6            0            2
#> [3,]           23            3           16           12            1
#> [4,]           24           17           19           11            3
#> [5,]            8           12           19            4           19
#> [6,]           14           17           22            5            2
#>      taxon53_lag2 taxon54_lag2 taxon55_lag2 taxon56_lag2 taxon57_lag2
#> [1,]           15           78           17            0            1
#> [2,]           12           10            4            3           25
#> [3,]           19           19            8            3           14
#> [4,]           10           32            8            5           17
#> [5,]            5           26           16           17            8
#> [6,]           13           10            3            9            6
#>      taxon58_lag2 taxon59_lag2 taxon60_lag2 taxon61_lag2 taxon62_lag2
#> [1,]           12            5           14           14           21
#> [2,]            8            3           21           56           10
#> [3,]            6            1            5           40           12
#> [4,]            8           14            9           17           12
#> [5,]            3           17           26            9            7
#> [6,]            8            9           10           52            4
#>      taxon63_lag2 taxon64_lag2 taxon65_lag2 taxon66_lag2 taxon67_lag2
#> [1,]           26           28           13           16            2
#> [2,]            2           15            7           29            1
#> [3,]           25           21           10           32            9
#> [4,]            2           13            9           11           21
#> [5,]            4           25           11            6            4
#> [6,]            6           29            4            8            6
#>      taxon68_lag2 taxon69_lag2 taxon70_lag2 taxon71_lag2 taxon72_lag2
#> [1,]           12           13           29            6           13
#> [2,]            0            8           40            5            1
#> [3,]           31            3           31            7            1
#> [4,]            9            3           10            4            5
#> [5,]           11            2            7           31           13
#> [6,]           15           26           13            6            4
#>      taxon73_lag2 taxon74_lag2 taxon75_lag2 taxon76_lag2 taxon77_lag2
#> [1,]            5           17            6            4            7
#> [2,]            7           22            1            0            2
#> [3,]           13           12            6            5           18
#> [4,]           16            2            3            4           16
#> [5,]           18           16            5           21            2
#> [6,]           17            9            8            0            5
#>      taxon78_lag2 taxon79_lag2 taxon80_lag2 taxon81_lag2 taxon82_lag2
#> [1,]            2            6           14            5           16
#> [2,]            6            0            4            7            4
#> [3,]           13            3            7           19            3
#> [4,]            3            8            0           12            5
#> [5,]           22           14           11            9           11
#> [6,]            7           11           15            4            5
#>      taxon83_lag2 taxon84_lag2 taxon85_lag2 taxon86_lag2 taxon87_lag2
#> [1,]            5            5           10           18            7
#> [2,]           12            3            6            6           15
#> [3,]           36           10            7           24            1
#> [4,]           16           12            2            6           11
#> [5,]            9           10            0           12            6
#> [6,]           13           28            0           18            9
#>      taxon88_lag2 taxon89_lag2 taxon90_lag2 taxon91_lag2 taxon92_lag2
#> [1,]           37            6           26           21            7
#> [2,]           10            7           12           18           13
#> [3,]           13           14            5            3            9
#> [4,]            5            3            1           13            8
#> [5,]           12            8           17            4           25
#> [6,]           11            4            6            6           20
#>      taxon93_lag2 taxon94_lag2 taxon95_lag2 taxon96_lag2 taxon97_lag2
#> [1,]           10           10            3           58           10
#> [2,]            4           23           32           23           41
#> [3,]            0            9           13           15            8
#> [4,]           31            6           31           20           16
#> [5,]           14            8            9            4           15
#> [6,]            4            6           31            9           13
#>      taxon98_lag2 taxon99_lag2 taxon100_lag2 taxon1_lag1 taxon2_lag1
#> [1,]           10           12             5           7          43
#> [2,]           17            4             0           5          19
#> [3,]           14            9             9           3           3
#> [4,]            0           10             1           6           7
#> [5,]           12            7            16           2          14
#> [6,]            1           18             8           7           2
#>      taxon3_lag1 taxon4_lag1 taxon5_lag1 taxon6_lag1 taxon7_lag1 taxon8_lag1
#> [1,]          40          12           2          20          42           3
#> [2,]          42          15           3          23          22          12
#> [3,]          25           7          13           1           6          19
#> [4,]           8          23           7          15          27          13
#> [5,]          13          13           4           2           9          14
#> [6,]          29          13           6          11          14          13
#>      taxon9_lag1 taxon10_lag1 taxon11_lag1 taxon12_lag1 taxon13_lag1
#> [1,]           8            9            5           53            8
#> [2,]           4           13            8           47            2
#> [3,]          14           15           21           22           21
#> [4,]          15           33           11           14           12
#> [5,]          22           19           13           14            4
#> [6,]          11            2           15            1            2
#>      taxon14_lag1 taxon15_lag1 taxon16_lag1 taxon17_lag1 taxon18_lag1
#> [1,]            4           18            6            7            5
#> [2,]           23           18            2            7            9
#> [3,]            5           18            8            6            4
#> [4,]            7           13            1           11            7
#> [5,]            6            7            3           11            5
#> [6,]           14            2           19            3            3
#>      taxon19_lag1 taxon20_lag1 taxon21_lag1 taxon22_lag1 taxon23_lag1
#> [1,]            2           27            3            2            1
#> [2,]            9           21            5           48            3
#> [3,]           10            0            4           13           21
#> [4,]            9            5            8            1            3
#> [5,]           24            7            1            6           12
#> [6,]            7           23           29           36           13
#>      taxon24_lag1 taxon25_lag1 taxon26_lag1 taxon27_lag1 taxon28_lag1
#> [1,]            7            8            5            2           13
#> [2,]           29            2           12            8            6
#> [3,]            1            7            7            6            5
#> [4,]            4            4            4           22           12
#> [5,]            6           21            8            5           10
#> [6,]            5           21           14           14            4
#>      taxon29_lag1 taxon30_lag1 taxon31_lag1 taxon32_lag1 taxon33_lag1
#> [1,]            4            5            1            0            5
#> [2,]            4           14            7            1            2
#> [3,]            0            5            4            8           10
#> [4,]            0            6           14            0            2
#> [5,]           10            6            2            2            7
#> [6,]           15            8            0           15            8
#>      taxon34_lag1 taxon35_lag1 taxon36_lag1 taxon37_lag1 taxon38_lag1
#> [1,]            2           13            1           14            9
#> [2,]            8            5            9            8            7
#> [3,]           10            1           11            4            7
#> [4,]           13           14           17           14            4
#> [5,]           15            4            5           21            4
#> [6,]            9            3            7            8            9
#>      taxon39_lag1 taxon40_lag1 taxon41_lag1 taxon42_lag1 taxon43_lag1
#> [1,]            4           11            8           56           30
#> [2,]           13            4            4            5           26
#> [3,]            1            7            8            7            2
#> [4,]            2            7           10           16           10
#> [5,]            7            6            2            8            3
#> [6,]            4           10            3            9           16
#>      taxon44_lag1 taxon45_lag1 taxon46_lag1 taxon47_lag1 taxon48_lag1
#> [1,]           22           13           79            1            7
#> [2,]           36           10           11            2           23
#> [3,]            9            7            6            7           24
#> [4,]           32           21            6            5            8
#> [5,]            0            8           30            9           14
#> [6,]           12            7            3            4            6
#>      taxon49_lag1 taxon50_lag1 taxon51_lag1 taxon52_lag1 taxon53_lag1
#> [1,]            9            6            0            2           12
#> [2,]            3           16           12            1           19
#> [3,]           17           19           11            3           10
#> [4,]           12           19            4           19            5
#> [5,]           17           22            5            2           13
#> [6,]           18            4            6            2            7
#>      taxon54_lag1 taxon55_lag1 taxon56_lag1 taxon57_lag1 taxon58_lag1
#> [1,]           10            4            3           25            8
#> [2,]           19            8            3           14            6
#> [3,]           32            8            5           17            8
#> [4,]           26           16           17            8            3
#> [5,]           10            3            9            6            8
#> [6,]           12            9            3            0            3
#>      taxon59_lag1 taxon60_lag1 taxon61_lag1 taxon62_lag1 taxon63_lag1
#> [1,]            3           21           56           10            2
#> [2,]            1            5           40           12           25
#> [3,]           14            9           17           12            2
#> [4,]           17           26            9            7            4
#> [5,]            9           10           52            4            6
#> [6,]            8           13            0            5            6
#>      taxon64_lag1 taxon65_lag1 taxon66_lag1 taxon67_lag1 taxon68_lag1
#> [1,]           15            7           29            1            0
#> [2,]           21           10           32            9           31
#> [3,]           13            9           11           21            9
#> [4,]           25           11            6            4           11
#> [5,]           29            4            8            6           15
#> [6,]            8           10            3            6            9
#>      taxon69_lag1 taxon70_lag1 taxon71_lag1 taxon72_lag1 taxon73_lag1
#> [1,]            8           40            5            1            7
#> [2,]            3           31            7            1           13
#> [3,]            3           10            4            5           16
#> [4,]            2            7           31           13           18
#> [5,]           26           13            6            4           17
#> [6,]           14            7            4            7            3
#>      taxon74_lag1 taxon75_lag1 taxon76_lag1 taxon77_lag1 taxon78_lag1
#> [1,]           22            1            0            2            6
#> [2,]           12            6            5           18           13
#> [3,]            2            3            4           16            3
#> [4,]           16            5           21            2           22
#> [5,]            9            8            0            5            7
#> [6,]           30           25           29            4            6
#>      taxon79_lag1 taxon80_lag1 taxon81_lag1 taxon82_lag1 taxon83_lag1
#> [1,]            0            4            7            4           12
#> [2,]            3            7           19            3           36
#> [3,]            8            0           12            5           16
#> [4,]           14           11            9           11            9
#> [5,]           11           15            4            5           13
#> [6,]           14           14           12            6           28
#>      taxon84_lag1 taxon85_lag1 taxon86_lag1 taxon87_lag1 taxon88_lag1
#> [1,]            3            6            6           15           10
#> [2,]           10            7           24            1           13
#> [3,]           12            2            6           11            5
#> [4,]           10            0           12            6           12
#> [5,]           28            0           18            9           11
#> [6,]            6            7            7           16           13
#>      taxon89_lag1 taxon90_lag1 taxon91_lag1 taxon92_lag1 taxon93_lag1
#> [1,]            7           12           18           13            4
#> [2,]           14            5            3            9            0
#> [3,]            3            1           13            8           31
#> [4,]            8           17            4           25           14
#> [5,]            4            6            6           20            4
#> [6,]           26           12            2           10           15
#>      taxon94_lag1 taxon95_lag1 taxon96_lag1 taxon97_lag1 taxon98_lag1
#> [1,]           23           32           23           41           17
#> [2,]            9           13           15            8           14
#> [3,]            6           31           20           16            0
#> [4,]            8            9            4           15           12
#> [5,]            6           31            9           13            1
#> [6,]            3           15           17            4            0
#>      taxon99_lag1 taxon100_lag1 intervention1_lag1 intervention1_lag0
#> [1,]            4             0                  0                  0
#> [2,]            9             9                  0                  0
#> [3,]           10             1                  0                  0
#> [4,]            7            16                  0                  0
#> [5,]           18             8                  0                  0
#> [6,]           11            11                  0                  0
head(patches$y)
#>      taxon1 taxon2 taxon3 taxon4 taxon5 taxon6 taxon7 taxon8 taxon9 taxon10
#> [1,]      3      3     25      7     13      1      6     19     14      15
#> [2,]      6      7      8     23      7     15     27     13     15      33
#> [3,]      2     14     13     13      4      2      9     14     22      19
#> [4,]      7      2     29     13      6     11     14     13     11       2
#> [5,]     11     18     10      7     12     22     16      9      5      19
#> [6,]      6      6      5      8     15     35      2     11     24      24
#>      taxon11 taxon12 taxon13 taxon14 taxon15 taxon16 taxon17 taxon18 taxon19
#> [1,]      21      22      21       5      18       8       6       4      10
#> [2,]      11      14      12       7      13       1      11       7       9
#> [3,]      13      14       4       6       7       3      11       5      24
#> [4,]      15       1       2      14       2      19       3       3       7
#> [5,]       3       8       7       4       9       9       6      32      12
#> [6,]       4       9       9       4      13      10      10       1       1
#>      taxon20 taxon21 taxon22 taxon23 taxon24 taxon25 taxon26 taxon27 taxon28
#> [1,]       0       4      13      21       1       7       7       6       5
#> [2,]       5       8       1       3       4       4       4      22      12
#> [3,]       7       1       6      12       6      21       8       5      10
#> [4,]      23      29      36      13       5      21      14      14       4
#> [5,]       4      12      13       3       3       3       6       6       1
#> [6,]       9      10      17      21       7      10      14       5      24
#>      taxon29 taxon30 taxon31 taxon32 taxon33 taxon34 taxon35 taxon36 taxon37
#> [1,]       0       5       4       8      10      10       1      11       4
#> [2,]       0       6      14       0       2      13      14      17      14
#> [3,]      10       6       2       2       7      15       4       5      21
#> [4,]      15       8       0      15       8       9       3       7       8
#> [5,]       3       5      12       5      13      36      20       8       7
#> [6,]       3       7      13      32       0      11      14      10      21
#>      taxon38 taxon39 taxon40 taxon41 taxon42 taxon43 taxon44 taxon45 taxon46
#> [1,]       7       1       7       8       7       2       9       7       6
#> [2,]       4       2       7      10      16      10      32      21       6
#> [3,]       4       7       6       2       8       3       0       8      30
#> [4,]       9       4      10       3       9      16      12       7       3
#> [5,]       7       8      19       7      10       9      10       6       4
#> [6,]       3       3       9      14       9       5      16       5       1
#>      taxon47 taxon48 taxon49 taxon50 taxon51 taxon52 taxon53 taxon54 taxon55
#> [1,]       7      24      17      19      11       3      10      32       8
#> [2,]       5       8      12      19       4      19       5      26      16
#> [3,]       9      14      17      22       5       2      13      10       3
#> [4,]       4       6      18       4       6       2       7      12       9
#> [5,]      17      24       5      11      11      15      22       3      19
#> [6,]      10      23      14       1       0      15       8       3       7
#>      taxon56 taxon57 taxon58 taxon59 taxon60 taxon61 taxon62 taxon63 taxon64
#> [1,]       5      17       8      14       9      17      12       2      13
#> [2,]      17       8       3      17      26       9       7       4      25
#> [3,]       9       6       8       9      10      52       4       6      29
#> [4,]       3       0       3       8      13       0       5       6       8
#> [5,]       6       7       5       6      12      13      18       6      12
#> [6,]       4       6       2       5      14      45       4       2       9
#>      taxon65 taxon66 taxon67 taxon68 taxon69 taxon70 taxon71 taxon72 taxon73
#> [1,]       9      11      21       9       3      10       4       5      16
#> [2,]      11       6       4      11       2       7      31      13      18
#> [3,]       4       8       6      15      26      13       6       4      17
#> [4,]      10       3       6       9      14       7       4       7       3
#> [5,]       4       1      19      11      10       5      12      13      12
#> [6,]      26       9       3       6       7      14       8       4       2
#>      taxon74 taxon75 taxon76 taxon77 taxon78 taxon79 taxon80 taxon81 taxon82
#> [1,]       2       3       4      16       3       8       0      12       5
#> [2,]      16       5      21       2      22      14      11       9      11
#> [3,]       9       8       0       5       7      11      15       4       5
#> [4,]      30      25      29       4       6      14      14      12       6
#> [5,]       7      16       9      18      17       3      20       5       9
#> [6,]       7       5       5       3       8      13       4       1      13
#>      taxon83 taxon84 taxon85 taxon86 taxon87 taxon88 taxon89 taxon90 taxon91
#> [1,]      16      12       2       6      11       5       3       1      13
#> [2,]       9      10       0      12       6      12       8      17       4
#> [3,]      13      28       0      18       9      11       4       6       6
#> [4,]      28       6       7       7      16      13      26      12       2
#> [5,]       2       6      19       6       3      41      17       9       7
#> [6,]       4      20      10      25       6       3       9       1      11
#>      taxon92 taxon93 taxon94 taxon95 taxon96 taxon97 taxon98 taxon99 taxon100
#> [1,]       8      31       6      31      20      16       0      10        1
#> [2,]      25      14       8       9       4      15      12       7       16
#> [3,]      20       4       6      31       9      13       1      18        8
#> [4,]      10      15       3      15      17       4       0      11       11
#> [5,]       5      16      16      21       3       5      30      17       13
#> [6,]      14      12       9       9      14       8       9       7       15
```
