# Build a `ts`

Build a `ts`

## Usage

``` r
ts_from_dfs(reads, interventions, metadata, subject_data = NULL)
```

## Examples

``` r
library(readr)
library(tibble)
subject <- read_figshare_csv("https://ndownloader.figshare.com/files/40275934")
#> Downloading file 40275934 from Figshare...
#> Rows: 20 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): subject, diet
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
interventions <- read_figshare_csv("https://ndownloader.figshare.com/files/40279171") |>
  column_to_rownames("sample")
#> Downloading file 40279171 from Figshare...
#> Rows: 236 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): sample
#> dbl (2): Plant, Animal
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
reads <- read_figshare_csv("https://ndownloader.figshare.com/files/40279108") |>
  column_to_rownames("sample")
#> Downloading file 40279108 from Figshare...
#> Rows: 236 Columns: 192
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr   (1): sample
#> dbl (191): Otu000001, Otu000002, Otu000003, Otu000004, Otu000005, Otu000006,...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
samples <- read_figshare_csv("https://ndownloader.figshare.com/files/40275943")
#> Downloading file 40275943 from Figshare...
#> Rows: 236 Columns: 4
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (3): sample, subject, condition
#> dbl (1): time
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
ts <- as.matrix(reads) |>
  ts_from_dfs(interventions, samples, subject)
#> Joining with `by = join_by(subject)`
ts
#> # A ts_inter object | 191 taxa | 20 subjects | 11.80 ± 3.57 timepoints:
#> 
#> Plant5:
#>             DD53 DD111  DD24  DD29  
#> Otu000001 10.218 9.149 8.392 9.399 …
#> Otu000002      0     0  5.15 8.454 …
#> Otu000003  6.413     0 3.509  3.57 …
#> Otu000004   7.46 6.997 6.855  8.04 …
#>                ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> Plant7:
#>            DD79  DD34  DD72 DD144  
#> Otu000001 9.311 9.368 8.052 9.734 …
#> Otu000002     0     0 4.838 1.844 …
#> Otu000003 7.201 7.575 5.674 7.127 …
#> Otu000004 9.375 9.591 8.313 9.435 …
#>               ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> Plant4:
#>            DD31 DD146  DD38 DD138  
#> Otu000001 8.176 8.512 8.284 8.466 …
#> Otu000002 4.545     0     0     0 …
#> Otu000003 8.471 8.226 7.681 6.527 …
#> Otu000004     0     0     0     0 …
#>               ⋮     ⋮     ⋮     ⋮ ⋱
#> 
#> and 17 more subjects.
```
