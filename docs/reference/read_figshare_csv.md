# Read a CSV file hosted on Figshare

Figshare's download redirector (`ndownloader`) can block a bare
`read_csv(url)` call because the `figshare.com/ndownloader` redirect
path returns a 202 Accepted response instead of a direct 302 → S3 hop.
This wrapper normalises both common URL forms to the
`ndownloader.figshare.com` subdomain (which redirects reliably to S3),
downloads the file with
[`download.file()`](https://rdrr.io/r/utils/download.file.html), then
hands the local copy to
[`readr::read_csv()`](https://readr.tidyverse.org/reference/read_delim.html).
Within an R session the file is kept in a temporary cache directory so
repeated calls for the same URL skip the network entirely.

## Usage

``` r
read_figshare_csv(url, ...)
```

## Arguments

- url:

  A Figshare `ndownloader` URL in either of the two common forms:

  - `"https://ndownloader.figshare.com/files/40275934"` (canonical)

  - `"https://figshare.com/ndownloader/files/40275934/subject.csv"`
    (also accepted; the hostname and optional trailing filename are
    normalised automatically).

- ...:

  Additional arguments forwarded to
  [`readr::read_csv()`](https://readr.tidyverse.org/reference/read_delim.html).

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
as returned by
[`readr::read_csv()`](https://readr.tidyverse.org/reference/read_delim.html).
