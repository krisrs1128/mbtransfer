#' Read a CSV file hosted on Figshare
#'
#' Figshare's download redirector (`ndownloader`) can block a bare
#' `read_csv(url)` call because the `figshare.com/ndownloader` redirect path
#' returns a 202 Accepted response instead of a direct 302 → S3 hop.
#' This wrapper normalises both common URL forms to the
#' `ndownloader.figshare.com` subdomain (which redirects reliably to S3),
#' downloads the file with [download.file()], then hands the local copy to
#' [readr::read_csv()].  Within an R session the file is kept in a
#' temporary cache directory so repeated calls for the same URL skip the
#' network entirely.
#'
#' @param url A Figshare `ndownloader` URL in either of the two common forms:
#'   * `"https://ndownloader.figshare.com/files/40275934"` (canonical)
#'   * `"https://figshare.com/ndownloader/files/40275934/subject.csv"`
#'     (also accepted; the hostname and optional trailing filename are
#'     normalised automatically).
#' @param ... Additional arguments forwarded to [readr::read_csv()].
#' @return A [tibble::tibble()] as returned by [readr::read_csv()].
#' @importFrom readr read_csv
#' @export
read_figshare_csv <- function(url, ...) {
  # Normalise figshare.com/ndownloader/files/{id}[/name] to the subdomain form
  # that actually returns a 302 redirect to S3.  The canonical form served by
  # the Figshare API is https://ndownloader.figshare.com/files/{id}.
  url <- sub(
    "^https?://figshare\\.com/ndownloader/files/([0-9]+).*",
    "https://ndownloader.figshare.com/files/\\1",
    url
  )

  cache_dir <- file.path(tempdir(), "mbtransfer_figshare_cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  # Use the numeric file ID (last path component) as the cache filename.
  file_id <- basename(url)
  # Recover a human-readable name from the original URL when available.
  dest <- file.path(cache_dir, file_id)
  if (!file.exists(dest)) {
    message("Downloading file ", file_id, " from Figshare...")
    status <- download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    if (status != 0 || !file.exists(dest)) {
      stop("download.file() failed for ", url, " (status ", status, ")")
    }
  }

  read_csv(dest, ...)
}

#' Build a  `ts`
#' @importFrom dplyr pull filter left_join distinct
#' @examples
#' library(readr)
#' library(tibble)
#' subject <- read_figshare_csv("https://ndownloader.figshare.com/files/40275934")
#' interventions <- read_figshare_csv("https://ndownloader.figshare.com/files/40279171") |>
#'   column_to_rownames("sample")
#' reads <- read_figshare_csv("https://ndownloader.figshare.com/files/40279108") |>
#'   column_to_rownames("sample")
#' samples <- read_figshare_csv("https://ndownloader.figshare.com/files/40275943")
#' ts <- as.matrix(reads) |>
#'   ts_from_dfs(interventions, samples, subject)
#' ts
#' @export
ts_from_dfs <- function(reads, interventions, metadata, subject_data = NULL) {
  subjects <- unique(metadata$subject)
  series <- list()

  # fill in each subject's intervention and abundance data
  for (i in seq_along(subjects)) {
    sample_ix <- metadata |>
      filter(subject == subjects[i]) |>
      pull(time, sample)

    x <- t(as.matrix(reads[names(sample_ix), ]))
    z <- t(as.matrix(interventions[names(sample_ix), ]))
    rownames(z) <- colnames(interventions)

    series[[i]] <- new(
      "ts_inter_single",
      values = x[, order(sample_ix)],
      interventions = z[, order(sample_ix), drop = FALSE],
      time = sort(sample_ix)
    )
  }

  # ensure subject and metadata order agree
  if (!is.null(subject_data)) {
    names(series) <- subjects
    subject_data <- metadata |>
      distinct(subject) |>
      left_join(subject_data)
  }

  ts <- new("ts_inter", series = series, subject_data = subject_data)
  names(ts) <- subjects
  ts
}

#' @export
ts_to_dfs <- function(ts) {
  if (is.null(names(ts))) {
    names(ts) <- seq_along(ts)
  }

  reads <- do.call(cbind, map(ts@series, ~ values(.)))
  interventions <- do.call(cbind, map(ts@series, ~ interventions(.))) |>
    t() |>
    as.data.frame() |>
    rownames_to_column("sample")
  metadata <- map_dfr(
    ts@series,
    ~ tibble(sample = colnames(values(.)), time = .@time),
    .id = "subject"
  )

  list(
    reads = reads,
    interventions = interventions,
    metadata = metadata,
    subject_data = subject_data(ts)
  )
}

#' Joined ts dfs
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom tidyr pivot_longer
#' @export
#' @examples
#' data(sim_ts)
#' pivoted <- pivot_ts(sim_ts)
#' head(pivoted)
#' @export
pivot_ts <- function(ts) {
  dfs <- ts_to_dfs(ts)

  reads <- data.frame(dfs$reads) |>
    rownames_to_column("taxon") |>
    pivot_longer(-taxon, names_to = "sample")

  reads |>
    left_join(dfs$metadata) |>
    left_join(dfs$interventions) |>
    left_join(dfs$subject)
}
