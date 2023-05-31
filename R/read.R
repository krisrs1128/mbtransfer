
#'Build a  `ts` 
#' @importFrom dplyr pull filter left_join distinct
#' @examples
#' library(tibble)
#' subject <- read_csv("https://figshare.com/ndownloader/files/40275934/subject.csv")
#' interventions <- read_csv("https://figshare.com/ndownloader/files/40279171/interventions.csv") |>
#' column_to_rownames("sample")
#' reads <- read_csv("https://figshare.com/ndownloader/files/40279108/reads.csv") |>
#'  column_to_rownames("sample")
#' samples <- read_csv("https://figshare.com/ndownloader/files/40275943/samples.csv")
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
      interventions = z[, order(sample_ix), drop=FALSE],
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
