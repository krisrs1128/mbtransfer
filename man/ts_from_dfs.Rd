% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{ts_from_dfs}
\alias{ts_from_dfs}
\title{Build a  \code{ts}}
\usage{
ts_from_dfs(reads, interventions, metadata, subject_data = NULL)
}
\description{
Build a  \code{ts}
}
\examples{
library(readr)
library(tibble)
subject <- read_csv("https://figshare.com/ndownloader/files/40275934/subject.csv")
interventions <- read_csv("https://figshare.com/ndownloader/files/40279171/interventions.csv") |>
  column_to_rownames("sample")
reads <- read_csv("https://figshare.com/ndownloader/files/40279108/reads.csv") |>
  column_to_rownames("sample")
samples <- read_csv("https://figshare.com/ndownloader/files/40275943/samples.csv")
ts <- as.matrix(reads) |>
  ts_from_dfs(interventions, samples, subject)
ts
}
