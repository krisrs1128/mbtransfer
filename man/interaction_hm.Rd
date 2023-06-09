% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis.R
\name{interaction_hm}
\alias{interaction_hm}
\title{Visualize Time Series Heatmaps}
\usage{
interaction_hm(values_df, taxa, condition = NULL, r = 0, ...)
}
\arguments{
\item{values_df}{A data.frame containing joined information from a \code{ts_inter}
object. Each row should be one count from one taxon/subject pair. The
output of \code{pivot_ts} reshape data to this format.}

\item{taxa}{A character vector of taxa to filter down to in the final
visualization.}

\item{condition}{A subject-level variable to use during faceting. Defaults to
no faceting.}

\item{r}{When ordering subjects, we consider averages after rounding
timepoints into bins. Large \code{r} corresponds to fine-grained binning,
negative \code{r} are more coarse bins. Defaults to 0. See \code{subject_order} for
details.}
}
\description{
This is a helper function to visualize time series for a subset of taxa
across all subjects. It is used in the raw data figures in the case studies
section of the accompanying manuscript.
}
\examples{
library(dplyr)
data(sim_ts)
pivoted <- sim_ts |>
  pivot_ts() |>
  mutate(v_pos = V1 > 0)
interaction_hm(pivoted, c("tax1", "tax2"), "v_pos")

pivoted |>
  group_by(taxon) |>
  mutate(value = rank(value) / n()) |>
  interaction_hm(c("tax1", "tax2"), "v_pos")
}
