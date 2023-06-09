% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patches.R
\name{predictor_names}
\alias{predictor_names}
\title{Create uniform column names for patchified df}
\usage{
predictor_names(x_dim, w_dim)
}
\description{
Create uniform column names for patchified df
}
\examples{
data(sim_ts)
predictor_names(c(2, 3), c(2, 2))
predictor_names(c(5, 2), c(1, 1))
}
