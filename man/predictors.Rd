% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patches.R
\name{predictors}
\alias{predictors}
\title{Generate a Predictor Matrix}
\usage{
predictors(ts_inter, lags, subject, interactions = NULL)
}
\arguments{
\item{ts_inter}{An object of class \code{ts_inter_single} over which to generate
sliding windows.}

\item{lags}{A vector specifying \code{P} and \code{Q} in the trained mbtransfer model.}

\item{subject}{A static data frame of subject-level variables. This will be
concatenated to time-varying intervention and taxonomic covariates when
making predictions. This is analogous to the training process.}
}
\description{
\code{patchify_df} is used to generate training data. At test time, we don't need
all the patches, but we do need to construct regressors for the immediate
past. This will be matched to the model's original training data and used to
generate predictions.
}
\examples{
data(sim_ts)
ts <- subset_values(sim_ts, 1:5)
predictors(ts[[1]], c(2, 2), NULL)
}
