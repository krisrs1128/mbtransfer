% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{model_predict_step}
\alias{model_predict_step}
\title{Predictions for a single timepoint and subject}
\usage{
model_predict_step(ts_inter, fit, lags, subject = NULL, interactions = NULL)
}
\arguments{
\item{ts_inter}{A new \code{ts_inter_single} object over which to perform
prediction. This method will make predictions for every timepoint that
appears in the \verb{@interventions} slot but not the \verb{@values}. This is assumed
to be a single subject from a larger \code{ts_inter} object.}

\item{fit}{An object of class \code{mbtransfer_model}, as generated using the
\code{mbtransfer} function. This includes trained boosting models for every
taxon, stored within the \verb{@parameters} slot.}

\item{lags}{A vector specifying \code{P} and \code{Q} in the trained mbtransfer model.}

\item{subject}{A static data frame of subject-level variables. This will be
concatenated to time-varying intervention and taxonomic covariates when
making predictions. This is analogous to the training process.}
}
\description{
This make predictions for all taxa for a single timestep ahead in one
subject. It loops over the trained boosting models for each taxon predicts a
single value for each.
}
