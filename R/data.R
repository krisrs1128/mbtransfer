#' Simulated Time Series
#'
#' These are data generated according to a negative binomial vector
#' autoregressive model. The specific generative mechanism is described in
#' Section 3.1 of "mbtransfer: Microbiome Intervention Analysis using Transfer
#' Functions and Mirror Statistics" by Sankaran and Jeganathan.
#'
#' @format A `ts_inter` class of 50 subjets with 30 timepoints each.
#'
#' @source This run is `sim_input_004.rda` from this collection,
#'   https://drive.google.com/drive/folders/19mDoKbT2qvO9AlRdn_wQKBbqNh8Cx22F as
#'   reshaped and interpolated by this script:
#'   https://github.com/krisrs1128/microbiome_interventions/blob/main/scripts/forecasting_metrics.Rmd
"sim_ts"
