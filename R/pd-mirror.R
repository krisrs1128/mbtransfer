
#' Mirror Statistics from a Single Split
#' 
#' This is based on the mirror statistic in Dai et al. (2022). The idea is that,
#' if a feature is null, the sign of the effect is as likely to be positive or
#' negative (this symmetry supports FDR estimation). If there is a real effect,
#' then the signs are more likely to agree (sign == 1 below) and the magnitude
#' should be large.
#' @param effects A list of arrays containing estimated partial dependence
#'   effects. The list indexes different splits. Within each list element, the
#'   expected dimensions are n_taxa x time_lag x random_split_index.
#' @export
#' @examples
#' effects <- matrix(rnorm(500), 250, 2)
#' m <- consistency_mirror(effects)
#' hist(m, 20)
#' 
#' # long tail on the right is the real effect
#' effects[1:5, ] <- runif(10, 2, 4)
#' m <- consistency_mirror(effects)
#' hist(m, 20)
consistency_mirror <- function(effects) {
  sgn <- sign(effects[, 1] * effects[, 2])
  magnitude <- rowMeans(abs(effects))
  sgn * magnitude
}

#' Compute Mirrors across Splits
#' The mirror of `consistency_mirror` works on effects derived from a single
#' split. In practice, we will want to have mirror statistics across a multiple
#' splits. This is a small wrapper of that function that computes mirrors for
#' effects available along a list.
#' @param effects A list of arrays containing estimated partial dependence
#'   effects. The list indexes different splits. Within each list element, the
#'   expected dimensions are n_taxa x time_lag x random_split_index.
#' @importFrom dplyr mutate bind_rows row_number
#' @importFrom tibble tibble
#' @examples
#' effects <- list()
#' effects[[1]] <- array(rnorm(1000), dim = c(250, 2, 2))
#' effects[[2]] <- array(rnorm(1000), dim = c(250, 2, 2))
#' m <- consistency_mirror_multisplit(effects)
#' @export
consistency_mirror_multisplit <- function(effects) {
  ms <- list()
  k <- 1
  for (s in seq_along(effects)) {
    for (lag in seq_len(dim(effects[[s]])[3])) {
      ms[[k]] <- tibble(
        m = consistency_mirror(effects[[s]][,, lag]),
        lag = lag,
        multisplit = s
      ) |>
        mutate(taxon = row_number())
      k <- k + 1
    }
  }

  bind_rows(ms)
}

#' Mirror Splits from Partial Dependence
#'
#' For ech split, we return an array of dimension n_taxa x 2 fits x n_lags. Each
#' entry contains the estimated counterfactual effect for that taxon and lag
#' combination across the two models (each fit on a different random sample of
#' data).
#'
#' @param ts An object of class `ts_inter` containing the time-varying
#'   microbiome community, environmental interventions, and static host features
#'   data. The columns for each element of the `values` matrix are expected to
#'   be consecutive timepoints in an changing community. `ts_from_dfs()` is a
#'   helper to create these objects from `data.frames` with the relevant
#'   information.
#' @param w0 One of the counterfactuals with which to compute partial dependence
#'   profiles. See `steps` or `pulses` for helpers in generating these
#'   counterfactuals. The procedure concatenates these counterfactuals to the
#'   end of the series and computes the difference in the forecasts.
#' @param w1 One of the counterfactuals with which to compute partial dependence
#'   profiles. See `steps` or `pulses` for helpers in generating these
#'   counterfactuals. The procedure concatenates these counterfactuals to the
#'   end of the series and computes the difference in the forecasts.
#' @param tr_fun A function that can be used to train over random splits. In the
#'   examples in this package, we use an anonymous function that fills our
#'   chosen hyperparameters in `mbtransfer`. For example \(x) mbtransfer(x, P =
#'   1, Q = 1) will fit a lag-1 transfer function model on all the random
#'   splits.
#' @export
#' @examples
#' data(sim_ts)
#' w0 <- cbind(sim_ts[[1]]@interventions, matrix(0, nrow = 1, ncol = 3))
#' w1 <- cbind(sim_ts[[1]]@interventions, matrix(1, nrow = 1, ncol = 3))
#' pd_splits(sim_ts, w0, w1, \(x) mbtransfer(x, 1, 1, nrounds = 10), n_splits = 1)
pd_splits <- function(ts, w0, w1, tr_fun, n_splits = 20, ...) {
  effects <- replicate(n_splits, array(dim = c(nrow(ts[[1]]), 2, ncol(w0))), simplify = FALSE)

  for (s in seq_len(n_splits)) {
    print(glue("Training models for split {s}/{n_splits}"))
    split_ix <- sample(length(ts), 0.5 * length(ts))
    ts_split <- list(ts[split_ix], ts[-split_ix])
    fits <- map(ts_split, tr_fun)

    for (i in seq_along(ts_split)) {
      effects[[s]][, i,] <- pd_effects(fits[[i]], ts_split[[i]], w0, w1, ...)
    }
  }

  effects
}

#' @importFrom purrr map
#' @importFrom abind abind
#' @export
pd_summary <- function(y0, y1, ix, summary_fun = mean) {
  map((y0 - y1)[, ix], ~ values(.)) |>
    abind(along = 3) |>
    apply(1:2, summary_fun)
}

#' Counterfactual Partial Dependence Effects
#' 
#' Our selection algorithm, `select_taxa` depends on estimates for
#' counterfactual partial dependence effects across all taxa and time lags. This
#' supports estimation of these effects for a single random data split pair.
#' @param fit An object of class `mbtransfer_model`, as generated using the
#'   `mbtransfer` function. This includes trained boosting models for every
#'   taxon, stored within the `@parameters` slot.
#' @param ts An object of class `ts_inter` containing the time-varying
#'   microbiome community, environmental interventions, and static host features
#'   data. The columns for each element of the `values` matrix are expected to
#'   be consecutive timepoints in an changing community. `ts_from_dfs()` is a
#'   helper to create these objects from `data.frames` with the relevant
#'   information.
#' @param w0 One of the counterfactuals with which to compute partial dependence
#'   profiles. See `steps` or `pulses` for helpers in generating these
#'   counterfactuals. The procedure concatenates these counterfactuals to the
#'   end of the series and computes the difference in the forecasts.
#' @param w1 One of the counterfactuals with which to compute partial dependence
#'   profiles. See `steps` or `pulses` for helpers in generating these
#'   counterfactuals. The procedure concatenates these counterfactuals to the
#'   end of the series and computes the difference in the forecasts.
#' @examples
#' data(sim_ts)
#' w0 <- cbind(sim_ts[[1]]@interventions, matrix(0, nrow = 1, ncol = 3))
#' w1 <- cbind(sim_ts[[1]]@interventions, matrix(1, nrow = 1, ncol = 3))
#' fit <- mbtransfer(sim_ts, 1, 1, nrounds = 10)
#' pd_effects(fit, sim_ts, w0, w1)
pd_effects <- function(fit, ts, w0, w1, n_sample = NULL, patch_len = 8, intervention_len = NULL) {
  if (is.null(n_sample)) {
    n_sample <- 0.5 * length(ts) * ncol(ts[[1]]) / patch_len
  }

  # sampled patches under two counterfactual interventions
  ts_star <- sample_ts(ts, n_sample, patch_len, intervention_len = ncol(w0)) |>
    counterfactual_ts(w0, w1)

  # compute difference in predictions across future ix
  y_hats <- map(ts_star, ~ predict(fit, .))
  ix <- seq(patch_len + 1, patch_len + ncol(w0))
  pd_summary(y_hats[[1]], y_hats[[2]], ix)
}

#' Significant Taxa using Mirror Statistics
#' 
#' This selects taxa through a data splitting procedure. It retrains an
#' mbtransfer model across random splits of the data. Partial dependence
#' profiles for pairs of splits are compared with one another -- if they agree
#' for a given feature, then that feature is considered more likely to have a
#' true effect on the response. This function only supports inference of the
#' effects of interventions on taxa responses, but the same principle could
#' apply to estimate significant relationships between taxa.
#' 
#' @param ts An object of class `ts_inter` containing the time-varying
#'   microbiome community, environmental interventions, and static host features
#'   data. The columns for each element of the `values` matrix are expected to
#'   be consecutive timepoints in an changing community. `ts_from_dfs()` is a
#'   helper to create these objects from `data.frames` with the relevant
#'   information.
#' @param w0 One of the counterfactuals with which to compute partial dependence
#'   profiles. See `steps` or `pulses` for helpers in generating these
#'   counterfactuals. The procedure concatenates these counterfactuals to the
#'   end of the series and computes the difference in the forecasts.
#' @param w1 One of the counterfactuals with which to compute partial dependence
#'   profiles. See `steps` or `pulses` for helpers in generating these
#'   counterfactuals. The procedure concatenates these counterfactuals to the
#'   end of the series and computes the difference in the forecasts.
#' @param tr_fun A function that can be used to train over random splits. In the
#'   examples in this package, we use an anonymous function that fills our
#'   chosen hyperparameters in `mbtransfer`. For example \(x) mbtransfer(x, P =
#'   1, Q = 1) will fit a lag-1 transfer function model on all the random
#'   splits.
#' @param qvalue The target False Discovery Rate. Defaults to 0.2.
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @examples
#' data(sim_ts)
#' w0 <- cbind(sim_ts[[1]]@interventions, matrix(0, nrow = 1, ncol = 3))
#' w1 <- cbind(sim_ts[[1]]@interventions, matrix(1, nrow = 1, ncol = 3))
#' select_taxa(sim_ts, w0, w1, \(x) mbtransfer(x, 1, 1, nrounds = 10), n_splits = 1)
#' @export
select_taxa <- function(ts, w0, w1, tr_fun, qvalue = 0.2, ...) {
  effects <- pd_splits(ts, w0, w1, tr_fun, ...)
  ms <- consistency_mirror_multisplit(effects)
  
  taxa <- ms |>
    select(multisplit, m, lag) %>%
    split(.$lag) %>%
    map(~ split(., .$multisplit) %>% map(~ pull(., m))) |>
    map(~ which(multiple_data_splitting(., q = qvalue))) |>
    map(~ taxa(ts)[.])
  
  list(taxa = taxa, ms = ms, effects = effects)
}
