
#' @export
consistency_mirror <- function(effects) {
  sgn <- sign(effects[, 1] * effects[, 2])
  magnitude <- rowMeans(abs(effects))
  sgn * magnitude
}

#' @importFrom dplyr mutate bind_rows row_number
#' @importFrom tibble tibble
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
#' @export
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

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
select_taxa <- function(ts, w0, w1, tr_fun, qvalue = 0.2, ...) {
  effects <- pd_splits(ts, w0, w1, tr_fun, ...) # immediate effect
  ms <- consistency_mirror_multisplit(effects)
  
  taxa <- ms |>
    select(multisplit, m, lag) %>%
    split(.$lag) %>%
    map(~ split(., .$multisplit) %>% map(~ pull(., m))) |>
    map(~ which(multiple_data_splitting(., q = qvalue))) |>
    map(~ taxa(ts)[.])
  
  list(taxa = taxa, ms = ms, effects = effects)
}
