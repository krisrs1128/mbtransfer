
#' Transfer Function Model for `ts_inter` objects
#' 
#' This is the main prediction function in the `mbtransfer` package. Given an
#' object of class `ts_inter` (see `ts_from_dfs()`), this will fit a collection
#' of linear gradient boosting-based transfer function models. The output is an
#' object of class `mbtransfer_model`. Each component boosting model is
#' contained in the `@parameters` slot, which is a list whose j^th element is
#' the model for the j^th taxon (row) within each `ts`'s values matrix.
#' 
#' @param ts_inter An object of class `ts_inter` containing the time-varying
#'   microbiome community, environmental interventions, and static host features
#'   data. The columns for each element of the `values` matrix are expected to
#'   be consecutive timepoints in an changing community. `ts_from_dfs()` is a
#'   helper to create these objects from `data.frames` with the relevant
#'   information.
#' @param P The number of historical community composition timepoints to
#'   consider when making predictions.
#' @param Q The number of historical intervention timepoints to consider when
#'   making predictions.
#' @param nrounds The maximum number of rounds for each taxon's gradient
#'   boosting model. Smaller values will lead to faster training, but at the
#'   risk of poorer fits. Defaults to 500.
#' @param early_stopping_rounds If the loss function does not improve after this
#'   many rounds, then the model is assumed to have converged and training is
#'   stopped. Defaults to 5.
#' @param verbose Should information about each gradient boosting model's
#'   performance be printed? Allowable values are 2 (all information), 1 (some
#'   information), and 0 (no information, default).
#' @param lambda The l2-regularization value in the linear gradient boosting
#'   model. Defaults to 1e-2.
#' @param alpha The l1-regularization value in the linear gradient boosting
#'   model. Defaults to 1e-2. This value generally leads to less sparse fits,
#'   which creates useful variation in potential downstream mirror statistics
#'   calculations.
#' @param eta The learning rate. Defaults to 0.05. This is slower than the
#'   default in xgboost (0.3) but has been found to improve stability when
#'   needing to train on taxa with a wide range of abundances.
#' @importFrom progress progress_bar
#' @importFrom glue glue
#' @importFrom xgboost xgboost
#' @export
mbtransfer <- function(ts_inter, P = 1, Q = 1, nrounds = 500,
                       early_stopping_rounds = 5, verbose = 0, lambda = 1e-2,
                       alpha = 1e-2, eta = 0.05, ...) {
  train_data <- patchify_df(ts_inter, P, Q)
  fit <- list()

  pb <- progress_bar$new(total = length(train_data$y), format = "[:bar] :percent ETA: :eta")
  for (j in seq_along(train_data$y)) {
    pb$tick()
    fit[[j]] <- xgboost(
      data = train_data$x, label = train_data$y[[j]], nrounds = nrounds,
      booster = "gblinear", lambda = lambda, alpha = alpha,
      early_stopping_rounds = early_stopping_rounds, verbose = verbose, ...
    )
  }

  hyper <- list(P = P, Q = Q, nrounds = nrounds, ...)
  new("mbtransfer_model", parameters = fit, method = "mbtransfer", hyper = hyper)
}

#' @export
mbtransfer_predict <- function(object, newdata) {
  lags <- time_lags(object@parameters[[1]])
  result <- list()
  subject <- subject_data(newdata) |>
    select(-subject) |>
    as.matrix()

  for (i in seq_along(newdata)) {
    result[[i]] <- model_predict_single(
      object@parameters,
      newdata[[i]],
      lags,
      subject[i,, drop = FALSE]
    )
  }

  names(result) <- names(newdata)
  new("ts_inter", series = result, subject_data = subject_data(newdata))
}

model_predict_single <- function(fit, ts_inter, lags, subject = NULL) {
  n_time <- ncol(ts_inter)
  w <- interventions(ts_inter)
  while(ncol(ts_inter) < ncol(w)) {
    ts_inter <- model_predict_step(ts_inter, fit, lags, subject)
  }

  colnames(values(ts_inter)) <- colnames(w)
  ts_inter
}

model_predict_step <- function(ts_inter, fit, lags, subject = NULL) {
  xz <- predictors(ts_inter, lags, subject)
  y_hat <- vector(length = nrow(ts_inter))
  for (j in seq_len(nrow(ts_inter))) {
    y_hat[j] <- predict(fit[[j]], xz)
  }

  values(ts_inter) <- cbind(values(ts_inter), y_hat)
  delta <- median(diff(ts_inter@time))
  ts_inter
}

#' @export
setClass(
  "mbtransfer_model",
  slots = c(
    parameters = "ANY",
    method = "character",
    hyper = "list"
  )
)

#' @export
setMethod("predict",  c(object = "mbtransfer_model"), mbtransfer_predict)
