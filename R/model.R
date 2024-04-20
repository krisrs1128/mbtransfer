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
#'   model. Defaults to 1e-1.
#' @param alpha The l1-regularization value in the linear gradient boosting
#'   model. Defaults to 1e-2. This value generally leads to less sparse fits,
#'   which creates useful variation in potential downstream mirror statistics
#'   calculations.
#' @param eta The learning rate. Defaults to 0.1. This is slower than the
#'   default in xgboost (0.3) but has been found to improve stability when
#'   needing to train on taxa with a wide range of abundances.
#' @importFrom progress progress_bar
#' @importFrom glue glue
#' @importFrom xgboost xgboost
#' @export
#' @examples
#' data(sim_ts)
#' fit <- mbtransfer(sim_ts)
#' fit@parameters[[1]]
mbtransfer <- function(ts_inter, P = 1, Q = 1, nrounds = 500,
                       early_stopping_rounds = 5, verbose = 0, lambda = 1e-1,
                       alpha = 1e-2, eta = 0.1, interactions = "search", nthread = -1, ...) {
  train_data <- patchify_df(ts_inter, P, Q, interactions)
  if (!is.null(train_data$interactions)) {
    train_data$x <- append_interactions(train_data$x, train_data$interactions)
  }

  fit <- list()
  pb <- progress_bar$new(total = length(train_data$y), format = "[:bar] :percent ETA: :eta")
  for (j in seq_along(train_data$y)) {
    pb$tick()
    fit[[j]] <- xgboost(
      data = train_data$x, label = train_data$y[[j]], nrounds = nrounds,
      booster = "gblinear", alpha = alpha, lambda = lambda,
      eta = eta, nthread = nthread,
      early_stopping_rounds = early_stopping_rounds, verbose = verbose, ...
    )
  }

  hyper <- list(P = P, Q = Q, nrounds = nrounds, ...)
  new("mbtransfer_model", parameters = fit, method = "mbtransfer", hyper = hyper, interactions = train_data$interactions)
}

#' Prediction method for `mbtransfer` objects
#'
#' @param object An object of class `mbtransfer_model`, as generated using the
#'   `mbtransfer` function. This includes trained boosting models for every
#'   taxon, stored within the `@parameters` slot.
#' @param newdata A new `ts_inter` object over which to perform prediction. This
#'   method will make predictions for every timepoint that appears in the
#'   `@interventions` slot but not the `@values`.
#' @export
#' @examples
#' data(sim_ts)
#' fit <- mbtransfer(sim_ts)
#' ts_subset <- subset_values(sim_ts, 1:25)
#' predict(fit, ts_subset)
mbtransfer_predict <- function(object, newdata) {
  lags <- unlist(object@hyper[c("P", "Q")])
  result <- list()
  subject <- subject_data(newdata) |>
    select(-subject) |>
    as.matrix()

  for (i in seq_along(newdata)) {
    result[[i]] <- model_predict_single(
      object@parameters,
      newdata[[i]],
      lags,
      subject[i, , drop = FALSE],
      object@interactions
    )
  }

  names(result) <- names(newdata)
  new("ts_inter", series = result, subject_data = subject_data(newdata))
}

#' Prediction for a Single Subject
#'
#' This loops over all timepoints for a single subject and makes predictions by
#' comparing the number of timepoints in @interventions and in @values. The gap
#' will be filled in one step at a time using `mbtransfer_predict_step()`.
#'
#' @param fit An object of class `mbtransfer_model`, as generated using the
#'   `mbtransfer` function. This includes trained boosting models for every
#'   taxon, stored within the `@parameters` slot.
#' @param ts_inter A new `ts_inter_single` object over which to perform
#'   prediction. This method will make predictions for every timepoint that
#'   appears in the `@interventions` slot but not the `@values`. This is assumed
#'   to be a single subject from a larger `ts_inter` object.
#' @param lags A vector specifying `P` and `Q` in the trained mbtransfer model.
#' @param subject A static data frame of subject-level variables. This will be
#'   concatenated to time-varying intervention and taxonomic covariates when
#'   making predictions. This is analogous to the training process.
model_predict_single <- function(fit, ts_inter, lags, subject = NULL, interactions = NULL) {
  n_time <- ncol(ts_inter)
  w <- interventions(ts_inter)
  while (ncol(ts_inter) < ncol(w)) {
    ts_inter <- model_predict_step(ts_inter, fit, lags, subject, interactions)
  }

  colnames(values(ts_inter)) <- colnames(w)
  ts_inter
}

#' Predictions for a single timepoint and subject
#'
#' This make predictions for all taxa for a single timestep ahead in one
#' subject. It loops over the trained boosting models for each taxon predicts a
#' single value for each.
#'
#' @param fit An object of class `mbtransfer_model`, as generated using the
#'   `mbtransfer` function. This includes trained boosting models for every
#'   taxon, stored within the `@parameters` slot.
#' @param ts_inter A new `ts_inter_single` object over which to perform
#'   prediction. This method will make predictions for every timepoint that
#'   appears in the `@interventions` slot but not the `@values`. This is assumed
#'   to be a single subject from a larger `ts_inter` object.
#' @param lags A vector specifying `P` and `Q` in the trained mbtransfer model.
#' @param subject A static data frame of subject-level variables. This will be
#'   concatenated to time-varying intervention and taxonomic covariates when
#'   making predictions. This is analogous to the training process.
model_predict_step <- function(ts_inter, fit, lags, subject = NULL, interactions = NULL) {
  xz <- predictors(ts_inter, lags, subject, interactions)
  y_hat <- vector(length = nrow(ts_inter))
  for (j in seq_len(nrow(ts_inter))) {
    y_hat[j] <- predict(fit[[j]], xz)
  }

  values(ts_inter) <- cbind(values(ts_inter), y_hat)
  ts_inter
}

#' Transfer Function Model Class
#' @slot parameters The model training result. For mbtransfer models, this is a
#'   list of all the boosting models that have been trained, one for each taxon.
#' @slot method The model type used for prediction.
#' @slot hyper A list. containing any hyperparameters used during model
#'   training.
#' @export
setClass(
  "mbtransfer_model",
  slots = c(
    parameters = "ANY",
    method = "character",
    hyper = "list",
    interactions = "ANY"
  )
)

#' Prediction method for mbtransfer models
#' @examples
#' data(sim_ts)
#' fit <- mbtransfer(sim_ts, 1, 1)
#' sim_sub <- subset_values(sim_ts, 1:26) # remove last timepoint
#' y_hat <- predict(fit, sim_sub)
#' plot(values(sim_ts[[1]])[, 27], values(y_hat[[1]])[, 27], xlab = "Truth", ylab = "Predicted")
#' @export
setMethod("predict", c(object = "mbtransfer_model"), mbtransfer_predict)
