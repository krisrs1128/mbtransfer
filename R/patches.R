#' Sliding Windows for a ts_inter_single object
#'
#' This creates sliding windows of intervention and community composition
#' features for a single subject. It returns a list giving the "patchified"
#' data, which can then be organized into matrices for prediction.
#'
#' @param ts_inter An object of class `ts_inter_single` over which to generate
#'   sliding windows.
#' @param p The number of time lags to. use in the sliding window for the
#'   microbiome features.
#' @param q The number of time lags to use in the sliding window for
#'   interventions.
#' @importFrom slider slide
#' @export
#' @examples
#' data(sim_ts)
#' patches <- patchify_single(sim_ts[[1]])
#' head(patches$x[[1]])
#' head(patches$w[[1]])
#' head(patches$y[[1]])
patchify_single <- function(ts_inter, p = 2, q = 3) {
  ix <- seq_len(ncol(ts_inter))
  x_indices <- slide(ix, ~., .before = p, .after = -1)
  w_indices <- slide(ix, ~., .before = q - 1, .after = 0)
  y_indices <- slide(ix, ~., .before = -1, .after = 1)

  # initialize result
  k <- 1
  data <- replicate(3, list())
  names(data) <- c("x", "y", "w")
  values_ <- values(ts_inter)
  interventions_ <- interventions(ts_inter)

  # extract x (taxa), w (intervention), and y (future taxa) patches
  for (i in seq_along(x_indices)) {
    if (length(x_indices[[i]]) == p &
      length(w_indices[[i]]) == q &
      length(y_indices[[i]]) == 1) {
      data$x[[k]] <- cbind(values_[, x_indices[[i]], drop = FALSE])
      data$w[[k]] <- interventions_[, w_indices[[i]], drop = FALSE]
      data$y[[k]] <- values_[, y_indices[[i]], drop = FALSE]
      k <- k + 1
    }
  }

  data
}

#' Convert `patchify_single` output into a `data.frame`
#'
#' This generates a data.frame of sliding window predictors and responses by
#' wrapping `patchify_single`. It is in a format that can directly be used by
#' different regression methods.
#'
#' @param ts_inter An object of class `ts_inter_single` over which to generate
#'   sliding windows.
#' @param p The number of time lags to. use in the sliding window for the
#'   microbiome features.
#' @param q The number of time lags to use in the sliding window for
#'   interventions.
#'
#' @examples
#' data(sim_ts)
#' patches <- patchify_single_df(sim_ts[[1]], 2, 2)
#' head(patches$x)
#' head(patches$y)
#' @importFrom glue glue
#' @export
patchify_single_df <- function(ts_inter, p, q) {
  data <- patchify_single(ts_inter, p, q)
  x <- data$x
  y <- data$y
  w <- data$w

  result <- list(
    x = matrix(nrow = length(x), ncol = length(x[[1]]) + length(w[[1]])),
    y = matrix(nrow = length(x), ncol = length(y[[1]]))
  )

  for (i in seq_along(x)) {
    xi <- matrix(x[[i]], nrow = 1)
    result$y[i, ] <- matrix(y[[i]], nrow = 1)
    wi <- matrix(w[[i]], nrow = 1)
    result$x[i, ] <- cbind(xi, wi)
  }

  colnames(result$x) <- predictor_names(dim(x[[1]]), dim(w[[1]]))
  colnames(result$y) <- glue("taxon{seq_len(nrow(y[[1]]))}")
  result
}

#' Search for Candidate Interactions
#' @importFrom insight check_if_installed
#' @importFrom xyz xyz_search
interaction_search <- function(x, y, interactions = "none", ...) {
  if (interactions == "none") {
    return(NULL)
  } else if (interactions == "search") {
    check_if_installed("xyz", "to search for candidate interaction effects. Please install using devtools::install_github('gathanei/xyz)")
    interactions <- t(xyz_search(x, rowMeans(y), binary = FALSE, ...)[[1]])
    return(interactions[interactions[, 1] != interactions[, 2], ])
  }
  stop("'interactions' must be one of \"none\" or \"search\".")
}

#' Add Interactions to Predictor Matrix
#' @param x The original matrix containing the source of the interaction terms.
#' @param interactions An n_interactions x 2 matrix containing interactions
#'   to add to the original matrix x.
#' @return The matrix x with new columns associated wtih each interaction.
append_interactions <- function(x, interactions = NULL) {
  new_cols <- list()
  for (k in 1:nrow(interactions)) {
    nm <- paste0(colnames(x)[interactions[k, ]], collapse = "*")
    new_cols[[nm]] <- apply(x[, interactions[k, ], drop = FALSE], 1, prod)
  }
  new_cols <- do.call(cbind, new_cols)
  return(cbind(x, new_cols))
}


#' Sliding Windows for a ts_inter object
#'
#' This creates sliding windows of intervention and community composition
#' features for all. subjects in a `ts_inter` object. It returns a list giving
#' the "patchified" data. The first component, `x`, gives the microbiome and
#' intervention features immediately preceding the values in the second
#' component, `y`. This is constructed by running `patchify_single_df` over all
#' subjects in the dataset.
#'
#' @param ts_inter An object of class `ts_inter` over which to generate
#'   sliding windows.
#' @param p The number of time lags to. use in the sliding window for the
#'   microbiome features.
#' @param q The number of time lags to use in the sliding window for
#'   interventions.
#' @importFrom purrr map_dfr
#' @importFrom tibble as_tibble
#' @importFrom dplyr select bind_cols
#' @examples
#' data(sim_ts)
#' result <- patchify_df(sim_ts)
#' lapply(result, head)
#' @export
patchify_df <- function(ts_inter, p = 2, q = 3, interactions = "none") {
  patches <- list()
  for (i in seq_along(ts_inter)) {
    patches[[i]] <- patchify_single_df(ts_inter[[i]], p, q)
    sdata <- subject_data(ts_inter)
    if (!is.null(sdata)) {
      patches[[i]]$x <- patches[[i]]$x |>
        bind_cols(select(sdata[i, ], -subject))
    }
  }

  x <- map_dfr(patches, ~ as_tibble(.$x)) |>
    as.matrix()

  y <- map_dfr(patches, ~ as_tibble(.$y))
  z <- interaction_search(x, y, interactions)
  list(x = x, y = y, interactions = z)
}

#' Create uniform column names for patchified df
#'
#' @examples
#' data(sim_ts)
#' predictor_names(c(2, 3), c(2, 2))
#' predictor_names(c(5, 2), c(1, 1))
#' @importFrom glue glue
#' @export
predictor_names <- function(x_dim, w_dim) {
  n1 <- rep(glue("taxon{seq_len(x_dim[1])}"), x_dim[2])
  n2 <- rep(glue("lag{seq(x_dim[2], 1)}"), each = x_dim[1])
  x_names <- paste(n1, n2, sep = "_")

  n1 <- rep(glue("intervention{seq_len(w_dim[1])}"), w_dim[2])
  n2 <- rep(glue("lag{seq(w_dim[2] - 1, 0)}"), each = w_dim[1])
  w_names <- paste(n1, n2, sep = "_")

  c(x_names, w_names)
}

#' Detect the time lags based on column names
#'
#' This is a helper function for deciding on P and Q in transfer function models
#' using just the names of an example set of covariates.
#' @importFrom utils strcapture
#' @examples
#' lag_from_names(c("taxon1_lag1", "taxon1_lag2", "taxon1_lag3"))
#' @export
lag_from_names <- function(names, group = "taxon") {
  names_ix <- names[grepl(group, names)]
  regex <- paste0(group, "([0-9]+_lag[0-9]+)")
  names_value <- strcapture(regex, names_ix, data.frame(chr = character()))
  gsub("[0-9]+_lag", "", names_value$chr) |>
    as.numeric() |>
    max()
}

#' Detect the time lags over which a model was trained
#'
#' This is a helper for determining the P, Q associated with a fitted mbtransfer
#' model.
#'
#' @examples
#' data(sim_ts)
#' fit <- mbtransfer(sim_ts)
#' time_lags(fit@parameters[[1]])
#' @export
time_lags <- function(fit) {
  # different names for gbm and lasso, resp.
  if (!is.null(fit$feature_names)) {
    inputs <- fit$feature_names
  } else {
    inputs <- rownames(fit$beta)
  }

  P <- lag_from_names(inputs, "taxon")
  Q <- lag_from_names(inputs, "intervention") + 1
  c(P, Q)
}

pad_lag <- function(x, lag) {
  if (ncol(x) < lag) {
    x <- cbind(matrix(0, nrow(x), lag - ncol(x)), x)
  }
  x
}

#' Generate a Predictor Matrix
#'
#' `patchify_df` is used to generate training data. At test time, we don't need
#' all the patches, but we do need to construct regressors for the immediate
#' past. This will be matched to the model's original training data and used to
#' generate predictions.
#'
#' @param ts_inter An object of class `ts_inter_single` over which to generate
#'   sliding windows.
#' @param lags A vector specifying `P` and `Q` in the trained mbtransfer model.
#' @param subject A static data frame of subject-level variables. This will be
#'   concatenated to time-varying intervention and taxonomic covariates when
#'   making predictions. This is analogous to the training process.
#' @importFrom purrr set_names
#' @examples
#' data(sim_ts)
#' ts <- subset_values(sim_ts, 1:5)
#' predictors(ts[[1]], c(2, 2), NULL)
#' @export
predictors <- function(ts_inter, lags, subject, interactions = NULL) {
  x <- values(ts_inter) |>
    pad_lag(lags[1])
  w <- interventions(ts_inter) |>
    pad_lag(lags[2])
  n_time <- ncol(x)

  x_prev <- x[, seq(n_time - lags[1] + 1, n_time), drop = FALSE]
  w_prev <- w[, seq(n_time - lags[2] + 2, n_time + 1), drop = FALSE]
  xw <- cbind(matrix(x_prev, nrow = 1), matrix(w_prev, nrow = 1)) |>
    as.data.frame() |>
    set_names(predictor_names(dim(x_prev), dim(w_prev))) |>
    as.matrix()

  if (!is.null(subject)) {
    xw <- cbind(xw, subject[rep(1, nrow(xw)), , drop = FALSE])
  }
  if (!is.null(interactions)) {
    xw <- append_interactions(xw, interactions)
  }

  xw
}
