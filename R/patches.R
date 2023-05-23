
#' @importFrom slider slide
patchify_single <- function(ts_inter, p = 2, q = 3) {
  ix <- seq_len(ncol(ts_inter))
  x_indices <- slide(ix, ~ ., .before = p, .after = -1)
  w_indices <- slide(ix, ~ ., .before = q - 1, .after = 0)
  y_indices <- slide(ix, ~ ., .before = -1, .after = 1)
  
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

#' @importFrom stringr str_c
patchify_single_df <- function(ts_inter, p, q) {
  data <- patchify_single(ts_inter, p, q)
  x <- data$x
  y <- data$y
  w <- data$w
  
  result <- list(
    x = matrix(nrow = length(x), ncol = 1 + length(x[[1]]) + length(w[[1]])),
    y = matrix(nrow = length(x), ncol = length(y[[1]]))
  )
  
  for (i in seq_along(x)) {
    xi <- matrix(x[[i]], nrow = 1)
    result$y[i, ] <- matrix(y[[i]], nrow = 1)
    wi <- matrix(w[[i]], nrow = 1)
    result$x[i, ] <- cbind(1, xi, wi)
  }
  
  colnames(result$x) <- predictor_names(dim(x[[1]]), dim(w[[1]]))
  colnames(result$y) <- str_c("taxon", seq_len(nrow(y[[1]])))
  result
}

#' @importFrom purrr map_dfr
#' @importFrom tibble as_tibble
#' @importFrom dplyr select bind_cols
#' @export
patchify_df <- function(ts_inter, p = 2, q = 3) {
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
  list(x = x, y = y)
}

predictor_names <- function(x_dim, w_dim) {
  n1 <- rep(str_c("taxon", seq_len(x_dim[1])), x_dim[2])
  n2 <- rep(str_c("lag", seq(x_dim[2], 1)), each = x_dim[1])
  x_names <- c("intercept", str_c(n1, "_", n2))
  
  n1 <- rep(str_c("intervention", seq_len(w_dim[1])), w_dim[2])
  n2 <- rep(str_c("lag", seq(w_dim[2] - 1, 0)), each = w_dim[1])
  w_names <- str_c(n1, "_", n2)
  
  c(x_names, w_names)
}

#' @importFrom stringr str_extract str_detect str_remove
lag_from_names <- function(names, group = "taxon") {
  names[str_detect(names, group)] |>
    str_extract("lag[0-9]+") |>
    str_remove("lag") |>
    as.numeric() |>
    max()
}

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
#' @importFrom purrr set_names
predictors <- function(ts_inter, lags, subject) {
  x <- values(ts_inter) |>
    pad_lag(lags[1])
  w <- interventions(ts_inter) |>
    pad_lag(lags[2])
  n_time <- ncol(x)
  
  x_prev <- x[, seq(n_time - lags[1] + 1, by = 1, length.out = lags[1]), drop = FALSE]
  w_prev <- w[, seq(n_time - lags[2] + 2, by = 1, length.out = lags[2]), drop = FALSE]
  
  xw <- cbind(intercept = 1, matrix(x_prev, nrow = 1), matrix(w_prev, nrow = 1)) |>
    as.data.frame() |>
    set_names(predictor_names(dim(x_prev), dim(w_prev))) |>
    as.matrix()
  
  if (!is.null(subject)) {
    xw  <- cbind(xw, subject[rep(1, nrow(xw)),, drop = FALSE])
  }
  
  xw
}
