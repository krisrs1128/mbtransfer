#' @importFrom glue glue
#' @export
interpolate <- function(ts_inter, delta = 1, method = "constant") {
  for (i in seq_along(ts_inter)) {
    ts_inter[[i]] <- interpolate_(ts_inter[[i]], delta, method)
    colnames(values(ts_inter[[i]])) <- glue("{names(ts_inter)[i]}_{colnames(values(ts_inter[[i]]))}")
    colnames(interventions(ts_inter[[i]])) <- glue("{names(ts_inter)[i]}_{colnames(interventions(ts_inter[[i]]))}")
  }
  ts_inter
}

approx_mat <- function(time, y_mat, times_out, method) {
  n <- length(times_out)
  new_y <- matrix(nrow = nrow(y_mat), ncol = n, dimnames = list(rownames(y_mat), seq_len(n)))

  for (i in seq_len(nrow(y_mat))) {
    new_y[i, ] <- approx(
      time,
      y_mat[i, ],
      times_out,
      method = method,
      ties = mean
    )$y
  }

  colnames(new_y) <- str_c("T", seq_len(ncol(new_y)))
  new_y
}

interpolate_ <- function(ts_inter_single, delta, method) {
  time <- ts_inter_single@time
  times_out <- seq(min(time), max(time), by = delta)
  ts_inter_single@time <- times_out

  v <- values(ts_inter_single)
  inter <- interventions(ts_inter_single)
  values(ts_inter_single) <- approx_mat(time, v, times_out, method)
  interventions(ts_inter_single) <- approx_mat(time, inter, times_out, method)
  ts_inter_single
}
