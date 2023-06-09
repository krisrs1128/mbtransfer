
#' ts is a single element of a ts_inter class
replace_inter_ <- function(ts, new_inter, start_ix = NULL) {
  inter <- interventions(ts)[, seq_len(start_ix - 1), drop = FALSE]
  interventions(ts) <- cbind(inter, new_inter)
  ts
}

#' @importFrom glue glue
#' @export
replace_inter <- function(ts, new_inter, start_ix = NULL) {
  if (length(start_ix) == 1) {
    start_ix <- rep(start_ix, length(ts))
  }
  
  for (i in seq_along(ts)) {
    inter_ <- new_inter
    colnames(inter_) <- glue("{names(ts)[i]}_{colnames(inter_)}")
    ts[[i]] <- replace_inter_(ts[[i]], inter_, start_ix[i])
  }
  
  ts
}

#' @export
replace_subject <- function(ts, new_subject) {
  subject <- subject_data(ts)
  shared_cols <- intersect(colnames(subject), colnames(new_subject))
  subject[, shared_cols] <- new_subject
  subject_data(ts) <- subject
  ts
}

#' Randomly subset a `ts_inter` object
#' 
#' There are times we want to randomly sample patches from across timepoints and
#' subjects. This is a helper function that does this random sampling adapted to
#' `ts_inter` objects, subsetting taxa and intervention timepoints
#' simultaneously.
#' 
#' @param ts A `ts_inter` object whose random patches we want to sample.
#' @param n The number of randomly sampled time series to return.
#' @param patch_len The length of the randomly sampled patches
#' @param intervention_len We will sample random patches padded by an extra
#'   intervention_len set of timepoints.
#' 
#' @examples
#' data(sim_ts)
#' sample_ts(sim_ts, 10, 4, 1)
#' @export
sample_ts <- function(ts, n, patch_len = 5, intervention_len = NULL) {
  # randomly subset series
  weights <- map_dbl(ts, ncol)
  weights <- weights / sum(weights)
  ix <- sample(seq_along(ts), n, replace = TRUE, prob = weights)
  ts_star <- ts[ix]
  
  # randomly subset windows
  for (i in seq_along(ts_star)) {
    start_ix <- sample(seq_len(ncol(ts_star[[i]]) - patch_len), 1)
    tmp <- ts_star[[i]]@time[seq(start_ix, start_ix + patch_len + intervention_len)]
    ts_star[[i]] <- ts_star[[i]][, seq(start_ix, start_ix + patch_len)]
    ts_star[[i]]@time <- tmp
  }
  
  ts_star
}

diff_ts_inter <- function(e1, e2) {
  result <- e1
  for (i in seq_along(e1)) {
    values(result[[i]]) <- values(e1[[i]]) - values(e2[[i]])
  }
  result
}

setMethod("-", c("ts_inter", "ts_inter"), diff_ts_inter)
