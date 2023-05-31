
###############################################################################
# Helpers to generate counterfactual ts_inter objects
###############################################################################

to_vector <- function(x) {
  if (length(x) == 1) {
    x <- c(x)
  }
  x
}

#' Hypothetical Step Interventions
#' @examples
#' library(mbtransfer)
#' steps(c("P1" = TRUE), 1:3, 2:3, 4)
#' 
#' @importFrom glue glue
#' @export
steps <- function(p_states, starts = 1, lengths = 1:3, L = 3, w_star = c(0, 1)) {
  w0 <- matrix(0, nrow = length(p_states), ncol = L)
  rownames(w0) <- names(p_states)
  colnames(w0) <- glue("t{seq_len(ncol(w0))}")
  starts <- to_vector(starts)
  lengths <- to_vector(lengths)
  active_p <- names(p_states[p_states])

  result <- list()
  k <- 1
  for (i in seq_along(w_star)) {
    for (s in starts) {
      for (l in lengths) {
        wi <- w0
        for (p in active_p) {
          wi[p, seq(s, min(s + l - 1, ncol(wi)))] <- w_star[i]
        }
        result[[k]] <- wi
        k <- k + 1
      }
    }
  }
  
  unique(result)
}

#' Hypothetical Pulse Interventions
#' 
#' Create counterfactual perturbations with specific 0/1 intervention and total
#' length patterns.
#' 
#' @param p_stats A list specifying the row names and whether they should
#'   include interventions. Any name set to TRUE will include intervention
#'   effects, while those with FALSE will not.
#' @param lags The time length of any interventions that we generate.
#' @param L The total length returned, including both intervention and non-intervention timepoints
#' @param w_star The unique values to include in and out of the intervention.
#'   Defaults to 1/0, respectively.
#' 
#' @examples
#' library(mbtransfer)
#' pulses(c("P1" = TRUE), 1, 4)
#' pulses(c("P1" = TRUE), 1:3, 4)
#" pulses(c("P1" = TRUE), 1:3, 4, seq(0, 1, .2))
#' @importFrom glue glue
#' @export
pulses <- function(p_states, lags = 1, L = 3, w_star = c(0, 1)) {
  w0 <- matrix(0, nrow = nrow(inter), ncol = L)
  rownames(w0) <- rownames(inter)
  colnames(w0) <- glue("Tn_{seq_len(ncol(w0))}")
  lengths <- to_vector(lengths)
  active_p <- names(p_states[p_states])
  
  result <- list()
  k <- 1
  for (i in seq_along(w_star)) {
    for (l in seq_along(lags)) {
      wi <- w0
      for (p in seq_along(active_p)) {
        wi[p, L - l + 1] <- w_star[i]
      }
      result[[k]] <- wi
      k <- k + 1
    }
  }
  
  result
}

#' Generate Counterfactual versions of a ts_inter object
#' @export
counterfactual_ts <- function(ts, w0, w1, start_ix = NULL) {
  if (is.null(start_ix)) {
    start_ix <- map_dbl(ts, ncol)
  } else if (length(start_ix) == 1) {
    start_ix <- rep(start_ix, length(ts))
  }
  
  for (i in seq_along(ts)) {
    times_ <- ts[[i]]@time
    ts[[i]] <- ts[[i]][, seq_len(start_ix[i] - 1)]
    ts[[i]]@time <- times_[seq_len(start_ix[i] + ncol(w0) - 1)]
  }

  ts0 <- replace_inter(ts, w0, start_ix)
  ts1 <- replace_inter(ts, w1, start_ix)
  list(ts0 = ts0, ts1 = ts1)
}
