setClassUnion("data.frameOrNull", members = c("data.frame", "NULL"))

#' An S4 class to represent a single subject's series
#'
#' Taxonomic abundances, binary/continuous interventions, and underlying
#' timepoints for a single subject are stored in this unified data structure. A
#' `ts_inter` object is list of these `ts_inter_single` classes (together wtih
#' static subject data).
#'
#' @slot values A matrix whose rows are taxa and whose columns are samples.
#' @slot time A vector of the timepoints associated with the samples in `values`.
#' @slot interventions A matrix whose rows are perturbations, columns are
#'   samples, and values are either binary interventions or continuous input
#'   series, representing the value of the exogenous influence.
#'
#' @examples
#' data(sim_ts)
#' sim_ts[[1]]
#' @export
setClass(
  "ts_inter_single",
  slots = c(
    values = "matrix",
    time = "numeric",
    interventions = "matrix"
  )
)

#' An S4 class representing a collection of time series under environmental shifts
#'
#' While `ts_inter_single` represents shifts for a single subject, `ts_inter`
#' includes shifts across a collection of subjects. It includes an additional
#' slot, `subject_data` to store non-time-varying subject metadata.
#'
#' @slot series A list of `ts_inter_single` objects, each representing
#'   interventions and compositional responses for a subject.
#' @slot subject_data An optional data.frame storing static host-level metadata
#'   associated with each series.
#' @export
#' @examples
#' data(sim_ts)
#' sim_ts
setClass(
  "ts_inter",
  slots = c(
    series = "list",
    subject_data = "data.frameOrNull"
  )
)

single_subset <- function(x, i = NULL, j = NULL, ..., drop = FALSE) {
  if (!missing(i)) {
    values <- x@values[i, j, drop = drop]
  } else {
    values <- x@values[, j, drop = drop]
  }

  time <- x@time[j]
  interventions <- x@interventions[, j, drop = drop]
  new("ts_inter_single", values = values, time = time, interventions = interventions)
}

multi_subset <- function(x, i = NULL, j = NULL, ..., drop = FALSE) {
  result <- list()
  for (k in seq_along(x)) {
    if (!missing(i)) {
      result[[k]] <- x[[k]][i, j]
    } else {
      result[[k]] <- x[[k]][, j]
    }
  }

  ts <- new("ts_inter", series = result, subject_data = x@subject_data)
  names(ts) <- names(x)
  ts
}

#' Subset values of a ts object
#'
#' This is a helper to filter the number of timepoints across all `values` slots
#' in a `ts_inter` class.
#'
#' @examples
#' data(sim_ts)
#' subset_values(sim_ts, 1:2)
#' subset_values(sim_ts, 4:8)
#' @export
subset_values <- function(ts, ix) {
  ts_missing <- ts
  for (i in seq_along(ts)) {
    values(ts_missing[[i]]) <- values(ts_missing[[i]][, ix])
  }
  ts_missing
}

setMethod("length", "ts_inter", function(x) length(x@series))
setMethod("nrow", "ts_inter_single", function(x) nrow(x@values))
setMethod("ncol", "ts_inter_single", function(x) ncol(x@values))
setMethod("colnames", "ts_inter_single", function(x) colnames(x@values))
setMethod("rownames", "ts_inter_single", function(x) rownames(x@values))
setMethod("dim", "ts_inter_single", function(x) dim(x@values))
setMethod("[", c("ts_inter_single", "numeric", "numeric"), single_subset)
setMethod("[", c("ts_inter_single", "numeric", "missing"), single_subset)
setMethod("[", c("ts_inter_single", "missing", "numeric"), single_subset)
setMethod("[", c("ts_inter", "numeric", "missing", "ANY"), function(x, i, j, ..., drop = TRUE) initialize(x, series = x@series[i], subject_data = x@subject_data[i, ]))
setMethod("[[", c("ts_inter", "numeric", "missing"), function(x, i, j, ...) x@series[[i]])
setMethod("[", c("ts_inter", "logical", "missing", "ANY"), function(x, i, j, ..., drop = TRUE) initialize(x, series = x@series[i], subject_data = x@subject_data[i, ]))
setMethod("[[", c("ts_inter", "logical", "missing"), function(x, i, j, ...) x@series[[i]])
setMethod("[", c("ts_inter", "missing", "numeric", "ANY"), multi_subset)
setMethod("[", c("ts_inter", "numeric", "numeric", "ANY"), multi_subset)

#' @export
setGeneric("values", function(x) standardGeneric("values"))
setMethod("values", "ts_inter_single", function(x) x@values)

#' @export
setGeneric("values<-", function(x, values) standardGeneric("values<-"))
setMethod("values<-", "ts_inter_single", function(x, values) {
  x@values <- values
  x
})

setMethod("[[<-", "ts_inter", function(x, i, j, value) {
  x@series[[i]] <- value
  x
})

#' @importFrom utils head
print_ts_inter <- function(object) {
  n_time <- unlist(lapply(object@series, function(x) ncol(values(x))))
  cat(sprintf(
    "# A ts_inter object | %d taxa | %d subjects | %.2f \U00B1 %.2f timepoints:\n",
    length(taxa(object)), length(object), mean(n_time), 1.9 * sd(n_time)
  ))
  for (i in seq_len(min(3, length(n_time)))) {
    n_col <- min(4, ncol(values(object[[i]])))

    cat(sprintf("\n%s:\n", names(object)[i]))
    v <- data.frame(round(values(object[[i]])[1:n_col, 1:n_col], 3))
    v <- cbind(v, " " = rep("\U2026", n_col))
    v <- rbind(v, " " = c(rep("\U22EE", n_col), "\U22F1"))
    print(v)
  }

  if (length(n_time) > 3) {
    cat(sprintf("\nand %d more subjects.", length(n_time) - 3))
  }
}

setMethod("show", "ts_inter", print_ts_inter)

print_ts_inter_single <- function(object) {
  cat(sprintf(
    "# A ts_inter_single object | %d taxa | %d timepoints:\n",
    length(nrow(values(object))), ncol(values(object))
  ))

  cat("taxa:\n")
  n_col <- min(4, ncol(values(object)))
  v <- data.frame(round(values(object)[1:n_col, 1:n_col], 3))
  v <- cbind(v, " " = rep("\U2026", n_col))
  v <- rbind(v, " " = c(rep("\U22EE", n_col), "\U22F1"))
  print(v)

  cat("interventions:\n")
  i <- interventions(object)[, 1:n_col, drop = FALSE]
  i <- as.data.frame(i)
  i <- cbind(i, " " = "\U2026")
  print(i)
}

setMethod("show", "ts_inter_single", print_ts_inter_single)


#' @export
setGeneric("interventions", function(x) standardGeneric("interventions"))
setMethod("interventions", "ts_inter_single", function(x) x@interventions)

#' @export
setGeneric("interventions<-", function(x, value) standardGeneric("interventions<-"))
setMethod("interventions<-", "ts_inter_single", function(x, value) {
  x@interventions <- value
  x
})

#' @export
setGeneric("times", function(x) standardGeneric("times"))
setMethod("times", "ts_inter_single", function(x) x@time)

#' @export
setGeneric("subject_data", function(x) standardGeneric("subject_data"))
setMethod("subject_data", "ts_inter", function(x) x@subject_data)

#' @export
setGeneric("subject_data<-", function(x, value) standardGeneric("subject_data<-"))
setMethod("subject_data<-", "ts_inter", function(x, value) {
  x@subject_data <- value
  x
})

#' @export
setGeneric("taxa", function(x) standardGeneric("taxa"))
setMethod("taxa", "ts_inter", function(x) {
  rownames(values(x[[1]]))
})


#' @export
setMethod("names", "ts_inter", function(x) names(x@series))

#' @export
setMethod("names<-", "ts_inter", function(x, value) {
  names(x@series) <- value
  x
})
