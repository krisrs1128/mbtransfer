setClassUnion("data.frameOrNull", members=c("data.frame", "NULL"))

#' @export
setClass(
  "ts_inter_single", 
  slots = c(
    values = "matrix",
    time = "numeric",
    interventions = "matrix"
  )
)

#' @export
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
setMethod("[", c("ts_inter", "numeric", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) initialize(x, series=x@series[i], subject_data = x@subject_data[i, ]))
setMethod("[[", c("ts_inter", "numeric", "missing"), function(x, i, j, ...) x@series[[i]])
setMethod("[", c("ts_inter", "logical", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) initialize(x, series=x@series[i], subject_data = x@subject_data[i, ]))
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
