
#' @importFrom dplyr filter select mutate
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider unite
subject_order <- function(values_df, taxa, r = 0) {
  values_wide <- values_df |>
    filter(taxon %in% taxa) |>
    select(subject, taxon, value, time) |>
    mutate(time = round(time, r)) |>
    unite("txtime", taxon, time) |>
    pivot_wider(names_from = txtime) |>
    column_to_rownames("subject") |>
    as.matrix()
  
  values_wide[is.na(values_wide)] <- 0
  rownames(values_wide)[hclust(dist(values_wide))$order]
}

#' Visualize Time Series Heatmaps
#' 
#' This is a helper function to visualize time series for a subset of taxa
#' across all subjects. It is used in the raw data figures in the case studies
#' section of the accompanying manuscript.
#' 
#' @param values_df A data.frame containing joined information from a `ts_inter`
#'   object. Each row should be one count from one taxon/subject pair. The
#'   output of `pivot_ts` reshape data to this format.
#' @param taxa A character vector of taxa to filter down to in the final
#'   visualization.
#' @param condition A subject-level variable to use during faceting. Defaults to
#'   no faceting.
#' @param r When ordering subjects, we consider averages after rounding
#'   timepoints into bins. Large `r` corresponds to fine-grained binning,
#'   negative `r` are more coarse bins. Defaults to 0. See `subject_order` for
#'   details.
#' 
#' @importFrom ggplot2 scale_color_gradient scale_fill_gradient
#' @examples
#' library(dplyr)
#' data(sim_ts)
#' pivoted <- sim_ts |>
#'   pivot_ts() |>
#'   mutate(v_pos = V1 > 0)
#' interaction_hm(pivoted, c("tax1", "tax2"), "v_pos")
#' 
#' pivoted |>
#'   group_by(taxon) |>
#'   mutate(value = rank(value) / n()) |>
#'   interaction_hm(c("tax1", "tax2"), "v_pos")
#' @export
interaction_hm <- function(values_df, taxa, condition = NULL, r = 0, ...) {
  p <- values_df |>
    filter(taxon %in% taxa) |>
    mutate(subject = factor(subject, subject_order(values_df, taxa, r))) |>
    ggplot() +
    geom_tile(aes(time, subject, fill = value, col = value), ...) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_gradient(low = "#eaf7f7", high = "#037F8C") +
    scale_fill_gradient(low = "#eaf7f7", high = "#037F8C") +
    theme(
      strip.text.y = element_text(angle = 0),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 1, fill = NA, color = "#545454"),
      panel.spacing = unit(0, "cm")
    )
  
  if (!is.null(condition)) {
    p <- p + 
      facet_grid(.data[[condition]] ~ reorder(taxon, -value), scales = "free", space = "free")
  }
  
  p
}

#' @importFrom ggplot2 facet_grid ggplot geom_tile scale_x_continuous
#'   scale_color_gradient scale_fill_gradient
#' @export
interaction_barcode <- function(values_df, taxa, condition = NULL, r = 0, ...) {
  p <- values_df |>
    filter(taxon %in% taxa) |>
    mutate(subject = factor(subject, levels = subject_order(values_df, taxa, r))) |>
    ggplot() +
    geom_tile(aes(time, taxon, fill = value, col = value), ...) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_gradient(low = "#eaf7f7", high = "#037F8C") +
    scale_fill_gradient(low = "#eaf7f7", high = "#037F8C") +
    theme(
      panel.spacing = unit(0, "line"),
      panel.border = element_rect(linewidth = 1, fill = NA, color = "#d3d3d3"),
      strip.text.y = element_text(angle = 0),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  if (!is.null(condition)) {
    p <- p + 
      #ggh4x::facet_nested(reorder(.data[[condition]], -value) + subject ~ ., scales = "free", space = "free")
      facet_grid(reorder(.data[[condition]], -value) + subject ~ ., scales = "free", space = "free")
  }
  
  p
}

#' @importFrom dplyr filter group_by summarise
#' @export
ribbon_data <- function(ts1, ts0, focus_taxa = NULL, delta = NULL, q_lower = 0.25, q_upper = 0.75) {
  result <- (ts1 - ts0) |>
    pivot_ts()
    
  if (!is.null(focus_taxa)) {
    result <- filter(result, taxon %in% focus_taxa)
  }
  if (!is.null(delta)) {
    result <- result |>
      mutate(time = delta * round(time / delta))
  }

  result |>
    group_by(taxon, time) |>
    summarise(
      q_lower = quantile(value, q_lower),
      median = median(value),
      q_upper = quantile(value, q_upper)
    )
}

#' @importFrom ggplot2 ggplot geom_ribbon aes geom_hline geom_line
#'   scale_x_continuous theme facet_wrap
#' @export
ribbon_plot <- function(rdata, group = NULL, reorder_var = NULL) {
  p <- ggplot(rdata, aes(x = time)) +
    geom_hline(yintercept = 0, size = .25, col = "#787878") +
    scale_x_continuous(expand = c(0, 0)) +
    theme(axis.text = element_text(size = 8))
  
  # filled or grey ribbons
  if (!is.null(group)) {
    p <- p + 
      geom_ribbon(aes(ymin = q_lower, ymax = q_upper, fill = .data[[group]], group = .data[[group]]), alpha = 0.6) +
      geom_line(aes(y = median, col = .data[[group]], group = .data[[group]]), size = 1) +
      scale_fill_brewer(palette = "Set2") +
      scale_color_brewer(palette = "Set2")
  } else {
    p <- p + 
      geom_ribbon(aes(ymin = q_lower, ymax = q_upper), fill = "#d3d3d3") +
      geom_line(aes(y = median), col = "#8e8e8e", size = 1)
  }
  
  # order facets (or not)
  if (!is.null(reorder_var)) {
    p <- p + 
      facet_wrap(~ reorder(taxon, .data[[reorder_var]]), scales = "free_y")
  } else if (n_distinct(rdata$taxon) > 1) {
    p <- p +
      facet_wrap(~ taxon, scales = "free_y")
  }
  
  p
}

#' @importFrom dplyr rename select mutate full_join left_join group_by ungroup
#'   summarise
#' @export
reshape_preds <- function(ts, ts_pred, n_quantile = 4, lag = 3) {
  ts_df <- pivot_ts(ts_pred) |>
    rename(y_hat = value) |>
    select(taxon, sample, time, y_hat) |>
    mutate(h = time - lag) |>
    full_join(pivot_ts(ts)) |>
    rename(y = value)
  
  taxa_totals <- ts_df |>
    group_by(taxon) |>
    summarise(total = median(y)) |>
    ungroup() |>
    mutate(quantile = cut(total, unique(quantile(total, 0:n_quantile/n_quantile)), include.lowest = TRUE))
  
  ts_df |>
    left_join(taxa_totals)
}