#' Plot individual longitudinal trajectories (spaghetti plot)
#'
#' @description
#' Draws one line per subject for a longitudinal biomarker, optionally
#' colored by the subject's eventual outcome, with a smoothed mean
#' trajectory per outcome group overlaid. With \code{align = "event"} the
#' horizontal axis becomes the time remaining until the subject's
#' event/censoring time (\code{survival_variable - time_variable}, so 0 is
#' the event time, at the left edge). This mirrors the backward joint
#' model's view of the biomarker as a function of time-to-event. Unlike
#' \code{\link{cmtPlot}}, which averages over subjects sharing one event
#' time, this shows every subject's raw history.
#'
#' @param data_plot_all A long-format \code{data.frame} with one row per
#' biomarker measurement.
#' @param bio_variable Name of the biomarker variable to plot.
#' @param time_variable Name of the measurement-time variable.
#' @param id_variable Name of the subject ID column. Default is \code{"id"}.
#' @param survival_variable Name of the time-to-event variable. Required when
#' \code{align = "event"}; otherwise ignored.
#' @param event_type_variable Name of the variable used to color subjects
#' (e.g. the event indicator or event type). Must be constant within each
#' subject. \code{NULL} (default) draws every subject in one color.
#' @param align Either \code{"baseline"} (default), plotting against
#' \code{time_variable}, or \code{"event"}, plotting against the time until
#' the event, \code{survival_variable - time_variable}.
#' @param n_subjects Optional number of subjects to draw, sampled at random
#' (call \code{set.seed()} first for a reproducible sample). \code{NULL}
#' (default) draws every subject. The smoothed means use the drawn subjects
#' only.
#' @param smooth Logical; overlay a LOESS mean trajectory per group. Default
#' \code{TRUE}.
#' @param alpha Transparency of the individual lines. Default \code{0.3}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(pbc3)
#'
#' # trajectories colored by death (1) vs alive/transplanted (0)
#' spaghettiPlot(pbc3, bio_variable = "serBilir", time_variable = "year",
#'   event_type_variable = "status2")
#'
#' # aligned at each subject's event/censoring time
#' spaghettiPlot(pbc3, bio_variable = "serBilir", time_variable = "year",
#'   survival_variable = "years", event_type_variable = "status2",
#'   align = "event")
#'
#' @export
spaghettiPlot <- function(data_plot_all, bio_variable, time_variable, id_variable = "id",
                          survival_variable = NULL, event_type_variable = NULL,
                          align = c("baseline", "event"), n_subjects = NULL,
                          smooth = TRUE, alpha = 0.3) {

  align <- match.arg(align)
  assert_data_frame(data_plot_all, "data_plot_all")
  assert_string(bio_variable, "bio_variable")
  assert_string(time_variable, "time_variable")
  assert_string(id_variable, "id_variable")
  vars <- c(bio_variable, time_variable, id_variable)
  if (align == "event") {
    if (is.null(survival_variable)) {
      stop("`survival_variable` must be supplied when `align = \"event\"`.", call. = FALSE)
    }
    assert_string(survival_variable, "survival_variable")
    vars <- c(vars, survival_variable)
  }
  if (!is.null(event_type_variable)) {
    assert_string(event_type_variable, "event_type_variable")
    vars <- c(vars, event_type_variable)
  }
  assert_vars_in_data(vars, data_plot_all, "spaghettiPlot() variables", "data_plot_all")
  if (!is.null(n_subjects)) {
    assert_positive_integer(n_subjects, "n_subjects")
  }
  assert_scalar_numeric(alpha, "alpha", positive = TRUE)

  d <- data_plot_all[!is.na(data_plot_all[[bio_variable]]) & !is.na(data_plot_all[[time_variable]]), ,
                     drop = FALSE]
  if (align == "event") {
    d <- d[!is.na(d[[survival_variable]]), , drop = FALSE]
  }
  if (nrow(d) == 0) {
    stop("`data_plot_all` has no rows with non-missing values of the plotted variables.", call. = FALSE)
  }

  ids <- unique(d[[id_variable]])
  if (!is.null(n_subjects) && n_subjects < length(ids)) {
    ids <- ids[sample.int(length(ids), n_subjects)]
    d <- d[d[[id_variable]] %in% ids, , drop = FALSE]
  }

  if (!is.null(event_type_variable)) {
    n_types <- tapply(d[[event_type_variable]], as.character(d[[id_variable]]),
                      function(v) length(unique(v)))
    if (any(n_types > 1)) {
      stop(sprintf("`%s` must be constant within each subject; it varies for %d subject(s).",
                   event_type_variable, sum(n_types > 1)), call. = FALSE)
    }
  }

  plot_df <- data.frame(
    id = as.character(d[[id_variable]]),
    x = if (align == "event") d[[survival_variable]] - d[[time_variable]] else d[[time_variable]],
    y = d[[bio_variable]],
    group = if (is.null(event_type_variable)) "All" else factor(d[[event_type_variable]], exclude = NULL)
  )
  plot_df <- plot_df[order(plot_df$id, plot_df$x), ]

  p <- ggplot(plot_df, aes(x = x, y = y, group = id))
  if (is.null(event_type_variable)) {
    p <- p +
      geom_line(alpha = alpha, color = "grey40") +
      geom_point(alpha = alpha, color = "grey40", size = 0.6)
    if (smooth) {
      p <- p + geom_smooth(aes(group = 1), method = "loess", formula = y ~ x,
                           color = "red", linewidth = 1.2)
    }
  } else {
    p <- p +
      geom_line(aes(color = group), alpha = alpha) +
      geom_point(aes(color = group), alpha = alpha, size = 0.6)
    if (smooth) {
      p <- p + geom_smooth(aes(group = group, color = group, fill = group), method = "loess",
                           formula = y ~ x, linewidth = 1.2)
    }
    p <- p + labs(color = event_type_variable, fill = event_type_variable)
  }

  if (align == "event") {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed") +
      xlab(sprintf("Time before event (%s - %s)", survival_variable, time_variable))
  } else {
    p <- p + xlab(time_variable)
  }
  p + ylab(bio_variable) + theme_bw()
}
