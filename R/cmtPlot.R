#' Plot conditional mean trajectories (CMT)
#'
#' @description This function generates the Conditional Mean Trajectories (CMT)
#' plot: the mean biomarker trajectory of the patients whose event happened in
#' the window \code{(condi_time2event - event_window, condi_time2event]}. The
#' time from the first measurement up to \code{condi_time2event} is cut into
#' consecutive intervals of length \code{interval_time}; each point is the mean
#' of those patients' biomarker values measured within one interval, plotted at
#' the interval's end, and a straight line is fitted through the points. With
#' competing risks there is one trajectory per event type.
#'
#' All times -- \code{time_variable}, \code{survival_variable},
#' \code{condi_time2event}, \code{event_window} and \code{interval_time} -- are
#' in the units of the data (years, months, days, ...), and the time axis shows
#' \code{time_variable} on that scale, labelled with its name.
#'
#' Only patients with an observed event enter: with competing risks, those
#' whose \code{event_type_variable} equals one of \code{event_type}; without,
#' those whose \code{status_variable} is not 0. Without competing risks and
#' with \code{status_variable = NULL}, patients are selected by their survival
#' time alone, so censored ones are included too, with a warning. An interval
#' without any measurement is left out rather than filled in.
#'
#' @param data_plot_all A \code{data.frame} that includes the biomarker used for plotting.
#' It is utilized to plot conditional mean trajectories (CMT).
#' @param condi_time2event Event time the patients are selected around. If
#' \code{NULL}, it defaults to the midpoint of the observed range of
#' \code{time_variable} in \code{data_plot_all}.
#' @param event_type_variable Competing risks variable indicator name. Set to NULL if
#' there are no competing risks.
#' @param event_type A vector containing the values of \code{event_type_variable}
#' for the event types to plot (at least two).
#' @param bio_variable Name of the biomarker variable used for plotting.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survival_variable Name of the time-to-event outcomes variable.
#' @param interval_time Length of the intervals the follow-up is cut into, in the
#' units of \code{time_variable}.
#' @param id_variable Name of the patient ID column in \code{data_plot_all}. Default is \code{"id"}.
#' @param status_variable Without competing risks, the name of the event
#' indicator (0 = censored), so that only patients with an event are used.
#' Ignored with competing risks.
#' @param event_window Width of the window of event times before
#' \code{condi_time2event} that a patient's event must fall in, in the units of
#' \code{survival_variable}. Default 1 (one year for data in years; use e.g. 12
#' for months or 365 for days).
#'
#' @return Conditional mean trajectories plot.
#'
#' @examples
#'
#' # example without competing risks
#'
#' data(pbc3)
#'
#' pbc.cmt <- cmtPlot(data_plot_all = pbc3, condi_time2event = 5,
#'    event_type_variable = NULL, event_type = NULL,
#'    bio_variable = "serBilir", time_variable = "year",
#'    survival_variable = "years",
#'    interval_time = 1/12, status_variable = "status3"
#' )
#'
#' pbc.cmt
#'
#' @examples
#'
#' # example with competing risks
#'
#' data(pbc3)
#'
#' pbc.cmt.cr <- cmtPlot(pbc3, condi_time2event = 5,
#'    event_type_variable = 'status4', event_type = c("0", "1"),
#'    bio_variable = "albumin", time_variable = "year",
#'    survival_variable = "years",
#'    interval_time = 1/4
#' )
#'
#' pbc.cmt.cr
#'
#' @export
cmtPlot = function(data_plot_all, condi_time2event, event_type_variable, event_type,
                   bio_variable, time_variable, survival_variable, interval_time = 1/12,
                   id_variable = "id", status_variable = NULL, event_window = 1){

  assert_data_frame(data_plot_all, "data_plot_all")
  assert_string(bio_variable, "bio_variable")
  assert_string(time_variable, "time_variable")
  assert_string(survival_variable, "survival_variable")
  assert_string(id_variable, "id_variable")
  assert_vars_in_data(c(bio_variable, time_variable, survival_variable, id_variable), data_plot_all,
                       "bio_variable/time_variable/survival_variable/id_variable", "data_plot_all")
  has_cr <- !is.null(event_type_variable)
  if (has_cr) {
    assert_string(event_type_variable, "event_type_variable")
    assert_vars_in_data(event_type_variable, data_plot_all, "event_type_variable", "data_plot_all")
    if (length(event_type) < 2) {
      stop("`event_type` must have at least 2 elements (one per competing event type) when `event_type_variable` is supplied.", call. = FALSE)
    }
  } else if (!is.null(status_variable)) {
    assert_string(status_variable, "status_variable")
    assert_vars_in_data(status_variable, data_plot_all, "status_variable", "data_plot_all")
  }
  if (!is.null(condi_time2event)) {
    assert_scalar_numeric(condi_time2event, "condi_time2event")
  }
  assert_scalar_numeric(interval_time, "interval_time", positive = TRUE)
  assert_scalar_numeric(event_window, "event_window", positive = TRUE)

  times <- data_plot_all[[time_variable]]
  if (is.null(condi_time2event)) {
    ### the midpoint of the observed follow-up (it used to be rounded up to a
    ### whole unit, which assumed the data were in years)
    condi_time2event <- (min(times, na.rm = TRUE) + max(times, na.rm = TRUE)) / 2
  }

  ### the patients whose event is in (condi_time2event - event_window,
  ### condi_time2event]; the window used to be a fixed 1 unit of time
  surv <- data_plot_all[[survival_variable]]
  in_window <- !is.na(surv) & surv > condi_time2event - event_window & surv <= condi_time2event
  groups <- if (has_cr) {
    lapply(event_type, function(k) in_window & as.character(data_plot_all[[event_type_variable]]) %in% as.character(k))
  } else if (!is.null(status_variable)) {
    status <- data_plot_all[[status_variable]]
    list(in_window & !is.na(status) & status != 0)
  } else {
    warning(paste0(
      "`status_variable` is not given, so patients are selected by their survival time alone and ",
      "censored ones are included; pass the event indicator (e.g. status_variable = \"status\") ",
      "to use only patients with an event."), call. = FALSE)
    list(in_window)
  }
  labels <- if (has_cr) paste0("Event", seq_along(event_type)) else "Event1"

  ### mean biomarker value in each interval (t - interval_time, t] up to
  ### condi_time2event, plotted at the interval's end. An interval without a
  ### measurement is left out: it used to be filled with the previous mean,
  ### and those copies were fitted by the line as if they were data.
  cmt_means <- function(rows) {
    d <- data_plot_all[rows & !is.na(data_plot_all[[bio_variable]]) & !is.na(times), , drop = FALSE]
    if (nrow(d) == 0) return(data.frame(time = numeric(0), probEvent = numeric(0)))
    ends <- seq(min(d[[time_variable]]), condi_time2event, by = interval_time)
    if (length(ends) == 0) ends <- condi_time2event
    if (max(d[[time_variable]]) > max(ends)) ends <- c(ends, max(ends) + interval_time)
    bin <- findInterval(d[[time_variable]], c(-Inf, ends), left.open = TRUE)
    m <- tapply(d[[bio_variable]], factor(bin, levels = seq_along(ends)), mean)
    out <- data.frame(time = ends, probEvent = as.vector(m))
    out[!is.na(out$probEvent), , drop = FALSE]
  }
  cmt_data <- lapply(groups, cmt_means)
  n_patients <- vapply(groups, function(g) length(unique(data_plot_all[[id_variable]][g])), integer(1))
  if (all(n_patients == 0)) {
    stop(sprintf("No patient has an event between %g and %g (condi_time2event - event_window and condi_time2event).",
                 condi_time2event - event_window, condi_time2event), call. = FALSE)
  }

  colours <- c("Event1" = "red", "Event2" = "black")
  Plot_p1 <- ggplot()
  for (k in seq_along(cmt_data)) {
    if (nrow(cmt_data[[k]]) == 0) next
    dk <- cmt_data[[k]]
    dk$line <- labels[k]
    Plot_p1 <- Plot_p1 +
      geom_point(data = dk, aes(x = time, y = probEvent, color = line), size = 3) +
      geom_smooth(data = dk, aes(x = time, y = probEvent, color = line), method = "lm", formula = y ~ x)
  }
  Plot_p1 +
    scale_color_manual(name = "Lines", values = colours[labels]) +
    ylab("Biomarker value") + xlab(time_variable) +
    theme_bw(base_size = 22) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position = if (has_cr) "right" else "none")
}
