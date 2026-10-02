#' Plot cumulative incidence functions
#'
#' @description
#' Plots the nonparametric cumulative incidence of each event type over
#' time, estimated by the Aalen--Johansen estimator (via
#' \code{\link[survival]{survfit}} on a multi-state \code{Surv} outcome).
#' With a single event type this equals one minus the Kaplan--Meier
#' survival curve; with competing risks each event type gets its own curve,
#' correctly accounting for the other event types instead of treating them
#' as censoring. An optional \code{group_variable} draws separate curves per
#' group (e.g. treatment arm).
#'
#' @param data_plot_all A \code{data.frame} with the time-to-event outcome.
#' It may be in long format (several rows per subject, as for the
#' longitudinal sub-model); only the first row of each subject is used, and
#' \code{survival_variable}, \code{event_type_variable} and
#' \code{group_variable} must be constant within each subject.
#' @param survival_variable Name of the time-to-event variable.
#' @param event_type_variable Name of the variable holding each subject's
#' event type, or censoring status.
#' @param id_variable Name of the subject ID column. Default is \code{"id"}.
#' @param censor_value The value of \code{event_type_variable} that marks a
#' censored subject. Default \code{0}. Use \code{NA} when censored subjects
#' have a missing event type, as with \code{pbc3$status4}. Every other value
#' is treated as an event type.
#' @param event_labels Optional legend labels for the event types: either a
#' character vector named by the event type values (e.g.
#' \code{c("1" = "Death", "2" = "Transplant")}), or an unnamed one in the
#' sorted order of the event type values. Default uses the values themselves.
#' @param group_variable Optional name of a subject-level variable to stratify
#' the curves by, shown by line type (and also by color when there is a
#' single event type). \code{NULL} (default) for no grouping.
#' @param conf_int Logical; draw pointwise 95\% confidence bands. Default
#' \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(pbc3)
#'
#' # competing risks: death (1) and transplantation (2); 0 = censored
#' cifPlot(pbc3, survival_variable = "years", event_type_variable = "status5",
#'   event_labels = c("1" = "Death", "2" = "Transplant"))
#'
#' # the same events coded as in status4 (censored subjects are NA)
#' cifPlot(pbc3, survival_variable = "years", event_type_variable = "status4",
#'   censor_value = NA, event_labels = c("0" = "Death", "1" = "Transplant"))
#'
#' # single event type, by treatment arm
#' cifPlot(pbc3, survival_variable = "years", event_type_variable = "status2",
#'   group_variable = "drug", conf_int = FALSE)
#'
#' @export
cifPlot <- function(data_plot_all, survival_variable, event_type_variable, id_variable = "id",
                    censor_value = 0, event_labels = NULL, group_variable = NULL,
                    conf_int = TRUE) {

  assert_data_frame(data_plot_all, "data_plot_all")
  assert_string(survival_variable, "survival_variable")
  assert_string(event_type_variable, "event_type_variable")
  assert_string(id_variable, "id_variable")
  vars <- c(survival_variable, event_type_variable, id_variable)
  if (!is.null(group_variable)) {
    assert_string(group_variable, "group_variable")
    vars <- c(vars, group_variable)
  }
  assert_vars_in_data(vars, data_plot_all, "cifPlot() variables", "data_plot_all")
  if (length(censor_value) != 1) {
    stop("`censor_value` must be a single value (or NA).", call. = FALSE)
  }

  ### one row per subject; the outcome must not vary within a subject
  id_chr <- as.character(data_plot_all[[id_variable]])
  for (v in setdiff(vars, id_variable)) {
    n_values <- tapply(data_plot_all[[v]], id_chr, function(x) length(unique(x)))
    if (any(n_values > 1)) {
      stop(sprintf("`%s` must be constant within each subject; it varies for %d subject(s).",
                   v, sum(n_values > 1)), call. = FALSE)
    }
  }
  d <- data_plot_all[!duplicated(id_chr), , drop = FALSE]

  status <- d[[event_type_variable]]
  if (is.na(censor_value)) {
    censored <- is.na(status)
  } else {
    if (anyNA(status)) {
      warning(sprintf("Dropping %d subject(s) with a missing `%s`; set `censor_value = NA` if NA means censored.",
                      sum(is.na(status)), event_type_variable), call. = FALSE)
      d <- d[!is.na(status), , drop = FALSE]
      status <- status[!is.na(status)]
    }
    censored <- as.character(status) == as.character(censor_value)
  }
  keep <- !is.na(d[[survival_variable]])
  if (!is.null(group_variable)) {
    keep <- keep & !is.na(d[[group_variable]])
  }
  d <- d[keep, , drop = FALSE]
  status <- status[keep]
  censored <- censored[keep]
  if (nrow(d) == 0) {
    stop("`data_plot_all` has no subjects with a non-missing outcome.", call. = FALSE)
  }

  event_values <- sort(unique(status[!censored]))
  if (length(event_values) == 0) {
    stop(sprintf("No events found: every subject is censored according to `censor_value` (%s).",
                 format(censor_value)), call. = FALSE)
  }
  event_values <- as.character(event_values)
  event_names <- cif_event_labels(event_values, event_labels)

  state <- factor(ifelse(censored, "(censored)", as.character(status)),
                  levels = c("(censored)", event_values))
  groups <- if (is.null(group_variable)) factor(rep("All", nrow(d))) else factor(d[[group_variable]])

  curves <- list()
  for (g in levels(groups)) {
    in_g <- groups == g
    fit <- survfit(Surv(d[[survival_variable]][in_g], state[in_g]) ~ 1)
    for (k in seq_along(event_values)) {
      col <- match(event_values[k], fit$states)
      curves[[length(curves) + 1]] <- data.frame(
        time = c(0, fit$time),
        cif = c(0, fit$pstate[, col]),
        lower = c(0, if (is.null(fit$lower)) NA else fit$lower[, col]),
        upper = c(0, if (is.null(fit$upper)) NA else fit$upper[, col]),
        event = event_names[k],
        group = g
      )
    }
  }
  curves <- do.call(rbind, curves)
  ### before an event type's first event its estimate is exactly 0 and
  ### survfit() gives no limits; use the estimate so the band stays a step
  no_limits <- is.na(curves$lower) | is.na(curves$upper)
  curves$lower[no_limits & curves$cif == 0] <- 0
  curves$upper[no_limits & curves$cif == 0] <- 0
  curves$event <- factor(curves$event, levels = event_names)
  curves$group <- factor(curves$group, levels = levels(groups))
  curves$curve <- interaction(curves$event, curves$group)

  ### color by event type; with one event type and a grouping variable,
  ### color by group instead so the groups' curves and bands stand apart
  color_by_group <- length(event_values) == 1 && !is.null(group_variable)
  curves$hue <- if (color_by_group) curves$group else curves$event
  hue_title <- if (color_by_group) group_variable else "Event type"

  p <- ggplot(curves, aes(x = time, y = cif, color = hue, group = curve))
  if (conf_int) {
    p <- p + geom_ribbon(data = cif_step_ribbon(curves),
                         aes(x = time, ymin = lower, ymax = upper, fill = hue, group = curve),
                         inherit.aes = FALSE, alpha = 0.15)
  }
  if (is.null(group_variable)) {
    p <- p + geom_step(linewidth = 1)
  } else {
    p <- p + geom_step(aes(linetype = group), linewidth = 1) + labs(linetype = group_variable)
  }
  if (length(event_values) == 1 && is.null(group_variable) && is.null(event_labels)) {
    ### a single, unlabeled event type: the color legend adds nothing
    p <- p + guides(color = "none", fill = "none")
  }
  p + labs(color = hue_title, fill = hue_title) +
    scale_y_continuous(limits = c(0, 1)) +
    xlab(survival_variable) + ylab("Cumulative incidence") +
    theme_bw()
}

#' Legend labels for the event types in cifPlot()
#'
#' @param event_values Sorted event type values, as character.
#' @param event_labels The user's \code{event_labels} argument.
#' @return A character vector of labels, one per event type.
#' @keywords internal
cif_event_labels <- function(event_values, event_labels) {
  if (is.null(event_labels)) {
    return(event_values)
  }
  if (!is.character(event_labels)) {
    stop("`event_labels` must be a character vector.", call. = FALSE)
  }
  if (is.null(names(event_labels))) {
    if (length(event_labels) != length(event_values)) {
      stop(sprintf("An unnamed `event_labels` must have one label per event type (%d: %s).",
                   length(event_values), paste(event_values, collapse = ", ")), call. = FALSE)
    }
    return(unname(event_labels))
  }
  missing_values <- setdiff(event_values, names(event_labels))
  if (length(missing_values) > 0) {
    stop(sprintf("`event_labels` has no label for event type(s): %s.",
                 paste(missing_values, collapse = ", ")), call. = FALSE)
  }
  unname(event_labels[event_values])
}

#' Expand cumulative incidence confidence limits into a step-shaped ribbon
#'
#' @description \code{geom_ribbon()} joins points with straight lines; to
#' follow the step curves, each jump time gets two rows, one with the
#' previous limits and one with the new limits.
#' @param curves The curve \code{data.frame} built in \code{cifPlot()}.
#' @return A \code{data.frame} with columns \code{time}, \code{lower},
#' \code{upper}, \code{hue}, \code{curve}.
#' @keywords internal
cif_step_ribbon <- function(curves) {
  pieces <- lapply(split(curves, curves$curve, drop = TRUE), function(cv) {
    n <- nrow(cv)
    idx_time <- c(1, rep(seq_len(n)[-1], each = 2))
    idx_val <- c(1, as.vector(rbind(seq_len(n - 1), seq_len(n)[-1])))
    data.frame(time = cv$time[idx_time], lower = cv$lower[idx_val], upper = cv$upper[idx_val],
               hue = cv$hue[1], curve = cv$curve[1])
  })
  out <- do.call(rbind, pieces)
  out[!is.na(out$lower) & !is.na(out$upper), , drop = FALSE]
}
