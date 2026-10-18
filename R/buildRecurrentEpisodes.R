#' Reshape recurrent-event data into gap-time episodes for the backward joint model
#'
#' @description
#' \strong{Experimental.} Prepares recurrent-event data for the backward joint
#' model by a gap-time (renewal) reformulation: it splits each subject's
#' event history into \emph{episodes} -- the interval from one restart (the
#' previous event, or the discharge from the stay it began) to the next event,
#' with the clock reset to \code{0} at the start of each -- and re-clocks any
#' accompanying longitudinal measurements onto the same per-episode timeline.
#' The result plugs into the existing \code{survivalSub()}/
#' \code{longitudinalSub()}/\code{predictRisk()} pipeline, unchanged, to
#' predict for a subject currently at risk in an open episode whether their
#' \emph{next} event (within a horizon) is a recurrence or a terminal event;
#' \code{\link{predictRecurrentCount}} builds on the same data to predict the
#' \emph{number} of recurrences in a window. Unlike the single-event and
#' competing-risks methods the rest of the package rests on (Shen and Li,
#' 2021; Li et al., 2026 -- see \code{citation("BJM")}), this extension is
#' unpublished, and has been checked only in simulation studies, not on real
#' data (see the "Recurrent events" article on the package website).
#'
#' @details
#' \strong{What counts as an episode.} Episode 1 runs from time \code{0} to
#' a subject's first row in \code{data_events}; episode \code{k > 1} runs
#' from the end of episode \code{k - 1} to that subject's \code{k}-th row.
#' A subject contributes no further episodes after a row with
#' \code{status \%in\% c(0, 2)} (administrative censoring or a terminal
#' event); any later rows for that \code{id} are ignored, with a message
#' saying how many. \code{status == 1}
#' rows (recurrences) always start a new episode.
#'
#' \strong{Feeding the result into \code{survivalSub()}.} Fit the gap-time
#' recurrence hazard and the recurrence-vs-terminal-event type model exactly
#' as \code{survivalSub()} already fits a marginal hazard and a competing-risks
#' type model, using \code{episode_data} in place of one-row-per-subject
#' data, \code{Surv(gap, event_any) ~ covariates} in place of
#' \code{form_marginal_surv}, and \code{event_type ~ gap + covariates} in
#' place of \code{form_conditional_cr} (\code{event_type} is \code{0} for a
#' recurrence, \code{1} for a terminal event -- so \code{predictRisk()}'s
#' \code{risk_prob_1}/\code{risk_prob_2} come out as the recurrence/terminal
#' risk respectively). Administratively-censored episodes get
#' \code{event_any = 0} and \code{event_type = NA}, exactly like
#' \code{status3 = 0}/\code{status4 = NA} for a censored subject in
#' \code{\link{pbc3}}. The episodes of one subject are then treated as
#' independent: the coefficients are still consistent for the gap-time
#' model, but the standard errors that \code{summary()} reports ignore the
#' correlation between a subject's episodes and are too small.
#'
#' \strong{Feeding the result into \code{longitudinalSub()}/\code{predictRisk()}.}
#' \code{gap} plays the same dual role \code{years} plays for \link{pbc3}: it
#' is both the \code{Surv()} time in \code{episode_data} and the backward
#' covariate swept over candidate values inside \code{predictRisk()}, so it
#' must appear in \code{long_sub_fixed} (e.g. \code{biomarker ~
#' time_in_episode + gap + gap * time_in_episode + covariates}), and the
#' training data passed to \code{longitudinalSub()} must be filtered to
#' \code{event_any == 1} rows of \code{long_data} -- only completed episodes
#' have an observed \code{gap} -- exactly as the package's own examples
#' filter \code{pbc3} to \code{status3 == 1} before fitting
#' \code{longitudinalSub()}. \code{time_in_episode} (elapsed time since the
#' current episode began) plays the role \code{year} plays for
#' \code{pbc3}. To predict from "now" for a subject currently at risk in an
#' open episode, add a trailing \code{status = 0} row to \code{data_events}
#' at the current data-cutoff time, and call \code{predictRisk()} with
#' \code{prediction_time} set to that subject's elapsed
#' \code{time_in_episode} at the cutoff.
#'
#' \strong{Random effects: group by episode, not by patient.} Fit the biomarker
#' model with the random effects grouped by \code{episode_id}
#' (\code{random = ~ time_in_episode | episode_id}), so that each episode is
#' its own "subject", as it already is in the survival sub-model. A patient's
#' biomarker level generally differs from one episode to the next (the
#' underlying trajectory keeps moving), so a single random effect per patient
#' is the wrong unit; in a simulation in which the biomarker's current level
#' drives the risk, grouping by patient gave a clearly worse AUC than grouping
#' by episode. \code{predictRisk()} and \code{predictRecurrentCount()} then
#' take one \code{episode_id} per prediction, so give each prediction its own
#' \code{id}.
#'
#' \strong{The open episode's \code{gap} is \code{NA}.} For an episode that
#' has not ended (\code{event_any == 0}) -- typically the one a subject is
#' currently at risk in -- the \code{gap} column of \code{long_data} is
#' \code{NA}: its final length is exactly what is being predicted.
#' \code{predictRisk()} reads a missing survival-time variable as "event time
#' not yet known, still at risk", so that episode's rows can be passed
#' straight to it, with no placeholder to set or override (\code{episode_data}
#' keeps the real elapsed time of such an episode, which is what the Cox fit
#' needs as a censoring time). Training rows are unaffected, since
#' \code{longitudinalSub()} is fit on \code{event_any == 1} rows only. Note that
#' \code{predictRisk()} also expects the outcome columns it was fit with
#' (the time and event-type variables) to exist in the prediction data;
#' \code{long_data} already carries \code{event_type} (\code{NA} for an open
#' episode).
#'
#' \strong{Episodes that start at discharge (\code{restart_variable},
#' \code{lookback}).} By default an episode starts at the previous event
#' (time \code{0} for the first). For hospitalization data that is not the
#' at-risk clock: after a readmission the patient is in hospital, not at risk,
#' until discharge. Give \code{data_events} a numeric column -- passed as
#' \code{restart_variable} -- holding, for each recurrence row, the time the
#' process restarts (the discharge ending that stay); the next episode then
#' starts there and its gap is measured from it. Time \code{0} is the first
#' episode's start (e.g. the index discharge). Biomarkers measured during the
#' stay, or before the index discharge, precede the episode they predict:
#' \code{lookback} (default \code{0}, none) keeps measurements up to that
#' long before an episode's start -- never reaching back past the previous
#' event, so a stay's own labs are not mixed into the earlier at-risk period --
#' and gives them \emph{negative} \code{time_in_episode}. Use \code{Inf} to keep
#' the whole stay (and, for the first episode, everything before time
#' \code{0}). Negative times are fine for \code{longitudinalSub()} and
#' \code{predictRisk()}.
#'
#' \strong{Death, or the end of follow-up, during a hospital stay.} With
#' \code{restart_variable}, a terminal event or censoring row whose time is at
#' or before the restart time of the readmission before it happened while the
#' patient was in hospital, not at risk. That readmission then ends the
#' subject's process: no further episode is opened, the row is dropped with a
#' message, and so a death in hospital is not counted as an out-of-hospital
#' terminal event (\code{predictRecurrentCount()}'s \code{p_death} is the
#' probability of dying after discharge). Record such a stay's restart time
#' as the time of death or of the end of follow-up, or later, when the
#' discharge never happened. A data cutoff during a stay gives no open
#' episode to predict for.
#'
#' \strong{Covariates that change over time.} Each episode's covariates come
#' from the row of \code{data_events} that ends it, as in counting-process
#' data, where each interval carries the covariate values in force during it.
#' A covariate that changes over time must therefore be recorded on each row
#' as its value at the start of that row's episode (age at the discharge
#' that started it, the number of earlier admissions, the length of the
#' previous stay), and the open episode's row in the same way. A value
#' measured at the event itself -- age at readmission, say -- depends on how
#' long the episode lasted, which is what the survival sub-model explains,
#' and so biases it; a value measured at the time of the prediction is not
#' what the model was fit with. Biomarkers that change within an episode
#' belong in \code{data_long}. Covariates built from earlier episodes can
#' also be added to \code{episode_data} afterwards, e.g. the length of the
#' previous gap:
#' \preformatted{ed <- episodes$episode_data
#' ed$previous_gap <- ave(ed$gap, ed$id, FUN = function(g) c(0, g[-length(g)]))}
#'
#' @param data_events One row per subject per event-or-censoring occurrence,
#'   in calendar time (i.e. not yet re-clocked to gap time -- this function
#'   does that). Must contain \code{id_variable}, \code{time_variable}
#'   (numeric, strictly increasing within each subject), and
#'   \code{status_variable} (\code{0} = censored/end of follow-up, \code{1} =
#'   recurrent event, \code{2} = terminal event -- the same 0/1/2 convention
#'   \link{pbc3}'s \code{status5} already uses). Any other columns are
#'   covariates, carried through to \code{episode_data} from the row that
#'   ends each episode: so a covariate that changes over time must hold, on
#'   each row, its value at the \emph{start} of the episode that row ends
#'   (see Details), not its value at the event.
#' Subjects are matched across \code{data_events} and \code{data_long} by the
#' printed value of their id, so an id that is a number in one table and a
#' factor or character in the other still matches; \code{id} keeps the type
#' it has in \code{data_events}.
#'
#' @param data_long Optional repeated-measurement data to re-clock onto the
#'   same per-episode timeline, in the same format \code{longitudinalSub()}/
#'   \code{predictRisk()} accept elsewhere in this package: a single bare
#'   \code{data.frame}, or a list of \code{data.frame}s (one per
#'   longitudinal outcome). Must contain \code{id_variable} and
#'   \code{time_variable_long}. Rows whose \code{time_variable_long} falls
#'   outside every episode -- before \code{0} (beyond \code{lookback}), in a
#'   hospital stay \code{lookback} does not reach, or after the subject's
#'   process has already ended (per \code{data_events}) -- are dropped, with
#'   a message reporting how many. A
#'   row whose \code{time_variable_long} falls exactly on the shared
#'   boundary between two episodes (e.g. a biomarker drawn at the same
#'   visit an event is recorded) is genuinely ambiguous -- it is at once
#'   "the last measurement of the ending episode" and "the first
#'   (\code{time_in_episode = 0}) measurement of the next one" -- so it is
#'   kept \emph{twice}, once per episode, rather than silently assigned to
#'   only one (which would systematically discard every fresh-episode
#'   baseline measurement that happens to coincide with the previous
#'   episode's event time, and bias the fitted backward-covariate effect); a
#'   message says how many rows this applies to.
#' @param id_variable,time_variable,status_variable Column names in
#'   \code{data_events} (defaults \code{"id"}, \code{"time"}, \code{"status"}).
#' @param time_variable_long Column name of the measurement time in
#'   \code{data_long}; defaults to \code{time_variable} (i.e. the same name
#'   as in \code{data_events}) when \code{NULL}.
#' @param episode_variable,gap_variable,time_in_episode_variable,event_any_variable,event_type_variable,episode_id_variable
#'   Column names to create in the output (defaults \code{"episode"},
#'   \code{"gap"}, \code{"time_in_episode"}, \code{"event_any"},
#'   \code{"event_type"}, \code{"episode_id"}). \code{episode_id} is
#'   \code{"<id>:<episode>"} (a factor), one value per episode: the grouping variable
#'   to use for the biomarker model's random effects (see Details). Change these if they collide with an existing
#'   column of \code{data_events}/\code{data_long} -- collisions are an
#'   error, not a silent overwrite.
#' @param restart_variable Optional name of a numeric column of
#'   \code{data_events}: for each recurrence (\code{status == 1}) row, the time
#'   the process restarts after it (for hospitalizations, the discharge ending
#'   that stay; must not be earlier than the row's \code{time_variable}). The
#'   next episode starts there instead of at the event time. \code{NULL}
#'   (default) restarts at the event time itself. A missing value on a
#'   recurrence row that starts another episode is taken as no stay, with a
#'   warning; a column with no value at all is allowed. A death or the end of
#'   follow-up at or before the restart time happened during the stay: see
#'   Details. Ignored on rows that end the process. Not carried into
#'   \code{episode_data}.
#' @param lookback A single non-negative number (\code{Inf} allowed): how far
#'   before an episode's start its measurements are still kept, as negative
#'   \code{time_in_episode}, but never from before the previous event. Default
#'   \code{0} keeps only measurements from the start on.
#'
#' @return A list with elements:
#' \describe{
#'   \item{episode_data}{One row per subject per episode:
#'   \code{id_variable}, \code{episode_variable} (episode number, starting
#'   at 1), \code{gap_variable} (that episode's gap time),
#'   \code{event_any_variable} (\code{0}/\code{1}: did the episode end in an
#'   event of either type, i.e. is \code{gap} an observed event time rather
#'   than an administrative censoring time), \code{event_type_variable}
#'   (\code{0} = recurrence, \code{1} = terminal event, \code{NA} when
#'   \code{event_any == 0}), and every other column of \code{data_events}.}
#'   \item{long_data}{\code{NULL} if \code{data_long} was not supplied;
#'   otherwise \code{data_long}, in the same bare-\code{data.frame}-or-list
#'   shape it was supplied in, with \code{episode_variable},
#'   \code{gap_variable}, \code{time_in_episode_variable},
#'   \code{event_any_variable}, and \code{event_type_variable} columns added
#'   (and out-of-range rows dropped, as described above). Any other column
#'   of \code{episode_data} (baseline or event-time-varying covariates, e.g.
#'   \code{age}/\code{sex}) is also carried onto matching rows here, the way
#'   every row of \link{pbc3}'s long-format table already repeats
#'   \code{age}/\code{sex}/\code{years}/\code{status3} -- but only into a
#'   column \code{data_long} does not already have of that name.}
#' }
#'
#' @examples
#' ## A minimal example: reshaping alone (no model fitting).
#' ## Subject 1: two recurrences, then a terminal event; subject 2: censored.
#' data_events <- data.frame(
#'   id     = c(1, 1, 1, 2),
#'   time   = c(1.2, 2.5, 4.0, 3.1),
#'   status = c(1,   1,   2,   0),
#'   age    = c(58,  58,  58,  61)
#' )
#' episodes <- buildRecurrentEpisodes(data_events)
#' episodes$episode_data
#'
#' ## Hospitalizations: each episode starts at DISCHARGE, not at the readmission.
#' ## `time` = readmission (or death/censoring) day; `discharge` = day the stay
#' ## it started ends. Labs drawn in the stay (and before the index discharge at
#' ## day 0) are kept with negative time_in_episode via lookback = Inf.
#' events <- data.frame(id = 1, time = c(20, 60, 90), status = c(1, 1, 0),
#'                      discharge = c(26, 66, NA))
#' labs <- data.frame(id = 1, day = c(-3, 0, 5, 21, 26, 40, 61, 66, 80))
#' hf <- buildRecurrentEpisodes(events, labs, time_variable_long = "day",
#'                              restart_variable = "discharge", lookback = Inf)
#' hf$episode_data   # gaps measured from each discharge: 20, 34, 24
#' hf$long_data[, c("day", "episode", "time_in_episode", "gap")]
#'
#' \donttest{
#' ## A full pipeline on simulated recurrent-event + longitudinal data:
#' ## build episodes -> survivalSub() -> longitudinalSub() -> predictRisk() and
#' ## predictRecurrentCount().
#' set.seed(1)
#' n <- 60
#' data_events <- do.call(rbind, lapply(seq_len(n), function(id) {
#'   k <- sample(1:3, 1)
#'   gaps <- rexp(k, rate = 1)
#'   data.frame(
#'     id = id, time = cumsum(gaps),
#'     status = c(rep(1, k - 1), sample(c(0, 2), 1, prob = c(0.6, 0.4))),
#'     age = round(50 + stats::rnorm(1, sd = 8)), sex = id %% 2
#'   )
#' }))
#'
#' ## subject 9999 is still at risk right now, 0.7 time units into their
#' ## 2nd episode (1 prior recurrence, no resolution yet) -- this is who
#' ## we will predict for below.
#' data_events <- rbind(data_events, data.frame(
#'   id = 9999, time = c(0.8, 1.6), status = c(1, 0), age = 55, sex = 1
#' ))
#'
#' data_long <- do.call(rbind, lapply(c(seq_len(n), 9999), function(id) {
#'   maxt <- max(data_events$time[data_events$id == id]) + 0.2
#'   tt <- seq(0, maxt, by = 0.25)
#'   data.frame(id = id, year = tt, biomarker = rnorm(length(tt)))
#' }))
#'
#' built <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year")
#'
#' ## gap-time recurrence hazard + recurrence-vs-terminal-event type model
#' survival_fit_recur <- survivalSub(
#'   built$episode_data,
#'   Surv(gap, event_any) ~ age + sex + episode,
#'   event_type ~ gap + age + sex + episode
#' )
#'
#' ## longitudinal sub-model, trained on completed episodes only, with the
#' ## random effects grouped by episode (see Details)
#' long_fit_recur <- longitudinalSub(
#'   built$long_data[built$long_data$event_any == 1, ],
#'   biomarker ~ time_in_episode + age + sex + gap + gap * time_in_episode,
#'   ~ time_in_episode | episode_id
#' )
#'
#' ## subject 9999's open episode: its gap is NA, which predictRisk() reads as
#' ## "still at risk" -- nothing to override
#' open <- built$long_data[built$long_data$id == 9999 & built$long_data$event_any == 0, ]
#' landmark <- max(open$time_in_episode)
#'
#' ## the risk that their *next* event, within 1 time unit, is a recurrence /
#' ## a terminal event
#' risk <- predictRisk(
#'   open, long_fit_recur, survival_fit_recur,
#'   prediction_time = landmark, horizon = 1, time_variable = "time_in_episode",
#'   survival_variable_all = list(), survival_trans_function = list(),
#'   bandcount1 = 10, bandcount2 = 10
#' )
#' risk$risk_prob_1 # P(next event within horizon is a recurrence)
#' risk$risk_prob_2 # P(next event within horizon is the terminal event)
#'
#' ## the number of recurrences over the next 3 time units
#' predictRecurrentCount(
#'   open, long_fit_recur, survival_fit_recur,
#'   prediction_time = landmark, horizon = 3, time_variable = "time_in_episode",
#'   n_sim = 100, seed = 1
#' )
#' }
#'
#' @export
buildRecurrentEpisodes <- function(data_events,
                                    data_long = NULL,
                                    id_variable = "id",
                                    time_variable = "time",
                                    status_variable = "status",
                                    time_variable_long = NULL,
                                    episode_variable = "episode",
                                    gap_variable = "gap",
                                    time_in_episode_variable = "time_in_episode",
                                    event_any_variable = "event_any",
                                    event_type_variable = "event_type",
                                    restart_variable = NULL,
                                    lookback = 0,
                                    episode_id_variable = "episode_id") {

  assert_data_frame(data_events, "data_events")
  for (nm in c("id_variable", "time_variable", "status_variable", "episode_variable",
               "gap_variable", "time_in_episode_variable", "event_any_variable",
               "event_type_variable", "episode_id_variable")) {
    assert_string(get(nm), nm)
  }
  assert_vars_in_data(c(id_variable, time_variable, status_variable), data_events,
                       "id_variable/time_variable/status_variable", "data_events")

  new_cols <- c(episode_variable, gap_variable, event_any_variable, event_type_variable, episode_id_variable)
  clash <- intersect(new_cols, names(data_events))
  if (length(clash) > 0) {
    stop(sprintf(
      "`data_events` already has column(s) %s, which collide with buildRecurrentEpisodes()'s output column names. Rename those columns first, or pass different *_variable argument(s).",
      paste(sprintf("`%s`", clash), collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.null(restart_variable)) {
    assert_string(restart_variable, "restart_variable")
    assert_vars_in_data(restart_variable, data_events, "restart_variable", "data_events")
    ### a column with no restart time at all is logical NA in R
    if (is.logical(data_events[[restart_variable]]) && all(is.na(data_events[[restart_variable]]))) {
      data_events[[restart_variable]] <- as.numeric(data_events[[restart_variable]])
    }
    if (!is.numeric(data_events[[restart_variable]])) {
      stop(sprintf("`%s` must be a numeric column.", restart_variable), call. = FALSE)
    }
  }
  if (!is.numeric(lookback) || length(lookback) != 1 || is.na(lookback) || lookback < 0) {
    stop("`lookback` must be a single non-negative number (`Inf` allowed).", call. = FALSE)
  }

  status_vals <- data_events[[status_variable]]
  if (anyNA(status_vals) || !all(status_vals %in% c(0, 1, 2))) {
    stop(sprintf(
      "`%s` must only contain 0 (censored, no event), 1 (recurrent event), or 2 (terminal event), with no missing values.",
      status_variable
    ), call. = FALSE)
  }
  if (!is.numeric(data_events[[time_variable]]) || anyNA(data_events[[time_variable]])) {
    stop(sprintf("`%s` must be a numeric column with no missing values.", time_variable), call. = FALSE)
  }
  if (anyNA(data_events[[id_variable]])) {
    stop(sprintf("`%s` must not contain missing values.", id_variable), call. = FALSE)
  }

  ord <- order(data_events[[id_variable]], data_events[[time_variable]])
  data_events <- data_events[ord, , drop = FALSE]
  rownames(data_events) <- NULL

  baseline_cols <- setdiff(names(data_events),
                           c(id_variable, time_variable, status_variable, restart_variable))
  id_key <- id_label(data_events[[id_variable]])
  rows_by_id <- split(seq_len(nrow(data_events)), factor(id_key, levels = unique(id_key)))
  time_all <- data_events[[time_variable]]
  status_all <- as.numeric(as.character(data_events[[status_variable]]))
  restart_all <- if (is.null(restart_variable)) rep(NA_real_, nrow(data_events)) else data_events[[restart_variable]]

  # An episode runs from t_start (time 0 for the first; afterwards the previous
  # event's `restart_variable` if given -- e.g. discharge from the stay a
  # readmission started -- else the previous event time itself) to t_end (that
  # row's time). A patient's episodes stop at the first row that is not a
  # recurrence (censoring or a terminal event).
  pieces <- vector("list", length(rows_by_id))
  n_ignored <- 0L
  missing_restart <- character(0)
  died_in_stay <- ended_in_stay <- character(0)
  for (g in seq_along(rows_by_id)) {
    idx <- rows_by_id[[g]]
    st <- status_all[idx]; tt <- time_all[idx]; rs <- restart_all[idx]
    n_ep <- if (any(st != 1)) which(st != 1)[1] else length(idx)
    n_ignored <- n_ignored + (length(idx) - n_ep)
    ### a death, or the end of follow-up, during the stay a readmission
    ### began (at or before its restart time): the patient was never back at
    ### risk, so the process ends at that readmission (reported below)
    if (!is.null(restart_variable) && n_ep > 1 && st[n_ep] != 1 &&
        !is.na(rs[n_ep - 1]) && tt[n_ep] <= rs[n_ep - 1]) {
      if (st[n_ep] == 2) died_in_stay <- c(died_in_stay, id_key[idx[1]]) else
        ended_in_stay <- c(ended_in_stay, id_key[idx[1]])
      n_ep <- n_ep - 1
    }
    keep <- seq_len(n_ep)
    t_end <- tt[keep]
    is_rec <- st[keep] == 1
    bad <- which(is_rec & !is.na(rs[keep]) & rs[keep] < t_end)
    if (length(bad) > 0) {
      r <- bad[1]
      stop(sprintf("`%s` (%s) is earlier than the event time it belongs to (%s) for id = %s, row %d.",
                   restart_variable, format(rs[r]), format(t_end[r]), id_key[idx[1]], r), call. = FALSE)
    }
    ### a recurrence that starts another episode but has no restart time:
    ### the stay is taken to be empty (reported below)
    if (!is.null(restart_variable) && n_ep > 1) {
      gaps <- which(is_rec[-n_ep] & is.na(rs[keep][-n_ep]))
      missing_restart <- c(missing_restart, rep(id_key[idx[1]], length(gaps)))
    }
    next_start <- ifelse(is.na(rs[keep]), t_end, rs[keep])
    t_start <- c(0, next_start[-n_ep])
    bad <- which(!is.finite(t_end) | t_end <= t_start)
    if (length(bad) > 0) {
      stop(sprintf(
        "`%s` must be strictly increasing within each `%s` -- and later than 0 for the first row, and later than the previous row's `restart_variable` if one is given; id = %s violates this at row %d.",
        time_variable, id_variable, id_key[idx[1]], bad[1]
      ), call. = FALSE)
    }
    pieces[[g]] <- list(idx = idx[keep], ep = keep, t_start = t_start, t_end = t_end,
                        t_prev = c(-Inf, t_end[-n_ep]))
  }
  if (length(missing_restart) > 0) {
    warning(sprintf(paste0(
      "`%s` is missing on %d recurrence row(s) followed by another episode (id %s). The next episode ",
      "then starts at the event time itself, so the time until the missing restart counts as time at ",
      "risk. Fill it in, or set it to the event time for a recurrence without a stay; for a death, or ",
      "the end of follow-up, during the stay, set it to that time or later."),
      restart_variable, length(missing_restart),
      paste(first_n(unique(missing_restart), 5), collapse = ", ")), call. = FALSE)
  }
  if (length(died_in_stay) + length(ended_in_stay) > 0) {
    message(sprintf(paste0(
      "buildRecurrentEpisodes(): %d subject(s) died and %d had their follow-up end during a hospital ",
      "stay (at or before the `%s` of their last readmission). They were not back at risk, so their ",
      "episodes end at that readmission: there is no open episode for them, and a death in hospital ",
      "is not an out-of-hospital terminal event."),
      length(died_in_stay), length(ended_in_stay), restart_variable))
  }
  if (n_ignored > 0) {
    message(sprintf(
      "buildRecurrentEpisodes(): ignored %d row(s) of `data_events` that come after a subject's censoring or terminal event (the recurrent-event process had already ended).",
      n_ignored))
  }
  collect <- function(field) unlist(lapply(pieces, `[[`, field), use.names = FALSE)
  use_idx <- collect("idx"); ep_no <- as.integer(collect("ep"))
  t_start <- collect("t_start"); t_end <- collect("t_end"); t_prev <- collect("t_prev")
  st_used <- status_all[use_idx]

  episode_data <- data_events[use_idx, baseline_cols, drop = FALSE]
  rownames(episode_data) <- NULL
  episode_data[[id_variable]] <- data_events[[id_variable]][use_idx]
  episode_data[[episode_variable]] <- ep_no
  episode_data[[episode_id_variable]] <- factor(paste(id_label(episode_data[[id_variable]]), ep_no, sep = ":"))
  episode_data[[gap_variable]] <- t_end - t_start
  episode_data[[event_any_variable]] <- as.integer(st_used != 0)
  episode_data[[event_type_variable]] <- ifelse(st_used == 1, 0L, ifelse(st_used == 2, 1L, NA_integer_))
  episode_data <- episode_data[, c(id_variable, episode_variable, episode_id_variable, gap_variable,
                                    event_any_variable, event_type_variable, baseline_cols),
                                drop = FALSE]

  long_data <- NULL
  if (!is.null(data_long)) {
    long_is_list <- is.list(data_long) && !is.data.frame(data_long)
    data_long_list <- if (long_is_list) data_long else list(data_long)

    tv_long <- if (is.null(time_variable_long)) time_variable else time_variable_long
    assert_string(tv_long, "time_variable_long")

    # For every episode (aligned with the rows of episode_data): the earliest time
    # whose measurements count toward it (t_lo). With `lookback = 0` (the
    # default) that is the episode's own start; a larger `lookback` also takes in
    # measurements up to that long before it -- but never from before the
    # previous event, which belong to the previous episode's own at-risk period
    # -- giving them negative time_in_episode. Between a readmission and the
    # discharge that restarts the process (the stay) that window is exactly the
    # stay. An id's episodes are contiguous, which the matching below relies on.
    t_lo <- pmax(t_start - lookback, t_prev)
    gap_all <- t_end - t_start
    ep_key <- id_label(episode_data[[id_variable]])
    uid <- unique(ep_key)
    first_pos <- match(uid, ep_key)
    n_ep_by_id <- tabulate(match(ep_key, uid), nbins = length(uid))
    # Baseline/event-time covariates already carried onto episode_data (e.g.
    # age, sex) get carried onto long_data too, exactly as every row of
    # pbc3's long-format table already repeats age/sex/years/status3 --
    # but only into columns data_long does not already have, so a genuinely
    # different variable that happens to share a name is left alone.
    carry_cols <- setdiff(names(episode_data),
                           c(id_variable, episode_variable, gap_variable,
                             event_any_variable, event_type_variable, episode_id_variable))

    augment_one <- function(df) {
      assert_data_frame(df, "data_long")
      assert_vars_in_data(c(id_variable, tv_long), df, "id_variable/time_variable_long", "data_long")
      clash2 <- intersect(c(episode_variable, gap_variable, time_in_episode_variable,
                             event_any_variable, event_type_variable, episode_id_variable), names(df))
      if (length(clash2) > 0) {
        stop(sprintf(
          "`data_long` already has column(s) %s, which collide with buildRecurrentEpisodes()'s output column names. Rename those columns first, or pass different *_variable argument(s).",
          paste(sprintf("`%s`", clash2), collapse = ", ")
        ), call. = FALSE)
      }

      carry_new <- setdiff(carry_cols, names(df))
      n <- nrow(df)

      # Which episode(s) each row's timestamp falls into, using a CLOSED
      # [t_lo, t_end] interval for every episode. A timestamp strictly
      # inside exactly one episode matches only that episode, as expected.
      # A timestamp exactly AT a shared boundary between episode k (ending
      # there) and episode k + 1 (starting there) matches *both* -- this is
      # a genuine ambiguity, not an error to resolve one way: that instant
      # is simultaneously "the last measurement of the ending episode" and
      # "the first measurement (time_in_episode = 0) of the next one", and
      # both readings are informative for the backward-time trick, so both
      # rows are kept (i.e. the input row is duplicated, once per episode
      # it borders). Silently keeping only one (e.g. always the ending
      # episode) would systematically discard every fresh-start measurement
      # taken exactly when a new episode begins -- plausible in real data
      # whenever a biomarker is recorded at the same visit an event is
      # recorded -- and bias the fitted backward covariate's effect instead
      # of just losing a documented edge case.
      blk <- match(id_label(df[[id_variable]]), uid)               # which subject's block of episodes
      n_blk <- ifelse(is.na(blk), 0L, n_ep_by_id[blk])
      cand_row <- rep(seq_len(n), n_blk)
      cand_pos <- rep(first_pos[blk], n_blk) + sequence(n_blk) - 1L
      tm <- df[[tv_long]][cand_row]
      hit <- which(!is.na(tm) & tm >= t_lo[cand_pos] & tm <= t_end[cand_pos])
      row_idx <- cand_row[hit]; ep_pos <- cand_pos[hit]

      dropped <- n - length(unique(row_idx))
      if (dropped > 0) {
        message(sprintf(
          "buildRecurrentEpisodes(): dropped %d row(s) of `data_long` whose %s fell outside every modeled episode (before time 0 or the `lookback` window, inside a hospital stay not covered by `lookback`, or after the subject's recurrent-event process had already ended).",
          dropped, tv_long
        ))
      }
      doubled <- length(row_idx) - length(unique(row_idx))
      if (doubled > 0) {
        message(sprintf(
          "buildRecurrentEpisodes(): %d row(s) of `data_long` fall exactly on the boundary between two episodes and are kept in both (once as the end of one, once as the start of the next).",
          doubled
        ))
      }

      out <- df[row_idx, , drop = FALSE]
      rownames(out) <- NULL
      gp <- gap_all[ep_pos]
      gp[episode_data[[event_any_variable]][ep_pos] != 1] <- NA_real_     # an unfinished episode has no length yet
      out[[episode_variable]] <- ep_no[ep_pos]
      out[[episode_id_variable]] <- factor(paste(id_label(out[[id_variable]]), ep_no[ep_pos], sep = ":"))
      out[[gap_variable]] <- gp
      out[[time_in_episode_variable]] <- out[[tv_long]] - t_start[ep_pos]
      out[[event_any_variable]] <- episode_data[[event_any_variable]][ep_pos]
      out[[event_type_variable]] <- episode_data[[event_type_variable]][ep_pos]
      if (length(carry_new) > 0) {
        carried <- episode_data[ep_pos, carry_new, drop = FALSE]
        rownames(carried) <- NULL
        out[carry_new] <- carried
      }
      out
    }

    out_list <- lapply(data_long_list, augment_one)
    long_data <- if (long_is_list) out_list else out_list[[1]]
  }

  list(episode_data = episode_data, long_data = long_data)
}

## A subject's id as the string used to match rows across tables and to build
## episode ids: a number is written out in full (as.character(1e5) is "1e+05"),
## and a factor, character or number that print the same match each other.
id_label <- function(x) {
  if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE, digits = 15) else as.character(x)
}
