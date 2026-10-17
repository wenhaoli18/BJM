#' Reshape recurrent-event data into gap-time episodes for the backward joint model
#'
#' @description
#' \strong{Experimental.} \code{survivalSub()}/\code{longitudinalSub()}/
#' \code{predictRisk()} were developed and validated (Shen and Li, 2021;
#' Li et al., 2026 -- see \code{citation("BJM")}) for a single terminal event,
#' optionally with competing risks: every subject contributes exactly one
#' row of survival data and one backward time-to-event covariate.
#' \code{buildRecurrentEpisodes()} extends the same backward-time trick to
#' \emph{recurrent} events by a gap-time (renewal) reformulation: it splits
#' each subject's event history into episodes (the intervals between
#' consecutive events, gap-time-clock reset to \code{0} at the start of each
#' episode), and re-clocks any accompanying longitudinal measurements onto
#' the same per-episode timeline. The output plugs directly into the
#' existing \code{survivalSub()}/\code{longitudinalSub()}/\code{predictRisk()}
#' pipeline -- no changes to those functions are needed -- to predict, for a
#' subject currently at risk within an ongoing episode, the probability that
#' their \emph{next} event (within a horizon) is another recurrence versus a
#' terminal event. This is unpublished, unvalidated-by-simulation
#' methodology (unlike the single-event and competing-risks cases cited
#' above); treat its output accordingly, and note that it only predicts the
#' next event, not the total number of future recurrences -- see the package
#' README/vignette discussion of the gap-time extension for why the latter
#' needs materially more machinery.
#'
#' @details
#' \strong{What counts as an episode.} Episode 1 runs from time \code{0} to
#' a subject's first row in \code{data_events}; episode \code{k > 1} runs
#' from the end of episode \code{k - 1} to that subject's \code{k}-th row.
#' A subject contributes no further episodes after a row with
#' \code{status \%in\% c(0, 2)} (administrative censoring or a terminal
#' event); any later rows for that \code{id} are ignored. \code{status == 1}
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
#' \code{\link{pbc3}}.
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
#' @param data_events One row per subject per event-or-censoring occurrence,
#'   in calendar time (i.e. not yet re-clocked to gap time -- this function
#'   does that). Must contain \code{id_variable}, \code{time_variable}
#'   (numeric, strictly increasing within each subject), and
#'   \code{status_variable} (\code{0} = censored/end of follow-up, \code{1} =
#'   recurrent event, \code{2} = terminal event -- the same 0/1/2 convention
#'   \link{pbc3}'s \code{status5} already uses). Any other columns
#'   (baseline or event-time-varying covariates) are carried through to
#'   \code{episode_data}, taken from the same row as the event they end.
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
#'   episode's event time, and bias the fitted backward-covariate effect).
#' @param id_variable,time_variable,status_variable Column names in
#'   \code{data_events} (defaults \code{"id"}, \code{"time"}, \code{"status"}).
#' @param time_variable_long Column name of the measurement time in
#'   \code{data_long}; defaults to \code{time_variable} (i.e. the same name
#'   as in \code{data_events}) when \code{NULL}.
#' @param episode_variable,gap_variable,time_in_episode_variable,event_any_variable,event_type_variable
#'   Column names to create in the output (defaults \code{"episode"},
#'   \code{"gap"}, \code{"time_in_episode"}, \code{"event_any"},
#'   \code{"event_type"}). Change these if they collide with an existing
#'   column of \code{data_events}/\code{data_long} -- collisions are an
#'   error, not a silent overwrite.
#' @param restart_variable Optional name of a numeric column of
#'   \code{data_events}: for each recurrence (\code{status == 1}) row, the time
#'   the process restarts after it (for hospitalizations, the discharge ending
#'   that stay; must not be earlier than the row's \code{time_variable}). The
#'   next episode starts there instead of at the event time. \code{NULL}
#'   (default) restarts at the event time itself. Ignored on rows that end the
#'   process. Not carried into \code{episode_data}.
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
#' data_events <- data.frame(
#'   id     = c(1, 1, 1, 2, 2),
#'   time   = c(1.2, 2.5, 4.0, 3.1, 3.1 + 1e-9), # subject 2's 2nd row just illustrates ordering
#'   status = c(1,   1,   2,   0,   0),
#'   age    = c(58,  58,  58,  61,  61)
#' )
#' data_events <- data_events[-5, ] # keep one row per subject/episode boundary above
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
#' ## build episodes -> survivalSub() -> longitudinalSub() -> predictRisk().
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
#' # longitudinalSub()'s id-grouping needs a factor id, exactly as pbc3$id is
#' # already a factor -- this is required regardless of how the data was
#' # reshaped, not something specific to buildRecurrentEpisodes().
#' built$long_data$id <- factor(built$long_data$id)
#'
#' ## gap-time recurrence hazard + recurrence-vs-terminal-event type model
#' survival_fit_recur <- survivalSub(
#'   built$episode_data,
#'   Surv(gap, event_any) ~ age + sex,
#'   event_type ~ gap + age + sex
#' )
#'
#' ## longitudinal sub-model, trained on completed episodes only
#' long_fit_recur <- longitudinalSub(
#'   built$long_data[built$long_data$event_any == 1, ],
#'   biomarker ~ time_in_episode + age + sex + gap + gap * time_in_episode,
#'   ~ time_in_episode | id
#' )
#'
#' ## predict, for subject 9999, the risk that their *next* event (within 1
#' ## time unit) is a recurrence vs. a terminal event, using their history in
#' ## the episode they are currently at risk in
#' data_predict <- built$long_data[built$long_data$id == 9999, ]
#' cur_episode <- max(data_predict$episode)
#' landmark <- max(data_predict$time_in_episode[data_predict$episode == cur_episode])
#' data_predict <- data_predict[data_predict$episode == cur_episode, ]
#' ## (the open episode's gap is NA, which predictRisk() reads as "still at
#' ## risk" -- nothing to override)
#'
#' risk <- predictRisk(
#'   data_predict,
#'   long_fit_recur, survival_fit_recur,
#'   prediction_time = landmark, horizon = 1, time_variable = "time_in_episode",
#'   survival_variable_all = list(), survival_trans_function = list(),
#'   bandcount1 = 10, bandcount2 = 10
#' )
#' risk$risk_prob_1 # P(next event within horizon is a recurrence)
#' risk$risk_prob_2 # P(next event within horizon is the terminal event)
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
                                    lookback = 0) {

  assert_data_frame(data_events, "data_events")
  for (nm in c("id_variable", "time_variable", "status_variable", "episode_variable",
               "gap_variable", "time_in_episode_variable", "event_any_variable",
               "event_type_variable")) {
    assert_string(get(nm), nm)
  }
  assert_vars_in_data(c(id_variable, time_variable, status_variable), data_events,
                       "id_variable/time_variable/status_variable", "data_events")

  new_cols <- c(episode_variable, gap_variable, event_any_variable, event_type_variable)
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

  baseline_cols <- setdiff(names(data_events),
                           c(id_variable, time_variable, status_variable, restart_variable))
  ids <- unique(data_events[[id_variable]])

  episode_rows <- vector("list", length(ids))
  names(episode_rows) <- as.character(ids)
  bounds_rows <- episode_rows

  for (pid in ids) {
    di <- data_events[data_events[[id_variable]] == pid, , drop = FALSE]
    # an episode runs from t_start (time 0 for the first; afterwards the
    # previous event's `restart_variable` if given -- e.g. discharge from the
    # stay a readmission started -- else the previous event time itself) to
    # t_end (this row's time)
    t_start <- 0
    t_prev_event <- -Inf
    rows <- list(); brows <- list()
    k <- 0L
    for (r in seq_len(nrow(di))) {
      st <- di[[status_variable]][r]
      tt <- di[[time_variable]][r]
      if (!is.finite(tt) || tt <= t_start) {
        stop(sprintf(
          "`%s` must be strictly increasing within each `%s` -- and later than 0 for the first row, and later than the previous row's `restart_variable` if one is given; id = %s violates this at row %d.",
          time_variable, id_variable, pid, r
        ), call. = FALSE)
      }
      k <- k + 1L
      row_out <- di[r, baseline_cols, drop = FALSE]
      row_out[[id_variable]] <- pid
      row_out[[episode_variable]] <- k
      row_out[[gap_variable]] <- tt - t_start
      row_out[[event_any_variable]] <- as.integer(st != 0)
      row_out[[event_type_variable]] <- if (st == 1) 0L else if (st == 2) 1L else NA_integer_
      rows[[k]] <- row_out
      brows[[k]] <- data.frame(episode = k, t_start = t_start, t_end = tt, t_prev_event = t_prev_event)
      if (st != 1) break
      t_prev_event <- tt
      t_start <- tt
      if (!is.null(restart_variable) && !is.na(di[[restart_variable]][r])) {
        rs <- di[[restart_variable]][r]
        if (rs < tt) {
          stop(sprintf("`%s` (%s) is earlier than the event time it belongs to (%s) for id = %s, row %d.",
                       restart_variable, format(rs), format(tt), pid, r), call. = FALSE)
        }
        t_start <- rs
      }
    }
    episode_rows[[as.character(pid)]] <- do.call(rbind, rows)
    bounds_rows[[as.character(pid)]] <- do.call(rbind, brows)
  }
  episode_data <- do.call(rbind, episode_rows)
  rownames(episode_data) <- NULL
  episode_data <- episode_data[, c(id_variable, episode_variable, gap_variable,
                                    event_any_variable, event_type_variable, baseline_cols),
                                drop = FALSE]

  long_data <- NULL
  if (!is.null(data_long)) {
    long_is_list <- is.list(data_long) && !is.data.frame(data_long)
    data_long_list <- if (long_is_list) data_long else list(data_long)

    tv_long <- if (is.null(time_variable_long)) time_variable else time_variable_long
    assert_string(tv_long, "time_variable_long")

    # per subject: each episode's start, end, and the earliest time whose
    # measurements count toward it (t_lo). With `lookback = 0` (the default)
    # t_lo is the episode's own start, i.e. only measurements from t_start on;
    # a larger `lookback` also takes in measurements up to that long before
    # t_start -- but never from before the previous event, which belong to the
    # previous episode's own at-risk period -- giving them negative
    # time_in_episode. Between a readmission and the discharge that restarts
    # the process (the stay) that window is exactly the stay.
    boundaries <- lapply(bounds_rows, function(b) {
      data.frame(episode = b$episode, t_start = b$t_start, t_end = b$t_end,
                 t_lo = pmax(b$t_start - lookback, b$t_prev_event),
                 gap = b$t_end - b$t_start)
    })
    ep_lookup <- split(episode_data, episode_data[[id_variable]])
    # Baseline/event-time covariates already carried onto episode_data (e.g.
    # age, sex) get carried onto long_data too, exactly as every row of
    # pbc3's long-format table already repeats age/sex/years/status3 --
    # but only into columns data_long does not already have, so a genuinely
    # different variable that happens to share a name is left alone.
    carry_cols <- setdiff(names(episode_data),
                           c(id_variable, episode_variable, gap_variable,
                             event_any_variable, event_type_variable))

    augment_one <- function(df) {
      assert_data_frame(df, "data_long")
      assert_vars_in_data(c(id_variable, tv_long), df, "id_variable/time_variable_long", "data_long")
      clash2 <- intersect(c(episode_variable, gap_variable, time_in_episode_variable,
                             event_any_variable, event_type_variable), names(df))
      if (length(clash2) > 0) {
        stop(sprintf(
          "`data_long` already has column(s) %s, which collide with buildRecurrentEpisodes()'s output column names. Rename those columns first, or pass different *_variable argument(s).",
          paste(sprintf("`%s`", clash2), collapse = ", ")
        ), call. = FALSE)
      }

      carry_new <- setdiff(carry_cols, names(df))
      n <- nrow(df)

      # Which episode(s) each row's timestamp falls into, using a CLOSED
      # [t_start, t_end] interval for every episode. A timestamp strictly
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
      row_idx <- integer(0); hit_idx <- integer(0)
      for (r in seq_len(n)) {
        pid <- as.character(df[[id_variable]][r])
        tm <- df[[tv_long]][r]
        b <- boundaries[[pid]]
        if (is.null(b) || is.na(tm)) next
        hits <- which(tm >= b$t_lo & tm <= b$t_end)
        if (length(hits) == 0) next
        row_idx <- c(row_idx, rep(r, length(hits)))
        hit_idx <- c(hit_idx, hits)
      }

      dropped <- n - length(unique(row_idx))
      if (dropped > 0) {
        message(sprintf(
          "buildRecurrentEpisodes(): dropped %d row(s) of `data_long` whose %s fell outside every modeled episode (before time 0 or the `lookback` window, inside a hospital stay not covered by `lookback`, or after the subject's recurrent-event process had already ended).",
          dropped, tv_long
        ))
      }

      out <- df[row_idx, , drop = FALSE]
      rownames(out) <- NULL
      m <- nrow(out)
      ep <- rep(NA_integer_, m); gp <- rep(NA_real_, m); tie <- rep(NA_real_, m)
      ev_any <- rep(NA_integer_, m); ev_type <- rep(NA_integer_, m)
      carried <- if (length(carry_new) > 0) {
        tmpl <- episode_data[rep(1, max(m, 1)), carry_new, drop = FALSE][seq_len(m), , drop = FALSE]
        tmpl[] <- NA
        rownames(tmpl) <- NULL
        tmpl
      } else {
        NULL
      }

      for (j in seq_len(m)) {
        pid <- as.character(out[[id_variable]][j])
        b <- boundaries[[pid]]
        hit <- hit_idx[j]
        ep[j] <- b$episode[hit]
        tie[j] <- out[[tv_long]][j] - b$t_start[hit]
        erow <- ep_lookup[[pid]]
        erow <- erow[erow[[episode_variable]] == b$episode[hit], , drop = FALSE]
        ev_any[j] <- erow[[event_any_variable]][1]
        ev_type[j] <- erow[[event_type_variable]][1]
        # An episode that has not ended (event_any == 0) has no observed
        # length yet, so its gap is NA here -- exactly what predictRisk()
        # reads as "event time not yet known, still at risk". (episode_data
        # keeps the real elapsed time for the Cox fit.)
        gp[j] <- if (ev_any[j] == 1) b$gap[hit] else NA_real_
        if (!is.null(carried)) carried[j, ] <- erow[1, carry_new]
      }

      out[[episode_variable]] <- ep
      out[[gap_variable]] <- gp
      out[[time_in_episode_variable]] <- tie
      out[[event_any_variable]] <- ev_any
      out[[event_type_variable]] <- ev_type
      if (!is.null(carried)) out[carry_new] <- carried
      out
    }

    out_list <- lapply(data_long_list, augment_one)
    long_data <- if (long_is_list) out_list else out_list[[1]]
  }

  list(episode_data = episode_data, long_data = long_data)
}
