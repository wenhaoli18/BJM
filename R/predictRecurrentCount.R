#' Predict the number of recurrent events in a window, by rolling episodes forward
#'
#' @description
#' \strong{Experimental.} For a patient at the start of (or partway through)
#' an episode of a gap-time recurrent-event model fit with
#' \code{\link{buildRecurrentEpisodes}}, predicts the distribution of the
#' \emph{number} of recurrences in \code{(prediction_time, prediction_time +
#' horizon]}, together with the probability of dying in it, by Monte Carlo:
#' every patient is simulated forward \code{n_sim} times, episode after
#' episode, until the window closes or the patient dies. \code{predictRisk()}
#' answers "is the \emph{next} event a recurrence or a terminal event, and
#' when"; this answers "how many recurrences".
#'
#' @details
#' Each simulated path goes as follows.
#' \enumerate{
#'   \item The \emph{current} (open) episode's event time and type are drawn
#'   with \code{\link{simulateTrajectory}}, from their posterior given the
#'   patient's biomarker history in that episode -- the same distribution
#'   \code{\link{predictRisk}} integrates, so the first step uses the
#'   biomarkers exactly as \code{predictRisk()} does.
#'   \item If the event is a terminal event inside the window, the path
#'   ends (the patient dies). If it is a recurrence inside the window, the
#'   count goes up by one; a stay length is then drawn (from
#'   \code{stay_lengths}) and the next episode starts at discharge. If the
#'   event (or the discharge) falls after the window closes, the path ends.
#'   \item Each \emph{later} episode is drawn the same way, but with no
#'   biomarker history -- the patient has not been observed in it yet -- so
#'   from the model's prior given the patient's covariates, with the episode
#'   counter advanced by \code{next_covariates}. Repeat until no path is
#'   still running, or \code{max_stages} episodes have been drawn.
#' }
#' \strong{Limitation.} Only the first step uses the patient's biomarkers.
#' Later episodes depend on the patient only through their covariates (and
#' the episode counter), so any patient-level tendency to recur that is not
#' in those covariates -- a frailty -- is not carried forward, and the
#' predicted spread of the counts (e.g. the chance of two or more
#' recurrences) will be too small when such a tendency is present. Putting
#' the number of previous episodes, or the length of the previous gap, among
#' the survival sub-model's covariates (and updating them in
#' \code{next_covariates}) recovers some of it.
#'
#' The stay length is not part of the backward joint model, so it has to be
#' supplied. The default, \code{0}, means a recurrence takes no time (a
#' tumour recurrence, say); for hospitalizations pass observed stay lengths
#' (they are resampled with replacement).
#'
#' @param data_predict The rows of the patients' \emph{open} episode, as
#'   \code{predictRisk()} takes them (for example the \code{long_data}
#'   rows of \code{buildRecurrentEpisodes()} for each patient's last
#'   episode, whose \code{gap} is \code{NA}): one \code{data.frame}, used
#'   for every biomarker (a list of data frames, one per biomarker, as
#'   \code{predictRisk()} also accepts, is not supported here), one patient
#'   per \code{id}: more than one episode for the same patient is an error.
#'   A patient whose \code{gap} is not \code{NA} and lies below
#'   \code{prediction_time} counts as no longer at risk, and is an error too
#'   rather than a silent zero. The patient's last row supplies the covariates of later
#'   episodes.
#' @param long_fit_all,survival_fit_all Fits from \code{\link{longitudinalSub}}
#'   and \code{\link{survivalSub}}, on episode data.
#' @param prediction_time The time (on the episode clock, as in
#'   \code{predictRisk()}) the prediction is made at: \code{0} at discharge.
#' @param horizon Length of the window after \code{prediction_time}, in the
#'   same time unit.
#' @param time_variable Name of the time variable of the longitudinal
#'   sub-model (e.g. \code{"time_in_episode"}).
#' @param survival_variable_all,survival_trans_function As in
#'   \code{\link{predictRisk}}.
#' @param episode_variable Name of the episode-counter covariate (default
#'   \code{"episode"}); used only by the default \code{next_covariates}, and
#'   ignored if the data has no such column.
#' @param next_covariates Optional \code{function(row, k)} returning the
#'   one-row covariate data frame of the episode \code{k} episodes after the
#'   open one (\code{k = 1, 2, ...}), given the patient's last row of
#'   \code{data_predict} as \code{row}. The default adds \code{k} to
#'   \code{episode_variable}, but not beyond the largest episode number in the
#'   data the survival sub-model was fit on: past it the fitted effect of the
#'   counter is an extrapolation (a hazard that grows with each earlier
#'   recurrence would be multiplied a hundredfold or more by episode 30), and
#'   the simulated paths would run on through ever shorter episodes.
#' @param stay_lengths A numeric vector of non-negative stay lengths (the
#'   time between a recurrence and the discharge that restarts the process),
#'   resampled with replacement. \code{0} (default): none.
#' @param n_sim Paths simulated per patient.
#' @param max_stages Most episodes drawn per path (including the first);
#'   paths still running at the cap are counted as they stand, with a
#'   warning.
#' @param bandcount2 Event-time grid size passed to
#'   \code{\link{simulateTrajectory}}.
#' @param seed Optional integer seed; the caller's random number stream is
#'   left unchanged.
#'
#' @return An object of class \code{"predictRecurrentCount.BJM"}, a list with
#'   \describe{
#'   \item{summary}{One row per patient: \code{mean_count} (expected number
#'   of recurrences in the window), \code{p_count_0}, \code{p_count_1},
#'   \code{p_count_2}, \code{p_count_3plus}, \code{p_ge1}, \code{p_ge2}, and
#'   \code{p_death} (probability of dying in the window).}
#'   \item{counts}{The raw draws: a patients x \code{n_sim} matrix of
#'   recurrence counts.}
#'   \item{died}{A patients x \code{n_sim} logical matrix: died in the window.}
#'   \item{horizon, prediction_time, n_sim}{As supplied.}
#'   \item{stages}{Most episodes drawn on any path.}
#'   \item{capped}{\code{TRUE} if some path was still running when
#'   \code{max_stages} stopped it (a warning, giving their share, is also
#'   given): its later recurrences are not counted, so the counts are
#'   slightly low.}
#'   \item{capped_fraction}{The share of all simulated paths stopped by
#'   \code{max_stages} (\code{0} if none).}
#'   }
#'
#' @seealso \code{\link{buildRecurrentEpisodes}}, \code{\link{predictRisk}},
#'   \code{\link{simulateTrajectory}}
#' @export
predictRecurrentCount <- function(data_predict, long_fit_all, survival_fit_all,
                                  prediction_time, horizon, time_variable,
                                  survival_variable_all = list(),
                                  survival_trans_function = list(),
                                  episode_variable = "episode",
                                  next_covariates = NULL,
                                  stay_lengths = 0,
                                  n_sim = 200, max_stages = 10,
                                  bandcount2 = 100, seed = NULL) {
  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  assert_data_frame(data_predict, "data_predict")
  assert_scalar_numeric(prediction_time, "prediction_time")
  assert_scalar_numeric(horizon, "horizon", positive = TRUE)
  assert_string(time_variable, "time_variable")
  assert_string(episode_variable, "episode_variable")
  assert_positive_integer(n_sim, "n_sim")
  assert_positive_integer(max_stages, "max_stages")
  if (!is.numeric(stay_lengths) || length(stay_lengths) == 0 || anyNA(stay_lengths) ||
      any(stay_lengths < 0)) {
    stop("`stay_lengths` must be a non-empty numeric vector of non-negative values.", call. = FALSE)
  }
  if (!is.null(next_covariates) && !is.function(next_covariates)) {
    stop("`next_covariates` must be NULL or a function(row, k).", call. = FALSE)
  }
  if (!is.null(seed)) local_r_seed(seed)

  id <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  assert_vars_in_data(c(id, time_variable), data_predict, "the id / time variables", "data_predict")
  id_chr <- as.character(data_predict[[id]])
  if (anyNA(id_chr)) stop(sprintf("`%s` in `data_predict` must not contain missing values.", id), call. = FALSE)
  pids <- unique(id_chr)
  n_p <- length(pids)
  if (episode_variable %in% names(data_predict)) {
    n_episodes <- tapply(data_predict[[episode_variable]], id_chr, function(x) length(unique(x)))
    if (any(n_episodes > 1)) {
      stop(sprintf(paste0(
        "`data_predict` has more than one `%s` for %d patient(s) (e.g. %s). It must hold only each ",
        "patient's open episode, one `%s` per prediction -- predictions for the same patient at ",
        "several discharges need separate ids."),
        episode_variable, sum(n_episodes > 1), paste(first_n(names(n_episodes)[n_episodes > 1], 3), collapse = ", "), id),
        call. = FALSE)
    }
  }
  has_cr <- length(survival_fit_all$form_conditional_cr) != 0
  type_variable <- if (has_cr) all.vars(survival_fit_all$form_conditional_cr[[2]]) else NULL
  response_names <- vapply(long_fit_all$long_sub_fixed, function(f) all.vars(f)[1], character(1))
  survival_variable <- survival_time_variable(survival_fit_all)

  if (is.null(next_covariates)) {
    next_covariates <- default_next_covariates(episode_variable,
                                               training_episode_max(survival_fit_all, episode_variable))
  }

  draw_stage <- function(df, pred_time) {
    sim <- simulateTrajectory(df, long_fit_all, survival_fit_all,
                              prediction_time = pred_time, times = pred_time,
                              time_variable = time_variable,
                              survival_variable_all = survival_variable_all,
                              survival_trans_function = survival_trans_function,
                              n_sim = n_sim, bandcount2 = bandcount2,
                              max_event_time = Inf, truncate = FALSE)
    row_i <- match(as.character(sim[[id]]), pids)
    ev <- ty <- matrix(NA_real_, n_p, n_sim)
    ev[cbind(row_i, sim$sim)] <- sim$event_time
    ty[cbind(row_i, sim$sim)] <- if (has_cr) sim[[type_variable]] else 0
    list(event_time = ev, type = ty)
  }

  ## stage 1: the open episode, given its biomarkers
  s1 <- draw_stage(data_predict, prediction_time)
  no_draw <- pids[rowSums(!is.na(s1$event_time)) == 0]
  if (length(no_draw) > 0) {
    stop(sprintf(paste0(
      "No prediction could be made for %d of %d patient(s) (e.g. %s). The usual cause is that their `%s` ",
      "is not NA and lies below `prediction_time`, which marks them as no longer at risk: for an open ",
      "episode it should be NA (buildRecurrentEpisodes() leaves it so)."),
      length(no_draw), n_p, paste(first_n(no_draw, 3), collapse = ", "), survival_variable), call. = FALSE)
  }
  ev_time <- s1$event_time - prediction_time      # time from the prediction time
  ev_type <- s1$type
  active <- !is.na(ev_time)
  count <- matrix(0L, n_p, n_sim)
  died <- matrix(FALSE, n_p, n_sim)
  last_rows <- data_predict[tapply(seq_len(nrow(data_predict)), id_chr, max)[pids], , drop = FALSE]

  stage <- 1L
  capped <- FALSE
  repeat {
    in_win <- active & ev_time <= horizon
    died[in_win & ev_type == 1] <- TRUE
    recur <- in_win & ev_type == 0
    count[recur] <- count[recur] + 1L
    # stay lengths: index with sample.int (sample(x) on a length-1 x draws from 1:x)
    stay <- matrix(stay_lengths[sample.int(length(stay_lengths), n_p * n_sim, replace = TRUE)],
                   n_p, n_sim)
    start_next <- ev_time + stay
    active <- recur & start_next <= horizon
    if (!any(active, na.rm = TRUE)) break
    if (stage >= max_stages) { capped <- TRUE; break }
    stage <- stage + 1L

    ## a later episode: no biomarker history, covariates advanced by one more episode
    who <- which(rowSums(active) > 0)
    new_rows <- do.call(rbind, lapply(who, function(i) {
      row <- next_covariates(last_rows[i, , drop = FALSE], stage - 1L)
      for (v in c(response_names, survival_variable, type_variable)) row[[v]] <- NA_real_
      row[[time_variable]] <- 0
      row
    }))
    s2 <- draw_stage(new_rows, 0)
    if (any(active & is.na(s2$event_time))) {
      stop("A later episode could not be simulated for some patient(s); check that `data_predict`'s covariates ",
           "(including the episode counter) are complete.", call. = FALSE)
    }
    ev_time <- ifelse(active, start_next + s2$event_time, Inf)
    ev_type <- s2$type
    active[is.na(ev_time)] <- FALSE
  }
  capped_fraction <- if (capped) mean(active, na.rm = TRUE) else 0
  if (capped) {
    warning(sprintf(paste0(
      "%d of %d simulated paths (%s%%) were still running after `max_stages` = %d episodes; their ",
      "remaining recurrences are not counted. Raise `max_stages` if that share matters."),
      sum(active, na.rm = TRUE), length(active), format(signif(100 * capped_fraction, 2)), max_stages),
      call. = FALSE)
  }

  summ <- data.frame(
    id = pids,
    mean_count = rowMeans(count),
    p_count_0 = rowMeans(count == 0L), p_count_1 = rowMeans(count == 1L),
    p_count_2 = rowMeans(count == 2L), p_count_3plus = rowMeans(count >= 3L),
    p_ge1 = rowMeans(count >= 1L), p_ge2 = rowMeans(count >= 2L),
    p_death = rowMeans(died), stringsAsFactors = FALSE)
  names(summ)[1] <- id
  rownames(count) <- rownames(died) <- pids
  out <- list(summary = summ, counts = count, died = died, horizon = horizon,
              prediction_time = prediction_time, n_sim = n_sim, stages = stage, capped = capped,
              capped_fraction = capped_fraction)
  class(out) <- "predictRecurrentCount.BJM"
  out
}

## first n elements (base R only: utils is not imported)
first_n <- function(x, n) x[seq_len(min(n, length(x)))]

#' Default covariate update for the episodes after the open one
#'
#' @description Adds \code{k} to the episode counter, but never past
#' \code{cap}: the largest episode number in the data the survival
#' sub-model was fit on. Beyond it the model's effect of the counter is an
#' extrapolation, and a hazard that rises with each earlier recurrence, run
#' out to episode 30, multiplies the hazard many hundredfold, so paths drag
#' on through ever shorter episodes and the counts are inflated.
#' @keywords internal
default_next_covariates <- function(episode_variable, cap = Inf) {
  function(row, k) {
    if (episode_variable %in% names(row)) row[[episode_variable]] <- min(row[[episode_variable]] + k, cap)
    row
  }
}

#' Largest value of the episode counter the survival sub-model was fit on
#' @keywords internal
training_episode_max <- function(survival_fit_all, episode_variable) {
  x <- survival_fit_all$coxph_fit$x
  if (!is.null(x) && episode_variable %in% colnames(x)) max(x[, episode_variable]) else Inf
}

#' @export
print.predictRecurrentCount.BJM <- function(x, ...) {
  cat(sprintf("Predicted recurrences in (%g, %g], %d patient(s), %d simulated paths each\n",
              x$prediction_time, x$prediction_time + x$horizon, nrow(x$summary), x$n_sim))
  cat(sprintf("mean expected count %.3f; P(>=1) %.3f; P(>=2) %.3f; P(death) %.3f\n\n",
              mean(x$summary$mean_count), mean(x$summary$p_ge1), mean(x$summary$p_ge2),
              mean(x$summary$p_death)))
  print(x$summary[seq_len(min(6, nrow(x$summary))), , drop = FALSE], row.names = FALSE, digits = 3)
  if (nrow(x$summary) > 6) cat(sprintf("... %d more patient(s)\n", nrow(x$summary) - 6))
  invisible(x)
}
