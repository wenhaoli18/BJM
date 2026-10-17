# buildRecurrentEpisodes() reshapes a calendar-time recurrent-event history
# (data_events) -- plus, optionally, repeated-measurement data (data_long) --
# into gap-time episodes for the gap-time/renewal extension of the backward
# joint model (see ?buildRecurrentEpisodes).

test_that("buildRecurrentEpisodes computes gap times and stops after censoring/terminal events", {
  # id 1: two recurrences then a terminal event (3 episodes).
  # id 2: administratively censored at their very first row (1 episode).
  # id 3: one recurrence then censored (2 episodes); a spurious 3rd row
  # after censoring must be ignored.
  data_events <- data.frame(
    id     = c(1, 1, 1, 2, 3, 3, 3),
    time   = c(1.0, 2.5, 4.0, 3.0, 0.5, 1.5, 9.0),
    status = c(1,   1,   2,   0,   1,   0,   2),
    age    = c(58,  58,  58,  61,  70,  70,  70)
  )

  expect_message(out <- buildRecurrentEpisodes(data_events), "ignored 1 row")
  ed <- out$episode_data

  expect_null(out$long_data)
  expect_equal(nrow(ed), 6) # 3 + 1 + 2, id 3's 3rd row dropped
  expect_setequal(names(ed), c("id", "episode", "episode_id", "gap", "event_any", "event_type", "age"))

  id1 <- ed[ed$id == 1, ]
  id1 <- id1[order(id1$episode), ]
  expect_equal(id1$episode, 1:3)
  expect_equal(id1$gap, c(1.0, 1.5, 1.5))
  expect_equal(id1$event_any, c(1L, 1L, 1L))
  expect_equal(id1$event_type, c(0L, 0L, 1L)) # recurrence, recurrence, terminal

  id2 <- ed[ed$id == 2, ]
  expect_equal(id2$episode, 1)
  expect_equal(id2$gap, 3.0)
  expect_equal(id2$event_any, 0L)
  expect_true(is.na(id2$event_type))

  id3 <- ed[ed$id == 3, ]
  id3 <- id3[order(id3$episode), ]
  expect_equal(nrow(id3), 2) # the row after censoring is dropped
  expect_equal(id3$gap, c(0.5, 1.0))
  expect_equal(id3$event_any, c(1L, 0L))
  expect_equal(id3$event_type, c(0L, NA_integer_))

  # baseline covariate carried straight through from data_events
  expect_equal(unique(id1$age), 58)
})

test_that("buildRecurrentEpisodes sorts unsorted data_events by id then time first", {
  # id 1's rows are given out of chronological order: (time = 2, censor) is
  # listed before (time = 1, recurrence). Sorted correctly, the recurrence
  # at time 1 comes first and the subject is still at risk for a 2nd
  # episode, ending in censoring at time 2; processed in the given
  # (unsorted) order instead, the censoring row would be seen first and
  # stop the process after a single episode.
  data_events <- data.frame(
    id     = c(2, 1, 1, 2),
    time   = c(3.0, 2.0, 1.0, 1.0),
    status = c(0,   0,   1,   0)
  )
  expect_message(out <- buildRecurrentEpisodes(data_events), "ignored 1 row")
  id1 <- out$episode_data[out$episode_data$id == 1, ]
  id1 <- id1[order(id1$episode), ]
  expect_equal(nrow(id1), 2)
  expect_equal(id1$gap, c(1.0, 1.0))
  expect_equal(id1$event_any, c(1L, 0L))
})

test_that("buildRecurrentEpisodes rejects invalid status/time/id columns", {
  expect_error(
    buildRecurrentEpisodes(data.frame(id = 1, time = 1, status = 3)),
    "must only contain 0"
  )
  expect_error(
    buildRecurrentEpisodes(data.frame(id = 1, time = 1, status = NA)),
    "must only contain 0"
  )
  expect_error(
    buildRecurrentEpisodes(data.frame(id = 1, time = "a", status = 0)),
    "must be a numeric column"
  )
  expect_error(
    buildRecurrentEpisodes(data.frame(id = NA, time = 1, status = 0)),
    "must not contain missing values"
  )
})

test_that("buildRecurrentEpisodes rejects duplicate event times within a subject", {
  # data_events is sorted by (id, time) before processing, so two distinct
  # times for the same id always come out in increasing order; the only way
  # to actually violate strict monotonicity afterwards is a tie.
  expect_error(
    buildRecurrentEpisodes(data.frame(id = c(1, 1), time = c(2, 2), status = c(1, 0))),
    "must be strictly increasing"
  )
})

test_that("buildRecurrentEpisodes errors on output-column collisions instead of silently overwriting", {
  expect_error(
    buildRecurrentEpisodes(data.frame(id = 1, time = 1, status = 1, gap = 99)),
    "collide with buildRecurrentEpisodes"
  )
})

test_that("buildRecurrentEpisodes re-clocks data_long onto per-episode gap time", {
  data_events <- data.frame(
    id = c(1, 1), time = c(1.0, 3.0), status = c(1, 2), age = c(58, 58)
  )
  data_long <- data.frame(
    id = 1, year = c(0, 0.5, 1.0, 1.5, 2.9, 3.0, 3.5)
  )

  msgs <- capture_messages(out <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year"))
  expect_match(msgs, "dropped 1 row", all = FALSE)
  expect_match(msgs, "boundary between two episodes", all = FALSE)
  ld <- out$long_data

  # the final row (year = 3.5) is after the process ended (2 episodes only).
  # year = 1.0 sits exactly on the shared boundary between episode 1 (ending
  # there) and episode 2 (starting there) and is genuinely ambiguous, so it
  # is kept *twice* -- once per episode -- rather than arbitrarily assigned
  # to only one of them (see ?buildRecurrentEpisodes).
  expect_equal(nrow(ld), 7)
  expect_equal(sort(ld$year), sort(c(0, 0.5, 1.0, 1.0, 1.5, 2.9, 3.0)))

  ep1 <- ld[ld$episode == 1, ]
  expect_equal(sort(ep1$year), c(0, 0.5, 1.0))
  expect_true(all(ep1$gap == 1.0))
  expect_equal(ep1$time_in_episode, ep1$year) # episode 1 starts at t = 0
  expect_true(all(ep1$event_any == 1))
  expect_true(all(ep1$event_type == 0)) # episode 1 ended in a recurrence

  ep2 <- ld[ld$episode == 2, ]
  expect_equal(sort(ep2$year), c(1.0, 1.5, 2.9, 3.0))
  expect_true(all(ep2$gap == 2.0))
  expect_equal(ep2$time_in_episode, ep2$year - 1.0) # episode 2 starts at t = 1
  expect_true(all(ep2$event_type == 1)) # episode 2 ended in the terminal event

  # the duplicated boundary row itself: time_in_episode = 1.0 (the full gap)
  # as episode 1's last point, and 0 (fresh start) as episode 2's first
  boundary_rows <- ld[ld$year == 1.0, ]
  expect_equal(sort(boundary_rows$time_in_episode), c(0, 1.0))

  # baseline covariate (age) carried from episode_data onto every long_data row
  expect_true(all(ld$age == 58))
})

test_that("buildRecurrentEpisodes drops out-of-range data_long rows with a message", {
  data_events <- data.frame(id = 1, time = 2.0, status = 2)
  data_long <- data.frame(id = 1, year = c(-1, 0, 1, 2, 3))

  expect_message(
    out <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year"),
    "dropped 2 row"
  )
  expect_equal(out$long_data$year, c(0, 1, 2))
})

test_that("buildRecurrentEpisodes leaves an existing same-named data_long column alone", {
  data_events <- data.frame(id = 1, time = 2.0, status = 2, age = 70)
  data_long <- data.frame(id = 1, year = c(0, 1), age = c(999, 999))

  out <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year")
  # data_long already had its own `age` column: not overwritten with
  # episode_data's age (70)
  expect_equal(out$long_data$age, c(999, 999))
})

test_that("buildRecurrentEpisodes accepts a list of data_long data.frames (one per longitudinal outcome)", {
  data_events <- data.frame(id = 1, time = 2.0, status = 1)
  data_events <- rbind(data_events, data.frame(id = 1, time = 4.0, status = 0))
  data_long_list <- list(
    data.frame(id = 1, year = c(0, 1)),
    data.frame(id = 1, year = c(0, 1, 3))
  )

  out <- buildRecurrentEpisodes(data_events, data_long_list, time_variable_long = "year")
  expect_true(is.list(out$long_data) && !is.data.frame(out$long_data))
  expect_length(out$long_data, 2)
  expect_equal(nrow(out$long_data[[1]]), 2)
  expect_equal(nrow(out$long_data[[2]]), 3)
})

test_that("buildRecurrentEpisodes errors on missing required columns", {
  expect_error(
    buildRecurrentEpisodes(data.frame(time = 1, status = 0)),
    "not found in `data_events`"
  )
  data_events <- data.frame(id = 1, time = 1, status = 1)
  data_events <- rbind(data_events, data.frame(id = 1, time = 2, status = 0))
  expect_error(
    buildRecurrentEpisodes(data_events, data.frame(id = 1, day = 0)),
    "not found in `data_long`"
  )
})

test_that("the gap-time episode_data/long_data plug directly into survivalSub()/longitudinalSub()/predictRisk()", {
  set.seed(42)
  n <- 40
  data_events <- do.call(rbind, lapply(seq_len(n), function(id) {
    k <- sample(1:3, 1)
    gaps <- rexp(k, rate = 1.2)
    data.frame(
      id = id, time = cumsum(gaps),
      status = c(rep(1, k - 1), sample(c(0, 2), 1, prob = c(0.6, 0.4))),
      age = round(50 + stats::rnorm(1, sd = 8)), sex = id %% 2
    )
  }))
  extra_id <- 9999
  data_events <- rbind(data_events, data.frame(
    id = extra_id, time = c(0.8, 1.6), status = c(1, 0), age = 55, sex = 1
  ))

  data_long <- do.call(rbind, lapply(c(seq_len(n), extra_id), function(id) {
    maxt <- max(data_events$time[data_events$id == id]) + 0.2
    tt <- seq(0, maxt, by = 0.25)
    data.frame(id = id, year = tt, biomarker = stats::rnorm(length(tt)))
  }))

  expect_message(built <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year"),
                 "dropped 38 row")
  built$long_data$id <- factor(built$long_data$id)

  survival_fit_recur <- survivalSub(
    built$episode_data,
    Surv(gap, event_any) ~ age + sex,
    event_type ~ gap + age + sex
  )
  expect_s3_class(survival_fit_recur, "survivalSub.BJM")

  long_fit_recur <- longitudinalSub(
    built$long_data[built$long_data$event_any == 1, ],
    biomarker ~ time_in_episode + age + sex + gap + gap * time_in_episode,
    ~ time_in_episode | id
  )
  expect_s3_class(long_fit_recur, "longitudinalSub.BJM")

  data_predict <- built$long_data[built$long_data$id == extra_id, ]
  cur_episode <- max(data_predict$episode)
  landmark <- max(data_predict$time_in_episode[data_predict$episode == cur_episode])

  risk <- predictRisk(
    data_predict[data_predict$episode == cur_episode, ],
    long_fit_recur, survival_fit_recur,
    prediction_time = landmark, horizon = 1, time_variable = "time_in_episode",
    survival_variable_all = list(), survival_trans_function = list(),
    bandcount1 = 10, bandcount2 = 10
  )
  expect_s3_class(risk, "predictRisk.BJM")
  expect_true(is.numeric(risk$risk_prob_1) && risk$risk_prob_1 >= 0 && risk$risk_prob_1 <= 1)
  expect_true(is.numeric(risk$risk_prob_2) && risk$risk_prob_2 >= 0 && risk$risk_prob_2 <= 1)
})

# ---- gap-time episodes that start at discharge (restart_variable/lookback) ----

test_that("the open episode's gap is NA in long_data, but episode_data keeps the elapsed time", {
  data_events <- data.frame(id = c(1, 1), time = c(1.0, 2.5), status = c(1, 0))
  data_long <- data.frame(id = 1, year = c(0.5, 1.5, 2.0))

  out <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year")

  expect_equal(out$episode_data$gap, c(1.0, 1.5))     # Cox fit needs the censoring time
  ld <- out$long_data
  expect_true(all(ld$gap[ld$episode == 1] == 1.0))    # resolved episode: real gap
  expect_true(all(is.na(ld$gap[ld$episode == 2])))    # open episode: unknown
  expect_true(all(ld$event_any[ld$episode == 2] == 0))
})

test_that("restart_variable starts the next episode at discharge and keeps stay labs as negative time", {
  # readmitted day 20, discharged day 26; readmitted day 60, discharged day 66;
  # censored day 90
  events <- data.frame(id = 1, time = c(20, 60, 90), status = c(1, 1, 0),
                       discharge = c(26, 66, NA))
  labs <- data.frame(id = 1, day = c(-3, 0, 5, 21, 26, 40, 61, 66, 80))

  out <- buildRecurrentEpisodes(events, labs, time_variable_long = "day",
                                restart_variable = "discharge", lookback = Inf)
  ed <- out$episode_data
  expect_equal(ed$gap, c(20, 34, 24))                 # 20-0, 60-26, 90-66
  expect_false("discharge" %in% names(ed))            # restart column not carried

  ld <- out$long_data
  ep1 <- ld[ld$episode == 1, ]
  expect_equal(sort(ep1$time_in_episode), c(-3, 0, 5))   # index-stay lab at -3, discharge 0, day 5
  # day 21 is inside the stay after the day-20 readmission: belongs to episode 2 as -5
  ep2 <- ld[ld$episode == 2, ]
  expect_equal(sort(ep2$time_in_episode), sort(c(21, 26, 40) - 26))
  ep3 <- ld[ld$episode == 3, ]
  expect_equal(sort(ep3$time_in_episode), sort(c(61, 66, 80) - 66))
})

test_that("lookback = 0 (default) drops stay labs and restart_variable alone changes only the clock", {
  events <- data.frame(id = 1, time = c(20, 60), status = c(1, 0), discharge = c(26, NA))
  labs <- data.frame(id = 1, day = c(5, 21, 26, 40))

  expect_message(out <- buildRecurrentEpisodes(events, labs, time_variable_long = "day",
                                               restart_variable = "discharge"),   # lookback = 0
                 "dropped 1 row")
  expect_equal(out$episode_data$gap, c(20, 34))
  # day 21 is in the stay (20, 26): not kept without lookback; day 26 = start of episode 2
  expect_equal(sort(out$long_data$day), c(5, 26, 40))
  expect_equal(out$long_data$time_in_episode[out$long_data$day == 26], 0)
})

test_that("a missing restart time on a recurrence is reported and taken as no stay", {
  events <- data.frame(id = 1, time = c(20, 60, 90), status = c(1, 1, 0), discharge = c(26, NA, NA))
  expect_warning(out <- buildRecurrentEpisodes(events, restart_variable = "discharge"),
                 "missing on 1 recurrence row")
  expect_equal(out$episode_data$gap, c(20, 34, 30))   # episode 3 starts at the day-60 readmission
  # the last row ends the process, so its missing restart is not reported
  expect_no_warning(buildRecurrentEpisodes(events[-3, ], restart_variable = "discharge"))
})

test_that("lookback never reaches back past the previous event", {
  events <- data.frame(id = 1, time = c(20, 60), status = c(1, 0), discharge = c(26, NA))
  labs <- data.frame(id = 1, day = c(10, 21))   # day 10 is episode 1's own at-risk time
  out <- buildRecurrentEpisodes(events, labs, time_variable_long = "day",
                                restart_variable = "discharge", lookback = 100)
  ep2 <- out$long_data[out$long_data$episode == 2, ]
  expect_equal(ep2$day, 21)       # day 10 (before the day-20 readmission) is NOT episode 2's
  expect_true(10 %in% out$long_data$day[out$long_data$episode == 1])
})

test_that("restart earlier than its own event, bad lookback, or a non-numeric restart column are rejected", {
  expect_error(
    buildRecurrentEpisodes(data.frame(id = c(1, 1), time = c(20, 60), status = c(1, 0), r = c(10, NA)),
                           restart_variable = "r"),
    "earlier than the event time")
  expect_error(
    buildRecurrentEpisodes(data.frame(id = 1, time = 5, status = 0), lookback = -1),
    "`lookback` must be")
  expect_error(
    buildRecurrentEpisodes(data.frame(id = 1, time = 5, status = 0, r = "a"), restart_variable = "r"),
    "must be a numeric column")
  # next event before the restart (readmitted before being discharged)
  expect_error(
    buildRecurrentEpisodes(data.frame(id = c(1, 1), time = c(20, 24), status = c(1, 0), r = c(26, NA)),
                           restart_variable = "r"),
    "must be strictly increasing")
})

test_that("episode_id identifies each episode, in episode_data and long_data", {
  events <- data.frame(id = c(7, 7, 9), time = c(1.0, 2.5, 3.0), status = c(1, 0, 2))
  labs <- data.frame(id = c(7, 7, 7, 9), day = c(0.2, 1.4, 2.0, 1.0))
  out <- buildRecurrentEpisodes(events, labs, time_variable_long = "day")
  expect_equal(as.character(out$episode_data$episode_id), c("7:1", "7:2", "9:1"))
  expect_s3_class(out$episode_data$episode_id, "factor")
  expect_s3_class(out$long_data$episode_id, "factor")
  expect_equal(as.character(out$long_data$episode_id),
               paste(out$long_data$id, out$long_data$episode, sep = ":"))
  expect_equal(length(unique(out$long_data$episode_id)), 3)
  # a column of that name already there is an error, like the other output columns
  expect_error(buildRecurrentEpisodes(cbind(events, episode_id = 1)), "collide")
  renamed <- buildRecurrentEpisodes(events, labs, time_variable_long = "day", episode_id_variable = "unit")
  expect_true("unit" %in% names(renamed$long_data) && !"episode_id" %in% names(renamed$long_data))
})
