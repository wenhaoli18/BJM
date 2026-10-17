# predictRecurrentCount() rolls gap-time episodes forward by Monte Carlo.

make_count_fixture <- function(n = 250, seed = 11) {
  set.seed(seed)
  ev <- list(); lab <- list()
  for (id in seq_len(n)) {
    age <- rnorm(1, 70, 8); sex <- rbinom(1, 1, 0.5); age_c <- (age - 70) / 10
    u0 <- rnorm(1, 0, 0.5); s <- 0; k <- 0
    repeat {
      k <- k + 1
      G <- (-log(runif(1)) / (0.35 * exp(0.3 * age_c + 0.2 * (k - 1))))^(1 / 1.2)
      type <- rbinom(1, 1, plogis(-2 + 0.5 * G))               # 1 = death
      tv <- c(s - 0.1, s, s + 0.4 * (1:8)); tv <- tv[tv < s + G]
      tie <- tv - s
      u1 <- rnorm(2, 0, 0.3)                                    # episode-level random slopes
      y <- rnorm(length(tv), 5 + u0 + 0.4 * G - 0.5 * G * tie + u1[1] * tie, 0.4)
      lab[[length(lab) + 1]] <- data.frame(id = id, t = tv, bio = y,
        bio2 = 2 + rnorm(1, 0, 0.4) + 0.3 * G + u1[2] * tie + rnorm(length(tv), 0, 0.3), age = age, sex = sex)
      end <- s + G
      if (type == 1) {
        ev[[length(ev) + 1]] <- data.frame(id = id, time = end, status = 2L, disc = NA_real_, age = age, sex = sex)
        break
      }
      ev[[length(ev) + 1]] <- data.frame(id = id, time = end, status = 1L, disc = end + 0.3, age = age, sex = sex)
      if (k >= 4) break
      s <- end + 0.3
    }
  }
  ev <- do.call(rbind, ev); lab <- do.call(rbind, lab)
  b <- buildRecurrentEpisodes(ev, lab, time_variable_long = "t", restart_variable = "disc", lookback = Inf)
  surv_fit <- survivalSub(b$episode_data, Surv(gap, event_any) ~ age + sex + episode,
                          event_type ~ gap + age + sex + episode)
  trn <- b$long_data[b$long_data$event_any == 1, ]
  # random effects by episode (the recommended grouping)
  long_fit <- longitudinalSub(trn, bio ~ time_in_episode + age + sex + gap + gap:time_in_episode,
                              ~ time_in_episode | episode_id)
  long_fit2 <- longitudinalSub(list(trn, trn),
    list(bio ~ time_in_episode + age + sex + gap + gap:time_in_episode,
         bio2 ~ time_in_episode + age + sex + gap + gap:time_in_episode),
    list(~ time_in_episode | episode_id, ~ time_in_episode | episode_id))
  # episodes still going at 0.5, seen only up to 0.5, treated as open (gap and type unknown)
  ld <- b$long_data
  mid <- ld[!is.na(ld$gap) & ld$gap > 0.5 & ld$time_in_episode <= 0.5, ]
  mid <- mid[ave(mid$time_in_episode, as.character(mid$episode_id), FUN = max) > 0.2, ]
  mid <- droplevels(mid[mid$episode_id %in% head(unique(as.character(mid$episode_id)), 30), ])
  mid$gap <- NA_real_; mid$event_type <- NA_integer_; mid$event_any <- 0L
  # new patients at their first discharge (open first episode, labs up to time 0)
  set.seed(seed + 1)
  m <- 25
  new_ev <- data.frame(id = 1000 + seq_len(m), time = 0.001, status = 0L, disc = NA_real_,
                       age = rnorm(m, 70, 8), sex = rbinom(m, 1, 0.5))
  new_lab <- do.call(rbind, lapply(seq_len(m), function(i) {
    u0 <- rnorm(1, 0, 0.5)
    y <- rnorm(2, 5 + u0 + 0.4 * 0.8, 0.4)
    data.frame(id = new_ev$id[i], t = c(-0.1, 0), bio = y, bio2 = 2 + rnorm(1, 0, 0.4) + 0.3 * 0.8 + rnorm(2, 0, 0.3),
               age = new_ev$age[i], sex = new_ev$sex[i])
  }))
  bn <- suppressMessages(buildRecurrentEpisodes(new_ev, new_lab, time_variable_long = "t",
                                                restart_variable = "disc", lookback = Inf))
  list(long_fit = long_fit, long_fit2 = long_fit2, surv_fit = surv_fit, open = bn$long_data,
       mid = mid, stays = rep(0.3, 5))
}

fx <- suppressMessages(make_count_fixture())

predict_counts <- function(..., horizon = 3) {
  predictRecurrentCount(fx$open, fx$long_fit, fx$surv_fit, prediction_time = 0, horizon = horizon,
                        time_variable = "time_in_episode", ...)
}

test_that("predictRecurrentCount returns well-formed probabilities", {
  res <- predict_counts(n_sim = 150, stay_lengths = fx$stays, seed = 1)
  expect_s3_class(res, "predictRecurrentCount.BJM")
  s <- res$summary
  expect_equal(nrow(s), 25)
  expect_true(all(c("mean_count", "p_count_0", "p_count_1", "p_count_2", "p_count_3plus",
                    "p_ge1", "p_ge2", "p_death") %in% names(s)))
  expect_equal(rowSums(s[c("p_count_0", "p_count_1", "p_count_2", "p_count_3plus")]), rep(1, 25),
               tolerance = 1e-12)
  expect_true(all(s$p_ge1 >= s$p_ge2))
  expect_true(all(unlist(s[-1]) >= 0 & unlist(s[c("p_count_0", "p_ge1", "p_ge2", "p_death")]) <= 1))
  expect_equal(dim(res$counts), c(25, 150))
  expect_equal(s$mean_count, unname(rowMeans(res$counts)))
  expect_output(print(res), "Predicted recurrences")
})

test_that("with an endless stay only the first episode counts and it matches predictRisk()", {
  res <- predict_counts(n_sim = 1500, stay_lengths = 1e6, seed = 2)
  expect_true(all(res$counts <= 1))
  risk <- predictRisk(fx$open, fx$long_fit, fx$surv_fit, prediction_time = 0, horizon = 3,
                      time_variable = "time_in_episode", survival_variable_all = list(),
                      survival_trans_function = list(), bandcount1 = 30, bandcount2 = 60)
  expect_lt(mean(abs(res$summary$p_ge1 - unname(risk$risk_prob_1))), 0.03)
  expect_lt(mean(abs(res$summary$p_death - unname(risk$risk_prob_2))), 0.03)
})

test_that("a longer window gives more recurrences, and later episodes add to the first", {
  short <- predict_counts(horizon = 1, n_sim = 300, stay_lengths = fx$stays, seed = 3, max_stages = 40)
  long  <- predict_counts(horizon = 6, n_sim = 300, stay_lengths = fx$stays, seed = 3, max_stages = 40)
  expect_gt(mean(long$summary$mean_count), mean(short$summary$mean_count))
  expect_gt(mean(long$summary$p_ge2), mean(short$summary$p_ge2))
  expect_gt(long$stages, 1)
  one_episode <- predict_counts(horizon = 6, n_sim = 300, stay_lengths = 1e6, seed = 3)
  expect_gt(mean(long$summary$mean_count), mean(one_episode$summary$mean_count))
})

test_that("max_stages stops runaway paths with a warning and flags the result", {
  expect_warning(res <- predict_counts(horizon = 6, n_sim = 100, stay_lengths = fx$stays, seed = 6,
                                       max_stages = 2), "max_stages")
  expect_true(res$capped)
  ok <- predict_counts(horizon = 1, n_sim = 100, stay_lengths = fx$stays, seed = 6, max_stages = 40)
  expect_false(ok$capped)
})

test_that("seed makes draws reproducible and leaves the caller's random stream alone", {
  set.seed(99); before <- .Random.seed
  a <- predict_counts(n_sim = 80, stay_lengths = fx$stays, seed = 5)
  expect_identical(.Random.seed, before)
  b <- predict_counts(n_sim = 80, stay_lengths = fx$stays, seed = 5)
  expect_identical(a$counts, b$counts)
})

test_that("next_covariates is called for each later episode, with the episode offset", {
  seen <- integer(0)
  hook <- function(row, k) { seen <<- c(seen, k); row$episode <- row$episode + k; row }
  predict_counts(horizon = 6, n_sim = 100, stay_lengths = fx$stays, seed = 4, next_covariates = hook,
                 max_stages = 40)
  expect_true(length(seen) > 0)
  expect_equal(min(seen), 1)
})

test_that("predictRecurrentCount rejects bad inputs", {
  expect_error(predict_counts(stay_lengths = -1), "stay_lengths")
  expect_error(predict_counts(stay_lengths = numeric(0)), "stay_lengths")
  expect_error(predict_counts(horizon = 0), "`horizon`")
  expect_error(predict_counts(next_covariates = 3), "next_covariates")
  expect_error(predictRecurrentCount(fx$open, fx$surv_fit, fx$long_fit, 0, 1, "time_in_episode"),
               "must be the output of")
})

test_that("the default episode update never goes past the largest episode the model was fit on", {
  cap <- training_episode_max(fx$surv_fit, "episode")
  expect_true(is.finite(cap) && cap >= 1)
  upd <- default_next_covariates("episode", cap)
  row <- data.frame(episode = 2, age = 70)
  expect_equal(upd(row, 1)$episode, min(3, cap))
  expect_equal(upd(row, 500)$episode, cap)               # clamped, not 502
  expect_equal(upd(row, 1)$age, 70)                      # other covariates untouched
  # a model without the counter has nothing to clamp to
  expect_equal(default_next_covariates("episode", Inf)(row, 5)$episode, 7)
  expect_equal(upd(data.frame(age = 70), 3)$age, 70)     # no counter column: left alone
})

test_that("a patient with more than one episode in data_predict is refused rather than merged", {
  two <- fx$open[fx$open$id == unique(fx$open$id)[1], ]
  two2 <- two; two2$episode <- two2$episode + 1
  expect_error(predictRecurrentCount(rbind(fx$open, two2), fx$long_fit, fx$surv_fit, prediction_time = 0,
                                     horizon = 1, time_variable = "time_in_episode"),
               "more than one `episode`")
})

test_that("a patient who is not at risk at prediction_time is an error, not a silent zero", {
  bad <- fx$open
  bad$gap[bad$id == unique(bad$id)[1]] <- -1     # a known event time below prediction_time = 0
  expect_error(predictRecurrentCount(bad, fx$long_fit, fx$surv_fit, prediction_time = 0, horizon = 1,
                                     time_variable = "time_in_episode", n_sim = 20),
               "No prediction could be made for 1 of 25")
})

test_that("a landmark partway through the episode (prediction_time > 0) agrees with predictRisk()", {
  skip_if(nrow(fx$mid) == 0, "no episode in the fixture lasted past 0.5")
  res <- predictRecurrentCount(fx$mid, fx$long_fit, fx$surv_fit, prediction_time = 0.5, horizon = 1,
                               time_variable = "time_in_episode", stay_lengths = 1e6, n_sim = 2000, seed = 8)
  risk <- predictRisk(fx$mid, fx$long_fit, fx$surv_fit, prediction_time = 0.5, horizon = 1,
                      time_variable = "time_in_episode", survival_variable_all = list(),
                      survival_trans_function = list(), bandcount1 = 30, bandcount2 = 60)
  key <- as.character(res$summary[[1]])
  expect_lt(mean(abs(res$summary$p_ge1 - unname(risk$risk_prob_1[key]))), 0.03)
  expect_lt(mean(abs(res$summary$p_death - unname(risk$risk_prob_2[key]))), 0.03)
})

test_that("two biomarkers measured at the same times work (one data frame serves both)", {
  res <- predictRecurrentCount(fx$open, fx$long_fit2, fx$surv_fit, prediction_time = 0, horizon = 2,
                               time_variable = "time_in_episode", stay_lengths = fx$stays, n_sim = 60, seed = 2,
                               max_stages = 40)
  expect_equal(nrow(res$summary), 25)
  expect_true(all(is.finite(res$summary$mean_count)))
})
