## Number of recurrences in a window: predictRecurrentCount() vs baselines.
##
## Same simulated hospitalization process as recurrent-discharge-calibration.R
## (episodes start at discharge; plain and frailty data-generating processes).
## Every discharge a test patient reaches is a landmark; the answer key is the
## TRUE number of readmissions in (discharge, discharge + horizon], counted from
## the uncensored simulated process, and whether the patient died in it.
##
## Methods:
##   BJM-roll      predictRecurrentCount(): first episode given the patient's
##                 biomarkers, later episodes from the model's prior given covariates
##   BJM-roll-nobio  the same with the biomarkers withheld, so the first episode is
##                 also drawn from the prior: the value of the biomarkers
##   Poisson       log-linear Poisson regression of the count on age, sex, #previous
##                 readmissions and the lab at discharge (training episodes with a
##                 full horizon of potential follow-up, so the count is complete)
## Metrics per seed and subset: mean predicted vs observed count, Spearman
## correlation of predicted mean and observed count, AUC and Brier for
## "at least 1" and "at least 2" recurrences, and the mean death probability.
## Run from the package root:  Rscript validation/recurrent-count-calibration.R

source("validation/recurrent-discharge-calibration.R")   # sim_patient(), auc(), calib()
max_ep <- 12                   # long windows: do not cap the number of episodes

seeds    <- 1:6
n_train  <- 1000
n_test   <- 600                # test patients (each contributes every discharge it reaches)
horizon  <- 6                  # months
n_sim    <- 200
n_cores  <- max(1, min(4, parallel::detectCores() - 1))

count_truth <- function(p, kk, h) {
  tr <- p$truth; s0 <- tr$s[kk]; end <- tr$s + tr$G
  idx <- seq(kk, nrow(tr))
  in_win <- end[idx] <= s0 + h
  c(count = sum(in_win & tr$type[idx] == 0), death = as.integer(any(in_win & tr$type[idx] == 1)))
}

one_run <- function(job) {
  frailty_sd <- job$frailty_sd; seed <- job$seed
  set.seed(seed * 100 + (frailty_sd > 0) + 7)

  tr <- lapply(seq_len(n_train), sim_patient, frailty_sd = frailty_sd)
  ev_tr  <- do.call(rbind, lapply(tr, `[[`, "events"))
  lab_tr <- do.call(rbind, lapply(tr, `[[`, "long"))
  b <- buildRecurrentEpisodes(ev_tr, lab_tr, time_variable_long = "t",
                              restart_variable = "discharge", lookback = Inf)
  b$long_data$id <- factor(b$long_data$id)
  surv_fit <- survivalSub(b$episode_data, Surv(gap, event_any) ~ age + sex + episode,
                          event_type ~ gap + age + sex + episode)
  long_fit <- longitudinalSub(b$long_data[b$long_data$event_any == 1, ],
                              bio ~ time_in_episode + age + sex + gap + gap:time_in_episode,
                              ~ time_in_episode | id)
  stays <- ev_tr$discharge[ev_tr$status == 1] - ev_tr$time[ev_tr$status == 1]
  stays <- stays[!is.na(stays) & stays >= 0]

  ## Poisson baseline: counts in the horizon after each training discharge,
  ## from the patient's later readmission times, on episodes with full follow-up
  ed <- b$episode_data
  at0 <- b$long_data[b$long_data$time_in_episode == 0, ]
  ed$last <- at0$bio[match(paste(ed$id, ed$episode), paste(at0$id, at0$episode))]
  ed$n_in <- NA_integer_
  by_id <- split(ev_tr[ev_tr$status == 1, ], ev_tr$id[ev_tr$status == 1])
  for (i in seq_len(nrow(ed))) {
    r <- by_id[[as.character(ed$id[i])]]
    ed$n_in[i] <- if (is.null(r)) 0L else sum(r$time > ed$s_start[i] & r$time <= ed$s_start[i] + horizon)
  }
  ed <- ed[(ed$fu_end - ed$s_start) >= horizon & !is.na(ed$last), ]
  pois <- glm(n_in ~ age + sex + episode + last, data = ed, family = poisson)

  ## test landmarks
  te <- lapply(seq_len(n_test), function(i) sim_patient(200000 + i, frailty_sd, censor = FALSE))
  hist_ev <- list(); hist_lab <- list(); key <- list(); m <- 0
  for (i in seq_along(te)) {
    p <- te[[i]]
    for (kk in seq_len(nrow(p$truth))) {
      m <- m + 1
      lid <- p$truth$id[kk] * 100 + kk
      sk <- p$truth$s[kk]
      evs <- p$events[seq_len(kk - 1), , drop = FALSE]
      cur <- data.frame(id = p$truth$id[kk], time = sk + 0.001, status = 0L, discharge = NA_real_,
                        age = p$events$age[1], sex = p$events$sex[1], fu_end = Inf, s_start = sk)
      he <- rbind(evs, cur); he$id <- lid
      hl <- p$long[p$long$t <= sk + 1e-9, , drop = FALSE]; hl$id <- lid
      hist_ev[[m]] <- he; hist_lab[[m]] <- hl
      tc <- count_truth(p, kk, horizon)
      key[[m]] <- data.frame(id = lid, k = kk, n_true = tc[["count"]], died = tc[["death"]])
    }
  }
  hist_ev <- do.call(rbind, hist_ev); hist_lab <- do.call(rbind, hist_lab); key <- do.call(rbind, key)
  bt <- buildRecurrentEpisodes(hist_ev, hist_lab, time_variable_long = "t",
                               restart_variable = "discharge", lookback = Inf)
  cur_rows <- bt$long_data[paste(bt$long_data$id, bt$long_data$episode) %in% paste(key$id, key$k), ]
  cur_rows$id <- factor(cur_rows$id)

  roll <- predictRecurrentCount(cur_rows, long_fit, surv_fit, prediction_time = 0, horizon = horizon,
                                time_variable = "time_in_episode", stay_lengths = stays,
                                n_sim = n_sim, bandcount2 = 100, max_stages = 30)
  nobio_rows <- cur_rows; nobio_rows$bio <- NA_real_
  roll0 <- predictRecurrentCount(nobio_rows, long_fit, surv_fit, prediction_time = 0, horizon = horizon,
                                 time_variable = "time_in_episode", stay_lengths = stays,
                                 n_sim = n_sim, bandcount2 = 100, max_stages = 30)

  capped <- roll$capped
  ids <- roll$summary$id
  key <- key[match(as.integer(ids), key$id), ]
  at0t <- cur_rows[cur_rows$time_in_episode == 0, ]
  te_d <- data.frame(age = at0t$age[match(ids, as.character(at0t$id))],
                     sex = at0t$sex[match(ids, as.character(at0t$id))],
                     episode = key$k, last = at0t$bio[match(ids, as.character(at0t$id))])
  mu <- predict(pois, te_d, type = "response")
  pred <- list(
    `BJM-roll` = list(mean = roll$summary$mean_count, ge1 = roll$summary$p_ge1,
                      ge2 = roll$summary$p_ge2, death = roll$summary$p_death),
    `BJM-roll-nobio` = list(mean = roll0$summary$mean_count, ge1 = roll0$summary$p_ge1,
                            ge2 = roll0$summary$p_ge2, death = roll0$summary$p_death),
    Poisson = list(mean = mu, ge1 = 1 - ppois(0, mu), ge2 = 1 - ppois(1, mu), death = rep(NA_real_, length(mu))))

  subsets <- list(all = rep(TRUE, nrow(key)), `1st discharge` = key$k == 1, `2nd+ discharge` = key$k >= 2)
  rows <- list()
  for (mth in names(pred)) for (sn in names(subsets)) {
    sel <- subsets[[sn]]; pr <- pred[[mth]]
    y <- key$n_true[sel]
    pe <- function(p, yy) { r <- calib(p[sel], yy[sel]); r[c("mean_pred", "mean_obs", "slope", "brier", "auc")] }
    g1 <- pe(pr$ge1, key$n_true >= 1); g2 <- pe(pr$ge2, key$n_true >= 2)
    rows[[length(rows) + 1]] <- data.frame(
      dgp = if (frailty_sd > 0) "frailty" else "plain", seed = seed, method = mth, subset = sn,
      n = sum(sel), mean_count_pred = mean(pr$mean[sel]), mean_count_obs = mean(y),
      spearman = suppressWarnings(cor(pr$mean[sel], y, method = "spearman")),
      auc_ge1 = g1[["auc"]], auc_ge2 = g2[["auc"]], brier_ge2 = g2[["brier"]],
      slope_ge2 = g2[["slope"]], p_ge2_pred = g2[["mean_pred"]], p_ge2_obs = g2[["mean_obs"]],
      death_pred = if (all(is.na(pr$death))) NA_real_ else mean(pr$death[sel]),
      death_obs = mean(key$died[sel]), capped = capped)
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0) {
  jobs <- do.call(c, lapply(c(0, 0.7), function(fs) lapply(seeds, function(s) list(frailty_sd = fs, seed = s))))
  t0 <- Sys.time()
  res <- parallel::mclapply(jobs, function(j) tryCatch(one_run(j), error = function(e) {
    message("job dgp=", j$frailty_sd, " seed=", j$seed, " failed: ", conditionMessage(e)); NULL
  }), mc.cores = n_cores)
  res <- do.call(rbind, res)
  saveRDS(res, "validation/recurrent_count_results.rds")
  cat(sprintf("\n%d/%d jobs completed in %.1f min (%d cores)\n", length(unique(paste(res$dgp, res$seed))),
              length(jobs), as.numeric(difftime(Sys.time(), t0, units = "mins")), n_cores))
  num <- c("n", "mean_count_pred", "mean_count_obs", "spearman", "auc_ge1", "auc_ge2", "brier_ge2",
           "slope_ge2", "p_ge2_pred", "p_ge2_obs", "death_pred", "death_obs")
  agg <- aggregate(res[num], by = res[c("dgp", "subset", "method")], FUN = function(x) mean(x, na.rm = TRUE))
  agg[num] <- round(agg[num], 3)
  options(width = 220)
  print(agg[order(agg$dgp, agg$subset, agg$method), ], row.names = FALSE)
}
