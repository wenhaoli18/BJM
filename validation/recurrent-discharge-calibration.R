## Scenario B: SAME patient, SEVERAL hospitalizations -- every discharge is a
## starting point. Each discharge begins an episode (clock 0 = discharge) that
## ends at the next readmission or death. Predict, at discharge, whether the
## next event within 30 days is a readmission or a death.
##
## Uses buildRecurrentEpisodes(restart_variable = , lookback = Inf), then the
## unmodified survivalSub()/longitudinalSub()/predictRisk().
##
## Test design: EVERY discharge a test patient actually reaches is a landmark
## (each with its own pseudo-id). Whether a patient reaches the k-th discharge
## depends only on events before it, so there is no selection on the outcome;
## the answer key is the true, uncensored length/type of that episode.
## (An earlier version drew ONE discharge per patient uniformly from those the
## patient eventually reaches -- but how many discharges a patient reaches
## depends on how the later episodes end, e.g. death ends the process, so the
## "1st discharge" group was enriched for patients who die in episode 1.
## Administrative censoring only shapes what the training data can see.)
##
## Two data-generating processes:
##   "plain"   hazard depends on age, sex, number of prior hospitalizations --
##             exactly what the fitted models contain (correctly specified)
##   "frailty" additionally a subject-level frailty that multiplies the hazard
##             of EVERY episode and shifts the biomarker level -- BJM's survival
##             sub-model has no frailty, so this is misspecified for BJM (the
##             model class frailtypack is built for)
## Baselines: LR-last = logistic on age, sex, #prior hospitalizations, and the
## lab value at discharge; LR-cov = same without the lab.
## FP-JFM = frailtypack's joint frailty model (recurrent readmissions + death,
## gamma frailty, gap-time scale, spline hazards with fixed smoothing), predicted
## with frailtypack::prediction() at t = 0 (the discharge) given the patient's
## number of previous readmissions. It does NOT use the biomarker, so it answers
## "what does a frailty model on the event history alone achieve". Note its two
## outputs are not defined exactly like BJM's: its recurrent-event probability is
## P(a readmission in the window | alive), not "the next event is a readmission"
## (death is not treated as competing inside the window), and its death
## probability is P(death in the window | alive, history).
## Run from the package root:  Rscript validation/recurrent-discharge-calibration.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
have_fp <- requireNamespace("frailtypack", quietly = TRUE)
if (have_fp) suppressMessages(library(frailtypack))

seeds     <- 1:10
n_train   <- 1000        # patients
n_test    <- 1000        # test PATIENTS (each contributes every discharge it reaches)
horizon   <- 1           # months (~30 days)
max_ep    <- 6
n_cores   <- max(1, min(4, parallel::detectCores() - 1))

## ---- data-generating process -------------------------------------------

sim_patient <- function(id, frailty_sd, censor = TRUE) {
  age <- rnorm(1, 70, 8); sex <- rbinom(1, 1, 0.5); age_c <- (age - 70) / 10
  z  <- if (frailty_sd > 0) rnorm(1, 0, frailty_sd) else 0
  u0 <- rnorm(1, 0, 0.5) + 0.8 * z; u1 <- rnorm(1, 0, 0.3)
  C  <- if (censor) runif(1, 6, 12) else Inf
  stay <- function() max(0.08, rexp(1, 1 / 0.25))   # months in hospital

  s <- 0; prevL <- stay(); k <- 0
  ev <- list(); lab <- list(); truth <- list()
  repeat {
    k <- k + 1
    lin <- 0.3 * age_c + 0.3 * sex + 0.25 * (k - 1) + z
    G <- (-log(runif(1)) / (0.25 * exp(lin)))^(1 / 1.2)
    type <- rbinom(1, 1, plogis(-2 + 0.5 * G + 0.3 * age_c + 0.3 * (k - 1)))  # 1 = death
    end <- s + G
    cens <- end > C
    stop_obs <- if (cens) C else end

    # labs for this episode: two in the preceding stay, one at discharge, then
    # a visit every 0.5 months while still at risk
    tvec <- c(s - prevL * c(0.7, 0.3), s)
    j <- 1
    while (s + 0.5 * j < stop_obs) { tvec <- c(tvec, s + 0.5 * j); j <- j + 1 }
    tie <- tvec - s
    mu <- 5 + u0 + (-0.3 + u1) * tie + 0.4 * G - 0.5 * G * tie + 0.05 * age_c + 0.2 * sex
    lab[[k]] <- data.frame(id = id, t = tvec, bio = rnorm(length(tvec), mu, 0.4),
                           age = age, sex = sex)
    truth[[k]] <- data.frame(id = id, k = k, s = s, G = G, type = type)

    # fu_end (administrative end of follow-up, known for everyone) and s_start
    # (this episode's start) ride along as extra columns, only so the baselines
    # can keep episodes with a full horizon of potential follow-up
    if (cens) { ev[[k]] <- data.frame(id = id, time = C, status = 0L, discharge = NA_real_,
                                      age = age, sex = sex, fu_end = C, s_start = s); break }
    if (type == 1) { ev[[k]] <- data.frame(id = id, time = end, status = 2L, discharge = NA_real_,
                                           age = age, sex = sex, fu_end = C, s_start = s); break }
    L <- stay(); s_next <- end + L
    ev[[k]] <- data.frame(id = id, time = end, status = 1L, discharge = s_next, age = age, sex = sex,
                          fu_end = C, s_start = s)
    if (s_next >= C || k >= max_ep) break      # follow-up ends inside this stay
    prevL <- L; s <- s_next
  }
  list(events = do.call(rbind, ev), long = do.call(rbind, lab), truth = do.call(rbind, truth))
}

## ---- metrics -------------------------------------------------------------

auc <- function(p, y) {
  if (sum(y) == 0 || sum(!y) == 0) return(NA_real_)
  r <- rank(p); (sum(r[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
calib <- function(p, y) {
  p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
  cf <- tryCatch(coef(glm(y ~ qlogis(p), family = binomial)), error = function(e) c(NA, NA))
  c(n = length(y), mean_pred = mean(p), mean_obs = mean(y), intercept = unname(cf[1]),
    slope = unname(cf[2]), brier = mean((p - y)^2), auc = auc(p, y))
}

## ---- one (dgp, seed) -------------------------------------------------------

one_run <- function(job) {
  frailty_sd <- job$frailty_sd; seed <- job$seed
  set.seed(seed * 100 + (frailty_sd > 0))

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

  ed0 <- b$episode_data
  ## baselines: training episodes whose outcome within the horizon is known
  ed <- b$episode_data
  at0 <- b$long_data[b$long_data$time_in_episode == 0, ]
  ed$last <- at0$bio[match(paste(ed$id, ed$episode), paste(at0$id, at0$episode))]
  # keep episodes with a FULL horizon of potential follow-up (administrative
  # end of follow-up is known for everyone): their outcome within the horizon is
  # then known for certain. (Keeping "events + long-enough censored" instead
  # would select on the outcome and inflate the event rate.)
  known <- (ed$fu_end - ed$s_start) >= horizon
  ed <- ed[known & !is.na(ed$last), ]
  hit <- ed$event_any == 1 & ed$gap <= horizon
  ed$y_readm <- hit & ed$event_type == 0
  ed$y_death <- hit & ed$event_type == 1

  ## test landmarks: patient reaches discharges 1..K (uncensored); pick one uniformly
  te <- lapply(seq_len(n_test), function(i) sim_patient(200000 + i, frailty_sd, censor = FALSE))
  hist_ev <- list(); hist_lab <- list(); key <- list(); m <- 0
  for (i in seq_along(te)) {
    p <- te[[i]]
    for (kk in seq_len(nrow(p$truth))) {          # every discharge this patient reaches
      m <- m + 1
      lid <- p$truth$id[kk] * 10 + kk             # unique pseudo-id per landmark
      sk <- p$truth$s[kk]
      evs <- p$events[seq_len(kk - 1), , drop = FALSE]
      cur <- data.frame(id = p$truth$id[kk], time = sk + 0.001, status = 0L, discharge = NA_real_,
                        age = p$events$age[1], sex = p$events$sex[1], fu_end = Inf, s_start = sk)
      he <- rbind(evs, cur); he$id <- lid
      hl <- p$long[p$long$t <= sk + 1e-9, , drop = FALSE]; hl$id <- lid
      hist_ev[[m]] <- he; hist_lab[[m]] <- hl
      key[[m]] <- data.frame(id = lid, k = kk, G = p$truth$G[kk], type = p$truth$type[kk])
    }
  }
  hist_ev <- do.call(rbind, hist_ev); hist_lab <- do.call(rbind, hist_lab); key <- do.call(rbind, key)
  bt <- buildRecurrentEpisodes(hist_ev, hist_lab, time_variable_long = "t",
                               restart_variable = "discharge", lookback = Inf)
  cur_rows <- bt$long_data[paste(bt$long_data$id, bt$long_data$episode) %in% paste(key$id, key$k), ]
  cur_rows$id <- factor(cur_rows$id)

  risk <- predictRisk(cur_rows, long_fit, surv_fit, prediction_time = 0, horizon = horizon,
                      time_variable = "time_in_episode", survival_variable_all = list(),
                      survival_trans_function = list(), bandcount1 = 20, bandcount2 = 20)
  ids <- names(risk$risk_prob_1)
  key <- key[match(as.integer(ids), key$id), ]
  y_readm <- key$G <= horizon & key$type == 0
  y_death <- key$G <= horizon & key$type == 1

  at0t <- cur_rows[cur_rows$time_in_episode == 0, ]
  te_d <- data.frame(age = at0t$age[match(ids, as.character(at0t$id))],
                     sex = at0t$sex[match(ids, as.character(at0t$id))],
                     episode = key$k,
                     last = at0t$bio[match(ids, as.character(at0t$id))])
  fit_lr <- function(y, f) predict(glm(f, data = cbind(ed, y = ed[[y]]), family = binomial), te_d, type = "response")
  fp_pred <- NULL
  if (have_fp) fp_pred <- tryCatch({
    d <- data.frame(id = ed0$id, time = ed0$gap,
                    event = as.integer(ed0$event_any == 1 & !is.na(ed0$event_type) & ed0$event_type == 0),
                    death = as.integer(!is.na(ed0$event_type) & ed0$event_type == 1),
                    age_c = (ed0$age - 70) / 10, sex = ed0$sex)
    sink("/dev/null")
    fp <- try(suppressWarnings(frailtyPenal(
      Surv(time, event) ~ cluster(id) + age_c + sex + terminal(death),
      formula.terminalEvent = ~ age_c + sex, data = d, hazard = "Splines", n.knots = 6,
      kappa = c(10, 10), recurrentAG = FALSE)), silent = TRUE)
    sink()
    if (inherits(fp, "try-error")) stop("frailtypack fit failed")
    # history at each landmark: its previous readmissions' gaps (event = 1), or
    # a single event = 0 row for a patient with none
    ebt <- bt$episode_data
    prev <- ebt[ebt$episode < ave(ebt$episode, ebt$id, FUN = max), ]
    lm_age <- ebt$age[match(key$id, ebt$id)]; lm_sex <- ebt$sex[match(key$id, ebt$id)]
    hist_rows <- data.frame(time = prev$gap, event = 1L, id = prev$id,
                            age_c = (prev$age - 70) / 10, sex = prev$sex)
    none <- setdiff(key$id, prev$id)
    none_rows <- data.frame(time = rep(0.01, length(none)), event = 0L, id = none,
                            age_c = (lm_age[match(none, key$id)] - 70) / 10, sex = lm_sex[match(none, key$id)])
    dp <- rbind(hist_rows, none_rows)
    pr <- suppressWarnings(prediction(fp, dp, t = 0, window = horizon, event = "Both"))
    ord <- match(key$id, as.integer(as.character(pr$group)))
    list(readm = pr$pred1_rec[ord, 1], death = pr$pred1[ord, 1])
  }, error = function(e) { try(sink(), silent = TRUE); message("frailtypack failed: ", conditionMessage(e)); NULL })

  pred <- list(
    BJM       = list(readm = unname(risk$risk_prob_1), death = unname(risk$risk_prob_2)),
    `LR-cov`  = list(readm = fit_lr("y_readm", y ~ age + sex + episode),
                     death = fit_lr("y_death", y ~ age + sex + episode)),
    `LR-last` = list(readm = fit_lr("y_readm", y ~ age + sex + episode + last),
                     death = fit_lr("y_death", y ~ age + sex + episode + last)))
  if (!is.null(fp_pred)) pred[["FP-JFM"]] <- fp_pred

  subsets <- list(all = rep(TRUE, nrow(key)), `1st discharge` = key$k == 1, `2nd+ discharge` = key$k >= 2)
  rows <- list()
  for (m in names(pred)) for (o in c("readmission", "death")) for (sn in names(subsets)) {
    sel <- subsets[[sn]]
    pp <- if (o == "readmission") pred[[m]]$readm else pred[[m]]$death
    yy <- if (o == "readmission") y_readm else y_death
    rows[[length(rows) + 1]] <- data.frame(
      dgp = if (frailty_sd > 0) "frailty" else "plain", seed = seed, method = m, outcome = o,
      subset = sn, t(calib(pp[sel], yy[sel])))
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
  saveRDS(res, "validation/recurrent_discharge_results.rds")
  cat(sprintf("\n%d/%d jobs completed in %.1f min (%d cores)\n", length(unique(paste(res$dgp, res$seed))),
              length(jobs), as.numeric(difftime(Sys.time(), t0, units = "mins")), n_cores))

  agg <- aggregate(cbind(n, mean_pred, mean_obs, intercept, slope, brier, auc) ~ dgp + outcome + subset + method,
                   data = res, FUN = function(x) mean(x, na.rm = TRUE))
  sds <- aggregate(cbind(slope, auc) ~ dgp + outcome + subset + method, data = res,
                   FUN = function(x) sd(x, na.rm = TRUE))
  names(sds)[5:6] <- c("slope_sd", "auc_sd")
  agg <- merge(agg, sds)
  agg[, 5:ncol(agg)] <- round(agg[, 5:ncol(agg)], 3)
  options(width = 200)
  print(agg[order(agg$dgp, agg$outcome, agg$subset, agg$method), ], row.names = FALSE)
}
