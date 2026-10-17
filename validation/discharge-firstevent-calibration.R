## Scenario A: post-discharge "readmission vs death within 30 days", ONE index
## discharge per patient, with the current, unmodified BJM pipeline.
##
## Time origin = discharge (months). Biomarker is measured before discharge, so
## its times are negative. A new patient's event time is NA (unknown), which
## predictRisk() treats as "still at risk". No recurrent-event machinery is
## involved: this is ordinary BJM with competing risks (first readmission vs
## death), and it checks that the pipeline is accurate in this setting and how
## much it adds over a simple baseline.
##
## Test design (the point of this rewrite): every test patient is scored at
## discharge (time 0), so there is no landmark selection, and the answer key is
## each patient's TRUE, uncensored event time/type -- administrative censoring
## only affects what the training data can see, never the answer key.
##
## Methods compared on the same test patients:
##   BJM       predictRisk() on the fitted survivalSub()/longitudinalSub()
##   LR-cov    logistic regression on age + sex only
##   LR-last   logistic regression on age + sex + the biomarker's last value
## Metrics, per seed: mean predicted vs observed, calibration intercept/slope
## (logit scale; 0/1 is perfect), Brier score, AUC.
## Run from the package root:  Rscript validation/discharge-firstevent-calibration.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

seeds     <- 1:10
n_train   <- 1000
n_test    <- 2000
horizon   <- 1          # 1 month ~ 30 days
cens_time <- 6          # administrative censoring (training only)
n_cores   <- max(1, min(4, parallel::detectCores() - 1))

sim_patients <- function(n, id0 = 0) {
  id <- id0 + seq_len(n)
  age <- rnorm(n, 70, 8); sex <- rbinom(n, 1, 0.5); age_c <- (age - 70) / 10
  Tev <- (-log(runif(n)) / (0.25 * exp(0.3 * age_c + 0.3 * sex)))^(1 / 1.2)
  type <- rbinom(n, 1, plogis(-2 + 0.5 * Tev + 0.3 * age_c))   # 0 readmission, 1 death
  list(id = id, age = age, sex = sex, Tev = Tev, type = type,
       u0 = rnorm(n, 0, 0.5), u1 = rnorm(n, 0, 0.3))
}

bio_rows <- function(p, tvec = c(-0.4, -0.2, 0)) {   # pre-discharge: negative times
  do.call(rbind, lapply(seq_along(p$id), function(i) {
    age_c <- (p$age[i] - 70) / 10
    mu <- 5 + p$u0[i] + (-0.3 + p$u1[i]) * tvec + 0.4 * p$Tev[i] - 0.5 * p$Tev[i] * tvec +
      0.05 * age_c + 0.2 * p$sex[i]
    data.frame(id = p$id[i], t = tvec, bio = rnorm(length(tvec), mu, 0.4),
               age = p$age[i], sex = p$sex[i])
  }))
}

auc <- function(p, y) {
  if (sum(y) == 0 || sum(!y) == 0) return(NA_real_)
  r <- rank(p); (sum(r[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
calib <- function(p, y) {
  p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
  cf <- tryCatch(coef(glm(y ~ qlogis(p), family = binomial)), error = function(e) c(NA, NA))
  c(mean_pred = mean(p), mean_obs = mean(y), intercept = unname(cf[1]), slope = unname(cf[2]),
    brier = mean((p - y)^2), auc = auc(p, y))
}

one_seed <- function(seed) {
  set.seed(seed)
  tr <- sim_patients(n_train)
  tr$anyev <- as.integer(tr$Tev <= cens_time); tr$Tobs <- pmin(tr$Tev, cens_time)
  surv_df <- data.frame(id = tr$id, Tobs = tr$Tobs, anyev = tr$anyev,
                        type = ifelse(tr$anyev == 1, tr$type, NA), age = tr$age, sex = tr$sex)
  long_tr <- bio_rows(tr)
  long_tr$Tobs  <- tr$Tobs[match(long_tr$id, tr$id)]
  long_tr$anyev <- tr$anyev[match(long_tr$id, tr$id)]
  long_tr$id <- factor(long_tr$id)

  surv_fit <- survivalSub(surv_df, Surv(Tobs, anyev) ~ age + sex, type ~ Tobs + age + sex)
  long_fit <- longitudinalSub(long_tr[long_tr$anyev == 1, ],
                              bio ~ t + age + sex + Tobs + Tobs:t, ~ t | id)

  te <- sim_patients(n_test, id0 = 100000)
  long_te <- bio_rows(te); long_te$Tobs <- NA_real_; long_te$id <- factor(long_te$id)
  # the event-type column must exist in the prediction data (predictRisk() builds
  # the type model's design matrix from it); NA = not yet known
  long_te$type <- NA_integer_
  risk <- predictRisk(long_te, long_fit, surv_fit, prediction_time = 0, horizon = horizon,
                      time_variable = "t", survival_variable_all = list(),
                      survival_trans_function = list(), bandcount1 = 20, bandcount2 = 20)
  y_readm <- te$Tev <= horizon & te$type == 0
  y_death <- te$Tev <= horizon & te$type == 1

  # baselines: training outcomes within the horizon are fully observed
  # (censoring at cens_time > horizon), so plain logistic regression applies
  tr_last <- long_tr[long_tr$t == 0, ]
  tr_d <- data.frame(age = tr$age, sex = tr$sex, last = tr_last$bio[match(tr$id, as.integer(as.character(tr_last$id)))],
                     readm = tr$Tev <= horizon & tr$type == 0, death = tr$Tev <= horizon & tr$type == 1)
  te_last <- long_te[long_te$t == 0, ]
  te_d <- data.frame(age = te$age, sex = te$sex,
                     last = te_last$bio[match(te$id, as.integer(as.character(te_last$id)))])
  pr <- function(f, y) predict(glm(f, data = cbind(tr_d, y = tr_d[[y]]), family = binomial), te_d, type = "response")
  pred <- list(
    BJM       = list(readm = risk$risk_prob_1, death = risk$risk_prob_2),
    `LR-cov`  = list(readm = pr(y ~ age + sex, "readm"),        death = pr(y ~ age + sex, "death")),
    `LR-last` = list(readm = pr(y ~ age + sex + last, "readm"), death = pr(y ~ age + sex + last, "death")))
  out <- list()
  for (m in names(pred)) {
    out[[paste(m, "readmission")]] <- calib(unname(pred[[m]]$readm), y_readm)
    out[[paste(m, "death")]]       <- calib(unname(pred[[m]]$death), y_death)
  }
  res <- do.call(rbind, out)
  data.frame(seed = seed, method = sub(" .*", "", rownames(res)),
             outcome = sub("^[^ ]+ ", "", rownames(res)), res, row.names = NULL)
}

t0 <- Sys.time()
all <- parallel::mclapply(seeds, function(s) tryCatch(one_seed(s), error = function(e) {
  message("seed ", s, " failed: ", conditionMessage(e)); NULL }), mc.cores = n_cores)
all <- do.call(rbind, all)
saveRDS(all, "validation/discharge_firstevent_results.rds")
cat(sprintf("\n%d/%d seeds completed in %.1f min (%d cores)\n",
            length(unique(all$seed)), length(seeds), as.numeric(difftime(Sys.time(), t0, units = "mins")), n_cores))

agg <- aggregate(cbind(mean_pred, mean_obs, intercept, slope, brier, auc) ~ method + outcome,
                 data = all, FUN = function(x) mean(x, na.rm = TRUE))
sds <- aggregate(cbind(slope, auc) ~ method + outcome, data = all, FUN = function(x) sd(x, na.rm = TRUE))
names(sds)[3:4] <- c("slope_sd", "auc_sd")
agg <- merge(agg, sds)
agg[, 3:ncol(agg)] <- round(agg[, 3:ncol(agg)], 3)
cat("\nMean over seeds (slope_sd / auc_sd = spread across seeds):\n")
print(agg[order(agg$outcome, agg$method), ], row.names = FALSE)
