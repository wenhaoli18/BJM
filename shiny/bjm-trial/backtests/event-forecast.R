# Backtest: forecasting the date of the D-th event in an event-driven trial.
#
# Each replicate splits pbc3 at random into a "historical" half (complete
# follow-up, used to fit models, as a sponsor would use a natural-history
# cohort or earlier trials) and a "trial" half. Trial patients get random
# entry dates over an accrual period. At the calendar time when a share of
# the D target events has been observed, each method forecasts when the
# D-th trial event will occur, and the forecast is compared with the date
# it actually occurred.
#
# Methods:
#   bjm       BJM fit on the historical half; each at-risk trial patient's
#             event time is drawn given their biomarker history so far
#   wb_hist   Weibull fit on the historical half, no biomarkers
#   wb_blind  Weibull fit on the blinded trial data so far (usual practice)
#
# All methods share the same dropout model (exponential, from the losses
# observed in the trial so far) and the same model of future enrollment.
#
# Usage: Rscript event-forecast.R <n_reps> <n_sim> <out.rds>
# The app's table: Rscript event-forecast.R 50 500 bt50.rds, then
#   Rscript summarize-event-forecast.R bt50.rds

suppressMessages({library(BJM); library(survival)})
args <- commandArgs(TRUE)
n_reps <- as.integer(args[1]); n_sim <- as.integer(args[2]); out_file <- args[3]

data(pbc3)
bio_vars <- c("serBilir", "albumin", "prothrombin")
base <- pbc3[!duplicated(pbc3$id), c("id", "years", "status3", "age", "sex")]
base$id <- as.character(base$id)
visits <- pbc3[, c("id", "year", "age", "sex", bio_vars, "years", "status3")]
visits$id <- as.character(visits$id)

accrual <- 3                        # years
frac_at_cut <- c(0.3, 0.5, 0.7)     # share of D observed at the cut
# The longitudinal sub-models' event-time term is capped at T_cap years. A
# term linear in T is extrapolated past the last event time (13.9 years in
# pbc3), which overstates the risk of low- and mid-risk patients
# (backtests/calibration.R); the cap keeps BJM calibrated, and its ranking.
T_cap <- 8
trans <- list(survival_variable_all = list("Tc"),
              survival_trans_function = list(function(x) pmin(x, T_cap)))
long_fixed <- setNames(lapply(bio_vars, function(b)
  as.formula(paste(b, "~ year + age + sex + Tc + Tc * year"))), bio_vars)
long_random <- setNames(rep(list(~ year | id), length(bio_vars)), bio_vars)
methods <- c("bjm", "wb_hist", "wb_blind")

summarize <- function(dates, actual) {
  q <- quantile(dates, c(0.05, 0.5, 0.95), names = FALSE)
  c(median = q[2], lo = q[1], hi = q[3], err = q[2] - actual,
    covered = actual >= q[1] && actual <= q[3], width = q[3] - q[1],
    p_never = mean(!is.finite(dates)))
}

# Weibull sampler, conditional on being event-free at t0
weibull_sampler <- function(fit) {
  k <- 1 / fit$scale; lambda <- exp(coef(fit))
  function(t0, n) lambda * ((t0 / lambda)^k - log(runif(n)))^(1 / k)
}

one_cut <- function(trial, cut, actual, D, fits) {
  enrolled <- trial[trial$entry < cut, ]
  n_future <- sum(trial$entry >= cut)
  enrolled$t_cut <- cut - enrolled$entry
  enrolled$t_obs <- pmin(enrolled$years, enrolled$t_cut)
  enrolled$d_obs <- as.integer(enrolled$status3 == 1 & enrolled$years <= enrolled$t_cut)
  # censored before the cut: lost to follow-up, no further events
  enrolled$lost <- enrolled$status3 == 0 & enrolled$years < enrolled$t_cut
  at_risk <- enrolled[enrolled$d_obs == 0 & !enrolled$lost, ]
  observed_dates <- with(enrolled[enrolled$d_obs == 1, ], entry + t_obs)
  drop_rate <- max(sum(enrolled$lost), 0.5) / sum(enrolled$t_obs)
  rate_in <- nrow(enrolled) / min(cut, accrual)

  ## BJM: n_sim x patient matrices of event times (Inf = event after follow-up)
  draw_bjm <- function(history, t0) {
    s <- simulateTrajectory(history, fits$long, fits$surv, prediction_time = t0,
                            times = t0, time_variable = "year",
                            trans$survival_variable_all, trans$survival_trans_function,
                            n_sim = n_sim, seed = sample.int(1e6, 1))
    ifelse(s$event == 1, s$event_time, Inf)   # one row per draw (times = t0)
  }
  bjm_risk <- vapply(seq_len(nrow(at_risk)), function(i) {
    p <- at_risk[i, ]
    h <- visits[visits$id == p$id & visits$year <= p$t_cut, c("id", "year", "age", "sex", bio_vars)]
    h$years <- NA
    h$Tc <- NA
    draw_bjm(h, p$t_cut)
  }, numeric(n_sim))
  bjm_future <- vapply(seq_len(n_future), function(j) {
    donor <- enrolled[sample.int(nrow(enrolled), 1), ]   # covariates not known yet
    h <- data.frame(id = paste0("new", j), year = 0, age = donor$age, sex = donor$sex,
                    serBilir = NA, albumin = NA, prothrombin = NA, years = NA, Tc = NA)
    draw_bjm(h, 0)
  }, numeric(n_sim))

  rw_hist <- weibull_sampler(fits$weibull)
  rw_blind <- weibull_sampler(survreg(Surv(t_obs, d_obs) ~ 1, data = enrolled, dist = "weibull"))

  dates <- sapply(methods, function(m) vapply(seq_len(n_sim), function(s) {
    t_risk <- switch(m, bjm = bjm_risk[s, ],
                     wb_hist = rw_hist(at_risk$t_cut, nrow(at_risk)),
                     wb_blind = rw_blind(at_risk$t_cut, nrow(at_risk)))
    ev_risk <- at_risk$entry + t_risk
    ev_risk[ev_risk > cut + rexp(nrow(at_risk), drop_rate)] <- Inf
    entry_new <- cut + cumsum(rexp(n_future, rate_in))
    t_new <- switch(m, bjm = bjm_future[s, ],
                    wb_hist = rw_hist(0, n_future), wb_blind = rw_blind(0, n_future))
    t_new[t_new > rexp(n_future, drop_rate)] <- Inf
    all <- sort(c(observed_dates, ev_risk, entry_new + t_new))
    if (length(all) >= D) all[D] else Inf
  }, numeric(1)))
  t(apply(dates, 2, summarize, actual = actual))
}

results <- list()
for (r in seq_len(n_reps)) {
  set.seed(1000 + r)
  hist_ids <- sample(base$id, nrow(base) / 2)
  hist <- base[base$id %in% hist_ids, ]
  trial <- base[!base$id %in% hist_ids, ]
  trial$entry <- runif(nrow(trial), 0, accrual)
  D <- round(0.8 * sum(trial$status3 == 1))
  ev_dates <- sort(with(trial[trial$status3 == 1, ], entry + years))
  actual <- ev_dates[D]

  hv <- visits[visits$id %in% hist_ids, ]

  hv$Tc <- pmin(hv$years, T_cap)
  fits <- list(
    surv = survivalSub(hist, Surv(years, status3) ~ age + sex, NULL),
    long = suppressWarnings(longitudinalSub(hv[hv$status3 == 1, ], long_fixed, long_random)),
    weibull = survreg(Surv(years, status3) ~ 1, data = hist, dist = "weibull"))

  for (f in frac_at_cut) {
    cut <- ev_dates[round(f * D)] + 1e-6
    t0 <- Sys.time()
    res <- one_cut(trial, cut, actual, D, fits)
    secs <- as.numeric(Sys.time() - t0, units = "secs")
    results[[length(results) + 1]] <- data.frame(
      rep = r, frac = f, cut = cut, actual = actual, D = D, method = rownames(res), res, secs = secs)
    message(sprintf("rep %d, %.0f%% of D (cut %.2f): actual %.2f | %s (%.0fs)", r, 100 * f, cut, actual,
                    paste(sprintf("%s %.2f [%.2f, %.2f]", rownames(res), res[, "median"], res[, "lo"],
                                  res[, "hi"]), collapse = " | "), secs))
  }
}
results <- do.call(rbind, results)
saveRDS(list(results = results, accrual = accrual, n_sim = n_sim), out_file)
