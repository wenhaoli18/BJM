# 06. Landmark AUC and Brier score for interval-censored fits.
# The true values use every at-risk subject's true event time (case:
# s < T <= s + 2, control: T > s + 2; Brier over all subjects at risk).
# Estimates: IPCW over the subjects whose status is certain (as published by
# Yang et al., 2026, and normalised by the sum of the weights) and the
# model-based version (performancePlot()'s default). Also with a
# misspecified survival sub-model (x left out).
source(file.path("simulations", "interval-censoring", "R", "setup.R"))
N_REPS <- n_reps_for(40)
landmarks <- c(2, 4)
horizon <- 2

auc_of <- function(r, case, ctrl) {
  if (!any(case) || !any(ctrl)) return(NA_real_)
  mean(outer(r[case], r[ctrl], function(a, b) (a > b) + 0.5 * (a == b)))
}

one <- function(r, gap, surv_formula) {
  d <- simulate_ic(300, r, "weibull", gap)
  fit <- suppressWarnings(fitIntervalBJM(d$surv, d$long, surv_formula, "Tev", LONG_FIXED, LONG_RANDOM, "year",
                                         n_burnin = 3, n_imputations = 1, seed = r))
  ev <- merge(d$long, d$surv[, c("id", "L", "R")], by = "id")
  ev$type <- NULL
  preds <- landmark_predictions(list(ev), fit$long_fit_all, fit$survival_fit_all, landmarks, horizon,
                                "year", NULL, NULL, 20, 40, posterior = TRUE)
  do.call(rbind, lapply(preds, function(lp) {
    o <- lp$outcome
    risk <- lp$risks[[1]]
    s <- lp$landmark
    e <- s + horizon
    Tt <- d$surv$trueT[match(o$id, d$surv$id)]
    case <- Tt > s & Tt <= e
    ctrl <- Tt > e
    ip <- performance_metrics_interval_ipcw(risk, o$L, o$R, o$followup, o$cause, 1, s, horizon)
    certain <- (is.finite(o$R) & o$L >= s & o$R <= e) | o$L >= e
    # the same IPCW weights, for the weight-normalised Brier score
    cens <- survfit(Surv(o$followup, as.numeric(!is.finite(o$R))) ~ 1)
    G <- stats::stepfun(cens$time, c(1, cens$surv))
    det <- is.finite(o$R) & o$L >= s & o$R <= e
    w <- numeric(length(risk))
    w[det] <- 1 / G(o$R[det] - 1e-8)
    w[o$L >= e] <- 1 / G(e)
    probs <- t(vapply(seq_len(nrow(o)), function(i)
      window_status_probs(lp$posterior[[o$id[i]]], o$L[i], o$R[i], o$cause[i], s, e, 1), numeric(2)))
    mb <- performance_metrics_interval_model(risk, probs[, 1], probs[, 2])
    truth <- c(AUC = auc_of(risk, case, ctrl), Brier = mean((case - risk)^2))
    rbind(truth = truth,
          ipcw = c(ip$auc, ip$brier) - truth,
          `ipcw, weight-normalised` = c(ip$auc, sum(w * (det - risk)^2) / sum(w)) - truth,
          model = c(mb$auc, mb$brier) - truth,
          `share certain` = mean(certain))
  }))
}

designs <- list(
  list(label = "visits every 1-3 yr", gap = c(1, 3), formula = IC_FORMULA),
  list(label = "visits every 0.25-0.75 yr", gap = c(0.25, 0.75), formula = IC_FORMULA),
  list(label = "visits every 1-3 yr, survival sub-model without x", gap = c(1, 3),
       formula = Surv(L, R, type = "interval2") ~ 1))
text <- sprintf(paste0("06-performance-measures: %d replicates, n = 300, Weibull, %g-yr window; ",
                       "'truth' rows give the true value, the others the estimate minus the truth; mean (MC SE)"),
                N_REPS, horizon)
results <- list()
for (des in designs) {
  reps <- run_reps(function(r) one(r, des$gap, des$formula))
  results[[des$label]] <- reps
  for (j in seq_along(landmarks)) {
    rows <- (j - 1) * 5 + 1:5
    text <- c(text, capture_table(sprintf("\n== %s, landmark %g", des$label, landmarks[j]),
                                  summarise_reps(lapply(reps, function(m) m[rows, ]))))
  }
}
save_result("06-performance-measures", results, text)
