# 04. The hazard before the first visit.
# Before the earliest interval endpoint the baseline hazard is extrapolated,
# not estimated, and events before the first visit are imputed from that
# extrapolation. How much does it matter, and does anchoring the spline's
# lower tail to the Weibull fit (icphFit(lower_tail = "weibull")) help?
# Paired difference from the true-T fit (Tev), and the relative error of
# H0(1) (before the first visit in the 2-6 yr design).
source(file.path("simulations", "interval-censoring", "R", "setup.R"))
N_REPS <- n_reps_for(60)

scenarios <- expand.grid(hazard = c("weibull", "gompertz"), gap = c("1-3", "2-6"), stringsAsFactors = FALSE)
text <- sprintf("04-early-hazard: %d replicates, n = 300; mean (MC SE)", N_REPS)
results <- list()
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  gap <- as.numeric(strsplit(sc$gap, "-")[[1]])
  reps <- run_reps(function(r) {
    d <- simulate_ic(300, r, sc$hazard, gap)
    oracle <- fit_fixef(d, true_T(d))[["Tev"]]
    fits <- list(piecewise = survivalSub(d$surv, IC_FORMULA, NULL, event_time = "Tev", baseline = "piecewise"),
                 `spline, own tail` = survivalSub(d$surv, IC_FORMULA, NULL, event_time = "Tev"))
    fits$`spline, Weibull tail` <- fits$`spline, own tail`
    fits$`spline, Weibull tail`$ic_fit <- icphFit(IC_FORMULA, d$surv, "Tev", "spline", 3, lower_tail = "weibull")
    ev <- observed_events(d)
    rbind(sapply(fits, function(f) mean(sem_chain(d, f)[6:15, "Tev"]) - oracle),
          sapply(fits, function(f) cumulative_baseline_at(f$ic_fit$cum_basehaz, 1) / H0_TRUE[[sc$hazard]](1) - 1),
          share_before_first_visit = mean(ev$L == 0))
  })
  reps <- lapply(reps, function(m) { rownames(m) <- c("Tev difference", "H0(1) rel. error", "events before 1st visit"); m })
  results[[paste(sc$hazard, sc$gap)]] <- reps
  text <- c(text, capture_table(sprintf("\n== %s baseline, visits every %s yr", sc$hazard, sc$gap),
                                summarise_reps(reps)),
            sprintf("SD over replicates of the Tev difference: %s",
                    paste(names(reps[[1]][1, ]), round(apply(sapply(reps, function(m) m[1, ]), 1, stats::sd), 4),
                          sep = " ", collapse = "; ")))
}
save_result("04-early-hazard", results, text)
