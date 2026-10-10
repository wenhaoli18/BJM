# 02. Longitudinal sub-model f(Y | T) when T is interval-censored.
# Paired difference from the fit with the true event times ("oracle"), for
#   midpoint     T set to the interval midpoint
#   spline       fitIntervalBJM() with the default spline baseline
#   piecewise    fitIntervalBJM(baseline = "piecewise")
#   true f(T)    the same stochastic EM, imputing from the true f(T)
# True values: Tev -0.2, (Intercept) 1.
source(file.path("simulations", "interval-censoring", "R", "setup.R"))
N_REPS <- n_reps_for(60)
keep <- c("(Intercept)", "Tev")

scenarios <- expand.grid(hazard = c("weibull", "gompertz"), gap = c("1-3", "2-6"), stringsAsFactors = FALSE)
text <- sprintf("02-imputation-vs-midpoint: %d replicates, n = 300; mean (Monte Carlo SE) of the difference from the true-T fit", N_REPS)
results <- list()
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  gap <- as.numeric(strsplit(sc$gap, "-")[[1]])
  reps <- run_reps(function(r) {
    d <- simulate_ic(300, r, sc$hazard, gap)
    oracle <- fit_fixef(d, true_T(d))[keep]
    est <- rbind(
      midpoint = fit_fixef(d, midpoint_T(d))[keep],
      spline = fit_interval_bjm(d, r)$fixef[keep],
      piecewise = fit_interval_bjm(d, r, baseline = "piecewise")$fixef[keep],
      `true f(T)` = colMeans(sem_chain(d, true_survival_fit(d, sc$hazard))[6:15, keep]))
    sweep(est, 2, oracle)
  })
  results[[paste(sc$hazard, sc$gap)]] <- reps
  text <- c(text, capture_table(sprintf("\n== %s baseline, visits every %s yr", sc$hazard, sc$gap),
                                summarise_reps(reps)))
}
save_result("02-imputation-vs-midpoint", results, text)
