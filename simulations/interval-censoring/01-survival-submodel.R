# 01. Interval-censored survival sub-model on its own.
# Does survivalSub(Surv(L, R, type = "interval2") ~ x) recover the
# coefficient (with honest standard errors) and the baseline cumulative
# hazard, with the spline and the piecewise-constant baseline?
source(file.path("simulations", "interval-censoring", "R", "setup.R"))
N_REPS <- n_reps_for(100)
times <- c(1, 3, 6)

scenarios <- expand.grid(hazard = c("weibull", "gompertz"), gap = c("1-3", "2-6"), stringsAsFactors = FALSE)
text <- sprintf("01-survival-submodel: %d replicates, n = 300, true coefficient 0.5", N_REPS)
results <- list()
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  gap <- as.numeric(strsplit(sc$gap, "-")[[1]])
  reps <- run_reps(function(r) {
    d <- simulate_ic(300, r, sc$hazard, gap)
    sapply(c(spline = "spline", piecewise = "piecewise"), function(b) {
      ic <- survivalSub(d$surv, IC_FORMULA, NULL, event_time = "Tev", baseline = b)$ic_fit
      c(coef = unname(ic$coefficients), se = sqrt(ic$var[1, 1]),
        stats::setNames(cumulative_baseline_at(ic$cum_basehaz, times) / H0_TRUE[[sc$hazard]](times) - 1,
                        paste0("H0_relerr_t", times)))
    })
  })
  a <- simplify2array(reps)
  tab <- rbind(mean_coef = apply(a["coef", , ], 1, mean),
               empirical_sd = apply(a["coef", , ], 1, stats::sd),
               mean_model_se = apply(a["se", , ], 1, mean),
               apply(a[grep("H0_relerr", dimnames(a)[[1]]), , , drop = FALSE], 1:2, mean))
  results[[paste(sc$hazard, sc$gap)]] <- a
  text <- c(text, capture_table(sprintf("\n== %s baseline, visits every %s yr (H0 rows: mean relative error)",
                                        sc$hazard, sc$gap), round(tab, 3)))
}
save_result("01-survival-submodel", results, text)
