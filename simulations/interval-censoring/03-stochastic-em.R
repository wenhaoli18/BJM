# 03. Properties of the stochastic EM in fitIntervalBJM().
#  (a) Burn-in: how fast does the chain leave the midpoint start?
#  (b) Is the imputation distribution right? Imputing from the true f(T)
#      AND the true f(Y | T) parameters ("full oracle"), with no EM, must be
#      unbiased; compared with the stochastic EM using the true f(T).
#  (c) Monte Carlo EM: K draws per iteration stacked into one fit (K = 1 is
#      the stochastic EM of fitIntervalBJM()).
# Gompertz baseline, visits every 2-6 years: the design where a small
# stochastic-EM bias was seen. Paired difference from the true-T fit, Tev.
source(file.path("simulations", "interval-censoring", "R", "setup.R"))
N_REPS <- n_reps_for(200)
hazard <- "gompertz"
gap <- c(2, 6)

reps <- run_reps(function(r) {
  d <- simulate_ic(300, r, hazard, gap)
  oracle <- fit_fixef(d, true_T(d))[["Tev"]]
  sf_true <- true_survival_fit(d, hazard)

  # (a) burn-in, estimated survival sub-model
  sf_est <- survivalSub(d$surv, IC_FORMULA, NULL, event_time = "Tev")
  burn <- sem_chain(d, sf_est, n_iter = 20)[, "Tev"] - oracle

  # (b) full-oracle imputation: true f(T) and true f(Y | T) parameters
  bounds <- subject_intervals(d$surv, IC_FORMULA, "id")
  bounds <- bounds[is.finite(bounds$R), ]
  long1 <- d$long[as.character(d$long$id) %in% bounds$id, ]
  lt <- suppressWarnings(longitudinalSub(list(transform(long1, Tev = unname(true_T(d)[as.character(id)]))),
                                         LONG_FIXED, LONG_RANDOM))
  lt$lfit[[1]]$coefficients$fixed[] <- c(1, 0.3, -0.2, 0.4)
  lt$lfit[[1]]$sigma <- 0.3
  lt$Sigma_fit[] <- diag(c(0.5, 0.1)^2)
  full_oracle <- mean(sapply(1:8, function(m) {
    Tm <- imputeEventTime(list(long1), lt, sf_true, bounds, "year", n_grid = 40)
    fit_fixef(d, Tm)[["Tev"]]
  })) - oracle

  # (c) stochastic EM (K = 1) and Monte Carlo EM (K = 5, 20), true f(T)
  em <- sapply(c(K1 = 1, K5 = 5, K20 = 20), function(K) {
    mean(sem_chain(d, sf_true, n_iter = 12, K = K)[5:12, "Tev"]) - oracle
  })
  list(burn = burn, other = c(full_oracle = full_oracle, em))
})

burn <- summarise_reps(lapply(reps, `[[`, "burn"))
names(burn) <- paste0("iter", seq_along(burn))
text <- c(sprintf("03-stochastic-em: %d replicates, n = 300, Gompertz, visits every 2-6 yr; Tev minus the true-T fit, mean (MC SE)", N_REPS),
          capture_table("\n(a) stochastic EM from the midpoint start, by iteration", burn),
          capture_table("(b, c) full-oracle imputation vs EM with the true f(T), K draws per iteration",
                        summarise_reps(lapply(reps, `[[`, "other"))))
save_result("03-stochastic-em", reps, text)
