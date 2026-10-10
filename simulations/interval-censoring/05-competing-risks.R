# 05. Competing risks, with the event type known once the event is detected.
# Paired difference from the fits with the true event times, for the
# event-type model's T coefficient (true 0.3) and the longitudinal
# sub-model's T and D coefficients (true -0.2 and 0.5): midpoint vs
# fitIntervalBJM(form_conditional_cr = type ~ Tev + x).
source(file.path("simulations", "interval-censoring", "R", "setup.R"))
N_REPS <- n_reps_for(60)

coefs <- function(d, T_values) {
  ev <- observed_events(d)
  ev$Tev <- unname(T_values[as.character(ev$id)])
  g <- stats::coef(glm(type ~ Tev + x, ev, family = binomial))[["Tev"]]
  f <- fit_fixef(d, T_values, LONG_FIXED_CR)
  c(glm_T = g, long_T = f[["Tev"]], long_D = f[["type"]])
}

text <- sprintf("05-competing-risks: %d replicates, n = 300, Weibull; mean (MC SE) of the difference from the true-T fit", N_REPS)
results <- list()
for (gap in list(c(1, 3), c(2, 6))) {
  reps <- run_reps(function(r) {
    d <- simulate_ic(300, r, "weibull", gap, competing = TRUE)
    oracle <- coefs(d, true_T(d))
    fit <- fit_interval_bjm(d, r, LONG_FIXED_CR, form_conditional_cr = type ~ Tev + x)$fit
    imputed <- colMeans(t(sapply(seq_along(fit$long_fit_all_list), function(m) {
      c(glm_T = stats::coef(fit$survival_fit_all_list[[m]]$glm_fit)[["Tev"]],
        long_T = nlme::fixef(fit$long_fit_all_list[[m]]$lfit[[1]])[["Tev"]],
        long_D = nlme::fixef(fit$long_fit_all_list[[m]]$lfit[[1]])[["type"]])
    })))
    rbind(midpoint = coefs(d, midpoint_T(d)) - oracle, fitIntervalBJM = imputed - oracle)
  })
  key <- paste(gap, collapse = "-")
  results[[key]] <- reps
  text <- c(text, capture_table(sprintf("\n== visits every %s yr", key), summarise_reps(reps)))
}
save_result("05-competing-risks", results, text)
