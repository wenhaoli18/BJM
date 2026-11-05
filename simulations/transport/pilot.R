.libPaths(c("rlib", .libPaths())); suppressMessages(library(BJM)); data(pbc3)
args <- commandArgs(TRUE); k_true <- as.numeric(args[1]); reps <- as.integer(args[2]); seed0 <- as.integer(args[3])
first <- pbc3[!duplicated(pbc3$id), ]
st <- survivalTrans(c(1,3,5,7)); SV <- st$survival_variable_all; SF <- st$survival_trans_function
fixed <- list(serBilir = serBilir ~ year + age + sex + (years) + (years)*year,
              albumin  = albumin  ~ year + age + sex + (years) + (years)*year)
rand <- list(serBilir = ~ year | id, albumin = ~ year | id)
## generating model: BJM fitted on the real pbc3
src_true <- survivalSub(first, Surv(years, status3) ~ age + sex, NULL)
long_true <- longitudinalSub(list(pbc3[pbc3$status3==1,], pbc3[pbc3$status3==1,]), fixed, rand)
maxfu <- max(first$years)
shifted <- function(fit, k) {   # true model of a population with baseline hazard x k
  bh <- survival::basehaz(fit$coxph_fit, centered = FALSE); bh$hazard <- bh$hazard * k
  fit$transport <- list(method = "truth", cum_basehaz = bh[, c("hazard","time")], lp_slope = 1,
                        max_followup = maxfu); fit }
gen <- function(n, k) {
  b <- first[sample(nrow(first), n, replace = TRUE), c("age", "sex")]
  b$id <- seq_len(n); b$year <- 0; b$serBilir <- NA_real_; b$albumin <- NA_real_; b$years <- NA_real_
  sim <- simulateTrajectory(list(b, b), long_true, shifted(src_true, k), prediction_time = 0,
                            times = seq(0, 14, by = 0.5), time_variable = "year",
                            survival_variable_all = SV, survival_trans_function = SF, n_sim = 1, bandcount2 = 200)
  sim <- as.data.frame(sim)
  sim <- merge(sim, b[, c("id", "age", "sex")], by = "id")
  C <- setNames(runif(n, 4, 14), seq_len(n))
  sim$true_T <- sim$event_time; sim$true_event <- sim$event
  sim$years <- pmin(sim$true_T, C[as.character(sim$id)])
  sim$status3 <- as.numeric(sim$true_event == 1 & sim$true_T <= C[as.character(sim$id)])
  keep <- sim$year < sim$years & !is.na(sim$serBilir) & !is.na(sim$albumin) & (sim$year == 0 | runif(nrow(sim)) < 0.5)
  sim <- sim[keep, ]
  sim[order(sim$id, sim$year), c("id","year","age","sex","serBilir","albumin","years","status3","true_T","true_event")]
}
H <- 2
auc <- function(p, y) { r <- rank(p); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1*(n1+1)/2)/(n1*n0) }
res <- list()
for (rep in seq_len(reps)) {
  set.seed(seed0 + rep)
  S <- gen(400, 1); Tg <- gen(400, k_true)
  S1 <- S[!duplicated(S$id), ]; T1 <- Tg[!duplicated(Tg$id), ]
  src_hat <- survivalSub(S1, Surv(years, status3) ~ age + sex, NULL)
  long_hat <- longitudinalSub(list(S[S$status3==1,], S[S$status3==1,]), fixed, rand)
  lab <- T1[sample(nrow(T1), 100), ]
  hide <- function(d) { d$years <- NA; d$status3 <- NA; d }   # target outcomes unknown
  t0 <- Sys.time()
  fits <- list(
    source_only  = src_hat,
    unlabeled_noref = transportBJM(src_hat, "unlabeled", data_predict_all = hide(Tg), long_fit_all = long_hat,
                       time_variable = "year", survival_variable_all = SV, survival_trans_function = SF),
    unlabeled    = transportBJM(src_hat, "unlabeled", data_predict_all = hide(Tg), long_fit_all = long_hat,
                       time_variable = "year", survival_variable_all = SV, survival_trans_function = SF,
                       data_source = hide(S)),
    baseline100  = transportBJM(src_hat, "baseline", data_target = lab),
    recal100     = transportBJM(src_hat, "recalibrate", data_target = lab),
    local_refit400 = survivalSub(T1, Surv(years, status3) ~ age + sex, NULL))
  t_unl <- Sys.time() - t0
  ## evaluation: target patients in follow-up at L, history up to L, true outcome in (L, L+H]
  for (L in c(0, 1, 3)) {
  at_risk <- T1$id[T1$years > L]
  dp <- Tg[Tg$id %in% at_risk & Tg$year <= L, ]; dp$years <- L + 1e-6  # outcome not used
  y1 <- Tg[!duplicated(Tg$id), ]; y1 <- y1[match(unique(dp$id), y1$id), ]
  y <- y1$true_event == 1 & y1$true_T <= L + H
  for (m in names(fits)) {
    r <- predictRisk(list(dp, dp), long_hat, fits[[m]], prediction_time = L, horizon = H, time_variable = "year",
                     survival_variable_all = SV, survival_trans_function = SF, bandcount1 = 10, bandcount2 = 30)$risk_prob_1
    r <- r[as.character(y1$id)]
    a <- fits[[m]]$transport$alpha
    res[[length(res)+1]] <- data.frame(k = k_true, rep = rep, landmark = L, method = m, n_eval = length(y),
      observed = mean(y), predicted = mean(r), OE = mean(y)/mean(r), brier = mean((r - y)^2), auc = auc(r, y),
      HR_hat = if (is.null(a)) NA else exp(a), alpha_se = if (is.null(a)) NA else fits[[m]]$transport$alpha_se)
  }
  }
  saveRDS(do.call(rbind, res), sprintf("sim/res_k%g_s%d.rds", k_true, seed0))
  cat(sprintf("k=%g rep %d done (transport fits %.1fs)\n", k_true, rep, as.numeric(t_unl, units = "secs")))
}
out <- do.call(rbind, res)
saveRDS(out, sprintf("sim/res_k%g_s%d.rds", k_true, seed0))
