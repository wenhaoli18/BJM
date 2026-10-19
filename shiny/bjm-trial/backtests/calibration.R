# Calibration of BJM on pbc3, by quintile of predicted risk: the 6-year
# risk from entry predicted by BJM (fit on all of pbc3) against what
# happened (Kaplan-Meier), and the C-index.
# Usage: Rscript calibration.R
suppressMessages({library(BJM); library(survival)})
data(pbc3); pbc3$id <- as.character(pbc3$id)
base <- pbc3[!duplicated(pbc3$id), c("id", "years", "status3", "age", "sex")]
# The longitudinal sub-models' event-time term is capped at T_cap years. A
# term linear in T is extrapolated past the last event time (13.9 years in
# pbc3), which overstates the risk of low- and mid-risk patients
# (backtests/calibration.R); the cap keeps BJM calibrated, and its ranking.
T_cap <- 8
trans <- list(survival_variable_all = list("Tc"),
              survival_trans_function = list(function(x) pmin(x, T_cap)))
bv <- c("serBilir", "albumin", "prothrombin")
lf <- setNames(lapply(bv, function(b) as.formula(paste(b, "~ year + age + sex + Tc + Tc * year"))), bv)
lr <- setNames(rep(list(~ year | id), 3), bv)
fs <- survivalSub(base, Surv(years, status3) ~ age + sex, NULL)
fl <- suppressWarnings(longitudinalSub(transform(pbc3[pbc3$status3 == 1, ], Tc = pmin(years, T_cap)), lf, lr))
H <- 6
pr <- vapply(base$id, function(i) { h <- pbc3[pbc3$id == i & pbc3$year <= 0, c("id","year","age","sex",bv)]; h$years <- NA; h$Tc <- NA
  unname(predictRisk(h, fl, fs, prediction_time = 0, horizon = H, time_variable = "year", trans$survival_variable_all, trans$survival_trans_function, bandcount1 = 10, bandcount2 = 10)$risk_prob_1) }, numeric(1))
g <- cut(rank(pr), 5, labels = paste0("Q", 1:5))
res <- do.call(rbind, lapply(split(seq_along(pr), g), function(ix) {
  km <- survfit(Surv(pmin(base$years[ix], H), base$status3[ix] == 1 & base$years[ix] <= H) ~ 1)
  data.frame(n = length(ix), mean_pred = mean(pr[ix]), observed = 1 - summary(km, times = H, extend = TRUE)$surv,
             share_censored_early = mean(base$status3[ix] == 0 & base$years[ix] < H))
}))
print(round(res, 3))
cat("C-index:", round(concordance(Surv(pmin(base$years, H), base$status3 == 1 & base$years <= H) ~ I(-pr))$concordance, 3), "\n")
