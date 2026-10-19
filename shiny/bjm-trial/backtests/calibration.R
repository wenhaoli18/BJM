# Calibration of BJM on pbc3, by quintile of predicted risk: the 6-year
# risk from entry predicted by BJM (fit on all of pbc3) against what
# happened (Kaplan-Meier), and the C-index.
# Usage: Rscript calibration.R
suppressMessages({library(BJM); library(survival)})
data(pbc3); pbc3$id <- as.character(pbc3$id)
base <- pbc3[!duplicated(pbc3$id), c("id", "years", "status3", "age", "sex")]
trans <- survivalTrans(c(1, 3, 5, 7))
bv <- c("serBilir", "albumin", "prothrombin")
lf <- setNames(lapply(bv, function(b) as.formula(paste(b, "~ year + age + sex + years + years * year"))), bv)
lr <- setNames(rep(list(~ year | id), 3), bv)
fs <- survivalSub(base, Surv(years, status3) ~ age + sex, NULL)
fl <- suppressWarnings(longitudinalSub(pbc3[pbc3$status3 == 1, ], lf, lr))
H <- 6
pr <- vapply(base$id, function(i) { h <- pbc3[pbc3$id == i & pbc3$year <= 0, c("id","year","age","sex",bv)]; h$years <- NA
  unname(predictRisk(h, fl, fs, prediction_time = 0, horizon = H, time_variable = "year", trans$survival_variable_all, trans$survival_trans_function, bandcount1 = 10, bandcount2 = 10)$risk_prob_1) }, numeric(1))
g <- cut(rank(pr), 5, labels = paste0("Q", 1:5))
res <- do.call(rbind, lapply(split(seq_along(pr), g), function(ix) {
  km <- survfit(Surv(pmin(base$years[ix], H), base$status3[ix] == 1 & base$years[ix] <= H) ~ 1)
  data.frame(n = length(ix), mean_pred = mean(pr[ix]), observed = 1 - summary(km, times = H, extend = TRUE)$surv,
             share_censored_early = mean(base$status3[ix] == 0 & base$years[ix] < H))
}))
print(round(res, 3))
cat("C-index:", round(concordance(Surv(pmin(base$years, H), base$status3 == 1 & base$years <= H) ~ I(-pr))$concordance, 3), "\n")
