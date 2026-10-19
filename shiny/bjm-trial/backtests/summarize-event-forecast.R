# Summarize a backtest: per method and cut, accuracy of the forecast date of
# the D-th event (in months) and calibration of the 90% interval.
# Usage: Rscript summarize-event-forecast.R <results.rds>
x <- readRDS(commandArgs(TRUE)[1])$results
x$abs_err <- abs(x$err)
s <- do.call(rbind, lapply(split(x, list(x$method, x$frac)), function(d) data.frame(
  method = d$method[1], at = paste0(100 * d$frac[1], "% of D"),
  MAE_months = round(12 * mean(d$abs_err), 1),
  bias_months = round(12 * mean(d$err), 1),
  coverage_90 = round(mean(d$covered), 2),
  width_months = round(12 * mean(d$width), 1),
  never = round(mean(d$p_never), 3))))
s <- s[order(s$at, match(s$method, c("bjm", "wb_hist", "wb_blind"))), ]
rownames(s) <- NULL
cat("replicates:", length(unique(x$rep)), " median D:", median(x$D),
    " median months from cut to D-th event:",
    round(12 * median(x$actual - x$cut)), "\n\n")
print(s)
# how often BJM is closer than each Weibull method
w <- reshape(x[, c("rep", "frac", "method", "abs_err")], idvar = c("rep", "frac"),
             timevar = "method", direction = "wide")
cat("\nshare of forecasts where BJM is closer than:\n")
print(round(sapply(split(w, w$frac), function(d) c(
  wb_hist = mean(d$abs_err.bjm < d$abs_err.wb_hist),
  wb_blind = mean(d$abs_err.bjm < d$abs_err.wb_blind))), 2))
