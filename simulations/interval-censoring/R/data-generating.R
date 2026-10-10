# Data-generating model shared by every interval-censoring simulation.
#
# Backward joint model truth:
#   covariate      x ~ N(0, 1)
#   event time     T | x proportional hazards, exp(0.5 x), with baseline
#                  "weibull":  H0(t) = 0.05 t^1.5          (increasing hazard)
#                  "gompertz": H0(t) = 0.08 (exp(0.25 t) - 1)
#   event type     D | T, x ~ Bernoulli(plogis(-1 + 0.3 T + 0.5 x))   (competing = TRUE)
#   biomarker      y(t) | T, D = 1 + 0.3 t - 0.2 T + 0.5 D + 0.4 x + b0 + b1 t + e
#                  b0 ~ N(0, 0.5^2), b1 ~ N(0, slope_sd^2), e ~ N(0, 0.3^2)
#                  (the 0.5 D term only with competing = TRUE)
#   visits         0, then gaps ~ U(gap[1], gap[2]), until the end of follow-up
#                  C ~ U(4, 15); the biomarker is measured and the event
#                  checked at every visit
#   observation    event in (L, R]: L = last negative visit, R = the visit
#                  that detected it; R = NA if not detected by the last visit
#                  (right-censored at L). The event type is known once the
#                  event is detected.

H0_TRUE <- list(
  weibull = function(t) 0.05 * t^1.5,
  gompertz = function(t) 0.08 * (exp(0.25 * t) - 1)
)

draw_event_time <- function(hazard, lp) {
  e <- stats::rexp(length(lp)) / exp(lp)
  switch(hazard,
         weibull = (e / 0.05)^(1 / 1.5),
         gompertz = log(1 + e / 0.08) / 0.25)
}

#' @return list(surv = one row per subject with id, x, L, R, type (observed
#'   event type), trueT, trueD; long = one row per visit with id, year, x,
#'   type, y)
simulate_ic <- function(n, seed, hazard = c("weibull", "gompertz"), gap = c(1, 3),
                        competing = FALSE, slope_sd = 0.1) {
  hazard <- match.arg(hazard)
  set.seed(seed)
  x <- stats::rnorm(n)
  T <- draw_event_time(hazard, 0.5 * x)
  D <- if (competing) stats::rbinom(n, 1, stats::plogis(-1 + 0.3 * T + 0.5 * x)) else rep(NA_real_, n)
  C <- stats::runif(n, 4, 15)
  surv <- data.frame(id = seq_len(n), x = x, L = NA_real_, R = NA_real_, type = NA_real_,
                     trueT = T, trueD = D)
  long <- vector("list", n)
  for (i in seq_len(n)) {
    v <- c(0, cumsum(stats::runif(40, gap[1], gap[2])))
    v <- v[v < C[i]]
    if (T[i] > max(v)) {
      surv$L[i] <- max(v)
      obs <- v
    } else {
      k <- findInterval(T[i], v)
      surv$L[i] <- v[k]
      surv$R[i] <- v[k + 1]
      surv$type[i] <- D[i]
      obs <- v[v <= T[i]]
    }
    b <- stats::rnorm(2, 0, c(0.5, slope_sd))
    y <- 1 + 0.3 * obs - 0.2 * T[i] + (if (competing) 0.5 * D[i] else 0) + 0.4 * x[i] +
      b[1] + b[2] * obs + stats::rnorm(length(obs), 0, 0.3)
    long[[i]] <- data.frame(id = i, year = obs, x = x[i], type = surv$type[i], y = y)
  }
  list(surv = surv, long = do.call(rbind, long))
}
