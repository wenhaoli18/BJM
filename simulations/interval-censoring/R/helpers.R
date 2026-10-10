# Shared machinery for the interval-censoring simulations. Sourced by every
# numbered script (via setup.R); uses BJM's internal functions, so the
# package is loaded with pkgload::load_all().

IC_FORMULA <- Surv(L, R, type = "interval2") ~ x
LONG_FIXED <- y ~ year + Tev + x
LONG_FIXED_CR <- y ~ year + Tev + type + x
LONG_RANDOM <- ~ year | id

#' Run `fun(rep)` for rep = 1..n_reps in parallel (forked; macOS/Linux). Each
#' replicate sets its own seed from `rep`, so results do not depend on the
#' number of cores or scheduling. Failed replicates are reported and dropped.
run_reps <- function(fun, n_reps = N_REPS, cores = N_CORES) {
  out <- parallel::mclapply(seq_len(n_reps), function(r) {
    tryCatch(fun(r), error = function(e) structure(conditionMessage(e), class = "rep_error"))
  }, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(out, inherits, logical(1), "rep_error")
  if (any(failed)) {
    message(sprintf("%d of %d replicates failed: %s", sum(failed), n_reps,
                    paste(unique(unlist(out[failed])), collapse = " | ")))
  }
  out[!failed]
}

#' Fixed effects of the longitudinal sub-model fit with given event times
#' (named by subject id), on the subjects with a detected event.
fit_fixef <- function(d, T_values, formula = LONG_FIXED) {
  long <- d$long[as.character(d$long$id) %in% names(T_values), ]
  long$Tev <- unname(T_values[as.character(long$id)])
  fit <- suppressWarnings(longitudinalSub(list(long), formula, LONG_RANDOM))
  nlme::fixef(fit$lfit[[1]])
}

observed_events <- function(d) d$surv[!is.na(d$surv$R), ]
true_T <- function(d) { ev <- observed_events(d); stats::setNames(ev$trueT, ev$id) }
midpoint_T <- function(d) { ev <- observed_events(d); stats::setNames((ev$L + ev$R) / 2, ev$id) }

#' Interval-censored survival sub-model whose baseline and coefficient are
#' replaced by the truth: imputation then uses the true f(T).
true_survival_fit <- function(d, hazard) {
  fit <- survivalSub(d$surv, IC_FORMULA, NULL, event_time = "Tev")
  fit$ic_fit$coefficients[] <- 0.5
  fit$ic_fit$cum_basehaz$hazard <- H0_TRUE[[hazard]](fit$ic_fit$cum_basehaz$time)
  fit
}

#' The stochastic-EM loop of fitIntervalBJM(), with a given survival
#' sub-model and K draws per iteration stacked into one longitudinal fit
#' (K = 1: stochastic EM, as fitIntervalBJM(); K > 1: Monte Carlo EM).
#' Returns the coefficient trace (iterations x fixed effects).
sem_chain <- function(d, survival_fit_all, n_iter = 15, K = 1, n_grid = 40, formula = LONG_FIXED) {
  bounds <- subject_intervals(d$surv, IC_FORMULA, "id")
  bounds <- bounds[is.finite(bounds$R), ]
  long1 <- d$long[as.character(d$long$id) %in% bounds$id, ]
  fit_stack <- function(Tmat) {
    Tmat <- as.matrix(Tmat)
    stacked <- do.call(rbind, lapply(seq_len(ncol(Tmat)), function(k) {
      dk <- long1
      dk$Tev <- unname(Tmat[as.character(dk$id), k])
      dk$id <- if (ncol(Tmat) > 1) paste(dk$id, k, sep = "_") else dk$id
      dk
    }))
    suppressWarnings(longitudinalSub(list(stacked), formula, LONG_RANDOM))
  }
  lf <- fit_stack(stats::setNames((bounds$L + bounds$R) / 2, bounds$id))
  trace <- NULL
  for (it in seq_len(n_iter)) {
    for (attempt in 1:5) {
      Tm <- imputeEventTime(list(long1), lf, survival_fit_all, bounds, "year", n_grid = n_grid, n_draws = K)
      new_fit <- try(fit_stack(Tm), silent = TRUE)
      if (!inherits(new_fit, "try-error")) break
    }
    if (inherits(new_fit, "try-error")) stop("longitudinal fit failed on 5 successive draws")
    lf <- new_fit
    trace <- rbind(trace, nlme::fixef(lf$lfit[[1]]))
  }
  trace
}

#' Coefficients from fitIntervalBJM(), averaged over its saved imputations.
fit_interval_bjm <- function(d, seed, formula = LONG_FIXED, form_conditional_cr = NULL, ...) {
  fit <- suppressWarnings(fitIntervalBJM(d$surv, d$long, IC_FORMULA, "Tev", formula, LONG_RANDOM, "year",
                                         form_conditional_cr = form_conditional_cr,
                                         n_burnin = 5, n_imputations = 5, seed = seed, ...))
  list(fit = fit,
       fixef = colMeans(do.call(rbind, lapply(fit$long_fit_all_list, function(l) nlme::fixef(l$lfit[[1]])))))
}

#' Mean and Monte Carlo SE over replicates of a list of equally shaped
#' numeric matrices/vectors, formatted "mean (se)".
summarise_reps <- function(results, digits = 4) {
  a <- simplify2array(results)
  if (is.null(dim(a))) a <- matrix(a, nrow = 1)
  last <- length(dim(a))
  m <- apply(a, seq_len(last - 1), mean, na.rm = TRUE)
  se <- apply(a, seq_len(last - 1), stats::sd, na.rm = TRUE) / sqrt(dim(a)[last])
  fmt <- sprintf(paste0("%+.", digits, "f (%.", digits, "f)"), m, se)
  if (is.null(dim(m))) noquote(stats::setNames(fmt, names(m))) else
    noquote(matrix(fmt, nrow(m), dimnames = dimnames(m)))
}

#' Save a result object and a printed summary under results/.
save_result <- function(name, results, summary_text) {
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(results = results, n_reps = N_REPS, git_commit = GIT_COMMIT, date = Sys.time()),
          file.path(RESULTS_DIR, paste0(name, ".rds")))
  writeLines(summary_text, file.path(RESULTS_DIR, paste0(name, ".txt")))
  cat(summary_text, sep = "\n")
}

capture_table <- function(title, x) c(title, utils::capture.output(print(x)), "")
