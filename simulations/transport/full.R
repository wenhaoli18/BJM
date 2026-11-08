## Full simulation study for transportBJM()
## Usage: Rscript sim/full.R <first_rep> <last_rep>
## Each replicate runs every scenario and is saved to sim/full/rep_<r>.rds.
.libPaths(c("rlib", .libPaths())); suppressMessages(library(BJM)); data(pbc3)
args <- commandArgs(TRUE); reps <- seq(as.integer(args[1]), as.integer(args[2]))
dir.create("sim/full", showWarnings = FALSE)

first <- pbc3[!duplicated(pbc3$id), ]
st <- survivalTrans(c(1, 3, 5, 7)); SV <- st$survival_variable_all; SF <- st$survival_trans_function
fixed <- list(serBilir = serBilir ~ year + age + sex + (years) + (years) * year,
              albumin  = albumin  ~ year + age + sex + (years) + (years) * year)
rand <- list(serBilir = ~ year | id, albumin = ~ year | id)
H <- 2; LANDMARKS <- c(0, 2); N_S <- 400; N_EVAL <- 500

## ---- generating models (BJM fitted on the real pbc3) ------------------------
long_true <- longitudinalSub(list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ]), fixed, rand)
surv_true <- list(
  noCR = survivalSub(first, Surv(years, status3) ~ age + sex, NULL),
  CR   = survivalSub(first, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex))
maxfu <- max(first$years)

## a survival sub-model whose baseline cumulative hazard is `cbh * exp(alpha)`
with_alpha <- function(fit, alpha, cbh = BJM:::survival_cum_basehaz(fit)[, c("hazard", "time")],
                       max_followup = BJM:::survival_max_followup(fit)) {
  fit$transport <- NULL
  cbh$hazard <- cbh$hazard * exp(alpha)
  fit$transport <- list(method = "manual", cum_basehaz = cbh, lp_slope = 1, max_followup = max_followup)
  class(fit) <- c("transportBJM", "survivalSub.BJM")
  fit
}
## target longitudinal model with a weaker link between trajectory and event time
weaken_link <- function(long, factor) {
  for (i in seq_along(long$lfit)) {
    b <- long$lfit[[i]]$coefficients$fixed
    j <- grepl("years", names(b))
    long$lfit[[i]]$coefficients$fixed[j] <- b[j] * factor
  }
  long
}

## ---- data generation ---------------------------------------------------------
## n subjects from (surv, long): visits every 0.5 years kept with prob `visit`
## (baseline always), uniform administrative censoring on (4, 14) if `censor`,
## `assay` added to serBilir (a lab calibration offset)
gen <- function(n, surv, long, visit = 0.5, assay = 0, censor = TRUE) {
  b <- first[sample(nrow(first), n, replace = TRUE), c("age", "sex")]
  b$id <- seq_len(n); b$year <- 0; b$serBilir <- NA_real_; b$albumin <- NA_real_; b$years <- NA_real_
  sim <- as.data.frame(simulateTrajectory(list(b, b), long, surv, prediction_time = 0,
                                          times = seq(0, 14, by = 0.5), time_variable = "year",
                                          survival_variable_all = SV, survival_trans_function = SF,
                                          n_sim = 1, bandcount2 = 200))
  sim <- merge(sim, b[, c("id", "age", "sex")], by = "id")
  C <- if (censor) runif(n, 4, 14) else rep(Inf, n)
  C <- C[as.integer(sim$id)]
  sim$true_T <- sim$event_time; sim$true_event <- sim$event
  sim$years <- pmin(sim$true_T, C)
  sim$status3 <- as.numeric(sim$true_event == 1 & sim$true_T <= C)
  has_cr <- "status4" %in% names(sim)
  if (has_cr) {
    sim$true_type <- sim$status4
    sim$status4 <- ifelse(sim$status3 == 1, sim$true_type, NA)
  }
  keep <- sim$year < sim$years & !is.na(sim$serBilir) & !is.na(sim$albumin) &
    (sim$year == 0 | runif(nrow(sim)) < visit)
  sim <- sim[keep, ]
  sim$serBilir <- sim$serBilir + assay
  cols <- c("id", "year", "age", "sex", "serBilir", "albumin", "years", "status3",
            if (has_cr) c("status4", "true_type"), "true_T", "true_event")
  sim[order(sim$id, sim$year), cols]
}
hide <- function(d) { d$years <- NA; d$status3 <- NA; if ("status4" %in% names(d)) d$status4 <- NA; d }
one_row <- function(d) d[!duplicated(d$id), ]

## ---- scenarios -------------------------------------------------------------
sc <- function(name, k, n_T, world = "noCR", visit = 0.5, assay = 0, link = 1)
  data.frame(name, k, n_T, world, visit, assay, link)
scenarios <- rbind(
  sc("A k0.5 n100", 0.5, 100), sc("A k0.5 n400", 0.5, 400),
  sc("A k1 n100", 1, 100),     sc("A k1 n400", 1, 400),
  sc("A k2 n100", 2, 100),     sc("A k2 n400", 2, 400),
  sc("B assay k1", 1, 400, assay = 0.3), sc("B assay k2", 2, 400, assay = 0.3),
  sc("B link k1", 1, 400, link = 0.5),   sc("B link k2", 2, 400, link = 0.5),
  sc("C sparse k2", 2, 400, visit = 0.2), sc("C dense k2", 2, 400, visit = 1),
  sc("D CR k2", 2, 400, world = "CR"))

auc <- function(p, y) {
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(p); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
## labeled proportional shift of the baseline (calibration in the large):
## Poisson MLE of alpha with the source linear predictor and baseline as offset
citl_alpha <- function(src_hat, d1) {
  H0 <- BJM:::cumulative_baseline_at(BJM:::survival_cum_basehaz(src_hat), d1$years)
  e <- sum(H0 * exp(BJM:::survival_lp(src_hat, d1)))
  list(alpha = log(sum(d1$status3) / e), se = 1 / sqrt(sum(d1$status3)))
}
total_risk <- function(r) if (is.null(r$risk_prob_2)) r$risk_prob_1 else r$risk_prob_1 + r$risk_prob_2

## ---- one replicate -------------------------------------------------------------
run_rep <- function(rep) {
  out <- list()
  for (world in unique(scenarios$world)) {
    set.seed(10000 * rep + match(world, unique(scenarios$world)))
    S <- gen(N_S, surv_true[[world]], long_true)
    S1 <- one_row(S)
    cr_form <- if (world == "CR") status4 ~ years + age + sex else NULL
    src_hat <- survivalSub(S1, Surv(years, status3) ~ age + sex, cr_form)
    long_hat <- longitudinalSub(list(S[S$status3 == 1, ], S[S$status3 == 1, ]), fixed, rand)
    a_S <- BJM:::unlabeled_alpha(src_hat, hide(S), long_hat, "year", SV, SF, NULL, 200, NULL,
                                 c(-4, 4), TRUE)
    for (s in which(scenarios$world == world)) {
      p <- scenarios[s, ]
      set.seed(10000 * rep + 100 + s)
      truth_surv <- with_alpha(surv_true[[world]], log(p$k), max_followup = maxfu)
      truth_long <- if (p$link != 1) weaken_link(long_true, p$link) else long_true
      Tg <- gen(p$n_T, truth_surv, truth_long, p$visit, p$assay)
      T1 <- one_row(Tg)
      E <- gen(N_EVAL, truth_surv, truth_long, p$visit, p$assay, censor = FALSE)

      ## ---- methods
      a_T <- BJM:::unlabeled_alpha(src_hat, hide(Tg), long_hat, "year", SV, SF, NULL, 200, NULL,
                                   c(-4, 4), TRUE)
      est <- list(unlabeled_noref = c(a_T$alpha, a_T$alpha_se),
                  unlabeled = c(a_T$alpha - a_S$alpha, sqrt(a_T$alpha_se^2 + a_S$alpha_se^2)))
      fits <- list(source_only = src_hat,
                   unlabeled_noref = with_alpha(src_hat, est$unlabeled_noref[1]),
                   unlabeled = with_alpha(src_hat, est$unlabeled[1]))
      for (nL in unique(pmin(c(50, 200), p$n_T))) {
        lab <- T1[sample(nrow(T1), nL), ]
        ci <- citl_alpha(src_hat, lab)
        est[[paste0("citl", nL)]] <- c(ci$alpha, ci$se)
        fits[[paste0("citl", nL)]] <- with_alpha(src_hat, ci$alpha)
        fits[[paste0("baseline", nL)]] <- transportBJM(src_hat, "baseline", data_target = lab)
        fits[[paste0("recal", nL)]] <- transportBJM(src_hat, "recalibrate", data_target = lab)
      }
      fits$local_refit <- survivalSub(T1, Surv(years, status3) ~ age + sex, cr_form)

      ## ---- evaluation on an independent, uncensored target cohort
      E1 <- one_row(E)
      for (L in LANDMARKS) {
        ids <- E1$id[E1$true_T > L]
        dp <- E[E$id %in% ids & E$year <= L, ]; dp$years <- NA
        y1 <- E1[match(unique(dp$id), E1$id), ]
        y <- y1$true_event == 1 & y1$true_T <= L + H
        pr <- function(surv, long, d) {
          r <- predictRisk(list(d, d), long, surv, prediction_time = L, horizon = H,
                           time_variable = "year", survival_variable_all = SV,
                           survival_trans_function = SF, bandcount1 = 10, bandcount2 = 30)
          total_risk(r)[as.character(y1$id)]
        }
        dp_true <- dp; dp_true$serBilir <- dp_true$serBilir - p$assay   # truth knows the assay offset
        r_or <- pr(truth_surv, truth_long, dp_true)
        for (m in names(fits)) {
          r <- pr(fits[[m]], long_hat, dp)
          e <- est[[m]]
          out[[length(out) + 1]] <- data.frame(
            rep = rep, scenario = p$name, k = p$k, landmark = L, method = m, n_eval = length(y),
            observed = mean(y), predicted = mean(r), oracle = mean(r_or),
            rmse_oracle = sqrt(mean((r - r_or)^2)), brier = mean((r - y)^2),
            brier_oracle = mean((r_or - y)^2), auc = auc(r, y),
            alpha_hat = if (is.null(e)) NA else e[1], alpha_se = if (is.null(e)) NA else e[2],
            slope = if (!is.null(fits[[m]]$transport$lp_slope)) fits[[m]]$transport$lp_slope else NA)
        }
      }
    }
  }
  do.call(rbind, out)
}

for (rep in reps) {
  t0 <- Sys.time()
  res <- tryCatch(run_rep(rep), error = function(e) { message("rep ", rep, " failed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) saveRDS(res, sprintf("sim/full/rep_%03d.rds", rep))
  cat(sprintf("rep %d done in %.0fs\n", rep, as.numeric(Sys.time() - t0, units = "secs")))
}
