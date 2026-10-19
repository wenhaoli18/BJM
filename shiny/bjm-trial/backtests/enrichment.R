# Prognostic enrichment on pbc3, held-out evaluation.
#
# Each replicate splits pbc3 at random into a "historical" half (used to fit
# the models) and a "candidate" half (the patients screened for a trial).
# Candidates are screened for S years (their visits up to S) and, if still
# event-free and in follow-up at S, ranked by their predicted risk of an
# event in (S, S + H]. Enrolling the top share q, the event probability of
# the enrolled patients over the trial's follow-up is estimated by
# Kaplan-Meier from what actually happened to them, and turned into the
# number of patients a trial needs.
#
# Risk models compared:
#   baseline  Cox model on age, sex and the biomarkers at entry (year 0)
#   landmark  Cox model from S on, on age, sex and the last biomarker
#             values up to S (a landmark model)
#   bjm       BJM, from the whole biomarker history up to S
#
# Usage: Rscript enrichment.R <n_reps> <out.rds> [screening years S, default 1]
# The app's table: Rscript enrichment.R 50 en_S1.rds 1; Rscript enrichment.R 50 en_S2.rds 2

suppressMessages({library(BJM); library(survival)})
args <- commandArgs(TRUE)
n_reps <- as.integer(args[1]); out_file <- args[2]

data(pbc3)
bio_vars <- c("serBilir", "albumin", "prothrombin")
pbc3$id <- as.character(pbc3$id)
base <- pbc3[!duplicated(pbc3$id), c("id", "years", "status3", "age", "sex")]
# The longitudinal sub-models' event-time term is capped at T_cap years. A
# term linear in T is extrapolated past the last event time (13.9 years in
# pbc3), which overstates the risk of low- and mid-risk patients
# (backtests/calibration.R); the cap keeps BJM calibrated, and its ranking.
T_cap <- 8
trans <- list(survival_variable_all = list("Tc"),
              survival_trans_function = list(function(x) pmin(x, T_cap)))
long_fixed <- setNames(lapply(bio_vars, function(b)
  as.formula(paste(b, "~ year + age + sex + Tc + Tc * year"))), bio_vars)
long_random <- setNames(rep(list(~ year | id), length(bio_vars)), bio_vars)

S <- if (length(args) >= 3) as.numeric(args[3]) else 1; H <- 3
shares <- c(1, 0.75, 0.5, 1/3, 0.25)

# per patient: biomarkers at entry, and the last value up to S
first_vals <- function(ids) {
  v <- pbc3[pbc3$id %in% ids & pbc3$year == 0, c("id", bio_vars)]
  names(v)[-1] <- paste0(bio_vars, "_0"); v
}
last_vals <- function(ids) {
  v <- pbc3[pbc3$id %in% ids & pbc3$year <= S, c("id", "year", bio_vars)]
  do.call(rbind, lapply(split(v, v$id), function(d) {
    d <- d[order(d$year), ]
    data.frame(id = d$id[1], t(sapply(bio_vars, function(b) {
      x <- d[[b]][!is.na(d[[b]])]; if (length(x)) x[length(x)] else NA })))
  }))
}

# Cox risk of an event in (from, from + H] given event-free at `from`
cox_risk <- function(fit, newdata, from, H) {
  sf <- survfit(fit, newdata = newdata)
  s_at <- function(t) summary(sf, times = t, extend = TRUE)$surv
  s0 <- if (from > 0) s_at(from) else rep(1, nrow(newdata))
  1 - as.numeric(s_at(from + H)) / as.numeric(s0)
}

# event probability over (S, S + H] among the enrolled, by Kaplan-Meier
km_prob <- function(d) {
  t <- pmin(d$years - S, H); e <- d$status3 == 1 & d$years - S <= H
  km <- survfit(Surv(t, e) ~ 1)
  1 - summary(km, times = H, extend = TRUE)$surv
}

# patients needed: Schoenfeld events for a hazard ratio, over the pooled
# event probability of the two arms (1:1)
patients_needed <- function(p, hr = 0.7, alpha = 0.05, power = 0.8) {
  d <- 4 * (qnorm(1 - alpha / 2) + qnorm(power))^2 / log(hr)^2
  d / ((p + 1 - (1 - p)^hr) / 2)
}

results <- list()
for (r in seq_len(n_reps)) {
  set.seed(2000 + r)
  hist_ids <- sample(base$id, nrow(base) / 2)
  hist <- base[base$id %in% hist_ids, ]
  cand <- base[!base$id %in% hist_ids & base$years > S, ]   # event-free and followed at S

  # baseline Cox
  hb <- merge(hist, first_vals(hist$id))
  fit_b <- coxph(Surv(years, status3) ~ age + sex + serBilir_0 + albumin_0 + prothrombin_0, data = hb)
  cb <- merge(cand, first_vals(cand$id))
  # landmark Cox, from S on
  hl <- merge(hist[hist$years > S, ], last_vals(hist$id[hist$years > S]))
  hl$t_lm <- hl$years - S
  fit_l <- coxph(Surv(t_lm, status3) ~ age + sex + serBilir + albumin + prothrombin, data = hl)
  cl <- merge(cand, last_vals(cand$id))
  # BJM
  hv <- pbc3[pbc3$id %in% hist_ids & pbc3$status3 == 1, ]
  hv$Tc <- pmin(hv$years, T_cap)
  fit_ss <- survivalSub(hist, Surv(years, status3) ~ age + sex, NULL)
  fit_ls <- suppressWarnings(longitudinalSub(hv, long_fixed, long_random))

  risk <- data.frame(id = cand$id)
  risk$baseline <- cox_risk(fit_b, cb[match(risk$id, cb$id), ], S, H)
  risk$landmark <- cox_risk(fit_l, cl[match(risk$id, cl$id), ], 0, H)
  risk$bjm <- vapply(risk$id, function(i) {
    h <- pbc3[pbc3$id == i & pbc3$year <= S, c("id", "year", "age", "sex", bio_vars)]
    h$years <- NA
    h$Tc <- NA
    unname(predictRisk(h, fit_ls, fit_ss, prediction_time = S, horizon = H,
                       time_variable = "year", trans$survival_variable_all,
                       trans$survival_trans_function, bandcount1 = 10,
                       bandcount2 = 10)$risk_prob_1)
  }, numeric(1))
  risk <- merge(risk, cand[, c("id", "years", "status3")])

  for (m in c("baseline", "landmark", "bjm")) {
    ok <- !is.na(risk[[m]])
    c_index <- concordance(Surv(pmin(years - S, H), status3 == 1 & years - S <= H) ~ I(-risk[[m]]),
                           data = risk[ok, ])$concordance
    for (q in shares) {
      cutoff <- quantile(risk[[m]][ok], 1 - q, names = FALSE)
      enrolled <- risk[ok & risk[[m]] >= cutoff, ]
      p <- km_prob(enrolled)
      results[[length(results) + 1]] <- data.frame(
        rep = r, method = m, share = q, n_cand = sum(ok), c_index = c_index,
        p_event = p, n_needed = patients_needed(p), n_screened = patients_needed(p) / q)
    }
  }
  message(sprintf("rep %d done", r))
}
results <- do.call(rbind, results)
saveRDS(list(results = results, S = S, H = H), out_file)

# summary
agg <- aggregate(cbind(c_index, p_event, n_needed, n_screened) ~ method + share,
                 data = results, FUN = mean)
agg <- agg[order(-agg$share, match(agg$method, c("baseline", "landmark", "bjm"))), ]
agg[, -(1:2)] <- lapply(agg[, -(1:2)], function(x) round(x, 3))
print(agg, row.names = FALSE)
