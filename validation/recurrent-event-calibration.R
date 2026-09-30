## Simulation calibration check for buildRecurrentEpisodes()'s Tier-1
## "next event" risk predictions.
##
## This is NOT a testthat test and is not part of the installed package
## (it lives outside R/, so R CMD check never sees it). It answers one
## narrow, specific question: under a data-generating process for which
## the gap-time Cox + LMM sub-models are *correctly specified* (right
## functional form for the hazard, the recurrence-vs-terminal split, and
## the biomarker mean), does predictRisk() -- run on data reshaped by
## buildRecurrentEpisodes() -- give predicted probabilities that match
## reality? If calibration fails here, that is a real bug (in
## buildRecurrentEpisodes()'s reshaping, or in how the existing engine
## responds to gap-time input). It does NOT test robustness to model
## misspecification, which is a separate, harder question about how well
## this extension would hold up on real data.
##
## Design: simulate recurrent-event + longitudinal histories from a known
## DGP (Weibull-PH gap times; logistic recurrence-vs-terminal split;
## linear mixed-effects biomarker whose within-episode trajectory shape
## depends on the eventual gap length -- exactly the functional form
## long_sub_fixed asks for). Fit survivalSub()/longitudinalSub() on a
## training cohort. For an independent test cohort, landmark each subject
## at a random point *inside* one of their own (fully resolved) episodes,
## keep only their history up to that landmark, and predict the risk that
## their *next* event (within a fixed horizon) is a recurrence / a
## terminal event. Compare predicted probabilities to the realized
## (ground-truth, since we generated the data) outcome via binned
## calibration plots, a calibration slope/intercept, and the Brier score.

devtools::load_all(".", quiet = TRUE)

set.seed(20260930)

## ---- 1. True data-generating process ----------------------------------

params <- list(
  kshape = 1.3, lambda0 = 0.7, beta_age = 0.10, beta_sex = 0.30,   # gap hazard
  alpha0 = -1.5, alpha1 = 0.40, alpha_age = 0.30, alpha_sex = 0.20, # P(terminal | gap)
  b0 = 5, b1 = -0.30, b2 = 0.50, b3 = -0.80, b_age = 0.05, b_sex = 0.20, # biomarker mean
  sigma_u0 = 0.6, sigma_u1 = 0.3, sigma_eps = 0.4,                  # biomarker variance
  max_episodes = 8
)

simulate_subject <- function(id, params) {
  age <- rnorm(1, 50, 8)
  sex <- rbinom(1, 1, 0.5)
  age_c <- (age - 50) / 10
  admin_cutoff <- runif(1, 3, 6)

  u0 <- rnorm(1, 0, params$sigma_u0)
  u1 <- rnorm(1, 0, params$sigma_u1)

  t_cal <- 0
  events <- vector("list", params$max_episodes)
  long_rows <- vector("list", params$max_episodes)
  k <- 0L

  repeat {
    k <- k + 1L
    lin_gap <- params$beta_age * age_c + params$beta_sex * sex
    U <- runif(1)
    g_true <- (-log(U) / (params$lambda0 * exp(lin_gap)))^(1 / params$kshape)

    censored_now <- (t_cal + g_true) > admin_cutoff
    if (censored_now) {
      g_obs <- admin_cutoff - t_cal
      status <- 0L
    } else {
      g_obs <- g_true
      lin_type <- params$alpha0 + params$alpha1 * g_obs + params$alpha_age * age_c + params$alpha_sex * sex
      is_terminal <- rbinom(1, 1, plogis(lin_type))
      status <- if (is_terminal) 2L else 1L
    }
    t_end <- t_cal + g_obs

    events[[k]] <- data.frame(id = id, time = t_end, status = status, age = age, sex = sex)

    ## Measurements start a tiny bit *after* the episode actually begins --
    ## never exactly at the prior event's timestamp -- matching how real
    ## data is collected (a biomarker drawn at the same visit an event is
    ## recorded is still drawn a nonzero instant after it, down to whatever
    ## precision the timestamps carry). Landing exactly on a boundary would
    ## make every non-final episode's transition point coincide exactly
    ## with calendar time, which buildRecurrentEpisodes() then (correctly,
    ## for genuinely ambiguous real data, see ?buildRecurrentEpisodes)
    ## treats as belonging to *both* adjacent episodes -- duplicating that
    ## single simulated biomarker draw under two different, only one of
    ## which is how it was actually generated here. That is a sensible
    ## real-data default, but it is not what this calibration check is
    ## trying to measure, so the DGP simply avoids ever creating an exact
    ## coincidence in the first place.
    tie_seq <- if (g_obs > 0.05) seq(0.05, g_obs, by = 0.2) else numeric(0)
    if (length(tie_seq) > 0) {
      mean_y <- params$b0 + u0 + (params$b1 + u1) * tie_seq + params$b2 * g_obs +
        params$b3 * g_obs * tie_seq + params$b_age * age_c + params$b_sex * sex
      y <- rnorm(length(tie_seq), mean_y, params$sigma_eps)
      long_rows[[k]] <- data.frame(
        id = id, year = t_cal + tie_seq, biomarker = y,
        .episode = k, .gap_true = g_obs, .tie = tie_seq, .status = status
      )
    } else {
      long_rows[[k]] <- data.frame(
        id = integer(0), year = numeric(0), biomarker = numeric(0),
        .episode = integer(0), .gap_true = numeric(0), .tie = numeric(0), .status = integer(0)
      )
    }

    t_cal <- t_end
    if (status != 1L || k >= params$max_episodes) break
  }

  list(
    events = do.call(rbind, events[seq_len(k)]),
    long = do.call(rbind, long_rows[seq_len(k)]),
    age = age, sex = sex, n_episodes = k
  )
}

## ---- 2. Training cohort: fit the sub-models ----------------------------

n_train <- 800
train_sims <- lapply(seq_len(n_train), simulate_subject, params = params)

train_events <- do.call(rbind, lapply(train_sims, `[[`, "events"))
train_long_full <- do.call(rbind, lapply(train_sims, `[[`, "long"))
train_long <- train_long_full[c("id", "year", "biomarker")]

built_train <- buildRecurrentEpisodes(train_events, train_long, time_variable_long = "year")
built_train$long_data$id <- factor(built_train$long_data$id)

cat(sprintf(
  "Training cohort: %d subjects, %d episodes (%d recurrences, %d terminal, %d admin-censored)\n",
  n_train, nrow(built_train$episode_data),
  sum(built_train$episode_data$event_type == 0, na.rm = TRUE),
  sum(built_train$episode_data$event_type == 1, na.rm = TRUE),
  sum(built_train$episode_data$event_any == 0)
))

survival_fit <- survivalSub(
  built_train$episode_data,
  Surv(gap, event_any) ~ age + sex,
  event_type ~ gap + age + sex
)
long_fit <- longitudinalSub(
  built_train$long_data[built_train$long_data$event_any == 1, ],
  biomarker ~ time_in_episode + age + sex + gap + gap * time_in_episode,
  ~ time_in_episode | id
)

cat("\nFitted gap-time Cox coefficients (truth: age = 0.10, sex = 0.30):\n")
print(coef(survival_fit$coxph_fit))
cat("\nFitted recurrence-vs-terminal glm coefficients (truth: intercept -1.5, gap 0.4, age 0.3, sex 0.2):\n")
print(coef(survival_fit$glm_fit))
cat("\nFitted biomarker fixed effects (truth: intercept 5, time_in_episode -0.3, gap 0.5, interaction -0.8):\n")
print(nlme::fixef(long_fit$lfit[[1]]))

## ---- 3. Test cohort: landmark mid-episode, predict the next event -----

n_test <- 500
horizon <- 0.8

test_cases <- list()
attempts <- 0L
while (length(test_cases) < n_test && attempts < n_test * 5L) {
  attempts <- attempts + 1L
  sim <- simulate_subject(id = 100000 + attempts, params = params)
  resolved <- sim$long[sim$long$.status %in% c(1, 2), ]
  if (nrow(resolved) == 0) next
  ep_choices <- unique(resolved$.episode)
  # sample(x, 1) on a length-1 x samples from 1:x, not the value x itself --
  # guard against that classic gotcha explicitly.
  k_star <- if (length(ep_choices) == 1) ep_choices else sample(ep_choices, 1)
  ep_rows <- sim$long[sim$long$.episode == k_star, ]
  g_true <- ep_rows$.gap_true[1]
  true_status <- ep_rows$.status[1]
  if (g_true < 0.4) next # too short an episode to get a meaningful mid-episode landmark
  landmark <- runif(1, 0.2 * g_true, 0.8 * g_true)

  ep_start_cal <- min(ep_rows$year) # calendar time this episode started
  ## observed history: every completed episode strictly before k_star,
  ## unchanged, plus episode k_star truncated to the landmark with a
  ## trailing status = 0 ("still at risk as of now") row -- exactly the
  ## documented buildRecurrentEpisodes()/predictRisk() usage pattern for a
  ## subject who has not yet resolved their current episode.
  prior_events <- sim$events[sim$events$id == sim$events$id[1] & seq_len(nrow(sim$events)) < k_star, ]
  cur_censor_row <- data.frame(
    id = sim$events$id[1], time = ep_start_cal + landmark, status = 0L,
    age = sim$events$age[1], sex = sim$events$sex[1]
  )
  obs_events <- rbind(prior_events, cur_censor_row)

  obs_long <- sim$long[sim$long$.episode < k_star | (sim$long$.episode == k_star & sim$long$.tie <= landmark), ]
  obs_long <- obs_long[c("id", "year", "biomarker")]

  test_cases[[length(test_cases) + 1]] <- list(
    events = obs_events, long = obs_long, landmark = landmark,
    g_true = g_true, true_status = true_status,
    remaining = g_true - landmark
  )
}
cat(sprintf("\nBuilt %d test landmark cases (from %d simulated subjects)\n", length(test_cases), attempts))

pred1 <- pred2 <- numeric(length(test_cases))
outcome_recur <- outcome_terminal <- logical(length(test_cases))

## predictRisk()'s own upper_bound = 2 * max(gap column of data_predict)
## (see predictRisk.R). For a still-OPEN episode, the "gap" value
## buildRecurrentEpisodes() records is just elapsed time so far -- which,
## used directly, makes upper_bound ~= 2 * landmark, an *upper bound on the
## denominator's integration range* that has no relationship to the true
## support of the gap distribution and can badly starve the "still at risk"
## denominator integral. It is only ever used for (a) the subset_at_risk()
## eligibility check gap >= prediction_time and (b) sizing upper_bound --
## never as the actual value being predicted -- so it is safe to override
## it to something comfortably larger than the training cohort's observed
## gap range once buildRecurrentEpisodes() has done its reshaping.
gap_ceiling <- 1.2 * max(built_train$episode_data$gap)

for (i in seq_along(test_cases)) {
  tc <- test_cases[[i]]
  built_i <- buildRecurrentEpisodes(tc$events, tc$long, time_variable_long = "year")
  built_i$long_data$id <- factor(built_i$long_data$id)
  cur_ep <- max(built_i$long_data$episode)
  data_predict <- built_i$long_data[built_i$long_data$episode == cur_ep, ]
  data_predict$gap <- gap_ceiling

  risk <- tryCatch(
    predictRisk(
      data_predict, long_fit, survival_fit,
      prediction_time = tc$landmark, horizon = horizon, time_variable = "time_in_episode",
      survival_variable_all = list(), survival_trans_function = list(),
      bandcount1 = 30, bandcount2 = 30
    ),
    error = function(e) NULL
  )
  if (is.null(risk)) {
    pred1[i] <- pred2[i] <- NA
  } else {
    pred1[i] <- risk$risk_prob_1
    pred2[i] <- risk$risk_prob_2
  }
  outcome_recur[i] <- (tc$true_status == 1) && (tc$remaining <= horizon)
  outcome_terminal[i] <- (tc$true_status == 2) && (tc$remaining <= horizon)

  if (i %% 50 == 0) cat(sprintf("  ...%d/%d predictions done\n", i, length(test_cases)))
}

ok <- !is.na(pred1) & !is.na(pred2)
cat(sprintf("\n%d/%d predictions succeeded (%d failed numerically and were dropped)\n",
            sum(ok), length(ok), sum(!ok)))

pred1 <- pred1[ok]; pred2 <- pred2[ok]
outcome_recur <- outcome_recur[ok]; outcome_terminal <- outcome_terminal[ok]

saveRDS(list(pred1 = pred1, pred2 = pred2, outcome_recur = outcome_recur,
             outcome_terminal = outcome_terminal),
        "validation/calibration_raw.rds")
cat(sprintf("\npred1 range: [%.4f, %.4f], %% at/near 0 or 1: %.1f%%\n",
            min(pred1), max(pred1), 100 * mean(pred1 < 0.01 | pred1 > 0.99)))
cat(sprintf("pred2 range: [%.4f, %.4f], %% at/near 0 or 1: %.1f%%\n",
            min(pred2), max(pred2), 100 * mean(pred2 < 0.01 | pred2 > 0.99)))

## ---- 4. Calibration summary ---------------------------------------------

calibration_table <- function(pred, outcome, n_bins = 10) {
  breaks <- unique(quantile(pred, probs = seq(0, 1, length.out = n_bins + 1)))
  if (length(breaks) < 3) {
    breaks <- c(min(pred) - 1e-9, max(pred) + 1e-9)
  }
  bin <- cut(pred, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  n_bins_actual <- length(breaks) - 1
  data.frame(
    bin = seq_len(n_bins_actual),
    n = as.integer(table(factor(bin, levels = seq_len(n_bins_actual)))),
    mean_predicted = tapply(pred, bin, mean),
    observed_rate = tapply(outcome, bin, mean)
  )
}

brier_score <- function(pred, outcome) mean((pred - outcome)^2)

calibration_intercept_slope <- function(pred, outcome) {
  eps <- 1e-6
  p <- pmin(pmax(pred, eps), 1 - eps)
  logit_p <- log(p / (1 - p))
  fit <- glm(outcome ~ logit_p, family = binomial())
  coef(fit)
}

cal1 <- calibration_table(pred1, outcome_recur)
cal2 <- calibration_table(pred2, outcome_terminal)

cat("\n== Calibration: P(next event is a recurrence within horizon) ==\n")
print(cal1, row.names = FALSE)
cat(sprintf("Brier score: %.4f | mean predicted: %.4f | observed rate: %.4f\n",
            brier_score(pred1, outcome_recur), mean(pred1), mean(outcome_recur)))
cs1 <- calibration_intercept_slope(pred1, outcome_recur)
cat(sprintf("Calibration intercept/slope (on logit scale; 0/1 = perfect): %.3f / %.3f\n", cs1[1], cs1[2]))

cat("\n== Calibration: P(next event is the terminal event within horizon) ==\n")
print(cal2, row.names = FALSE)
cat(sprintf("Brier score: %.4f | mean predicted: %.4f | observed rate: %.4f\n",
            brier_score(pred2, outcome_terminal), mean(pred2), mean(outcome_terminal)))
cs2 <- calibration_intercept_slope(pred2, outcome_terminal)
cat(sprintf("Calibration intercept/slope (on logit scale; 0/1 = perfect): %.3f / %.3f\n", cs2[1], cs2[2]))

saveRDS(list(pred1 = pred1, pred2 = pred2, outcome_recur = outcome_recur,
             outcome_terminal = outcome_terminal, cal1 = cal1, cal2 = cal2),
        "validation/calibration_results.rds")

## ---- 5. Plot -------------------------------------------------------------

png("validation/calibration_plot.png", width = 1000, height = 500, res = 130)
par(mfrow = c(1, 2))
plot(cal1$mean_predicted, cal1$observed_rate, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Mean predicted risk_prob_1 (bin)", ylab = "Observed recurrence rate (bin)",
     main = "Calibration: next event = recurrence", pch = 19, col = "steelblue")
abline(0, 1, lty = 2, col = "gray40")
plot(cal2$mean_predicted, cal2$observed_rate, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Mean predicted risk_prob_2 (bin)", ylab = "Observed terminal-event rate (bin)",
     main = "Calibration: next event = terminal", pch = 19, col = "firebrick")
abline(0, 1, lty = 2, col = "gray40")
dev.off()

cat("\nSaved validation/calibration_plot.png and validation/calibration_results.rds\n")
