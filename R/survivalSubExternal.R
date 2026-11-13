#' Survival sub-model that borrows a long-term tail from an external cohort
#'
#' @description \strong{Experimental.} Predicting beyond the last follow-up
#' of the study needs the baseline hazard there, which \code{survivalSub()}
#' can only extrapolate (linearly in the cumulative hazard). When an
#' external cohort with longer follow-up is available -- an earlier
#' generation of the same device, a registry, an older trial -- its
#' baseline hazard can supply the shape of the tail.
#' \code{survivalSubExternal()} fits one Cox model to both cohorts with a
#' generation indicator (1 for the study, 0 for the external cohort):
#' \deqn{h(t \mid x, g) = h_0(t) \exp(\beta^\top x + \gamma g),}
#' so the study's hazard beyond its own follow-up is the external one
#' scaled by \eqn{e^\gamma}. Only the survival sub-model uses the external
#' cohort; the longitudinal sub-model is fit on the study as usual, and the
#' external cohort needs no biomarkers.
#'
#' The key assumption is proportional hazards between the cohorts, also
#' beyond the study's follow-up, where it cannot be checked. Within the
#' overlap it is tested with \code{\link[survival]{cox.zph}} (stored in
#' \code{$external$ph_test}, with a warning when \eqn{p < 0.05}); the test has
#' little power with short study follow-up, so a non-significant result is
#' weak reassurance, and a sensitivity analysis on the tail is advisable.
#' The hazard ratio itself, estimated from the overlap only, can be noisy
#' when the study is short and has few events, and the proportional hazards
#' assumption carries that noise to the whole tail. \code{tail_hr} fixes it
#' instead (the other coefficients are then re-estimated given it), for a
#' sensitivity analysis over plausible values, e.g. \code{c(0.8, 1, 1.25)},
#' or to use prior knowledge, e.g. that the generations are equivalent.
#'
#' With \code{splice = TRUE} (the default) the external cohort only supplies
#' the tail: up to the study's last follow-up \eqn{\tau} the baseline
#' cumulative hazard is the Breslow estimate from the study alone (at the
#' coefficients of the joint fit), and beyond \eqn{\tau} it continues with the
#' increments of the joint fit's baseline, i.e. the external shape scaled by
#' \eqn{e^\gamma}. Within the study's follow-up predictions then rest on the
#' study's own events. With \code{splice = FALSE} the joint fit's baseline is
#' used throughout, so the external cohort (scaled by \eqn{e^\gamma}) also
#' shapes the study's hazard within its follow-up, which matters when the
#' cohorts are not proportional there.
#'
#' Predictions are for the study population: when the prediction data have
#' no \code{generation} column it is taken as 1.
#'
#' @param data_internal,data_external One row per subject, for the study and
#'   for the external cohort. Both need the variables of
#'   \code{form_marginal_surv} (and of \code{form_conditional_cr}, if given).
#' @param form_marginal_surv,form_conditional_cr As in
#'   \code{\link{survivalSub}}, without the generation term, which is added.
#' @param generation Name of the generation indicator column to create (it
#'   must not already exist).
#' @param tail_hr \code{NULL} (default) to estimate the hazard ratio of the
#'   study against the external cohort, or a positive number to fix it.
#' @param splice If \code{TRUE} (default), use the study's own baseline up to
#'   its last follow-up and the external shape only beyond it (see Details).
#'   Not available for a stratified Cox model.
#' @return A \code{survivalSub.BJM} fit, usable wherever \code{survivalSub()}'s
#'   is (e.g. \code{predictRisk()}, \code{fitImputedBJM(survival_fit_all = )}),
#'   with an \code{external} element: the generation variable, the cohort
#'   sizes and follow-up, the hazard ratio of the study against the external
#'   cohort (\code{hazard_ratio}, the one used; \code{hazard_ratio_estimated},
#'   from the data; \code{tail_hr_fixed}) and the PH test.
#' @examples
#' \donttest{
#' data(pbc3)
#' surv <- pbc3[!duplicated(pbc3$id), ]
#' # pretend the study had 4 years of follow-up and the rest of pbc3 is external
#' id <- as.numeric(as.character(surv$id))
#' study <- surv[id %% 2 == 0, ]
#' study$status3[study$years > 4] <- 0
#' study$years <- pmin(study$years, 4)
#' external <- surv[id %% 2 == 1, ]
#' fit <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, NULL)
#' fit$external
#' # sensitivity to the hazard ratio between the cohorts
#' fits <- lapply(c(0.8, 1, 1.25), function(hr)
#'   survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, NULL, tail_hr = hr))
#' }
#' @export
survivalSubExternal <- function(data_internal, data_external, form_marginal_surv,
                                form_conditional_cr = NULL, generation = "generation",
                                tail_hr = NULL, splice = TRUE) {
  assert_data_frame(data_internal, "data_internal")
  assert_data_frame(data_external, "data_external")
  assert_string(generation, "generation")
  if (!is.null(tail_hr) && !(is.numeric(tail_hr) && length(tail_hr) == 1 && is.finite(tail_hr) && tail_hr > 0)) {
    stop("`tail_hr` must be NULL or a single positive number.", call. = FALSE)
  }
  if (!is.logical(splice) || length(splice) != 1 || is.na(splice)) {
    stop("`splice` must be TRUE or FALSE.", call. = FALSE)
  }
  if (generation %in% c(names(data_internal), names(data_external))) {
    stop(sprintf("`generation` (\"%s\") is already a column; choose another name.", generation),
         call. = FALSE)
  }
  vars <- unique(c(all.vars(form_marginal_surv), all.vars(form_conditional_cr)))
  assert_vars_in_data(vars, data_internal, "the survival formulas", "data_internal")
  assert_vars_in_data(vars, data_external, "the survival formulas", "data_external")
  if (is_interval_surv(form_marginal_surv, data_internal)) {
    stop("survivalSubExternal() supports right-censored outcomes only.", call. = FALSE)
  }
  a <- data_internal[vars]
  b <- data_external[vars]
  a[[generation]] <- 1
  b[[generation]] <- 0
  stacked <- rbind(a, b)
  form <- stats::update(form_marginal_surv, stats::as.formula(paste(". ~ . +", generation)))
  environment(form) <- environment(form_marginal_surv)
  fit <- survivalSub(stacked, form, form_conditional_cr)
  hr_estimated <- unname(exp(stats::coef(fit$coxph_fit)[generation]))
  if (!is.null(tail_hr)) fit$coxph_fit <- cox_fixed_generation(fit$coxph_fit, stacked, form, generation, tail_hr)

  y_int <- eval(form_marginal_surv[[2]], data_internal, environment(form_marginal_surv))
  y_ext <- eval(form_marginal_surv[[2]], data_external, environment(form_marginal_surv))
  tau <- max(y_int[, "time"])
  # PH between the cohorts on the overlap, where both are observed
  overlap <- stacked
  y_all <- eval(form_marginal_surv[[2]], stacked, environment(form_marginal_surv))
  overlap$.bjm_time <- pmin(y_all[, "time"], tau)
  overlap$.bjm_status <- y_all[, "status"] * (y_all[, "time"] <= tau)
  rhs <- paste(c(attr(stats::terms(form_marginal_surv), "term.labels"), generation), collapse = " + ")
  ph_test <- tryCatch({
    z <- survival::cox.zph(survival::coxph(stats::as.formula(paste("Surv(.bjm_time, .bjm_status) ~", rhs)),
                                           data = overlap))
    z$table[generation, ]
  }, error = function(e) NULL)
  if (!is.null(ph_test) && ph_test[["p"]] < 0.05) {
    warning(sprintf(paste0(
      "The hazards of the study and the external cohort do not look proportional within the study's ",
      "follow-up (cox.zph p = %.3g), so the borrowed tail may be off."), ph_test[["p"]]), call. = FALSE)
  }
  spliced <- NULL
  if (splice) {
    if (length(survival::untangle.specials(stats::terms(fit$coxph_fit), "strata")$vars) > 0) {
      stop("`splice = TRUE` is not available for a stratified Cox model.", call. = FALSE)
    }
    fit$external <- list(generation = generation)
    spliced <- splice_baseline(fit, data_internal, y_int, tau)
  }
  fit$external <- list(generation = generation,
                       n_internal = nrow(data_internal), n_external = nrow(data_external),
                       max_followup_internal = tau, max_followup_external = max(y_ext[, "time"]),
                       hazard_ratio = unname(exp(stats::coef(fit$coxph_fit)[generation])),
                       hazard_ratio_estimated = hr_estimated, tail_hr_fixed = !is.null(tail_hr),
                       ph_test = ph_test, splice = splice, cum_basehaz = spliced)
  fit
}

#' Baseline cumulative hazard spliced at the study's last follow-up
#'
#' @description For \code{survivalSubExternal(splice = TRUE)}: the Breslow
#' estimate from the study alone, at the joint fit's coefficients (with the
#' study's generation term), up to \code{tau}; beyond it, the joint fit's
#' baseline increments. Both are on the scale of the joint fit's uncentered
#' linear predictor, so \code{survival_lp()} applies unchanged.
#' @keywords internal
splice_baseline <- function(fit, data_internal, y_int, tau) {
  lp <- survival_lp(fit, data_internal)
  time <- y_int[, "time"]
  status <- y_int[, "status"]
  ok <- !is.na(lp) & !is.na(time) & !is.na(status)
  time <- time[ok]; status <- status[ok]; w <- exp(lp[ok])
  times <- sort(unique(time))
  d <- vapply(times, function(t) sum(status[time == t]), numeric(1))
  at_risk <- vapply(times, function(t) sum(w[time >= t]), numeric(1))
  internal <- data.frame(hazard = cumsum(d / at_risk), time = times)
  joint <- survival::basehaz(fit$coxph_fit, centered = FALSE)[, c("hazard", "time")]
  h_tau <- stats::approx(c(0, joint$time), c(0, joint$hazard), xout = tau, method = "constant", rule = 2)$y
  beyond <- joint[joint$time > tau, , drop = FALSE]
  beyond$hazard <- internal$hazard[nrow(internal)] + beyond$hazard - h_tau
  rbind(internal, beyond)
}

#' Cox model with the generation coefficient fixed
#'
#' @description For \code{survivalSubExternal(tail_hr = )}: the other
#' coefficients are estimated with the generation term as an offset, then
#' the model is refit at all coefficients with no iterations, so that the
#' result is an ordinary \code{coxph} fit (its prediction and Breslow
#' baseline then need no special handling of the offset).
#' @keywords internal
cox_fixed_generation <- function(free_fit, stacked, form, generation, tail_hr) {
  rhs <- attr(stats::terms(form), "term.labels")
  others <- setdiff(rhs, generation)
  off <- sprintf("offset(%.17g * %s)", log(tail_hr), generation)
  form_off <- stats::reformulate(c(others, off), response = form[[2]], env = environment(form))
  beta <- stats::coef(free_fit)
  beta[] <- 0
  beta[generation] <- log(tail_hr)
  if (length(others) > 0) {
    off_fit <- coxph(form_off, data = stacked, ties = "breslow")
    beta[names(stats::coef(off_fit))] <- stats::coef(off_fit)
  }
  suppressWarnings(coxph(form, data = stacked, ties = "breslow", x = TRUE, y = TRUE, init = beta,
                         control = survival::coxph.control(iter.max = 0)))
}

#' Prediction data for an external-tail fit
#'
#' @description Adds the generation indicator (1, the study) to
#' \code{newdata} when a \code{survivalSubExternal()} fit needs it and it is
#' missing.
#' @keywords internal
with_generation <- function(survival_fit_all, newdata) {
  g <- survival_fit_all$external$generation
  if (!is.null(g) && !g %in% names(newdata)) newdata[[g]] <- rep(1, nrow(newdata))
  newdata
}
