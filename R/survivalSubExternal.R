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
#' @return A \code{survivalSub.BJM} fit, usable wherever \code{survivalSub()}'s
#'   is (e.g. \code{predictRisk()}, \code{fitImputedBJM(survival_fit_all = )}),
#'   with an \code{external} element: the generation variable, the cohort
#'   sizes and follow-up, the hazard ratio of the study against the external
#'   cohort and the PH test.
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
#' }
#' @export
survivalSubExternal <- function(data_internal, data_external, form_marginal_surv,
                                form_conditional_cr = NULL, generation = "generation") {
  assert_data_frame(data_internal, "data_internal")
  assert_data_frame(data_external, "data_external")
  assert_string(generation, "generation")
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
  fit$external <- list(generation = generation,
                       n_internal = nrow(data_internal), n_external = nrow(data_external),
                       max_followup_internal = tau, max_followup_external = max(y_ext[, "time"]),
                       hazard_ratio = unname(exp(stats::coef(fit$coxph_fit)[generation])),
                       ph_test = ph_test)
  fit
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
