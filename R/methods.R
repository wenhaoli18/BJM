# -- Internal formatting helper (not exported) ----------------------------------
.format_survivalSub <- function(x, digits = 4, extended = FALSE) {
  if (is_interval_censored(x)) return(.format_survivalSub_ic(x, digits, extended))

  sep_line  <- paste(rep("=", 65), collapse = "")
  dash_line <- paste(rep("-", 65), collapse = "")

  coxph_fit   <- x$coxph_fit
  cox_formula <- x$form_marginal_surv
  glm_fit     <- x$glm_fit
  glm_form    <- x$form_conditional_cr
  has_cr      <- !is.null(glm_fit)
  s_cox       <- summary(coxph_fit)

  # -- Header ------------------------------------------------------------
  cat("\n")
  cat("Call:\n")
  cat(sprintf("survivalSub(form_marginal_surv = %s",
              deparse(cox_formula, width.cutoff = 50L)))
  if (has_cr)
    cat(sprintf(",\n            form_conditional_cr = %s",
                deparse(glm_form, width.cutoff = 50L)))
  cat(")\n")

  # -- Data descriptives -------------------------------------------------
  cat("\n")
  cat("Data Descriptives:\n")
  cat(sprintf("  Number of subjects        : %d\n", coxph_fit$n))
  cat(sprintf("  Number of events          : %d\n", coxph_fit$nevent))
  if (has_cr) {
    n_ev   <- nrow(glm_fit$data)
    n_cr   <- sum(glm_fit$y)
    cat(sprintf("  Cause-1 events (CR model) : %d\n",
                as.integer(n_cr)))
    cat(sprintf("  Cause-2 events (CR model) : %d\n",
                as.integer(n_ev - n_cr)))
  }

  # -- Marginal Survival Sub-model ---------------------------------------
  cat("\n", sep_line, "\n", sep = "")
  cat(" Marginal Survival Sub-model  [Cox PH]\n")
  cat(dash_line, "\n", sep = "")
  cat(" Formula: ")
  print(cox_formula)
  cat("\n")

  coef_cox <- s_cox$coefficients
  tab_cox  <- coef_cox[, c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"),
                        drop = FALSE]
  colnames(tab_cox) <- c("Coef", "exp(Coef)", "SE", "z", "p-value")
  stats::printCoefmat(tab_cox,
               digits       = digits,
               P.values     = TRUE,
               has.Pvalue   = TRUE,
               signif.stars = getOption("show.signif.stars"),
               cs.ind       = 1:3,
               tst.ind      = 4)

  cat("\n")
  cat(sprintf("  n = %d,  events = %d\n", coxph_fit$n, coxph_fit$nevent))
  cat(sprintf("  Concordance       = %.3f  (se = %.4f)\n",
              s_cox$concordance["C"], s_cox$concordance["se(C)"]))
  cat(sprintf("  Likelihood ratio  = %.2f  on %d df,  p = %.4g\n",
              s_cox$logtest["test"], s_cox$logtest["df"],
              s_cox$logtest["pvalue"]))
  cat(sprintf("  Wald test         = %.2f  on %d df,  p = %.4g\n",
              s_cox$waldtest["test"], s_cox$waldtest["df"],
              s_cox$waldtest["pvalue"]))
  cat(sprintf("  Score (logrank)   = %.2f  on %d df,  p = %.4g\n",
              s_cox$sctest["test"], s_cox$sctest["df"],
              s_cox$sctest["pvalue"]))

  # -- Extended: baseline hazard summary ---------------------------------
  if (extended) {
    bh <- survival::basehaz(coxph_fit, centered = FALSE)
    cat(sprintf("\n  Baseline cumulative hazard range: [%.4f, %.4f]\n",
                min(bh$hazard), max(bh$hazard)))
    cat(sprintf("  Time range: [%.4f, %.4f]\n",
                min(bh$time), max(bh$time)))
  }

  # -- Conditional Competing-Risks Sub-model -----------------------------
  if (has_cr) {
    cat("\n", sep_line, "\n", sep = "")
    cat(" Conditional Competing-Risks Sub-model  [Logistic GLM]\n")
    cat(dash_line, "\n", sep = "")
    cat(" Formula: ")
    print(glm_form)
    cat("\n")

    s_glm   <- summary(glm_fit)
    tab_glm <- s_glm$coefficients
    colnames(tab_glm) <- c("Coef", "SE", "z", "p-value")
    stats::printCoefmat(tab_glm,
                 digits       = digits,
                 P.values     = TRUE,
                 has.Pvalue   = TRUE,
                 signif.stars = getOption("show.signif.stars"),
                 cs.ind       = 1:2,
                 tst.ind      = 3)

    cat("\n")
    cat(sprintf("  Null deviance     : %.2f  on %d df\n",
                glm_fit$null.deviance, glm_fit$df.null))
    cat(sprintf("  Residual deviance : %.2f  on %d df\n",
                stats::deviance(glm_fit), glm_fit$df.residual))
    cat(sprintf("  AIC : %.2f\n", stats::AIC(glm_fit)))

    if (extended) {
      cat(sprintf("  BIC : %.2f\n", stats::BIC(glm_fit)))
      cat(sprintf("  McFadden R2 : %.4f\n",
                  1 - glm_fit$deviance / glm_fit$null.deviance))
    }
  }

  cat(sep_line, "\n\n", sep = "")
}


# -- Interval-censored survival sub-model (not exported) -------------------------
icph_coef_table <- function(ic_fit) {
  b <- ic_fit$coefficients
  se <- sqrt(diag(ic_fit$var))
  z <- b / se
  cbind(Coef = b, `exp(Coef)` = exp(b), SE = se, z = z,
        `p-value` = 2 * stats::pnorm(-abs(z)))
}

.format_survivalSub_ic <- function(x, digits = 4, extended = FALSE) {
  sep_line  <- paste(rep("=", 65), collapse = "")
  dash_line <- paste(rep("-", 65), collapse = "")
  ic <- x$ic_fit

  cat("\nCall:\n")
  cat(sprintf("survivalSub(form_marginal_surv = %s,\n            event_time = \"%s\")\n",
              deparse(x$form_marginal_surv, width.cutoff = 50L), ic$event_time))
  cat("\nData Descriptives:\n")
  cat(sprintf("  Number of subjects        : %d\n", ic$n))
  cat(sprintf("  Interval-censored         : %d\n", ic$n_interval))
  cat(sprintf("  Exactly observed          : %d\n", ic$n_exact))
  cat(sprintf("  Right-censored            : %d\n", ic$n_right))

  cat("\n", sep_line, "\n", sep = "")
  cat(sprintf(" Marginal Survival Sub-model  [PH, %s baseline, interval-censored]\n",
              if (ic$baseline == "spline") sprintf("spline (df = %d)", ic$df) else "piecewise-constant"))
  cat(dash_line, "\n", sep = "")
  cat(" Formula: ")
  print(x$form_marginal_surv)
  cat("\n")
  if (length(ic$coefficients) > 0) {
    stats::printCoefmat(icph_coef_table(ic), digits = digits, P.values = TRUE, has.Pvalue = TRUE,
                        signif.stars = getOption("show.signif.stars"), cs.ind = 1:3, tst.ind = 4)
  }
  cat(sprintf("\n  Log-likelihood    = %.2f\n", ic$loglik))
  if (extended && ic$baseline == "spline") {
    cat("  Spline knots (time scale):", paste(signif(ic$base$knot_times, digits), collapse = ", "), "\n")
  }
  if (extended && ic$baseline == "piecewise") {
    finite_cuts <- ic$cuts[is.finite(ic$cuts)]
    cat("  Baseline hazard by piece:\n")
    print(data.frame(from = finite_cuts, to = c(finite_cuts[-1], Inf), hazard = signif(ic$lambda, digits)),
          row.names = FALSE)
  }
  cat(sep_line, "\n\n", sep = "")
}


#' Print method for \code{survivalSub.BJM} objects
#'
#' Automatically called when you type \code{survival_fit_all} or
#' \code{print(survival_fit_all)} at the console. Displays a JMbayes2-style
#' formatted summary of the survival sub-model.
#'
#' @param x A \code{survivalSub.BJM} object returned by \code{\link{survivalSub}}.
#' @param digits Number of significant digits. Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' data(pbc3)
#' data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
#' form_marginal_surv  <- Surv(years, status3) ~ age + sex
#' form_conditional_cr <- status4 ~ years + age + sex
#' survival_fit_all  <- survivalSub(data_survival_fitting,
#'                                  form_marginal_surv, form_conditional_cr)
#' survival_fit_all   # triggers print.survivalSub.BJM automatically
#'
#' @keywords internal
#' @export
print.survivalSub.BJM <- function(x, digits = 4, ...) {
  .format_survivalSub(x, digits = digits, extended = FALSE)
  invisible(x)
}


#' Summary method for \code{survivalSub.BJM} objects
#'
#' Called via \code{summary(survival_fit_all)}. Returns (and prints) an
#' extended summary including baseline hazard range, BIC, and McFadden R2
#' for the competing-risks GLM.
#'
#' @param object A \code{survivalSub.BJM} object returned by \code{\link{survivalSub}}.
#' @param digits Number of significant digits. Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list with components \code{cox_summary} and
#'   (if competing risks) \code{glm_summary}.
#'
#' @examples
#' data(pbc3)
#' data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
#' form_marginal_surv  <- Surv(years, status3) ~ age + sex
#' form_conditional_cr <- status4 ~ years + age + sex
#' survival_fit_all  <- survivalSub(data_survival_fitting,
#'                                  form_marginal_surv, form_conditional_cr)
#' summary(survival_fit_all)
#'
#' @keywords internal
#' @export
summary.survivalSub.BJM <- function(object, digits = 4, ...) {
  .format_survivalSub(object, digits = digits, extended = TRUE)

  out <- list(
    cox_summary = if (!is_interval_censored(object)) summary(object$coxph_fit),
    glm_summary = if (!is.null(object$glm_fit)) summary(object$glm_fit) else NULL
  )
  if (is_interval_censored(object)) out$ic_coefficients <- icph_coef_table(object$ic_fit)
  class(out) <- "summary.survivalSub.BJM"
  invisible(out)
}


# -- Internal formatting helper for longitudinalSub (not exported) --------------
.format_longitudinalSub <- function(x, digits = 4, extended = FALSE) {

  sep_line  <- paste(rep("=", 65), collapse = "")
  dash_line <- paste(rep("-", 65), collapse = "")

  lfit         <- x$lfit
  Sigma_fit    <- x$Sigma_fit
  long_sub_fixed <- x$long_sub_fixed
  M            <- length(lfit)
  ### biomarker_type is only present on fits from longitudinalSubCopula();
  ### absent (NULL) means every biomarker is continuous, i.e. every lfit[[m]]
  ### is an nlme::lme object, exactly as before biomarker_type existed.
  biomarker_type <- x$biomarker_type
  if (is.null(biomarker_type)) biomarker_type <- rep("continuous", M)

  cat("\n")
  cat("Call:\n")
  cat(sprintf("longitudinalSub(M = %d longitudinal outcome%s)\n",
              M, if (M > 1) "s" else ""))

  if (extended) {
    cat("\nData Descriptives:\n")
    for (m in seq_len(M)) {
      nm <- if (!is.null(names(long_sub_fixed)[m]) && names(long_sub_fixed)[m] != "")
              names(long_sub_fixed)[m] else paste0("Outcome ", m)
      if (biomarker_type[m] == "continuous") {
        n_subj <- lfit[[m]]$dims$ngrps[1]
        n_obs  <- lfit[[m]]$dims$N
      } else {
        n_subj <- unname(lfit[[m]]$dims$nlev.re[1])
        n_obs  <- lfit[[m]]$dims$nobs
      }
      cat(sprintf("  [%d] %-15s  subjects = %d,  observations = %d\n",
                  m, nm, n_subj, n_obs))
    }
  }

  cat("\n", sep_line, "\n", sep = "")
  cat(" Longitudinal Sub-model(s)\n")

  for (m in seq_len(M)) {
    fit_m  <- lfit[[m]]
    form_m <- long_sub_fixed[[m]]
    nm     <- if (!is.null(names(long_sub_fixed)[m]) && names(long_sub_fixed)[m] != "")
                names(long_sub_fixed)[m] else paste0("Outcome ", m)

    cat(dash_line, "\n", sep = "")
    cat(sprintf(" [%d] %s (%s)\n", m, nm, biomarker_type[m]))
    cat(" Formula: ")
    print(form_m)
    cat("\n")

    if (biomarker_type[m] == "continuous") {
      fe   <- nlme::fixef(fit_m)
      se   <- sqrt(diag(fit_m$varFix))
      tval <- fe / se
      pval <- 2 * stats::pt(abs(tval), df = fit_m$fixDF$terms[1], lower.tail = FALSE)
      tab_fe <- cbind(Value = fe, SE = se, t = tval, "p-value" = pval)
      stats::printCoefmat(tab_fe,
                   digits       = digits,
                   P.values     = TRUE,
                   has.Pvalue   = TRUE,
                   signif.stars = getOption("show.signif.stars"),
                   cs.ind       = 1:2,
                   tst.ind      = 3)

      cat("\n")
      cat(sprintf("  Residual std. error (sigma): %.4f\n", fit_m$sigma))
      cat(sprintf("  n (subjects) = %d,  N (observations) = %d\n",
                  fit_m$dims$ngrps[1], fit_m$dims$N))
      cat(sprintf("  Log-likelihood: %.2f\n", stats::logLik(fit_m)[1]))
      cat(sprintf("  AIC: %.2f,  BIC: %.2f\n", stats::AIC(fit_m), stats::BIC(fit_m)))

      if (extended) {
        vc <- nlme::getVarCov(fit_m)
        cat("\n  Random-effects variance-covariance (subject level):\n")
        print(round(vc, digits))
      }
    } else {
      ### Ordinal biomarker: fit_m is an ordinal::clmm object (probit
      ### cumulative link mixed model). Its fixed effects (fit_m$beta) omit
      ### an intercept by construction -- absorbed into the cumulative
      ### thresholds (fit_m$alpha) instead -- so thresholds are reported
      ### alongside, not folded into the same coefficient table.
      fe   <- fit_m$beta
      se   <- sqrt(diag(stats::vcov(fit_m))[names(fe)])
      zval <- fe / se
      pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)
      tab_fe <- cbind(Value = fe, SE = se, z = zval, "p-value" = pval)
      stats::printCoefmat(tab_fe,
                   digits       = digits,
                   P.values     = TRUE,
                   has.Pvalue   = TRUE,
                   signif.stars = getOption("show.signif.stars"),
                   cs.ind       = 1:2,
                   tst.ind      = 3)

      cat("\n  Cumulative thresholds (probit link):\n")
      print(round(fit_m$alpha, digits))

      cat("\n")
      cat(sprintf("  n (subjects) = %d,  N (observations) = %d\n",
                  unname(fit_m$dims$nlev.re[1]), fit_m$dims$nobs))
      cat(sprintf("  Log-likelihood: %.2f\n", stats::logLik(fit_m)[1]))
      cat(sprintf("  AIC: %.2f,  BIC: %.2f\n", stats::AIC(fit_m), stats::BIC(fit_m)))

      if (extended) {
        vc <- nlme::VarCorr(fit_m)[[1]]
        cat("\n  Random-effects variance-covariance (subject level):\n")
        print(round(vc, digits))
      }
    }
    cat("\n")
  }

  cat(sep_line, "\n", sep = "")
  cat(" Multivariate Random-Effects Covariance Matrix (D)\n")
  cat(dash_line, "\n", sep = "")
  print(round(Sigma_fit, digits))

  if (extended) {
    cat("\n  Correlations:\n")
    D <- Sigma_fit
    corr_mat <- diag(1 / sqrt(diag(D))) %*% D %*% diag(1 / sqrt(diag(D)))
    dimnames(corr_mat) <- dimnames(D)
    print(round(corr_mat, digits))
  }

  cat(sep_line, "\n\n", sep = "")
}


#' Print method for \code{longitudinalSub.BJM} objects
#'
#' Automatically called when you type \code{long_fit_all} at the console.
#'
#' @param x A \code{longitudinalSub.BJM} object returned by \code{\link{longitudinalSub}}.
#' @param digits Number of significant digits. Default is 4.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @keywords internal
#' @export
print.longitudinalSub.BJM <- function(x, digits = 4, ...) {
  .format_longitudinalSub(x, digits = digits, extended = FALSE)
  invisible(x)
}


#' Summary method for \code{longitudinalSub.BJM} objects
#'
#' Like \code{print} but adds per-outcome random-effects variance components
#' and the full correlation matrix of D.
#'
#' @param object A \code{longitudinalSub.BJM} object.
#' @param digits Number of significant digits. Default is 4.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns a list of per-outcome \code{summary.lme} objects.
#' @keywords internal
#' @export
summary.longitudinalSub.BJM <- function(object, digits = 4, ...) {
  .format_longitudinalSub(object, digits = digits, extended = TRUE)
  out <- lapply(object$lfit, summary)
  class(out) <- "summary.longitudinalSub.BJM"
  invisible(out)
}


# -- Internal convenience wrapper: prints both sub-model fits together (not
# exported; long_fit_all and survival_fit_all already auto-print themselves
# via print.longitudinalSub.BJM/print.survivalSub.BJM, so this saves callers
# who have both objects on hand one line -- not a distinct public capability) --
printBJM <- function(long_fit_all, survival_fit_all, digits = 4) {
  cat("\n")
  cat("Backward Joint Model (BJM) - Model Summary\n")
  print.longitudinalSub.BJM(long_fit_all,  digits = digits)
  print.survivalSub.BJM(survival_fit_all,  digits = digits)
  invisible(list(long_fit_all = long_fit_all,
                 survival_fit_all = survival_fit_all))
}


# -- Internal formatting helper for predictRisk (not exported) ------------
.format_predictRisk <- function(x, digits = 4,
                                      prediction_time = NULL,
                                      horizon = NULL,
                                      subject_ids = NULL,
                                      extended = FALSE) {

  sep_line  <- paste(rep("=", 65), collapse = "")
  dash_line <- paste(rep("-", 65), collapse = "")

  risk0  <- x$risk_prob_1
  risk1  <- x$risk_prob_2
  has_cr <- !is.null(risk1)
  n_subj <- length(risk0)
  ids    <- if (!is.null(subject_ids)) as.character(subject_ids) else
            if (!is.null(names(risk0))) names(risk0) else
              paste0("S", seq_len(n_subj))

  cat("\n", sep_line, "\n", sep = "")
  cat(" Dynamic Prediction - Event Risk\n")
  cat(dash_line, "\n", sep = "")
  if (!is.null(prediction_time))
    cat(sprintf("  Prediction time   : %g\n", prediction_time))
  if (!is.null(horizon))
    cat(sprintf("  Prediction horizon: %g\n", horizon))
  if (!is.null(prediction_time) && !is.null(horizon))
    cat(sprintf("  Risk window       : (%g, %g]\n",
                prediction_time, prediction_time + horizon))
  cat(sprintf("  Competing risks   : %s\n", ifelse(has_cr, "Yes", "No")))
  cat(sprintf("  Subjects          : %d\n", n_subj))
  cat(dash_line, "\n\n", sep = "")

  if (!has_cr) {
    tab <- data.frame(Subject   = ids,
                      Risk_Prob = round(risk0, digits),
                      stringsAsFactors = FALSE)
    colnames(tab) <- c("Subject", "Risk Prob")
    print(tab, row.names = FALSE, right = TRUE)
  } else {
    tab <- data.frame(Subject = ids,
                      Risk0   = round(risk0, digits),
                      Risk1   = round(risk1, digits),
                      Total   = round(risk0 + risk1, digits),
                      stringsAsFactors = FALSE)
    colnames(tab) <- c("Subject", "Cause 1 Risk", "Cause 2 Risk", "Total Risk")
    print(tab, row.names = FALSE, right = TRUE)
  }

  if (extended) {
    cat("\n")
    cat(dash_line, "\n", sep = "")
    cat(" Summary Statistics\n")
    cat(dash_line, "\n", sep = "")
    if (!has_cr) {
      cat(sprintf("  Mean : %.4f\n", mean(risk0, na.rm = TRUE)))
      cat(sprintf("  SD   : %.4f\n", stats::sd(risk0,   na.rm = TRUE)))
      cat(sprintf("  Min / Max : %.4f / %.4f\n",
                  min(risk0, na.rm = TRUE), max(risk0, na.rm = TRUE)))
    } else {
      cat(sprintf("  Mean (Cause 1 / Cause 2 / Total): %.4f / %.4f / %.4f\n",
                  mean(risk0, na.rm = TRUE),
                  mean(risk1, na.rm = TRUE),
                  mean(risk0 + risk1, na.rm = TRUE)))
      cat(sprintf("  SD   (Cause 1 / Cause 2 / Total): %.4f / %.4f / %.4f\n",
                  stats::sd(risk0, na.rm = TRUE),
                  stats::sd(risk1, na.rm = TRUE),
                  stats::sd(risk0 + risk1, na.rm = TRUE)))
    }
  }

  cat("\n", sep_line, "\n\n", sep = "")
}


#' Print method for \code{predictRisk.BJM} objects
#'
#' Automatically called when you type the result of \code{predictRisk()}
#' at the console.
#'
#' @param x A \code{predictRisk.BJM} object.
#' @param prediction_time Landmark time (for display). Default \code{NULL}.
#' @param horizon Prediction horizon (for display). Default \code{NULL}.
#' @param subject_ids Optional subject ID labels.
#' @param digits Decimal places. Default 4.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @keywords internal
#' @export
print.predictRisk.BJM <- function(x, prediction_time = NULL,
                                        horizon = NULL,
                                        subject_ids = NULL,
                                        digits = 4, ...) {
  .format_predictRisk(x, digits = digits,
                             prediction_time = prediction_time,
                             horizon = horizon,
                             subject_ids = subject_ids,
                             extended = FALSE)
  invisible(x)
}


#' Summary method for \code{predictRisk.BJM} objects
#'
#' Like \code{print} but also shows mean, SD, and range of predicted risks.
#'
#' @param object A \code{predictRisk.BJM} object.
#' @param prediction_time Landmark time (for display). Default \code{NULL}.
#' @param horizon Prediction horizon (for display). Default \code{NULL}.
#' @param subject_ids Optional subject ID labels.
#' @param digits Decimal places. Default 4.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{object}.
#' @keywords internal
#' @export
summary.predictRisk.BJM <- function(object, prediction_time = NULL,
                                          horizon = NULL,
                                          subject_ids = NULL,
                                          digits = 4, ...) {
  .format_predictRisk(object, digits = digits,
                             prediction_time = prediction_time,
                             horizon = horizon,
                             subject_ids = subject_ids,
                             extended = TRUE)
  invisible(object)
}


# -- Internal formatting helper for predictLongitudinal (not exported) ----------
.format_predictLongitudinal <- function(x, digits = 4,
                                        bio_i = NULL,
                                        long_fit_all = NULL,
                                        prediction_time = NULL,
                                        horizon = NULL,
                                        subject_ids = NULL,
                                        extended = FALSE) {

  sep_line  <- paste(rep("=", 65), collapse = "")
  dash_line <- paste(rep("-", 65), collapse = "")

  Y_predict <- x$Y_predict
  Y_density <- x$Y_density
  Y_all     <- x$Y_all
  n_subj    <- length(Y_predict)
  ids       <- if (!is.null(subject_ids)) as.character(subject_ids) else
               if (!is.null(names(Y_predict))) names(Y_predict) else
                 paste0("S", seq_len(n_subj))

  bio_name <- if (!is.null(bio_i) && !is.null(long_fit_all)) {
    tryCatch(as.character(formula(long_fit_all$long_sub_fixed[[bio_i]])[[2]]),
             error = function(e) paste0("Biomarker ", bio_i))
  } else if (!is.null(bio_i)) {
    paste0("Biomarker ", bio_i)
  } else "Biomarker"

  cat("\n", sep_line, "\n", sep = "")
  cat(sprintf(" Dynamic Prediction - Future %s\n", bio_name))
  cat(dash_line, "\n", sep = "")
  if (!is.null(prediction_time))
    cat(sprintf("  Prediction time   : %g\n", prediction_time))
  if (!is.null(horizon))
    cat(sprintf("  Prediction horizon: %g\n", horizon))
  if (!is.null(prediction_time) && !is.null(horizon))
    cat(sprintf("  Predicted at      : t = %g\n", prediction_time + horizon))
  cat(sprintf("  Y grid range      : [%.4g, %.4g]  (%d points)\n",
              min(Y_all), max(Y_all), length(Y_all)))
  cat(sprintf("  Subjects          : %d\n", n_subj))
  cat(dash_line, "\n\n", sep = "")

  get_stats <- function(i) {
    dens <- Y_density[, i]
    w    <- dens / (sum(dens) + 1e-300)
    E_Y  <- sum(w * Y_all)
    SD_Y <- sqrt(max(sum(w * Y_all^2) - E_Y^2, 0))
    cw   <- cumsum(w)
    lo   <- Y_all[which(cw >= 0.025)[1]]
    hi   <- Y_all[which(cw >= 0.975)[1]]
    c(Mean = E_Y, SD = SD_Y, Lower = lo, Upper = hi)
  }

  stats_mat <- t(sapply(seq_len(n_subj), get_stats))

  tab <- data.frame(
    Subject = ids,
    MAP     = round(Y_predict,              digits),
    Mean    = round(stats_mat[, "Mean"],    digits),
    SD      = round(stats_mat[, "SD"],      digits),
    L95     = round(stats_mat[, "Lower"],   digits),
    U95     = round(stats_mat[, "Upper"],   digits),
    stringsAsFactors = FALSE
  )
  colnames(tab) <- c("Subject",
                     sprintf("MAP(%s)", bio_name),
                     "Post.Mean", "Post.SD", "2.5%", "97.5%")
  print(tab, row.names = FALSE, right = TRUE)

  cat(sprintf("\n  MAP = mode of the predictive density (Maximum A Posteriori)\n"))
  cat(sprintf("  Posterior summaries from density grid of %d points.\n",
              length(Y_all)))

  if (extended) {
    cat("\n")
    cat(dash_line, "\n", sep = "")
    cat(" Predictive Distribution Summary (across subjects)\n")
    cat(dash_line, "\n", sep = "")
    cat(sprintf("  MAP       - Mean (SD): %.4f (%.4f),  Range: [%.4f, %.4f]\n",
                mean(Y_predict, na.rm = TRUE), stats::sd(Y_predict, na.rm = TRUE),
                min(Y_predict,  na.rm = TRUE), max(Y_predict, na.rm = TRUE)))
    cat(sprintf("  Post.Mean - Mean (SD): %.4f (%.4f)\n",
                mean(stats_mat[, "Mean"]), stats::sd(stats_mat[, "Mean"])))
    cat(sprintf("  Post.SD   - Mean (SD): %.4f (%.4f)\n",
                mean(stats_mat[, "SD"]),   stats::sd(stats_mat[, "SD"])))
    cat(sprintf("  95%% CI width - Mean: %.4f\n",
                mean(stats_mat[, "Upper"] - stats_mat[, "Lower"])))
  }

  cat("\n", sep_line, "\n\n", sep = "")
}


#' Print method for \code{predictLongitudinal.BJM} objects
#'
#' Automatically called when you type the result of
#' \code{\link{predictLongitudinal}()} (with a single \code{bio_i}) at the
#' console.
#'
#' @param x A \code{predictLongitudinal.BJM} object.
#' @param bio_i Biomarker index (for label lookup). Default \code{NULL}.
#' @param long_fit_all \code{longitudinalSub.BJM} object for name lookup.
#' @param prediction_time Landmark time (for display). Default \code{NULL}.
#' @param horizon Prediction horizon (for display). Default \code{NULL}.
#' @param subject_ids Optional subject ID labels.
#' @param digits Decimal places. Default 4.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @keywords internal
#' @export
print.predictLongitudinal.BJM <- function(x, bio_i = NULL,
                                          long_fit_all = NULL,
                                          prediction_time = NULL,
                                          horizon = NULL,
                                          subject_ids = NULL,
                                          digits = 4, ...) {
  .format_predictLongitudinal(x, digits = digits,
                              bio_i = bio_i,
                              long_fit_all = long_fit_all,
                              prediction_time = prediction_time,
                              horizon = horizon,
                              subject_ids = subject_ids,
                              extended = FALSE)
  invisible(x)
}


#' Summary method for \code{predictLongitudinal.BJM} objects
#'
#' Like \code{print} but also shows distribution-level summaries across subjects.
#'
#' @param object A \code{predictLongitudinal.BJM} object.
#' @param bio_i Biomarker index (for label lookup). Default \code{NULL}.
#' @param long_fit_all \code{longitudinalSub.BJM} object for name lookup.
#' @param prediction_time Landmark time (for display). Default \code{NULL}.
#' @param horizon Prediction horizon (for display). Default \code{NULL}.
#' @param subject_ids Optional subject ID labels.
#' @param digits Decimal places. Default 4.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{object}.
#' @keywords internal
#' @export
summary.predictLongitudinal.BJM <- function(object, bio_i = NULL,
                                            long_fit_all = NULL,
                                            prediction_time = NULL,
                                            horizon = NULL,
                                            subject_ids = NULL,
                                            digits = 4, ...) {
  .format_predictLongitudinal(object, digits = digits,
                              bio_i = bio_i,
                              long_fit_all = long_fit_all,
                              prediction_time = prediction_time,
                              horizon = horizon,
                              subject_ids = subject_ids,
                              extended = TRUE)
  invisible(object)
}
