#' Pool longitudinal sub-model fits across multiple imputations with Rubin's rules
#'
#' @description When \code{\link{imputeLongitudinal}} is run with
#' \code{impute = "multiple"}, it returns \code{n_imputations} independently
#' completed versions of \code{data_fit_all} in \code{data_fit_all_list},
#' instead of a single completed dataset filled with the across-draw mean
#' (\code{impute = "single"}). Fitting \code{\link{longitudinalSub}} once per
#' completion and combining the resulting fixed-effect estimates with
#' Rubin's rules (Rubin, 1987) -- rather than averaging into one completed
#' dataset up front -- is the standard way multiple imputation propagates
#' the extra uncertainty from not knowing the true missing values into the
#' final standard errors; \code{impute = "single"} does not do this (see
#' the "Optional: imputing interrupted follow-up before fitting" vignette
#' section). \code{poolLongitudinalSub()} performs that combination step.
#'
#' @details For each fixed-effect coefficient, across the \code{m} fits in
#' \code{long_fit_all_list}:
#' \itemize{
#'   \item the pooled estimate is the mean of the \code{m} per-completion
#'   estimates (\eqn{\bar{Q}});
#'   \item the pooled variance (\eqn{T}) adds the average within-imputation
#'   variance (\eqn{\bar{U}}, from each fit's own standard error) to the
#'   between-imputation variance (\eqn{B}, the sample variance of the
#'   \code{m} estimates) inflated by a factor of \eqn{1 + 1/m}:
#'   \eqn{T = \bar{U} + (1 + 1/m) B};
#'   \item the pooled degrees of freedom use the Barnard & Rubin (1999)
#'   adjustment, which interpolates between the classical Rubin (1987)
#'   degrees of freedom (based only on \code{m} and the fraction of missing
#'   information) and the complete-data degrees of freedom reported by each
#'   \code{nlme::lme()} fit. This adjustment matters most at the small
#'   sample sizes typical of this package's use case, where the classical,
#'   unadjusted degrees of freedom can otherwise come out larger than the
#'   complete-data degrees of freedom, which is not sensible.
#' }
#' A large fraction of missing information (\code{fmi}, close to 1) for a
#' coefficient means the missing data cost a lot of precision for that
#' coefficient specifically, even if the biomarker's overall missingness
#' rate is modest.
#'
#' @param long_fit_all_list A list of two or more \code{longitudinalSub.BJM}
#' objects: the result of calling \code{\link{longitudinalSub}} once on
#' each completed dataset in \code{imputeLongitudinal(..., impute =
#' "multiple")$data_fit_all_list}, using the same \code{long_sub_fixed}/
#' \code{long_sub_random} for every completion.
#'
#' @return An object of class \code{"poolLongitudinalSub.BJM"}, a named list
#' with elements:
#' \describe{
#'   \item{pooled}{A list, one \code{data.frame} per longitudinal outcome,
#'   with one row per fixed-effect coefficient and columns \code{estimate},
#'   \code{std_error}, \code{df}, \code{statistic}, \code{p_value},
#'   \code{riv} (relative increase in variance due to missingness), and
#'   \code{fmi} (fraction of missing information).}
#'   \item{m}{The number of completed-data fits pooled.}
#'   \item{long_sub_fixed}{The \code{long_sub_fixed} used by every fit (taken
#'   from the first element of \code{long_fit_all_list}).}
#'   \item{long_fit_all}{A \code{longitudinalSub.BJM} object usable for
#'   prediction (e.g. as \code{long_fit_all} in \code{\link{predictRisk}}):
#'   the first fit, with each biomarker's fixed effects replaced by their
#'   pooled estimates, its residual variance by the average across
#'   completions, and \code{Sigma_fit} by the average of the \code{m} fits'
#'   \code{Sigma_fit} (Rubin's rules point estimates). Averaging the \code{m}
#'   fits' predictions instead is a common alternative.}
#' }
#'
#' @references Rubin, D. B. (1987). \emph{Multiple Imputation for
#' Nonresponse in Surveys}. John Wiley & Sons.
#' @references Barnard, J. and Rubin, D. B. (1999). Small-Sample Degrees of
#' Freedom with Multiple Imputation. \emph{Biometrika}, 86(4):948-955.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE)) {
#'   data(pbc3)
#'   data_fit_all <- pbc3[pbc3$status3 == 1, ]
#'
#'   set.seed(1)
#'   n <- nrow(data_fit_all)
#'   data_fit_all$serBilir[sample.int(n, floor(0.1 * n))] <- NA
#'   data_fit_all$albumin[sample.int(n, floor(0.1 * n))] <- NA
#'
#'   long_sub_fixed <- list(
#'     "long1" = serBilir ~ year + age + sex + years,
#'     "long2" = albumin ~ year + age + sex + years)
#'   long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
#'
#'   imputed <- imputeLongitudinal(data_fit_all, long_sub_fixed, long_sub_random,
#'                                  time_variable = "year", n_imputations = 5,
#'                                  impute = "multiple", epochs = 50, seed = 1)
#'
#'   long_fit_all_list <- lapply(imputed$data_fit_all_list, function(d) {
#'     longitudinalSub(d, long_sub_fixed, long_sub_random)
#'   })
#'   pooled <- poolLongitudinalSub(long_fit_all_list)
#'   pooled
#' }
#' }
#'
#' @export
poolLongitudinalSub <- function(long_fit_all_list) {
  assert_poolable_longitudinal_fits(long_fit_all_list)

  m <- length(long_fit_all_list)
  M <- length(long_fit_all_list[[1]]$lfit)
  biomarker_names <- names(long_fit_all_list[[1]]$long_sub_fixed)
  if (is.null(biomarker_names) || any(biomarker_names == "")) {
    biomarker_names <- paste0("Outcome ", seq_len(M))
  }

  pooled <- vector("list", M)
  for (k in seq_len(M)) {
    fits_k <- lapply(long_fit_all_list, function(f) f$lfit[[k]])
    coef_names <- names(nlme::fixef(fits_k[[1]]))

    ### vapply() drops the matrix dimension (returning a plain length-m
    ### vector instead of a 1 x m matrix) when there is exactly one
    ### coefficient -- e.g. an intercept-only long_sub_fixed -- so the
    ### shape is forced explicitly rather than relying on vapply's default
    n_coef <- length(coef_names)
    est_mat <- matrix(vapply(fits_k, function(f) nlme::fixef(f), numeric(n_coef)),
                       nrow = n_coef, dimnames = list(coef_names, NULL))
    se_mat <- matrix(vapply(fits_k, function(f) sqrt(diag(f$varFix)), numeric(n_coef)),
                      nrow = n_coef, dimnames = list(coef_names, NULL))
    ### complete-data df per coefficient: lme() gives each fixed effect its
    ### own df (fixDF$X) -- e.g. a between-subject covariate such as age has
    ### far fewer than the within-subject time terms. The intercept's df
    ### used to be applied to every coefficient, overstating the df (and
    ### understating the p-value) of between-subject effects.
    dfcom_mat <- matrix(vapply(fits_k, function(f) {
      d <- f$fixDF$X
      if (is.null(d)) rep(f$fixDF$terms[1], n_coef) else as.numeric(d[coef_names])
    }, numeric(n_coef)), nrow = n_coef)

    ### the m completions from a single imputeLongitudinal(impute =
    ### "multiple") call fill the same missing cells (never drop/add rows),
    ### so every completed dataset has the same complete-data df; a
    ### mismatch usually means these fits came from differently-shaped
    ### data (e.g. a caller mixing in an unrelated fit) -- use the smallest
    ### (most conservative) value and warn, rather than silently picking one
    dfcom <- apply(dfcom_mat, 1, min)
    if (any(apply(dfcom_mat, 1, function(d) any(d != d[1])))) {
      warning(sprintf(paste0(
        "Outcome %d ('%s'): complete-data degrees of freedom differ across ",
        "`long_fit_all_list` (%s); using the smallest (most conservative) ",
        "value per coefficient. This is expected only if the completed datasets do not ",
        "all have the same number of retained subjects/observations."
      ), k, biomarker_names[k], paste(unique(c(dfcom_mat)), collapse = ", ")), call. = FALSE)
    }

    rows <- lapply(seq_along(coef_names), function(j) {
      rubin_pool_scalar(est_mat[j, ], se_mat[j, ], dfcom = dfcom[j])
    })
    tab <- do.call(rbind, rows)
    rownames(tab) <- coef_names
    pooled[[k]] <- tab
  }
  names(pooled) <- biomarker_names

  out <- list(pooled = pooled, m = m,
              long_sub_fixed = long_fit_all_list[[1]]$long_sub_fixed,
              long_fit_all = pooled_longitudinal_fit(long_fit_all_list, pooled))
  class(out) <- "poolLongitudinalSub.BJM"
  out
}

#' Combine one coefficient's per-completion estimates with Rubin's rules
#'
#' @description Internal helper for \code{\link{poolLongitudinalSub}}:
#' implements the scalar-parameter pooling rules of Rubin (1987, Ch. 3) and
#' the Barnard & Rubin (1999) small-sample degrees-of-freedom adjustment,
#' for a single fixed-effect coefficient's \code{m} per-completion
#' estimates/standard errors.
#' @param estimates Length-\code{m} numeric vector of per-completion point
#' estimates.
#' @param std_errors Length-\code{m} numeric vector of per-completion
#' standard errors (the square root of the within-imputation variance).
#' @param dfcom The complete-data degrees of freedom (assumed common across
#' completions; see \code{\link{poolLongitudinalSub}}).
#' @keywords internal
rubin_pool_scalar <- function(estimates, std_errors, dfcom) {
  m <- length(estimates)
  qbar <- mean(estimates)
  ubar <- mean(std_errors^2)
  b <- stats::var(estimates)
  t_var <- ubar + (1 + 1 / m) * b
  se <- sqrt(t_var)

  ### relative increase in variance / proportion of variance attributable
  ### to missingness; guarded against ubar == 0 (e.g. a coefficient that
  ### happens to be estimated with zero within-imputation variance)
  riv <- if (ubar > 0) (1 + 1 / m) * b / ubar else Inf
  lambda <- (1 + 1 / m) * b / t_var

  df <- barnard_rubin_df(m, lambda, dfcom)
  ### Barnard-Rubin FMI uses the pooled (adjusted) df, as in mice::pool();
  ### it used the complete-data df dfcom before
  fmi <- if (is.finite(riv)) (riv + 2 / (df + 3)) / (riv + 1) else 1

  statistic <- qbar / se
  p_value <- 2 * stats::pt(-abs(statistic), df = df)

  data.frame(estimate = qbar, std_error = se, df = df,
             statistic = statistic, p_value = p_value,
             riv = riv, fmi = fmi)
}

#' Barnard & Rubin (1999) small-sample pooled degrees of freedom
#'
#' @description Internal helper for \code{\link{rubin_pool_scalar}}.
#' Interpolates between the classical Rubin (1987) degrees of freedom
#' (\code{df_old}, which ignores the complete-data sample size and can
#' therefore exceed it) and the complete-data degrees of freedom
#' (\code{dfcom}), weighted by how much of the total variance is due to
#' missingness (\code{lambda}).
#' @keywords internal
barnard_rubin_df <- function(m, lambda, dfcom) {
  ### clamp away from exactly 0 (a coefficient with identical estimates
  ### across every completion) to avoid a division by zero in df_old; this
  ### mirrors the standard practical safeguard used by other Rubin's-rules
  ### implementations, and simply caps df_old very high rather than at Inf
  lambda <- max(lambda, 1e-4)
  df_old <- (m - 1) / lambda^2
  if (!is.finite(dfcom)) return(df_old)
  df_obs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
  df_old * df_obs / (df_old + df_obs)
}

#' Print method for \code{poolLongitudinalSub.BJM} objects
#'
#' Automatically called when you type \code{pooled} or
#' \code{print(pooled)} at the console. Displays a JMbayes2-style
#' formatted summary of the Rubin's-rules-pooled longitudinal sub-model(s).
#'
#' @param x A \code{poolLongitudinalSub.BJM} object returned by
#' \code{\link{poolLongitudinalSub}}.
#' @param digits Number of significant digits. Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE)) {
#'   data(pbc3)
#'   data_fit_all <- pbc3[pbc3$status3 == 1, ]
#'
#'   set.seed(1)
#'   n <- nrow(data_fit_all)
#'   data_fit_all$serBilir[sample.int(n, floor(0.1 * n))] <- NA
#'   data_fit_all$albumin[sample.int(n, floor(0.1 * n))] <- NA
#'
#'   long_sub_fixed <- list(
#'     "long1" = serBilir ~ year + age + sex + years,
#'     "long2" = albumin ~ year + age + sex + years)
#'   long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
#'
#'   imputed <- imputeLongitudinal(data_fit_all, long_sub_fixed, long_sub_random,
#'                                  time_variable = "year", n_imputations = 5,
#'                                  impute = "multiple", epochs = 50, seed = 1)
#'
#'   long_fit_all_list <- lapply(imputed$data_fit_all_list, function(d) {
#'     longitudinalSub(d, long_sub_fixed, long_sub_random)
#'   })
#'   pooled <- poolLongitudinalSub(long_fit_all_list)
#'   pooled   # triggers print.poolLongitudinalSub.BJM automatically
#' }
#' }
#'
#' @keywords internal
#' @export
print.poolLongitudinalSub.BJM <- function(x, digits = 4, ...) {
  sep_line  <- paste(rep("=", 65), collapse = "")
  dash_line <- paste(rep("-", 65), collapse = "")
  M <- length(x$pooled)

  cat("\n")
  cat("Call:\n")
  cat(sprintf("poolLongitudinalSub(m = %d completed-data fits, M = %d longitudinal outcome%s)\n",
              x$m, M, if (M > 1) "s" else ""))

  cat("\n", sep_line, "\n", sep = "")
  cat(" Rubin's-Rules-Pooled Longitudinal Sub-model(s)\n")

  for (k in seq_len(M)) {
    nm  <- names(x$pooled)[k]
    tab <- x$pooled[[k]]

    cat(dash_line, "\n", sep = "")
    cat(sprintf(" [%d] %s\n\n", k, nm))

    coef_tab <- as.matrix(tab[, c("estimate", "std_error", "df", "statistic", "p_value")])
    colnames(coef_tab) <- c("Value", "SE", "DF", "t", "p-value")
    stats::printCoefmat(coef_tab, digits = digits, P.values = TRUE, has.Pvalue = TRUE,
                         signif.stars = getOption("show.signif.stars"),
                         cs.ind = 1:2, tst.ind = 4)

    cat("\n  Fraction of missing information (FMI), by coefficient:\n  ")
    cat(paste(sprintf("%s = %.3f", rownames(tab), tab$fmi), collapse = ",  "), "\n")
    cat("\n")
  }

  cat(sep_line, "\n", sep = "")
  invisible(x)
}

#' A prediction-ready fit from Rubin's-rules point estimates
#'
#' @description Helper for \code{poolLongitudinalSub()}: see its
#' \code{long_fit_all} return element.
#' @param long_fit_all_list The completed-data fits.
#' @param pooled \code{poolLongitudinalSub()}'s per-outcome pooled tables.
#' @return A \code{longitudinalSub.BJM} object.
#' @keywords internal
pooled_longitudinal_fit <- function(long_fit_all_list, pooled) {
  fit <- long_fit_all_list[[1]]
  for (k in seq_along(fit$lfit)) {
    est <- pooled[[k]]$estimate
    names(est) <- rownames(pooled[[k]])
    fit$lfit[[k]]$coefficients$fixed <- est
    sigmas <- vapply(long_fit_all_list, function(f) {
      s <- f$lfit[[k]]$sigma
      if (is.null(s)) NA_real_ else s
    }, numeric(1))
    if (all(is.finite(sigmas))) fit$lfit[[k]]$sigma <- sqrt(mean(sigmas^2))
  }
  Sigmas <- lapply(long_fit_all_list, `[[`, "Sigma_fit")
  if (!any(vapply(Sigmas, is.null, logical(1)))) {
    fit$Sigma_fit <- Reduce(`+`, Sigmas) / length(Sigmas)
  }
  fit
}
