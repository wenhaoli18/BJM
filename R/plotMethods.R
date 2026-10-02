#' Plot a fitted longitudinal sub-model
#'
#' @description
#' Diagnostic plots for the per-biomarker mixed models and the shared
#' random-effects covariance of a \code{\link{longitudinalSub}} fit, one
#' panel per biomarker:
#' \describe{
#'   \item{\code{"residuals"}}{Standardized (Pearson) residuals against
#'   fitted values, with a LOESS smoother; a trend or a funnel shape
#'   suggests a misspecified mean or non-constant variance.}
#'   \item{\code{"qq"}}{Normal Q-Q plot of the standardized residuals.}
#'   \item{\code{"ranef"}}{Normal Q-Q plot of each predicted random effect
#'   (the model assumes they are normal).}
#'   \item{\code{"corr"}}{Heat map of the correlation matrix of the
#'   multivariate random-effects covariance \code{Sigma_fit}, i.e. how the
#'   subject-level deviations of the different biomarkers move together.}
#' }
#' Residual plots are drawn for continuous biomarkers only; an ordinal
#' biomarker (fit by \code{ordinal::clmm()}) has no residuals on the
#' response scale and is left out of them with a message.
#'
#' @param x A \code{longitudinalSub.BJM} object returned by
#' \code{\link{longitudinalSub}}.
#' @param which Which plot to draw: \code{"residuals"} (default),
#' \code{"qq"}, \code{"ranef"} or \code{"corr"}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' data(pbc3)
#' long_sub_fixed <- list("long1" = serBilir ~ year + age + sex + years,
#'                        "long2" = albumin ~ year + age + sex + years)
#' long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
#' data_fit <- pbc3[pbc3$status3 == 1, ]
#' long_fit_all <- longitudinalSub(list(data_fit, data_fit), long_sub_fixed, long_sub_random)
#'
#' plot(long_fit_all)
#' plot(long_fit_all, which = "qq")
#' plot(long_fit_all, which = "ranef")
#' plot(long_fit_all, which = "corr")
#' }
#'
#' @export
plot.longitudinalSub.BJM <- function(x, which = c("residuals", "qq", "ranef", "corr"), ...) {
  which <- match.arg(which)
  lfit <- x$lfit
  M <- length(lfit)
  bio_names <- unname(vapply(x$long_sub_fixed, function(f) deparse(formula(f)[[2]]), character(1)))
  ### biomarker_type is only present on fits with an ordinal biomarker
  biomarker_type <- x$biomarker_type
  if (is.null(biomarker_type)) biomarker_type <- rep("continuous", M)

  if (which %in% c("residuals", "qq")) {
    continuous <- which(biomarker_type == "continuous")
    if (length(continuous) == 0) {
      stop("Residual plots need at least one continuous biomarker; every biomarker here is ordinal.",
           call. = FALSE)
    }
    if (length(continuous) < M) {
      message("Leaving out ordinal biomarker(s) without response-scale residuals: ",
              paste(bio_names[-continuous], collapse = ", "), ".")
    }
    res <- do.call(rbind, lapply(continuous, function(m) {
      data.frame(biomarker = bio_names[m],
                 fitted = unname(stats::fitted(lfit[[m]])),
                 resid = unname(stats::residuals(lfit[[m]], type = "pearson")))
    }))
    res$biomarker <- factor(res$biomarker, levels = bio_names[continuous])

    if (which == "residuals") {
      p <- ggplot(res, aes(x = fitted, y = resid)) +
        geom_point(alpha = 0.3, size = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "red") +
        facet_wrap(~ biomarker, scales = "free") +
        xlab("Fitted value") + ylab("Standardized residual")
    } else {
      p <- ggplot(res, aes(sample = resid)) +
        geom_qq(alpha = 0.3, size = 0.8) + geom_qq_line(color = "red") +
        facet_wrap(~ biomarker, scales = "free") +
        xlab("Theoretical normal quantile") + ylab("Standardized residual")
    }
    return(p + theme_bw())
  }

  if (which == "ranef") {
    re <- do.call(rbind, lapply(seq_len(M), function(m) {
      b <- nlme::ranef(lfit[[m]])
      ### lme gives a data.frame; clmm a list with one data.frame per grouping factor
      if (!is.data.frame(b)) b <- b[[1]]
      do.call(rbind, lapply(names(b), function(eff) {
        data.frame(panel = paste0(bio_names[m], ": ", eff), value = b[[eff]])
      }))
    }))
    re$panel <- factor(re$panel, levels = unique(re$panel))
    p <- ggplot(re, aes(sample = value)) +
      geom_qq(alpha = 0.5, size = 0.8) + geom_qq_line(color = "red") +
      facet_wrap(~ panel, scales = "free") +
      xlab("Theoretical normal quantile") + ylab("Predicted random effect")
    return(p + theme_bw())
  }

  ### which == "corr": label each row/column "<biomarker>: <effect>"; the
  ### fitted names end in "_<m>", the biomarker's position
  corr <- stats::cov2cor(as.matrix(x$Sigma_fit))
  raw <- rownames(corr)
  if (is.null(raw)) raw <- paste0("b", seq_len(nrow(corr)), "_", seq_len(nrow(corr)))
  m_idx <- suppressWarnings(as.integer(sub(".*_", "", raw)))
  labels <- ifelse(is.na(m_idx) | m_idx > M, raw,
                   paste0(bio_names[pmin(m_idx, M)], ": ", sub("_[0-9]+$", "", raw)))
  cd <- data.frame(row = factor(rep(labels, times = ncol(corr)), levels = rev(labels)),
                   col = factor(rep(labels, each = nrow(corr)), levels = labels),
                   value = as.vector(corr))
  ggplot(cd, aes(x = col, y = row, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Correlation") +
    xlab(NULL) + ylab(NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}


#' Plot a fitted survival sub-model
#'
#' @description
#' \describe{
#'   \item{\code{"forest"}}{Forest plot of the hazard ratios (with 95\%
#'   confidence intervals) of the Cox marginal survival model and, when the
#'   fit has competing risks, of the odds ratios (with 95\% Wald confidence
#'   intervals) of the logistic event-type model, on a log scale with a
#'   reference line at 1.}
#'   \item{\code{"basehaz"}}{The Cox model's baseline cumulative hazard
#'   (covariates at 0, not centered), one step curve per stratum.}
#' }
#'
#' @param x A \code{survivalSub.BJM} object returned by \code{\link{survivalSub}}.
#' @param which Which plot to draw: \code{"forest"} (default) or \code{"basehaz"}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(pbc3)
#' data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
#' survival_fit_all <- survivalSub(data_survival_fitting,
#'                                 Surv(years, status3) ~ age + sex,
#'                                 status4 ~ years + age + sex)
#' plot(survival_fit_all)
#' plot(survival_fit_all, which = "basehaz")
#'
#' @export
plot.survivalSub.BJM <- function(x, which = c("forest", "basehaz"), ...) {
  which <- match.arg(which)

  if (which == "basehaz") {
    bh <- basehaz(x$coxph_fit, centered = FALSE)
    bh <- rbind(if (is.null(bh$strata)) data.frame(time = 0, hazard = 0) else
                  data.frame(time = 0, hazard = 0, strata = unique(bh$strata)),
                bh[, intersect(c("time", "hazard", "strata"), names(bh)), drop = FALSE])
    p <- if (is.null(bh$strata)) {
      ggplot(bh, aes(x = time, y = hazard)) + geom_step(linewidth = 1)
    } else {
      ggplot(bh, aes(x = time, y = hazard, color = strata)) + geom_step(linewidth = 1)
    }
    return(p + xlab(all.vars(x$form_marginal_surv[[2]])[1]) +
             ylab("Baseline cumulative hazard") + theme_bw())
  }

  ci <- summary(x$coxph_fit)$conf.int
  est <- data.frame(model = "Survival model: hazard ratio", term = rownames(ci),
                    estimate = ci[, "exp(coef)"], lower = ci[, "lower .95"], upper = ci[, "upper .95"])
  if (!is.null(x$glm_fit)) {
    b <- stats::coef(x$glm_fit)
    se <- sqrt(diag(stats::vcov(x$glm_fit)))
    keep <- names(b) != "(Intercept)"
    est <- rbind(est, data.frame(
      model = "Competing-risks model: odds ratio", term = names(b)[keep],
      estimate = exp(b[keep]), lower = exp(b[keep] - stats::qnorm(0.975) * se[keep]),
      upper = exp(b[keep] + stats::qnorm(0.975) * se[keep])))
  }
  est$model <- factor(est$model, levels = unique(est$model))
  est$term <- factor(est$term, levels = rev(unique(est$term)))

  ggplot(est, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
    geom_point(size = 2.5) +
    scale_x_log10() +
    facet_wrap(~ model, ncol = 1, scales = "free") +
    xlab("Estimate (95% CI, log scale)") + ylab(NULL) +
    theme_bw()
}


#' Plot predicted biomarker distributions
#'
#' @description
#' Plots the predicted distribution of the biomarker at
#' \code{prediction_time + horizon} returned by \code{\link{predictLongitudinal}},
#' one curve per patient, with a dashed line at the patient's point
#' prediction (the most likely value, \code{Y_predict}). For an ordinal
#' biomarker, plots each patient's predicted category probabilities as bars
#' instead.
#'
#' @param x A \code{predictLongitudinal.BJM} object, i.e. the result of
#' \code{predictLongitudinal()} for a single biomarker. For several
#' biomarkers, plot one element of the returned list, e.g.
#' \code{plot(result[["serBilir"]])}.
#' @param subject Patients to plot, as ids (the names of \code{x$Y_predict})
#' or positions. \code{NULL} (default) plots every patient, or the first six
#' when there are more.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' data(pbc3)
#' survival_fit_all <- survivalSub(pbc3[!duplicated(pbc3$id), ],
#'                                 Surv(years, status3) ~ age + sex, NULL)
#' long_sub_fixed <- list("long1" = serBilir ~ year + age + sex + years)
#' long_sub_random <- list("long1" = ~ year | id)
#' long_fit_all <- longitudinalSub(list(pbc3[pbc3$status3 == 1, ]),
#'                                 long_sub_fixed, long_sub_random)
#'
#' data_predict <- pbc3[pbc3$id %in% c(2, 4, 5) & pbc3$year <= 3, ]
#' data_predict$years <- NA
#' trans <- survivalTrans(c(1, 3, 5, 7))
#' pred <- predictLongitudinal(data_predict, long_fit_all, survival_fit_all,
#'   prediction_time = 3, horizon = 1, time_variable = "year",
#'   trans$survival_variable_all, trans$survival_trans_function, bio_i = 1)
#' plot(pred)
#' }
#'
#' @export
plot.predictLongitudinal.BJM <- function(x, subject = NULL, ...) {
  ids <- names(x$Y_predict)
  if (is.null(ids)) ids <- as.character(seq_along(x$Y_predict))
  if (length(ids) == 0) {
    stop("`x` has no patients to plot (none were at risk at the prediction time).", call. = FALSE)
  }
  if (is.null(subject)) {
    keep <- seq_len(min(6, length(ids)))
    if (length(ids) > 6) {
      message(sprintf("Plotting the first 6 of %d patients; choose others with `subject`.", length(ids)))
    }
  } else if (is.numeric(subject)) {
    if (any(subject < 1 | subject > length(ids) | subject != round(subject))) {
      stop(sprintf("`subject` positions must be whole numbers between 1 and %d.", length(ids)), call. = FALSE)
    }
    keep <- subject
  } else {
    keep <- match(as.character(subject), ids)
    if (anyNA(keep)) {
      stop("Unknown patient id(s) in `subject`: ",
           paste(subject[is.na(keep)], collapse = ", "), ".", call. = FALSE)
    }
  }

  category_labels <- attr(x$Y_all, "category_labels")
  dens <- data.frame(
    value = rep(as.vector(x$Y_all), times = length(keep)),
    density = as.vector(x$Y_density[, keep, drop = FALSE]),
    subject = factor(rep(ids[keep], each = length(x$Y_all)), levels = ids[keep])
  )

  if (!is.null(category_labels)) {
    dens$category <- factor(category_labels[dens$value], levels = category_labels)
    levels(dens$subject) <- paste("Patient", levels(dens$subject))
    return(ggplot(dens, aes(x = category, y = density)) +
             geom_col(fill = "grey50") +
             facet_wrap(~ subject) +
             xlab("Predicted category") + ylab("Predicted probability") +
             theme_bw())
  }

  point <- data.frame(subject = factor(ids[keep], levels = ids[keep]),
                      value = as.vector(x$Y_predict[keep]))
  ggplot(dens, aes(x = value, y = density, color = subject)) +
    geom_line(linewidth = 1) +
    geom_vline(data = point, aes(xintercept = value, color = subject), linetype = "dashed") +
    labs(color = "Patient") +
    xlab("Predicted biomarker value") + ylab("Predicted density") +
    theme_bw()
}
