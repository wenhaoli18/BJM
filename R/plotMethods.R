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


#' Plot simulated futures
#'
#' @description
#' Plots the draws of \code{\link{simulateTrajectory}}, pooling every draw of
#' the patients chosen by \code{id} (all of them by default: one patient's
#' draws picture that patient's predictive distribution, one draw each of a
#' synthetic cohort pictures the cohort):
#' \describe{
#'   \item{\code{"trajectory"}}{One panel per biomarker. For a continuous
#'   biomarker, \code{n_paths} of the simulated trajectories as thin lines,
#'   with the median and the central \code{level} interval of the draws at
#'   each time, and the measurements conditioned on as points. For an
#'   ordinal biomarker, the share of the draws in each category at each
#'   time. With \code{truncate = TRUE} in \code{simulateTrajectory()} the
#'   summaries at a time describe the draws still event-free then.}
#'   \item{\code{"event"}}{The cumulative incidence of the event (of each
#'   event type, with competing risks) after \code{prediction_time},
#'   computed from the drawn event times.}
#' }
#'
#' @param x A \code{simulateTrajectory.BJM} object returned by
#'   \code{\link{simulateTrajectory}}.
#' @param which \code{"trajectory"} (default) or \code{"event"}.
#' @param id Patient id(s) whose draws to plot; \code{NULL} (default) for
#'   all.
#' @param n_paths Number of simulated trajectories drawn as lines (the
#'   first ones; the draws are independent, so they are a random sample).
#' @param level Coverage of the interval drawn around the median.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' data(pbc3)
#' survival_fit_all <- survivalSub(pbc3[!duplicated(pbc3$id), ],
#'                                 Surv(years, status3) ~ age + sex, NULL)
#' long_fit_all <- longitudinalSub(list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ]),
#'                                 list(serBilir ~ year + age + sex + years,
#'                                      albumin ~ year + age + sex + years),
#'                                 list(~ year | id, ~ year | id))
#' history <- pbc3[pbc3$id == 2 & pbc3$year <= 3, ]
#' sims <- simulateTrajectory(history, long_fit_all, survival_fit_all,
#'                            prediction_time = 3, times = seq(3, 8, by = 0.5),
#'                            time_variable = "year", survival_variable_all = list(),
#'                            survival_trans_function = list(), n_sim = 200, seed = 1)
#' plot(sims)
#' plot(sims, which = "event")
#' }
#'
#' @export
plot.simulateTrajectory.BJM <- function(x, which = c("trajectory", "event"), id = NULL,
                                        n_paths = 30, level = 0.9, ...) {
  which <- match.arg(which)
  assert_positive_integer(n_paths, "n_paths")
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number between 0 and 1.", call. = FALSE)
  }
  vars <- attr(x, "variables")
  prediction_time <- attr(x, "prediction_time")
  max_event_time <- attr(x, "max_event_time")
  history <- attr(x, "history")
  sims <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  sims[[vars$id]] <- as.character(sims[[vars$id]])
  if (!is.null(id)) {
    unknown <- setdiff(as.character(id), sims[[vars$id]])
    if (length(unknown) > 0) {
      stop(sprintf("No draws for id(s) %s.", paste(unknown, collapse = ", ")), call. = FALSE)
    }
    sims <- sims[sims[[vars$id]] %in% as.character(id), , drop = FALSE]
    history <- history[as.character(history$id) %in% as.character(id), , drop = FALSE]
  }
  sims$path <- paste(sims[[vars$id]], sims$sim)

  if (which == "event") {
    draws <- sims[!duplicated(sims$path), , drop = FALSE]
    type <- if (is.null(vars$event_type)) rep("Event", nrow(draws)) else
      paste0(vars$event_type, " = ", draws[[vars$event_type]])
    has_event <- draws$event == 1
    end <- if (is.finite(max_event_time)) max_event_time else max(draws$event_time)
    grid <- sort(unique(c(prediction_time, draws$event_time[has_event], end)))
    cif <- do.call(rbind, lapply(sort(unique(type[has_event])), function(tp) {
      times <- draws$event_time[has_event & type == tp]
      data.frame(time = grid, cif = vapply(grid, function(t) sum(times <= t), numeric(1)) / nrow(draws),
                 curve = tp)
    }))
    p <- ggplot(cif, aes(x = time, y = cif)) +
      scale_y_continuous(limits = c(0, 1)) +
      xlab(vars$time) + ylab("Cumulative incidence") + theme_bw()
    p <- if (is.null(vars$event_type)) p + geom_step(linewidth = 1) else
      p + geom_step(aes(color = curve), linewidth = 1) + labs(color = NULL)
    if (is.finite(max_event_time)) {
      p <- p + labs(caption = sprintf("%d draws; %.1f%% event-free through %s = %s",
                                      nrow(draws), 100 * mean(!has_event), vars$time,
                                      format(signif(max_event_time, 3))))
    }
    return(p)
  }

  ordinal <- vapply(vars$biomarkers, function(b) is.factor(sims[[b]]), logical(1))
  shown_paths <- unique(sims$path)
  shown_paths <- shown_paths[seq_len(min(n_paths, length(shown_paths)))]
  alpha_tail <- (1 - level) / 2

  continuous_long <- do.call(rbind, lapply(vars$biomarkers[!ordinal], function(b)
    data.frame(path = sims$path, time = sims[[vars$time]], biomarker = b, value = sims[[b]])))
  continuous_long <- continuous_long[!is.na(continuous_long$value), , drop = FALSE]
  bands <- NULL
  if (!is.null(continuous_long) && nrow(continuous_long) > 0) {
    bands <- do.call(rbind, lapply(split(continuous_long, list(continuous_long$biomarker, continuous_long$time),
                                         drop = TRUE), function(d) {
      q <- stats::quantile(d$value, c(alpha_tail, 0.5, 1 - alpha_tail), names = FALSE)
      data.frame(time = d$time[1], biomarker = d$biomarker[1], lower = q[1], value = q[2], upper = q[3])
    }))
  }

  shares <- do.call(rbind, lapply(vars$biomarkers[ordinal], function(b) {
    d <- sims[!is.na(sims[[b]]), , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    counts <- table(d[[vars$time]], d[[b]])
    data.frame(time = as.numeric(rownames(counts))[row(counts)], biomarker = b,
               category = factor(colnames(counts)[col(counts)], levels = levels(sims[[b]])),
               value = as.vector(counts / rowSums(counts)))
  }))

  p <- ggplot() + geom_vline(xintercept = prediction_time, linetype = "dashed", color = "grey50")
  if (!is.null(bands)) {
    p <- p +
      geom_line(data = continuous_long[continuous_long$path %in% shown_paths, , drop = FALSE],
                aes(x = time, y = value, group = path), color = "grey60", alpha = 0.4, linewidth = 0.3) +
      geom_ribbon(data = bands, aes(x = time, ymin = lower, ymax = upper), fill = "#2166AC", alpha = 0.2) +
      geom_line(data = bands, aes(x = time, y = value), color = "#2166AC", linewidth = 1)
  }
  history <- history[history$biomarker %in% vars$biomarkers[!ordinal], , drop = FALSE]
  if (nrow(history) > 0) {
    p <- p + geom_point(data = history, aes(x = time, y = value), size = 1.8)
  }
  if (!is.null(shares)) {
    width <- if (length(unique(shares$time)) > 1) 0.8 * min(diff(sort(unique(shares$time)))) else 0.8
    p <- p + geom_col(data = shares, aes(x = time, y = value, fill = category), width = width) +
      labs(fill = NULL)
  }
  ### panels in the fit's biomarker order, not alphabetical
  for (layer in seq_along(p$layers)) {
    d <- p$layers[[layer]]$data
    if (is.data.frame(d) && "biomarker" %in% names(d)) {
      p$layers[[layer]]$data$biomarker <- factor(d$biomarker, levels = vars$biomarkers)
    }
  }
  p + facet_wrap(~ biomarker, scales = "free_y") +
    xlab(vars$time) +
    ylab(if (any(ordinal)) "Value, or share of draws in each category" else "Value") +
    labs(caption = sprintf("Line and band: median and %g%% interval of the draws%s; grey: %d draws.",
                           100 * level, if (any(is.na(sims[vars$biomarkers]))) " still event-free" else "",
                           length(shown_paths))) +
    theme_bw()
}
