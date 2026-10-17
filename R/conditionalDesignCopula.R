#' Build the per-patient longitudinal design pieces for a mixed
#' continuous/ordinal (Gaussian-copula) joint model
#'
#' @description Copula-aware counterpart to \code{build_conditional_design()}
#' (see \code{R/conditionalDesign.R}), for use whenever
#' \code{long_fit_all$biomarker_type} contains at least one \code{"ordinal"}
#' biomarker (see \code{longitudinalSubCopula()}). Continuous biomarkers'
#' rows are handled identically to \code{build_conditional_design()}: their
#' exact observed value is recorded. For an ordinal biomarker's rows, the
#' exact underlying latent Gaussian score is unobserved -- only the pair of
#' cumulative-link thresholds (\code{long_fit_all$thresholds[[m]]}) bracketing
#' the observed category is known (the classical Albert & Chib probit
#' data-augmentation representation, same convention as
#' \code{longitudinalSubVarCopula()} in \code{R/copulaLongitudinal.R}: latent
#' score \code{Z = eta + epsilon}, \code{epsilon ~ N(0, 1)}, category
#' \code{k} iff \code{full_alpha[k] < Z <= full_alpha[k+1]} with
#' \code{full_alpha = c(-Inf, alpha_m, Inf)}, thresholds used directly with
#' no sign flip).
#'
#' The random-effects design (\code{A_i}) and \code{Sigma_all} construction
#' is unchanged in form from \code{build_conditional_design()} -- marker type
#' only affects \code{sigma.longitudinal} (already fixed at \code{1} for
#' ordinal markers by the caller, the probit identification constraint; see
#' \code{longitudinalSubCopula()}), which flows in exactly the same way for
#' both continuous and ordinal markers.
#'
#' @return A list with \code{y_all_vec} (stacked observed values; \code{NA}
#' at ordinal-marker rows), \code{alpha_lower_vec}/\code{alpha_upper_vec}
#' (stacked threshold bounds; \code{NA} at continuous-marker rows),
#' \code{row_marker_type} (stacked per-row \code{"continuous"}/\code{"ordinal"}
#' label), \code{parameter_matrix} (block-diagonal fixed-effect coefficients,
#' \code{lfit[[i]]$coefficients$fixed} for continuous markers,
#' \code{lfit[[i]]$beta} -- no intercept -- for ordinal markers), and
#' \code{Sigma_all}.
#' @keywords internal
build_conditional_design_copula <- function(rep_num_i_list, data_num_i_list, lfit, Sigma,
                                             sigma.longitudinal, time_variable, n_longitudinal,
                                             biomarker_type, long_sub_fixed, thresholds,
                                             long_sub_random) {
  n_total <- sum(sapply(data_num_i_list, nrow))
  y_all_vec <- rep(NA_real_, n_total)
  alpha_lower_vec <- rep(NA_real_, n_total)
  alpha_upper_vec <- rep(NA_real_, n_total)
  row_marker_type <- character(n_total)

  length_y <- rep(0, n_longitudinal + 1)
  for (i in 1:n_longitudinal) {
    ni <- nrow(data_num_i_list[[i]])
    length_y[i + 1] <- length_y[i] + ni
    idx <- (length_y[i] + 1):length_y[i + 1]
    row_marker_type[idx] <- biomarker_type[i]

    if (biomarker_type[i] == "continuous") {
      longname <- as.character(formula(lfit[[i]]))[2]
      y_all_vec[idx] <- unlist(data_num_i_list[[i]][, longname])
    } else {
      ### ordinal marker: only the observed category is known, so record the
      ### pair of thresholds bracketing its (unobserved) latent Gaussian
      ### score instead of a point value -- see description above.
      resp_name <- all.vars(long_sub_fixed[[i]])[1]
      y_cat <- as.numeric(data_num_i_list[[i]][[resp_name]])
      full_alpha <- c(-Inf, thresholds[[i]], Inf)
      alpha_lower_vec[idx] <- full_alpha[y_cat]
      alpha_upper_vec[idx] <- full_alpha[y_cat + 1]
    }
  }

  ### regression parameter matrix, block-diagonal per marker (mirrors
  ### build_conditional_design()'s parameter_matrix); ordinal markers use
  ### lfit[[i]]$beta, which has no intercept column -- the per-category
  ### thresholds play that role instead -- exactly matching the
  ### names(lfit[[i]]$beta) column subsetting longitudinalSubCopula() uses
  ### to build the matching design matrix at fitting time.
  coef_list <- lapply(seq_len(n_longitudinal), function(i) {
    if (biomarker_type[i] == "continuous") lfit[[i]]$coefficients$fixed else lfit[[i]]$beta
  })
  n_lfit_total <- sum(sapply(coef_list, length))
  parameter_matrix <- matrix(0, nrow = n_lfit_total, ncol = n_longitudinal)
  length_p <- rep(0, n_longitudinal + 1)
  for (i in 1:n_longitudinal) {
    length_p[i + 1] <- length_p[i] + length(coef_list[[i]])
    parameter_matrix[(length_p[i] + 1):length_p[i + 1], i] <- coef_list[[i]]
  }

  ### random-effects design (A_i) / Sigma_all -- identical in form to
  ### build_conditional_design(); see random_effects_design().
  A_i <- random_effects_design(data_num_i_list, long_sub_random, Sigma)

  Sigma_vector <- c()
  for (i in 1:n_longitudinal) {
    Sigma_vector <- c(Sigma_vector, rep(sigma.longitudinal[i]^2, dim(data_num_i_list[[i]])[1]))
  }
  ### see build_conditional_design() for why the explicit length is passed
  ### to diag() rather than relying on diag(Sigma_vector) alone.
  Sigma_all <- A_i %*% Sigma %*% t(A_i) + diag(Sigma_vector, length(Sigma_vector))

  list(y_all_vec = y_all_vec,
       alpha_lower_vec = alpha_lower_vec,
       alpha_upper_vec = alpha_upper_vec,
       row_marker_type = row_marker_type,
       parameter_matrix = parameter_matrix,
       Sigma_all = Sigma_all)
}

#' Joint density x probability for a mixed continuous/ordinal stacked vector
#'
#' @description Given one patient's full stacked mean vector and covariance
#' matrix (across all jointly-fit biomarkers' repeated measurements),
#' computes the exact Gaussian density at the observed values for the
#' continuous-marker rows, times the Gaussian-copula probability that the
#' unobserved ordinal-marker latent scores fall in the box implied by their
#' observed categories' cumulative-link thresholds -- conditional on the
#' continuous rows' observed values (via the standard multivariate-normal
#' conditioning/Schur-complement formula). This is the quantity
#' \code{conditionalYT()}/\code{conditionalYDT()} compute directly as a
#' single Gaussian density when every marker is continuous (see the
#' quadratic-form comment in \code{R/conditionalYT.R}); setting
#' \code{oo_idx = integer(0)} here reduces to exactly that computation
#' (the full \code{(2*pi)^{-n/2}} normalizing constant is kept, not
#' dropped, even though it would cancel in \code{predictRisk()}'s own
#' ratio -- \code{dynamicPredictionBio()} instead divides this function's
#' output by \code{mvtnorm::dmvnorm()}-normalized values from
#' \code{conditionalYTBioCopula()}/\code{conditionalYDTBioCopula()}, so the
#' constant must be correct in absolute terms, not just consistent within
#' one call, for \code{Y_density} to remain a genuine, correctly-scaled
#' conditional density -- see \code{R/conditionalYTBio.R}'s analogous
#' all-continuous computation, which keeps its own constant for the same
#' reason).
#' Returns the \strong{log} of that density/probability, so that it neither
#' overflows nor underflows (see \code{predictRisk()} for how it is
#' exponentiated); \code{-Inf} if the ordinal box probability is 0.
#' @param Sigma_all Full stacked covariance matrix for this patient (as
#' returned by \code{build_conditional_design_copula()}).
#' @param mu_full Full stacked mean vector, same row order as \code{Sigma_all}.
#' @param y_all_vec Full stacked observed-value vector; entries at
#' \code{oo_idx} are ignored.
#' @param alpha_lower_vec,alpha_upper_vec Full stacked threshold-bound
#' vectors; entries at \code{cc_idx} are ignored.
#' @param cc_idx,oo_idx Integer index vectors (into the stacked row order)
#' of the continuous- and ordinal-marker rows respectively.
#' @keywords internal
mixed_density_prob_copula <- function(Sigma_all, mu_full, y_all_vec, alpha_lower_vec, alpha_upper_vec,
                                       cc_idx, oo_idx) {
  mixed_density_copula_prep(Sigma_all, cc_idx, oo_idx)(mu_full, y_all_vec, alpha_lower_vec,
                                                       alpha_upper_vec)
}

#' \code{mixed_density_prob_copula()} for many means under one covariance
#'
#' @description Factors the parts of \code{mixed_density_prob_copula()}
#' that depend only on the covariance -- the inverse and log-determinant of
#' the continuous block, and the conditional covariance of the ordinal
#' block -- once, and returns a function of the mean, observed values and
#' threshold bounds that does the rest with the same arithmetic. Used to
#' evaluate a patient's density at every integration grid point.
#'
#' @inheritParams mixed_density_prob_copula
#' @return A function \code{(mu_full, y_all_vec, alpha_lower_vec,
#'   alpha_upper_vec)} returning the log density.
#' @keywords internal
mixed_density_copula_prep <- function(Sigma_all, cc_idx, oo_idx) {
  if (length(cc_idx) > 0) {
    Sigma_cc <- Sigma_all[cc_idx, cc_idx, drop = FALSE]
    Sigma_cc_solve <- solve(Sigma_cc)
    ### on the log scale: det() and exp() over/underflow with many
    ### observations or large-scale biomarkers
    log_const_c <- -length(cc_idx) / 2 * log(2 * pi) -
      0.5 * as.numeric(determinant(Sigma_cc, logarithm = TRUE)$modulus)
  }
  if (length(oo_idx) > 0) {
    Sigma_oo <- Sigma_all[oo_idx, oo_idx, drop = FALSE]
    if (length(cc_idx) > 0) {
      Sigma_oc <- Sigma_all[oo_idx, cc_idx, drop = FALSE]
      Sigma_co <- Sigma_all[cc_idx, oo_idx, drop = FALSE]
      gain_o <- Sigma_oc %*% Sigma_cc_solve
      cond_cov_o <- as.matrix(Sigma_oo - Sigma_oc %*% Sigma_cc_solve %*% Sigma_co)
    } else {
      cond_cov_o <- as.matrix(Sigma_oo)
    }
  }

  function(mu_full, y_all_vec, alpha_lower_vec, alpha_upper_vec) {
    if (length(cc_idx) > 0) {
      resid_c <- y_all_vec[cc_idx] - mu_full[cc_idx]
      quad_c <- as.numeric(t(resid_c) %*% Sigma_cc_solve %*% resid_c)
      log_dens_c <- log_const_c - 0.5 * quad_c
    } else {
      log_dens_c <- 0
    }
    if (length(oo_idx) == 0) {
      return(log_dens_c)
    }
    mu_o <- mu_full[oo_idx]
    cond_mean_o <- if (length(cc_idx) > 0) mu_o + as.numeric(gain_o %*% resid_c) else mu_o
    prob_o <- as.numeric(mvtnorm::pmvnorm(lower = alpha_lower_vec[oo_idx], upper = alpha_upper_vec[oo_idx],
                                           mean = cond_mean_o, sigma = cond_cov_o))
    log_dens_c + log(prob_o)
  }
}

#' A biomarker's mean at every grid point, continuous or ordinal
#'
#' @description For the copula density functions: \code{fixed_mean_grid()}
#' for a continuous biomarker; for an ordinal one, the linear predictor
#' \eqn{X(t)\beta} of its \code{clmm()} fit (no intercept; the thresholds
#' take its place) from \code{long_sub_fixed[[i]]}, with all grid points
#' stacked into one \code{model.matrix()} call. An ordinal formula with a
#' term whose values depend on the other rows (\code{poly()},
#' \code{splines::ns()}/\code{bs()}, \code{factor()}, \code{scale()},
#' \code{cut()}, which \code{longitudinalSub()} warns about) is evaluated
#' one grid point at a time, as before, so that its value does not change.
#'
#' @return A \code{nrow(data_i)} x \code{length(l_i)} matrix.
#' @keywords internal
copula_mean_grid <- function(i, data_i, long_fit_all, l_i, survival_variable, survival_variable_all,
                             survival_trans_function, fixed_values = list()) {
  lfit_i <- long_fit_all$lfit[[i]]
  if (long_fit_all$biomarker_type[i] == "continuous") {
    xlev_i <- if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
    return(fixed_mean_grid(data_i, lfit_i, xlev_i, l_i, survival_variable, survival_variable_all,
                           survival_trans_function, fixed_values = fixed_values))
  }
  formula_i <- long_fit_all$long_sub_fixed[[i]]
  beta <- lfit_i$beta
  n <- nrow(data_i)
  K <- length(l_i)
  for (name in names(fixed_values)) data_i[[name]] <- fixed_values[[name]]
  if (length(unsafe_formula_calls(formula_i, c("poly", "ns", "bs", "factor", "scale", "cut"))) > 0) {
    out <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
      d <- set_survival_columns(data_i, survival_variable, l_i[k], survival_variable_all,
                                survival_trans_function)
      X <- model.matrix(formula_i, d)[, names(beta), drop = FALSE]
      mu <- c(X %*% beta)
      out[, k] <- c(mu, rep(NA_real_, n - length(mu)))
    }
    return(out)
  }
  trans <- lapply(seq_along(survival_variable_all), function(surv_i)
    vapply(l_i, function(l) apply_survival_trans(survival_trans_function[[surv_i]], l, surv_i),
           numeric(1)))
  stacked <- data_i[rep(seq_len(n), times = K), , drop = FALSE]
  stacked[[survival_variable]] <- rep(l_i, each = n)
  for (surv_i in seq_along(survival_variable_all)) {
    stacked[[survival_variable_all[[surv_i]]]] <- rep(trans[[surv_i]], each = n)
  }
  mf <- model.frame(formula_i, stacked, na.action = stats::na.pass)
  X <- model.matrix(formula_i, mf)[, names(beta), drop = FALSE]
  matrix(X %*% beta, n, K)
}

#' A patient's copula log densities at every grid point
#'
#' @description Shared by \code{conditionalYTCopula()},
#' \code{conditionalYDTCopula()}, \code{conditionalYTBioCopula()} and
#' \code{conditionalYDTBioCopula()}. The patient's covariance, observed
#' values and threshold bounds do not depend on the event time, so the
#' design (\code{build_conditional_design_copula()}) is built once per
#' candidate value of the predicted biomarker instead of once per grid
#' point, and the means at all grid points come from
#' \code{copula_mean_grid()}. The densities are then evaluated in the same
#' order as before -- grid point, candidate value, event type 1 then 0 --
#' so the Monte Carlo draws of \code{mvtnorm::pmvnorm()} are the same.
#'
#' @param data_num_i_list,rep_num_i_list The patient's rows per biomarker.
#' @param sigma.longitudinal Residual SDs (1 for an ordinal biomarker).
#' @param event_type_variable \code{NULL} without competing risks;
#'   otherwise the event-type column, evaluated at 1 and at 0.
#' @param bio NULL, or a list with \code{index}, \code{name}, \code{row}
#'   and \code{values}: the predicted biomarker, its response column, the
#'   row of the predicted measurement and its candidate values.
#' @return A list with one element per event type (1, then 0; just one
#'   without competing risks), each a \code{length(l_i)} x (number of
#'   candidate values, 1 without \code{bio}) matrix of log densities.
#' @keywords internal
copula_logdens_grid <- function(data_num_i_list, rep_num_i_list, long_fit_all, l_i,
                                survival_variable, time_variable, survival_variable_all,
                                survival_trans_function, sigma.longitudinal,
                                event_type_variable = NULL, bio = NULL) {
  lfit <- long_fit_all$lfit
  biomarker_type <- long_fit_all$biomarker_type
  long_sub_fixed <- long_fit_all$long_sub_fixed
  n_longitudinal <- length(lfit)

  for (i in seq_len(n_longitudinal)) {
    assert_uses_event_time(formula(long_sub_fixed[[i]]), i, survival_variable, survival_variable_all)
    data_num_i_list[[i]] <- set_survival_columns(data_num_i_list[[i]], survival_variable, l_i[1],
                                                 survival_variable_all, survival_trans_function)
    if (biomarker_type[i] == "continuous") {
      ### NA in nlme outcome (longitudinal biomarkers), replace with 999
      ### -- same convention as conditionalYT().
      resp_i <- as.character(formula(long_sub_fixed[[i]])[[2]])
      data_num_i_list[[i]][[resp_i]][is.na(data_num_i_list[[i]][[resp_i]])] <- 999
    }
  }

  event_types <- if (is.null(event_type_variable)) list(list()) else
    list(stats::setNames(list(1), event_type_variable), stats::setNames(list(0), event_type_variable))
  means <- lapply(event_types, function(fixed_values) {
    do.call(rbind, lapply(seq_len(n_longitudinal), function(i)
      copula_mean_grid(i, data_num_i_list[[i]], long_fit_all, l_i, survival_variable,
                       survival_variable_all, survival_trans_function, fixed_values)))
  })

  values <- if (is.null(bio)) list(NULL) else as.list(bio$values)
  designs <- lapply(values, function(value) {
    d <- data_num_i_list
    if (!is.null(bio)) d[[bio$index]][[bio$name]][bio$row] <- value
    build_conditional_design_copula(rep_num_i_list, d, lfit, long_fit_all$Sigma_fit,
                                    sigma.longitudinal, time_variable, n_longitudinal,
                                    biomarker_type, long_sub_fixed, long_fit_all$thresholds,
                                    long_fit_all$long_sub_random)
  })
  cc_idx <- which(designs[[1]]$row_marker_type == "continuous")
  oo_idx <- which(designs[[1]]$row_marker_type == "ordinal")
  dens <- mixed_density_copula_prep(designs[[1]]$Sigma_all, cc_idx, oo_idx)

  out <- lapply(event_types, function(e) matrix(NA_real_, length(l_i), length(values)))
  for (it in seq_along(l_i)) {
    for (v in seq_along(values)) {
      for (e in seq_along(event_types)) {
        out[[e]][it, v] <- dens(means[[e]][, it], designs[[v]]$y_all_vec,
                                designs[[v]]$alpha_lower_vec, designs[[v]]$alpha_upper_vec)
      }
    }
  }
  out
}
