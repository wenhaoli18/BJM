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
#' \code{impute_latent_ordinal()} in \code{R/copulaLongitudinal.R}: latent
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
                                             biomarker_type, long_sub_fixed, thresholds) {
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
  ### build_conditional_design(); see that function for the random
  ### intercept-only vs. intercept+slope branching logic.
  A_i_ <- list()
  if (dim(Sigma)[1] == n_longitudinal) {
    for (i in 1:n_longitudinal) A_i_[[i]] <- rbind(rep_num_i_list[[i]])
  } else {
    for (i in 1:n_longitudinal) A_i_[[i]] <- rbind(rep_num_i_list[[i]], unlist(data_num_i_list[[i]][time_variable]))
  }

  A_i <- matrix(0, nrow = sum(sapply(A_i_, ncol)), ncol = sum(sapply(A_i_, nrow)))
  length_A <- rep(0, n_longitudinal + 1)
  for (i in 1:n_longitudinal) {
    length_A[i + 1] <- length_A[i] + dim(A_i_[[i]])[2]
    if (dim(Sigma)[1] == n_longitudinal) {
      A_i[(length_A[i] + 1):length_A[i + 1], i] <- t(A_i_[[i]])
    } else {
      A_i[(length_A[i] + 1):length_A[i + 1], (2 * i - 1):(2 * i)] <- t(A_i_[[i]])
    }
  }

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

#' Build one biomarker's transposed design matrix at prediction time
#' (Gaussian-copula-aware)
#'
#' @description Copula-aware counterpart to the \code{model.matrix()} call
#' inline in \code{conditionalYT()}/\code{conditionalYDT()}. Continuous
#' markers use the cached \code{terms}/\code{xlevels}/\code{contrasts} from
#' the \code{nlme::lme()} fit, exactly as those functions do. Ordinal
#' markers instead reuse \code{longitudinalSubCopula()}'s own construction
#' verbatim (\code{model.matrix(long_sub_fixed[[i]], data)} subset to
#' \code{names(lfit[[i]]$beta)}), since an \code{ordinal::clmm} fit's
#' \code{$terms} does not carry the no-intercept convention
#' \code{lfit[[i]]$beta} expects.
#' @return A matrix with one row per fixed-effect coefficient and one column
#' per observation (row of \code{data_i}).
#' @keywords internal
build_LME_indi_matrix_copula <- function(i, data_i, lfit, long_fit_all, long_sub_fixed) {
  if (long_fit_all$biomarker_type[i] == "continuous") {
    terms_model <- lfit[[i]]$terms
    xlev_i <- if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
    mf_i <- model.frame(terms_model, data_i, xlev = xlev_i)
    out <- t(model.matrix(terms_model, mf_i, contrasts.arg = lfit[[i]]$contrasts))
  } else {
    beta_names_i <- names(lfit[[i]]$beta)
    full_mm <- model.matrix(long_sub_fixed[[i]], data_i)
    out <- t(full_mm[, beta_names_i, drop = FALSE])
  }
  if (dim(out)[2] != dim(data_i)[1]) {
    out <- cbind(out, matrix(NA, dim(out)[1], dim(data_i)[1] - dim(out)[2]))
  }
  out
}

#' Combine per-marker transposed design matrices into one block-diagonal
#' matrix
#'
#' @description Shared helper: replicates the \code{cbind}/\code{rbind}
#' block-diagonal assembly duplicated inline in \code{conditionalYT.R}/
#' \code{conditionalYDT.R}, so the copula prediction functions can reuse it
#' rather than re-copy it a third time.
#' @keywords internal
build_block_diagonal_matrix <- function(matrix_list) {
  rows <- list()
  for (i in seq_along(matrix_list)) {
    columns <- list()
    for (j in seq_along(matrix_list)) {
      if (i == j) {
        columns[[j]] <- matrix_list[[i]]
      } else {
        columns[[j]] <- matrix(0, nrow = nrow(matrix_list[[i]]), ncol = ncol(matrix_list[[j]]))
      }
    }
    rows[[i]] <- do.call(cbind, columns)
  }
  do.call(rbind, rows)
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
  if (length(cc_idx) > 0) {
    Sigma_cc <- Sigma_all[cc_idx, cc_idx, drop = FALSE]
    mu_c <- mu_full[cc_idx]
    y_c <- y_all_vec[cc_idx]
    resid_c <- y_c - mu_c
    Sigma_cc_solve <- solve(Sigma_cc)
    quad_c <- as.numeric(t(resid_c) %*% Sigma_cc_solve %*% resid_c)
    dens_c <- (2 * pi)^(-length(cc_idx) / 2) * det(Sigma_cc)^(-0.5) * exp(-0.5 * quad_c)
  } else {
    dens_c <- 1
  }

  if (length(oo_idx) == 0) {
    return(dens_c)
  }

  mu_o <- mu_full[oo_idx]
  Sigma_oo <- Sigma_all[oo_idx, oo_idx, drop = FALSE]
  if (length(cc_idx) > 0) {
    Sigma_oc <- Sigma_all[oo_idx, cc_idx, drop = FALSE]
    Sigma_co <- Sigma_all[cc_idx, oo_idx, drop = FALSE]
    cond_mean_o <- mu_o + as.numeric(Sigma_oc %*% Sigma_cc_solve %*% resid_c)
    cond_cov_o <- Sigma_oo - Sigma_oc %*% Sigma_cc_solve %*% Sigma_co
  } else {
    cond_mean_o <- mu_o
    cond_cov_o <- Sigma_oo
  }

  lower <- alpha_lower_vec[oo_idx]
  upper <- alpha_upper_vec[oo_idx]
  prob_o <- as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper,
                                         mean = cond_mean_o, sigma = as.matrix(cond_cov_o)))
  dens_c * prob_o
}
