#' Filter longitudinal data down to a single patient
#'
#' @description Shared helper for \code{conditionalYT} and \code{conditionalYDT}:
#' extracts each biomarker's rows for one patient ID, replicating the shared
#' data frame across biomarkers when only one was supplied.
#'
#' @return A list with \code{rep_num_i_list} and \code{data_num_i_list}, or
#' \code{NULL} if any biomarker has zero rows for this patient (caller should
#' skip the patient in that case).
#' @keywords internal
select_patient_longitudinal_data <- function(data.long, num, num_i, n_longitudinal, time_variable) {
  rep_num_i_list <- list()
  data_num_i_list <- list()

  for (i in 1:n_longitudinal) {
    df <- data.long[[i]]
    selected_data <- df[which(as.character(df[[num]]) == as.character(num_i)), , drop = FALSE]
    rep_num_i_list[[i]] <- rep(1, length(unlist(selected_data[time_variable])))
    data_num_i_list[[i]] <- selected_data
  }
  if (length(data.long) == 1) {
    for (i in 1:n_longitudinal) {
      rep_num_i_list[[i]] <- rep_num_i_list[[1]]
      data_num_i_list[[i]] <- data_num_i_list[[1]]
    }
  }

  if (any(sapply(data_num_i_list, nrow) == 0)) {
    return(NULL)
  }

  list(rep_num_i_list = rep_num_i_list, data_num_i_list = data_num_i_list)
}

#' Random-effects design matrix for one patient, across all biomarkers
#'
#' @description Shared helper for \code{build_conditional_design()},
#' \code{build_conditional_design_copula()}, and \code{process_variance()}:
#' builds each biomarker's random-effects design from the left-hand side of
#' its own \code{long_sub_random} formula (e.g. \code{~ 1}, \code{~ year},
#' \code{~ year + I(year^2)}) and stacks them block-diagonally, matching the
#' block-diagonal \code{Sigma_fit} built by \code{longitudinalSub()}.
#' Previously this assumed every biomarker had either a random intercept
#' only, or a random intercept plus a slope on \code{time_variable}, so a
#' mix of the two (or any other random-effects formula) failed with
#' "non-conformable arguments".
#'
#' @param data_num_i_list One patient's data, one data frame per biomarker.
#' @param long_sub_random The list of random-effects formulas.
#' @param Sigma The random-effects covariance matrix (\code{Sigma_fit}).
#' @return A matrix with one row per observation (stacked across
#' biomarkers) and one column per random effect.
#' @keywords internal
random_effects_design <- function(data_num_i_list, long_sub_random, Sigma) {
  blocks <- lapply(seq_along(data_num_i_list), function(i) {
    ffk <- nlme::splitFormula(long_sub_random[[i]], "|")[[1]]
    mf <- model.frame(ffk, data_num_i_list[[i]], na.action = stats::na.pass)
    model.matrix(ffk, mf)
  })
  A_i <- as.matrix(Matrix::bdiag(blocks))
  if (ncol(A_i) != nrow(Sigma)) {
    stop(sprintf(
      "The random-effects formulas in long_fit_all$long_sub_random give %d random effect(s), but long_fit_all$Sigma_fit is %d x %d.",
      ncol(A_i), nrow(Sigma), ncol(Sigma)), call. = FALSE)
  }
  A_i
}

#' Build the per-patient longitudinal design matrices
#'
#' @description Shared helper for \code{conditionalYT} and \code{conditionalYDT}:
#' constructs the longitudinal outcome matrix, fixed-effect parameter matrix,
#' and random-effect variance-covariance pieces used by the conditional
#' density quadratic form, for a single patient.
#'
#' @return A list with \code{longitudinal_all_matrix}, \code{parameter_matrix},
#' and \code{cov_fac}, the patient's factored covariance from
#' \code{cov_factor()}.
#' @keywords internal
build_conditional_design <- function(rep_num_i_list, data_num_i_list, lfit, Sigma,
                                      sigma.longitudinal, time_variable, n_longitudinal,
                                      long_sub_random) {
  ####Initialize longitudinal matrix for all biomarkers
  longitudinal_all_matrix <- matrix(0, nrow = sum(sapply(data_num_i_list, nrow)), ncol = n_longitudinal)
  ####Constructing the longitudinal matrix for all biomarkers
  length_y = rep(0, n_longitudinal + 1)
  for (i in 1:n_longitudinal) {
    longname = as.character(formula(lfit[[i]]))[2]
    length_y[i + 1] = length_y[i] + length(c(unlist(data_num_i_list[[i]][, longname])))
    longitudinal_all_matrix[c((length_y[i] + 1):length_y[i + 1]), i] <-
      unlist(data_num_i_list[[i]][, longname])
  }

  ####constructing the regression parameters' matrix
  n_lfit_total = 0
  for (i in seq_len(length(lfit))) {
    n_lfit_total = n_lfit_total + length(lfit[[i]]$coefficients$fixed)
  }
  ####Initialize the regression parameters' matrix
  parameter_matrix <- matrix(0, nrow = n_lfit_total, ncol = length(lfit))
  length_p = rep(0, n_longitudinal + 1)
  for (i in 1:n_longitudinal) {
    length_p[i + 1] = length_p[i] + length(lfit[[i]]$coefficients$fixed)
    parameter_matrix[c((length_p[i] + 1):length_p[i + 1]), i] <- lfit[[i]]$coefficients$fixed
  }

  A_i <- random_effects_design(data_num_i_list, long_sub_random, Sigma)

  Sigma_vector = c()
  for (i in 1:n_longitudinal) {
    Sigma_vector = c(Sigma_vector, rep(sigma.longitudinal[i]^2, dim(data_num_i_list[[i]])[1]))
  }
  ### Sigma_all = A_i %*% Sigma %*% t(A_i) + diag(Sigma_vector), factored
  ### once per patient without necessarily forming it (see cov_factor()).
  ### The log-determinant is kept on the log scale: det() overflows to Inf
  ### (and the density to 0) with many observations or large-scale markers.
  cov_fac <- cov_factor(A_i, Sigma, Sigma_vector)

  list(longitudinal_all_matrix = longitudinal_all_matrix,
       parameter_matrix = parameter_matrix,
       cov_fac = cov_fac)
}

#' Filter longitudinal data down to a single patient, substituting the
#' candidate biomarker prediction value
#'
#' @description Shared helper for \code{conditionalYTBio} and
#' \code{conditionalYDTBio}: extracts each biomarker's rows for one patient
#' ID, and for the biomarker being predicted (\code{bio_i}), appends a row at
#' \code{time_new} with each candidate value in \code{Y_all} substituted in
#' turn.
#'
#' @return A list with \code{rep_num_i_list}, \code{data_num_i_list}, and
#' \code{Y_select_all} (a matrix of the substituted biomarker values, one
#' column per element of \code{Y_all}), or \code{NULL} if any biomarker has
#' zero rows for this patient (caller should skip the patient in that case).
#' @keywords internal
select_patient_longitudinal_data_bio <- function(data.long, num, num_i, n_longitudinal, time_variable,
                                                  bio_i, time_new, Y_all, long_fit_all) {
  rep_num_i_list <- list()
  data_num_i_list <- list()
  Y_select_all <- NULL

  for (i in 1:n_longitudinal) {
    df <- data.long[[i]]
    selected_data <- df[which(as.character(df[[num]]) == as.character(num_i)), , drop = FALSE]

    if (i == bio_i) {
      selected_data <- rbind(selected_data, selected_data[nrow(selected_data), ])
      selected_data[time_variable][nrow(selected_data), ] = time_new
      bio_i_name = as.character(formula(long_fit_all$long_sub_fixed[[i]])[[2]])
      Y_select_all = c()
      for (Y_new in Y_all) {
        selected_data[bio_i_name][nrow(selected_data), ] = Y_new
        Y_select_all = cbind(Y_select_all, unlist(selected_data[bio_i_name]))
      }
    }
    rep_num_i_list[[i]] <- rep(1, length(unlist(selected_data[time_variable])))
    data_num_i_list[[i]] <- selected_data
  }
  if (length(data.long) == 1) {
    for (i in 1:n_longitudinal) {
      rep_num_i_list[[i]] <- rep_num_i_list[[1]]
      data_num_i_list[[i]] <- data_num_i_list[[1]]
    }
  }

  if (any(sapply(data_num_i_list, nrow) == 0)) {
    return(NULL)
  }

  list(rep_num_i_list = rep_num_i_list, data_num_i_list = data_num_i_list, Y_select_all = Y_select_all)
}

#' Build the longitudinal outcome matrix for all candidate biomarker values
#'
#' @description Shared helper for \code{conditionalYTBio} and
#' \code{conditionalYDTBio}: for each candidate value in \code{Y_all},
#' assembles the row of observed longitudinal outcomes across biomarkers,
#' substituting the candidate value for the biomarker being predicted.
#'
#' @return A matrix with \code{length(Y_all)} rows, one per candidate value.
#' @keywords internal
build_longitudinal_matrix_bio <- function(data_num_i_list, lfit, bio_i, Y_select_all, n_longitudinal, Y_all) {
  longitudinal_all_matrix <- c()
  for (Y_i in 1:length(Y_all)) {
    longitudinal_all_matrix_tran <- c()
    for (i in 1:n_longitudinal) {
      longname = as.character(formula(lfit[[i]]))[2]
      if (i == bio_i) {
        longitudinal_all_matrix_tran = c(longitudinal_all_matrix_tran, Y_select_all[, Y_i])
      } else {
        longitudinal_all_matrix_tran = c(longitudinal_all_matrix_tran, unlist(data_num_i_list[[i]][, longname]))
      }
    }
    longitudinal_all_matrix = rbind(longitudinal_all_matrix, longitudinal_all_matrix_tran)
  }
  longitudinal_all_matrix
}

#' Set a patient's survival-time columns to one integration grid point
#'
#' @description Sets the survival-time column and each
#' \code{survival_variable_all} column (through its
#' \code{survival_trans_function}) to their values at event time
#' \code{l}.
#'
#' @return \code{df}, updated.
#' @keywords internal
set_survival_columns <- function(df, survival_variable, l, survival_variable_all,
                                 survival_trans_function) {
  df[[survival_variable]] <- l
  for (surv_i in seq_along(survival_variable_all)) {
    df[[survival_variable_all[[surv_i]]]] <- apply_survival_trans(survival_trans_function[[surv_i]], l, surv_i)
  }
  df
}


#' A biomarker's fixed-effects mean at every integration grid point
#'
#' @description The longitudinal mean \eqn{X(t)\beta} depends on the event
#' time \eqn{t} only through the survival-time column and its
#' \code{survival_variable_all} transformations. This evaluates it for all
#' rows of \code{data_i} (all patients) at all grid points \code{l_i} with
#' one \code{model.frame()}/\code{model.matrix()} call on the rows stacked
#' once per grid point, instead of one call per patient and grid point,
#' whose fixed overhead dominated the run time of the prediction functions.
#' The terms, factor levels and contrasts are those of the fit, as in
#' \code{predict.lme()}, so the design matrix is the one the model was fit
#' with, including data-dependent bases such as \code{poly()} and
#' \code{splines::ns()}. Transformed terms such as \code{log(years)} are
#' evaluated at each grid point. Rows are processed in chunks of at most
#' \code{max_rows} stacked rows.
#'
#' @param data_i Prediction rows for the biomarker.
#' @param lfit_i The biomarker's fitted model (\code{terms},
#'   \code{contrasts}, \code{coefficients$fixed}).
#' @param xlev_i Its factor levels, or \code{NULL}.
#' @param l_i Grid of event times.
#' @param fixed_values Named list of further columns to set (e.g. the event
#'   type).
#' @param beta Fixed effects; by default those of \code{lfit_i}.
#' @param max_rows Largest number of stacked rows per chunk.
#' @return A \code{nrow(data_i)} x \code{length(l_i)} matrix; a row with a
#'   missing covariate is \code{NA}.
#' @keywords internal
fixed_mean_grid <- function(data_i, lfit_i, xlev_i, l_i, survival_variable,
                            survival_variable_all, survival_trans_function,
                            fixed_values = list(), beta = lfit_i$coefficients$fixed,
                            max_rows = 2e5) {
  n <- nrow(data_i)
  K <- length(l_i)
  out <- matrix(NA_real_, n, K)
  if (n == 0 || K == 0) return(out)
  trans <- lapply(seq_along(survival_variable_all), function(surv_i)
    vapply(l_i, function(l) apply_survival_trans(survival_trans_function[[surv_i]], l, surv_i),
           numeric(1)))
  for (name in names(fixed_values)) data_i[[name]] <- fixed_values[[name]]
  per_chunk <- max(1, floor(max_rows / n))
  for (start in seq(1, K, by = per_chunk)) {
    ks <- start:min(K, start + per_chunk - 1)
    stacked <- data_i[rep(seq_len(n), times = length(ks)), , drop = FALSE]
    stacked[[survival_variable]] <- rep(l_i[ks], each = n)
    for (surv_i in seq_along(survival_variable_all)) {
      stacked[[survival_variable_all[[surv_i]]]] <- rep(trans[[surv_i]][ks], each = n)
    }
    mf <- model.frame(lfit_i$terms, stacked, xlev = xlev_i, na.action = stats::na.pass)
    X <- model.matrix(lfit_i$terms, mf, contrasts.arg = lfit_i$contrasts)
    out[, ks] <- matrix(X %*% beta, n, length(ks))
  }
  out
}

#' Check that a longitudinal sub-model uses the event time
#'
#' @param formula_i The biomarker's fixed-effects formula.
#' @return Invisibly \code{TRUE}; an error naming biomarker \code{i}
#'   otherwise.
#' @keywords internal
assert_uses_event_time <- function(formula_i, i, survival_variable, survival_variable_all) {
  if (!any(c(survival_variable, unlist(survival_variable_all)) %in% all.vars(formula_i))) {
    stop(sprintf("The longitudinal sub-model for biomarker %d uses neither the survival-time variable `%s` nor any of its transformations (survival_variable_all); the backward joint model needs the event time as a covariate, directly or transformed, in every long_sub_fixed formula.", i, survival_variable), call. = FALSE)
  }
  invisible(TRUE)
}

#' One patient's stacked fixed-effects means at every grid point
#'
#' @description \code{fixed_mean_grid()} for each biomarker's rows of one
#' patient (\code{data_num_i_list}), stacked in biomarker order as in
#' \code{build_conditional_design()}.
#'
#' @return A matrix, one row per stacked observation and one column per
#'   grid point in \code{l_i}.
#' @keywords internal
patient_mean_grid <- function(data_num_i_list, long_fit_all, l_i, survival_variable,
                              survival_variable_all, survival_trans_function,
                              fixed_values = list()) {
  lfit <- long_fit_all$lfit
  do.call(rbind, lapply(seq_along(lfit), function(i) {
    xlev_i <- if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
    fixed_mean_grid(data_num_i_list[[i]], lfit[[i]], xlev_i, l_i, survival_variable,
                    survival_variable_all, survival_trans_function, fixed_values = fixed_values)
  }))
}

#' Log densities of every candidate outcome vector at every grid point
#'
#' @description For \code{conditionalYTBio()}/\code{conditionalYDTBio()}:
#' the log density of each row of \code{x} (one per candidate value of the
#' predicted biomarker) with mean \code{mu[, k]}, for every grid point
#' \code{k}, in one \code{cov_logdens()} call.
#'
#' @param fac The patient's factored covariance (\code{cov_factor()}).
#' @param x Matrix, one row per candidate value, one column per stacked
#'   observation.
#' @param mu Matrix of means, one row per stacked observation, one column
#'   per grid point.
#' @return A \code{nrow(x)} x \code{ncol(mu)} matrix.
#' @keywords internal
cov_logdens_grid <- function(fac, x, mu) {
  if (is.vector(x)) x <- matrix(x, ncol = length(x))
  if (ncol(x) != fac$N || nrow(mu) != fac$N) stop("x, mean and the covariance have non-conforming size")
  n_y <- nrow(x)
  K <- ncol(mu)
  E <- t(x)[, rep(seq_len(n_y), times = K), drop = FALSE] - mu[, rep(seq_len(K), each = n_y), drop = FALSE]
  matrix(cov_logdens(fac, E), n_y, K)
}
