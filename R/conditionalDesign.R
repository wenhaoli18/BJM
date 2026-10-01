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
#' \code{Sigma_all}, \code{log_det_Var_cov_estep}, \code{Sigma_all_solve}, and
#' \code{long_sigma_long}.
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
  ### diag(Sigma_vector) alone is unsafe when Sigma_vector has length 1
  ### (e.g. a patient with a single longitudinal observation to condition
  ### on): diag() then treats that single number as a matrix *dimension*
  ### and returns an n x n identity matrix instead of the intended 1 x 1
  ### diagonal matrix, corrupting Sigma_all's dimensions downstream (a
  ### classic base R diag() gotcha -- see ?diag). Passing the length
  ### explicitly avoids the ambiguity for every length, including 1.
  Sigma_all = A_i %*% Sigma %*% t(A_i) + diag(Sigma_vector, length(Sigma_vector))

  ### log-determinant: det() itself overflows to Inf (and the density to 0)
  ### once there are many observations or biomarkers on a large scale
  log_det_Var_cov_estep = as.numeric(determinant(2 * pi * Sigma_all, logarithm = TRUE)$modulus)
  Sigma_all_solve = solve(Sigma_all)

  ### t(Y) %*% Sigma %*% Y
  long_sigma_long = t(longitudinal_all_matrix) %*% Sigma_all_solve %*% longitudinal_all_matrix

  list(longitudinal_all_matrix = longitudinal_all_matrix,
       parameter_matrix = parameter_matrix,
       Sigma_all = Sigma_all,
       log_det_Var_cov_estep = log_det_Var_cov_estep,
       Sigma_all_solve = Sigma_all_solve,
       long_sigma_long = long_sigma_long)
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

#' Can a model frame's survival-time columns be overwritten in place?
#'
#' @description \code{conditionalYDT()} and \code{conditionalYDTBio()}
#' build each biomarker's model frame once per patient and, for every
#' integration grid point, overwrite its survival-time columns in place
#' instead of rebuilding it. That is only valid when the survival time (and
#' every \code{survival_variable_all} column) enters the model frame as a
#' bare column: a model-frame column such as \code{log(years)},
#' \code{I(years^2)} or \code{poly(years, 2)} is named by its expression, so
#' it was never found and kept its value at the first grid point for the
#' whole integral. This reports whether every model-frame variable that
#' refers to one of \code{time_vars} is one of those bare columns.
#'
#' @param terms_model The fitted model's \code{terms}.
#' @param time_vars Names of the survival-time column and its user-supplied
#'   transformations.
#' @return \code{TRUE} if in-place updating is valid.
#' @keywords internal
time_columns_bare <- function(terms_model, time_vars) {
  variables <- as.list(attr(terms_model, "variables"))[-1]
  response <- attr(terms_model, "response")
  if (response > 0) variables <- variables[-response]
  all(vapply(variables, function(v) {
    !any(all.vars(v) %in% time_vars) || (is.name(v) && as.character(v) %in% time_vars)
  }, logical(1)))
}

#' Set a patient's survival-time columns to one integration grid point
#'
#' @description Sets the survival-time column and each
#' \code{survival_variable_all} column (through its
#' \code{survival_trans_function}) to their values at event time
#' \code{l}. With \code{only_existing = TRUE} (for a model frame updated in
#' place, see \code{time_columns_bare()}) columns not already present are
#' left out.
#'
#' @return \code{df}, updated.
#' @keywords internal
set_survival_columns <- function(df, survival_variable, l, survival_variable_all,
                                 survival_trans_function, only_existing = FALSE) {
  if (!only_existing || survival_variable %in% names(df)) df[[survival_variable]] <- l
  for (surv_i in seq_along(survival_variable_all)) {
    svar <- survival_variable_all[[surv_i]]
    if (!only_existing || svar %in% names(df)) {
      df[[svar]] <- apply_survival_trans(survival_trans_function[[surv_i]], l, surv_i)
    }
  }
  df
}

#' A biomarker's model frame at one integration grid point
#'
#' @description Overwrites the survival-time columns of the per-patient
#' template \code{mf} in place when \code{in_place} (see
#' \code{time_columns_bare()}); otherwise rebuilds the model frame from
#' \code{data} with those columns set, so that transformed terms such as
#' \code{log(years)} are re-evaluated at \code{l}.
#'
#' @return The model frame.
#' @keywords internal
survival_model_frame_at <- function(mf, data, in_place, terms_model, xlev, survival_variable, l,
                                    survival_variable_all, survival_trans_function) {
  if (in_place) {
    return(set_survival_columns(mf, survival_variable, l, survival_variable_all,
                                survival_trans_function, only_existing = TRUE))
  }
  data <- set_survival_columns(data, survival_variable, l, survival_variable_all,
                               survival_trans_function)
  model.frame(terms_model, data, xlev = xlev)
}
