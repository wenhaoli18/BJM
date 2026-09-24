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
    selected_data <- df[df[num] == num_i, ]
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

#' Build the per-patient longitudinal design matrices
#'
#' @description Shared helper for \code{conditionalYT} and \code{conditionalYDT}:
#' constructs the longitudinal outcome matrix, fixed-effect parameter matrix,
#' and random-effect variance-covariance pieces used by the conditional
#' density quadratic form, for a single patient.
#'
#' @return A list with \code{longitudinal_all_matrix}, \code{parameter_matrix},
#' \code{Sigma_all}, \code{det_Var_cov_estep}, \code{Sigma_all_solve}, and
#' \code{long_sigma_long}.
#' @keywords internal
build_conditional_design <- function(rep_num_i_list, data_num_i_list, lfit, Sigma,
                                      sigma.longitudinal, time_variable, n_longitudinal) {
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

  A_i_ = list()
  ###random intercept or slope, depend on variance-covariance matrix
  if (dim(Sigma)[1] == n_longitudinal) {
    ###random intercept
    for (i in 1:n_longitudinal) {
      A_i_[[i]] = rbind(rep_num_i_list[[i]])
    }
  } else {
    ###random slope
    for (i in 1:n_longitudinal) {
      A_i_[[i]] = rbind(rep_num_i_list[[i]], unlist(data_num_i_list[[i]][time_variable]))
    }
  }

  A_i <- matrix(0, nrow = sum(sapply(A_i_, ncol)), ncol = sum(sapply(A_i_, nrow)))
  length_A = rep(0, n_longitudinal + 1)
  for (i in 1:n_longitudinal) {
    length_A[i + 1] = length_A[i] + dim(A_i_[[i]])[2]
    if (dim(Sigma)[1] == n_longitudinal) {
      ###random intercept
      A_i[c((length_A[i] + 1):length_A[i + 1]), i] <- t(A_i_[[i]])
    } else {
      ###random slope
      A_i[c((length_A[i] + 1):length_A[i + 1]), (2 * i - 1):(2 * i)] <- t(A_i_[[i]])
    }
  }

  Sigma_vector = c()
  for (i in 1:n_longitudinal) {
    Sigma_vector = c(Sigma_vector, rep(sigma.longitudinal[i]^2, dim(data_num_i_list[[i]])[1]))
  }
  Sigma_all = A_i %*% Sigma %*% t(A_i) + diag(Sigma_vector)

  det_Var_cov_estep = det(2 * pi * Sigma_all)
  Sigma_all_solve = solve(Sigma_all)

  ### t(Y) %*% Sigma %*% Y
  long_sigma_long = t(longitudinal_all_matrix) %*% Sigma_all_solve %*% longitudinal_all_matrix

  list(longitudinal_all_matrix = longitudinal_all_matrix,
       parameter_matrix = parameter_matrix,
       Sigma_all = Sigma_all,
       det_Var_cov_estep = det_Var_cov_estep,
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
    selected_data <- df[df[num] == num_i, ]

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
