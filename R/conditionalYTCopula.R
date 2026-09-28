#' Conditional distribution of Y|T for a mixed continuous/ordinal
#' (Gaussian-copula) joint model, if no competing risk
#'
#' @description Copula-aware counterpart to \code{conditionalYT()}, used by
#' \code{predictRisk()} whenever \code{long_fit_all$biomarker_type}
#' contains at least one \code{"ordinal"} biomarker (see
#' \code{longitudinalSubCopula()}). Continuous biomarkers contribute an exact
#' Gaussian density factor, exactly as \code{conditionalYT()} computes for
#' every biomarker; ordinal biomarkers instead contribute a Gaussian-copula
#' box probability (via \code{mixed_density_prob_copula()}), since only the
#' cumulative-link category -- not the exact underlying latent score -- is
#' observed for them. All-continuous fits are unaffected: they are still
#' routed to \code{conditionalYT()} by \code{predictRisk()}, this
#' function is never called for them, and \code{mixed_density_prob_copula()}
#' reduces to (a constant multiple of) \code{conditionalYT()}'s own
#' computation when there are no ordinal markers -- see that function's
#' documentation.
#'
#' @inheritParams conditionalYT
#' @keywords internal
conditionalYTCopula <- function(data_predict_all, long_fit_all, l_i, survival_variable,
                                 time_variable, survival_variable_all, survival_trans_function) {
  lfit <- long_fit_all$lfit
  Sigma <- long_fit_all$Sigma_fit
  biomarker_type <- long_fit_all$biomarker_type
  thresholds <- long_fit_all$thresholds
  long_sub_fixed <- long_fit_all$long_sub_fixed
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]

  n_longitudinal <- length(lfit)
  sigma.longitudinal <- vapply(seq_len(n_longitudinal), function(i) {
    if (biomarker_type[i] == "continuous") lfit[[i]]$sigma else 1
  }, numeric(1))

  data.long <- data_predict_all
  if (!is.list(data.long) || is.data.frame(data.long)) {
    data.long <- list(data.long)
    data.long <- rep(data.long, each = n_longitudinal)
  }

  f_Y_T_D_w1 <- matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))

  iii <- 0
  for (num_i in unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])) {
    iii <- iii + 1

    patient_data <- select_patient_longitudinal_data(data.long, num, num_i, n_longitudinal, time_variable)
    if (is.null(patient_data)) next
    rep_num_i_list <- patient_data$rep_num_i_list
    data_num_i_list <- patient_data$data_num_i_list

    for (it in seq_along(l_i)) {
      data_it <- data_num_i_list
      for (i in 1:n_longitudinal) {
        data_it[[i]][survival_variable] <- l_i[it]
        if (length(survival_variable_all) != 0) {
          for (surv_i in seq_along(survival_variable_all)) {
            data_it[[i]][survival_variable_all[[surv_i]]] <-
              apply_survival_trans(survival_trans_function[[surv_i]], l_i[it], surv_i)
          }
        }

        all_variables <- all.vars(formula(long_sub_fixed[[i]]))
        if (!(survival_variable %in% all_variables)) {
          stop("Error: Condition is false. Please add survival variable to linear mixed model.")
        }

        if (biomarker_type[i] == "continuous") {
          ### NA in nlme outcome (longitudinal biomarkers), replace with 999
          ### -- same convention as conditionalYT().
          resp_i <- as.character(formula(long_sub_fixed[[i]])[[2]])
          data_it[[i]][[resp_i]][is.na(data_it[[i]][[resp_i]])] <- 999
        }
      }

      design <- build_conditional_design_copula(rep_num_i_list, data_it, lfit, Sigma,
                                                  sigma.longitudinal, time_variable, n_longitudinal,
                                                  biomarker_type, long_sub_fixed, thresholds)

      LME_indi_matrix <- lapply(seq_len(n_longitudinal), function(i) {
        build_LME_indi_matrix_copula(i, data_it[[i]], lfit, long_fit_all, long_sub_fixed)
      })
      LME_all_matrix <- build_block_diagonal_matrix(LME_indi_matrix)

      mean_all_matrix <- t(LME_all_matrix) %*% design$parameter_matrix
      mu_full <- rowSums(mean_all_matrix)

      cc_idx <- which(design$row_marker_type == "continuous")
      oo_idx <- which(design$row_marker_type == "ordinal")

      f_Y_T_D_w1[it, iii] <- mixed_density_prob_copula(design$Sigma_all, mu_full, design$y_all_vec,
                                                          design$alpha_lower_vec, design$alpha_upper_vec,
                                                          cc_idx, oo_idx)
    }
  }
  return(f_Y_T_D = list(f_Y_T_D_w1))
}
