#' Conditional distribution of Y|D,T for a mixed continuous/ordinal
#' (Gaussian-copula) joint model, if with competing risk -- biomarker-value
#' prediction
#'
#' @description Copula-aware counterpart to \code{conditionalYDTBio()}, used
#' by \code{dynamicPredictionBio()} whenever \code{long_fit_all$biomarker_type}
#' contains at least one \code{"ordinal"} biomarker and a competing-risk
#' \code{survival_fit_all} is supplied. Both a continuous and an
#' \strong{ordinal} \code{bio_i} target are supported: for an ordinal
#' \code{bio_i}, \code{Y_all} is expected to already be the vector of
#' candidate \emph{category labels} (the fitted factor's \code{levels()}, in
#' threshold order) rather than a numeric grid -- see
#' \code{compute_bio_marker_step()}, which builds that vector and translates
#' the result back into integer category codes for its caller. See
#' \code{conditionalYTBioCopula()} for why no branching is actually needed
#' here -- \code{build_conditional_design_copula()} already looks up each
#' row of biomarker \code{bio_i} generically via
#' \code{long_fit_all$biomarker_type[bio_i]}, whether that row is historical
#' or the candidate row assigned below -- and for why this function
#' evaluates \code{mvtnorm::pmvnorm()} once per patient/\code{l_i}/
#' \code{Y_all} candidate/event-type branch, rather than one vectorized
#' \code{dmvnorm()} call per patient/\code{l_i} as \code{conditionalYDTBio()}
#' does.
#'
#' @inheritParams conditionalYDTBio
#' @keywords internal
conditionalYDTBioCopula <- function(Y_all, time_new, bio_i, data_predict_all, long_fit_all,
                                     survival_fit_all, l_i, survival_variable, time_variable,
                                     survival_variable_all, survival_trans_function) {
  lfit <- long_fit_all$lfit
  Sigma <- long_fit_all$Sigma_fit
  biomarker_type <- long_fit_all$biomarker_type
  thresholds <- long_fit_all$thresholds
  long_sub_fixed <- long_fit_all$long_sub_fixed
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  event_type_variable <- as.character(formula(survival_fit_all$form_conditional_cr)[[2]])

  n_longitudinal <- length(lfit)
  sigma.longitudinal <- vapply(seq_len(n_longitudinal), function(i) {
    if (biomarker_type[i] == "continuous") lfit[[i]]$sigma else 1
  }, numeric(1))

  data.long <- data_predict_all
  if (!is.list(data.long) || is.data.frame(data.long)) {
    data.long <- list(data.long)
    data.long <- rep(data.long, each = n_longitudinal)
  }

  bio_i_name <- as.character(formula(long_sub_fixed[[bio_i]])[[2]])

  f_Y_T_D_w1 <- list(); f_Y_T_D_w0 <- list()
  for (Y_i in seq_along(Y_all)) {
    f_Y_T_D_w1[[Y_i]] <- matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
    f_Y_T_D_w0[[Y_i]] <- matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
  }

  iii <- 0
  for (num_i in unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])) {
    iii <- iii + 1

    patient_data <- select_patient_longitudinal_data_bio(data.long, num, num_i, n_longitudinal, time_variable,
                                                           bio_i, time_new, Y_all, long_fit_all)
    if (is.null(patient_data)) next
    rep_num_i_list <- patient_data$rep_num_i_list
    data_num_i_list <- patient_data$data_num_i_list
    n_rows_bio_i <- nrow(data_num_i_list[[bio_i]])

    for (it in seq_along(l_i)) {
      data_it_1 <- data_num_i_list
      data_it_0 <- data_num_i_list
      for (i in 1:n_longitudinal) {
        data_it_1[[i]][survival_variable] <- l_i[it]
        data_it_0[[i]][survival_variable] <- l_i[it]
        if (length(survival_variable_all) != 0) {
          for (surv_i in seq_along(survival_variable_all)) {
            trans_val <- apply_survival_trans(survival_trans_function[[surv_i]], l_i[it], surv_i)
            data_it_1[[i]][survival_variable_all[[surv_i]]] <- trans_val
            data_it_0[[i]][survival_variable_all[[surv_i]]] <- trans_val
          }
        }
        data_it_1[[i]][event_type_variable] <- 1
        data_it_0[[i]][event_type_variable] <- 0

        all_variables <- all.vars(formula(long_sub_fixed[[i]]))
        if (!(survival_variable %in% all_variables)) {
          stop("Error: Condition is false. Please add survival variable to linear mixed model.")
        }

        if (biomarker_type[i] == "continuous") {
          resp_i <- as.character(formula(long_sub_fixed[[i]])[[2]])
          data_it_1[[i]][[resp_i]][is.na(data_it_1[[i]][[resp_i]])] <- 999
          data_it_0[[i]][[resp_i]][is.na(data_it_0[[i]][[resp_i]])] <- 999
        }
      }

      for (Y_i in seq_along(Y_all)) {
        data_it_1_Y <- data_it_1
        data_it_0_Y <- data_it_0
        data_it_1_Y[[bio_i]][[bio_i_name]][n_rows_bio_i] <- Y_all[Y_i]
        data_it_0_Y[[bio_i]][[bio_i_name]][n_rows_bio_i] <- Y_all[Y_i]

        ### the random-effects design (A_i)/Sigma_all does not depend on the
        ### event-type covariate value, so it is identical for the w0/w1
        ### branches -- only the mean vector differs (see
        ### build_conditional_design_copula()).
        design_1 <- build_conditional_design_copula(rep_num_i_list, data_it_1_Y, lfit, Sigma,
                                                      sigma.longitudinal, time_variable, n_longitudinal,
                                                      biomarker_type, long_sub_fixed, thresholds,
                                                  long_fit_all$long_sub_random)
        design_0 <- build_conditional_design_copula(rep_num_i_list, data_it_0_Y, lfit, Sigma,
                                                      sigma.longitudinal, time_variable, n_longitudinal,
                                                      biomarker_type, long_sub_fixed, thresholds,
                                                  long_fit_all$long_sub_random)

        LME_indi_matrix_1 <- lapply(seq_len(n_longitudinal), function(i) {
          build_LME_indi_matrix_copula(i, data_it_1_Y[[i]], lfit, long_fit_all, long_sub_fixed)
        })
        LME_indi_matrix_0 <- lapply(seq_len(n_longitudinal), function(i) {
          build_LME_indi_matrix_copula(i, data_it_0_Y[[i]], lfit, long_fit_all, long_sub_fixed)
        })
        LME_all_matrix_1 <- build_block_diagonal_matrix(LME_indi_matrix_1)
        LME_all_matrix_0 <- build_block_diagonal_matrix(LME_indi_matrix_0)

        mu_full_1 <- rowSums(t(LME_all_matrix_1) %*% design_1$parameter_matrix)
        mu_full_0 <- rowSums(t(LME_all_matrix_0) %*% design_0$parameter_matrix)

        cc_idx <- which(design_1$row_marker_type == "continuous")
        oo_idx <- which(design_1$row_marker_type == "ordinal")

        f_Y_T_D_w1[[Y_i]][it, iii] <- mixed_density_prob_copula(design_1$Sigma_all, mu_full_1, design_1$y_all_vec,
                                                                    design_1$alpha_lower_vec, design_1$alpha_upper_vec,
                                                                    cc_idx, oo_idx)
        f_Y_T_D_w0[[Y_i]][it, iii] <- mixed_density_prob_copula(design_0$Sigma_all, mu_full_0, design_0$y_all_vec,
                                                                    design_0$alpha_lower_vec, design_0$alpha_upper_vec,
                                                                    cc_idx, oo_idx)
      }
    }
  }
  return(f_Y_T_D = list(f_Y_T_D_w0, f_Y_T_D_w1))
}
