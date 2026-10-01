#' Conditional distribution of Y|T for a mixed continuous/ordinal
#' (Gaussian-copula) joint model, if no competing risk -- biomarker-value
#' prediction
#'
#' @description Copula-aware counterpart to \code{conditionalYTBio()}, used
#' by \code{dynamicPredictionBio()} whenever \code{long_fit_all$biomarker_type}
#' contains at least one \code{"ordinal"} biomarker (see
#' \code{longitudinalSubCopula()}). Both a continuous and an
#' \strong{ordinal} \code{bio_i} target are supported: for an ordinal
#' \code{bio_i}, \code{Y_all} is expected to already be the vector of
#' candidate \emph{category labels} (the fitted factor's \code{levels()}, in
#' threshold order) rather than a numeric grid -- see
#' \code{compute_bio_marker_step()}, which builds that vector and translates
#' the result back into integer category codes for its caller. No branching
#' is actually needed here: \code{build_conditional_design_copula()} already
#' looks up each row of biomarker \code{bio_i} (historical \strong{and} the
#' candidate row assigned below) generically via
#' \code{long_fit_all$biomarker_type[bio_i]}, bracketing an ordinal row's
#' latent score between its category's cumulative-link thresholds instead of
#' matching it to an exact value -- exactly the same mechanism already used
#' for every \emph{other} ordinal biomarker's observed history in the joint
#' density. Assigning a candidate label into \code{data_it_Y}'s ordinal
#' factor column (below) therefore evaluates the same
#' \code{mixed_density_prob_copula()} box probability as any other ordinal
#' row would, with no separate code path required.
#'
#' Unlike \code{conditionalYTBio()} -- which can evaluate
#' \code{mvtnorm::dmvnorm()} once per \code{l_i} across every \code{Y_all}
#' candidate in one vectorized call, because the multivariate-normal density
#' does not need re-normalizing per candidate -- here each candidate value in
#' \code{Y_all} changes the conditional block's mean/covariance (see
#' \code{mixed_density_prob_copula()}) and therefore requires its own
#' \code{mvtnorm::pmvnorm()} evaluation. This makes the copula path in this
#' function \code{O(patients * length(l_i) * length(Y_all))} Monte-Carlo box
#' probability evaluations, rather than one \code{dmvnorm()} call per
#' patient/\code{l_i} -- keep \code{bandcount3} modest for mixed fits with a
#' continuous \code{bio_i} (an ordinal \code{bio_i}'s \code{Y_all} length is
#' fixed at its category count, not controlled by \code{bandcount3} at all).
#'
#' @inheritParams conditionalYTBio
#' @keywords internal
conditionalYTBioCopula <- function(Y_all, time_new, bio_i, data_predict_all, long_fit_all, l_i,
                                    survival_variable, time_variable, survival_variable_all,
                                    survival_trans_function) {
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

  bio_i_name <- as.character(formula(long_sub_fixed[[bio_i]])[[2]])

  f_Y_T_D_w1 <- list()
  for (Y_i in seq_along(Y_all)) {
    f_Y_T_D_w1[[Y_i]] <- matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
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
          resp_i <- as.character(formula(long_sub_fixed[[i]])[[2]])
          data_it[[i]][[resp_i]][is.na(data_it[[i]][[resp_i]])] <- 999
        }
      }

      for (Y_i in seq_along(Y_all)) {
        data_it_Y <- data_it
        data_it_Y[[bio_i]][[bio_i_name]][n_rows_bio_i] <- Y_all[Y_i]

        design <- build_conditional_design_copula(rep_num_i_list, data_it_Y, lfit, Sigma,
                                                    sigma.longitudinal, time_variable, n_longitudinal,
                                                    biomarker_type, long_sub_fixed, thresholds,
                                                  long_fit_all$long_sub_random)

        LME_indi_matrix <- lapply(seq_len(n_longitudinal), function(i) {
          build_LME_indi_matrix_copula(i, data_it_Y[[i]], lfit, long_fit_all, long_sub_fixed)
        })
        LME_all_matrix <- build_block_diagonal_matrix(LME_indi_matrix)

        mu_full <- rowSums(t(LME_all_matrix) %*% design$parameter_matrix)

        cc_idx <- which(design$row_marker_type == "continuous")
        oo_idx <- which(design$row_marker_type == "ordinal")

        f_Y_T_D_w1[[Y_i]][it, iii] <- mixed_density_prob_copula(design$Sigma_all, mu_full, design$y_all_vec,
                                                                    design$alpha_lower_vec, design$alpha_upper_vec,
                                                                    cc_idx, oo_idx)
      }
    }
  }
  return(f_Y_T_D = list(f_Y_T_D_w1))
}
