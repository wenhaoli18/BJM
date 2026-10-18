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

  ### one patient at a time, possibly on several cores (see patient_lapply())
  patient_ids <- unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])
  results <- patient_lapply(patient_ids, function(num_i) {
    patient_data <- select_patient_longitudinal_data(data.long, num, num_i, n_longitudinal, time_variable)
    if (is.null(patient_data)) return(NULL)
    ### all grid points at once (see copula_logdens_grid())
    copula_logdens_grid(patient_data$data_num_i_list, patient_data$rep_num_i_list, long_fit_all, l_i,
                        survival_variable, time_variable, survival_variable_all,
                        survival_trans_function, sigma.longitudinal)[[1]][, 1]
  })
  for (iii in seq_along(results)) {
    if (!is.null(results[[iii]])) f_Y_T_D_w1[, iii] <- results[[iii]]
  }
  return(f_Y_T_D = list(f_Y_T_D_w1))
}
