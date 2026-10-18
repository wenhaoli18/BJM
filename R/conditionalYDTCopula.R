#' Conditional distribution of Y|D,T for a mixed continuous/ordinal
#' (Gaussian-copula) joint model, if with competing risk
#'
#' @description Copula-aware counterpart to \code{conditionalYDT()}, used by
#' \code{predictRisk()} whenever \code{long_fit_all$biomarker_type}
#' contains at least one \code{"ordinal"} biomarker (see
#' \code{longitudinalSubCopula()}). See \code{conditionalYTCopula()} for the
#' mixed continuous/ordinal density/probability construction shared by both
#' event-type branches (\code{w0}/\code{w1}) here.
#'
#' @inheritParams conditionalYDT
#' @keywords internal
conditionalYDTCopula <- function(data_predict_all, long_fit_all, survival_fit_all,
                                  l_i, survival_variable,
                                  time_variable, survival_variable_all, survival_trans_function) {
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

  f_Y_T_D_w1 <- matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
  f_Y_T_D_w0 <- matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))

  ### one patient at a time, possibly on several cores (see patient_lapply())
  patient_ids <- unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])
  results <- patient_lapply(patient_ids, function(num_i) {
    patient_data <- select_patient_longitudinal_data(data.long, num, num_i, n_longitudinal, time_variable)
    if (is.null(patient_data)) return(NULL)
    ### all grid points and both event types at once (see copula_logdens_grid())
    copula_logdens_grid(patient_data$data_num_i_list, patient_data$rep_num_i_list, long_fit_all, l_i,
                        survival_variable, time_variable, survival_variable_all,
                        survival_trans_function, sigma.longitudinal,
                        event_type_variable = event_type_variable)
  })
  for (iii in seq_along(results)) {
    if (is.null(results[[iii]])) next
    f_Y_T_D_w1[, iii] <- results[[iii]][[1]][, 1]
    f_Y_T_D_w0[, iii] <- results[[iii]][[2]][, 1]
  }
  return(f_Y_T_D = list(f_Y_T_D_w0, f_Y_T_D_w1))
}
