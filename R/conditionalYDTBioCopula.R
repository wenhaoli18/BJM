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

    ### all grid points, candidate values and both event types at once (see
    ### copula_logdens_grid())
    log_dens <- copula_logdens_grid(data_num_i_list, rep_num_i_list, long_fit_all, l_i,
                                    survival_variable, time_variable, survival_variable_all,
                                    survival_trans_function, sigma.longitudinal,
                                    event_type_variable = event_type_variable,
                                    bio = list(index = bio_i, name = bio_i_name, row = n_rows_bio_i,
                                               values = Y_all))
    for (Y_i in seq_along(Y_all)) {
      f_Y_T_D_w1[[Y_i]][, iii] <- log_dens[[1]][, Y_i]
      f_Y_T_D_w0[[Y_i]][, iii] <- log_dens[[2]][, Y_i]
    }
  }
  return(f_Y_T_D = list(f_Y_T_D_w0, f_Y_T_D_w1))
}
