#' Restrict prediction data to patients still at risk
#'
#' @description Shared helper for \code{dynamicPrediction} and
#' \code{dynamicPredictionBio}: drops rows whose survival-time variable is
#' below \code{prediction_time} from every biomarker's data frame.
#'
#' @return \code{data_predict_all}, filtered in place per element.
#' @keywords internal
subset_at_risk <- function(data_predict_all, survival_variable, prediction_time) {
  for (i in seq_len(length(data_predict_all))) {
    data_predict_all[[i]] = data_predict_all[[i]][data_predict_all[[i]][survival_variable] >= prediction_time, ]
  }
  data_predict_all
}

#' Build the prediction-to-infinity integration grid and marginal survival
#'
#' @description Shared helper for \code{dynamicPrediction} and
#' \code{dynamicPredictionBio}: builds the numerical-integration grid from
#' \code{prediction_time} out to \code{upper_bound}, and evaluates the
#' marginal survival function \code{S(T)} over it.
#'
#' @return A list with \code{predict.time.infinity},
#' \code{predict.time.infinity.1}, and \code{S_T_all_infinity}.
#' @keywords internal
prepare_infinity_grid <- function(data_predict_all, long_fit_all, survival_fit_all,
                                   prediction_time, upper_bound, bandcount2) {
  bandwidth2 = (upper_bound - prediction_time) / bandcount2
  predict.time.infinity = seq(prediction_time, upper_bound, bandwidth2)
  predict.time.infinity.1 = seq(prediction_time - bandwidth2 / 2, upper_bound + bandwidth2 / 2, bandwidth2)

  S_T_all_infinity = marginalT(data_predict_all, long_fit_all, survival_fit_all,
                                l_i = predict.time.infinity.1, upper_bound)

  list(predict.time.infinity = predict.time.infinity,
       predict.time.infinity.1 = predict.time.infinity.1,
       S_T_all_infinity = S_T_all_infinity)
}

#' Normalize a risk-probability ratio into a valid probability
#'
#' @description Shared helper for \code{dynamicPrediction} and
#' \code{dynamicPredictionBio}: divides summed predicted-event mass by
#' summed total mass and clamps the result to \code{[0, 1]}.
#'
#' @return A numeric vector of risk probabilities in \code{[0, 1]}.
#' @keywords internal
clamp_risk_prob <- function(numerator_sum, denominator_sum) {
  risk.prob <- numerator_sum / denominator_sum
  risk.prob[risk.prob > 1] = 1
  risk.prob[risk.prob < 0] = 0
  risk.prob
}
