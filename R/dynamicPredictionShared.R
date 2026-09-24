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

#' Compare two prediction results' plain numeric-vector fields
#'
#' @description Shared comparison logic for \code{checkBandcountConvergence()}
#' and the \code{"auto"} bandcount support in \code{dynamicPrediction()}/
#' \code{dynamicPredictionBio()}: only plain numeric vectors (no \code{dim})
#' that have the same length in both results are compared. This naturally
#' skips fields whose *size* is itself controlled by the bandcount being
#' varied (e.g. \code{dynamicPredictionBio()}'s \code{Y_density} matrix and
#' \code{Y_all} grid, whose resolution is exactly what \code{bandcount3}
#' sets), while still comparing the actual per-patient estimates derived
#' from them (\code{risk_prob_1}/\code{risk_prob_2}, \code{Y_predict}).
#' @return A list with \code{max} (the largest relative change across all
#' comparable fields, or \code{NA} if none were comparable) and
#' \code{by_field} (a named numeric vector, one entry per comparable
#' field).
#' @keywords internal
max_relative_diff <- function(result_a, result_b) {
  is_comparable <- function(x) is.numeric(x) && is.null(dim(x))
  common_fields <- intersect(
    names(result_a)[vapply(result_a, is_comparable, logical(1))],
    names(result_b)[vapply(result_b, is_comparable, logical(1))]
  )
  by_field <- c()
  for (field in common_fields) {
    a <- result_a[[field]]
    b <- result_b[[field]]
    if (length(a) == 0 || length(a) != length(b)) next
    by_field[field] <- max(abs(a - b) / pmax(abs(a), 1e-8))
  }
  list(max = if (length(by_field) == 0) NA_real_ else max(by_field), by_field = by_field)
}

#' Starting values for "auto" bandcount doubling
#'
#' @description The built-in starting point \code{auto_tune_bandcount()}
#' doubles from for each bandcount argument. \code{bandcount1}/
#' \code{bandcount2}/\code{bandcount3} used to default to fixed numbers
#' (\code{10}, \code{40}, and \code{300} respectively, across
#' \code{dynamicPrediction()}/\code{dynamicPredictionBio()}); those same
#' numbers are reused here as the starting point for auto-tuning, so that
#' the first call \code{auto_tune_bandcount()} makes matches what a caller
#' relying on the old fixed defaults would have gotten. This cannot instead
#' be read off \code{formals(predict_fun)}, because that default is now the
#' literal string \code{"auto"} itself.
#' @keywords internal
bandcount_auto_start <- c(bandcount1 = 10, bandcount2 = 40, bandcount3 = 300)

#' Auto-select "auto" bandcount arguments by doubling until convergence
#'
#' @description Shared implementation backing \code{bandcount1}/
#' \code{bandcount2}/\code{bandcount3 = "auto"} support in
#' \code{dynamicPrediction()}/\code{dynamicPredictionBio()}, and the
#' bandcount pre-resolution done once, up front, by \code{predictPlot()}/
#' \code{riskPlot()} (so their internal horizon/landmark loops do not repeat
#' the auto-tuning search on every iteration).
#'
#' Starts every argument named in \code{auto_names} at its entry in
#' \code{bandcount_auto_start()}, doubles all of them together each round,
#' and compares consecutive results with \code{max_relative_diff()} until
#' the largest relative change drops below \code{tol}, or \code{max_rounds}
#' extra doublings have been tried -- a hard cap, so this never loops
#' indefinitely: at most \code{max_rounds + 1} calls to \code{predict_fun}
#' (the default \code{max_rounds = 2} means at most 3 calls). If the cap is
#' hit without converging, a warning is issued and the result/bandcount at
#' the largest value tried is returned anyway, rather than erroring, so
#' automated pipelines are not interrupted.
#'
#' @param predict_fun \code{dynamicPrediction} or \code{dynamicPredictionBio}.
#' @param args A named list of all of \code{predict_fun}'s arguments
#' (typically \code{as.list(environment())} captured right after argument
#' validation, before any other local variables are created).
#' @param auto_names Character vector naming which element(s) of \code{args}
#' to auto-tune (e.g. \code{"bandcount1"}, or \code{c("bandcount1", "bandcount2")}).
#' @return A list with \code{result} (\code{predict_fun}'s return value at
#' the resolved bandcount) and \code{bandcount} (a named list of the
#' resolved numeric bandcount value(s), one per element of \code{auto_names}).
#' @keywords internal
auto_tune_bandcount <- function(predict_fun, args, auto_names, tol = 0.01, max_rounds = 2, multiplier = 2) {
  for (n in auto_names) {
    args[[n]] <- unname(bandcount_auto_start[n])
  }

  prev_result <- do.call(predict_fun, args)
  round_i <- 0
  repeat {
    scaled_args <- args
    for (n in auto_names) scaled_args[[n]] <- scaled_args[[n]] * multiplier
    new_result <- do.call(predict_fun, scaled_args)
    round_i <- round_i + 1

    comparison <- max_relative_diff(prev_result, new_result)
    converged <- !is.na(comparison$max) && comparison$max < tol

    args <- scaled_args
    prev_result <- new_result

    if (converged || round_i >= max_rounds) {
      if (!converged) {
        warning(sprintf(paste0(
          "Auto-selected bandcount (%s) had not converged (max relative change %s) after ",
          "%d doubling(s) from the default; returning the result at the largest value tried ",
          "(%s). Pass an explicit, larger bandcount if you need tighter convergence, or use ",
          "checkBandcountConvergence() to investigate further."
        ), paste(auto_names, collapse = "/"),
           if (is.na(comparison$max)) "NA" else format(comparison$max, digits = 3),
           round_i,
           paste(sprintf("%s = %s", auto_names, unlist(args[auto_names])), collapse = ", ")
        ), call. = FALSE)
      }
      break
    }
  }
  list(result = prev_result, bandcount = args[auto_names])
}
