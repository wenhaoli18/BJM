#' Check whether bandcount1/bandcount2/bandcount3 are large enough
#'
#' @description \code{dynamicPrediction()}/\code{dynamicPredictionBio()}
#' approximate the integrals behind the predicted risk probabilities (and,
#' for \code{dynamicPredictionBio()}, the predicted biomarker density) with a
#' finite grid, controlled by \code{bandcount1}/\code{bandcount2}/
#' \code{bandcount3}. There is no universal correct value: too few grid
#' points silently bias the answer, and too many just cost more time, and
#' the right value depends on the data (e.g. how wide the follow-up range
#' is). Rather than guess, or auto-loop until some tolerance is met --
#' which multiplies runtime unpredictably, especially for
#' \code{dynamicPredictionBio()}'s nested \code{bandcount2} x
#' \code{bandcount3} grid -- this runs the prediction once at the
#' bandcount value(s) you supply, once more with those value(s) scaled up,
#' and reports the largest relative change between the two, so you can see
#' whether you have already converged or need to increase the checked
#' bandcount(s) and re-run. It costs exactly 2 calls to \code{predict_fun},
#' regardless of how many times you invoke it.
#'
#' @param predict_fun The prediction function to check: \code{dynamicPrediction}
#' or \code{dynamicPredictionBio} themselves (not a string, and not
#' \code{predictPlot()}/\code{riskPlot()}, which return a plot rather than
#' the underlying numeric predictions this function compares).
#' @param ... Arguments to forward to \code{predict_fun}, exactly as you
#' would call it directly, except for the bandcount argument(s) being
#' checked, which are supplied separately via \code{bandcount_args}. Any
#' bandcount argument of \code{predict_fun} not named in
#' \code{bandcount_args} must still be supplied here (it is held fixed at
#' that value for both calls).
#' @param bandcount_args A named list giving the bandcount value(s) to
#' check, e.g. \code{list(bandcount1 = 10, bandcount2 = 40)}. Every element
#' is scaled by \code{multiplier} for the second call. To isolate which
#' bandcount is driving instability, check one at a time (e.g.
#' \code{list(bandcount2 = 40)}, with \code{bandcount3} passed as a fixed
#' value via \code{...}), rather than checking all of them together.
#' @param multiplier Factor the checked bandcount value(s) are scaled by
#' for the second call. Must be greater than 1. Defaults to \code{2}
#' (doubling).
#' @param tol Relative-change tolerance below which the result is reported
#' as converged. Defaults to \code{0.01} (1%).
#'
#' @return An object of class \code{"checkBandcountConvergence.BJM"}, a
#' named list with elements:
#' \describe{
#'   \item{base_bandcount}{The bandcount value(s) checked, as supplied in
#'   \code{bandcount_args}.}
#'   \item{scaled_bandcount}{\code{base_bandcount} scaled by
#'   \code{multiplier}.}
#'   \item{tol}{The relative-change tolerance used to decide \code{converged}.}
#'   \item{by_field}{A named numeric vector giving the largest relative
#'   change, per comparable output field (e.g. \code{risk_prob_1},
#'   \code{Y_predict}), between the base and scaled call.}
#'   \item{max_rel_diff}{The largest value in \code{by_field}, or \code{NA}
#'   if there was no comparable output (e.g. \code{horizon <= 0}).}
#'   \item{converged}{\code{TRUE} if \code{max_rel_diff < tol}, \code{FALSE}
#'   if not, \code{NA} if there was nothing to compare.}
#'   \item{base_result}{The full result of calling \code{predict_fun} at
#'   \code{base_bandcount}.}
#'   \item{scaled_result}{The full result of calling \code{predict_fun} at
#'   \code{scaled_bandcount}.}
#' }
#' Printing the object summarizes these fields.
#'
#' @examples
#'
#' \donttest{
#' data(pbc3)
#'
#' data_survival_fitting = pbc3[!duplicated(pbc3$id), ]
#' survival_fit_all = survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
#'
#' long_sub_fixed = serBilir ~ year + age + sex + (years) + (years) * year
#' long_sub_random = ~ year | id
#' long_fit_all = longitudinalSub(pbc3[pbc3$status3 == 1, ], long_sub_fixed, long_sub_random)
#'
#' trans = survivalTrans(c(1, 3, 5, 7))
#'
#' data_predict_all = pbc3[pbc3$id == 2 & pbc3$year <= 3, ]
#'
#' # Check bandcount1/bandcount2 together: are they both already large enough?
#' check = checkBandcountConvergence(
#'   dynamicPrediction, data_predict_all, long_fit_all, survival_fit_all,
#'   prediction_time = 3, horizon = 3, time_variable = "year",
#'   trans$survival_variable_all, trans$survival_trans_function,
#'   bandcount_args = list(bandcount1 = 10, bandcount2 = 10)
#' )
#' print(check)
#' }
#'
#' @export
checkBandcountConvergence <- function(predict_fun, ..., bandcount_args, multiplier = 2, tol = 0.01) {
  is_dynamic_prediction <- identical(predict_fun, dynamicPrediction)
  is_dynamic_prediction_bio <- identical(predict_fun, dynamicPredictionBio)
  if (!is_dynamic_prediction && !is_dynamic_prediction_bio) {
    stop(paste(
      "`predict_fun` must be `dynamicPrediction` or `dynamicPredictionBio`",
      "(the function itself, not a string, and not `predictPlot`/`riskPlot`,",
      "which return a plot rather than the underlying numeric predictions",
      "this function compares)."
    ), call. = FALSE)
  }
  predict_fun_name <- if (is_dynamic_prediction) "dynamicPrediction" else "dynamicPredictionBio"

  if (!is.list(bandcount_args) || length(bandcount_args) == 0 ||
      is.null(names(bandcount_args)) || any(names(bandcount_args) == "")) {
    stop("`bandcount_args` must be a non-empty named list, e.g. list(bandcount1 = 10, bandcount2 = 40).", call. = FALSE)
  }
  valid_names <- grep("^bandcount", names(formals(predict_fun)), value = TRUE)
  bad_names <- setdiff(names(bandcount_args), valid_names)
  if (length(bad_names) > 0) {
    stop(sprintf(
      "`bandcount_args` names %s are not bandcount argument(s) of %s(); valid options are: %s.",
      paste(sprintf("`%s`", bad_names), collapse = ", "), predict_fun_name, paste(valid_names, collapse = ", ")
    ), call. = FALSE)
  }
  not_positive_scalar <- !vapply(bandcount_args, function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0, logical(1))
  if (any(not_positive_scalar)) {
    stop(sprintf("`bandcount_args$%s` must be a single positive number.",
                 names(bandcount_args)[which(not_positive_scalar)[1]]), call. = FALSE)
  }
  assert_scalar_numeric(multiplier, "multiplier", positive = TRUE)
  if (multiplier <= 1) {
    stop("`multiplier` must be greater than 1, so the scaled call actually checks a larger bandcount.", call. = FALSE)
  }
  assert_scalar_numeric(tol, "tol", positive = TRUE)

  extra_args <- list(...)
  scaled_bandcount_args <- lapply(bandcount_args, function(x) x * multiplier)

  base_result <- do.call(predict_fun, c(extra_args, bandcount_args))
  scaled_result <- do.call(predict_fun, c(extra_args, scaled_bandcount_args))

  # See max_relative_diff() for which fields are comparable and why.
  comparison <- max_relative_diff(base_result, scaled_result)
  by_field <- comparison$by_field
  max_rel_diff <- comparison$max
  converged <- if (is.na(max_rel_diff)) NA else (max_rel_diff < tol)

  out <- list(
    base_bandcount = bandcount_args,
    scaled_bandcount = scaled_bandcount_args,
    tol = tol,
    by_field = by_field,
    max_rel_diff = max_rel_diff,
    converged = converged,
    base_result = base_result,
    scaled_result = scaled_result
  )
  class(out) <- "checkBandcountConvergence.BJM"
  out
}

#' @export
print.checkBandcountConvergence.BJM <- function(x, ...) {
  format_bandcount <- function(b) paste(names(b), unlist(b), sep = " = ", collapse = ", ")
  cat("Bandcount convergence check\n")
  cat(sprintf("  base:   %s\n", format_bandcount(x$base_bandcount)))
  cat(sprintf("  scaled: %s\n", format_bandcount(x$scaled_bandcount)))
  if (length(x$by_field) > 0) {
    for (field in names(x$by_field)) {
      cat(sprintf("  max relative change in %s: %s\n", field, format(x$by_field[[field]], digits = 3)))
    }
  }
  if (is.na(x$converged)) {
    cat("  Could not compare (no comparable numeric output was returned, e.g. horizon <= 0).\n")
  } else if (x$converged) {
    cat(sprintf(
      "  Converged: max relative change (%s) is below tol (%s); the base bandcount value(s) look adequate.\n",
      format(x$max_rel_diff, digits = 3), format(x$tol, digits = 3)
    ))
  } else {
    cat(sprintf(
      "  NOT converged: max relative change (%s) exceeds tol (%s); consider increasing the checked bandcount(s) and re-running.\n",
      format(x$max_rel_diff, digits = 3), format(x$tol, digits = 3)
    ))
  }
  invisible(x)
}
