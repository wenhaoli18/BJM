#' Assert that an object is a non-empty data.frame
#'
#' @description Shared input-validation helper: raises a clear error instead
#' of letting a malformed argument fail deep inside model-fitting code with
#' a cryptic message.
#' @keywords internal
assert_data_frame <- function(x, arg_name) {
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a data.frame, not %s.", arg_name, class(x)[1]), call. = FALSE)
  }
  if (nrow(x) == 0) {
    stop(sprintf("`%s` has zero rows.", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that an object is a list of data.frame objects of the expected length
#'
#' @description Shared input-validation helper for the \code{data_fit_all}/
#' \code{data_predict_all} arguments, which must be a list with one
#' data.frame per longitudinal outcome. When \code{allow_bare_df = TRUE}, a
#' single bare data.frame is also accepted, matching the auto-repeat
#' convenience some functions apply for a single data.frame shared across
#' every outcome.
#' @keywords internal
assert_data_list <- function(x, arg_name, n_expected, allow_bare_df = FALSE) {
  if (is.data.frame(x)) {
    if (allow_bare_df) {
      return(invisible(TRUE))
    }
    stop(sprintf(
      "`%s` must be a list of data.frame objects (one per longitudinal outcome), not a bare data.frame. Wrap it, e.g. list(%s).",
      arg_name, arg_name
    ), call. = FALSE)
  }
  if (!is.list(x) || length(x) == 0) {
    stop(sprintf("`%s` must be a list of data.frame objects, not %s.", arg_name, class(x)[1]), call. = FALSE)
  }
  not_df <- !vapply(x, is.data.frame, logical(1))
  if (any(not_df)) {
    stop(sprintf("Every element of `%s` must be a data.frame; element %d is not.",
                 arg_name, which(not_df)[1]), call. = FALSE)
  }
  if (length(x) != n_expected) {
    stop(sprintf(
      "`%s` has %d element(s), but there %s %d longitudinal outcome(s); `%s` must have exactly %d element(s) (one per outcome), or be a single bare data.frame reused for every outcome.",
      arg_name, length(x), if (n_expected == 1) "is" else "are", n_expected, arg_name, n_expected
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that every element of a list is a formula
#'
#' @description Shared input-validation helper for \code{long_sub_fixed}/
#' \code{long_sub_random}-style arguments, after they have been normalized
#' to a list (a bare formula is wrapped in \code{list()} by the caller
#' before calling this).
#' @keywords internal
assert_all_formulas <- function(x, arg_name) {
  if (!is.list(x) || length(x) == 0) {
    stop(sprintf("`%s` must be a formula or a non-empty list of formulas.", arg_name), call. = FALSE)
  }
  not_formula <- !vapply(x, inherits, logical(1), what = "formula")
  if (any(not_formula)) {
    stop(sprintf("Every element of `%s` must be a formula; element %d is not.",
                 arg_name, which(not_formula)[1]), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a formula's variables are all present in a data.frame
#'
#' @description Shared input-validation helper: catches missing columns
#' before they surface as an opaque error from deep inside \code{coxph()},
#' \code{lme()}, or \code{model.matrix()}.
#' @keywords internal
assert_vars_in_data <- function(vars, data, source_name, data_name) {
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variable(s) used in `%s` not found in `%s`: %s.",
                 source_name, data_name, paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that an object is the output of a specific BJM fitting function
#'
#' @description Shared input-validation helper: checks the S3 class tag
#' attached by \code{survivalSub()}/\code{longitudinalSub()}, so passing the
#' wrong object (or the arguments in the wrong order) fails immediately with
#' a clear message instead of deep inside the prediction code.
#' @keywords internal
assert_class <- function(x, expected_class, arg_name, expected_source) {
  if (!inherits(x, expected_class)) {
    stop(sprintf("`%s` must be the output of %s(), not %s.",
                 arg_name, expected_source, class(x)[1]), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that an object is a single, non-missing numeric value
#'
#' @keywords internal
assert_scalar_numeric <- function(x, arg_name, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("`%s` must be a single numeric value.", arg_name), call. = FALSE)
  }
  if (positive && x <= 0) {
    stop(sprintf("`%s` must be a positive number.", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that an object is a single, non-missing character string
#'
#' @keywords internal
assert_string <- function(x, arg_name) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("`%s` must be a single character string.", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that an index is a valid, in-range biomarker position
#'
#' @description Shared input-validation helper for \code{bio_i}: catches an
#' out-of-range or non-integer value before it becomes a "subscript out of
#' bounds" error from indexing into \code{data_predict_all}/\code{lfit}.
#' @keywords internal
assert_index <- function(x, max_value, arg_name, context) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x)) {
    stop(sprintf("`%s` must be a single integer.", arg_name), call. = FALSE)
  }
  if (x < 1 || x > max_value) {
    stop(sprintf("`%s` must be between 1 and %d (the number of %s), not %s.",
                 arg_name, max_value, context, x), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that survival_variable_all/survival_trans_function are consistent
#'
#' @description Shared input-validation helper for \code{dynamicPrediction()}
#' and \code{dynamicPredictionBio()}: the two arguments must have matching
#' length, and every transform must be a function.
#' @keywords internal
assert_survival_trans <- function(survival_variable_all, survival_trans_function) {
  if (length(survival_variable_all) != length(survival_trans_function)) {
    stop(sprintf(
      "`survival_variable_all` has %d element(s) but `survival_trans_function` has %d; they must have the same length.",
      length(survival_variable_all), length(survival_trans_function)
    ), call. = FALSE)
  }
  if (length(survival_trans_function) > 0) {
    not_fun <- !vapply(survival_trans_function, is.function, logical(1))
    if (any(not_fun)) {
      stop(sprintf("Every element of `survival_trans_function` must be a function; element %d is not.",
                   which(not_fun)[1]), call. = FALSE)
    }
  }
  invisible(TRUE)
}
