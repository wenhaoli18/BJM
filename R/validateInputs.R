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

#' Assert that a bandcount argument is a positive number or "auto"
#'
#' @description Shared input-validation helper for the \code{bandcount1}/
#' \code{bandcount2}/\code{bandcount3} arguments of \code{predictRisk()},
#' \code{dynamicPredictionBio()}, \code{predictPlot()}, and \code{riskPlot()},
#' which now accept either an explicit positive number (the original
#' behavior) or the literal string \code{"auto"} to have the value chosen
#' automatically (see \code{auto_tune_bandcount()}).
#' @keywords internal
assert_bandcount <- function(x, arg_name) {
  if (identical(x, "auto")) {
    return(invisible(TRUE))
  }
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive number, or the string \"auto\".", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that an object is a single positive integer
#'
#' @description Shared input-validation helper for \code{n_cores}: catches a
#' non-integer, zero, negative, or non-scalar value before it reaches
#' \code{parallel::mclapply()}'s own (less informative) \code{mc.cores}
#' validation.
#' @keywords internal
assert_positive_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x) || x < 1) {
    stop(sprintf("`%s` must be a single positive integer.", arg_name), call. = FALSE)
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

#' Warn about formula terms whose basis is recomputed from whatever data
#' they are given
#'
#' @description \code{poly()} (in its default orthogonal mode),
#' \code{splines::ns()}/\code{splines::bs()}, and \code{factor()} compute
#' their basis/contrasts from whatever data is passed to \code{model.matrix()}.
#' For a continuous biomarker's \code{long_sub_fixed} formula this is
#' handled: the basis \code{lme()} was fit with (its \code{terms}'
#' \code{predvars}) and the training factor levels are reused both when
#' estimating \code{Sigma_fit} and at prediction time, so such terms are not
#' flagged there. They are still rebuilt from the data at hand -- a single
#' patient's rows at prediction time, the retained subjects when estimating
#' \code{Sigma_fit} -- for \code{long_sub_random} formulas and for an
#' ordinal biomarker's \code{long_sub_fixed} formula, which is where this
#' warns. \code{poly(..., raw = TRUE)}, \code{I(x^2)}, \code{log()},
#' \code{sqrt()}, and similar terms that do not depend on the surrounding
#' data are never flagged.
#' @keywords internal
warn_unsafe_formula_terms <- function(formula_list, arg_name) {
  find_unsafe_calls <- function(expr) {
    hits <- character(0)
    if (is.call(expr)) {
      fname <- as.character(expr[[1]])
      fname <- fname[length(fname)]
      if (fname %in% c("poly", "ns", "bs", "factor")) {
        raw_arg <- tryCatch(eval(expr[["raw"]]), error = function(e) FALSE)
        if (!(fname == "poly" && isTRUE(raw_arg))) {
          hits <- c(hits, deparse(expr))
        }
      }
      for (a in as.list(expr)[-1]) {
        hits <- c(hits, find_unsafe_calls(a))
      }
    }
    hits
  }

  for (i in seq_along(formula_list)) {
    hits <- unique(find_unsafe_calls(formula_list[[i]]))
    if (length(hits) > 0) {
      warning(sprintf(
        "`%s[[%d]]` uses %s. Its basis/contrasts depend on the data it is computed from, and for this formula BJM rebuilds the design matrix from the data at hand (a single patient's rows at prediction time) instead of reusing the training-time basis -- so this can silently produce incorrect predictions, or fail outright when there are too few distinct values. Prefer poly(..., raw = TRUE), I(x^2), log(), sqrt(), or other terms that do not depend on the surrounding data.",
        arg_name, i, paste(hits, collapse = ", ")
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Assert that survival_variable_all/survival_trans_function are consistent
#'
#' @description Shared input-validation helper for \code{predictRisk()},
#' \code{dynamicPredictionBio()}, \code{predictPlot()}, and \code{riskPlot()}:
#' the two arguments must have matching length, and every transform must be
#' a function. When \code{probe_value} is supplied, every transform is also
#' test-called once on it, and must return a single, finite, non-missing
#' numeric value. Without this, a transform that throws an error, or
#' returns a character value, a length != 1 vector, or a non-finite value
#' (e.g. \code{log(x)} evaluated at \code{x <= 0}), would only surface deep
#' inside the per-patient prediction grid built by \code{conditionalYT()}/
#' \code{conditionalYDT()}/\code{conditionalYTBio()}/\code{conditionalYDTBio()}
#' -- as a cryptic error, or, worse, as silently corrupted data with no
#' error at all. The probe is a single call per transform, so it is cheap
#' even though the same transform is later called many times inside the
#' prediction grid.
#' @keywords internal
assert_survival_trans <- function(survival_variable_all, survival_trans_function, probe_value = NULL) {
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
    if (!is.null(probe_value)) {
      for (i in seq_along(survival_trans_function)) {
        value <- tryCatch(
          survival_trans_function[[i]](probe_value),
          error = function(e) {
            stop(sprintf(paste0(
              "`survival_trans_function[[%d]]` failed when called on a representative ",
              "time value (%s): %s. Every element of `survival_trans_function` must be a ",
              "function that accepts a single numeric time value and returns a single, ",
              "finite numeric value."
            ), i, format(probe_value), conditionMessage(e)), call. = FALSE)
          }
        )
        ok <- is.numeric(value) && length(value) == 1 && is.finite(value)
        if (!ok) {
          stop(sprintf(paste0(
            "`survival_trans_function[[%d]]` must return a single, finite numeric value; ",
            "calling it on a representative time value (%s) returned a %s of length %d instead."
          ), i, format(probe_value), class(value)[1], length(value)), call. = FALSE)
        }
      }
    }
  }
  invisible(TRUE)
}

#' Assert that an optional (Suggests-only) package is installed
#'
#' @description Shared input-validation helper for functions that depend on
#' a package listed in \code{Suggests} rather than \code{Imports} (so that
#' most users -- and CRAN's own checks -- are unaffected by a large,
#' optional dependency they never call). Every call site into that package
#' must still be fully namespaced (\code{pkg::fun()}), never
#' \code{@importFrom}; this helper just turns a missing package into a
#' clear, actionable error instead of a cryptic "could not find function"
#' failure deep inside the calling code.
#' @keywords internal
assert_package_installed <- function(pkg, context) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "%s requires the '%s' package, which is not installed. Install it with install.packages(\"%s\")%s.",
      context, pkg, pkg,
      if (pkg == "torch") ' and then run torch::install_torch() once to download LibTorch' else ""
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Call a single survival_trans_function element and validate its output
#'
#' @description \code{assert_survival_trans()} only probes each transform
#' once, at a single representative time value, before the prediction grid
#' runs. That catches a transform that is broken everywhere (wrong return
#' type/length, or throws), but not one that only misbehaves away from the
#' probe point -- e.g. \code{log(x - 10)}, which is fine near the probe value
#' but returns \code{NaN} once the internal prediction grid (which can range
#' up to \code{2 * max(observed survival time)}) reaches \code{x <= 10}. This
#' helper wraps every actual call site inside \code{conditionalYT()}/
#' \code{conditionalYTBio()}/\code{conditionalYDT()}/\code{conditionalYDTBio()}
#' so a bad value is caught immediately, with a clear error, instead of
#' silently corrupting a data.frame column or surfacing later as a cryptic
#' "replacement has ... rows, data has ..." error.
#' @keywords internal
apply_survival_trans <- function(fun, x, surv_i) {
  value <- fun(x)
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value)) {
    stop(sprintf(paste0(
      "`survival_trans_function[[%d]]` must return a single, finite numeric value; ",
      "calling it on %s returned a %s of length %d instead."
    ), surv_i, format(x), class(value)[1], length(value)), call. = FALSE)
  }
  value
}

#' Assert that biomarker_type is a valid per-biomarker type vector
#'
#' @description Shared input-validation helper for the \code{biomarker_type}
#' argument of \code{\link{longitudinalSub}}: when supplied explicitly (as
#' opposed to being auto-detected from each biomarker's response column),
#' it must be a character vector with exactly one entry per longitudinal
#' outcome, and every entry must be one of \code{"continuous"} or
#' \code{"ordinal"}.
#' @keywords internal
assert_biomarker_type <- function(biomarker_type, M, valid_types = c("continuous", "ordinal")) {
  if (!is.character(biomarker_type) || length(biomarker_type) != M || anyNA(biomarker_type)) {
    stop(sprintf(
      "`biomarker_type` must be a character vector of length %d (one entry per longitudinal outcome, in the same order as `long_sub_fixed`), with no missing values.",
      M
    ), call. = FALSE)
  }
  bad <- setdiff(biomarker_type, valid_types)
  if (length(bad) > 0) {
    stop(sprintf(
      "`biomarker_type` must only contain %s; found invalid value(s): %s.",
      paste(sprintf('"%s"', valid_types), collapse = ", "),
      paste(sprintf('"%s"', unique(bad)), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a list of longitudinalSub.BJM fits is poolable with Rubin's rules
#'
#' @description Shared input-validation helper for
#' \code{\link{poolLongitudinalSub}}: checks that \code{long_fit_all_list}
#' is a list of at least two \code{longitudinalSub.BJM} objects (the S3
#' class \code{\link{longitudinalSub}} attaches to its return value), and
#' that every one of them was fit for the same set of biomarkers with the
#' same fixed-effect coefficients -- otherwise "the mean of the estimates"
#' would be averaging unrelated quantities across completions.
#' @keywords internal
assert_poolable_longitudinal_fits <- function(long_fit_all_list, arg_name = "long_fit_all_list") {
  if (!is.list(long_fit_all_list) || is.data.frame(long_fit_all_list)) {
    stop(sprintf(
      "`%s` must be a list of `longitudinalSub.BJM` objects (fit `longitudinalSub()` once per completed dataset in `imputeLongitudinal(..., impute = \"multiple\")`'s `data_fit_all_list`), not %s.",
      arg_name, class(long_fit_all_list)[1]
    ), call. = FALSE)
  }
  if (length(long_fit_all_list) < 2) {
    stop(sprintf(
      "`%s` has %d element(s); Rubin's rules pooling needs at least 2 completed-data fits (use `n_imputations >= 2` and `impute = \"multiple\"` in `imputeLongitudinal()`).",
      arg_name, length(long_fit_all_list)
    ), call. = FALSE)
  }
  for (i in seq_along(long_fit_all_list)) {
    assert_class(long_fit_all_list[[i]], "longitudinalSub.BJM",
                 sprintf("%s[[%d]]", arg_name, i), "longitudinalSub")
  }

  M1 <- length(long_fit_all_list[[1]]$lfit)
  coef_names_1 <- lapply(long_fit_all_list[[1]]$lfit, function(f) names(nlme::fixef(f)))
  for (i in seq_along(long_fit_all_list)[-1]) {
    fit_i <- long_fit_all_list[[i]]
    if (length(fit_i$lfit) != M1) {
      stop(sprintf(
        "`%s[[%d]]` has %d longitudinal outcome(s), but `%s[[1]]` has %d; every fit must be for the same set of biomarkers (fit `longitudinalSub()` with the same `long_sub_fixed`/`long_sub_random` on each completed dataset).",
        arg_name, i, length(fit_i$lfit), arg_name, M1
      ), call. = FALSE)
    }
    coef_names_i <- lapply(fit_i$lfit, function(f) names(nlme::fixef(f)))
    if (!identical(coef_names_i, coef_names_1)) {
      stop(sprintf(
        "`%s[[%d]]` has different fixed-effect coefficient names than `%s[[1]]`; every fit must use the same `long_sub_fixed` formula(s).",
        arg_name, i, arg_name
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Check that the survival time enters the event-type model linearly
#'
#' @description Used by \code{survivalSub()}: at prediction time the
#' event-type model \code{form_conditional_cr} is evaluated at every point of
#' the survival-time integration grid by multiplying the survival time's
#' coefficient by the grid value. That is only correct when the survival
#' time appears as a plain main effect (\code{status ~ time + ...}); a
#' transformation (\code{log(time)}, \code{I(time^2)}) or an interaction
#' (\code{time:age}) would silently be evaluated at the wrong value, so it is
#' rejected here. Leaving the survival time out entirely is allowed (event
#' type then does not depend on event time).
#'
#' @param form_conditional_cr The event-type formula.
#' @param survival_variable Name of the survival-time variable.
#' @keywords internal
assert_linear_time_term <- function(form_conditional_cr, survival_variable) {
  tt <- stats::terms(form_conditional_cr)
  factors <- attr(tt, "factors")
  if (length(factors) == 0) return(invisible(TRUE))
  uses_time <- vapply(rownames(factors),
                      function(r) survival_variable %in% all.vars(str2lang(r)), logical(1))
  uses_time[1] <- FALSE  # the response row
  if (!any(uses_time)) return(invisible(TRUE))
  time_terms <- colnames(factors)[colSums(factors[uses_time, , drop = FALSE] > 0) > 0]
  bad <- setdiff(time_terms, survival_variable)
  if (length(bad) > 0) {
    stop(sprintf(paste0(
      "`form_conditional_cr` uses the survival time `%s` in %s. It may only enter as a ",
      "plain main effect (e.g. status ~ %s + age), since prediction evaluates it ",
      "linearly at each point of the integration grid."),
      survival_variable, paste(sprintf("`%s`", bad), collapse = ", "), survival_variable), call. = FALSE)
  }
  invisible(TRUE)
}
