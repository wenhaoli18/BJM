#' Bio_i-independent shared computation for dynamicPredictionBio()
#'
#' @description Internal helper factoring out the part of
#' \code{dynamicPredictionBio()}'s pipeline that does not depend on which
#' biomarker (\code{bio_i}) is being predicted: restricting to at-risk
#' patients, building the survival-side integration grid out to
#' \code{upper_bound}, dispatching to the Gaussian-copula conditional-density
#' variants when needed (see \code{predictRisk()}), and evaluating the
#' denominator conditional density (\code{conditionalYT}/\code{conditionalYDT}).
#' \code{dynamicPredictionBioAll()} computes this once and reuses it across
#' every requested biomarker, instead of recomputing it once per biomarker
#' (including, under the copula path, its \code{mvtnorm::pmvnorm()} calls)
#' the way calling \code{dynamicPredictionBio()} once per biomarker would.
#'
#' @param data_predict_all See \code{\link{dynamicPredictionBio}}.
#' @param long_fit_all See \code{\link{dynamicPredictionBio}}.
#' @param survival_fit_all See \code{\link{dynamicPredictionBio}}.
#' @param prediction_time See \code{\link{dynamicPredictionBio}}.
#' @param time_variable See \code{\link{dynamicPredictionBio}}.
#' @param survival_variable_all See \code{\link{dynamicPredictionBio}}.
#' @param survival_trans_function See \code{\link{dynamicPredictionBio}}.
#' @param bandcount2 See \code{\link{dynamicPredictionBio}}.
#' @return A list with the pieces \code{compute_bio_marker_step()} needs:
#' \code{data_predict_all} (at-risk-filtered), \code{survival_variable},
#' \code{predict.time.infinity}, \code{S_T_all_infinity}, \code{has_cr},
#' \code{D_T_all_infinity} (\code{NULL} if no competing risk),
#' \code{f_y_D_all_infinity}, and the dispatched
#' \code{conditionalYTBio_fun}/\code{conditionalYDTBio_fun} functions.
#' @keywords internal
compute_bio_shared_step <- function(data_predict_all, long_fit_all, survival_fit_all,
                                     prediction_time, time_variable,
                                     survival_variable_all, survival_trans_function,
                                     bandcount2) {
  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2]

  ## at risk sample
  data_predict_all = subset_at_risk(data_predict_all, survival_variable, prediction_time)

  upper_bound = integration_upper_bound(data_predict_all, long_fit_all, survival_fit_all,
                                        prediction_time)

  infinity_grid <- prepare_infinity_grid(data_predict_all, long_fit_all, survival_fit_all,
                                          prediction_time, upper_bound, bandcount2)
  predict.time.infinity = infinity_grid$predict.time.infinity
  S_T_all_infinity = infinity_grid$S_T_all_infinity

  # See predictRisk() for the rationale: whenever long_fit_all was
  # fit via the Gaussian-copula path (at least one ordinal biomarker),
  # dispatch to the copula-aware conditional-density variants instead of
  # the all-continuous originals. This is a strict no-op for an
  # all-continuous fit (long_fit_all$biomarker_type is NULL there).
  use_copula = !is.null(long_fit_all$biomarker_type) && any(long_fit_all$biomarker_type == "ordinal")
  conditionalYT_fun = if (use_copula) conditionalYTCopula else conditionalYT
  conditionalYDT_fun = if (use_copula) conditionalYDTCopula else conditionalYDT
  conditionalYTBio_fun = if (use_copula) conditionalYTBioCopula else conditionalYTBio
  conditionalYDTBio_fun = if (use_copula) conditionalYDTBioCopula else conditionalYDTBio

  has_cr = length(survival_fit_all$form_conditional_cr) != 0
  D_T_all_infinity = NULL
  if (has_cr) {
    #conditional probability D|T
    D_T_all_infinity = conditionalDT(data_predict_all, long_fit_all, survival_fit_all,
                                      l_i = predict.time.infinity)
    #conditional probability Y|D,T
    f_y_D_all_infinity = conditionalYDT_fun(data_predict_all, long_fit_all, survival_fit_all,
                                             l_i = predict.time.infinity, survival_variable,
                                             time_variable, survival_variable_all,
                                             survival_trans_function)
  } else {
    #conditional probability Y|T
    f_y_D_all_infinity = conditionalYT_fun(data_predict_all, long_fit_all, l_i = predict.time.infinity,
                                            survival_variable, time_variable, survival_variable_all,
                                            survival_trans_function)
  }

  ### per-patient shift for the log densities (see patient_log_shift());
  ### the same shift is applied to every candidate value's numerator in
  ### compute_bio_marker_step(), so it cancels in the predicted density
  log_shift = if (has_cr) patient_log_shift(f_y_D_all_infinity[[1]], f_y_D_all_infinity[[2]]) else
    patient_log_shift(f_y_D_all_infinity[[1]])

  list(data_predict_all = data_predict_all, survival_variable = survival_variable,
       predict.time.infinity = predict.time.infinity, S_T_all_infinity = S_T_all_infinity,
       has_cr = has_cr, D_T_all_infinity = D_T_all_infinity,
       f_y_D_all_infinity = f_y_D_all_infinity, log_shift = log_shift,
       conditionalYTBio_fun = conditionalYTBio_fun, conditionalYDTBio_fun = conditionalYDTBio_fun)
}

#' Bio_i-specific computation for dynamicPredictionBio()
#'
#' @description Internal helper for the part of \code{dynamicPredictionBio()}'s
#' pipeline that does depend on which biomarker is being predicted: building
#' the candidate-value grid \code{Y_all}, evaluating the numerator
#' conditional density on it, and assembling/normalizing \code{Y_density}
#' and \code{Y_predict}. Takes the output of \code{compute_bio_shared_step()}
#' so the bio_i-independent pieces are not recomputed.
#'
#' @param shared Output of \code{compute_bio_shared_step()}.
#' @param bio_i See \code{\link{dynamicPredictionBio}}.
#' @param long_fit_all See \code{\link{dynamicPredictionBio}}.
#' @param survival_fit_all See \code{\link{dynamicPredictionBio}}.
#' @param prediction_time See \code{\link{dynamicPredictionBio}}.
#' @param horizon See \code{\link{dynamicPredictionBio}}.
#' @param time_variable See \code{\link{dynamicPredictionBio}}.
#' @param survival_variable_all See \code{\link{dynamicPredictionBio}}.
#' @param survival_trans_function See \code{\link{dynamicPredictionBio}}.
#' @param bandcount3 See \code{\link{dynamicPredictionBio}}.
#' @return A list with \code{Y_predict}, \code{Y_density}, \code{Y_all}
#' (matching \code{dynamicPredictionBio()}'s return value fields). For an
#' \strong{ordinal} \code{bio_i}, \code{Y_all}/\code{Y_predict} hold integer
#' category codes (\code{1:K}, in threshold order) rather than a numeric
#' grid -- see Details below -- and \code{Y_all} additionally carries a
#' \code{"category_labels"} attribute with the matching level-label strings.
#' @keywords internal
compute_bio_marker_step <- function(shared, bio_i, long_fit_all, survival_fit_all,
                                     prediction_time, horizon, time_variable,
                                     survival_variable_all, survival_trans_function, bandcount3) {
  bio_i_name <- as.character(formula(long_fit_all$long_sub_fixed[[bio_i]])[[2]])
  is_ordinal_target <- !is.null(long_fit_all$biomarker_type) &&
    long_fit_all$biomarker_type[bio_i] == "ordinal"

  if (is_ordinal_target) {
    # An ordinal bio_i's candidate grid is fixed by its category count, not
    # controlled by bandcount3 at all (see conditionalYTBioCopula()).
    # Y_query holds the actual level-label strings needed to assign into the
    # candidate row's factor column (see
    # select_patient_longitudinal_data_bio()/conditionalYTBioCopula()), while
    # Y_all -- what is exposed to the caller -- holds the matching integer
    # category codes (1:K, in threshold order), so that downstream
    # numeric-only consumers (.format_dynamicPredictionBio(),
    # max_relative_diff(), predictPlot()) need no changes.
    Y_labels <- levels(shared$data_predict_all[[bio_i]][[bio_i_name]])
    Y_query <- Y_labels
    Y_all <- seq_along(Y_labels)
  }

  ### Y_density at the candidate values Y_query (Y_all holds the matching
  ### numeric values), one row per candidate value, one column per patient
  density_on_grid <- function(Y_query, Y_all) {
    if (shared$has_cr) {
      f_y_D_all_predict = shared$conditionalYDTBio_fun(Y_query, time_new = prediction_time + horizon,
                                            bio_i, shared$data_predict_all, long_fit_all,
                                            survival_fit_all,
                                            l_i = shared$predict.time.infinity, shared$survival_variable,
                                            time_variable, survival_variable_all,
                                            survival_trans_function)

      # Y_density is filled row-by-row into a pre-allocated matrix rather than
      # grown with rbind() inside the loop: rbind()-in-a-loop reallocates and
      # copies the whole growing matrix on every iteration (O(bandcount3^2)
      # copies total), which becomes non-negligible at the large end of
      # bandcount3's auto-tuned range. The matrix is allocated on the first
      # iteration once the per-patient row length is known, so this makes no
      # assumption about that length elsewhere.
      ### log densities, exponentiated after the per-patient shift computed
      ### in compute_bio_shared_step(); the denominator does not depend on
      ### the candidate value, so it is computed once
      T.surv.infinity.0 = t(exp_shifted(shared$f_y_D_all_infinity[[1]], shared$log_shift) *
                              shared$D_T_all_infinity[[1]] * shared$S_T_all_infinity)
      T.surv.infinity.1 = t(exp_shifted(shared$f_y_D_all_infinity[[2]], shared$log_shift) *
                              shared$D_T_all_infinity[[2]] * shared$S_T_all_infinity)
      denominator = rowSums(T.surv.infinity.1 + T.surv.infinity.0)
      Y_density = NULL
      for (Y_i in seq_len(length(Y_all))) {
        T.surv.predict.0 = t(exp_shifted(f_y_D_all_predict[[1]][[Y_i]], shared$log_shift) *
                               shared$D_T_all_infinity[[1]] * shared$S_T_all_infinity)
        T.surv.predict.1 = t(exp_shifted(f_y_D_all_predict[[2]][[Y_i]], shared$log_shift) *
                               shared$D_T_all_infinity[[2]] * shared$S_T_all_infinity)

        ### a density, not a probability: it may exceed 1 (a narrow predictive
        ### distribution), so it must not go through clamp_risk_prob()
        Y_density_row = pmax(rowSums(T.surv.predict.1) / denominator, 0) +
          pmax(rowSums(T.surv.predict.0) / denominator, 0)
        if (is.null(Y_density)) Y_density = matrix(NA_real_, length(Y_all), length(Y_density_row))
        Y_density[Y_i, ] = Y_density_row
      }
    } else { #without competing risk

      f_y_D_all_predict = shared$conditionalYTBio_fun(Y_query, time_new = prediction_time + horizon,
                                            bio_i, shared$data_predict_all, long_fit_all,
                                            l_i = shared$predict.time.infinity, shared$survival_variable,
                                            time_variable, survival_variable_all,
                                            survival_trans_function)

      # See the competing-risk branch above for why Y_density is filled into a
      # pre-allocated matrix instead of grown with rbind() in the loop.
      ### shifted log densities, as in the competing-risk branch; this also
      ### replaces a "+ 1e-20" in the denominator, which swamped it whenever
      ### the densities were small and drove the predicted density towards 0
      denominator = rowSums(t(exp_shifted(shared$f_y_D_all_infinity[[1]], shared$log_shift) *
                                shared$S_T_all_infinity))
      Y_density = NULL
      for (Y_i in seq_len(length(Y_all))) {

        T.surv.predict.0 = t(exp_shifted(f_y_D_all_predict[[1]][[Y_i]], shared$log_shift) *
                               shared$S_T_all_infinity)

        ### a density, not a probability: see the competing-risk branch
        Y_density_row = pmax(rowSums(T.surv.predict.0) / denominator, 0)
        if (is.null(Y_density)) Y_density = matrix(NA_real_, length(Y_all), length(Y_density_row))
        Y_density[Y_i, ] = Y_density_row

      }
    }
    Y_density
  }

  if (is_ordinal_target) {
    Y_density <- density_on_grid(Y_query, Y_all)
  } else {
    Y_all <- continuous_value_grid(shared, bio_i, bio_i_name, long_fit_all, bandcount3,
                                   density_fun = function(y) density_on_grid(y, y))
    Y_density <- density_on_grid(Y_all, Y_all)
  }

  Y_predict = vapply(seq_len(dim(Y_density)[2]), function(i)
    density_mode(Y_all, Y_density[, i], refine = !is_ordinal_target), numeric(1))

  ### name each patient's prediction by id (see prediction_patient_ids())
  patient_ids <- prediction_patient_ids(shared$data_predict_all, long_fit_all)
  if (length(patient_ids) == length(Y_predict)) {
    names(Y_predict) <- patient_ids
    colnames(Y_density) <- patient_ids
  }

  if (is_ordinal_target) attr(Y_all, "category_labels") <- Y_labels

  list(Y_predict = Y_predict, Y_density = Y_density, Y_all = Y_all)
}

#' Predict a future biomarker value from fitted sub-models, for a single
#' biomarker
#'
#' @description
#' \strong{Internal single-biomarker engine} behind
#' \code{\link{predictLongitudinal}} -- call \code{predictLongitudinal()}
#' directly instead (it dispatches here automatically when \code{bio_i}
#' names exactly one biomarker, and to \code{\link{dynamicPredictionBioAll}}
#' otherwise). Kept as a separate internal function -- rather than folded
#' into \code{predictLongitudinal()} -- because \code{\link{predictPlot}}
#' and \code{\link{checkBandcountConvergence}} also call it directly for a
#' single biomarker at a time.
#'
#' Companion to \code{\link{predictRisk}}, using the same fitted
#' longitudinal (\code{\link{longitudinalSub}}) and survival
#' (\code{\link{survivalSub}}) sub-models, but instead of an event-risk
#' probability this returns a predictive density for a future value of one
#' chosen biomarker (\code{bio_i}) at \code{prediction_time + horizon},
#' conditional on the subject's observed longitudinal history up to
#' \code{prediction_time} and on being event-free at \code{prediction_time}.
#' The density (\code{Y_density}, evaluated over a candidate-value grid
#' \code{Y_all}) is obtained by integrating the biomarker's predictive
#' distribution against the survival sub-model's hazard, using the subject's
#' empirical-Bayes random-effects update from their observed history; its
#' mode (\code{Y_predict}) is reported as the point prediction. As in
#' \code{\link{predictRisk}}, the survival-side integrals are
#' evaluated on numerical grids controlled by \code{bandcount2}, and the
#' density itself is evaluated on a grid controlled by \code{bandcount3};
#' see Details.
#'
#' \code{bio_i} may refer to either a \strong{continuous} or an
#' \strong{ordinal} biomarker (see \code{\link{longitudinalSubCopula}});
#' for an ordinal \code{bio_i}, \code{Y_all} and \code{Y_predict} hold
#' integer \emph{category codes} (\code{1:K}, in the fitted factor's
#' \code{levels()}/threshold order) rather than a numeric grid, with
#' \code{Y_all} additionally carrying a \code{"category_labels"} attribute
#' giving the matching level-label strings; \code{bandcount3} is ignored in
#' that case, since the candidate grid is fixed at the biomarker's category
#' count (see Details).
#'
#' The prediction is conditional on the longitudinal history observed up to
#' \code{prediction_time}: rows of \code{data_predict_all} whose
#' \code{time_variable} is later than \code{prediction_time} are dropped,
#' with a warning, before predicting.
#'
#' @param bio_i Biomarker used to do prediction. May be continuous or
#' ordinal (see Details).
#' @param data_predict_all This involves a collection of \code{data.frame} objects for
#' dynamic prediction, each corresponding to a distinct longitudinal outcome.
#' These data frames should contain the variables specified in \code{long_sub_fixed}
#' and \code{long_sub_random}. Utilizing a list structure
#' allows for the incorporation of multiple longitudinal outcomes,
#' each potentially following different measurement protocols.
#' In instances where all longitudinal outcomes are recorded at identical
#' time points across patients, a singular \code{data.frame} object may
#' be used in a \code{list}. Alternatively, a single bare \code{data.frame}
#' (not wrapped in a list) may be supplied directly; it is then reused for
#' every longitudinal outcome. It is presumed that each data frame is
#' structured in a long format.
#'
#' @param long_fit_all Outputs from the model fitting process using the \code{nlme} package,
#' encompassing the results and parameters obtained from the analysis.
#' @param survival_fit_all Results and parameters generated from the model fitting
#' procedure, utilizing the \code{coxph} function. These outputs include the comprehensive
#' findings and variables derived from the analysis.
#' @param prediction_time Time used to make the prediction
#' @param horizon Prediction horizon
#' @param time_variable The name of time variable in linear mixed model.
#' @param survival_variable_all The name of the transformed time-to-event outcomes variable.
#' @param survival_trans_function The transformation function used for time-to-event outcomes,
#' in the order of \code{survival_variable_all}.
#' @param bandcount2 The number of grid points spanning
#' \code{[prediction_time, upper_bound]}, where \code{upper_bound} is set
#' internally as the earliest time by which every at-risk patient's
#' model-based probability of still being event-free (given event-free at
#' \code{prediction_time}) has dropped below \code{1e-4}; this approximates
#' integrating out to infinity for the denominator that normalizes the predicted density. A wider follow-up range
#' needs a larger \code{bandcount2} to keep the grid spacing comparable.
#' Defaults to \code{"auto"} (see Details).
#' @param bandcount3 The number of points in the candidate-biomarker-value
#' grid (\code{Y_all}) used to build the predicted density
#' (\code{Y_density}) and locate its mode (\code{Y_predict}). This controls
#' the resolution of the density curve, not a time integral; increase it if
#' the density looks jagged or \code{Y_predict} jumps erratically between
#' nearby grid points. Defaults to \code{"auto"} (see Details). Ignored when
#' \code{bio_i} is ordinal (see Details) -- that candidate grid is always
#' the biomarker's fixed category count.
#'
#' @details
#' There is no universal correct value for \code{bandcount2}/\code{bandcount3}:
#' as a practical check, double both and confirm the results barely change;
#' if they do, keep doubling. By default (\code{bandcount2 = "auto"},
#' \code{bandcount3 = "auto"}), this doubling check is done for you:
#' starting from small built-in values, both are doubled together, and the
#' result is compared to the previous round, until the largest relative
#' change in \code{Y_predict} drops below 1%, or 2 doublings have been
#' tried (so at most 3 calls' worth of work). If it still has not converged
#' by then, a warning reports this and the result at the largest value
#' tried is returned anyway (not an error), so this never silently loops
#' for an unbounded amount of time. Pass an explicit number for either
#' argument to skip auto-tuning it and use a fixed value instead (as in
#' previous package versions), or call \code{checkBandcountConvergence()}
#' directly for more control over the tolerance and doubling count. See
#' also \code{vignette("BJM-intro", package = "BJM")} for a worked example.
#'
#' @return An object of class \code{"dynamicPredictionBio.BJM"}, a named list with elements:
#' \describe{
#'   \item{Y_predict}{A vector, one entry per at-risk patient (named by patient id), giving the MAP (most likely) predicted
#'   value of biomarker \code{bio_i} at \code{prediction_time + horizon}. For an ordinal
#'   \code{bio_i}, this is an integer category code (see Details), not a raw value.}
#'   \item{Y_density}{A probability matrix whose rows correspond to the candidate biomarker
#'   values in \code{Y_all} and whose columns correspond to individual patients; each entry
#'   is the dynamically predicted density of the biomarker taking that value. For an ordinal
#'   \code{bio_i}, each row is instead that category's predicted probability.}
#'   \item{Y_all}{The grid of candidate biomarker values used to build \code{Y_density}. For an
#'   ordinal \code{bio_i}, this is the vector of integer category codes \code{1:K}, carrying a
#'   \code{"category_labels"} attribute with the matching level-label strings (see Details).}
#' }
#'
#' @keywords internal
dynamicPredictionBio = function(bio_i, data_predict_all, long_fit_all, survival_fit_all,
                                 prediction_time, horizon, time_variable,
                                 survival_variable_all, survival_trans_function,
                                 bandcount2 = "auto", bandcount3 = "auto"){

  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  assert_index(bio_i, length(long_fit_all$lfit), "bio_i", "longitudinal outcomes in long_fit_all")
  assert_data_list(data_predict_all, "data_predict_all", length(long_fit_all$lfit), allow_bare_df = TRUE)
  if (!is.list(data_predict_all) || is.data.frame(data_predict_all)) {
    data_predict_all <- rep(list(data_predict_all), each = length(long_fit_all$lfit))
  }
  assert_scalar_numeric(prediction_time, "prediction_time")
  assert_scalar_numeric(horizon, "horizon")
  assert_string(time_variable, "time_variable")
  assert_bandcount(bandcount2, "bandcount2")
  assert_bandcount(bandcount3, "bandcount3")
  assert_survival_trans(survival_variable_all, survival_trans_function, probe_value = prediction_time)
  data_predict_all <- drop_after_prediction_time(data_predict_all, time_variable, prediction_time)
  data_predict_all <- align_ordinal_levels(data_predict_all, long_fit_all)
  data_predict_all <- drop_missing_longitudinal(data_predict_all, long_fit_all,
                                               outcome_variables(survival_fit_all),
                                               survival_variable_all)

  # bandcount2/bandcount3 = "auto" (the default): see predictRisk()
  # for the rationale; same doubling-until-stable check, applied here to
  # Y_predict instead of the risk probabilities. `call_args` is captured
  # before `auto_names` is computed, so it is not swept up into `call_args`
  # too (see predictRisk() for the same pattern). For an ordinal
  # bio_i, the candidate grid is fixed at its category count and is not
  # controlled by bandcount3 at all (see compute_bio_marker_step()), so
  # bandcount3 is excluded from auto-tuning in that case even when left at
  # its "auto" default.
  call_args <- as.list(environment())
  bio_i_is_ordinal <- !is.null(long_fit_all$biomarker_type) &&
    long_fit_all$biomarker_type[bio_i] == "ordinal"
  auto_names <- c("bandcount2", "bandcount3")[c(identical(bandcount2, "auto"),
                                                  identical(bandcount3, "auto") && !bio_i_is_ordinal)]
  if (length(auto_names) > 0) {
    return(auto_tune_bandcount(dynamicPredictionBio, call_args, auto_names)$result)
  }

  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2] #survival_variable = "fuyrs"
  for (i in seq_along(data_predict_all)) {
    assert_vars_in_data(time_variable, data_predict_all[[i]],
                         "time_variable", sprintf("data_predict_all[[%d]]", i))
    assert_vars_in_data(survival_variable, data_predict_all[[i]],
                         "the survival-time variable used to fit survival_fit_all",
                         sprintf("data_predict_all[[%d]]", i))
  }

  # See compute_bio_shared_step()/compute_bio_marker_step() for what each
  # half computes; dynamicPredictionBioAll() calls the same two helpers
  # directly so the shared step can be reused across several biomarkers
  # instead of being recomputed once per call the way this does.
  shared <- compute_bio_shared_step(data_predict_all, long_fit_all, survival_fit_all,
                                     prediction_time, time_variable, survival_variable_all,
                                     survival_trans_function, bandcount2)
  marker <- compute_bio_marker_step(shared, bio_i, long_fit_all, survival_fit_all,
                                     prediction_time, horizon, time_variable,
                                     survival_variable_all, survival_trans_function, bandcount3)

  out <- list(Y_predict = marker$Y_predict, Y_density = marker$Y_density, Y_all = marker$Y_all)
  class(out) <- "dynamicPredictionBio.BJM"
  return(out)
}

#' Candidate-value grid for a continuous biomarker
#'
#' @description Helper for \code{compute_bio_marker_step()}: the grid of
#' candidate values (\code{Y_all}) a continuous biomarker's predictive
#' density is tabulated on, shared by all patients predicted. Built in two
#' passes:
#' \enumerate{
#'   \item A coarse grid over a deliberately wide range -- 5 spans either
#'   side of the patients' observed values, the span being their range but
#'   at least the biomarker's SD in the training data -- with a step of the
#'   biomarker's residual SD \eqn{\sigma}, at least 30 and at most 100
#'   points. A predicted future measurement includes its measurement error,
#'   so every predictive density has an SD of at least \eqn{\sigma}: within
#'   half a step of its mode it is still above 88\% of its peak, far above
#'   \code{rel_tol}, so no patient's density falls between coarse points.
#'   \item \code{bandcount3} + 1 equally spaced points over the part of the
#'   coarse grid where some patient's density exceeds \code{rel_tol} times
#'   that patient's maximum, widened by \eqn{2\sigma} (and at least one
#'   coarse step) on either side.
#' }
#' The wide range used to be the final grid, so most of the
#' \code{bandcount3} points (about 90\% in a simulated example) lay where
#' every density was practically 0. If the coarse pass finds no usable
#' density, the wide range is returned with \code{bandcount3} + 1 points,
#' as before.
#'
#' @param shared Output of \code{compute_bio_shared_step()}.
#' @param bio_i Index of the biomarker.
#' @param bio_i_name Its response variable name.
#' @param long_fit_all Output of \code{longitudinalSub()}.
#' @param bandcount3 Number of intervals of the final grid.
#' @param density_fun Function of a vector of candidate values returning the
#'   density matrix (one row per value, one column per patient).
#' @param rel_tol Relative density below which a value is outside the grid.
#' @return A numeric vector of candidate values.
#' @keywords internal
continuous_value_grid <- function(shared, bio_i, bio_i_name, long_fit_all, bandcount3,
                                  density_fun, rel_tol = 1e-6) {
  observed <- shared$data_predict_all[[bio_i]][[bio_i_name]]
  Y_upper <- max(observed, na.rm = TRUE)
  Y_lower <- min(observed, na.rm = TRUE)
  ### with a single observation (or all values equal) the range is 0 and the
  ### grid would collapse to one point, hence at least the training SD
  span <- max(Y_upper - Y_lower,
              stats::sd(nlme::getResponse(long_fit_all$lfit[[bio_i]]), na.rm = TRUE),
              na.rm = TRUE)
  if (!is.finite(span) || span <= 0) span <- max(abs(Y_upper), 1)
  wide <- c(Y_lower - 5 * span, Y_upper + 5 * span)
  wide_grid <- seq(wide[1], wide[2], length.out = bandcount3 + 1)

  sigma <- long_fit_all$lfit[[bio_i]]$sigma
  if (is.null(sigma) || !is.finite(sigma) || sigma <= 0) sigma <- span / 10
  n_coarse <- min(max(ceiling(diff(wide) / sigma), 30), 100)
  coarse <- seq(wide[1], wide[2], length.out = n_coarse)
  step <- coarse[2] - coarse[1]

  d0 <- as.matrix(density_fun(coarse))
  inside <- apply(d0, 2, function(d) {
    top <- suppressWarnings(max(d, na.rm = TRUE))
    if (!is.finite(top) || top <= 0) return(rep(FALSE, length(d)))
    !is.na(d) & d > rel_tol * top
  })
  inside <- matrix(inside, nrow = length(coarse))
  if (!any(inside)) return(wide_grid)

  pad <- max(2 * sigma, step)
  range_in <- range(coarse[rowSums(inside) > 0]) + c(-pad, pad)
  seq(range_in[1], range_in[2], length.out = bandcount3 + 1)
}

#' Mode of a density tabulated on a grid
#'
#' @description Helper for \code{compute_bio_marker_step()}: the grid point
#' where \code{density} peaks, refined (when \code{refine = TRUE}) by fitting
#' a parabola through the log-density at that point and its two neighbours.
#' For a smooth, near-normal predictive density this recovers the mode to
#' well within one grid step, instead of snapping to the grid: the
#' unrefined mode moved in whole grid steps as \code{bandcount3} changed,
#' which for a predicted value near 0 is a large relative change, so
#' \code{"auto"} tuning of \code{bandcount3} reported non-convergence even on
#' a fine grid. No refinement at the grid's ends, for an ordinal biomarker
#' (whose grid is its categories), or if a neighbour's density is 0.
#'
#' @param Y_all The (equally spaced, for continuous biomarkers) grid.
#' @param density The density at each grid point.
#' @param refine Whether to refine between grid points.
#' @return A single number, or \code{NA} if \code{density} has no finite maximum.
#' @keywords internal
density_mode <- function(Y_all, density, refine = TRUE) {
  if (!any(is.finite(density))) return(NA_real_)
  k <- which.max(density)
  if (!refine || k == 1 || k == length(density)) return(Y_all[k])
  y <- density[(k - 1):(k + 1)]
  if (any(!is.finite(y)) || any(y <= 0)) return(Y_all[k])
  ly <- log(y)
  curvature <- ly[1] - 2 * ly[2] + ly[3]
  if (curvature >= 0) return(Y_all[k])
  offset <- 0.5 * (ly[1] - ly[3]) / curvature
  Y_all[k] + offset * (Y_all[k + 1] - Y_all[k])
}

