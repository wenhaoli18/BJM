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

  upper_bound = 2 * max(data_predict_all[[1]][survival_variable])

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

  list(data_predict_all = data_predict_all, survival_variable = survival_variable,
       predict.time.infinity = predict.time.infinity, S_T_all_infinity = S_T_all_infinity,
       has_cr = has_cr, D_T_all_infinity = D_T_all_infinity,
       f_y_D_all_infinity = f_y_D_all_infinity,
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
  } else {
    Y_upper = max(shared$data_predict_all[[bio_i]][bio_i_name], na.rm = TRUE)
    Y_lower = min(shared$data_predict_all[[bio_i]][bio_i_name], na.rm = TRUE)
    Y_all = seq(Y_lower - 5 * (Y_upper - Y_lower), Y_upper + 5 * (Y_upper - Y_lower),
                11 * (Y_upper - Y_lower) / bandcount3) #seq(0, 100, 5)
    Y_query <- Y_all
  }

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
    Y_density = NULL
    for (Y_i in seq_len(length(Y_all))) {
      T.surv.predict.0 = t(f_y_D_all_predict[[1]][[Y_i]] * shared$D_T_all_infinity[[1]] * shared$S_T_all_infinity)
      T.surv.predict.1 = t(f_y_D_all_predict[[2]][[Y_i]] * shared$D_T_all_infinity[[2]] * shared$S_T_all_infinity)
      T.surv.infinity.0 = t(shared$f_y_D_all_infinity[[1]] * shared$D_T_all_infinity[[1]] * shared$S_T_all_infinity)
      T.surv.infinity.1 = t(shared$f_y_D_all_infinity[[2]] * shared$D_T_all_infinity[[2]] * shared$S_T_all_infinity)

      risk.prob.1 = clamp_risk_prob(rowSums(T.surv.predict.1), rowSums(T.surv.infinity.1 + T.surv.infinity.0))
      risk.prob.0 = clamp_risk_prob(rowSums(T.surv.predict.0), rowSums(T.surv.infinity.1 + T.surv.infinity.0))
      Y_density_row = risk.prob.1 + risk.prob.0
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
    Y_density = NULL
    for (Y_i in seq_len(length(Y_all))) {

      T.surv.predict.0 = t(f_y_D_all_predict[[1]][[Y_i]] * shared$S_T_all_infinity)
      T.surv.infinity.0 = t(shared$f_y_D_all_infinity[[1]] * shared$S_T_all_infinity)

      risk.prob.0 = clamp_risk_prob(rowSums(T.surv.predict.0), rowSums(T.surv.infinity.0 + 1e-20))
      if (is.null(Y_density)) Y_density = matrix(NA_real_, length(Y_all), length(risk.prob.0))
      Y_density[Y_i, ] = risk.prob.0

    }
  }

  Y_predict = c()
  for (i in 1:dim(Y_density)[2]) {
    if (length(Y_all[which.max(Y_density[, i])]) == 0) {
      Y_predict = c(Y_predict, NA)
    } else {
      Y_predict = c(Y_predict, Y_all[which.max(Y_density[, i])])
    }
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
#' The time values in the prediction data subset must be less than the
#' specified \code{prediction_time} which is the prediction time. The time points for
#' longitudinal repeated measurements must not surpass the prediction time.
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
#' internally to twice the longest observed survival/censoring time among
#' at-risk patients; this approximates integrating out to infinity for the
#' denominator that normalizes the predicted density. A wider follow-up range
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
#'   \item{Y_predict}{A vector, one entry per patient, giving the MAP (most likely) predicted
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
#' @examples
#'
#' \donttest{
#' data(pbc3)
#'
#' data_survival_fitting =  pbc3[!duplicated(pbc3$id), ]
#'
#' form_marginal_surv = Surv(years, status3) ~ age + sex
#' form_conditional_cr = NULL
#'
#' survival_fit_all = survivalSub(data_survival_fitting, form_marginal_surv,
#'                                form_conditional_cr)
#'
#' long_sub_fixed = list(
#'   "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,
#'   "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,
#'   "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,
#'   "long4" = alkaline ~ year + age + sex + (years) + (years) * year,
#'   "long5" = SGOT ~ year + age + sex + (years) + (years) * year,
#'   "long6" = platelets ~ year + age + sex + (years)  + (years) * year)
#'
#' long_sub_random =list(
#'   "long1" =  ~ year| id,
#'   "long2" =  ~ year| id,
#'   "long3" =  ~ year| id,
#'   "long4" =  ~ year| id,
#'   "long5" =  ~ year| id,
#'   "long6" =  ~ year| id)
#'
#' survival_variable_all = list(
#'   "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
#' )
#'
#' survival_trans_function = list(
#'   fun1 = function(x){abs(x - 1)},
#'   fun2 = function(x){abs(x - 3)},
#'   fun3 = function(x){abs(x - 5)},
#'   fun4 = function(x){abs(x - 7)}
#' )
#'
#' # Complete case analysis
#' data_fit_all = list()
#' for(i in seq_len(length(long_sub_fixed))){
#'   data_fit_all[[i]] = pbc3[pbc3$status3 == 1, ]
#' }
#'
#' # fitting longitudinal submodel
#' long_fit_all = longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)
#'
#' i_PID = 2
#' data.raw.predict.1 = pbc3[pbc3$id == i_PID, ]
#'
#' data_predict_all = list()
#' for(i in seq_len(length(long_sub_fixed))){
#'   data_predict_all[[i]] = data.raw.predict.1[data.raw.predict.1$year <= 3,]
#' }
#'
#' Y_predict = dynamicPredictionBio(bio_i = 1, data_predict_all, long_fit_all,
#'                                  survival_fit_all, prediction_time = 3,
#'                                  horizon = 3, time_variable = "year",
#'                                  survival_variable_all, survival_trans_function,
#'                                  bandcount2 = 40, bandcount3 = 400)
#'
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

