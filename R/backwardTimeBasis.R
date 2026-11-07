#' Add a backward-time basis to longitudinal sub-model formulas
#'
#' @description The backward joint model describes a biomarker given the
#' event time \eqn{T}. \code{\link{survivalTrans}()} adds flexibility in
#' \eqn{T} itself, through terms \eqn{|T - c|} that are the same at every
#' visit of a patient. Biomarkers often change instead with the time left
#' until the event, the \emph{backward time} \eqn{u = T - s}, where \eqn{s}
#' is the visit time: for example, a marker that stays flat for years and
#' then accelerates in the last year or two before the event. A term in
#' \eqn{T}, \eqn{s} and \eqn{T \times s} can only describe such a trajectory
#' as linear in \eqn{u}.
#'
#' \code{backwardTimeBasis()} adds the piecewise-linear basis
#' \eqn{|u - c| = |T - s - c|}, one term per cut point \eqn{c}, to each
#' formula. The terms are written as \code{I(abs(years - year - c))}, so they
#' are ordinary columns of the fixed-effects design. They need no
#' \code{survival_variable_all}/\code{survival_trans_function} entries:
#' \code{\link{predictRisk}()}, \code{\link{predictLongitudinal}()} and
#' \code{\link{simulateTrajectory}()} evaluate them at every integration grid
#' point of the event time from the formula itself, using each row's visit
#' time.
#'
#' @param long_sub_fixed A fixed-effects formula, or a list of them (one per
#'   biomarker, as passed to \code{\link{longitudinalSub}()}). Each formula
#'   should already use the event time, e.g.
#'   \code{serBilir ~ year + age + sex + years}.
#' @param cut_points A non-empty numeric vector of cut points on the
#'   backward-time scale, in the units of \code{survival_variable}. Knots
#'   where the trajectory may bend, e.g. \code{c(0.5, 1, 2, 4)} for changes
#'   concentrated in the last few years before the event.
#' @param survival_variable Name of the event-time column (\code{"years"} in
#'   \code{pbc3}).
#' @param time_variable Name of the visit-time column (\code{"year"} in
#'   \code{pbc3}).
#'
#' @return \code{long_sub_fixed} with the terms
#'   \code{I(abs(survival_variable - time_variable - c))} appended to the
#'   right-hand side of each formula. A single formula in gives a single
#'   formula out, and a list keeps its names. Each formula keeps its
#'   environment.
#'
#' @seealso \code{\link{survivalTrans}()} for a basis in the event time
#'   itself; \code{vignette("nonlinear-terms", package = "BJM")}.
#'
#' @examples
#' backwardTimeBasis(serBilir ~ year + age + sex + years, c(1, 3))
#'
#' long_sub_fixed <- backwardTimeBasis(
#'   list("serBilir" = serBilir ~ year + age + sex + years,
#'        "albumin"  = albumin  ~ year + age + sex + years),
#'   cut_points = c(1, 3, 5)
#' )
#' long_sub_fixed
#'
#' \donttest{
#' data(pbc3)
#' data_fit_all <- list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])
#' long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed,
#'                                 list("serBilir" = ~ year | id, "albumin" = ~ year | id))
#' }
#' @export
backwardTimeBasis <- function(long_sub_fixed, cut_points, survival_variable = "years",
                              time_variable = "year") {
  assert_string(survival_variable, "survival_variable")
  assert_string(time_variable, "time_variable")
  if (identical(survival_variable, time_variable)) {
    stop("`survival_variable` and `time_variable` must name two different columns.", call. = FALSE)
  }
  if (!is.numeric(cut_points) || length(cut_points) == 0 || anyNA(cut_points) ||
      any(!is.finite(cut_points))) {
    stop("`cut_points` must be a non-empty numeric vector of finite values.", call. = FALSE)
  }
  if (anyDuplicated(cut_points)) {
    stop("`cut_points` must not contain duplicates (they would give identical design columns).",
         call. = FALSE)
  }

  single <- inherits(long_sub_fixed, "formula")
  formulas <- if (single) list(long_sub_fixed) else long_sub_fixed
  if (!is.list(formulas) || length(formulas) == 0 ||
      !all(vapply(formulas, inherits, logical(1), what = "formula"))) {
    stop("`long_sub_fixed` must be a formula or a non-empty list of formulas.", call. = FALSE)
  }

  T_sym <- as.name(survival_variable)
  s_sym <- as.name(time_variable)
  terms_new <- lapply(cut_points, function(k) bquote(I(abs(.(T_sym) - .(s_sym) - .(k)))))
  ### deparse1() needs R >= 4.0; the package supports R >= 3.5
  term_label <- function(tm) paste(deparse(tm, width.cutoff = 500L), collapse = " ")

  out <- lapply(seq_along(formulas), function(i) {
    f <- formulas[[i]]
    if (length(f) != 3) {
      stop(sprintf("`long_sub_fixed[[%d]]` must be a two-sided formula (biomarker ~ covariates).", i),
           call. = FALSE)
    }
    present <- vapply(terms_new, term_label, character(1)) %in%
      attr(stats::terms(f), "term.labels")
    if (any(present)) {
      stop(sprintf("`long_sub_fixed[[%d]]` already contains the backward-time term(s) %s.", i,
                   paste(vapply(terms_new[present], term_label, character(1)), collapse = ", ")),
           call. = FALSE)
    }
    ### edit the right-hand side in place, so the formula keeps its class,
    ### its environment and the order of its existing terms
    rhs <- f[[3]]
    for (tm in terms_new) rhs <- call("+", rhs, tm)
    f[[3]] <- rhs
    f
  })
  names(out) <- names(formulas)
  if (single) out[[1]] else out
}
