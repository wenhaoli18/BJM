#' Build a survival-time transform basis from cut points
#'
#' @description \code{dynamicPrediction()}, \code{dynamicPredictionBio()},
#' \code{predictPlot()}, and \code{riskPlot()} all take a pair of arguments,
#' \code{survival_variable_all}/\code{survival_trans_function}, that describe
#' transformed basis variables of the (remaining) survival time; these can
#' optionally be referenced in \code{long_sub_fixed} formulas to let the
#' longitudinal sub-model depend flexibly on time-to-event. In every example
#' in this package, that pair follows the same convention: variables named
#' \code{"Tyears1"}, \code{"Tyears2"}, ... , each defined as the absolute
#' distance from a fixed cut point (\code{function(x) abs(x - k)}).
#' \code{survivalTrans()} builds exactly that pair from a plain vector of cut
#' points, so you do not have to hand-write two parallel lists of matching
#' names and closures. You remain free to construct
#' \code{survival_variable_all}/\code{survival_trans_function} by hand for any
#' other transform.
#'
#' @param cut_points A non-empty numeric vector of cut points, one per
#' transformed basis variable.
#' @param prefix Prefix used for the generated variable names in
#' \code{survival_variable_all} (\code{"Tyears"} by default, giving
#' \code{"Tyears1"}, \code{"Tyears2"}, ...).
#'
#' @return A named list with elements \code{survival_variable_all} and
#' \code{survival_trans_function}, in the format expected by
#' \code{dynamicPrediction()}, \code{dynamicPredictionBio()},
#' \code{predictPlot()}, and \code{riskPlot()}.
#'
#' @examples
#' trans <- survivalTrans(c(1, 3, 5, 7))
#' trans$survival_variable_all
#' trans$survival_trans_function[[1]](2)
#'
#' # equivalent to hand-writing:
#' survival_variable_all <- list("Tyears1", "Tyears2", "Tyears3", "Tyears4")
#' survival_trans_function <- list(
#'   fun1 = function(x) abs(x - 1),
#'   fun2 = function(x) abs(x - 3),
#'   fun3 = function(x) abs(x - 5),
#'   fun4 = function(x) abs(x - 7)
#' )
#'
#' @export
survivalTrans <- function(cut_points, prefix = "Tyears") {
  assert_string(prefix, "prefix")
  if (!is.numeric(cut_points) || length(cut_points) == 0 || anyNA(cut_points)) {
    stop("`cut_points` must be a non-empty numeric vector with no missing values.", call. = FALSE)
  }

  survival_variable_all <- as.list(paste0(prefix, seq_along(cut_points)))
  survival_trans_function <- stats::setNames(
    lapply(cut_points, function(k) {
      force(k)
      function(x) abs(x - k)
    }),
    paste0("fun", seq_along(cut_points))
  )

  list(survival_variable_all = survival_variable_all,
       survival_trans_function = survival_trans_function)
}
