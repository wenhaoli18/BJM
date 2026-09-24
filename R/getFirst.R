#' @keywords internal
getFirst <- function(input_data) {
  if (is.list(input_data)) {
    return(input_data[[1]])
  } else {
    return(input_data)
  }
}
