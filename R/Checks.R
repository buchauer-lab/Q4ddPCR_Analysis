#' several functions that check the correctness of data frames, vectors, input,
#' and output format etc.

#' Check if is named vector
#'
#' Check if an object is a named vector
#' @param obj The object to check
#' @return Boolean
is_named_vector <- function(obj) {
  # Check if the object is an atomic vector
  is_atomic_vector <- is.atomic(obj) && !is.null(obj)

  # Check if the object has non-empty names
  has_names <- !is.null(names(obj)) && all(nzchar(names(obj)))

  # Return TRUE only if both conditions are satisfied
  is_atomic_vector && has_names
}

#' Check if df is dtQC
#'
#' Check if data frame object fulfills requirements for dtQC.
#' @param dtQC The data frame to be tested.
#' @return Boolean
is_dtQC <- function(dtQC) {
  # check if dtQC is a data frame
  if (!(is.data.frame(dtQC))) {
    return(FALSE)
  }
  # check if all required columns are present
  if (!all(c(
    "Well", "Sample description 1", "DyeName(s)", "Target",
    "Conc(copies/uL)", "Accepted Droplets", "Positives",
    "Negatives", "Threshold", "Total positives"
  ) %in% colnames(dtQC))) {
    return(FALSE)
  }
  return(TRUE)
}

#' Check if is grouped data
#'
#' Check if an object is a grouped data df
#' @param grouped_data The object to check
#' @return Boolean
is_grouped_data <- function(grouped_data) {
  # check if it is a data frame
  if (!is.data.frame(grouped_data)) {
    return(FALSE)
  }
  # check if all required columns are present
  if (!all(c("Sample description 1", "Dilution Factor", "shared_Wells") %in% colnames(dtQC))) {
    return(FALSE)
  }
  return(TRUE)
}

# is_household_table <- function(){
#
# }

# is_output_format <- function(){
#
# }
