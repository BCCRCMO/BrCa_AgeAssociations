#' Sample count condition for MS analyses
#'
#' The number of non-missing samples per treatment group determine whether it
#' can be used in subsequent analyses.
#'
#' The Variance Condition (VC) (default) requires at least two non-missing
#' samples in each of the treatment groups. A rudimentary estimate of variance
#' is possible when this phenomenon exists. The Minimum Reasonable Condition
#' (MRC) requires all treatment groups to have at least one non-missing sample
#' and one treatment group to have at least two non-missing samples. Equal
#' variance assumption is used in the case of MRC. Finally, the Minimum
#' Condition (MC) requires at least one non-missing sample in each of the
#' treatment groups. MC is not generally suggested in practice.
#'
#' Additionally, if there are controls, they are handled by MC.
#'
#' @param data The data matrix with rows as observations, columns as variables
#' @param treatment A character vector indicating how different treatment groups
#'   are coded. Regular expression is used on the column names to find all
#'   relevant treatment variables.
#' @param control Optionally specify a character string indicating the control
#'   group if it exists. Regular expression is used on the column names to find
#'   all relevant control variables.
#' @param condition Either one of "VC" (default), "MRC", or "MC". See details.
#' @return A logical vector of length equal to the number of rows in
#'   `data`, indicating which observations are kept and which are to be
#'   removed.
#' @family Mass Spectrometry functions
#' @author Derek Chiu
#' @export
#' @examples
#' # Example data matrix with some missing sample values
#' g <- c("X5", "X19")
#' r <- c("1", "2", "3")
#' col.nm <- sort(apply(expand.grid(g, r), 1, paste, collapse = "_"))
#' A <- matrix(c(seq_len(36)), ncol = 6, byrow = FALSE, dimnames = list(NULL, col.nm))
#' A[c(1, 21, 26, 7, 17, 8, 23, 4, 16, 19, 5, 25, 12, 18, 24)] <- NA
#'
#' # Two treatment experimental design
#' ms_condition(A, treatment = c("X19", "X5"), condition = "VC")
#' ms_condition(A, treatment = c("X19", "X5"), condition = "MRC")
#' ms_condition(A, treatment = c("X19", "X5"), condition = "MC")
#'
#' # Control ends in "_1", Treatments ends in "_2", "_3"
#' ms_condition(A, treatment = c("_2", "_3"), control = "_1", condition = "VC")
#' ms_condition(A, treatment = c("_2", "_3"), control = "_1", condition = "MRC")
#' ms_condition(A, treatment = c("_2", "_3"), control = "_1", condition = "MC")
ms_condition <- function(data, treatment, control = NULL,
                         condition = c("VC", "MRC", "MC")) {
  # Specify condition for non-missing treatments
  ind.trt <- ms_condition_indicator(data, treatment, condition = condition)
  if (!is.null(control)) {
    # At least 1 control sample
    ind.ctrl <- ms_condition_indicator(data, control, condition = "MC")
    return(ind.ctrl & ind.trt)
  } else {
    return(ind.trt)
  }
}

#' Returns the logical vector based on the condition
#' @noRd
ms_condition_indicator <- function(data, cols,
                                   condition = c("VC", "MRC", "MC")) {
  cond <- match.arg(condition)
  dat <- vapply(cols, n_notNA, data = data, double(nrow(data)))
  ind <- apply(dat, 1, ms_condition_choice, cond = cond)
  return(ind)
}

#' Logical predicates for each condition choice
#' @noRd
ms_condition_choice <- function(x, cond) {
  switch(cond,
         VC = all(x >= 2),
         MRC = all(x >= 1) & any(x >= 2),
         MC = all(x >= 1))
}

#' Number of elements not NA in each row of data for columns in cols
#' @noRd
n_notNA <- function(data, cols) {
  dat <- data[, grep(cols, colnames(data))]
  return(rowSums(!is.na(dat)))
}
