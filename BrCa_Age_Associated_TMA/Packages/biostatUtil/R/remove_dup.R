#' Remove duplicates
#'
#' Remove duplicates for specified columns of a data frame
#'
#' In Mass Spec data, there are occasionally duplicated entries that need to be
#' removed before further analysis. Duplication is indicated by the
#' `Quan.Info` and `PSM.Ambiguity` columns. `remove_dup` removes
#' duplicates for certain columns, then collapses repeated information into a
#' single row.
#'
#' This function is intended to be used after a call to
#' [dplyr::group_by()] such that the removal of duplicates is
#' performed within each group of unique protein IDs (e.g.
#' `Reporter.Quan.Result.ID`).
#'
#' @param x data frame
#' @param cols character vector of column names from `x` to remove
#'   duplicates
#' @return A data frame with potentially fewer rows than `x` after
#'   duplicated entries have been removed and repeated information has been
#'   collapsed.
#' @author Derek Chiu
#' @export
remove_dup <- function(x, cols) {
  idxp <- toupper(x$Quan.Info) == "UNIQUE" | toupper(x$PSM.Ambiguity) == "SELECTED"
  # If no duplicates, take only the first element
  if (any(idxp)) {
    odf <- x[which(idxp)[1], ]
  } else {
    # Otherwise take first element after collapsing info from non-missing rows
    odf <- mutate_at(x, vars(one_of(cols)),
                     funs(collapse_var(.[. != "" & !is.na(.)])), collapse = " | ")[1, ]
  }
  return(odf)
}
