#' Name cutpoint variables
#'
#' Create a name for cutpoint variables based on location of cut and number of
#' groups.
#'
#' The naming system is based on how [Hmisc::cut2()] cuts variables. Used in
#' [build_cuts()] for naming the new cutpoint variables.
#'
#' @param x vector to split by certain cutpoints
#' @param cuts where to cut `x`
#'
#' @return A character string representing the cutpoint variable name.
#' @author Derek Chiu
#' @seealso [Hmisc::cut2()], [build_cuts()]
#' @export
#'
#' @examples
#' set.seed(1108)
#' x <- sample(0:4, size = 1000, replace = TRUE)
#' name_cuts(x, c(1, 4))
#' name_cuts(x, c(2, 4))
name_cuts <- function(x, cuts) {
  . <- NULL
  levs <- sort(unique(x))
  n <- switch(as.character(length(cuts)),
              "2" = "b", "3" = "t", "4" = "qd", "5" = "qn")
  v.ind <- match(cuts[-length(cuts)], levs)
  cut.name <- rep("", length(levs)) %>%
    magrittr::inset(v.ind, "v") %>%
    paste0(levs) %>%
    paste(collapse = "") %>%
    paste0(n, .)
  return(cut.name)
}
