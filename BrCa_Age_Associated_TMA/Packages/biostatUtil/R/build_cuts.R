#' Build cutpoint variables
#'
#' Transforms an ordinal variable into anywhere from two to five groups for
#' cutpoint analysis.
#'
#' @param x vector to split by certain cutpoints
#' @param n either "b" for binarization (2 groups), "t" for trinarization (3
#'   groups), "q" for quads (4 groups), or "qn" for quints (5 groups)
#' @param var.prefix variable name prefix
#' @param list if `TRUE`, the variables are returned as a list.
#'
#' @return By default, a tibble of cutpoint variables built from a categorical
#'   biomarker. The number of columns correspond to all the ways the biomarker
#'   could be cut into `n` bins. Each column name starts with a "b", "t",
#'   "qd", or "qn" for "binarization", "trinarization", "quads", or "quints",
#'   respectively, with the levels being compared separated by "v". If
#'   `list = FALSE`, each cutpoint variable is an element of a list.
#'
#' @author Derek Chiu
#' @export
#'
#' @examples
#' set.seed(1108)
#' x <- sample(0:4, size = 1000, replace = TRUE)
#' build_cuts(x, n = "b")
#' build_cuts(x, n = "t")
#' build_cuts(x, n = "t", var.prefix = "PHGDH")
#' str(build_cuts(x, n = "qd", list = TRUE))
build_cuts <- function(x, n = c("b", "t", "qd", "qn"), var.prefix = NULL,
                       list = FALSE) {
  ulevs <- sort(unique(x[x > min(x, na.rm = TRUE)]))
  ng <- switch(match.arg(n), b = 2, t = 3, qd = 4, qn = 5)
  assertthat::assert_that(length(ulevs) >= ng - 1)
  res <- utils::combn(ulevs, ng - 1, simplify = FALSE) %>%
    purrr::map(c, max(x, na.rm = TRUE)) %>%
    magrittr::set_names(purrr::map_chr(., name_cuts, x = x)) %>%
    purrr::map_df(Hmisc::cut2, x = x) %>%
    magrittr::set_names(ifelse(rep(is.null(var.prefix), ncol(.)), names(.),
                               paste0(var.prefix, "_", names(.))))
  if (list) res <- as.list(res)
  return(res)
}
