#' Rounding of Small Numbers
#'
#' Rounding of a number smaller than specified precision doesn't coerce to 0.
#'
#' This function is useful when we have small p-values and don't want to show
#' the scientific notation, or coercion to 0. Instead we show an upper bound.
#' For example, if a p-value is 2e-05 and we want to round to 3 digits, the
#' function will return "< 0.001".
#'
#' @param x a numeric vector or matrix
#' @param method use either [round()] or [signif()] as
#'   rounding method
#' @param digits integer indicating number of decimal places to round to
#' @param sci if `TRUE`, scientific notation is used
#'
#' @return If precision of number is larger than desired rounding, the default
#'   `round` is used. Otherwise, we provide an upper bound instead of
#'   coercion to 0.
#' @author Derek Chiu
#' @export
#'
#' @examples
#' # Vector inputs
#' round_small(2e-04)
#' round_small(5e-04)
#' round_small(6e-04)
#'
#' # Matrix input
#' set.seed(12)
#' x <- matrix(rexp(25, 3), nrow = 5)
#' round_small(x, digits = 1)
round_small <- function(x, method = c("round", "signif"), digits = 3,
                        sci = FALSE) {
  method <- match.arg(method)
  x[] <- simplify2array(lapply(as.list(x), round_s, method, digits, sci))
  return(x)
}

#' Base function for rounding small numbers
#' @noRd
round_s <- function(x, method, digits, sci) {
  if (is.na(x)) {
    return(NA)
  } else {
    if (x <= 5 * 10 ^ -(digits + 1)) {
      return(paste("<", format(10 ^ -digits, scientific = sci)))
    } else {
      f <- switch(method, round = round, signif = signif)
      return(f(x, digits))
    }
  }
}
