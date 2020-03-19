#' Mean and bootstrap confidence interval
#'
#' Finds mean of a vector, with bootstrapped confidence bounds.
#'
#' Takes a numeric vector and resamples with replacement `num.boot` times.
#' If the input vector has any `NA` entries, include the argument
#' `na.rm = TRUE` from `mean`.
#'
#' @param x a vector (or matrix)
#' @param num.boot number of bootstrap samples. Defaults to 1000.
#' @param conf.level confidence level.
#' @param seed random seed for resampling
#' @param ... additional arguments to `mean`
#' @return A list with elements
#' \item{obs.mean}{mean}
#' \item{ci}{bootstraped confidence interval}
#' \item{n}{number of bootstrap samples used}
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' ## Vectors
#' set.seed(344)
#' bootMean(rnorm(100, 5, 3))
#' bootMean(rnorm(100, 5, 3), num.boot = 500)
#' bootMean(rnorm(100, 4, 3), conf.level = 0.90)
#'
#' ## Missing Values
#' set.seed(344)
#' s <- ifelse(rexp(100, 0.5) < 1, NA, rexp(100, 0.5))
#' bootMean(s)  # doesn't work
#' bootMean(s, na.rm = TRUE)
bootMean <- function(x, num.boot = 1000, conf.level = 0.95, seed = 12, ...) {
  set.seed(seed)
  obs.mean <- mean(x, ...)
  ci <- replicate(num.boot, mean(sample(x, replace = TRUE), na.rm = TRUE)) %>%
    sort() %>%
    magrittr::extract(c(floor(num.boot * (1 - conf.level) / 2),
                        ceiling(num.boot * (1 - (1 - conf.level) / 2))))
  list(obs.mean = obs.mean, ci = ci, n = length(x))
}
