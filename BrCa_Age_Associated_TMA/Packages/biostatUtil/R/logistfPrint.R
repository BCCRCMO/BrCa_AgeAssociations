#' Print summaries from logistf
#'
#' Print nice summaries from objects returned by `logistf`.
#'
#' There is a lot of raw output from `logistf`, and it is not
#' easy to extract the coefficient table. This function provides
#' a convenient wrapper to return summaries from `logistf`.
#'
#' @param fit a fit object returned by `logistf`
#' @param digits number of digits to round to
#' @return A data frame with a row for each predictor and the
#' following columns
#' \item{coef}{log odds ratio}
#' \item{exp(coef)}{odds ratio}
#' \item{lower}{lower confidence bound of odds ratio for specified level}
#' \item{upper}{upper confidence bound of odds ratio for specified level}
#' \item{p}{p-value}
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' library(logistf)
#' data(sex2)
#' fit <- logistf(case ~ age + oc + vic + vicl + vis + dia, data = sex2,
#' alpha = 0.1)
#' summary(fit)
#'
#' ## Streamlined summary
#' logistfPrint(fit)
logistfPrint <- function(fit, digits = 3) {
  conf.level <- 1 - fit$alpha
  colnames <- c("coef", "exp(coef)",
                paste("lower", conf.level), paste("upper", conf.level), "p")
  output <- cbind(fit$coefficients, exp(fit$coefficients),
                  exp(fit$ci.lower), exp(fit$ci.upper), fit$prob)
  dimnames(output)[[2]] <- colnames
  return(round(output, digits))
}
