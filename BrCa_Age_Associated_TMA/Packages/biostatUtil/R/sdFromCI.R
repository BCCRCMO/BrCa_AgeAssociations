#' Standard deviation of log hazard ratio
#'
#' Compute the standard deviation of a log hazard ratio when only given its
#' repoted confidence limits.
#'
#' Extract standard deviation from a confidence limit and obtain a corresponding
#' p-value, which can be compared to the reported p-value. The null p-values are
#' calculated from the Chi-squared 1 distribution.
#' @param HR reported hazard ratio
#' @param lower.limit reported lower confidence limit
#' @param upper.limit reported upper confidence limit
#' @param alpha significance level
#' @return A list with elements
#' \item{sd}{standard deviation of log HR}
#' \item{pval}{estimated p-value}
#' @author Aline Talhouk, Derek Chiu
#' @references http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1920534/
#' @export
#' @examples
#' library(meta)
#' studlab <- c("PORTEC", "Leuven", "TCGA", "Billingsley")
#' effects <- log(c(0.43, 0.18, 0.12, 0.37))
#' se_effects <- c(sdFromCI(0.43, 0.13, 1.37)$sd,
#'                 sdFromCI(0.18, 0.01, 3.11)$sd,
#'                 sdFromCI(0.12, 0.01, 2.11)$sd,
#'                 sdFromCI(0.37, 0.09, 1.54)$sd)
#' metagen(effects, se_effects, studlab = studlab, sm = "HR", comb.fixed = TRUE)
sdFromCI <- function(HR, lower.limit, upper.limit, alpha = 0.05) {
  sdlnHR <- sqrt(((log(upper.limit) - log(lower.limit)) /
                    (2 * stats::qnorm(1 - alpha / 2))) ^ 2)
  estPval <- stats::pchisq((log(HR) / sdlnHR) ^ 2, 1, lower.tail = FALSE)
  return(list(sd = sdlnHR, Pval = estPval))
}
