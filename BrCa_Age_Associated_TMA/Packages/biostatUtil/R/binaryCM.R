#' Confusion matrix summaries
#'
#' Calculates summaries from cross-tabulated reference and prediction labels for
#' a two-class variable.
#'
#' Given two dichotomous variables summarized in a confusion matrix, this
#' function provides performance summaries. The accuracy, sensitivity,
#' specificity, positive predictive value (PPV), negative predictive value
#' (NPV), and the kappa statistic, along with their bootstrapped confidence
#' intervals are returned.
#'
#' Note that the classes given in `x` and `y` must be binary.
#'
#' @param x a vector of reference classes
#' @param y a vector of predicted classes
#' @param pcond a string or value to be considered positive condition; defaults
#'   to 1.
#' @param seed random seed for bootstrapping
#' @param num.boot number of times to bootstrap. Defaults to 1000.
#' @param conf.level confidence level. Defaults to 95\%.
#' @param digits number of digits to round summaries to
#' @param method method for obtaining confidence intervals for binomial
#'   probabilities. See `Hmisc::binconf` for details.
#' @param verbose logical; if `TRUE`, outputs are printed to the screen
#' @return A confusion matrix for the predicted and reference classes. Then
#' the estimated statistics along with bootstrapped confidence intervals. A
#' list with the following elements
#' \item{CM}{The confusion matrix, whose columns are the predicted conditions
#' and its rows are the true conditions}
#' \item{Accuracy}{Accuracy point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{Sensitivity}{Sensitivity point estimate, lower bound and upper bound
#' for bootstrapped CI}
#' \item{Specificity}{Specificity point estimate, lower bound and upper bound
#' for bootstrapped CI}
#' \item{PPV}{PPV point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{NPV}{NPV point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{kappa}{kappa point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{table}{a data frame that contains all 6 of the estimated statistics
#' along with confidence intervals}
#' @family confusion matrix functions
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' ### 95% CI from 1000 bootstraped samples
#' set.seed(547)
#' n <- 80
#' x <- rbinom(n, size = 1, prob = 0.6)
#' y <- rbinom(n, size = 1, prob = 0.4)
#' binaryCM(x, y)
#'
#' ### 90% CI from 500 bootstrapped samples
#' binaryCM(x, y, num.boot = 500, conf.level = 0.90)
#'
#' ### Round to 2 digits
#' binaryCM(x, y, digits = 2)
binaryCM <- function(x, y, seed = 20, num.boot = 1000, pcond = 1,
                     conf.level = 0.95, digits = 4, method = "wilson",
                     verbose = FALSE) {
  CM0 <- table(Reference = x, Prediction = y)
  if (dim(CM0)[1] > 2 || dim(CM0)[2] > 2) {
    stop("This function only works for binary classes!")
  }
  if (verbose) {
    cat("Confusion Matrix", "\n")
    print(CM0)
  }
  # Check the position of the positive condition and move it to [1,1] in the
  # table
  xppos <- match(pcond, rownames(CM0))
  newxpos <- c(xppos, which(!(c(1, 2) %in% xppos)))
  yppos <- match(pcond, colnames(CM0))
  newypos <- c(yppos, which(!(c(1, 2) %in% yppos)))
  CM <- CM0[newxpos, newypos]
  # Compute marginal totals
  c1 <- CM[1, 1] + CM[2, 1]
  c2 <- CM[1, 2] + CM[2, 2]
  r1 <- CM[1, 1] + CM[1, 2]
  r2 <- CM[2, 1] + CM[2, 2]
  d0 <- CM[1, 1] + CM[2, 2]
  # Calculate statistics
  successes <- list(d0, CM[1, 1], CM[2, 2], CM[1, 1], CM[2, 2])
  observations <- list(c1 + c2, r1, r2, c1, c2)
  kappa <- kappaBootCI(x, y, seed, num.boot, conf.level)
  ci.lo <- paste0((1 - conf.level) / 2 * 100, "%")
  ci.hi <- paste0((1 - (1 - conf.level) / 2) * 100, "%")
  stats <- purrr::map2(successes, observations, Hmisc::binconf,
                       alpha = 1 - conf.level, method = method) %>%
    c(list(matrix(kappa, nrow = 1, dimnames = list("")))) %>%
    purrr::map(round, digits) %>%
    purrr::set_names(c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV",
                       "kappa")) %>%
    purrr::map(~ .x %>% magrittr::set_colnames(c("PointEst", ci.lo, ci.hi)))
  # Optionally print stats to screen
  if (verbose) {
    purrr::map2_chr(names(stats), stats,
                    ~ paste0("\n", .x, ": ", printCI(.y))) %>%
      c("\n\n") %>%
      cat()
  }
  # Result table
  table <- stats %>%
    purrr::invoke(rbind, .) %>%
    magrittr::set_rownames(names(stats)) %>%
    magrittr::set_colnames(c("Point Estimate", "Lower CI", "Upper CI"))
  c(dplyr::lst(CM), stats, dplyr::lst(table))
}
