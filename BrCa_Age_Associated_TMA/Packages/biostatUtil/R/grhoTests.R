#' Summary of Survival Curve differences using G-rho tests
#'
#' Runs the log-rank test and Breslow test to test for difference between two or
#' more survival curves.
#'
#' The log-rank test corresponds to `rho = 0` and the Breslow test
#' corresponds to `rho = 1` in `survdiff`. This function emulates the
#' "Overall Comparisons" table output from SPSS.
#'
#' @param formula formula expression of the form `Surv(time, status) ~
#'   predictors`
#' @param data data frame where variables from `formula` originate
#' @param digits number of significant digits to retain
#' @return The Chi-Square statistic, degrees of freedom, and p-value are given
#'   for both G-rho tests.
#' @author Derek Chiu
#' @export
#'
#' @examples
#' library(survival)
#' grhoTests(Surv(futime, fustat) ~ rx, data = ovarian)
grhoTests <- function(formula, data, digits = 4) {
  tab <- purrr::map(0:1, ~ {
    mod <- survdiff(formula, data, rho = .x)
    c <- mod$chisq
    d <- length(mod$n) - 1
    p <- stats::pchisq(c, d, lower.tail = FALSE)
    data.frame(c, d, p)
  }) %>%
    purrr::invoke(rbind, .) %>%
    signif(digits) %>%
    magrittr::set_rownames(c("Log Rank (Mantel-Cox)",
                             "Breslow (Generalized Wilcoxon)")) %>%
    magrittr::set_colnames(c("Chi-Square", "df", "Sig."))
  return(tab)
}
