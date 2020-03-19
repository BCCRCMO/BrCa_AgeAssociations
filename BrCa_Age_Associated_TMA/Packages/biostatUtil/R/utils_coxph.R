#' Nice output from Cox regression object
#'
#' Nicely formatted output from a Cox regression object, either `coxph` or
#' `coxphf`.
#'
#' For objects of class `coxphf`, the calculation of the number of events
#' used in the fit is slightly different.
#'
#' @param object a model fit object returned from `coxph` or `coxphf`
#' @param coefnames a vector of labels for the coefficient names returned by the
#'   fit. `NULL` by default, uses original coefficient names.
#' @param conf.level confidence level for hazard ratio. Defaults to 95\%.
#' @param digits number of digits to round to
#' @return A matrix with a row for each coefficient, and a column for each of
#'   the following elements
#' \item{n}{the number of observations used in the fit.}
#' \item{events}{the number of events used in the fit.}
#' \item{coef}{log hazard ratio}
#' \item{se}{standard error of `coef`}
#' \item{Z-Score}{Wald test statistic}
#' \item{P-value}{p-value for `Z-Score`}
#' \item{HR}{hazard ratio}
#' \item{Lower CI Limit}{Defaults to `2.5 \%`}
#' \item{Upper CI Limit}{Defaults to `97.5 \%`}
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' library(survival)
#' test1 <- list(time = c(4, 3, 1, 1, 2, 2, 3),
#' status = c(1, 1, 1, 0, 1, 1, 0),
#' x = c(0, 2, 1, 1, 1, 0, 0),
#' sex = c(0, 0, 0, 0, 1, 1, 1))
#'
#' # Stratified
#' mod1 <- coxph(Surv(time, status) ~ x + strata(sex), test1)
#' coxphOut(mod1)
#'
#' # Not stratified
#' mod2 <- coxph(Surv(time, status) ~ x + sex, test1)
#' coxphOut(mod2, coefnames = c("x", "gender"))
coxphOut <- function(object, coefnames = NULL, conf.level = 0.95,
                     digits = 2) {
  cox <- object
  n <- cox$n
  events <- ifelse("coxphf" %in% class(cox), sum(cox$y[, "status"]), cox$nevent)
  coef <- cox$coef
  se <- sqrt(diag(cox$var))
  z <- coef / se
  p <- 1 - stats::pchisq(z ^ 2, 1)
  HR <- exp(coef)
  CI <- exp(stats::confint(cox, level = conf.level))
  tmp <- cbind(n, events, coef, se, "Z-Score" = z, "P-value" = p, HR, CI)
  if (!is.null(coefnames))
    rownames(tmp) <- coefnames
  return(round(tmp, digits))
}

#' Univariate cox proprtional hazards model
#'
#' Concatenates hazard ratios and confidence limits for every covariate in a Cox
#' model.
#'
#' @param mod model fit object, returned from either `coxph` or
#'   `coxphf`.
#' @param digits number of digits to round
#' @return Hazard ratio and 95/% confidence interval
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' library(survival)
#' library(coxphf)
#' # One predictor
#' test1 <- list(
#'   time = c(4, 3, 1, 1, 2, 2, 3),
#'   status = c(1, 1, 1, 0, 1, 1, 0),
#'   x = c(0, 2, 1, 1, 1, 0, 0),
#'   sex = c(0, 0, 0, 0, 1, 1, 1)
#' )
#' mod <- coxph(Surv(time, status) ~ x + strata(sex), test1)
#' Xunivcoxph(mod)
#'
#' # Multiple predictors
#' bladder1 <- bladder[bladder$enum < 5, ]
#' mod <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) +
#' cluster(id), bladder1)
#' Xunivcoxph(mod, digits = 2)
#'
#' # Firth's correction
#' test2 <- data.frame(list(
#'   start = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
#'   stop = c(2, 3, 6, 7, 8, 9, 9, 9, 14, 17),
#'   event = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
#'   x    = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0)
#' ))
#' mod <- coxphf(formula = Surv(start, stop, event) ~ x, pl = FALSE,
#' data = test2)
#' Xunivcoxph(mod)
Xunivcoxph <- function(mod, digits = 3) {
  UseMethod("Xunivcoxph")
}

#' @export
Xunivcoxph.coxph <- function(mod, digits = 3) {
  mod %>%
    broom::tidy(exponentiate = TRUE) %>%
    magrittr::extract(c("estimate", "conf.low", "conf.high")) %>%
    format_hr_ci(digits) %>%
    paste("HR", .)
}

#' @export
Xunivcoxph.coxphf <- function(mod, digits = 3) {
  mod %>%
    magrittr::extract(c("coefficients", "ci.lower", "ci.upper")) %>%
    purrr::map_at("coefficients", exp) %>%
    format_hr_ci(digits) %>%
    paste("HR(F)", .)
}

#' @noRd
format_hr_ci <- function(stats, digits, labels = TRUE, method = c("Inf", "Sci")
                         ) {
  stats %>%
    purrr::map(round, digits) %>%
    purrr::map(sprintf, fmt = paste0("%.", digits, "f")) %>%
    unname() %>%
    purrr::pmap_chr(paste_hr_ci, labels = labels, method = method)
}

#' @noRd
paste_hr_ci <- function(hr, ci.lo, ci.hi, labels = TRUE,
                        method = c("Inf", "Sci")) {
  method <- match.arg(method)
  hr_ci <- list(hr, ci.lo, ci.hi) %>%
    as.numeric() %>%
    purrr::map_if(
      .p = ~ is.na(.) || . > 1000,
      .f = ~ switch(method,
                    `Inf` = "Inf",
                    Sci = format(., digits = 3, scientific = TRUE))
    )
  paste0(hr_ci[[1]], ifelse(labels, " (95% CI: ", " ("),
         hr_ci[[2]], "-", hr_ci[[3]], ")")
}

#' Print Cox model output
#'
#' Prints a Cox model output in a nice HTML table
#' @param cox is the matrix `summary(x)$coef`, where x is an object of
#'   class `coxph`
#' @param Capt a caption
#' @author Aline Talhouk
#' @export
printCoxMod <- function(cox, Capt) {
  TAB <- htmlTable::htmlTable(cox,
                              rowlabel = "Predictors",
                              caption = Capt,
                              ctable = TRUE)
  pander::pander(TAB, style = "rmarkdown")
}

#' Create a survival formula from time, status, event code and terms strings
#' @noRd
surv_formula <- function(time, status, event, terms) {
  stats::as.formula(paste0("Surv(", time, ", ", status, " == '", event, "') ~ ",
                           paste(terms, collapse = " + ")))
}

#' Round p-values according to specifications for coxph summaries
#' @noRd
round_pval <- function(pvalue, round.small = TRUE, scientific = FALSE,
                       digits = 4) {
  if (round.small) {
    p <- round_small(pvalue, method = "round", digits = digits,
                     sci = scientific)
    if (grepl("<", p)) {
      p
    } else {
      sprintf(paste0("%.", digits, "f"), p)
    }
  } else {
    sprintf(paste0("%.", digits, "f"), round(pvalue, digits = digits))
  }
}
