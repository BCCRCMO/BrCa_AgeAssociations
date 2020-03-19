#' Nicer outputs from coxph model fits
#'
#' Formats the coxph (or coxphf) model outputs into a prettier display.
#'
#' The original fit call is returned, along with other useful statistics. If
#' Firth is called, we first determine if the percentage of censored cases is
#' large enough to use Firth.
#'
#' @param input.formula a formula object, with the response on the left of a
#'   `~` operator, and the terms on the right. The response must be a
#'   survival object as returned by the `Surv` function.
#' @param input.d data.frame containing the variables in `input.formula`
#' @param ref.grp a named array indicating reference group(s) for any variables
#'   in the cox model.  this parameter will be ignored if not applicable, e.g.
#'   for continuous variable
#' @param use.firth if set to -1, Firth's correction is applied using
#' `coxphf`. The default value is 1, which uses `coxph`.
#' @param check.ph if `TRUE`, checks the proportional hazard assumption.
#' @param plot.ph if `TRUE` (default), the PH residual plot is graphed.
#' @param ph.test.plot.filename character string for filename of residual plot,
#'   including file extension. If `NULL` (default), plots to console.
#'   Otherwise plot is saved to disk.
#' @param ... additional arguments to `coxph` or `coxphf`
#' @return A list with elements
#' \item{output}{A character matrix showing the hazard ratio, confidence bounds,
#' and p-value for each coefficient. If `check.ph = TRUE`, the last column
#' shows the results of the PH test, otherwise it reads "NOT CALCULATED" for
#' every coefficient.}
#' \item{fit}{The fit output call returned from `coxph` or `coxphf`}
#' \item{n}{The number of observations used in the fit}
#' \item{nevent}{The number of events used in the fit}
#' \item{ph.test}{The rho, Chi-squared statistic, and p-value for each coefficient.
#' Not shown if `check.ph = FALSE`}
#' \item{used.firth}{logical; if `TRUE`, Firth's correction was applied
#' using `coxphf`}
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' library(survival)
#' # Base output
#' test1 <- list(time = c(4, 3, 1, 1, 2, 2, 3),
#' status = c(1, 1, 1, 0, 1, 1, 0),
#' x = c(0, 2, 1, 1, 1, 0, 0),
#' sex = c(0, 0, 0, 0, 1, 1, 1))
#' coxph(Surv(time, status) ~ x + strata(sex), test1)
#'
#' # Pretty output
#' prettyCoxph(Surv(time, status) ~ x + strata(sex), test1)
#'
#' # x is now a factor
#' test1$x <- factor(test1$x)
#' prettyCoxph(Surv(time, status) ~ x + strata(sex), test1)
#'
#' # Releveled reference group
#' prettyCoxph(Surv(time, status) ~ x + strata(sex), test1, ref.grp =
#' setNames("2", "x"))
#'
#' # Use Firth's correction
#' prettyCoxph(Surv(time, status) ~ x + strata(sex), test1, use.firth = -1)
prettyCoxph <- function(input.formula, input.d, ref.grp = NULL, use.firth = 1,
                        check.ph = FALSE, plot.ph = TRUE,
                        ph.test.plot.filename = NULL, ...) {
  pos <- 1
  assign(".my.formula", stats::as.formula(
    paste0(paste(deparse(input.formula), collapse = ""))),
    envir = as.environment(pos))
  # Modify input.d if ref.grp is defined
  if (!is.null(ref.grp)) {
    for (var.name in names(ref.grp)) {
      if (var.name %in% all.vars(input.formula[[3]])) {
        # Only relevel if var.name in formula
        input.d[[var.name]] <- stats::relevel(factor(input.d[[var.name]]),
                                              ref = ref.grp[var.name])
      }
    }
  }

  assign(".my.data", input.d, envir = as.environment(pos))
  ok.to.use.firth <- ifelse(use.firth == -1, TRUE, FALSE)
  if (use.firth < 1 & use.firth > -1) {
    for (var.name in all.vars(.my.formula)[-c(1, 2)]) {
      if (is.factor(.my.data[[var.name]])) {
        fit <- survfit(stats::as.formula(paste(deparse(.my.formula[[2]]),
                                               "~", var.name)), data = .my.data)
        for (i in seq_len(nrow(fit))) {
          if (sum(fit[i]$n.censor) / fit[i]$n > use.firth) {
            ok.to.use.firth <- TRUE
            break
          }
        }
        if (ok.to.use.firth)
          break
      }
    }
  }
  fit <- coxph(.my.formula, .my.data, ...)
  .my.formula <- fit$formula

  if (ok.to.use.firth) {
    .my.data <- .my.data %>%
      as.data.frame() %>%
      magrittr::extract(apply(.[all.vars(.my.formula)], 1,
                              function(x) !any(is.na(x))), )
    fit.firth <- coxphf::coxphf(.my.formula, .my.data, maxit = 100)
    fit.firth$nevent <- sum(fit.firth$y[, "status"])
    fit.firth$ci.lower[is.nan(fit.firth$ci.lower)] <- NA
    fit.firth$ci.upper[is.nan(fit.firth$ci.upper)] <- NA
    # If becomes NaN, change to NA for better display/understanding in KM plot
  } else {
    fit.firth <- NA
  }

  if (check.ph) {
    ph.test <- cox.zph(fit)
    ph.check <- ph.test$table %>%
      magrittr::extract(rownames(.) != "GLOBAL", "p") %>%
      cbind("PH test" = .)
    if (plot.ph) {
      if (is.null(ph.test.plot.filename)) {
        graphics::plot(ph.test)
      } else {
        grDevices::pdf(ph.test.plot.filename)
        graphics::plot(ph.test)
        grDevices::dev.off()
      }
    }
  } else {
    ph.test <- NA
    ph.check <- "NOT CALCULATED"
  }

  if (ok.to.use.firth) {
    result <- fit.firth %>%
      magrittr::extract(c("coefficients", "ci.lower", "ci.upper", "prob")) %>%
      purrr::set_names(c("estimate", "conf.low", "conf.high", "p.value")) %>%
      purrr::map_at("estimate", exp) %>%
      as.data.frame()
  } else {
    result <- fit %>%
      broom::tidy(exponentiate = TRUE) %>%
      as.data.frame() %>%
      magrittr::set_rownames(.[["term"]]) %>%
      magrittr::extract(c("estimate", "conf.low", "conf.high", "p.value"))
  }
  result <- cbind(result, ph.check)
  list(output = result, fit = fit, fit.firth = fit.firth, n = fit$n,
       nevent = fit$nevent, ph.test = ph.test, used.firth = ok.to.use.firth)
}
