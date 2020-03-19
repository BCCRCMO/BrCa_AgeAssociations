#' Do Schoenfeld residual plot
#'
#' @param input.d input `data.frame`
#' @param input.formula survival formula to `Surv`
#' @param vars.to.plot variable(s) to plot
#' @param main title for each variable (i.e. main can be an array)
#' @author Samuel Leung
#' @export
#' @examples
#' # Base output
#' test1 <- list(time = c(4, 3, 1, 1, 2, 2, 3),
#' status = c(1, 1, 1, 0, 1, 1, 0),
#' x = c(0, 2, 1, 1, 1, 0, 0),
#' sex = c(0, 0, 0, 0, 1, 1, 1))
#'
#' # Pretty output
#' plotSchoenfeld(test1, Surv(time, status) ~ x + strata(sex), "x")
plotSchoenfeld <- function(input.d, input.formula, vars.to.plot = NULL,
                           main = NULL) {
  pos <- 1
  assign(".my.data", input.d, envir = as.environment(pos))
  assign(".my.formula", stats::as.formula(paste0("survival::", deparse(input.formula))),
         envir = as.environment(pos))
  fit <- coxph(get(".my.formula"), get(".my.data"))
  var.names <-
    all.vars(input.formula[[3]])  # there may be > 1 terms in formula
  for (var.name in var.names) {
    names(fit$coefficients) <- names(fit$coefficients) %>%
      purrr::map_chr(~ ifelse(var.name == .x, .x, sub(var.name, "", .x)))
  }
  fit.zph <- survival::cox.zph(fit)
  # determine which variable to plot
  if (is.null(vars.to.plot)) {
    graphics::plot(fit.zph, main = ifelse(is.null(main), "", main)) # plot all
  } else {
    for (i in seq_along(vars.to.plot)) {
      graphics::plot(
        fit.zph,
        var = which(names(fit$coefficients) == vars.to.plot[i]),
        main = ifelse(is.null(main), "", main[min(length(main), i)])
      )
    }
  }
  dplyr::lst(fit, fit.zph)
}
