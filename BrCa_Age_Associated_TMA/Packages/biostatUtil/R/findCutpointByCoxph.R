#' Find cutpoint by Cox model statistics
#'
#' Try to find a single cutpoint by maximizing the hazard ratio or minimizing
#' the AIC.
#'
#' The formula must be univariable.
#'
#' @param input.d input data, typically a matrix or data frame
#' @param surv.formula a formula of type `Surv(time, status) ~ x`, where
#'   `x` is the variable of interest.
#' @return best cutpoint as determined by different metrics. Returns `NA`
#'   if `surv.formula` is not univariable, the variable is not numeric, or
#'   there is no variation in the variable of interest.
#' @section Warning: `surv.formula` cannot be multivariable. For example,
#'   `Surv(time, status) ~ x + age` won't work but `Surv(time, status)
#'   ~ x` is fine.
#' @author Samuel Leung, Derek Chiu
#' @export
findCutpointByCoxph <- function(input.d, surv.formula) {
  if (length(surv.formula[[3]]) > 1) {
    stop("formula must be univariable: Surv(time, status) ~ variable.")
  }

  var.name <- deparse(surv.formula[[3]])
  if (!is.numeric(input.d[, var.name])) {
    stop(paste(var.name, "is not numeric."))
  }

  # find the min,max value of the variable
  possible.values <- sort(unique(input.d[, var.name]))
  if (length(possible.values) < 2) {
    stop(paste0("variable '", var.name, "' does not have any variation (i.e. all values equal) for cutpoint assessment."))
  }

  # try all possible cut point
  cut.points <- aics <- hrs <- c()
  for (cut.point in possible.values[-1]) {
    cox.model.fit <- NA
    tryCatch({
      cox.model.fit <- coxph(stats::as.formula(paste(deparse(surv.formula), ">=", cut.point)), data = input.d)
    }, error = function(e) {
      print(paste("ERROR (find.cut.point.by.coxph) unexpected error occured:", e))
    })
    # cox.model.fit must not be NA, otherwise, length(cox.model.fit) would be 1
    if (length(cox.model.fit) > 1) {
      # cox model with no warning/error
      cut.points <- c(cut.points, cut.point)
      aics <- c(aics, stats::extractAIC(cox.model.fit)[2])
      hrs <- c(hrs, exp(cox.model.fit$coefficients[[1]]))
    }
  }
  return(list("best.cut.point.by.aic" = cut.points[which.min(aics)],
              "best.cut.point.by.hr" = cut.points[which.max(hrs)],
              "cut.points" = cut.points,
              "aics" = aics,
              "hazard.ratios" = hrs))
}
