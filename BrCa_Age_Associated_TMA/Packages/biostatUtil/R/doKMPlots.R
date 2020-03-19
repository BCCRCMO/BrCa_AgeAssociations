#' Make Kaplan-Meier plots
#'
#' @param input.d `data.frame` containing data
#' @param time follow up time
#' @param status status indicator
#' @param var.name name of variable to make Kaplan-Meier plots on
#' @param var.description description for `var.name`
#' @param surv.type survival outcome. Either "os", "dss", or "pfs".
#' @param shading.colors colors for survival curves
#' @param line.name names for each survival curve
#' @param line.pattern line type for survival curves
#' @param legend logical; if `TRUE`, the legend is overlaid on the graph
#'   (instead of on the side).
#' @param cox.ref.group specify reference group for cox model i.e. hazard
#'   ratio(s)
#' @param use.firth Whether to use Firth's correction for plotting the curves
#' @param CI logical; if `TRUE`, will plot confidence bands
#' @param HR logical; if `TRUE`, will show hazard ratios
#' @param show.risk logical; if `TRUE`, will show the number of people at
#'   risk at each time of death beneath the plot
#' @param km.plot.ref.group specify KM plot reference group; "single" means a
#'   lump log-rank statistic
#' @param single.test.type test to use for survival curves. Defaults to
#'   "logrank".
#' @param use.ggkm if `TRUE`, will use function `ggkm` for plotting
#' @param ... additional arguments to other functions and methods
#' @return A Kaplan-Meier plot for the specified survival outcome split on the
#'   desired variable.
#' @author Samuel Leung, Derek Chiu
#' @export
doKMPlots <- function(input.d, time, status, var.name, var.description,
                      surv.type = c("os", "dss", "pfs"),
                      shading.colors = c("blue2", "red2", "deepskyblue",
                                         "indianred3"),
                      line.name = NULL, line.pattern = NULL, legend = FALSE,
                      cox.ref.group = NULL, use.firth = -1, CI = TRUE,
                      HR = TRUE, show.risk = TRUE, km.plot.ref.group = "single",
                      single.test.type = "logrank", use.ggkm = FALSE, ...) {

  input.d[var.name] <- droplevels(input.d[var.name])
  input.d[[time]] <- as.numeric(input.d[[time]])
  levs <- names(table(input.d[, var.name]))
  if (is.null(line.name))
    line.name <- levs
  if (is.null(shading.colors))
    shading.colors <- seq_along(levs)
  formula.obj <- stats::as.formula(paste0("Surv(", time, ", ", status, ") ~ ", var.name))
  main <- ifelse(is.na(var.description) | var.description == "", "",
                 paste0(var.description, " (", toupper(match.arg(surv.type)), ")"))

  if (use.ggkm) {
    pos <- 1
    assign("formula.obj", formula.obj, envir = as.environment(pos))
    assign("input.d", input.d, envir = as.environment(pos))
    sfit <- survival::survfit(formula.obj, data = input.d)

    ggkm(sfit, shading.colors = shading.colors, ystratalabs = line.name,
         line.pattern = line.pattern, legend = legend,
         cox.ref.grp = cox.ref.group, use.firth = use.firth, CI = CI, HR = HR,
         table = show.risk, main = main, ...)

  } else {
    plotKM(input.d = input.d, input.formula = formula.obj,
           line.color = shading.colors, line.name = line.name,
           line.pattern = line.pattern, show.test = km.plot.ref.group,
           single.test.type = single.test.type, main.text = main, ...)
  }
}
