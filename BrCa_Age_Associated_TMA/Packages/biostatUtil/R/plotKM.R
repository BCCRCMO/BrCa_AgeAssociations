#' Plot Kaplan-Meier curves
#'
#' KM plots with basic annotations.
#'
#' @param input.d input `data.frame`
#' @param input.formula survival formula to `Surv`
#' @param main.text plot title
#' @param xlab.text x-axis label
#' @param ylab.text y-axis label
#' @param line.name name of legend
#' @param line.color line colour of survival curves
#' @param line.pattern line pattern of survival curves
#' @param line.width line width of survival curves
#' @param show.test show single or the reference group value (for pairwise
#'   comparisons). If `"none"`, then no test is show.
#' @param single.test.type test to show if specified `show.test = "single"`.
#'   Possible choices are `"logrank"` (default), `"wilcoxon"`, `"taroneware"`,
#'   or `"all"`.
#' @param digits number of digits for p-value
#' @param obs.survyrs show the observed survival years survival rate on KM plot
#' @param legend.pos legend position keyword
#' @param file.name name of file to save plot to
#' @param file.width width of figure in saved file
#' @param file.height height of figure in saved file
#' @param grey.scale logical. If `TRUE`, the plot will be in grey scale.
#' @param show.single.test.pos position to show single test; defaults to 0.5 if
#'   `legend.pos = "top"`. Otherwise 0.1
#' @param xlabs dummy param to match ggkm so they don't get passed to ...
#' @param legend.xy dummy param to match ggkm so they don't get passed to ...
#' @param ylabs dummy param to match ggkm so they don't get passed to ...
#' @param timeby dummy param to match ggkm so they don't get passed to ...
#' @param ... additional arguments to `plotKMDetail`
#' @author Samuel Leung
#' @seealso [plotKMDetail()]
#' @export
plotKM <- function(input.d, input.formula,
                   main.text = "", xlab.text = "", ylab.text = "",
                   line.name, line.color, line.pattern = NULL,
                   line.width = NULL, show.test = "single",
                   single.test.type = "logrank", digits = 3,
                   obs.survyrs = 3, legend.pos = "bottomleft",
                   file.name = "no.file", file.width = 7, file.height = 7,
                   grey.scale = FALSE, show.single.test.pos = 0.1,
                   xlabs = NULL, legend.xy = NULL, ylabs = NULL, timeby = NULL,
                   ...) {

  # Calculate "obs.survyrs"-yrs survival
  fit.obj <- survfit(input.formula, data = input.d)
  summary.surv.fit <- summary(fit.obj, time = obs.survyrs, extend = TRUE)
  surv.time <- lapply(summary.surv.fit[c("surv", "lower", "upper")],
                      function(x) format(x * 100, digits = 3))
  ten.years.surv.95CI <- mapply(function(s, l, u)
    paste0(s, "% (", l, "% - ", u, "%)"),
    s = surv.time$surv, l = surv.time$lower, u = surv.time$upper)

  # Summary of survival object to end of followup to calculate event count
  summary.surv.fit.all <- summary.surv.fit[["table"]]
  event.count <- apply(summary.surv.fit.all[, c("events", "records")], 1,
                       paste, collapse = "/")

  # If legend is on top, the test position is at halfway point of plot
  if (legend.pos == "top") {
    show.single.test.pos <- 0.5
  }

  # Plot km
  output <- plotKMDetail(
    input.data = input.d, surv.formula = input.formula,
    main.text = main.text, xlab.text = xlab.text, ylab.text = ylab.text,
    line.name = line.name, line.color = line.color, line.pattern = line.pattern,
    line.width = line.width, show.test = show.test,
    single.test.type = single.test.type, round.digits.p.value = digits,
    obs.survyrs = obs.survyrs, ten.years.surv.95CI = ten.years.surv.95CI,
    event.count = event.count, legend.pos = legend.pos,
    file.name = file.name, file.width = file.width, file.height = file.height,
    grey.scale = grey.scale, show.single.test.pos = show.single.test.pos,
    mark.time = TRUE, # for plot.survfit to show censor time point
    ...)
  return(list("log.rank.p.values" = output$log.rank.p.values,
              "wilcox.p.values" = output$wilcox.p.values,
              "n" = sum(fit.obj$n),
              "nevent" = sum(fit.obj$n.event)))
}
