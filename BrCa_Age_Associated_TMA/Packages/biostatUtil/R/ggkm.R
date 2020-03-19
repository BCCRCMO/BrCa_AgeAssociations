#' Kaplan-Meier Plots using ggplot
#'
#' Produce nicely annotated KM plots using ggplot style.
#'
#' @param sfit an object of class `survfit` containing one or more survival
#'   curves
#' @param sfit2 an (optional) second object of class `survfit` to compare
#'   with `sfit`
#' @param table logical; if `TRUE` (default), the numbers at risk at each
#'   time of death is shown as a table underneath the plot
#' @param returns logical; if `TRUE` the plot is returned
#' @param marks logical; if `TRUE` (default), censoring marks are shown on
#'   survival curves
#' @param CI logical; if `TRUE` (default), confidence bands are drawn for
#'   survival curves the using cumulative hazard, or log(survival).
#' @param line.pattern linetype for survival curves
#' @param shading.colors vector of colours for each survival curve
#' @param main plot title
#' @param xlabs horizontal axis label
#' @param ylabs vertical axis label
#' @param xlims horizontal limits for plot
#' @param ylims vertical limits for plot
#' @param ystratalabs labels for the strata being compared in `survfit`
#' @param cox.ref.grp indicates reference group for the variable of interest in
#'   the cox model.  this parameter will be ignored if not applicable, e.g. for
#'   continuous variable
#' @param timeby length of time between consecutive time points spanning the
#'   entire range of follow-up. Defaults to 5.
#' @param pval logical; if `TRUE` (default), the logrank test p-value is
#'   shown on the plot
#' @param HR logical; if `TRUE` (default), the estimated hazard ratio and
#'   its 95\% confidence interval will be shown
#' @param use.firth Firth's method for Cox regression is used if the percentage
#'   of censored cases exceeds `use.firth`. Setting `use.firth = 1`
#'   (default) means Firth is never used, and `use.firth = -1` means Firth
#'   is always used.
#' @param legend logical; if `TRUE`, the legend is overlaid on the graph
#'   (instead of on the side).
#' @param legend.xy named vector specifying the x/y position of the legend
#' @param legend.direction layout of items in legends ("horizontal" (default) or
#'   "vertical")
#' @param line.y.increment how much y should be incremented for each line
#' @param digits number of digits to round: p-values digits=nunber of
#'   significant digits, HR digits=number of digits after decimal point NOT
#'   significant digits
#' @param ... additional arguments to other methods
#' @return A kaplan-meier plot with optional annotations for hazard ratios, log
#'   rank test p-values, and risk table counts for each stratum.
#'
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' library(survival)
#' sfit <- survfit(Surv(time, status) ~ sex, lung)
#' ggkm(sfit, timeby = 200, main = "Survival curves by sex")
ggkm <- function(sfit, sfit2 = NULL, table = TRUE, returns = TRUE, marks = TRUE,
                 CI = TRUE, line.pattern = NULL, shading.colors = NULL,
                 main = "Kaplan-Meier Plot", xlabs = "Time",
                 ylabs = "Survival Probability", xlims = NULL, ylims = NULL,
                 ystratalabs = NULL, cox.ref.grp = NULL, timeby = 5,
                 pval = TRUE, HR = TRUE, use.firth = 1, legend = FALSE,
                 legend.xy = NULL, legend.direction = "horizontal",
                 line.y.increment = 0.05, digits = 3, ...) {
  time <- surv <- lower <- upper <- n.censor <- n.risk <- n.event <-
    estimate <- conf.high <- conf.low <- NULL
  times <- seq.int(0, max(sfit$time), by = timeby)
  s1 <- levels(summary(sfit)$strata)
  s2 <- summary(sfit, censored = TRUE)$strata
  s3 <- summary(sfit, times = times, extend = TRUE)$strata

  # Specifying plot parameter defaults
  shading.colors <- shading.colors %||% c("blue2", "red2",
                                          "deepskyblue", "indianred3")
  xlims <- xlims %||% c(0, max(sfit$time))
  ylims <- ylims %||% c(0, 1)
  legend.xy <- legend.xy %||% c(0.8, 0.88)
  ystratalabs <- ystratalabs %||% gsub(".*=(.)", "\\1", names(sfit$strata))
  if (is.null(line.pattern) | length(line.pattern) == 1)
    line.pattern <- stats::setNames(rep(1, length(sfit$strata)), ystratalabs)
  if (!is.null(cox.ref.grp))
    names(cox.ref.grp) <- all.vars(sfit$call)[3]  # works for two-sided formulas

  # Left margins for km plot and risk table
  mleft <- left_margin(ystratalabs)

  # Data for KM plot
  .df <- sfit %>%
    broom::tidy() %>%
    select(time, n.risk, n.event, n.censor, surv = estimate,
           strata, upper = conf.high, lower = conf.low) %>%
    mutate(strata = factor(s2, labels = ystratalabs)) %>%
    bind_rows(
      data.frame(time = 0, surv = 1,
                 strata = factor(ystratalabs, levels = levels(.$strata)),
                 upper = 1, lower = 1), .)

  # KM plot
  p <- ggplot(.df, aes(time, surv, color = strata, fill = strata,
                       linetype = strata)) +
    geom_step(size = .7) +
    scale_colour_manual(values = shading.colors) +
    scale_fill_manual(values = shading.colors) +
    scale_linetype_manual(values = line.pattern) +
    scale_x_continuous(xlabs, breaks = times, limits = xlims) +
    scale_y_continuous(ylabs, limits = ylims) +
    theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.5),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.margin = grid::unit(c(0.5, 1, .5, mleft$margin.km),
                                   c("lines", "lines", "lines", "in")),
          legend.position = "none") +
    ggtitle(main)
  if (legend)  # Legend
    p <- p + theme(legend.position = legend.xy,
                   legend.key = element_rect(colour = NA),
                   legend.title = element_blank(),
                   legend.direction = legend.direction)
  if (CI)  # Confidence Bands
    p <- p + geom_ribbon(data = .df, aes(ymin = lower, ymax = upper),
                         alpha = 0.05, linetype = 0)
  if (marks)  # Censor Marks
    p <- p + geom_point(data = subset(.df, n.censor >= 1),
                        aes(x = time, y = surv), shape = "/", size = 4)

  # HR statistic (95% CI), log rank test p-value for sfit (or sfit2, if exists)
  if (pval) {
    fit <- sfit2 %||% sfit
    p <- summarize_km(fit = fit, p = p, digits = digits, HR = HR,
                      cox.ref.grp = cox.ref.grp, use.firth = use.firth,
                      ystratalabs = ystratalabs,
                      line.y.increment = line.y.increment)
  }

  # Create table graphic to include at-risk numbers, keep at-risk numbers
  # same order as appears in HR (do not reverse levels)
  if (table) {
    risk.data <- sfit %>%
      summary(times = times, extend = TRUE) %>%
      list() %>%
      purrr::map_df(`[`, c("strata", "time", "n.risk"))
    data.table <- ggplot(risk.data, aes(x = time, y = strata,
                                        label = format(n.risk, nsmall = 0))) +
      geom_text(size = 3.5) +
      scale_y_discrete(labels = ystratalabs) +
      scale_x_continuous("Numbers at risk", limits = xlims) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = shading.colors,
                                       face = "bold", hjust = 1),
            axis.title.x = element_text(size = 12, vjust = 1),
            legend.position = "none",
            plot.margin = grid::unit(c(0.5, 1.25, 0.5, mleft$margin.rt),
                                     "lines")) +
      labs(y = NULL)
    if (returns)
      gridExtra::grid.arrange(p, data.table, clip = FALSE, nrow = 2, ncol = 1,
                              heights = grid::unit(c(2, .5), "null"))
  } else {
    return(p)
  }
}

#' Get left margin distances for km plot and risk table
#' @noRd
left_margin <- function(labels) {
  max.nc <- max(nchar(labels))
  if (max.nc <= 4) {
    mleft.km <- 0
    mleft.rt <- 2.4
  } else {
    mleft.km <- graphics::strwidth(labels[which.max(nchar(labels))],
                                   units = "in")
    if (max.nc <= 5) {
      mleft.km <- mleft.km - (max.nc / 100 + 0.4)
    } else {
      mleft.km <- mleft.km - (max.nc / 100 + 0.5)
    }
    mleft.rt <- 0.5
  }
  return(list(margin.km = mleft.km,
              margin.rt = mleft.rt))
}

#' Numerical summaries of km fit: HR (95\% CI), Log rank test p-value
#' @noRd
summarize_km <- function(fit, p, digits, HR, cox.ref.grp,
                         use.firth, ystratalabs, line.y.increment) {
  f <- eval(fit$call$formula)
  d <- eval(fit$call$data)
  pvalue <- survdiff(f, d) %>%
    getPval() %>%
    round_small(method = "signif", digits = digits)
  if (is.numeric(pvalue)) {
    pvalsep <- " = "
    # the following line tries to figure out how to print p-values with specified digits
    # e.g. digits=2 -> "0.20" NOT "0.2"
    pvalue <- sprintf(paste0("%.", (max(digits, min(grep("[1-9]", strsplit(as.character(pvalue), "")[[1]]) - 3 + digits))), "f"), pvalue)
  } else {
    pvalsep <- " "
  }
  pvaltxt <- paste("Log Rank p", pvalue, sep = pvalsep)
  if (HR) {
    pretty.coxph.obj <- prettyCoxph(input.formula = f,
                                    input.d = d,
                                    ref.grp = cox.ref.grp,
                                    use.firth = use.firth)
    if (pretty.coxph.obj$used.firth) {
      coxm <- pretty.coxph.obj$fit.firth
    } else {
      coxm <- pretty.coxph.obj$fit
    }
    HRtxts <- Xunivcoxph(coxm, digits = digits)
    cox.strata.labs <- ystratalabs
    if (!is.null(cox.ref.grp)) {
      cox.strata.labs <- c(cox.ref.grp,
                           ystratalabs[ystratalabs != cox.ref.grp])
    }
    for (i in seq_along(HRtxts)) {
      HRtxt <- HRtxts[i]
      HRtxt <- paste0(HRtxt, " ~ ", cox.strata.labs[i + 1],
                        " vs. ", cox.strata.labs[1])
      p <- p + annotate("text", x = 0.2 * max(fit$time), hjust = 0,
                        y = 0.01 + line.y.increment * i, label = HRtxt,
                        size = 3)
    }
  }
  p <- p + annotate("text", x = 0.2 * max(fit$time), hjust = 0, y = 0.01,
                    label = pvaltxt, size = 3)
}
