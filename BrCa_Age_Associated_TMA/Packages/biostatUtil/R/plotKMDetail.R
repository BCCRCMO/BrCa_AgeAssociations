#' Plot detailed Kaplan-Meier curves
#'
#' KM plots with details of event counts.
#'
#' @param input.data input `data.frame`
#' @param surv.formula survival formula to `Surv`
#' @param main.text plot title
#' @param xlab.text x-axis label
#' @param ylab.text y-axis label
#' @param line.name name of legend
#' @param line.color line colour of survival curves
#' @param line.pattern line pattern of survival curves
#' @param line.width line width of survival curves
#' @param show.test show single or the reference group value (for pairwise
#'   comparisons). If `"none"`, then no test is show.
#' @param single.test.type test to show if specified `show.test =
#'   "single"`. Possible choices are `"logrank"` (default),
#'   `"wilcoxon"`, `"taroneware"`, or `"all"`.
#' @param round.digits.p.value number of digits for p-value
#' @param obs.survyrs show the observed survival years survival rate on KM plot
#' @param ten.years.surv.95CI show ten year survival 95\% confidence interval
#' @param event.count show the number of events at each time point
#' @param legend.pos legend position keyword
#' @param file.name name of file to save plot to
#' @param file.width width of figure in saved file
#' @param file.height height of figure in saved file
#' @param grey.scale logical. If `TRUE`, the plot will be in grey scale.
#' @param show.single.test.pos position to show single test; defaults to 0.5 if
#'   `legend.pos = "top"`. Otherwise 0.1
#' @param ... additional arguments to `plot`
#' @author Samuel Leung
#' @references
#'   http://courses.nus.edu.sg/course/stacar/internet/st3242/handouts/notes2.pdf
#' @seealso [plotKM()]
#' @export
plotKMDetail <- function(input.data, surv.formula,
                         main.text = "", xlab.text = "", ylab.text = "",
                         line.name, line.color, line.pattern = NULL,
                         line.width = NULL, show.test = "single",
                         single.test.type = "logrank",
                         round.digits.p.value = 4,
                         obs.survyrs, ten.years.surv.95CI, event.count,
                         legend.pos = "bottomleft",
                         file.name = "no.file",
                         file.width = 7, file.height = 7,
                         grey.scale = FALSE, show.single.test.pos, ...) {

  var.name <- deparse(surv.formula[[3]]) # this should be the biomarker name
  # the deparse() function is used to make sure var.name is a string

  log.rank.p.values    <- c()
  wilcox.p.values      <- c()
  tarone.ware.p.values <- c()

  fit <- survival::survfit(surv.formula, data = input.data)
  # do not generate a file if "no.file" is specified
  if (file.name != "no.file" & nchar(file.name) > 4) {
    file.ext <- tools::file_ext(file.name)
    if (file.ext == "pdf") {
      grDevices::cairo_pdf(filename = file.name, width = file.width, height = file.height)  # good for unicode character in e.g. line.name
    } else if (file.ext %in% c("wmf", "emf", "wmz", "emz")) {
      grDevices::png(filename = file.name, width = file.width, height = file.height)
    } else if (file.ext == "tiff") {
      grDevices::tiff(filename = file.name, width = file.width * 100, height = file.height * 100)
    } else {
      stop("Extension must be one of 'pdf', 'wmf', 'emf', 'wmz', 'emz', and 'tiff'.")
    }
  }
  # in case some strata do not have any cases
  which.strata.have.cases <- table(input.data[, var.name]) > 0
  # default line width
  if (is.null(line.width)) {
    line.width <- 1
  }
  if (grey.scale) {
    # gray scale plot
    if (is.null(line.pattern)) {
      line.pattern <- c(1:length(line.name))[which.strata.have.cases]
    }
    graphics::plot(fit, lty = line.pattern, lwd = line.width, main = main.text,
                   xlab = xlab.text, ylab = ylab.text, ...)
  } else {
    # color plot
    if (is.null(line.pattern)) {
      line.pattern <- 1
    }
    graphics::plot(fit, col = line.color[which.strata.have.cases],
                   lty = line.pattern, lwd = line.width, main = main.text,
                   xlab = xlab.text, ylab = ylab.text, ...)
  }
  # Legend 1
  if (legend.pos == "top") {
    x.pos <- diff(range(fit$time, na.rm = TRUE)) / 2
    y.pos <- 0.99   # top 1% ... since survival plot always starts at 100% survival
  } else {
    x.pos <- legend.pos
    y.pos <- NULL
  }
  l1 <- graphics::legend(
    x = x.pos, y = y.pos, legend = line.name,
    lty = line.pattern, lwd = line.width, box.lty = 0, cex = 0.8
  )

  # there seems to be need for the y-axis adjustment depending on the file.height ...
  dy <- 0.02 * (file.height - 7) / (12 - 7) # determined empirically
  if (legend.pos == "top") {
    y.pos <- l1$rect$top + dy
  } else {
    y.pos <- l1$rect$h - dy
  }
  # Legend 2 & 3
  l2 <- graphics::legend(
    x = l1$rect$w + l1$rect$left, y = y.pos,
    legend = ten.years.surv.95CI,
    title = paste0(obs.survyrs, "yr 95% CI"), title.col = 1,
    box.lty = 0, cex = 0.8
  )
  l3 <- graphics::legend(
    x = l1$rect$w + l1$rect$left + l2$rect$w, y = y.pos,
    legend = event.count,
    title = "Events/N", title.col = 1,
    box.lty = 0, cex = 0.8
  )
  graphics::box()
  if (show.test == "single") {
    log.rank.test     <- survival::survdiff(surv.formula, data = input.data,
                                            rho = 0)
    gehan.wilcox.test <- survival::survdiff(surv.formula, data = input.data,
                                            rho = 1)
    tarone.ware.test  <- survival::survdiff(surv.formula, data = input.data,
                                            rho = 0.5)
    p.value <- getPval(log.rank.test)
    log.rank.p.values <- p.value
    p.value <- round(p.value, digits = round.digits.p.value)
    gehan.wilcox.p.value <- getPval(gehan.wilcox.test)
    wilcox.p.values <- gehan.wilcox.p.value
    gehan.wilcox.p.value <- round(gehan.wilcox.p.value, digits = round.digits.p.value)
    tarone.ware.p.value <- getPval(tarone.ware.test)
    tarone.ware.p.values <- tarone.ware.p.value
    tarone.ware.p.value <- round(tarone.ware.p.value, digits = round.digits.p.value)
    graphics::text(
      x = l1$rect$w + l1$rect$left + l2$rect$w + 1.3 * l3$rect$w,
      y = show.single.test.pos, # position of the test statistics on plot
      paste0(
        ifelse(sum(single.test.type %in% c("logrank",   "all")) >= 1,
               paste0("Log-Rank p=",
                      sprintf(paste0("%.", round.digits.p.value, "f"),
                              p.value), "\n"), ""),
        ifelse(sum(single.test.type %in% c("wilcoxon",  "all")) >= 1,
               paste0("Wilcoxon p=",
                      sprintf(paste0("%.", round.digits.p.value, "f"),
                              gehan.wilcox.p.value), "\n"), ""),
        ifelse(sum(single.test.type %in% c("taroneware", "all")) >= 1,
               paste0("Tarone-Ware p=",
                      sprintf(paste0("%.", round.digits.p.value, "f"),
                              tarone.ware.p.value), "\n"), "")),
      adj = c(0, 0),
      cex = 0.8)
  } else if (show.test != "none") {
    # assume show.test shows the reference group index
    legend.txt <- c()
    value.names <- names(table(input.data[, var.name]))
    for (value.name in value.names) {
      if (value.name == show.test) {
        # this is the reference group
        legend.txt <- c(legend.txt, "reference group")
      } else {
        # construct data
        temp.d <- input.data[input.data[, var.name] == show.test |
                               input.data[, var.name] == value.name, ]
        if (sum(input.data[, var.name] == value.name, na.rm = TRUE) == 0) {
          # no case in this group
          p.value <- NA
          w.p.value <- NA
          t.p.value <- NA
        } else {
          # calculate log rank p-values
          p.value   <- getPval(survival::survdiff(surv.formula, data = temp.d, rho = 0))
          log.rank.p.values    <- c(log.rank.p.values,    p.value)
          p.value   <- round(p.value,  digits = round.digits.p.value)

          w.p.value <- getPval(survival::survdiff(surv.formula, data = temp.d, rho = 1))
          wilcox.p.values      <- c(wilcox.p.values,      w.p.value)
          w.p.value <- round(w.p.value, digits = round.digits.p.value)

          t.p.value <- getPval(survival::survdiff(surv.formula, data = temp.d, rho = 0.5))
          tarone.ware.p.values <- c(tarone.ware.p.values, t.p.value)
          t.p.value <- round(t.p.value, digits = round.digits.p.value)
        }
        new.txt <- paste0(
          ifelse("logrank"    %in% single.test.type, paste0(p.value,  " / "), ""),
          ifelse("wilcoxon"   %in% single.test.type, paste0(w.p.value, " / "), ""),
          ifelse("taroneware" %in% single.test.type,       t.p.value,     ""))
        if (endsWith(new.txt, " / ")) {
          new.txt <- substr(new.txt, 0, nchar(new.txt) - 3)
        }
        legend.txt <- c(legend.txt, new.txt)
      }
    }

    legend.title <- paste0(
      ifelse("logrank"    %in% single.test.type, "Log-Rank / ",  ""),
      ifelse("wilcoxon"   %in% single.test.type, "Wilcoxon / ",  ""),
      ifelse("taroneware" %in% single.test.type, "Tarone-Ware ", ""))
    if (endsWith(legend.title, " / ")) {
      legend.title <- substr(legend.title, 0, nchar(legend.title) - 2)
    }
    legend.title <- paste0(legend.title, "P-values")

    l4 <- graphics::legend(
      x = l1$rect$w + l2$rect$w + l3$rect$w, y = y.pos, #y=l1$rect$h,
      legend = legend.txt,
      #text.col=line.color,
      title = legend.title,
      title.col = 1,
      box.lty = 0,
      cex = 0.8
    )
  }
  if (file.name != "no.file") {
    # do not generate a file if "no.file" is specified
    grDevices::dev.off()
  }
  return(list(
    "log.rank.p.values" = log.rank.p.values,
    "wilcox.p.values" = wilcox.p.values
  ))
}
