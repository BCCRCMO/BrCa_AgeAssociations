#' Barplot with counts
#'
#' Generates barplot on variables and shows counts per group.
#'
#' The number of data points in each group is reported.
#' The number of scorable and missing data are also reported by number
#' and percentage in the barplot title.
#'
#' @param data data.frame with column names
#' @param var string of variable in `data` to graph on
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param title barplot title
#' @param digits number of digits to round to
#' @return Barplot of frequencies for each group
#' @author Samuel Leung, Derek Chiu
#' @note Function expects missing to be `NA`. Do not filter out missing data
#' as this function reports missing data frequencies.
#' @export
#' @examples
#' doBarplot(mtcars, "cyl", "title = Number of cylinders")
#' doBarplot(mtcars, "cyl", "Cylinders", title = "Number of cylinders")
doBarplot <- function(data, var, xlab = var, ylab = "Frequency", title = NULL,
                      digits = 3) {
  dat.var <- data[, var]
  counts <- apply(cbind(names(table(dat.var)), as.numeric(table(dat.var))),
                  1, function(x) paste0(x, collapse = "\nn="))
  n.s <- sum(!is.na(dat.var))
  n.m <- sum(is.na(dat.var))
  graphics::barplot(
    table(dat.var), names.arg = counts, xlab = xlab, ylab = ylab,
    main = paste0(title, "\n# scorable, missing: ", n.s, "(",
                  format(n.s / nrow(data) * 100, digits = digits), "%)",
                  ", ", n.m, "(",
                  format(n.m / nrow(data) * 100, digits = digits), "%)")
  )
}

#' Do a boxplot among subtypes
#'
#' A stripchart (jitter plot) of a biomarker and some categorical subtype is
#' superimposed on top of the boxplot
#' @param input.d input data.frame
#' @param data.description boxplot title description
#' @param biomarker.var.name biomarker variable name in `input.d` to graph
#'   on
#' @param biomarker.name boxplot y-axis label
#' @param subtype.var.name subtype variable name in `input.d` to graph on
#' @param subtype.name boxplot x-axis label
#' @param pch stripchart plot style
#' @param jitter amount of jitter in stripchart
#' @param digits number of digits to round to
#' @param ... additional arguments to `boxplot`
#' @return a boxplot of biomarker values across different subtypes
#' @note Expects subtype variable to be a factor. Also expects biomarker and
#'   subtype variable missing values to be `NA`
#' @author Samuel Leung
#' @export
#' @examples
#' doBoxplotAmongSubtypes(mtcars, "Boxplot of qsec vs. gear", "qsec", "QSEC",
#' "gear", "GEAR")
doBoxplotAmongSubtypes <- function(input.d, data.description,
                                   biomarker.var.name, biomarker.name,
                                   subtype.var.name, subtype.name, pch = 4,
                                   jitter = 0.1, digits = 2, ...) {
  temp.d <- input.d[(!is.na(input.d[, biomarker.var.name])) &
                      (!is.na(input.d[, subtype.var.name])), ]
  biomarker <- temp.d[, biomarker.var.name]
  if (is.factor(temp.d[, subtype.var.name])) {
    subtype <- droplevels(temp.d[, subtype.var.name])
  } else {
    subtype <- temp.d[, subtype.var.name]
  }
  xbar <- tapply(biomarker, subtype, bootMean)
  test.name <- "Kruskal-Wallis"
  p.value <- stats::kruskal.test(biomarker ~ subtype)$p.value
  if (length(names(table(subtype))) == 2) {
    test.name <- "Wilcoxon Rank Sum"
    p.value <- stats::wilcox.test(biomarker ~ subtype)$p.value
  }

  graphics::boxplot(
    biomarker ~ subtype,
    names = paste0(paste0(names(xbar), rep("\nn=", length(xbar))),
                   purrr::map_int(xbar, "n")),
    ylab = biomarker.name, xlab = subtype.name,
    main = paste0(data.description, "\n", test.name, " test P=",
                  format(p.value, digits = digits)),
    outline = FALSE, ...
  )
  graphics::stripchart(
    jitter(biomarker) ~ subtype, vertical = TRUE, method = "jitter",
    pch = pch, add = TRUE, jitter = jitter
  )
}

#' Do histogram with median
#'
#' Plot a histogram with the median, first quartile, and third quartile reported.
#'
#' Expects missing to be `NA`. Do not filter out missing data
#' as this function reports missing data counts.
#'
#' @param data data.frame with column names
#' @param var string of variable in `data` to graph on
#' @param xlab x-axis label
#' @param title histogram title
#' @param show.title logical. If `TRUE` (default), the title is shown
#' @param br breaks in histogram. By default, the number of bins is taking as
#' the ceiling of the square root of the number of rows in `data`
#' @param digits number of digits to round for median, Q1, Q3
#' @param score.lab label for non-missing cases
#' @param ... additional arguments to `hist`
#' @return A histogram with some annotated values.
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' doHist(mtcars, "mpg")
#' doHist(mtcars, "mpg", "MPG")
#' doHist(mtcars, "mpg", title = "Distribution of MPG")
doHist <- function(data, var, xlab = var, title = NULL, show.title = TRUE,
                   br = ceiling(sqrt(nrow(data))), digits = 3,
                   score.lab = "scorable", ...) {
  dat.var <- data[, var]
  qs <- stats::quantile(dat.var, na.rm = TRUE)
  n.s <- sum(!is.na(dat.var))
  n.m <- sum(is.na(dat.var))
  if (show.title) {
    title <- ifelse(!is.null(title), paste0(title, "\n"), "")
    main <- paste0(
      title, "Mean (Min, Q1, Median, Q3, Max): ",
      format(mean(dat.var, na.rm = TRUE), digits = digits), " ", "(",
      paste(format(qs, digits = digits), collapse = ", "), ")",
      "\n# ", score.lab, ", missing: ", n.s, "(",
      format(n.s / nrow(data) * 100, digits = digits), "%)", ", ", n.m, "(",
      format(n.m / nrow(data) * 100, digits = digits), "%)")
  } else {
    main <- ""
  }
  graphics::hist(dat.var, br = br, xlab = xlab, main = main, ...)
}

#' Do a jitterplot among subtypes
#'
#' Generates stripchart (jitter plot) of a biomarker split on a categorical subtype.
#'
#' Expects subtype variable to be a factor. Also, biomarker and subtype variables
#' have missing cases as `NA`.
#'
#' @param input.d input `data.frame`
#' @param data.description title description
#' @param biomarker.var.name variable name of biomarker to plot
#' @param biomarker.name x-axis label for biomarker name
#' @param subtype.var.name variable name of subtype to separate biomarker
#' @param subtype.name y-axis label for subtype name
#' @param pch jitterplot dot type. For example, 20 = small dot, 16 = big dot.
#' @param jitter scalar indicating amount of jitter (spread of points)
#' @param digits number of digits to round to
#' @param cex.axis scalar indicating amount of scaling for axes
#' @return A jitterplot shown across different subtypes
#' @author Samuel Leung
#' @export
#' @examples
#' doJitterplotAmongSubtypes(mtcars, "Boxplot of qsec vs. vs", "qsec", "QSEC",
#' "vs", "VS")
doJitterplotAmongSubtypes <- function(input.d, data.description,
                                      biomarker.var.name, biomarker.name,
                                      subtype.var.name, subtype.name,
                                      pch = ".", jitter = 0.05, digits = 3,
                                      cex.axis = 0.9) {

  temp.d <- input.d[(!is.na(input.d[, biomarker.var.name])) &
                      (!is.na(input.d[, subtype.var.name])), ]
  biomarker <- temp.d[, biomarker.var.name]
  subtype <- temp.d[, subtype.var.name]
  xbar <- tapply(biomarker, subtype, bootMean)
  test.name <- "Kruskal-Wallis"
  p.value <- stats::kruskal.test(biomarker ~ subtype)$p.value
  if (length(names(table(subtype))) == 2) {
    test.name <- "Wilcoxon Rank Sum"
    p.value <- stats::wilcox.test(biomarker ~ subtype)$p.value
  }
  graphics::par(mar = c(5.1, 4.1, 5.1, 2.1))
  graphics::stripchart(
    biomarker ~ subtype, method = "jitter", jitter = jitter,
    pch = pch, cex.axis = cex.axis, vert = TRUE,
    group.names = paste0(
      paste0(names(xbar), rep("\nn=", length(xbar))),
      purrr::map_int(xbar, "n")),
    ylab = biomarker.name, xlab = subtype.name,
    main = paste0(data.description, "\n", test.name, " test P=",
                  format(p.value, digits = digits))
  )
  graphics::arrows(
    seq_along(xbar), purrr::map_dbl(xbar, c(2, 1)),
    seq_along(xbar), purrr::map_dbl(xbar, c(2, 2)),
    angle = 90, code = 3, length = 0.1
  )
  graphics::points(
    purrr::map_dbl(xbar, "obs.mean"), pch = 4, type = "p", cex = 2
  )
}
