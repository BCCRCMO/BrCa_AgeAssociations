#' Fit a full linear model with subsetting at age of interest
#'
#' @param data data object
#' @param x age variable in `data`
#' @param y expression variable in `data`
#' @param age where to cut `x` for subsetted age analysis
#' @param fc fold change
#' @param fdr false discovery rate
#' @param main plot title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param xlim.pval horizontal limits for p-values
#' @param ylim.pval vertical limits for p-values
#' @param annotate vector of genes with an age-dependent trend that can be
#'   annotated
#' @param xlim.trend horizontal limits for age-dependent trend annotation
#' @param ylim.trend vertical limits for age-dependent trend annotation
#' @param conf.level confidence level
#' @param show.avg show average?
#' @param erbtxt ER binding annotation text
#'
#' @return linear model fit
#' @author Derek Chiu, Steven McKinney
#' @export
regression_plot <- function(data, x, y, age, fc = 1.25, fdr = 0.01,
                            main = NULL, xlab = NULL, ylab = NULL, xlim = NULL,
                            ylim = NULL, xlim.pval = NULL, ylim.pval = NULL,
                            annotate = NULL, xlim.trend = NULL,
                            ylim.trend = NULL, conf.level = NULL,
                            show.avg = FALSE, erbtxt = NULL) {
  # Note: Function has two expressions of error rate (fdr and conf.level).
  # FDR is only to place in expressions of age association to indicate the
  # FDR implemented across all the expression probes outside of this function.
  # conf.level is to control width of added confidence intervals within this function.
  #
  # Specify subsets in linear model
  s <- setNames(list(data[[x]] <= age, data[[x]] > age, NULL),
                c(glue::glue("[<={age}] p = "),
                  glue::glue("[>{age}] p = "),
                  "[All] p = "))

  # Set defaults for y-axis labels and limits based on range of expression
  if (min(data[[y]], na.rm = TRUE) < 0) {
    ylab <- ylab %||% "Normalized Expression (-4, 4)"
    ylim <- ylim %||% c(-5, 6)
    ylim.pval <- ylim.pval %||% c(5.6, 5.0, 4.4)
  } else {
    ylab <- ylab %||% "Raw Expression (4, 16)"
    ylim <- ylim %||% c(4, 17)
    ylim.pval <- ylim.pval %||% c(16.3, 15.5, 14.7)
  }
  ylim.trend <- ylim.trend %||% head(ylim.pval, 2)

  # Set defaults for x-axis labels and limits and erbtxt
  xlab <- xlab %||% "Age at diagnosis"
  xlim <- xlim %||% c(20, 100)
  xlim.pval <- xlim.pval %||% rep(20, 3)
  xlim.trend <- xlim.trend %||% 98
  main <- main %||% ""
  erbtxt <- erbtxt %||% ""

  # Extract lm p-values
  pval <- purrr::map(s, ~ lm(data[[y]] ~ data[[x]], subset = .) %>%
                       broom::tidy() %>%
                       magrittr::extract2(2, "p.value") %>%
                       format(digits = 5)) %>%
    paste0(names(.), .)

  # Run loess regression while keeping standard error metadata
  plx <- predict(loess(data[[y]] ~ data[[x]]), se = TRUE)

  # Plot the title, axis labels/limits first
  plot(data[[x]], data[[y]], main = glue::glue("{main} {y} {erbtxt}"),
       xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, type = "n")

  # Store order of x
  xord <- order(data[[x]])

  # Pointwise confidence interval bands for the loess regression
  # 1 - conf.level = Type I error rate
  if (!is.null(conf.level)) {
    me <- qt(1 - (1 - conf.level) / 2, plx$df) * plx$se.fit
    lower <- plx$fit - me
    upper <- plx$fit + me
    polygon(
      c(rev(data[[x]][xord]), (data[[x]][xord])),
      c(rev(upper[xord]), lower[xord]),
      col = col2rgba(col = "wheat3", alpha = 128),
      border = NA
    )
  }

  # Add scatterplot with alpha shading
  points(data[[x]], data[[y]], pch = 20, col = col2rgba(col = "black"))

  # Loess fitted curve
  lines(data[[x]][xord], plx$fit[xord], lwd = 3, col = "#AA1010")

  # Vertical line showing where age is dichotomized
  abline(v = age, lty = 2)

  # Placement of p-values
  text(x = xlim.pval, y = ylim.pval, labels = pval, adj = 0, cex = 0.85)

  # Whether the response variable (probe) had a significant trend in the
  # age-related data output "aroutdf"
  if (!is.null(annotate)) {
    targs <- list(x = xlim.trend, y = ylim.trend, adj = 1, cex = 0.85)
    probe_id <- ifelse(grepl("\\|", y), strsplit(y, split = "\\|")[[1]][2], y)
    if (probe_id %in% annotate) {
      labs <- c("Age-dependent trend:",
                glue::glue("|FC|>{fc} and adjPval<{fdr}"))
    } else {
      labs <- c("No detectable trend:",
                glue::glue("|FC|<{fc} or adjPval>{fdr}"))
    }
    purrr::invoke(text, targs, labels = labs)
  }

  # Horizontal line showing the mean expression
  if (show.avg) {
    abline(h = mean(data[[y]]), lty = 3, col = "grey20")
  }
}

#' Check mean enclosure
#'
#' Is the mean expression fully enclosed by the loess confidence bands?
#'
#' @inheritParams regression_plot
#' @return Returns `TRUE` only if for every set of upper and lower pointwise
#'   confidence limits, the mean expression is between these bounds, `FALSE` if
#'   the mean expression falls outside of the confidence limits for at least one
#'   point.
#' @author Derek Chiu
#' @export
#' @examples
#' is_mean_enclosed(mtcars, "wt", "qsec", conf.level = 0.99)
#' is_mean_enclosed(mtcars, "wt", "qsec", conf.level = 0.95)
is_mean_enclosed <- function(data, x, y, conf.level = 0.99) {
  plx <- predict(loess(data[[y]] ~ data[[x]]), se = TRUE)
  xord <- order(data[[x]])
  me <- qt(1 - (1 - conf.level) / 2, plx$df) * plx$se.fit
  lower <- plx$fit - me
  upper <- plx$fit + me
  all(purrr::map2_lgl(lower[xord],
                      upper[xord],
                      dplyr::between,
                      x = mean(data[[y]])))
}
