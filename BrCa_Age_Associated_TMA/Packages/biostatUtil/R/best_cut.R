#' Find the best cutpoint for a covariate
#'
#' Repeatedly finds cutpoints for an explanatory variable in a univariable Cox
#' PH model. Also plots survival curves for each cutpoint.
#'
#' Takes the cutpoint resulting in the lowest AIC. If the range of AIC values is
#' within `AIC.range` units, take the cutpoint that results in the two groups
#' having the most similar numbers of events and cases. The function can cut a
#' variable into anywhere from 2 to 5 groups.
#'
#' @param f formula object
#' @param d data frame
#' @param n number of groups to transform variable into. Options are "b" (two),
#'   "t" (three), "qd" (four), and "qn" (five)
#' @param AIC.range If range of AIC is within `AIC.range` units, the
#'   likelihood is too flat. We choose the best cutpoint using the alternative
#'   method.
#' @param nround number of digits to round AIC and p-value on plots
#' @param plot logical; If `TRUE`, shows the survival curves for each
#'   cutpoint in a facetted figure
#' @param filename file name for saving a png image of figure
#' @param nrow number of rows in facetted plot
#' @param ncol number of columns in facetted plot
#' @param title title for plot
#' @param ... additional arguments for `plot`
#'
#' @return A list with the following elements
#' \item{cuts}{vector of cutpoints considered}
#' \item{fits}{A list of `coxph` objects run for each cutpoint}
#' \item{results}{A table showing the likelihood ratio test p-value,
#'   log likelihood, and AIC for each cutpoint}
#' \item{opt.cut}{optimal cutpoint value}
#' \item{flat.lik}{If `TRUE`, the likelihood was too flat and the
#'   alternative method was used}
#' Additionally, if `plot = TRUE`, the function also returns KM survival
#' curves for each possible cutpoint.
#'
#' @author Derek Chiu
#' @export
best_cut <- function(f, d, n = c("b", "t", "qd", "qn"), AIC.range = 3,
                     nround = 3, plot = TRUE, filename = NULL,
                     nrow = NULL, ncol = NULL, title = "", ...) {
  cutpoints <- p.value.log <- logLik <- AIC <- NULL
  pos <- 1
  assign("f", f, envir = as.environment(pos))
  assign("d", d, envir = as.environment(pos))

  # Build cutpoints, cox and km fits, summarize results
  bins <- build_cuts(d[, all.vars(f)[3]], n = n, list = TRUE)
  cuts <- stringr::str_extract_all(names(bins), ".v") %>%
    purrr::map_chr(~ paste(gsub("v", "", .x), collapse = ", "))
  coxs <- purrr::map(bins, ~
    coxph(stats::as.formula(paste(deparse(f[[2]]), "~ .x")), d))
  diffs <- purrr::map(bins, ~
    survfit(stats::as.formula(paste(deparse(f[[2]]), "~ .x")), d))
  results <- coxs %>%
    purrr::map_df(broom::glance, .id = "bin.names") %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$bin.names) %>%
    dplyr::select(p.value.log, logLik, AIC) %>%
    cbind(cutpoints = cuts, .)
  p.vals <- signif(results$p.value.log, nround)
  AIC.vals <- round(results$AIC, nround)

  # Check for flat likelihood issue using range of AIC
  if (diff(range(results$AIC)) < AIC.range) {
    # Cutpoint that distributes the group size and number of events most evenly
    opt.ind <- purrr::map_dbl(diffs, ~ prod(summary(.x)$table[, "events"] /
                                              sum(summary(.x)$n.event))) %>%
      which.max()
    flat.lik <- TRUE
  } else {
    opt.ind <- which.min(results$AIC)
    flat.lik <- FALSE
  }
  # Add annotation to title indicating best cutpoint
  opt.cut <- results$cutpoints[[opt.ind]]
  best.ind <- rep("", length(cuts)) %>%
    magrittr::inset(opt.ind, " (Best)")
  titles <- paste0(title, " ", names(bins), best.ind)

  # Plot survival curves for every cutpoint in console/PNG file
  if (plot) {
    if (is.null(filename)) {
      graphics::par(mfrow = c(nrow, ncol))
      purrr::pwalk(list(diffs, titles, p.vals, AIC.vals), best_cut_plot, ...)
      graphics::par(mfrow = c(1, 1))
    } else {
      grDevices::png(filename, width = 8.5, height = 11, units = "in", res = 300)
      graphics::par(mfrow = c(nrow %||% 2, ncol %||% 2))
      purrr::pwalk(list(diffs, titles, p.vals, AIC.vals), best_cut_plot, ...)
      graphics::par(mfrow = c(1, 1))
      grDevices::dev.off()
    }
  }
  return(list(cuts = cuts, fits = coxs, results = results, opt.cut = opt.cut,
              flat.lik = flat.lik))
}

#' Plotting function for best_cut
#' @noRd
best_cut_plot <- function(x, title, pval = NULL, aic = NULL, lwd = 1,
                          cex = 0.75, ...) {
  graphics::plot(x, main = title, col = seq_along(x$strata), lwd = lwd, ...)
  graphics::legend("bottomleft", legend = stringr::str_split_fixed(
    names(x$strata), "=", 2)[, 2], col = seq_along(x$strata),
    lwd = lwd, cex = cex)
  if (!is.null(pval))
    graphics::mtext(paste("P =", pval), side = 1, line = -2, at = max(x$time),
          adj = 1, cex = cex)
  if (!is.null(aic))
    graphics::mtext(paste("AIC:", aic), side = 1, line = -1, at = max(x$time),
          adj = 1, cex = cex)
}
