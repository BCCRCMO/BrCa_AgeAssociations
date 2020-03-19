#' Plot Principal Components
#'
#' Plot the first three PCs and draw ellipses highlighting
#' differences in the levels of a factor
#' @param by factor to plot as ellipses on PCA plots
#' @param tdat data matrix to compute principal components
#' @return PCA plots for every combination of PC1, PC2, and PC3. The
#' percentage of variation contribution is shown in the axes labels.
#' @author Aline Talhouk, Derek Chiu
#' @export
#'
#' @examples
#' PCAplot(mtcars$cyl, mtcars)
#' PCAplot(mtcars$gear, mtcars)
PCAplot <- function(by, tdat) {
  var.name <- NULL
  p <- stats::prcomp(tdat, retx = TRUE, scale. = TRUE)
  var3 <- round(((p$sdev ^ 2 / sum(p$sdev ^ 2)) * 100)[seq_len(3)], 2)
  df <- data.frame(p$x[, seq_len(3)], var.name = factor(by))
  PCA_labs <- purrr::map_chr(seq_len(3),
                             ~ paste0("PC", .x, " (", var3[.x], "% Var)"))
  PCA_list <- purrr::map2(c(1, 1, 2), c(2, 3, 3),
                          ~ ggplot(df, aes(df[, .x], df[, .y],
                                           colour = var.name)) +
                            geom_point() +
                            stat_ellipse(geom = "polygon", level = 0.9,
                                         alpha = 0.1, size = 0.05,
                                         aes(fill = var.name)) +
                            theme_bw() +
                            theme(legend.title = element_blank()) +
                            labs(x = PCA_labs[.x], y = PCA_labs[.y]))
  multiplot(plotlist = PCA_list, cols = 1)
}
