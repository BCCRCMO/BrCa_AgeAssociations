#' Scatter-fit plot
#'
#' Raw expression vs. age with smoother and annotations
#'
#' @param x data object
#' @param age age of interest
#' @param title plot title. Defaults to gene name and shows ER binding status if
#'   present
#' @param caption plot caption
#'
#' @return A scatterplot with fit trend line added
#' @author Derek Chiu
#' @export
scatterfit_plot <- function(x, age = 60, conf.level = 0.95, title = NULL, caption = NULL) {
  title <- purrr::`%||%`(title, x$Gene)
  if (!is.null(x$ERBS)) title <- paste(title, x$ERBS)

  p <- ggplot(x, aes_(~Age, ~Expression)) +
    geom_point(size = 0.5, alpha = 0.3) +
    labs(x = "Age at diagnosis", y = "log2(Raw expression)",
         title = title, caption = caption) +
    ylim(c(4, 17)) + ##ylim(c(4, 23)) +
    scale_x_continuous(breaks = pretty(x$Age, n = 4)) +
    stat_smooth(method = "loess", colour = "#AA1010", fill = "wheat3",
                se = TRUE, level = conf.level) +
    geom_vline(xintercept = age, linetype = "dashed") +
    geom_hline(yintercept = mean(x$Expression), linetype = "dotted") +
    annotate("text", label = unique(x[[paste0("LE", age, "_pval")]]),
             x = min(x$Age), y = 16.5, hjust = 0, size = 3) + ## y = 21
    annotate("text", label = unique(x[[paste0("GT", age, "_pval")]]),
             x = min(x$Age), y = 15.8, hjust = 0, size = 3) + ## y = 20
    annotate("text", label = unique(x$AllAges_pval),
             x = min(x$Age), y = 15.1, hjust = 0, size = 3) + ## y = 19
    annotate("text", label = unique(x$Trend),
             x = max(x$Age), y = 16.5, hjust = 1, vjust = 1, size = 3) + ## y = 21
    theme_linedraw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.caption = element_text(face = "bold.italic", hjust = 0.5))
  return(p)
}
