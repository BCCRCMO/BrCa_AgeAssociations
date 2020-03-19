#' Multiple plots
#'
#' Place multiple ggplot objects on the same figure space
#'
#' If the `layout` is something like `matrix(c(1, 2, 3, 3), nrow = 2,
#' byrow = TRUE)`, then plot 1 will go in the upper left, 2 will go in the upper
#' right, and 3 will go all the way across the bottom.
#'
#' @param ... pass ggplot objects
#' @param plotlist pass a list of ggplot objects
#' @param cols number of columns in layout
#' @param layout a matrix specifying the layout. If present, `cols` is
#'   ignored.
#' @return A grid object made up of multiple ggplots
#' @author Aline Talhouk
#' @export
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols, nrow = ceiling(numPlots / cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location

    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}
