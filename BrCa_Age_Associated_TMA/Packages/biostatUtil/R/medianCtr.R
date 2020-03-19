#' Median center rows
#'
#' Takes a matrix and centers the rows by their median.
#'
#' The data matrix should have genes as rows and samples as columns. After
#' centering, each row should have median approximately 0.
#'
#' @param x a matrix
#' @return A matrix with median centered rows.
#' @author Samuel Leung
#' @export
#' @examples
#' x <- matrix(rnorm(200), nrow = 10)
#' medianCtr(x)
medianCtr <- function(x) {
  annAll <- dimnames(x)
  medians <- apply(x, 1, stats::median, na.rm = TRUE)
  x <- t(scale(t(x), center = medians, scale = FALSE))
  dimnames(x) <- annAll
  return(x)
}
