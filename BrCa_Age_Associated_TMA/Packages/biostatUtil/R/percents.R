#' Row and Column Percentages
#'
#' Calculate percentages in a table. `rowPercent` gives row percentages,
#' `colPercent` gives column percentages, and `rowColPercent` gives
#' both row and column percentages.
#'
#' Generates a table of row and/or column percentages given table `t`.
#' Using `pretty.text = TRUE` will add the \% sign to the percentages.
#'
#' @param t a matrix
#' @param pretty.text logical; if `TRUE`, will format the table into nice
#'   display
#' @param keep logical; if `TRUE`, the original table counts will be kept
#'   along with the percentages
#' @param digits number of digits to round to
#' @return A table with row-wise/column-wise percentages added. The percentages
#'   sum to 1 per row/column.
#' @author Aline Talhouk, Samuel Leung, Derek Chiu
#' @name percents
#' @rdname percents
#' @export
#' @examples
#' # Base outputs
#' A <- matrix(c(2, 3, 5, 10), nrow = 2, dimnames = list(c("Row1", "Row2"), c("Col1", "Col2")))
#' rowPercent(A)
#' rowPercent(A, keep = FALSE)
#' colPercent(A, pretty.text = TRUE)
#' colPercent(A, pretty.text = TRUE, keep = FALSE)
#' rowColPercent(A, digits = 2)
colPercent <- function(t, pretty.text = FALSE, keep = TRUE, digits = 4) {
  t <- as.matrix(t)
  if (is.null(rownames(t)))
    rownames(t) <- seq_len(nrow(t))
  pcts <- round(t / colSums(t)[col(t)], digits = digits)
  if (pretty.text)
    pcts <- apply(pcts * 100, 1:2,
                  function(x) ifelse(!is.nan(x), paste0(x, "%"), "-"))
  if (keep) {
    pcts <- gdata::interleave(t, pcts)
    rownames(pcts) <- paste0(rownames(pcts), rep(c("", " Col %"), nrow(t)))
  }
  return(pcts)
}

#' @rdname percents
#' @export
rowPercent <- function(t, pretty.text = FALSE, keep = TRUE, digits = 4) {
  t <- as.matrix(t)
  if (is.null(rownames(t)))
    rownames(t) <- seq_len(nrow(t))
  pcts <- round(t / apply(t, 1, sum), digits = digits)
  if (pretty.text)
    pcts <- apply(pcts * 100, 1:2,
                  function(x) ifelse(!is.nan(x), paste0(x, "%"), "-"))
  if (keep) {
    pcts <- gdata::interleave(t, pcts)
    rownames(pcts) <- paste0(rownames(pcts), rep(c("", " Row %"), nrow(t)))
  }
  return(pcts)
}

#' @rdname percents
#' @param ... additional arguments from `colPercent` and `rowPercent`
#' for `rowColPercent`, or additional arguments from non-HTML functions
#' to HTML functions.
#' @export
rowColPercent <- function(t, keep = TRUE, ...) {
  if (is.null(rownames(t)))
    rownames(t) <- seq_len(nrow(t))
  row.p <- rowPercent(t, keep = !keep, ...)
  col.p <- colPercent(t, keep = !keep, ...)
  if (keep) {
    result <- as.matrix(gdata::interleave(t, row.p, col.p)) %>%
      magrittr::set_rownames(paste0(rownames(.), rep(c("", " Row %", " Col %"), nrow(t))))
  } else {
    result <- as.matrix(gdata::interleave(row.p, col.p)) %>%
      magrittr::extract(grep("%", rownames(.)), )
  }
  return(result)
}
