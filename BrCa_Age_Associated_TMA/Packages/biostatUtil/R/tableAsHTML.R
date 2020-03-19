#' Table as HTML format
#'
#' Generate HTML code for a table output. Assumes table has more than one row.
#'
#' @inheritParams rowPercentAsHTML
#' @return A character string that can be parsed as HTML code to produce a nicer
#'   table output.
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' library(htmlTable)
#' A <- matrix(c(2, 3, 5, 10), nrow = 2, dimnames = list(c("Row1", "Row2"), c("Col1", "Col2")))
#' htmlTable(tableAsHTML(A, banded.rows = TRUE))
tableAsHTML <- function(
  t, row.names = NULL, column.names = NULL,
  html.table.border = 0, banded.rows = FALSE,
  col.odd = "none", col.even = "lightgrey", caption = NA) {

  th.style <- COL.TH.STYLE
  if (!is.null(row.names)) {
    rownames(t) <- row.names
  } else {
    row.names <- rownames(t)
  }
  if (!is.null(column.names)) {
    colnames(t) <- column.names
  } else {
    column.names <- colnames(t)
  }
  if (banded.rows) {
    row.col <- rep(paste0("background-color: ", c(col.odd, col.even)),
                   nrow(t) / 2)
  } else {
    row.col <- rep("", nrow(t))
  }
  result <- paste0(HTML(paste0(
    tags$caption(style = TABLE.CAPTION.STYLE,
                 ifelse(is.na(caption), "", addTableNumber(caption))),
    tags$tr(HTML(paste0(
      tags$th(style = th.style),
      paste(unlist(lapply(column.names, function(x)
        paste(tags$th(style = th.style, x)))), collapse = "")
    )))
  )))
  for (i in seq_len(nrow(t))) {
    result <- paste0(HTML(paste0(
      result,
      tags$tr(style = row.col[i], HTML(paste0(
        tags$th(row.names[i]),
        paste(unlist(lapply(t[i, ], function(x)
          paste(tags$td(x)))), collapse = "")
      )))
    )))
  }
  result <- paste0(tags$table(border = html.table.border, HTML(result)))
  return(result)
}
