#' Add table number
#'
#' Prepends a table counter label to a table caption.
#'
#' The table counter automatically increments based on the number of existing
#' tables that have already called  using `addTableNumber`.
#'
#' @param caption table caption
#' @param table.counter.str.default the default string preceding the table count
#' @return A caption with a table counter prepended to its text.
#' @author Samuel Leung
#' @export
#' @examples
#' addTableNumber("Overall Survival Hazard Ratios")
#' addTableNumber("Correlation Heatmap", table.counter.str.default = "Figure %s:
#' ")
addTableNumber <- function(caption, table.counter.str.default = "Table %s: ") {
  if (!is.na(caption)) {
    tc <- getOption("table_counter")
    if (is.numeric(tc) & length(tc) > 0)
      tc <- tc + 1
    else
      tc <- 1
    options(table_counter = tc)
    caption_template <- getOption("table_counter_str",
                                  table.counter.str.default)
    caption <- paste0(sprintf(caption_template, as.character(tc)), caption)
  }
  return(caption)
}
