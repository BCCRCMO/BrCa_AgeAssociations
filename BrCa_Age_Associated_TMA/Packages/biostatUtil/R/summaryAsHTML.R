#' Summary table in HTML format
#'
#' Generate summary table as an HTML table
#'
#' @param d assume `d` is an array of numbers
#' @return summary table with annotated HTML code
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' library(htmlTable)
#' set.seed(1)
#' x <- rnorm(100)
#' htmlTable(summaryAsHTML(x))
summaryAsHTML <- function(d) {
  col.th.style <- COL.TH.STYLE
  s.table <- summary(d)
  result <- paste(tags$table(
    tags$tr(HTML(paste(lapply(names(s.table), function(x)
      paste(tags$th(style = col.th.style, x))), collapse = "")
    )),
    tags$tr(HTML(paste(lapply(s.table, function(x)
      paste(tags$td(x))), collapse = "")
    ))
  ))
  return(result)
}
