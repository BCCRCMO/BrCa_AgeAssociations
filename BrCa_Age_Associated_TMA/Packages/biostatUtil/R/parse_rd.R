#' Parse all Rd files in package
#'
#' Parse all Rd files in a package and return metadata as a data frame.
#'
#' The function requires that the working directory is set to the package's root
#' directory.
#'
#' @param path file path to save table
#' @param tags sections to extract from Rd files. Missing entries labelled as
#'   `NA`
#'
#' @return a table with information on a package's documented objects (functions
#'   and data). Rows are objects, and columns are tags.
#' @author Derek Chiu
#' @export
#'
#' @examples
#' \dontrun{
#' tab <- parse_rd()
#' str(tab)
#' }
parse_rd <- function(path = NULL, tags = c("name", "title", "desc",
                                           "details")) {
  rd.files <- list.files("man", full.names = TRUE)
  info.table <- rd.files %>%
    purrr::set_names() %>%
    purrr::map(Rd2roxygen::parse_file) %>%
    purrr::map(~ replace(.x[tags], purrr::map_lgl(.x[tags], is.null), "")) %>%
    purrr::transpose() %>%
    purrr::map(unlist) %>%
    data.frame(stringsAsFactors = FALSE)
  if (!is.null(path)) readr::write_csv(info.table, path = path)
  info.table
}
