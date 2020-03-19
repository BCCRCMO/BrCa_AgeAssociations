#' Latest and Earliest Dates
#'
#' Finds the latest (max) and earliest (min) date given a vector of dates. The
#' functions ending in `Array` find the extremes for an array of dates as
#' opposed to a vector.
#'
#' The input vector should have dates formatted as YYYY-MM-DD. If `na.rm`
#' is not set to the default (`TRUE`) and `dates` has `NA`
#' values, then the function will also return `NA`.
#'
#' @param dates vector of dates
#' @param na.rm logical. Should NAs be removed?
#' @return The latest or earliest date from a vector of dates.
#' @author Samuel Leung, Derek Chiu
#' @name date_extremes
#' @export
#' @examples
#' ## No NA
#' t1 <- as.Date(c("2015-03-01", "2015-02-15", "2015-05-01"))
#' minDate(t1)
#' maxDate(t1)
#'
#' ## With NA
#' t2 <- as.Date(c("2015-03-01", "2015-02-15", NA, "2014-05-01"))
#' maxDate(t2)
#' minDate(t2)
#' minDate(t2, na.rm = FALSE)
maxDate <- function(dates, na.rm = TRUE) {
  return(max(dates, na.rm = na.rm))
}

#' @rdname date_extremes
#' @export
minDate <- function(dates, na.rm = TRUE) {
  return(min(dates, na.rm = na.rm))
}

#' @param t.arr array of date strings
#' @param date.format format of the array of dates
#' @param existing.missing.codes missing dates
#' @param return.missing.code what to return if there is a missing input
#' @param ... additional arguments to `formatDate`
#' @param sep date separator. Defaults to "/"
#' @rdname date_extremes
#' @export
#' @examples
#'
#' ## Array of dates
#' many.dates <- c("03/21/1992", "04/21/2013", "10/10/2015")
#' maxDateArray(many.dates)
#' minDateArray(many.dates)
#'
#' many.dates <- c("2009-03-01", "2010-01-12", "2015-01-11")
#' maxDateArray(many.dates, sep = "-")
#' minDateArray(many.dates, sep = "-")
#'
#' ties.dates <- c("2009-03-01", "2010-01-12", "2010-01-12")
#' maxDateArray(ties.dates, sep = "-")
#' minDateArray(ties.dates, sep = "-")
maxDateArray <- function(t.arr, date.format = "MM.DD.YYYY",
                         existing.missing.codes = NA,
                         return.missing.code = NA, sep = "/") {
  t.arr <- t.arr[!is.na(t.arr)]
  if (n_distinct(existing.missing.codes, na.rm = TRUE) > 0)
    t.arr <- t.arr[!(t.arr %in% existing.missing.codes)]
  if (length(t.arr) == 0)
    return(return.missing.code)
  max.index <- purrr::map(t.arr, ~ as.Date(
    cleanDate(.x, original.format = date.format, preferred.format = date.format,
              existing.missing.codes = existing.missing.codes,
              return.missing.code = return.missing.code, sep = sep),
    format = getFormat(.x, date.format), origin = DATE.ORIGIN)) %>%
    which.max()
  return(t.arr[max.index])
}

#' @rdname date_extremes
#' @export
minDateArray <- function(t.arr, date.format = "MM.DD.YYYY",
                         existing.missing.codes = NA,
                         return.missing.code = NA, sep = "/") {
  t.arr <- t.arr[!is.na(t.arr)]
  if (n_distinct(existing.missing.codes, na.rm = TRUE) > 0)
    t.arr <- t.arr[!(t.arr %in% existing.missing.codes)]
  if (length(t.arr) == 0)
    return(return.missing.code)
  min.index <- purrr::map(t.arr, ~ as.Date(
    cleanDate(.x, original.format = date.format, preferred.format = date.format,
              existing.missing.codes = existing.missing.codes,
              return.missing.code = return.missing.code, sep = sep),
    format = getFormat(.x, date.format), origin = DATE.ORIGIN)) %>%
    which.min()
  return(t.arr[min.index])
}
