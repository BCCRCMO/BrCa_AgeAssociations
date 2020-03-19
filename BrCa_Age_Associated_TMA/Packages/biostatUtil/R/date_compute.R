#' Date computations
#'
#' `addToDate` adds a period of time to a date string and returns a new
#' date string. `diffDate` takes the difference between two dates and
#' returns in the specified unit of time. `compareDate` compares the
#' difference between two dates.
#'
#' If `delta` is negative, then the returned date will be earlier than
#' `org.date`. The output date format will be the same as the input date
#' format.
#'
#' The unit of time to add to or return in can be in `"days"`,
#' `"weeks"`, `"months"`, or `"years"`. `compareDate` calls
#' `diffDate` and returns integer values specifying which date is earlier
#' (or if they are the same). `d1` should be later than `d2` so the
#' function returns nonnegative values.
#'
#' @param org.date original date
#' @param delta amount of time to add to `org.date`
#' @param date.format how to format the resulting date
#' @param units for `addToDate`, the unit of time for `delta`; for
#'   `diffDate`, the unit of time for which to take the difference.
#'   Defaults to "days".
#' @param existing.missing.codes missing dates
#' @param return.missing.code what to return if there is a missing input
#' @param sep date separator string
#' @return `addDate` returns the new date after adding `delta` to
#'   `org.date`.
#' @name date_compute
#' @family date computation functions
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' ## Adding to a date
#' addToDate("2014/07/08", 10, date.format = "YYYY.MM.DD")
#' addToDate("2014-07-08", 10, date.format = "YYYY.MM.DD", sep = "-")
#' addToDate("2014/07/08", 10, date.format = "YYYY.MM.DD", units = "months")
#' addToDate("2014/07/08", -10, date.format = "YYYY.MM.DD", units = "years")
#'
#' ## Date differences
#' # Later date comes first, subtracts the second date
#' diffDate("2003/03/21", "1992/01/27", date.format = "YYYY.MM.DD")
#'
#' # Otherwise negative
#' diffDate("1992/01/27", "2003/03/21", date.format = "YYYY.MM.DD")
#'
#' # Different separator
#' diffDate("2003-03-21", "1992-01-27", date.format = "YYYY.MM.DD", sep = "-")
#'
#' ## Date comparisons
#' compareDate("01/22/1949", "04/13/1950", date.format = "MM.DD.YYYY")
#' compareDate("04/13/1950", "04/13/1950", date.format = "MM.DD.YYYY")
#' compareDate("04/13/1959", "04/13/1950", date.format = "MM.DD.YYYY")
#' compareDate("01-22-1949", "04-13-1950", date.format = "MM.DD.YYYY", sep =
#' "-")
addToDate <- function(org.date, delta, date.format = "MM.DD.YYYY",
                      units = c("days", "weeks", "months", "years"),
                      existing.missing.codes = NA, return.missing.code = NA,
                      sep = "/") {
  if (is.na(org.date) | is.na(delta)) return(NA)
  if (length(unique(existing.missing.codes
                    [!is.na(existing.missing.codes)])) > 0 &
      (org.date %in% existing.missing.codes |
       delta %in% existing.missing.codes))
    return(return.missing.code)
  units <- match.arg(units)
  delta.time <- delta %>%
    as.numeric() %>%
    switch(units,
           days = .,
           weeks = . * 7,
           months = . * NUM.DAYS.IN.MONTH,
           years = . * NUM.DAYS.IN.YEAR)
  result <- cleanDate(org.date, date.format, date.format,
                      existing.missing.codes = existing.missing.codes,
                      return.missing.code = return.missing.code, sep = sep) %>%
    as.Date(format = getFormat(org.date, date.format), origin = DATE.ORIGIN) %>%
    magrittr::add(delta.time) %>%
    as.Date(origin = DATE.ORIGIN) %>%
    format(format = getFormat(org.date, date.format))
  return(result)
}

#' @param d1 later date
#' @param d2 earlier date
#' @inheritParams addToDate
#' @return `diffDate` returns the difference between two dates `d1 -
#'   d2` in the specified unit of time.
#' @rdname date_compute
#' @family date computation functions
#' @export
diffDate <- function(d1, d2, date.format = "MM.DD.YYYY",
                     units = c("days", "weeks", "months", "years"),
                     existing.missing.codes = NA, return.missing.code = NA,
                     sep = "/") {
  if (is.na(d1) | is.na(d2)) return(NA)
  if (n_distinct(existing.missing.codes, na.rm = TRUE) > 0 &
      any(c(d1, d2) %in% existing.missing.codes))
    return(return.missing.code)
  result <- c(d1, d2) %>%
    unname() %>%
    purrr::map(~ strptime(
      cleanDate(.x, date.format, date.format,
                existing.missing.codes = existing.missing.codes,
                return.missing.code = return.missing.code, sep = sep),
      format = getFormat(.x, date.format))) %>%
    purrr::invoke(difftime, .) %>%
    as.numeric()
  switch(match.arg(units),
         days = result,
         weeks = result / 7,
         months = result / NUM.DAYS.IN.MONTH,
         years = result / NUM.DAYS.IN.YEAR)
}

#' @inheritParams diffDate
#' @return `compareDate` returns 1 if `d1 > d2`, -1 if `d1 < d2`,
#'   and 0 if `d1 == d2`.
#' @rdname date_compute
#' @family date computation functions
#' @export
compareDate <- function(d1, d2, date.format = "MM.DD.YYYY",
                        existing.missing.codes = NA,
                        return.missing.code = NA, sep = "/") {
  difference <- diffDate(d1, d2, date.format = date.format, units = "days",
                         existing.missing.codes = existing.missing.codes,
                         return.missing.code = return.missing.code, sep = sep)
  if (is.na(difference))
    return(return.missing.code)
  if (difference > 0) {
    return(1)
  } else if (difference < 0) {
    return(-1)
  } else {
    return(0)
  }
}
