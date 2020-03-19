#' Change numeric to date
#'
#' Change a numeric value to a date object by specifying a date of origin.
#'
#' @param x a number that represents the number of days after `date.origin`
#' @param date.origin the date from which we count the number of days passed
#' @return A date object, converted from a numeric object.
#' @family date formatting functions
#' @author Samuel Leung
#' @export
#' @examples
#' numericToDate(10)
#' numericToDate(10, "2000-09-11")
numericToDate <- function(x, date.origin = DATE.ORIGIN) {
  return(as.Date(x, origin = date.origin))
}

#' Get date format from character text
#'
#' Get the POSIX standard date formats from character text formats.
#'
#' @param date character string of date
#' @param char.format character text format of date
#' @param sep character string separating `date`
#' @return A character string representing the POSIX standard date format
#'   equivalent of the string in `char.format`.
#' @family date formatting functions
#' @author Derek Chiu
#' @export
#' @examples
#' getFormat("12/09/1993", "MM.DD.YYYY")
#' getFormat("2005-09-13", "YYYY.MM.DD")
getFormat <- function(date, char.format = c("MM.DD.YYYY", "MMM.DD.YYYY",
                                            "DD.MM.YYYY", "DD.MMM.YYYY",
                                            "YYYY.MM.DD", "YYYY.MMM.DD"),
                      sep = "") {
  if (grepl("-", date))
    sep <- "-"
  else if (grepl("/", date))
    sep <- "/"
  else if (grepl("\\|", date))
    sep <- "|"
  else
    sep <- sep
  switch(match.arg(char.format),
         MM.DD.YYYY = paste0("%m", sep, "%d", sep, "%Y"),
         MMM.DD.YYYY = paste0("%b", sep, "%d", sep, "%Y"),
         DD.MM.YYYY = paste0("%d", sep, "%m", sep, "%Y"),
         DD.MMM.YYYY = paste0("%d", sep, "%b", sep, "%Y"),
         YYYY.MM.DD = paste0("%Y", sep, "%m", sep, "%d"),
         YYYY.MMM.DD = paste0("%Y", sep, "%b", sep, "%d"))
}

#' Format a date.
#'
#' Prints a date into a pretty format.
#'
#' Given the day, month, and year of a date, returns the date in a specific
#' format. The order and separating string can be modified using
#' `date.format` and `sep` respectively. Take note the order of the
#' arguments: day, month, and year. Only accepts "MM.DD.YYYY", "MMM.DD.YYYY",
#' "DD.MM.YYYY", "DD.MMM.YYYY", "YYYY.MM.DD", "YYYY.MMM.DD".
#'
#' @param d day of the month (1-31)
#' @param m month of the year (1-12)
#' @param y year of date
#' @param date.format how to format the date. Defaults to "month/day/year".
#' @param sep string used to separate `d`, `m`, and `y`. Defaults
#'   to "/".
#' @return A character string of a formatted date.
#' @family date formatting functions
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' formatDate(8, 7, 2011)
#' formatDate(8, 7, 2011, date.format = "YYYY.MM.DD")
#' formatDate(8, 7, 2011, date.format = "DD.MM.YYYY", sep = "-")
#' formatDate(10, 1, 2015, date.format = "MMM.DD.YYYY", sep = "-")
formatDate <- function(d, m, y, date.format = c("MM.DD.YYYY", "MMM.DD.YYYY",
                                                "DD.MM.YYYY", "DD.MMM.YYYY",
                                                "YYYY.MM.DD", "YYYY.MMM.DD"),
                       sep = "/") {
  d <- as.numeric(d)
  m <- as.numeric(m)
  y <- as.numeric(y)
  switch(
    match.arg(date.format),
    MM.DD.YYYY = paste(sprintf("%02d", m), sprintf("%02d", d), y, sep = sep),
    MMM.DD.YYYY = paste(month.abb[m], sprintf("%02d", d), y, sep = sep),
    DD.MM.YYYY = paste(sprintf("%02d", d), sprintf("%02d", m), y, sep = sep),
    DD.MMM.YYYY = paste(sprintf("%02d", d), month.abb[m], y, sep = sep),
    YYYY.MM.DD = paste(y, sprintf("%02d", m), sprintf("%02d", d), sep = sep),
    YYYY.MMM.DD = paste(y, month.abb[m], sprintf("%02d", d), sep = sep)
  )
}

#' Clean dates
#'
#' Clean a date and reformat it to another style.
#'
#' @param x a date string or a numeric representation of a date
#' (e.g. January 13th, 1991 would be 19910113)
#' @param original.format format of input `x`
#' @param preferred.format format to change `x` to
#' @param existing.missing.codes missing dates
#' @param return.missing.code what to return if there is a missing input
#' @param sep date separator. Defaults to "/"
#' @return A date string cleaned and formatted from the original (unformatted)
#' date
#' @family date formatting functions
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' cleanDate("09/11/1991", original.format = "MM.DD.YYYY", preferred.format = "DD.MM.YYYY")
#' cleanDate(11091991, original.format = "DD.MM.YYYY", preferred.format = "YYYY.MMM.DD")
#' cleanDate(11091991, original.format = "DD.MM.YYYY", preferred.format = "YYYY.MMM.DD", sep = "-")
cleanDate <- function(x, original.format, preferred.format,
                      existing.missing.codes = "", return.missing.code = NA,
                      sep = "/") {
  if (is.na(x))
    return(return.missing.code)
  x <- trimws(x)
  if (x %in% existing.missing.codes)
    return(x)
  date.comp <- strsplit(x, "/|-|\\|")[[1]]
  temp <- suppressWarnings(as.numeric(date.comp[1]))
  if (is.na(temp))
    return(return.missing.code)
  # yyyymmdd or ddmmyyyy or mmddyyyy
  if (temp > 1000000) {
    temp <- paste0(ifelse(temp < 10000000, "0", ""), temp)  # pad leading 0 for jan-sept; turn temp back to string
    if (original.format == "DD.MM.YYYY") {
      return(formatDate(substr(temp, 1, 2), substr(temp, 3, 4),
                        substr(temp, 5, 8), date.format = preferred.format,
                        sep = sep))
    } else if (original.format == "MM.DD.YYYY") {
      return(formatDate(substr(temp, 3, 4), substr(temp, 1, 2),
                        substr(temp, 5, 8), date.format = preferred.format,
                        sep = sep))
    } else if (original.format == "YYYY.MM.DD") {
      return(formatDate(substr(temp, 7, 8), substr(temp, 5, 6),
                        substr(temp, 1, 4), date.format = preferred.format,
                        sep = sep))
    } else {
      stop('ERROR (cleanDate): original.format must be one of
           "DD.MM.YYYY", "MM.DD.YYYY", or "YYYY.MM.DD".')
    }
    # must be YYYY/MM/DD
  } else if (temp > 31) {
    return(formatDate(date.comp[3], date.comp[2], date.comp[1],
                      date.format = preferred.format, sep = sep))
    # must be DD/MM/YYYY
  } else if (temp > 12) {
    return(formatDate(date.comp[1], date.comp[2], date.comp[3],
                      date.format = preferred.format, sep = sep))
    # first component <= 12 ... however, we are not sure if it refers to a day or month
  } else {
    temp <- as.numeric(date.comp[2]) # second component can either be a day or month
    # must be MM/DD/YYYY
    if (temp > 12) {
      return(formatDate(date.comp[2], date.comp[1], date.comp[3],
                        date.format = preferred.format, sep = sep))
      # BOTH first & second component <=12; assume original.format
    } else {
      if (original.format == "MM.DD.YYYY") {
        return(formatDate(date.comp[2], date.comp[1], date.comp[3],
                          date.format = preferred.format, sep = sep))
      } else if (original.format == "DD.MM.YYYY") {
        return(formatDate(date.comp[1], date.comp[2], date.comp[3],
                          date.format = preferred.format, sep = sep))
      } else {
        stop("ERROR (cleanDate): unknown original.format specified: ",
             original.format)
      }
    }
  }
}
