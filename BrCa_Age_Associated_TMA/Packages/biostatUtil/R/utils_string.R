#' Escape string for regular expression
#'
#' Escape `[`, `]`, `(`, and `)` for use in `grep`.
#' @param x a character vector
#' @return A character vector with opening and closing square brackets and
#'   parentheses escaped for use in `grep`.
#' @author Samuel Leung
#' @export
#' @examples
#' escapeForGrep("[index]")
#' escapeForGrep("(parentheses)")
escapeForGrep <- function(x){
  sub("\\[", "\\\\[",
      sub("\\]", "\\\\]",
          sub("\\(", "\\\\(",
              sub("\\)", "\\\\)", x)
          )
      )
  )
}

#' Find occurences of string within another string
#'
#' Returns a vector of the indices where a string occurs in another string
#'
#' If `b` is longer than `a`, `indexOf` returns `NA`, since
#' it is not possible for a longer string to occur in a shorter string.
#'
#' @param a string to be checked against
#' @param b string to check
#' @param ignore.case logical; if `TRUE`, case is ignored when performing
#'   the check
#' @return Indices where `b` occurs in `a`. Returns `NA` if there
#'   are no occurences.
#' @author Samuel Leung
#' @seealso [stringr::str_locate_all()]
#' @export
#' @examples
#' indexOf("derek", "e")
#' indexOf("Animals", "a")
#' indexOf("Animals", "A")
#' indexOf("Animals", "a", ignore.case = TRUE)
indexOf <- function(a, b, ignore.case = FALSE) {
  if (ignore.case) {
    a <- toupper(a)
    b <- toupper(b)
  }
  b.len <- nchar(b)
  a.len <- nchar(a)
  if (b.len > a.len)
    return(NA)
  a.arr <- strsplit(a, "")[[1]]
  indexes <- c()
  i <- 1
  while (i <= (a.len - b.len + 1)) {
    if (paste(a.arr[i:(i + b.len - 1)], collapse = "") == b) {
      indexes <- c(indexes, i)
      i <- i + b.len - 1
    }
    i <- i + 1
  }
  if (length(indexes) > 0)
    return(indexes)
  else
    return(NA)
}

#' Is the first letter uppercase?
#'
#' Returns a logical indicating whether the first letter of a character string
#' is uppercase.
#'
#' If the input is an empty string, the function returns `TRUE`.
#'
#' @param x character string
#' @return logical; if `TRUE`, the first letter of the input string is
#'   uppercase.
#' @author Samuel Leung, Derek Chiu
#' @export
#' @examples
#' isFirstLetterUpperCase("peanut butter")
#' isFirstLetterUpperCase("peanut Butter")
#' isFirstLetterUpperCase("Samuel butter")
#' isFirstLetterUpperCase("")
isFirstLetterUpperCase <- function(x) {
  return(sub("(.).*", "\\1", x) %in% c(LETTERS, ""))
}

#' Simple capitalization
#'
#' Capitalize the first letter of every word in a character string.
#'
#' To capitalize only the first word, use `first.only = TRUE`.
#'
#' @param x character string
#' @param first.only logical; if `TRUE`, only the first word will be
#'   capitalized
#' @return A character string with every word's first letter capitalized.
#' @author Samuel Leung
#' @references
#'   http://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
#'
#' @export
#' @examples
#' simpleCap("clear cell")
#' simpleCap("high grade serous carcinoma")
#' simpleCap("ovarian cancer", first.only = TRUE)
simpleCap <- function(x, first.only = FALSE) {
  if (!first.only)
    s <- strsplit(x, " ")[[1]]
  else
    s <- x
  return(paste0(toupper(substring(s, 1, 1)), substring(s, 2),
                collapse = " "))
}
