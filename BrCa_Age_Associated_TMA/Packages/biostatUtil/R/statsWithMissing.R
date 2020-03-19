#' Calculates statistics with specifying missing values
#'
#' Returns various statistics, with the ability to specify which values are
#' missing.
#'
#' The statistics supported include `min, max, mean, sum, prod`, as
#' indicated by each function's prefix. In addition, `ratioWithMissing`
#' calculates the ratio `x / y`.
#'
#' @param x,y input vector
#' @param missing.value missing values in `x`
#' @param return.missing.value character to return for missing values
#' @name statsWithMissing
#' @note `NAs` are ignored.
#' @author Samuel Leung, Derek Chiu
#' @examples
#' x <- c(10:1)
#' y <- c(1:10)
#' z <- c(10:1)
#' minWithMissing(z, c(1, 3))
#' maxWithMissing(z, c(9, 10))
#' meanWithMissing(z, c(3:7))
#' sumWithMissing(z, c(9, 10))
#' prodWithMissing(z, c(9, 10))
#' ratioWithMissing(x, y, c(1:2))
#'
#' ## All missing
#' minWithMissing(z, c(1:10))
#' minWithMissing(z, c(1:10), return.missing.value = "all missing")
#' maxWithMissing(z, c(1:10))
#' maxWithMissing(z, c(1:10), return.missing.value = "all missing")
#' meanWithMissing(z, c(1:10))
#' meanWithMissing(z, c(1:10), return.missing.value = "all missing")
#' prodWithMissing(z, c(1:10))
#' prodWithMissing(z, c(1:10), return.missing.value = "all missing")
NULL

#' @rdname statsWithMissing
#' @export
minWithMissing <- function(x, missing.value = -1, return.missing.value = -1) {
  x.missing <- x %in% missing.value
  if (sum(x.missing) == length(x) | n_missing(x) == length(x))
    return(return.missing.value)
  return(min(as.numeric(x[!x.missing]), na.rm = TRUE))
}

#' @rdname statsWithMissing
#' @export
maxWithMissing <- function(x, missing.value = -1, return.missing.value = -1) {
  x.missing <- x %in% missing.value
  if (sum(x.missing) == length(x) | n_missing(x) == length(x))
    return(return.missing.value)
  return(max(as.numeric(x[!x.missing]), na.rm = TRUE))
}

#' @rdname statsWithMissing
#' @export
meanWithMissing <- function(x, missing.value = -1, return.missing.value = -1) {
  x.missing <- x %in% missing.value
  if (sum(x.missing) == length(x) | n_missing(x) == length(x))
    return(return.missing.value)
  return(mean(as.numeric(x[!x.missing]), na.rm = TRUE))
}

#' @rdname statsWithMissing
#' @export
sumWithMissing <- function(x, missing.value = -1, return.missing.value = -1) {
  x.missing <- x %in% missing.value
  if (sum(x.missing) == length(x) | n_missing(x) == length(x))
    return(return.missing.value)
  return(sum(as.numeric(x[!x.missing]), na.rm = TRUE))
}

#' @rdname statsWithMissing
#' @export
prodWithMissing <- function(x, missing.value = -1, return.missing.value = -1) {
  x.missing <- x %in% missing.value
  if (sum(x.missing) == length(x) | n_missing(x) == length(x))
    return(return.missing.value)
  return(prod(as.numeric(x[!x.missing]), na.rm = TRUE))
}

#' @rdname statsWithMissing
#' @export
ratioWithMissing <- function(x, y, missing.value = -1, return.missing.value = -1) {
  x[x %in% missing.value] <- y[y %in% missing.value] <- NA
  ratio <- as.numeric(x) / as.numeric(y)
  ratio[is.na(ratio)] <- return.missing.value
  return(ratio)
}
