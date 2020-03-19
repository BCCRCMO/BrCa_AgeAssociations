#' Generate a Frequency Table
#'
#' A frequency table emulating the SPSS FREQ output is generated.
#'
#' @param x vector of values to show frequencies for
#' @param levels (optional) vector for order of levels in `x`
#' @param missing vector of levels for which we consider missing and don't count
#'   in the valid percentage
#' @param description (optional) description for each level of `x`. Must be
#'   same length and order as `levels`
#' @param round number of digits to round percentages to
#' @param plot logical; if `TRUE`, a barplot is produced.
#'
#' @return A data frame with the following columns
#' \item{Class}{Tells you which scores are valid and which are missing. Not shown if
#'   there are no missing values.}
#' \item{Score}{Different levels}
#' \item{Frequency}{Count for each score}
#' \item{Percent}{Percent of Frequency out of the grand total}
#' \item{Valid Percent}{Percent of Frequency out of the Valid scores. Not applicable if
#'   there are no missing values.}
#' \item{Cumulative Percent}{Accumulated Percent of Frequency out of the Valid Scores}
#' \item{Description}{If `description` is given, a description for each level}
#'
#' @author Derek Chiu
#' @seealso [freq()]
#' @export
#'
#' @examples
#' # Create vector of randomly reordered alphabet with various frequencies
#' # for each letter
#' set.seed(123)
#' n <- sample(10, length(letters), replace = TRUE)
#' x <- sample(rep.int(letters, times = n))
#' freqTable(x, plot = TRUE)
#'
#' # Treat vowels as missing
#' freqTable(x, missing = c("a", "e", "i", "o", "u"), round = 2)
freqTable <- function(x, levels = sort(unique(as.character(x))),
                      missing = NULL, description = NULL,
                      round = 1, plot = FALSE) {
  . <- Class <- Frequency <- Score <- `Valid Percent` <- NULL

  tab <- descr::freq(x, user.missing = missing, plot = plot) %>%
    as.data.frame() %>%
    cbind(Score = factor(rownames(.), c(levels, "Total")), .) %>%
    cbind(Class = factor(ifelse(is.na(.$`Valid Percent`), "Missing",
                                ifelse(grepl("Total", .$Score), "Total", "Valid")),
                         c("Valid", "Missing", "Total")), .) %>%
    rbind(c("Valid", "Total",
            sum(.$Frequency[!is.na(.$`Valid Percent`) & .$Class == "Valid"]),
            sum(.$Percent[!is.na(.$`Valid Percent`) & .$Class == "Valid"]), 100)) %>%
    arrange(Class, Score) %>%
    mutate_each(funs(as.numeric), matches("Frequency|Percent")) %>%
    mutate(Score = ifelse(Class == "Total" & Score == "Total", "",
                          ifelse(is.na(Score), "Total", as.character(Score))),
           `Valid Percent` = ifelse(Class == "Total", NA, `Valid Percent`),
           `Cumulative Percent` = ifelse(!is.na(`Valid Percent`) & Score != "Total",
                                         cumsum(Frequency) / max(Frequency[Class == "Valid"]) * 100, NA),
           Class = ifelse(duplicated(Class), "", as.character(Class))) %>%
    mutate_each(funs(round(., 1)), contains("Percent"))
  if (is.null(missing)) {
    tab <- tab %>%
      magrittr::extract(-which(.$Class == "Total"), ) %>%
      select(-matches("Class|Valid Percent"))
  }
  if (!is.null(description)) {
    assertthat::assert_that(length(levels) == length(description))
    desc.ord <- order(match(levels, with(tab, Score[!grepl("Total|^$", Score)])))
    desc <- append(description[desc.ord], "", which(tab$Score == "Total") - 1)
    if (!is.null(missing))
      desc <- c(desc, "")
    tab <- mutate(tab, Description = desc)
  }
  return(tab)
}
