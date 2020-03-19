# magrittr placeholder
globalVariables(".")

#' Collapse unique levels of variable into one string
#' @noRd
collapse_var <- function(x, collapse = " || ") {
  return(paste(unique(x), collapse = collapse))
}

#' Print confidence interval wrapper
#' @noRd
printCI <- function(z) {
  paste0(z[1], " (", z[2], " - ", z[3], ")")
}

#' Count number of missing elements
#' @noRd
n_missing <- function(x, na.rm = FALSE) {
  return(sum(is.na(x), na.rm = na.rm))
}

#' Missing Value Formatting
#'
#' Takes a numeric vector and replaces all missing codes with NA and returns a
#' factor if the variable is categorical or a numeric variable if it's numeric.
#'
#' @param y a vector.
#' @param type whether the variable is `"cat"` (categorical) or
#'   `"cont"` (continuous). Defaults to `"cat"`.
#' @param codes vector of missing codes to replace with `NA`
#' @return A categorical or numerical vector with all missing formatted as
#'   `NA`.
#' @author Aline Talhouk, Derek Chiu
#' @export
#'
#' @examples
#' y <- c(1:10, "Unk", 12)
#' formatNA(y)
formatNA <- function(y, type = c("cat", "cont"), codes = c("", "Unk", "N/A")) {
  y[y %in% c(codes, NA)] <- NA
  res <- switch(match.arg(type), cat = factor(y), cont = as.numeric(y))
  return(res)
}

#' Generate a legend
#'
#' Given a ggplot object, generates a legend
#'
#' @param a.gplot ggplot object
#' @return ggplot object with legend
#' @author Aline Talhouk
#' @export
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(purrr::map_chr(tmp$grobs, ~ .x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' Get the p-value
#' @param x an object from [survival::survdiff()]
#' @return the Chi-squared p-value
#' @references Christos Hatzis
#'   (https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html)
#' @export
getPval <- function(x) {
  return(stats::pchisq(x$chisq, length(x$n) - 1, lower.tail = FALSE))
}

#' Standard error of the mean
#'
#' @param x input vector
#' @param missing.value values that are missing
#' @param return.missing.value the value to return where there are missing values
#' @return The standard error of the mean of `x`
#' @author Samuel Leung
#' @references http://en.wikipedia.org/wiki/Standard_error
#' @export
sem <- function(x, missing.value = NA, return.missing.value = NA) {
  x <- x[!is.na(x)]
  if (!is.na(missing.value))
    x <- x[!x %in% missing.value]
  return(ifelse(length(x) == 0, return.missing.value,
                sqrt(stats::var(x) / length(x))))
}

## CONSTANTS ##
# Dates
MM.DD.YYYY <- "%m/%d/%Y"
DD.MM.YYYY <- "%d/%m/%Y"
DD_MMM_YY  <- "%d-%b-%y"
YYYY_MM_DD <- "%Y-%m-%d"
YYYYMMDD   <- "%Y%m%d"
DDMMYYYY   <- "%d%m%Y"
MMDDYYYY   <- "%m%d%Y"
DATE.ORIGIN <- as.Date("1970-01-01")
NUM.DAYS.IN.YEAR <- 365.241 #365.25
NUM.DAYS.IN.MONTH <- 30.5

# Styles
COL.TH.STYLE <- "border-bottom: 1px solid grey; border-top: 4px double grey; text-align: center; padding-right:10px; padding-right:10px;"
ROW.TH.STYLE <- "text-align: center; padding-right:10px; padding-right:10px;"
TABLE.CAPTION.STYLE <- "display: table-caption; text-align: left;"

ROW.TD.STYLE.FOR.MULTI.COX <- "border-bottom: 1px solid grey; text-align: center; padding-right:10px; padding-right:10px;"
ROW.TD.STYLE.FOR.MULTI.COX.ALIGN.TOP <- "border-bottom: 1px solid grey; text-align: center; vertical-align: text-top; padding-right:10px; padding-right:10px;"

# Values
VALUE.CODING.INIT.TREATMENT.NO <- "No treatment"
VALUE.CODING.INIT.TREATMENT.CHEMO.ONLY <- "Chemo only"
VALUE.CODING.INIT.TREATMENT.RT.ONLY <- "Radiation only"
VALUE.CODING.INIT.TREATMENT.VAG.BRACHY.ONLY <- "Vag Brachy only"
VALUE.CODING.INIT.TREATMENT.BOTH <- "Both"

# Events
OS.EVENT  <- "os.event"
OS.CENSOR <- "os.censor"
DSS.EVENT  <- "dss.event"
DSS.CENSOR <- "dss.censor"
RFS.EVENT  <- "rfs.event"
RFS.CENSOR <- "rfs.censor"

# Missing codes

# missing value code for values that are explicitily indicated as missing from
# data source e.g. "X" in grade
MISSING.EXPLICIT <- "N/A"
# missing because values was not found (e.g. in data files) but the value must
# exist somewhere.
MISSING.UNK <- "Unk"
# data point not mentioned in data file.
MISSING...NOT.FOUND.IN.DATA.FILE <- ""
# missing value code for values that are explicitily indicated as missing from
# data source e.g. "X" in grade
MISSING.BIOMARKER.EXPLICIT <- MISSING.UNK
# data point not mentioned in data file.
MISSING.BIOMARKER...NOT.FOUND.IN.DATA.FILE <- ""
# combined missing codes
ALL.MISSING.CODES <- unique(c(
  MISSING.EXPLICIT,
  MISSING...NOT.FOUND.IN.DATA.FILE,
  MISSING.UNK,
  MISSING.BIOMARKER.EXPLICIT,
  MISSING.BIOMARKER...NOT.FOUND.IN.DATA.FILE
))

# Labels
FIRTH.THRESHOLD <- 0.8 # percent of censor cases to use Firth in Cox model
FIRTH.CAPTION <- "<sup>F</sup>" # text to indicate values are Firth corrected

BCSS.TITLE <- "Breast cancer specific survival"
BCSS.XLAB  <- "Total follow-up (years)"
BCSS.YLAB  <- "Cumulative breast cancer specific survival (BCSS)"

DSS.TITLE <- "Disease specific survival (DSS)"
DSS.XLAB  <- BCSS.XLAB
DSS.YLAB  <- DSS.TITLE

OS.TITLE <- "Overall survival"
OS.XLAB <- DSS.XLAB
OS.YLAB <- OS.TITLE

RFS.TITLE <- "Any relapse-free survival"
RFS.XLAB <- paste(RFS.TITLE, "time")
RFS.YLAB <- RFS.TITLE

DRFS.TITLE <- "Distant relapse-free survival"
DRFS.XLAB <- paste(DRFS.TITLE, "time")
DRFS.YLAB <- DRFS.TITLE

LRFS.TITLE <- "Rocal relapse-free survival"
LRFS.XLAB <- paste(LRFS.TITLE, "time")
LRFS.YLAB <- LRFS.TITLE

RRFS.TITLE <- "regional relapse-free survival"
RRFS.XLAB <- paste(RRFS.TITLE, "time")
RRFS.YLAB <- RRFS.TITLE

LRRFS.TITLE <- "Locoregional relapse-free survival"
LRRFS.XLAB <- paste(LRRFS.TITLE, "time")
LRRFS.YLAB <- LRRFS.TITLE
