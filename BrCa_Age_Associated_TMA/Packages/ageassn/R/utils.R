# magrittr placeholder
globalVariables(".")

# global environment
pos <- 1

#' Turn colour character string to hexadecimal form with alpha transparency
#' @param col colour
#' @param alpha alpha transparency value. Ranges from 0 to 255.
#' @noRd
col2rgba <- function(col, alpha = 30) {
  rgb(t(col2rgb(col)), alpha = alpha, maxColorValue = 255)
}

#' value at index if TRUE, otherwise return real NA
#' @param value numeric vector of values
#' @param index logical vector same length as `value` used for indexing
#' @noRd
value_if <- function(value, index) {
  purrr::map2_dbl(value, index, ~ ifelse(.y, .x[.y], NA_real_))
}

#' Tidy gene names
#'
#' Extract gene names from ER binding gene lists.
#'
#' Gene names with the empty string `""` or `NA` are removed. Duplicate gene
#' names are removed and finally the gene names are sorted.
#'
#' @param data ER binding gene dataset
#' @param var variable name with gene symbols in `data`
#' @export
tidy_genes <- function(data, var) {
  data[[var]] %>%
    `[`(!(is.na(.) | . == "")) %>%
    unique() %>%
    sort()
}

#' Tidy contingency tables
#'
#' Cross tabulation with Chi-Squared Test and Fisher Exact Test results
#'
#' The Fisher Exact Test is recalculated when there are errors pertaining to an
#' insufficient workspace ("LDKEY is too small for this problem.")
#'
#' @param table contingency table. Cross tabulation is performed by
#'   [descr::CrossTable()].
#' @param workspace fisher exact test workspace. Increase as needed.
#' @export
tidy_ct <- function(table, workspace = 1e6) {
  ct <- descr::CrossTable(
    table,
    prop.r = TRUE,
    prop.c = TRUE,
    prop.t = FALSE,
    prop.chisq = FALSE,
    chisq = TRUE,
    fisher = TRUE,
    format = "SPSS"
  )
  if (all(is.na(ct$fisher.ts))) {
    tab <- ct$tab
    ct$fisher.ts <- stats::fisher.test(tab, workspace = workspace)
  }
  ct
}



#' Tidy Fisher Exact Test results
#'
#' Extract relevant FET results from object of class "htest" into a table.
#'
#' The estimated odds ratio statistic and corresponding confidence interval only
#' apply to the 2 by 2 case. The FET p-value is returned for contingency tables
#' of any dimension.
#'
#' @param object object of class "htest" with FET results
#' @param digits number of digits to round statistics
#' @export
tidy_fet <- function(object, digits = 2) {
  fet_table <- object %>%
    tibble::enframe(name = "Statistics", value = "Value") %>%
    dplyr::mutate(
      Value = purrr::map_if(.data$Value, is.numeric, signif, digits = digits)
    ) %>%
    head(-2)
  if ("estimate" %in% names(object)) {
    dplyr::mutate(
      fet_table,
      Value = purrr::map_at(.data$Value, 2,
                            ~ paste0("[", paste(., collapse = ", "), "]"))
    )
  } else {
    fet_table
  }
}

#' Convert strings for use in file names
#'
#' Certain characters are not advisable for use in file names. `str_decimal`
#' converts decimal numbers and `str_sign` converts positive and negative signs.
#'
#' Decimal numbers contain a period character ".", which is traditionally used
#' to separate the name and extension of a file. `str_decimal` replaces all
#' periods with "p" (for "period").
#'
#' The plus sign "+", minus sign "-", and division sign "/" are all not
#' recommended for use in file names. There are two syntaxes used for
#' replacement. In the short version, they are replaced by "p", "n", and "_",
#' respectively. In the long version, they are replaced by "pos" (for
#' "positive"), "neg" (for "negative"), and "", respectively.
#' @param string character string
#' @name str_file
#' @export
#' @examples
#' str_decimal(1.23)
#' str_sign("ER+/HER2-")
#' str_sign("ER+/HER2-", type = "long")
str_decimal <- function(string) {
  gsub("\\.", "p", as.character(string))
}

#' @inheritParams str_decimal
#' @param type use "short" or "long" type of conversion
#' @name str_file
#' @export
str_sign <- function(string, type = c("short", "long")) {
  type <- match.arg(type)
  pattern <- switch(
    type,
    short = c("\\+" = "p", "\\-" = "n", "\\/" = "_"),
    long = c("\\+" = "pos", "\\-" = "neg", "\\/" = "")
  )
  stringr::str_replace_all(string, pattern)
}

#' Turn a variable name into a symbol for tidy evaluation specification in
#' tcga_aroutdf
#' @noRd
sym_var <- function(age, direction, name) {
  rlang::sym(paste0(direction, age, "_", name))
}
