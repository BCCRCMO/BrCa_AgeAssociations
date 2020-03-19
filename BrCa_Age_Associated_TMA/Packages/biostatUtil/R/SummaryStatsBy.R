#' Generate summary statistics
#'
#' Splits the data by two variables, computing relevant statistics for each
#' variable.
#'
#' The `data` is split by two variables, `by1` and `by2`, and
#' statistics are computed for continuous variables. Statistics currently
#' supported include `mean, sd, median, IQR, range`, and the number of
#' missing cases. For factor variables, the counts, column and row percentages
#' are shown for each of the variable levels.
#'
#' Note that marginal statistics are also shown for `by1`, so the order in
#' which you split `data` matters.
#'
#' There are four print options for the output: `raw` gives the output as a
#' character matrix, `pandoc` gives a Pandoc-friendly output for Word and
#' PDF reports, `html` gives HTML supported output, and `long` is a
#' tidy version of `raw`.
#'
#' @param data The `data.frame` containing the data
#' @param by1 character string of splitting variable 1
#' @param by2 character string of splitting variable 2
#' @param var.names character vector of variables to compute statistics for
#' @param stats statistics to compute for continuous variables
#' @param digits number of digits to round to
#' @param format format to return the table in. Either "raw", "pandoc" (for Word
#'   and PDF), "html", or "long" format for graphing and data manipulation using
#'   raw values.
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' mtcars$vs <- as.factor(mtcars$vs); mtcars$am <- as.factor(mtcars$am)
#' SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = c("mpg", "vs",
#' "qsec", "am"))
#' SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = c("vs",
#' "qsec"))
#' SummaryStatsBy(mtcars, by1 = "cyl", by2 = "gear", var.names = c("mpg"))
SummaryStatsBy <- function(data, by1, by2, var.names,
                           stats = c("mean", "sd", "median", "IQR", "range",
                                     "missing"), digits = 3,
                           format = c("raw", "pandoc", "html", "long")) {
  assertthat::assert_that(n_distinct(data[, by1]) >= 2,
                          n_distinct(data[, by2]) >= 2)
  bys <- c(by1, by2)
  stats <- match.arg(stats, c("mean", "sd", "median", "IQR", "range",
                              "missing"), several.ok = TRUE)
  types <- data[, var.names, drop = FALSE] %>% purrr::map_chr(class)
  num.ind <- types %in% c("numeric", "integer")
  fac.ind <- types %in% c("factor", "character")
  if (!(num.ind || fac.ind)) {
    stop("Variables must be numeric, integer, factor, or character.")
  }
  num.long <- num.res <- fac.long <- fac.res <- NULL
  col.names <- data[, bys] %>%
    unique() %>%
    rbind(matrix(c(unique(data[, by1]), rep(NA, n_distinct(data[, by1]))),
                 ncol = 2, dimnames = list(NULL, bys))) %>%
    arrange_(by1, by2) %>%
    purrr::map2(names(.), ., paste, sep = "=") %>%
    purrr::pmap_chr(paste, sep = ", ") %>%
    ifelse(grepl("NA", .), gsub(",.*", "\\1", .), .)

  # Compute numerical summaries
  if (sum(num.ind) > 0) {
    num.var <- var.names[num.ind]
    num.dat <- data[, c(num.var, bys)]
    num.ord <- purrr::map(num.var, ~ paste0(.x, c("", paste0(".", stats)))) %>%
      unlist()
    if (all(c("mean", "sd") %in% stats)) {
      num.ord <- num.ord[-grep("sd", num.ord)]
    }
    num.val <- paste_formula(num.var, bys) %>%
      doBy::summaryBy(num.dat, FUN = contSumFunc, digits = digits,
                      stats = stats) %>%
      mutate_all(funs(as.character))
    num.val.tot <- paste_formula(num.var, by1) %>%
      doBy::summaryBy(num.dat, FUN = contSumFunc, digits = digits,
                      stats = stats) %>%
      mutate_all(funs(as.character))
    num.all <- list(num.val, num.val.tot)
    num.long <- num.all %>%
      dplyr::bind_rows() %>%
      wtl(by1, by2)
    nnum <- n_distinct(num.long[, bys])
    num.res <- matrix(num.long$value, ncol = nnum,
                      dimnames = list(grep(collapse_var(num.var, "|"),
                                           names(num.val), value = TRUE),
                                      col.names)) %>%
      rbind(matrix("", nrow = length(num.var), ncol = nnum,
                   dimnames = list(num.var))) %>%
      magrittr::extract(num.ord, ) %>%
      magrittr::set_rownames(num.ord %>%
                               ifelse(!grepl("\\.", .),
                                      paste0("**", ., "**"), .) %>%
                               gsub(".+\\.", "\\1", .))
  }

  # Compute factor summaries
  if (sum(fac.ind) > 0) {
    fac.var <- var.names[fac.ind]
    fac.dat <- data[, c(fac.var, bys)]
    fac.ord <- purrr::map2(fac.var, fac.dat[, fac.var, drop = FALSE],
                           ~ paste0(.x, c("", paste0(".", levels(.y))))) %>%
      unlist()
    fac.val <- fac.var %>%
      purrr::map(~ paste_formula(.x, bys)) %>%
      purrr::map(~ as.matrix(aggregate(.x, fac.dat, summary))) %>%
      Reduce(merge, .) %>%
      as.data.frame()
    fac.val.tot <- fac.var %>%
      purrr::map(~ paste_formula(.x, by1)) %>%
      purrr::map(~ as.matrix(aggregate(.x, fac.dat, summary))) %>%
      Reduce(merge, .) %>%
      as.data.frame()
    fac.all <- dplyr::bind_rows(list(fac.val, fac.val.tot))
    fac.long <- wtl(fac.all, by1, by2)
    fac.pct <- fac.val %>%
      select(-one_of(bys)) %>%
      t() %>%
      as.data.frame() %>%
      split(factor(gsub("\\..*", "\\1", rownames(.)), levels = fac.var)) %>%
      purrr::map(rowColPercent) %>%
      purrr::map_df(split_pcts, n = 3)
    fac.pct.tot <- fac.val.tot %>%
      split(.[, by1]) %>%
      purrr::map(~ .x %>%
                   select(-matches(by1)) %>%
                   t() %>%
                   as.data.frame() %>%
                   split(rep(seq_len(nrow(.) / 2), each = 2)) %>%
                   purrr::map(colPercent) %>%
                   do.call(rbind, .)) %>%
      do.call(rbind, .) %>%
      split(factor(gsub("\\..*", "\\1", rownames(.)), levels = fac.var)) %>%
      purrr::map_df(~ split_pcts(matrix(.x, nrow = 4), n = 2))
    fac.res <- cbind(fac.pct, fac.pct.tot) %>%
      as.matrix() %>%
      magrittr::set_rownames(utils::tail(names(fac.val), -2)) %>%
      rbind(matrix("", nrow = length(fac.var), ncol = nrow(fac.all),
                   dimnames = list(fac.var, NULL))) %>%
      magrittr::extract(fac.ord, order(fac.all[, by1])) %>%
      magrittr::set_rownames(stringr::str_replace_all(
        rownames(.),
        c(stats::setNames(c(rep("", length(fac.var))), paste0(fac.var, ".")),
          stats::setNames(paste0("**", fac.var, "**"), fac.var)))) %>%
      magrittr::set_colnames(col.names)
  }

  # Final results in each format
  final.res <- rbind(num.res, fac.res)
  ind <- grep("\\*", rownames(final.res))
  org.ord <- gsub("\\*\\*", "", rownames(final.res)[ind])
  final.reord <- final.res %>%
    magrittr::extract(order(unlist(
      purrr::map2(match(org.ord, var.names),
                  diff(c(ind, nrow(final.res) + 1)),
                  rep))), )
  final.html <- final.res %>%
    magrittr::set_rownames(stringr::str_replace_all(
      rownames(.), c("^\\*\\*" = "<b>", "\\*\\*$" = "</b>")))
  final.long <- rbind(num.long, fac.long) %>%
    magrittr::extract(order(match(.$var, var.names)), )
  switch(match.arg(format),
         raw = final.reord,
         pandoc = pander::pandoc.table.return(final.reord,
                                              emphasize.rownames = FALSE),
         html = htmlTable::htmlTable(final.html),
         long = final.long)
}

#' Split data by row according to variables before calculating percentages
#' @noRd
split_pcts <- function(x, n) {
  x %>%
    as_data_frame() %>%
    split(rep(seq_len(nrow(.) / n), each = n)) %>%
    purrr::map_df(., . %>% purrr::map_df(pandoc_pcts))
}

#' Formatting percentages for Pandoc
#' @noRd
pandoc_pcts <- function(char) {
  count <- as.integer(char[1])
  pcts <- as.numeric(char[-1]) * 100
  if (length(char) > 2)
    paste0(count, " (", pcts[1], "%, ", pcts[2], "%)")
  else
    paste0(count, " (", pcts, "%)")
}

#' Continuous summary functions
#' @noRd
contSumFunc <- function(x, digits, stats = c("mean", "sd", "median", "IQR",
                                             "range", "missing")) {
  stats.choices <- c("mean", "sd", "median", "IQR", "range", "missing")
  funs.arg <- match.arg(stats, stats.choices, several.ok = TRUE)
  if ("missing" %in% stats) funs.arg[match("missing", funs.arg)] <- "n_missing"
  all.stats <- purrr::map_chr(funs.arg, ~ {
    match_fun_null(x = x, .x, na.rm = TRUE) %>%
      round(., digits = digits) %>%
      as.character() %>%
      ifelse(length(.) > 1, collapse_var(., "-"), .)
  }) %>%
    magrittr::set_names(stats)
  if (all(c("mean", "sd") %in% stats)) {
    all.stats["mean"] <- collapse_var(all.stats[c("mean", "sd")], " &#177; ")
    all.stats <- all.stats[-match("sd", names(all.stats))]
  }
  all.stats
}

#' Apply function on every element of list
#' @noRd
match_fun_null <- function(x, FUN, ...) {
  do.call(FUN, c(list(x), ...))
}

#' Construct formula object from character strings of response and terms
#' @noRd
paste_formula <- function(response, terms) {
  stats::as.formula(paste(collapse_var(response, " + "), "~",
                          collapse_var(terms, " + ")))
}

#' Munge data from Wide To Long format
#' @noRd
wtl <- function(data, by1, by2) {
  data %>%
    tidyr::gather_("stat", "value",
                   grep(collapse_var(c(by1, by2), "|"), names(.),
                        invert = TRUE, value = TRUE)) %>%
    tidyr::separate_("stat", c("var", "stat"), "\\.") %>%
    arrange_(by1, by2)
}
