#' Tests for Independence in Contingency Tables
#'
#' The Pearson's Chi-Squared test, likelihood ratio (G test) of independence,
#' Fisher's Exact test, and linear-by-linear association test are performed on
#' the data matrix.
#'
#' A Pearson's Chi-Squared test Yate's Continuity Correction is applied in the
#' case of 2 by 2 tables.
#'
#' @param x an object of class `CrossTable` containing the contingency
#'   table
#' @param digits number of digits to round to
#'
#' @return A table with method name, test statistic, degrees of freedom, and
#'   p-value reported for each Chi-squared test.
#' @author Derek Chiu
#' @seealso [descr::CrossTable()]
#' @export
#'
#' @examples
#' # Example from documentation of CrossTable
#' library(descr)
#' data(esoph, package = "datasets")
#' ct <- CrossTable(esoph$alcgp, esoph$agegp, expected = TRUE,
#'                  chisq = FALSE, prop.chisq = FALSE,
#'                  dnn = c("Alcohol consumption", "Tobacco consumption"))
#' indepTests(ct)
#'
#' # Better example
#' set.seed(1108)
#' A <- rbinom(100, 3, 0.2)
#' B <- rbinom(100, 4, 0.8)
#' ct <- CrossTable(A, B)
#' indepTests(ct)
indepTests <- function(x, digits = 3) {
  . <- `P-Value` <- Test <- Value <- df <- NULL
  Pearson <- x$CST
  if (any(is.na(Pearson))) {
    Pearson.obj <- rep(NA, 3)
  } else {
    if (any(Pearson$expected < 1) | mean(Pearson$expected < 5) > 0.2) {
      Pearson.obj <- rep(NA, 3)
    } else {
      Pearson.obj <- c(Pearson$statistic, Pearson$parameter, Pearson$p.value)
    }
  }

  CC <- x$chisq.corr
  if (all(is.na(CC))) {
    CC.obj <- rep(NA, 3)
  } else {
    if (!all(is.na(CC)) & !all(is.na(Pearson.obj))) {
      CC.obj <- c(CC$statistic, CC$parameter, CC$p.value)
    } else {
      CC.obj <- rep(NA, 3)
    }
  }

  G.test <- tryCatch(Deducer::likelihood.test(x$tab),
                     error = function(e) return(NULL))
  if (!is.null(G.test)) {
    G.test.obj <- c(G.test$statistic, G.test$parameter, G.test$p.value)
  } else {
    G.test.obj <- rep(NA, 3)
  }

  Fisher <- x$fisher.ts
  if (all(is.na(Fisher))) {
    Fisher.obj <- rep(NA, 3)
  } else {
    Fisher.obj <- c(NA, NA, Fisher$p.value)
  }

  LBL <- tryCatch(coin::lbl_test(x$tab),
                  error = function(e) return(NULL))
  if (!is.null(LBL)) {
    LBL.obj <- c(coin::statistic(LBL), 1, coin::pvalue(LBL))
  } else {
    LBL.obj <- rep(NA, 3)
  }

  res <- data.frame(Pearson.obj, CC.obj, G.test.obj, Fisher.obj, LBL.obj) %>%
    t() %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("Value", "df", "P-Value")) %>%
    magrittr::set_rownames(c("Pearson Chi-Square",
                             "Continuity Correction",
                             "Likelihood Ratio",
                             "Fisher's Exact Test",
                             "Linear-by-Linear Association")) %>%
    mutate(Test = rownames(.)) %>%
    mutate_each(funs(round(., digits)), 1:2) %>%
    mutate(`P-Value` = round_small(`P-Value`, digits = digits)) %>%
    select(Test, Value, df, `P-Value`) %>%
    rbind(., c("N of Valid Cases", x$gt, "", ""))
  return(res)
}
