#' Top variables in mass spectrometry analysis
#'
#' Show top genes/peptides defined by significance level and absolute fold
#' change.
#'
#' The input `x` is the result matrix returned by `ms_summarize`. We
#' want to filter `x` such that only interesting variables remain: i.e.
#' those that show statistical significance in an overall test and a
#' scientifically relevant effect size. Variables of interest are summarized
#' either on the gene level or peptide level.
#'
#' By default, the level of statistical significance is set to 5\% for the
#' Benjamini-Hochberg adjusted omnibus test p-value. The minimum absolute fold
#' change for determining scientific relevance is 1.25 by default. These default
#' values can be modified for different studies or projects, but offer a general
#' measure of validity for users to start with.
#'
#' @param x data frame object returned by `ms_summarize`
#' @param alpha significance level
#' @param fc minimum absolute fold change
#' @inheritParams ms_summarize
#'
#' @return A data frame showing the top variables. If `level = "Gene"`, the
#'   return value is a 4 column data frame, showing the Gene, Accession, BH
#'   adjusted omnibus p-value, and absolute fold change columns from `x`.
#'   If `level = "Peptide"`, the return value is the same except the Gene
#'   and Accession columns are replaced with the AGDSM column from `x`.
#' @note The fold change criterion only needs to be satisfied for *one*
#'   comparison if the experiment has 3 or more sample groups. For example,
#'   suppose we have 1 control group and treatments A and B. We filter variables
#'   on the fold change criterion where the absolute fold change is greater than
#'   `fc` for *either* A vs. control or B vs. control.
#' @family Mass Spectrometry functions
#' @author Derek Chiu
#' @export
ms_top <- function(x, level = c("Gene", "Peptide"), alpha = 0.05, fc = 1.25,
                   path = NULL) {
  # Determine variable level
  Var <- switch(match.arg(level),
                Gene = c("Gene", "Accession"),
                Peptide = "AGDSM")

  # Find Adjusted Omnibus Pvalue and Absolute Fold Change column positions
  AOP <- grep("(?=.*adj)(?=.*omnibus)(?=.*p-*val)", names(x),
              ignore.case = TRUE, perl = TRUE, value = TRUE)
  AFC <- grep("(?=.*abs)(?=.*effect)", names(x),
              ignore.case = TRUE, perl = TRUE, value = TRUE)

  # Filter for AOP under statistical significance, AFC over scientific relevance
  # The use of `any` pertains to the @note in the documentation
  topMat <- x[x[AOP] < alpha & apply((x[AFC] > fc), 1, any), c(Var, AOP, AFC)]

  if (!is.null(path))
    readr::write_csv(topMat, path)
  return(topMat)
}
