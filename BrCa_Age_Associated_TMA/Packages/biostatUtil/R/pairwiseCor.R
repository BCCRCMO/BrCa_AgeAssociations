#' Pairwise Correlations
#'
#' Computes all pairwise correlations between the columns of a data frame
#'
#' @param x A data frame containing numeric variables of interest.
#' @return all pairwise absolute correlations, correlations, Pval, Adj P val by
#'   decreasing order of absolute correlations.
#'
#' @author Aline Talhouk, Derek Chiu
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- data.frame(matrix(rnorm(25), nrow = 5))
#' pairwiseCor(x)
pairwiseCor <- function(x) {
  if (!all(purrr::map_lgl(x, is.numeric)))
    stop("All columns of data matrix must be numeric")
  Cor <- AbsCor <- Pval <- AdjP <- NULL
  pairs <- utils::combn(names(x), 2) %>%
    magrittr::set_rownames(paste0("Variable", 1:2))
  pairwiseCorDF <- data.frame(Cor = apply(pairs, 2, function(df)
    stats::cor(x[, df]))[2, ]) %>%
    mutate(AbsCor = abs(Cor),
           Pval = purrr::map2_dbl(x[pairs[1, ]], x[pairs[2, ]],
                                  ~ stats::cor.test(.x, .y)$p.value),
           AdjP = stats::p.adjust(Pval, "fdr")) %>%
    select(AbsCor, Cor, Pval, AdjP) %>%
    round(4) %>%
    data.frame(t(pairs), ., stringsAsFactors = FALSE) %>%
    arrange(desc(AbsCor))
  return(pairwiseCorDF)
}
