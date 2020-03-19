#' Graphical Exploration of a Dataset
#'
#' Function to provide a graphical exploration of a dataset will print results to file.
#' @param datmat the data frame (only categorical and numerical variables will be analysed)
#' @author Aline Talhouk
#' @export
#' @examples
#' mtcars$vs <- as.factor(mtcars$vs)
#' mtcars$am <- as.factor(mtcars$am)
#' exploreData(mtcars)
#' file.remove("DataSummary.pdf")
exploreData <- function(datmat) {
  types <- unname(purrr::map_chr(datmat, class))
  fd <- datmat[, types %in% c("factor", "numeric", "integer")]
  type.fd <- unname(purrr::map_chr(fd, class))
  num.ind <- type.fd %in% c("numeric", "integer")
  fac.ind <- type.fd %in% c("factor")
  catvars <- colnames(fd)[fac.ind]
  grDevices::pdf("DataSummary.pdf")

  for (i in seq_along(catvars)) {
    x <-  fd[, catvars[i]]
    tx <-  table(x, useNA = "ifany")
    graphics::par(mfrow = c(2, 1), mar = c(3.1, 9.5, 4.1, 2.1))
    barplotSum(tx, catvars[i])
    mat <- data.matrix(cbind(tx, round(prop.table(tx) * 100, 1)))
    colnames(mat) <- c("Freq", "%")
    rownames(mat)[is.na(rownames(mat))] <- "Missing"
    tmat <- rbind(mat, apply(mat[, 1:2], 2, sum))
    PerformanceAnalytics::textplot(tmat, wrap = FALSE)
  }

  numvars <- colnames(fd)[num.ind]
  for (i in seq_along(numvars)) {
    graphics::par(mfrow = c(2, 1))
    x <- fd[, numvars[i]]
    boxplotSum(x, numvars[i])
    histSum(x)
  }
  grDevices::dev.off()
}
