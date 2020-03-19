#' Confusion matrix summaries
#'
#' Calculates summaries from cross-tabulated reference and prediction labels for
#' a multi-class variable.
#'
#' Given two multi-class variables summarized in a confusion matrix, this
#' function provides performance summaries. It provides overall accuracy with
#' confidence intervals, as well as per class accuracy, sensitivity,
#' specificity, positive predictive value (PPV), negative predictive value (NPV).
#' if variable entered is binary, it will automatically call binaryCM
#'
#' @inheritParams binaryCM
#' @return A confusion matrix for the predicted and reference classes. Then the
#'   estimated statistics along with bootstrapped confidence intervals. A
#' list with the following elements
#' \item{Accuracy}{Accuracy point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{Sensitivity}{Sensitivity point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{Specificity}{Specificity point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{PPV}{PPV point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{NPV}{NPV point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' \item{kappa}{kappa point estimate, lower bound and upper bound for
#' bootstrapped CI}
#' @family confusion matrix functions
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' ### 95% CI from 1000 bootstraped samples
#' set.seed(23)
#' k <- 3
#' (x <- factor(sample(1:k, 100, replace = TRUE, prob = c(0.15, 0.25, 0.6))))
#' (y <- factor(sample(1:k, 100, replace = TRUE, prob = c(0.05, 0.4, 0.65))))
#' prop.table(table(y))
#' multiClassCM(x, y)
#'
#' ### 90% CI from 500 bootstrapped samples
#' multiClassCM(x, y, num.boot = 500, conf.level = 0.90)
#'
#' ### Round to 2 digits
#' multiClassCM(x, y, digits = 2)
multiClassCM <- function(x, y, seed = 20, num.boot = 1000, conf.level = 0.95,
                         digits = 2, method = "wilson") {
  if (!identical(sort(unique(x)), sort(unique(y)))) {
    stop("levels should be the same in the reference and predicted classes")
  }
  CM <- table(Reference = x, Prediction = y)
  clm <- colSums(CM)
  rwm <- rowSums(CM)
  N <-  sum(CM)
  TP <- diag(CM)
  FP <- clm - TP
  FN <- rwm - TP
  TN <- N - (TP + FP + FN)
  sens <- TP / clm
  spec <- TN / (N - clm)
  BA <- (sens + spec) / 2
  # Overall
  cc <- round(caret::confusionMatrix(y, x)$overall, digits)
  ckappa <- round(kappaBootCI(x, y, seed, num.boot, conf.level), digits)
  overall <- rbind(printCI(cc[c("Accuracy", "AccuracyLower", "AccuracyUpper")]),
                   printCI(ckappa),
                   cc["AccuracyNull"],
                   cc["AccuracyPValue"]) %>%
    magrittr::set_rownames(c("Overall Accuracy", "Cohen's kappa",
                             "No Information Rate", "P-Value [Acc > NIR]")) %>%
    magrittr::set_colnames("Overall Concordance Statistics")
  # By class
  successes <- list(TP, TN, TP, TN, clm, TP, rwm, TP + TN)
  observations <- list(clm, N - clm, rwm, N - rwm, N, N, N, N)
  stats <- purrr::map2(successes, observations, Hmisc::binconf,
                       alpha = 1 - conf.level, method = method) %>%
    purrr::map(round, digits) %>%
    purrr::set_names(c("Sensitivity", "Specificity", "Pos Pred Value",
                       "Neg Pred Value", "Prevalence", "Detection Rate",
                       "Detection Prevalence", "Accuracy"))
  # Result table
  Average <- purrr::map(stats, ~ .x[, "PointEst"]) %>%
    purrr::map_dbl(~ round(mean(.x), digits = digits))
  ByClass <- purrr::map_df(stats, apply, 1, printCI) %>%
    t() %>%
    magrittr::set_colnames(colnames(CM))
  table <- cbind(Average, ByClass) %>%
    rbind(`Balanced Accuracy` = round(c(mean(BA), BA), digits))

  list(CM = stats::addmargins(CM), overall = overall, table = table)
}
