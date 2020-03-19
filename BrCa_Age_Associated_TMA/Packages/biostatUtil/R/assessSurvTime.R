#' Assess survival time
#'
#' Given range of survival times and censoring status, provides different
#' time-related summary statistics
#'
#' The observation time is defined as the median time in days between all
#' `T1` and `T2`. Censoring time is the median time in days between
#' all `T1` and `T2` for events only.
#'
#' @param T1 vector of start dates
#' @param T2 vector of end dates
#' @param status logical;
#' @return A list with elements
#' \item{Otime}{Observation time}
#' \item{Stime}{Censoring time}
#' \item{Etime}{Time to end of study}
#' \item{KFT}{Known Function Time}
#' \item{RevKM}{Reverse Kaplan-Meier Time}
#' @author Samuel Leung
#' @export
#' @examples
#' set.seed(3)
#' starts <- seq(as.Date("2000/1/1"), as.Date("2003/1/1"), by = "quarter")
#' ends <-  seq(as.Date("2003/1/1"), as.Date("2006/1/1"), by = "quarter")
#' statuses <- sample(0:1, 13, replace = TRUE)
#' assessSurvTime(starts, ends, statuses)
assessSurvTime <- function(T1, T2, status) {
  # in case there are any is.na(status)
  # T2 may be NA as well for rfs!!!
  non.missing.cases <- !is.na(status) & !is.na(T2)
  T1 <- T1[non.missing.cases]
  T2 <- T2[non.missing.cases]
  status <- status[non.missing.cases]
  Otime <- T2 - T1
  Stime <- T2[status] - T1[status]
  Etime <- max(T2) - T1
  SurvTime <- T2 - T1
  KFT <- SurvTime
  KFT[status] <- T2[status] - T1[status]
  rev.status <- rep(1, length(status))
  rev.status[status] <- 0
  Ftime <- survfit(Surv(as.numeric(SurvTime), rev.status) ~ 1)
  SumServ <- utils::read.table(textConnection(utils::capture.output(Ftime)),
                               skip = 2, header = TRUE)
  MedianTime <- list(
    Otime = as.numeric(round(stats::median(Otime, na.rm = TRUE) /
                               NUM.DAYS.IN.YEAR, 2)),
    Stime = as.numeric(round(stats::median(Stime, na.rm = TRUE) /
                               NUM.DAYS.IN.YEAR, 2)),
    Etime = as.numeric(round(stats::median(Etime, na.rm = TRUE) /
                               NUM.DAYS.IN.YEAR, 2)),
    KFT = as.numeric(round(stats::median(KFT, na.rm = TRUE) /
                             NUM.DAYS.IN.YEAR, 2)),
    RevKM = as.numeric(round(SumServ[, "median"] /
                               NUM.DAYS.IN.YEAR, 2)))
  return(MedianTime)
}
