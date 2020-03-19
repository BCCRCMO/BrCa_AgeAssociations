#' Generate lifetables for multiclass variables
#'
#' Specify vector of time endpoints and create a cohort life table for two or
#' more strata
#'
#' Essentially a wrapper around [KMsurv::lifetab()] that allows the
#' user to input a `survfit` object instead of vectors of raw values.
#'
#' @param obj An object of class `survfit`
#' @param ntimes number of time intervals
#' @param times A vector of endpoints of time intervals to show life table
#'   calculations. By default, these are `ntimes` evenly spaced out
#'   endpoints based on the full range of survival times.
#' @param nround number of digits to round table values
#' @param show.strata logical; if `TRUE` (default), the variable name is
#'   appended to the beginning of each stratum in the lifetable's `strata`
#'   column
#' @param strata.name column name for the different strata
#' @param summary logical; if `TRUE`, a case processing summary is shown
#'   with number of subjects, events, censored, and percent censored per
#'   stratum.
#' @return A table with the following columns:
#' \item{strata}{name of specific group in variable}
#' \item{times}{time interval}
#' \item{nsubs}{See [KMsurv::lifetab()]}
#' \item{nlost}{See [KMsurv::lifetab()]}
#' \item{nrisk}{See [KMsurv::lifetab()]}
#' \item{nevent}{See [KMsurv::lifetab()]}
#' \item{surv}{See [KMsurv::lifetab()]}
#' \item{pdf}{See [KMsurv::lifetab()]}
#' \item{hazard}{See [KMsurv::lifetab()]}
#' \item{se.surv}{See [KMsurv::lifetab()]}
#' \item{se.pdf}{See [KMsurv::lifetab()]}
#' \item{se.hazard}{See [KMsurv::lifetab()]}
#' @author Derek Chiu
#' @seealso [KMsurv::lifetab()]
#' @export
#' @examples
#' library(survival)
#' obj <- survfit(Surv(futime, fustat) ~ rx, data = ovarian)
#' lifetable(obj)
#' lifetable(obj, ntimes = 4, show.strata = FALSE)
#' lifetable(obj, ntimes = 4, times = c(200, 500, 800, 1000))
lifetable <- function(obj, ntimes = 3, times = NULL, nround = 3,
                      show.strata = TRUE, strata.name = "strata",
                      summary = FALSE) {
  nevent <- nlost <- nsubs <- plost <- NULL
  cuts <- cumsum(obj$strata)
  times <-
    times %||% round(stats::quantile(obj$time, 1 / ntimes * seq_len(ntimes)))
  if (ntimes > 1) {
    ind <- purrr::map(split_pos(obj$time, cuts),
                      ~ purrr::map_int(times, function(.y)
                        which.min(abs(.y - .x)))) %>%
      purrr::map2(.y = obj$strata, ~ magrittr::inset(.x, ntimes, .y))
  } else {
    ind <- obj$strata
  }
  cs <- purrr::map2(split_pos(obj$n.censor, cuts), ind,
                    ~ purrr::map_dbl(split_pos(.x, .y), sum))
  es <- purrr::map2(split_pos(obj$n.event, cuts), ind,
                    ~ purrr::map_dbl(split_pos(.x, .y), sum))
  if (show.strata)
    strata <- names(obj$strata)
  else
    strata <- gsub(".+=", "\\1", names(obj$strata))
  tab <- purrr::pmap(list(obj$n, cs, es), KMsurv::lifetab,
                     tis = c(0, times)) %>%
    purrr::map(round, nround) %>%
    purrr::map(~ cbind(times = rownames(.x), .)) %>%
    purrr::invoke(rbind, .) %>%
    cbind(strata = rep(strata, each = ntimes), .) %>%
    magrittr::set_rownames(NULL)
  if (summary) {
    tab <- tab %>%
      select(strata, nsubs, nevent, nlost) %>%
      mutate(strata = as.character(strata)) %>%
      rbind(c("Overall", colSums(.[-1]))) %>%
      mutate_at(vars(names(.)[-1]), as.numeric) %>%
      mutate(plost = paste0(sprintf("%.1f", nlost / nsubs * 100), "%")) %>%
      dplyr::rename(`Total N` = nsubs, `N of Events` = nevent,
                    `N of Censored` = nlost, `Percent Censored` = plost)
  }
  colnames(tab)[1] <- eval(strata.name)
  return(tab)
}

#' Split vector at a position for calculating cumulative sums
#' @noRd
split_pos <- function(x, pos) {
  unname(split(x, cumsum(seq_along(x) %in% (pos + 1))))
}
