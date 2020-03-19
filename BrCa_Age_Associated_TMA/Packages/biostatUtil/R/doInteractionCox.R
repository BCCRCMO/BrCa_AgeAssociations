#' Do interaction test with cox model
#' likelihood ratio test only - two terms only
#'
#' Only a LRT can be performed (2015-08-12). TODO: add penalized LRT?
#'
#' @inheritParams doCoxphGeneric
#' @return A list with the following elements
#' @note the order of variable names in var.names dictates when will the term
#'   added in stepwise likelihood ratio test of nested models i.e. given
#'   `var.names = c("A", "B", "C")`, likelihood ratio test of nested model will
#'   be:
#'
#' A
#' A + B
#' A + B + C
#' @author Samuel Leung
#' @export
doInteractionCox <- function(
    input.d, var.names, var.descriptions,
    show.var.detail = FALSE,
    show.group.name.for.bin.var = FALSE, var.ref.groups = NULL,
    var.names.surv.time = c("os.yrs", "dss.yrs", "rfs.yrs"),
    var.names.surv.status = c("os.sts", "dss.sts", "rfs.sts"),
    event.codes.surv = c("os.event", "dss.event", "rfs.event"),
    surv.descriptions = c("OS", "DSS", "PFS"),
    missing.codes = c("N/A", "", "Unk"),
    use.firth = 1, firth.caption = FIRTH.CAPTION,
    round.digits.p.value = 4, round.small = FALSE, scientific = FALSE,
    caption = NA, html.table.border = 0, banded.rows = FALSE,
    css.class.name.odd = "odd", css.class.name.even = "even",
    split.table = 300, ...) {

  # Constants
  kLocalConstantHrSepFlag <- "kLocalConstantHrSepFlag" # separates HR estimates
  col.th.style <- COL.TH.STYLE
  row.th.style <- ROW.TH.STYLE
  row.td.style.for.multi.cox <- ROW.TD.STYLE.FOR.MULTI.COX
  row.td.style.for.multi.cox.align.top <- ROW.TD.STYLE.FOR.MULTI.COX.ALIGN.TOP

  # Initial assertion checks
  num.surv.endpoints <- length(var.names.surv.time)
  assertthat::assert_that(num.surv.endpoints == length(var.names.surv.status),
                          num.surv.endpoints == length(event.codes.surv),
                          num.surv.endpoints == length(surv.descriptions))

  # Check interactions
  for (var.name in var.names) {
    if (!grepl(":", var.name)) {
      if (is.factor(input.d[, var.name]))
        input.d[, var.name] <- droplevels(input.d[, var.name])
    } else {
      # this must be an interaction term ... separate them and droplevels on
      # individual variables
      me.vars <- vapply(strsplit(var.name, ":")[[1]], stringr::str_trim,
                        FUN.VALUE = character(length(var.names)),
                        USE.NAMES = FALSE)
      for (me.var in me.vars) {
        input.d[, me.var] <- droplevels(input.d[, me.var])
        if (!me.var %in% var.names) {
          warning("main effect not in model: ", me.var)
        }
      }
    }
  }
  # Ensure survival times are numeric
  input.d <- input.d %>%
    dplyr::mutate_at(var.names.surv.time, as.numeric)

  # Setup default for variable reference groups and result matrix
  nvar <- length(var.names)
  var.ref.groups <- var.ref.groups %||% rep(NA, nvar)
  rn <- paste(var.names, rep(surv.descriptions, each = nvar), sep = "-")
  cn <- c("# of events / n", "Hazard Ratio (95% CI)", "LRT P-value")
  result.table <- matrix(NA_character_, nrow = nvar * num.surv.endpoints,
                         ncol = 3, dimnames = list(rn, cn))
  for (i in seq_along(var.names)) {
    x <- var.names[i]
    if (!grepl(":", x)) {
      input.d <- input.d %>%  # remove any cases with NA's or missing values
        dplyr::filter(!is.na(.[, x]) & !(.[, x] %in% missing.codes))
      # automatically set ref.group to lowest group if not specified
      if (is.factor(input.d[, x])) {
        input.d[, x] <- droplevels(input.d[, x])
        if (is.na(var.ref.groups[i]))
          var.ref.groups[i] <- names(table(input.d[, x]))[1]
      }
      if (is.na(var.ref.groups[i])) {
        input.d[, x] <- as.numeric(input.d[, x])
      } else {
        var.levels <- names(table(input.d[, x]))
        var.levels <- c(var.ref.groups[i], var.levels[-which(
          var.levels == var.ref.groups[i])])
        input.d[, x] <- factor(input.d[, x], levels = var.levels)
      }
    }
  }

  cox.stats.output <- list()
  for (j in seq_len(num.surv.endpoints)) {
    temp.d <- input.d %>%
      dplyr::filter(!is.na(.[, var.names.surv.status[[j]]]) &
                      !is.na(.[, var.names.surv.time[[j]]]))
    surv.formula <- paste0("Surv(", var.names.surv.time[j], ", ", var.names.surv.status[j], "=='", event.codes.surv[j], "'  ) ~")
    curr.var.names <- c()
    var.idx <- 0
    for (i in seq_along(var.names)) {
      curr.var.names <- c(curr.var.names, var.names[i])
      var.idx <- max(var.idx) + 1
      if (!is.na(var.ref.groups[i]))
        var.idx <- c(var.idx:(var.idx + dplyr::n_distinct(temp.d[, var.names[i]]) - 2))
      cox.stats <- prettyCoxph(stats::as.formula(paste(surv.formula, paste(curr.var.names, collapse = "+"))),
                               input.d = temp.d, use.firth = use.firth)
      e.n <- paste(cox.stats$nevent, "/", cox.stats$n)
      hr.ci <- cox.stats$output %>%
        magrittr::extract(var.idx, c("estimate", "conf.low", "conf.high")) %>%
        format_hr_ci(digits = 2, labels = FALSE, method = "Sci") %>%
        paste0(ifelse(cox.stats$used.firth, firth.caption, "")) %>%
        paste(collapse = kLocalConstantHrSepFlag)
      if (length(curr.var.names) == 1) {
        p.value <- stats::anova(cox.stats$fit)[2, "Pr(>|Chi|)"]
      } else {
        p.value <- stats::anova(
          coxph(stats::as.formula(paste(surv.formula, paste(curr.var.names[-length(curr.var.names)], collapse = "+"))), data = temp.d),
          cox.stats$fit)[2, "P(>|Chi|)"]
      }
      p.value <- round_pval(p.value, round.small = round.small,
                            scientific = scientific,
                            digits = round.digits.p.value)
      result.table[nvar * (j - 1) + i, ] <- c(e.n, hr.ci, p.value)
    }
    cox.stats.output[[surv.descriptions[j]]] <- cox.stats # only capture the full model i.e. the last one
  }

  ### generate word-friendly table via pander i.e. result.table.bamboo ... ###
  result.table.bamboo <- result.table
  result.table.ncol <- ncol(result.table)
  result.table.bamboo.base.indexes <- c() # base indexes for each survival end point in result.table.bamboo
  num.surv <- length(surv.descriptions) # number of survival end points
  num.var <- length(var.descriptions) # number of variables
  for (i in 1:num.surv) {
    result.table.bamboo.base.index <- 1 + (i - 1) * (num.var + 1)
    if (i == 1) {
      result.table.bamboo <- rbind(rep("", result.table.ncol), result.table.bamboo)
    } else {
      result.table.bamboo <- rbind(
          result.table.bamboo[1:(result.table.bamboo.base.index - 1), ],
          rep("", result.table.ncol),
          result.table.bamboo[result.table.bamboo.base.index:nrow(result.table.bamboo), ])
    }
    rownames(result.table.bamboo)[result.table.bamboo.base.index] <- paste0("**", surv.descriptions[i], "**")
    rownames(result.table.bamboo)[result.table.bamboo.base.index + c(1:num.var)] <- var.descriptions
    # want to show # of events only once for each surv endpoint
    result.table.bamboo[result.table.bamboo.base.index, 1] <- result.table.bamboo[result.table.bamboo.base.index + 1, 1]
    result.table.bamboo[result.table.bamboo.base.index + c(1:num.var), 1] <- ""
    result.table.bamboo.base.indexes <- c(result.table.bamboo.base.indexes, result.table.bamboo.base.index)
  }
  # want to add a column to describe different factor level for categorical
  # whenever reference group is specified
  if (sum(is.na(var.ref.groups)) != length(var.ref.groups)) {
    first.col.name <- colnames(result.table.bamboo)[1]
    result.table.bamboo <- cbind(result.table.bamboo[, 1], "", result.table.bamboo[, 2:3])
    colnames(result.table.bamboo)[1] <- first.col.name
    hr.col.index <- 3 # column with the hazard ratios
    for (i in 1:num.surv) {
      result.table.bamboo.base.index <- result.table.bamboo.base.indexes[i]
      rows.added <- 0
      for (var.count in 1:length(var.names)) {
        if (!is.na(var.ref.groups[var.count])) {
          ref.group <- var.ref.groups[var.count]
          other.groups <- names(table(input.d[, var.names[var.count]]))
          other.groups <- other.groups[other.groups != ref.group & !(other.groups %in% missing.codes)]
          num.other.groups <- length(other.groups)
          curr.base.index <- result.table.bamboo.base.index + (var.count - 1) + rows.added + 1
          if (num.other.groups > 1) {
            for (j in 1:(num.other.groups - 1)) {
              if (curr.base.index < nrow(result.table.bamboo)) {
                last.row.name <- rownames(result.table.bamboo)[nrow(result.table.bamboo)]
                result.table.bamboo <- rbind(
                    result.table.bamboo[1:curr.base.index, ],
                    rep("", ncol(result.table.bamboo)),
                    result.table.bamboo[(curr.base.index + 1):nrow(result.table.bamboo), ])
                rownames(result.table.bamboo)[nrow(result.table.bamboo)] <- last.row.name
              } else {
                result.table.bamboo <- rbind(
                    result.table.bamboo[1:curr.base.index, ],
                    rep("", ncol(result.table.bamboo))
                )
              }
              rows.added <- rows.added + 1
            }
          }
          if (num.other.groups > 1 | show.group.name.for.bin.var) {
            result.table.bamboo[curr.base.index:(curr.base.index + num.other.groups - 1), hr.col.index] <-
                strsplit(result.table.bamboo[curr.base.index, hr.col.index], kLocalConstantHrSepFlag)[[1]]
            result.table.bamboo[curr.base.index:(curr.base.index + num.other.groups - 1), hr.col.index - 1] <- other.groups
          }
        }
      }

      # need to update result.table.bamboo.base.indexes since we've added rows!!!
      if (i < num.surv) {
        result.table.bamboo.base.indexes[(i + 1):num.surv] <- result.table.bamboo.base.indexes[(i + 1):num.surv] + rows.added
      }
    }
  }

  # subscript ("<sup>|</sup>") and line break ("<br>") syntax for pandoc
  options("table_counter" = options()$table_counter - 1)
  result.table.bamboo <- result.table.bamboo %>%
    gsub(pattern = "<sup>|</sup>", replacement = "^", .) %>%
    pander::pandoc.table.return(., caption = ifelse(is.na(caption), "", paste0("*", addTableNumber(caption), "*")),
                                emphasize.rownames = FALSE,
                                split.table = split.table, ...) %>%
    gsub(pattern = kLocalConstantHrSepFlag, replacement = "; ", .) %>%
    gsub(pattern = "<br>", replacement = "\\\\\n", .)
  ### end of result.table.bamboo ###

  ### generate html table ... ###
  result.table.html <- result.table.bamboo # just set it the same as bamboo
  ### end of generate html table ###

  ### clean result.table ###
  result.table <- gsub(kLocalConstantHrSepFlag, ", ", result.table)
  ### end of clean result.table ###

  return(list("result.table" = result.table,
              "result.table.html" = result.table.html,
              "result.table.bamboo" = result.table.bamboo,
              "cox.stats" = cox.stats.output))
}
