#' Generate cohort characteristics
#' @param input.d The `data.frame` containing the data
#' @param marker.name The variable that you want to split into different columns
#' @param marker.description The description for the variable(s) to split
#' @param var.names The variables that you want the statistics for
#' @param is.var.continuous Vector of length equal to the length of var.names
#'   with 1 indicating a continuous variable and 0 otherwise (this should be
#'   inferred in the function)
#' @param var.descriptions Vector of strings to describe the variables as they
#'   are to appear in the table
#' @param marker.value.labels.tolower Indicator as to whether to put marker
#'   value labels to lower case
#' @param show.missing an indicator to whether to show missing values
#' @param show.missing.continuous if set to `FALSE` and `show.missing
#'   == FALSE`, will not show the number of missing cases for continuous
#'   variables. Otherwise, it shows the number of missing for continuous
#'   variables even if `show.missing == FALSE`.
#' @param do.droplevels drop categories of unobserved levels set to `TRUE`
#' @param show.percent defaults to "both" which shows both rows and columns
#'   other possible values: "column", "row".
#' @param stat.tests statistical test to perform. `NULL` indicates do not
#'   do test for all variables, `NA` indicates do not do test for specified
#'   variable. Tests: chisq, fisher, ttest, wilcox, kendall, spearman, pearson,
#'   kruskal, confusionMarkerAsRef, confusionVarAsRef
#' @param chisq.test.simulate.p.value Whether to simulate p-value for chi-square
#'   test. this parameter is ignored if chi-square is not used.  Default
#'   value=FALSE
#' @param stat.test.column.header The name to show on the header defaults to
#'   "association/correlation test"
#' @param round.digits.p.value The number of digits to round the P values
#' @param num.boot the number of bootstrap samples for any bootstrap method that
#'   may be used
#' @param missing.codes.highlight default to `NULL` this indicates whether
#'   we wanted the missing values broken down down or lumped together.
#' @param missing.codes a vector to indicate how missing values are coded,
#'   default is `c("N/A", "", "Unk")`
#' @param decimal number of decimal places to show for aggregate numbers such as
#'   proportions or averages; default to 0
#' @param caption caption to use for the Table
#' @param html.table.border the border type to use for html tables
#' @param banded.rows If `TRUE`, rows have alternating shading colour
#' @param css.class.name.odd Used to set the row colour for odd rows
#' @param css.class.name.even Used to set the row colour for even rows
#' @param custom.marker.labels labels of marker to show; default `NULL`
#'   means using existing value label of the marker
#' @param custom.total.label label of the "Total" column; default `NULL`
#'   means show "Total"
#' @param split.table number of chars per row before table is split.
#' @param ... additional arguments to `pander`
#' @return A table with statistics reported for multiple variables, such as
#'   mean, median, and range for continuous variables and proportions and
#'   percentages for categorical variables. Relevant association and correlation
#'   tests are performed as well.
#' @author Aline Talhouk
#' @export
#' @examples
#' dcc <- doCohortCharacteristics( input.d = mtcars, marker.name = "cyl",
#' marker.description = "cylinders", var.names = c("disp", "hp"),
#' var.descriptions = c("displacement", "horsepower"), is.var.continuous =
#' c(TRUE, TRUE), caption = "Some mtcars summaries")
#' htmlTable::htmlTable(dcc$result.table.html)
doCohortCharacteristics <- function(input.d, marker.name, marker.description,
                                    var.names, is.var.continuous,
                                    var.descriptions,
                                    marker.value.labels.tolower = TRUE,
                                    show.missing = TRUE,
                                    show.missing.continuous = TRUE,
                                    do.droplevels = TRUE,
                                    show.percent = "both",
                                    stat.tests = NULL,
                                    chisq.test.simulate.p.value = FALSE,
                                    stat.test.column.header =
                                      "association/correlation test",
                                    round.digits.p.value = 4,
                                    num.boot = 1000,
                                    missing.codes.highlight = NULL,
                                    missing.codes = c("N/A", "", "Unk"),
                                    decimal = 0, caption = NA,
                                    html.table.border = 0,
                                    banded.rows = FALSE,
                                    css.class.name.odd = "odd",
                                    css.class.name.even = "even",
                                    custom.marker.labels = NULL,
                                    custom.total.label = NULL,
                                    split.table = 200, # set default for pander
                                    ...) {
  kLocalConstantRowColPercentBeginFlag <- "kLocalConstantRowColPercentBeginFlag"
  kLocalConstantRowColPercentEndFlag <- "kLocalConstantRowColPercentEndFlag"
  kLocalConstantRowColPercentSepFlag <- "kLocalConstantRowColPercentSepFlag"
  kLocalConstantItalicBeginFlag <- "kLocalConstantItalicBeginFlag"
  kLocalConstantItalicEndFlag <- "kLocalConstantItalicEndFlag"
  kLocalConstantStatTestBeginFlag <- "kLocalConstantStatTestBeginFlag"
  kLocalConstantNewLineFlag <- "kLocalConstantNewLineFlag"
  col.th.style <- "border-bottom: 1px solid grey; border-top: 4px double grey; text-align: center; padding-right:10px; padding-right:10px;"
  row.th.style <- "text-align: left; padding-right:10px; padding-right:10px;"

  if (is.null(missing.codes.highlight)) {
    missing.codes.highlight <- list()
    for (var.name in var.names) {
      missing.codes.highlight[var.name] <- NA
    }
  }

  # input.d.no.missing DEFINED HERE!!!

  #Filter Missing Value for the marker but other missing values may still exist
  input.d.no.missing <- input.d[!is.na(input.d[, marker.name]) & !input.d[, marker.name] %in% missing.codes, ]

  # Droplevels
  if (is.factor(input.d.no.missing[, marker.name]) & do.droplevels) {
    input.d.no.missing[, marker.name] <- droplevels(input.d.no.missing[, marker.name])
  }
  # Formatting Labels
  marker.categories <- names(table(input.d.no.missing[, marker.name]))
  marker.categories <- marker.categories[!marker.categories %in% missing.codes]

  total.label <- ifelse(sum(sapply(var.descriptions, function(x){
    return(ifelse(isFirstLetterUpperCase(x), 1, 0))})) == length(var.descriptions), "Total", "total")

  if (!is.null(custom.total.label)) {
    total.label <- custom.total.label
  }

  if (marker.value.labels.tolower) {
    result.table.col.names <- c(total.label, paste(marker.description, tolower(marker.categories)))
  } else {
    result.table.col.names <- c(total.label, paste(marker.description, marker.categories))
  }

  if (!is.null(custom.marker.labels)) {
    custom.marker.labels.length <- length(custom.marker.labels)
    result.table.col.names.length <- length(result.table.col.names)
    if (result.table.col.names.length == (custom.marker.labels.length + 1)) {
      result.table.col.names[2:result.table.col.names.length] <- custom.marker.labels
    } else if (result.table.col.names.length < (custom.marker.labels.length + 1)) {
      result.table.col.names[2:result.table.col.names.length] <- custom.marker.labels[c(1:(result.table.col.names.length - 1))]
    } else {
      result.table.col.names[2:result.table.col.names.length] <- custom.marker.labels[
        rep(c(1:custom.marker.labels.length), floor(result.table.col.names.length / custom.marker.labels.length)),
        c(1:(result.table.col.names.length %% custom.marker.labels.length))]
    }
  }

  # If all items in var.descriptions starts with capital letter, first row header should be capital as well i.e. "Total"
  result.table.row.names <- total.label
  result.table <- c(paste0(nrow(input.d), " (100%)"),
                    sapply(marker.categories, function(x){
                      return(paste0(sum(input.d.no.missing[, marker.name] == x), " (",
                                   format(round(sum(input.d.no.missing[, marker.name] == x) / nrow(input.d.no.missing) * 100, decimal),
                                          nsmall = decimal), "%)"))
                    })
  )

  do.stats <- !is.null(stat.tests)
  stat.tests.results <- c()
  num.var <- length(var.names)

  # Loop over all variables
  for (i in 1:num.var) {
    #############################
    # NOTE: three data matrices here ...
    # input.d.no.missing - no missing marker only (may contain missing var)
    # input.d.no.missing.var.only - missing var only (may contain missing marker)
    # input.d.no.missing.var - no missing marker and no missing var
    #############################

    # Header row for each var
    result.table <- rbind(result.table, c("", rep("", length(marker.categories))))
    result.table.row.names <- c(result.table.row.names, var.descriptions[i])
    var.name <- var.names[i]
    var.description <- var.descriptions[i]

    num.missing.row.header.name <- c()
    if (sum(!is.na(missing.codes.highlight[[var.name]])) > 0) {
      num.missing.row.header.name <- c(num.missing.row.header.name, missing.codes.highlight[[var.name]]) # need double [[]] for missing.codes.highlight because its a list
    }
    num.missing.row.header.name <- c(num.missing.row.header.name, "missing")

    # input.d.no.missing.var.only DEFINED HERE!!!
    input.d.no.missing.var.only <- input.d[!is.na(input.d[, var.name]) & !input.d[, var.name] %in%
                                             missing.codes & !input.d[, var.name] %in% missing.codes.highlight[[var.name]], ]
    if (is.factor(input.d.no.missing.var.only[, var.name]) & do.droplevels) {
      input.d.no.missing.var.only[, var.name] <- droplevels(input.d.no.missing.var.only[, var.name])
    }

    # input.d.no.missing.var DEFINED HERE!!!
    input.d.no.missing.var <- input.d.no.missing.var.only[!is.na(input.d.no.missing.var.only[, marker.name]) &
                                                            !input.d.no.missing.var.only[, marker.name] %in% missing.codes, ]
    if (is.factor(input.d.no.missing.var[, marker.name]) & do.droplevels) {
      input.d.no.missing.var[, marker.name] <- droplevels(input.d.no.missing.var[, marker.name])
    }
    # If continuous variable
    if (is.var.continuous[i]) {
      input.d.no.missing.var.only[, var.name] <- as.numeric(input.d.no.missing.var.only[, var.name])
      input.d.no.missing.var[, var.name] <- as.numeric(input.d.no.missing.var[, var.name])

      # 4 rows: mean (+/- std dev) / median / IQR / number of missing
      var.row.names <- c("mean (SD)", "median", "IQR", "range")
      if (show.missing | show.missing.continuous) {
        var.row.names <- c(var.row.names, num.missing.row.header.name)
      }
      result.table.row.names <- c(result.table.row.names, var.row.names)

      # Statistical testing for continuous variables (ignore non-applicable tests!) ...
      stat.test.result <- NA
      if (do.stats) {
        switch(stat.tests[i],
               spearman = {
                 spearman.result <- stats::cor.test(input.d.no.missing.var[, var.name],
                                                    as.numeric(input.d.no.missing.var[, marker.name]),
                                                    method = "spearman")
                 stat.test.result <- paste0("Spearman correlation", kLocalConstantStatTestBeginFlag,
                                            "rho = ", round(spearman.result$estimate, 2),
                                            kLocalConstantNewLineFlag, "P = ",
                                            sprintf(paste0("%.", round.digits.p.value, "f"),
                                                    round(spearman.result$p.value, digits = round.digits.p.value)))
               },
               kruskal = {
                 kruskal.result <- stats::kruskal.test(input.d.no.missing.var[, var.name] ~
                                                         as.numeric(input.d.no.missing.var[, marker.name]))
                 stat.test.result <- paste0("Kruskal-Wallis rank sum test", kLocalConstantStatTestBeginFlag,
                                            kLocalConstantNewLineFlag, "P = ",
                                            sprintf(paste0("%.", round.digits.p.value, "f"),
                                                    round(kruskal.result$p.value, digits = round.digits.p.value)))
               },
               wilcox = {
                 wilcox.result <- stats::wilcox.test(input.d.no.missing.var[, var.name] ~
                                                       as.numeric(input.d.no.missing.var[, marker.name]))
                 stat.test.result <- paste0("Wilcoxon rank sum test", kLocalConstantStatTestBeginFlag, "P = ",
                                            sprintf(paste0("%.", round.digits.p.value, "f"),
                                                    round(wilcox.result$p.value, digits = round.digits.p.value)))
               }
        )
        stat.tests.results <- c(stat.tests.results, stat.test.result)
      }

      # Assemble Continuous Results
      result.table <- rbind(
        result.table, # Add the mean
        c(paste0(format(round(mean(input.d.no.missing.var.only[, var.name]), decimal), nsmall = decimal),
                " (", format(round(stats::sd(input.d.no.missing.var.only[, var.name]), decimal), nsmall = decimal), ")"),
          sapply(marker.categories, function(x) {
            temp.d <- input.d.no.missing.var[input.d.no.missing.var[, marker.name] == x, var.name]
            if (length(temp.d) == 0) {
              return(MISSING.EXPLICIT)
            } else {
              return(paste0(
                format(round(mean(temp.d), decimal), nsmall = decimal),
                " (", format(round(stats::sd(temp.d), decimal), nsmall = decimal), ")"))
            }
          })
        ),
        # Add the median
        c(format(round(stats::median(input.d.no.missing.var.only[, var.name]), decimal), nsmall = decimal),
          sapply(marker.categories, function(x) {
            temp.d <- input.d.no.missing.var[input.d.no.missing.var[, marker.name] == x, var.name]
            if (length(temp.d) == 0) {
              return(MISSING.EXPLICIT)
            } else {
              return(format(round(stats::median(temp.d), decimal), nsmall = decimal))
            }
          })
        ),
        # Add inter quartile range
        c(paste(format(round(stats::quantile(input.d.no.missing.var.only[, var.name], c(0.25, 0.75)),
                             decimal), nsmall = decimal), collapse = " to "),
          sapply(marker.categories, function(x) {
            temp.d <- input.d.no.missing.var[input.d.no.missing.var[, marker.name] == x, var.name]
            if (length(temp.d) == 0) {
              return(MISSING.EXPLICIT)
            } else {
              return(paste(format(round(stats::quantile(temp.d, c(0.25, 0.75)), decimal), nsmall = decimal), collapse = " to "))
            }
          })
        ),
        # Add range
        c(paste(format(round(range(input.d.no.missing.var.only[, var.name]),
                             decimal), nsmall = decimal), collapse = " to "),
          sapply(marker.categories, function(x) {
            temp.d <- input.d.no.missing.var[input.d.no.missing.var[, marker.name] == x, var.name]
            if (length(temp.d) == 0) {
              return(MISSING.EXPLICIT)
            } else {
              return(paste(format(round(range(temp.d), decimal), nsmall = decimal), collapse = " to "))
            }
          })
        )
      )
      # If categorical variable
    } else {
      var.categories <- names(table(input.d.no.missing.var.only[, var.name]))
      var.row.names <- var.categories
      if (show.missing) {
        var.row.names <- c(var.row.names, num.missing.row.header.name)
      }
      result.table.row.names <- c(result.table.row.names, var.row.names)

      # do statistical test for continuous variables (ignore non-applicable tests!) ...
      stat.test.result <- NA
      if (do.stats) {
        switch(stat.tests[i],
               kendall = {
                 kendall.result <- stats::cor.test(as.numeric(input.d.no.missing.var[, var.name]),
                                                   as.numeric(input.d.no.missing.var[, marker.name]), method = "kendall")
                 stat.test.result <- paste0("Kendall correlation", kLocalConstantStatTestBeginFlag,
                                            "tau = ", round(kendall.result$estimate, 2),
                                            kLocalConstantNewLineFlag, "P = ",
                                            sprintf(paste0("%.", round.digits.p.value, "f"),
                                                    round(kendall.result$p.value, digits = round.digits.p.value)))
               },
               chisq = {
                 chisq.result <- stats::chisq.test(table(input.d.no.missing.var[, var.name],
                                                         input.d.no.missing.var[, marker.name]), simulate.p.value = chisq.test.simulate.p.value)
                 stat.test.result <- paste0("Chi-square test", kLocalConstantStatTestBeginFlag,
                                            "P = ", sprintf(paste0("%.", round.digits.p.value, "f"),
                                                            round(chisq.result$p.value, digits = round.digits.p.value)))
               },
               fisher = {
                 fisher.result <- stats::fisher.test(table(input.d.no.missing.var[, var.name],
                                                           input.d.no.missing.var[, marker.name]), workspace = 2e6)
                 stat.test.result <- paste0("Fisher's exact test", kLocalConstantStatTestBeginFlag, "P = ",
                                            sprintf(paste0("%.", round.digits.p.value, "f"),
                                                    round(fisher.result$p.value, digits = round.digits.p.value)))
               },
               confusionMarkerAsRef = {
                 # confusion matrix, marker as the reference
                 # require both marker and var to be factor ...
                 # if not, just print err msg
                 if (!is.factor(input.d.no.missing.var[, var.name]) | !is.factor(input.d.no.missing.var[, marker.name])) {
                   stat.test.result <- "error: both marker and variable needs to be factor"
                 } else {
                   stat.test.result <- binaryCMAsHTML(as.numeric(input.d.no.missing.var[, marker.name]),
                                                      as.numeric(input.d.no.missing.var[, var.name]),
                                                      marker.description, round.digits.p.value,
                                                      num.boot = num.boot)
                 }
               },
               confusionVarAsRef = {
                 # confusion matrix, variable as the reference
                 # require both marker and var to be factor ...
                 # if not, just print err msg
                 if (!is.factor(input.d.no.missing.var[, var.name]) | !is.factor(input.d.no.missing.var[, marker.name])) {
                   stat.test.result <- "error: both marker and variable needs to be factor"
                 } else {
                   stat.test.result <- binaryCMAsHTML(as.numeric(input.d.no.missing.var[, var.name]),
                                                      as.numeric(input.d.no.missing.var[, marker.name]),
                                                      var.description, round.digits.p.value,
                                                      num.boot = num.boot)
                 }
               }
        )
        stat.tests.results <- c(stat.tests.results, stat.test.result)
      }

      for (var.category in var.categories) {
        total.value <- paste0(sum(input.d.no.missing.var.only[, var.name] == var.category), " (",
                              format(round(sum(input.d.no.missing.var.only[, var.name] == var.category) /
                                             nrow(input.d.no.missing.var.only) * 100, decimal), nsmall = decimal), "%)")
        result.table <- rbind(
          result.table,
          switch(show.percent,
                 row = {
                   c(total.value, sapply(marker.categories, function(x) {
                     return(paste0(sum(input.d.no.missing.var[, var.name] == var.category &
                                         input.d.no.missing.var[, marker.name] == x), " (",
                                   ifelse(sum(input.d.no.missing.var[, var.name] == var.category) > 0,
                                          format(round(sum(input.d.no.missing.var[, var.name] == var.category &
                                                             input.d.no.missing.var[, marker.name] == x) /
                                                         sum(input.d.no.missing.var[, var.name] == x) * 100, decimal), nsmall = decimal), 0),
                                   "%)"))
                   }))
                 },
                 column = {
                   c(total.value, sapply(marker.categories, function(x) {
                     return(paste0(sum(input.d.no.missing.var[, var.name] == var.category &
                                         input.d.no.missing.var[, marker.name] == x), " (",
                                   ifelse(sum(input.d.no.missing.var[, marker.name] == x) > 0,
                                          format(round(sum(input.d.no.missing.var[, var.name] == var.category &
                                                             input.d.no.missing.var[, marker.name] == x) /
                                                         sum(input.d.no.missing.var[, marker.name] == x) * 100, decimal), nsmall = decimal), 0),
                                   "%)"))
                   }))
                 },
                 both = {
                   c(total.value, sapply(marker.categories, function(x) {
                     return(paste0(sum(input.d.no.missing.var[, var.name] == var.category &
                                         input.d.no.missing.var[, marker.name] == x),
                                   kLocalConstantRowColPercentBeginFlag,
                                   ifelse(sum(input.d.no.missing.var[, var.name] == var.category) > 0,
                                          format(round(sum(input.d.no.missing.var[, var.name] == var.category &
                                                             input.d.no.missing.var[, marker.name] == x) /
                                                         sum(input.d.no.missing.var[, var.name] == var.category) * 100, decimal), nsmall = decimal), 0),
                                   "%", kLocalConstantRowColPercentSepFlag,
                                   ifelse(sum(input.d.no.missing.var[, marker.name] == x) > 0,
                                          format(round(sum(input.d.no.missing.var[, var.name] == var.category &
                                                             input.d.no.missing.var[, marker.name] == x) /
                                                         sum(input.d.no.missing.var[, marker.name] == x) * 100, decimal), nsmall = decimal), 0),
                                   "%", kLocalConstantRowColPercentEndFlag))
                   }))
                 }
          )
        )
      }
    }

    if (show.missing | is.var.continuous[i] & show.missing.continuous) {
      if (sum(!is.na(missing.codes.highlight[[var.name]])) > 0) {
        # there's some missing values we want to highlight ...
        for (missing.code in missing.codes.highlight[[var.name]]) {
          result.table <- rbind(result.table,
                                c(sum(!is.na(input.d[, var.name]) & input.d[, var.name] == missing.code), # number of missing
                                  sapply(marker.categories, function(x) {
                                    temp.d <- input.d.no.missing[input.d.no.missing[, marker.name] == x, var.name]
                                    return(sum(!is.na(temp.d) & temp.d == missing.code))
                                  })
                                )
          )
        }
      }
      result.table <- rbind(result.table,
                            c(sum(is.na(input.d[, var.name]) | input.d[, var.name] %in% missing.codes), # number of missing
                              sapply(marker.categories, function(x) {
                                temp.d <- input.d.no.missing[input.d.no.missing[, marker.name] == x, var.name]
                                return(sum(is.na(temp.d) | temp.d %in% missing.codes))
                              })
                            )
      )
    }
  }
  colnames(result.table) <- result.table.col.names
  rownames(result.table) <- result.table.row.names

  ### end of "data" generation ###

  ##################################
  ### generate html table ...	 ###
  result.table.html <- paste0("<table border=", html.table.border, ">",
                              ifelse(is.na(caption), "",
                                     paste0("<caption style='", TABLE.CAPTION.STYLE, "'>", caption, "</caption>")))
  result.table.html <- paste0(
    result.table.html,
    "<tr><th style='", col.th.style, "' colspan=2></th><th style='", col.th.style, "'>",
    paste(result.table.col.names, collapse = paste0("</th><th style='", col.th.style, "'>")),
    "</th>",
    ifelse(do.stats, paste0("<th style='", col.th.style, "'>", stat.test.column.header, "</th>"), ""),
    "</tr>")

  row.band.toggle <- TRUE
  var.count <- -1 # want to skip the header row which contains the total count
  for (i in 1:nrow(result.table)) {
    num.missing.row.header.name <- c()
    if (var.count > 0) {
      var.name <- var.names[var.count]
      if (sum(!is.na(missing.codes.highlight[[var.name]])) > 0) {
        num.missing.row.header.name <- c(num.missing.row.header.name, missing.codes.highlight[[var.name]]) # need double [[]] for missing.codes.highlight because its a list
      }
    }
    num.missing.row.header.name <- c(num.missing.row.header.name, "missing")

    # check to see if this is the start of a new category or 1st row ...
    var.header.row <- FALSE
    if (i == 1 | sum(result.table[i, ] == "", na.rm = TRUE) > 1) {
      # this must be the start of a first category since row contains >1 empty cells (there will be 1 empty cell for any var with no stats test)
      var.header.row <- TRUE
      var.count <- var.count + 1
      # varname not available for i==1 (the row for the total)
      if (i > 1) {
        var.name <- var.names[var.count]
      }
      row.band.toggle <- !row.band.toggle
      tr.class <- ifelse(banded.rows, paste0(" class='",
                                             ifelse(row.band.toggle, css.class.name.even, css.class.name.odd), "'"), "")
      result.table.html <- paste0(
        result.table.html,
        "<tr", tr.class, "><th style='", row.th.style, "' colspan=", 2, ">", result.table.row.names[i], "</th>")
    } else {
      if (num.var == 1) {
        row.band.toggle <- !row.band.toggle # always toggle when there's only one variable
      }
      tr.class <- ifelse(banded.rows, paste0(" class='", ifelse(row.band.toggle, css.class.name.even, css.class.name.odd), "'"), "")
      result.table.html <- paste0(
        result.table.html,
        "<tr", tr.class, "><th>&nbsp;&nbsp;&nbsp;&nbsp;</th><th style='", row.th.style, "'>",
        result.table.row.names[i],
        ifelse(!is.var.continuous[var.count] & show.percent == "both",
               ifelse(!result.table.row.names[i] %in% num.missing.row.header.name | !show.missing,
                      "<i><br>&nbsp;&nbsp;&nbsp;&nbsp;row%<br>&nbsp;&nbsp;&nbsp;&nbsp;col%</i>", ""),
               ""),
        "</th>")
    }

    result.table.html <- paste0(
      result.table.html,
      "<td>", paste(
        gsub(kLocalConstantRowColPercentBeginFlag, "<br><i>",
             gsub(kLocalConstantRowColPercentEndFlag, "</i>",
                  gsub(kLocalConstantRowColPercentSepFlag, "<br>", result.table[i, ]))),
        collapse = "</td><td>"
      ), "</td>",
      ifelse(
        (var.header.row & do.stats & var.count > 0) | i == 1,
        ifelse(i == 1,
               "<td></td>",
               paste0(
                 "<td rowspan=",
                 ifelse(
                   is.var.continuous[var.count],
                   ifelse(show.missing, 1 + sum(!is.na(missing.codes.highlight[[var.name]])), 0) + 4,
                   ifelse(show.missing, 1 + sum(!is.na(missing.codes.highlight[[var.name]])), 0) + 1 + as.numeric(ifelse(
                     is.factor(input.d.no.missing[, var.name]),
                     length(names(table(droplevels(input.d.no.missing.var.only[!is.na(input.d.no.missing.var.only[, var.name]) &
                                                                                 !input.d.no.missing.var.only[, var.name] %in%
                                                                                 c(missing.codes, missing.codes.highlight[[var.name]]), var.name])))),
                     length(names(table(input.d.no.missing.var.only[!is.na(input.d.no.missing.var.only[, var.name]) &
                                                                      !input.d.no.missing.var.only[, var.name] %in%
                                                                      c(missing.codes, missing.codes.highlight[[var.name]]), var.name])))
                   ))
                 ),
                 ">", gsub(kLocalConstantNewLineFlag, "<br>",
                           gsub(kLocalConstantStatTestBeginFlag, "<br>",
                                stat.tests.results[var.count])), "</td>")
        ),
        ""
      ),
      "</tr>")
  }

  result.table.html <- paste0(result.table.html, "</table>")
  ### end of generate html table ###
  ##################################
  options("table_counter" = options()$table_counter - 1)
  ##################################
  ### generate table for pander  ###
  result.table.bamboo <- result.table
  for (i in 1:nrow(result.table.bamboo)) {
    for (j in 1:ncol(result.table.bamboo)) {
      result.table.bamboo[i, j] <- gsub(kLocalConstantRowColPercentBeginFlag, "(*",
                                        gsub(kLocalConstantRowColPercentEndFlag, "*)",
                                             gsub(kLocalConstantRowColPercentSepFlag,
                                                  "*, *", result.table.bamboo[i, j])))
    }
  }
  num.col.in.result.table.bamboo <- ncol(result.table.bamboo)
  # add stat column
  if (do.stats) {
    result.table.bamboo <- cbind(result.table.bamboo, "")
    num.col.in.result.table.bamboo <- num.col.in.result.table.bamboo + 1
    colnames(result.table.bamboo)[num.col.in.result.table.bamboo] <- stat.test.column.header
  }
  var.count <- 1
  for (i in 1:nrow(result.table)) {
    var.header.row <- FALSE
    if (i == 1 | sum(result.table[i, ] == "", na.rm = TRUE) > 1) {
      # this must be the start of a first category since row contains >1 empty cells (there will be 1 empty cell for any var with no stats test)
      var.header.row <- TRUE
      rownames(result.table.bamboo)[i] <- paste0("**", rownames(result.table)[i], "**") # make it bold
      if (i > 1) {
        if (do.stats) {
          result.table.bamboo[i, num.col.in.result.table.bamboo] <-
            gsub(kLocalConstantNewLineFlag, ", ",
                 gsub(kLocalConstantStatTestBeginFlag, ": ",
                      stat.tests.results[var.count]))
        }
        var.count <- var.count + 1
      }
    }
  }
  result.table.bamboo <- pander::pandoc.table.return(result.table.bamboo, emphasize.rownames = FALSE,
                                                     split.table = split.table, ...,
                                                     caption = caption)

  ##############################################################################
  ### remove the indicator flags e.g. kLocalConstantRowColPercentEndFlag     ###
  ### from result.table so that it can be viewed in R by others without so   ###
  ### many clutter.                                                          ###
  ##############################################################################
  for (i in 1:nrow(result.table)) {
    for (j in 1:ncol(result.table)) {
      result.table[i, j] <- gsub(
        kLocalConstantRowColPercentBeginFlag, "(",
        gsub(kLocalConstantRowColPercentEndFlag, ")",
             gsub(kLocalConstantRowColPercentSepFlag, ", ",
                  result.table[i, j])))
    }
  }
  for (i in 1:length(stat.tests.results)) {
    stat.tests.results[i] <- gsub(kLocalConstantNewLineFlag, ", ",
                                  gsub(kLocalConstantStatTestBeginFlag, ": ",
                                       stat.tests.results[i]))
  }

  return(list("result.table" = result.table,
              "stat.tests.results" = stat.tests.results,
              "result.table.html" = result.table.html,
              "result.table.bamboo" =  result.table.bamboo))
}
