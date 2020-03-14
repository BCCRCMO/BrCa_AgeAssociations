smCrossTable <- 
function (x, y, digits = 3, max.width = 5, expected = FALSE, 
    prop.r = TRUE, prop.c = TRUE, prop.t = TRUE, prop.chisq = TRUE, 
    chisq = FALSE, fisher = FALSE, mcnemar = FALSE, resid = FALSE, 
    sresid = FALSE, asresid = FALSE, missing.include = FALSE, 
    format = c("SAS", "SPSS"), dnn = NULL, ...) 
{
  require("gmodels")
    format = match.arg(format)
    RowData <- deparse(substitute(x))
    if (!missing(y)) 
        ColData <- deparse(substitute(y))
    if (max.width < 1) 
        stop("max.width must be >= 1")
    vector.x <- FALSE
    if (expected) 
        chisq <- TRUE
    if (missing(y)) {
        if (is.null(dim(x))) {
            if (missing.include) 
                x <- factor(x, exclude = NULL)
            else x <- factor(x)
            t <- t(as.matrix(table(x)))
            vector.x <- TRUE
        }
        else if (length(dim(x) == 2)) {
            if (any(x < 0) || any(is.na(x))) 
                stop("all entries of x must be nonnegative and finite")
            if (is.null(names(dimnames(x)))) {
                RowData <- ""
                ColData <- ""
            }
            else {
                RowData <- names(dimnames(x))[1]
                ColData <- names(dimnames(x))[2]
            }
            if (is.null(rownames(x))) 
                rownames(x) <- paste("[", 1:nrow(x), ",]", sep = "")
            if (is.null(colnames(x))) 
                colnames(x) <- paste("[,", 1:ncol(x), "]", sep = "")
            t <- x
        }
        else stop("x must be either a vector or a 2 dimensional matrix, if y is not given")
    }
    else {
        if (length(x) != length(y)) 
            stop("x and y must have the same length")
        if (missing.include) {
            x <- factor(x, exclude = c())
            y <- factor(y, exclude = c())
        }
        else {
            x <- factor(x)
            y <- factor(y)
        }
        t <- table(x, y)
    }
    if (all(dim(t) >= 2)) {
        if (!is.null(dnn)) {
            if (length(dnn) != 2) 
                stop("dnn must have length of 2, one element for each table dimension")
            else {
                RowData <- dnn[1]
                ColData <- dnn[2]
            }
        }
    }
    if (any(dim(t) < 2)) {
        prop.c <- prop.r <- prop.chisq <- chisq <- expected <- fisher <- mcnemar <- FALSE
    }
    CPR <- prop.table(t, 1)
    CPC <- prop.table(t, 2)
    CPT <- prop.table(t)
    GT <- sum(t)
    RS <- rowSums(t)
    CS <- colSums(t)
    if (length(dim(x) == 2)) 
        TotalN <- GT
    else TotalN <- length(x)
    ColTotal <- "Column Total"
    RowTotal <- "Row Total"
    CWidth <- max(digits + 2, c(nchar(t), nchar(dimnames(t)[[2]]), 
        nchar(RS), nchar(CS), nchar(RowTotal)))
    RWidth <- max(c(nchar(dimnames(t)[[1]]), nchar(ColTotal)))
    if (exists("RowData")) 
        RWidth <- max(RWidth, nchar(RowData))
    RowSep <- paste(rep("-", CWidth + 2), collapse = "")
    RowSep1 <- paste(rep("-", RWidth + 1), collapse = "")
    SpaceSep1 <- paste(rep(" ", RWidth), collapse = "")
    SpaceSep2 <- paste(rep(" ", CWidth), collapse = "")
    FirstCol <- formatC(dimnames(t)[[1]], width = RWidth, format = "s")
    ColTotal <- formatC(ColTotal, width = RWidth, format = "s")
    RowTotal <- formatC(RowTotal, width = CWidth, format = "s")
    if (chisq) {
        if (all(dim(t) == 2)) 
            CSTc <- chisq.test(t, correct = TRUE, ...)
        CST <- chisq.test(t, correct = FALSE, ...)
    }
    else CST <- suppressWarnings(chisq.test(t, correct = FALSE))
    if (asresid & !vector.x) 
        ASR <- (CST$observed - CST$expected)/sqrt(CST$expected * 
            ((1 - RS/GT) %*% t(1 - CS/GT)))
    print.CrossTable.SAS <- function() {
        if (exists("RowData")) {
            cat(SpaceSep1, "|", ColData, "\n")
            cat(formatC(RowData, width = RWidth, format = "s"), 
                formatC(dimnames(t)[[2]], width = CWidth, format = "s"), 
                RowTotal, sep = " | ", collapse = "\n")
        }
        else cat(SpaceSep1, formatC(dimnames(t)[[2]], width = CWidth, 
            format = "s"), RowTotal, sep = " | ", collapse = "\n")
        cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", collapse = "\n")
        for (i in 1:nrow(t)) {
            cat(FirstCol[i], formatC(c(t[i, ], RS[i]), width = CWidth, 
                format = "d"), sep = " | ", collapse = "\n")
            if (expected) 
                cat(SpaceSep1, formatC(CST$expected[i, ], digits = digits, 
                  format = "f", width = CWidth), SpaceSep2, sep = " | ", 
                  collapse = "\n")
            if (prop.chisq) 
                cat(SpaceSep1, formatC((((CST$expected[i, ] - 
                  t[i, ])^2)/CST$expected[i, ]), width = CWidth, 
                  digits = digits, format = "f"), SpaceSep2, 
                  sep = " | ", collapse = "\n")
            if (prop.r) 
                cat(SpaceSep1, formatC(c(CPR[i, ], RS[i]/GT), 
                  width = CWidth, digits = digits, format = "f"), 
                  sep = " | ", collapse = "\n")
            if (prop.c) 
                cat(SpaceSep1, formatC(CPC[i, ], width = CWidth, 
                  digits = digits, format = "f"), SpaceSep2, 
                  sep = " | ", collapse = "\n")
            if (prop.t) 
                cat(SpaceSep1, formatC(CPT[i, ], width = CWidth, 
                  digits = digits, format = "f"), SpaceSep2, 
                  sep = " | ", collapse = "\n")
            cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", 
                collapse = "\n")
        }
        cat(ColTotal, formatC(c(CS, GT), width = CWidth, format = "d"), 
            sep = " | ", collapse = "\n")
        if (prop.c) 
            cat(SpaceSep1, formatC(CS/GT, width = CWidth, digits = digits, 
                format = "f"), SpaceSep2, sep = " | ", collapse = "\n")
        cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", collapse = "\n")
    }
    print.CrossTable.SPSS <- function() {
        if (exists("RowData")) {
            cat(SpaceSep1, "|", ColData, "\n")
            cat(cat(formatC(RowData, width = RWidth, format = "s"), 
                sep = " | ", collapse = ""), cat(formatC(dimnames(t)[[2]], 
                width = CWidth - 1, format = "s"), sep = "  | ", 
                collapse = ""), cat(RowTotal, sep = " | ", collapse = "\n"), 
                sep = "", collapse = "")
        }
        else cat(SpaceSep1, formatC(dimnames(t)[[2]], width = CWidth, 
            format = "s"), RowTotal, sep = " | ", collapse = "\n")
        cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", collapse = "\n")
        for (i in 1:nrow(t)) {
            cat(cat(FirstCol[i], sep = " | ", collapse = ""), 
                cat(formatC(c(t[i, ], RS[i]), width = CWidth - 
                  1, format = "d"), sep = "  | ", collapse = "\n"), 
                sep = "", collapse = "")
            if (expected) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC(CST$expected[i, ], digits = digits, 
                    format = "f", width = CWidth - 1), sep = "  | ", 
                    collapse = ""), cat(SpaceSep2, sep = " | ", 
                    collapse = "\n"), sep = "", collapse = "")
            if (prop.chisq) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC((((CST$expected[i, ] - t[i, ])^2)/CST$expected[i, 
                    ]), digits = digits, format = "f", width = CWidth - 
                    1), sep = "  | ", collapse = ""), cat(SpaceSep2, 
                    sep = " | ", collapse = "\n"), sep = "", 
                  collapse = "")
            if (prop.r) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC(c(CPR[i, ] * 100, 100 * RS[i]/GT), 
                    width = CWidth - 1, digits = digits, format = "f"), 
                    sep = "% | ", collapse = "\n"), sep = "", 
                  collapse = "")
            if (prop.c) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC(CPC[i, ] * 100, width = CWidth - 
                    1, digits = digits, format = "f"), sep = "% | ", 
                    collapse = ""), cat(SpaceSep2, sep = " | ", 
                    collapse = "\n"), sep = "", collapse = "")
            if (prop.t) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC(CPT[i, ] * 100, width = CWidth - 
                    1, digits = digits, format = "f"), sep = "% | ", 
                    collapse = ""), cat(SpaceSep2, sep = " | ", 
                    collapse = "\n"), sep = "", collapse = "")
            if (resid) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC(CST$observed[i, ] - CST$expected[i, 
                    ], digits = digits, format = "f", width = CWidth - 
                    1), sep = "  | ", collapse = ""), cat(SpaceSep2, 
                    sep = " | ", collapse = "\n"), sep = "", 
                  collapse = "")
            if (sresid) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC(CST$residual[i, ], digits = digits, 
                    format = "f", width = CWidth - 1), sep = "  | ", 
                    collapse = ""), cat(SpaceSep2, sep = " | ", 
                    collapse = "\n"), sep = "", collapse = "")
            if (asresid) 
                cat(cat(SpaceSep1, sep = " | ", collapse = ""), 
                  cat(formatC(ASR[i, ], digits = digits, format = "f", 
                    width = CWidth - 1), sep = "  | ", collapse = ""), 
                  cat(SpaceSep2, sep = " | ", collapse = "\n"), 
                  sep = "", collapse = "")
            cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", 
                collapse = "\n")
        }
        cat(cat(ColTotal, sep = " | ", collapse = ""), cat(formatC(c(CS, 
            GT), width = CWidth - 1, format = "d"), sep = "  | ", 
            collapse = "\n"), sep = "", collapse = "")
        if (prop.c) 
            cat(cat(SpaceSep1, sep = " | ", collapse = ""), cat(formatC(100 * 
                CS/GT, width = CWidth - 1, digits = digits, format = "f"), 
                sep = "% | ", collapse = ""), cat(SpaceSep2, 
                sep = " | ", collapse = "\n"), sep = "", collapes = "")
        cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", collapse = "\n")
    }
    print.CrossTable.vector.SAS <- function() {
        if (length(t) > max.width) {
            final.row <- length(t)%%max.width
            max <- length(t) - final.row
            start <- seq(1, max, max.width)
            end <- start + (max.width - 1)
            if (final.row > 0) {
                start <- c(start, end[length(end)] + 1)
                end <- c(end, end[length(end)] + final.row)
            }
        }
        else {
            start <- 1
            end <- length(t)
        }
        SpaceSep3 <- paste(SpaceSep2, " ", sep = "")
        for (i in 1:length(start)) {
            cat(SpaceSep2, formatC(dimnames(t)[[2]][start[i]:end[i]], 
                width = CWidth, format = "s"), sep = " | ", collapse = "\n")
            cat(SpaceSep3, rep(RowSep, (end[i] - start[i]) + 
                1), sep = "|", collapse = "\n")
            cat(SpaceSep2, formatC(t[, start[i]:end[i]], width = CWidth, 
                format = "d"), sep = " | ", collapse = "\n")
            cat(SpaceSep2, formatC(CPT[, start[i]:end[i]], width = CWidth, 
                digits = digits, format = "f"), sep = " | ", 
                collapse = "\n")
            cat(SpaceSep3, rep(RowSep, (end[i] - start[i]) + 
                1), sep = "|", collapse = "\n")
            cat("\n\n")
        }
    }
    print.CrossTable.vector.SPSS <- function() {
        if (length(t) > max.width) {
            final.row <- length(t)%%max.width
            max <- length(t) - final.row
            start <- seq(1, max, max.width)
            end <- start + (max.width - 1)
            if (final.row > 0) {
                start <- c(start, end[length(end)] + 1)
                end <- c(end, end[length(end)] + final.row)
            }
        }
        else {
            start <- 1
            end <- length(t)
        }
        SpaceSep3 <- paste(SpaceSep2, " ", sep = "")
        for (i in 1:length(start)) {
            cat(cat(SpaceSep2, sep = " | ", collapse = ""), cat(formatC(dimnames(t)[[2]][start[i]:end[i]], 
                width = CWidth - 1, format = "s"), sep = "  | ", 
                collapse = "\n"), sep = "", collapse = "")
            cat(SpaceSep3, rep(RowSep, (end[i] - start[i]) + 
                1), sep = "|", collapse = "\n")
            cat(cat(SpaceSep2, sep = " | ", collapse = ""), cat(formatC(t[, 
                start[i]:end[i]], width = CWidth - 1, format = "d"), 
                sep = "  | ", collapse = "\n"), sep = "", collapse = "")
            cat(cat(SpaceSep2, sep = " | ", collapse = ""), cat(formatC(CPT[, 
                start[i]:end[i]] * 100, width = CWidth - 1, digits = digits, 
                format = "f"), sep = "% | ", collapse = ""), 
                sep = "", collapse = "\n")
            cat(SpaceSep3, rep(RowSep, (end[i] - start[i]) + 
                1), sep = "|", collapse = "\n")
        }
        if (GT < TotalN) 
            cat("\nNumber of Missing Observations: ", TotalN - 
                GT, " (", 100 * (TotalN - GT)/TotalN, "%)\n", 
                sep = "")
    }
    print.statistics <- function() {
        if (chisq) {
            cat(rep("\n", 2))
            cat("Statistics for All Table Factors\n\n\n")
            cat(CST$method, "\n")
            cat("------------------------------------------------------------\n")
            cat("Chi^2 = ", CST$statistic, "    d.f. = ", CST$parameter, 
                "    p = ", CST$p.value, "\n\n")
            if (all(dim(t) == 2)) {
                cat(CSTc$method, "\n")
                cat("------------------------------------------------------------\n")
                cat("Chi^2 = ", CSTc$statistic, "    d.f. = ", 
                  CSTc$parameter, "    p = ", CSTc$p.value, "\n")
            }
        }
        if (mcnemar) {
            McN <- mcnemar.test(t, correct = FALSE)
            cat(rep("\n", 2))
            cat(McN$method, "\n")
            cat("------------------------------------------------------------\n")
            cat("Chi^2 = ", McN$statistic, "    d.f. = ", McN$parameter, 
                "    p = ", McN$p.value, "\n\n")
            if (all(dim(t) == 2)) {
                McNc <- mcnemar.test(t, correct = TRUE)
                cat(McNc$method, "\n")
                cat("------------------------------------------------------------\n")
                cat("Chi^2 = ", McNc$statistic, "    d.f. = ", 
                  McNc$parameter, "    p = ", McNc$p.value, "\n")
            }
        }
        if (fisher) {
            cat(rep("\n", 2))
            FTt <- fisher.test(t, alternative = "two.sided")
            if (all(dim(t) == 2)) {
                FTl <- fisher.test(t, alternative = "less")
                FTg <- fisher.test(t, alternative = "greater")
            }
            cat("Fisher's Exact Test for Count Data\n")
            cat("------------------------------------------------------------\n")
            if (all(dim(t) == 2)) {
                cat("Sample estimate odds ratio: ", FTt$estimate, 
                  "\n\n")
                cat("Alternative hypothesis: true odds ratio is not equal to 1\n")
                cat("p = ", FTt$p.value, "\n")
                cat("95% confidence interval: ", FTt$conf.int, 
                  "\n\n")
                cat("Alternative hypothesis: true odds ratio is less than 1\n")
                cat("p = ", FTl$p.value, "\n")
                cat("95% confidence interval: ", FTl$conf.int, 
                  "\n\n")
                cat("Alternative hypothesis: true odds ratio is greater than 1\n")
                cat("p = ", FTg$p.value, "\n")
                cat("95% confidence interval: ", FTg$conf.int, 
                  "\n\n")
            }
            else {
                cat("Alternative hypothesis: two.sided\n")
                cat("p = ", FTt$p.value, "\n")
            }
        }
        cat(rep("\n", 2))
        CT <- list(t = t, prop.row = CPR, prop.col = CPC, prop.tbl = CPT)
        if (any(chisq, fisher, mcnemar)) {
            if (all(dim(t) == 2)) {
                if (chisq) 
                  CT <- c(CT, list(chisq = CST, chisq.corr = CSTc))
                if (fisher) 
                  CT <- c(CT, list(fisher.ts = FTt, fisher.tl = FTl, 
                    fisher.gt = FTg))
                if (mcnemar) 
                  CT <- c(CT, list(mcnemar = McN, mcnemar.corr = McNc))
            }
            else {
                if (chisq) 
                  CT <- c(CT, list(chisq = CST))
                if (fisher) 
                  CT <- c(CT, list(fisher.ts = FTt))
                if (mcnemar) 
                  CT <- c(CT, list(mcnemar = McN))
            }
        }
        invisible(CT)
    }
    if (format == "SAS") {
        cat(rep("\n", 2))
        cat("   Cell Contents\n")
        cat("|-------------------------|\n")
        cat("|                       N |\n")
        if (expected) 
            cat("|              Expected N |\n")
        if (prop.chisq) 
            cat("| Chi-square contribution |\n")
        if (prop.r) 
            cat("|           N / Row Total |\n")
        if (prop.c) 
            cat("|           N / Col Total |\n")
        if (prop.t) 
            cat("|         N / Table Total |\n")
        cat("|-------------------------|\n")
        cat(rep("\n", 2))
        cat("Total Observations in Table: ", GT, "\n")
        cat(rep("\n", 2))
        if (!vector.x) 
            print.CrossTable.SAS()
        else print.CrossTable.vector.SAS()
        print.statistics()
    }
    else if (format == "SPSS") {
        cat("\n")
        cat("   Cell Contents\n")
        cat("|-------------------------|\n")
        cat("|                   Count |\n")
        if (!vector.x) {
            if (expected) 
                cat("|         Expected Values |\n")
            if (prop.chisq) 
                cat("| Chi-square contribution |\n")
            if (prop.r) 
                cat("|             Row Percent |\n")
            if (prop.c) 
                cat("|          Column Percent |\n")
            if (prop.t) 
                cat("|           Total Percent |\n")
            if (resid) 
                cat("|                Residual |\n")
            if (sresid) 
                cat("|            Std Residual |\n")
            if (asresid) 
                cat("|           Adj Std Resid |\n")
        }
        else cat("|             Row Percent |\n")
        cat("|-------------------------|\n")
        cat("\n")
        cat("Total Observations in Table: ", GT, "\n")
        cat("\n")
        if (!vector.x) 
            print.CrossTable.SPSS()
        else print.CrossTable.vector.SPSS()
        CTres <- print.statistics()
        if (any(dim(t) >= 2) & any(chisq, mcnemar, fisher)) {
            MinExpF = min(CST$expected)
            cat("       Minimum expected frequency:", MinExpF, 
                "\n")
            NMinExpF = length(CST$expected[which(CST$expected < 
                5)])
            if (NMinExpF > 0) {
                NCells = length(CST$expected)
                cat("Cells with Expected Frequency < 5: ", NMinExpF, 
                  " of ", NCells, " (", 100 * NMinExpF/NCells, 
                  "%)\n", sep = "")
            }
            cat("\n")
        }
        CTres
    }
    else stop("unknown format")
}
