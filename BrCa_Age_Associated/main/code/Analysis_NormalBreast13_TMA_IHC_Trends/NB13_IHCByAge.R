###


if ( !file.exists("Plots") ) { dir.create("Plots") }

pddf <-
  read.table(file = "../../NormalBreast13Creation/RandomizedLat4coresBefore2cores_NormalBreastTMA_NB13.csv",
             sep = ",", stringsAsFactors = FALSE, header = TRUE)
eddf <- read.table(file = "../../Scoring/EZH2/NB13_EZH2_v02.tsv", sep = "\t",
                   stringsAsFactors = FALSE, header = TRUE)

hddf <- read.table(file = "../../Scoring/H3K27me3/NB13_H3K27me3_v02.tsv", sep = "\t",
                   stringsAsFactors = FALSE, header = TRUE)
dim(pddf)
dim(eddf)
dim(hddf)
names(pddf)
names(eddf)

intersect(names(pddf), names(eddf))
addf <- merge(merge(pddf, eddf), hddf, by = "nb13_id", suffixes = c(".e", ".h"))
head(addf)

sddf <- addf[addf$isUnpaired == 1, ]

table(sddf$PctPosCells.e)

sddf$pp.e <- as.numeric(gsub("s","",sddf$PctPosCells.e))
with(sddf, table(PctPosCells.e, pp.e, useNA = "always"))
sddf$pp.h <- as.numeric(gsub("s","",sddf$PctPosCells.h))
with(sddf, table(PctPosCells.h, pp.h, useNA = "always"))


plot( sddf$Age, sddf$pp.h)
lines(supsmu( sddf$Age, sddf$pp.h, bass = 3))

plot( sddf$Age, sddf$pp.h)
lines(supsmu( sddf$Age, sddf$pp.h, bass = 10))

plot( sddf$Age, sddf$pp.e)
lines(supsmu( sddf$Age, sddf$pp.h, bass = 10))


sddf$h3k27me3_bgtm_v0n <- (1.0 * (sddf$pp.h > median(sddf$pp.h, na.rm = TRUE)))
sddf$h3k27me3_bgtm_v0f <- factor(as.character((1.0 * (sddf$pp.h > median(sddf$pp.h, na.rm = TRUE)))),
                                 levels = c("0", "1"),
                                 labels = c("H3K27me3 <= median", "H3K27me3 > median")
                                 )
set.seed(3117); sddf$jpp.h <- jitter(1.0 * (sddf$pp.h > median(sddf$pp.h, na.rm = TRUE)), factor = 0.1)

sddf$ezh2_bgtm_v1n <- (1.0 * (sddf$pp.e > median(sddf$pp.e, na.rm = TRUE)))
sddf$ezh2_bgtm_v1f <- factor(as.character((1.0 * (sddf$pp.e > median(sddf$pp.e, na.rm = TRUE)))),
                                 levels = c("0", "1"),
                                 labels = c("EZH2 <= median", "EZH2 > median")
                                 )
set.seed(7743); sddf$jpp.e <- jitter(1.0 * (sddf$pp.e > median(sddf$pp.e, na.rm = TRUE)), factor = 0.1)

plot( sddf$Age, sddf$jpp.h )
lines(supsmu( sddf$Age, sddf$jpp.h, bass = 10))

### Contingency table analysis:
sddf$age_at_diagnosis <- sddf$Age
fscvars <- c("ezh2_bgtm_v1f", "h3k27me3_bgtm_v0f")
bscvars <- c("ezh2_bgtm_v1n", "h3k27me3_bgtm_v0n")

cscvars <- c("pp.e", "pp.h")

source("smCrossTable.R")
binchisq <- function(x,
                     xdf = sddf,
                     agecat = c(10, 40, 55, 70, 100) ) {
  ##browser()
  xdf$agecat <- cut(xdf$age_at_diagnosis, breaks=agecat,
                    labels = c(paste("<=", agecat[2]),
                      paste(agecat[2], "-", agecat[3]),
                      paste(agecat[3], "-", agecat[4]),
                      paste(">", agecat[4])
                      ), right = FALSE,
                    include.lowest = TRUE)
  idxp <- !is.na(xdf[, x])
  
  cat(paste("\n\nNB13 TMA set:", x, "\n"))
  ## print(chiout <- chisq.test(xdf[idxp, x], xdf[idxp, "agecat"]))
  ## mtext(text = paste("Proportion", x), side = 2, outer = TRUE)
  ## print(chiout$observed)
  ## print(chiout$expected)
  ctout <- try(smCrossTable(xdf[idxp, "agecat"], xdf[idxp, x], expected = TRUE, prop.r = TRUE,
                        prop.c = TRUE, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE,
                        format = "SPSS", digits = 1,
                        dnn = c("Age",
                          paste(toupper(strsplit(x, split = "_")[[1]][1]), "(", x, ")"))) )
  if ( class(ctout) == "try-error" ) {
    ctout <- NULL
    feout <- NULL
  } else {
    feout <- if (min(ctout$chisq$expected) < 2) {
      fisher.test(xdf[idxp, "agecat"], xdf[idxp, x])
    } else { NULL }
  }
  if ( !is.null(feout) ) { print(feout) }
  return(list(
    ## chiout = chiout,
              ctout = ctout,
              feout = feout))
  
}

bsctableout <- sapply(bscvars, binchisq, simplify = FALSE, USE.NAMES = TRUE)
wbscpvals <- c(sapply(bsctableout, function(x) x[[1]]$chisq$p.value))
wbscpvalsadj <- p.adjust(wbscpvals, method = "BH")
bsctableoutadj <-
  lapply(seq(length(bsctableout)),
         function(x) {
           bsctableout[[x]]$ctout$chisq$p.adjust <- wbscpvalsadj[x]
           bsctableout[[x]]$p.adjust <- wbscpvalsadj[x]
           bsctableout[[x]]
         })
names(bsctableoutadj) <- names(bsctableout)


### Regression

binregfit <- function(x, xdf = sddf, ydf = NULL,
                      binchisqlist = NULL, ybinchisqlist = NULL,
                      tlbl = "") {
     ##browser()
  idxp <- !is.na(xdf[, x])
  if ( length( unique(xdf[idxp, x]) ) == 2 ) {

    ## No missing data df
    nmdf <- data.frame(xdf[idxp, ])

    nmdf$Agec <- rep(NA_character_, nrow(nmdf))

    nmdf$Agec[nmdf$age_at_diagnosis <= 40] <- "Age <= 40"
    nmdf$Agec[nmdf$age_at_diagnosis >  40] <- "Age >  40"

    nmdf$Agef <- factor(nmdf$Agec, levels = c("Age <= 40", "Age >  40"))
###     nmdf$Agem40 <- nmdf$age_at_diagnosis - 40.0
###     nmdf$Agem60 <- nmdf$age_at_diagnosis - 60.0
###     clmf <- eval(parse(text = paste(
###               "lm(", x, " ~ age_at_diagnosis + Agem40 + Agem60, data = nmdf)" ) ) )
###     mlmf <- eval(parse(text = paste( "lm(", x, " ~ 1, data = nmdf)" ) ) )
###     clmf <- eval(parse(text = paste(
###               "lm(", x, " ~ AgeLE40 + Age4160 + AgeGT60, data = nmdf)" ) ) )
    almf <- eval(parse(text = paste(
              "lm(", x, " ~ Agef * age_at_diagnosis, data = nmdf)" ) ) )
    almr <- eval(parse(text = paste(
              "lm(", x, " ~ age_at_diagnosis, data = nmdf)" ) ) )
    plot(xdf[idxp, "age_at_diagnosis"], jitter(xdf[idxp, x], factor = 0.1),
         ylab = paste("Proportion", x), xlab = "Age (years)", xlim = xlims, ylim = ylims)
    lines(supsmu(xdf[idxp, "age_at_diagnosis"], xdf[idxp, x], bass = 10), lwd = 3, col = "grey")
    predictvals <- c(15:38, 42:58, 62:88)
###     points(predictvals, predict(almf, data.frame(age_at_diagnosis = predictvals,
###                                            Agef = factor(c(rep("Age <= 40", sum(predictvals <= 40)),
###                                              rep("Age 41-60", sum((predictvals > 40) & (predictvals <= 60))),
###                                              rep("Age >  60", sum(predictvals > 60))),
###                                              levels = c("Age <= 40", "Age 41-60", "Age >  60") ))),
###            pch = "*")
    points(predictvals,
           predict(almf, data.frame(age_at_diagnosis = predictvals,
                                    Agef = factor(c(rep("Age <= 40", sum(predictvals <= 40)),
                                      rep("Age >  40", sum(predictvals > 40))),
                                      levels = c("Age <= 40", "Age >  40") ))),
           pch = "*")
###     points(predictvals, predict(clmf, data.frame(AgeLE40 = predictvals * (1.0 * (predictvals <= 40)),
###                                                  Age4160 = predictvals * (1.0 * ((predictvals > 40) &
###                                                    (predictvals <= 60) ) ),
###                                                  AgeGT60 = predictvals * (1.0 * (predictvals > 60))  ) ) )
    lines(c(10, 80), predict(almr, data.frame(age_at_diagnosis = c(10, 80))))
    aovalmf <- anova(almr, almf)
###     aovmlmf <- anova(mlmf, almf)
    
    ## abline(almr, xlim = c(20, 80), omd = c(0.2, 0.8, 0.2, 0.8))
    
    abline(h = 0.5, v = 50, lty = 2)
    xdfttl <- "NB13 TMA set"
    if (tlbl != "") xdfttl <- paste(xdfttl, "(", tlbl, ")" )
    title(main = xdfttl)
    if ( !is.null(binchisqlist) && !is.null(binchisqlist[[x]]) ) {
      ypos <- ifelse(binchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
      xpos <- 75
      text(xpos, ypos, paste("Single slope p =",
                               format(aovalmf$"Pr(>F)"[2], digits = 3)))
###       if ( !is.null( binchisqlist[[x]]$p.adjust ) ) {
        ypos <- ypos - 0.06
        xpos <- xpos + 7
        text(xpos, ypos, paste("Zero slope p =",
                               format(summary(almr)$coefficients[2, 4], digits = 3)))
###       }
    }
  } else {
    plot(0, 0, type = "n", axes = FALSE, frame.plot = TRUE, xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = "", ylab = "")
    text(-0.2, 0, labels = paste("Data only available for level ",
                    unique(xdf[idxp, x]), "\nfor variable", x ) )    
  }
  if ( !is.null(ydf) ) {
    idxp <- !is.na(ydf[, x])
    if ( length( unique(ydf[idxp, x]) ) == 2 ) {
      plot(ydf[idxp, "age_at_diagnosis"], jitter(ydf[idxp, x], factor = 0.1),
           ylab = paste("Proportion", x), xlab = "Age (years)", xlim = xlims, ylim = ylims)
      lines(supsmu(ydf[idxp, "age_at_diagnosis"], ydf[idxp, x], bass = 10), lwd = 3)
      abline(h = 0.5, v = 50, lty = 2)
      ydfttl <- "No such for MB09 Validation set"
      if (tlbl != "") ydfttl <- paste(ydfttl, "(", tlbl, ")" )
      title(main = ydfttl)
      if ( !is.null(ybinchisqlist) && !is.null(ybinchisqlist[[x]]) ) {
        ypos <- ifelse(ybinchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
        xpos <- 75
        text(xpos, ypos, paste(toupper(strsplit(x, split = "_")[[1]][1]),
                               "by AgeCat: p =",
                               format(ybinchisqlist[[x]]$ctout$chisq$p.value, digits = 2)))
        if ( !is.null( ybinchisqlist[[x]]$p.adjust ) ) {
          ypos <- ypos - 0.06
          xpos <- xpos + 7
          text(xpos, ypos, paste("p.adjust =",
                                 format(ybinchisqlist[[x]]$p.adjust, digits = 2)))
        }
      }
    } else {
      plot(0, 0, type = "n", axes = FALSE, frame.plot = TRUE, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "", ylab = "")
      text(-0.2, 0, labels = paste("Data only available for level ",
                      unique(xdf[idxp, x]), "\nfor variable", x ) )    
    }
  }
  ## mtext(text = paste("Proportion", x), side = 2, outer = TRUE) 
}

xlims <- c(10, 100)
ylims <- c(-0.05, 1.05)

#binregfit(bscvars[2], binchisqlist = bsctableout)

pdf("./Plots/BiomarkerRatesPiecewiseByAge_NB13_v01.pdf", width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 0, 0), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binregfit, binchisqlist = bsctableoutadj )
dev.off()




### Randomization p-value:
### If no trend, ages exchangeable.  Randomize ages.  Calculate regression slope.
### How often does observed slope exceed randomization slope?



