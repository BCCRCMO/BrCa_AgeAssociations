###


if ( !file.exists("Plots") ) { dir.create("Plots") }

### BigSeries data to produce age histograms comparing age distributions
### between BigSeries and NormalBreast13 TMAs
wmadf <- read.table("c:/kilroy/Projects/TMA/02-008_BigSeries/Analysis/BuildData/02-008_BigSeries_Whole.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "")
names(wmadf)[1] <- "bseries_id"


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
addf <- merge(merge(pddf, eddf),
              hddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(hddf))],
              by = "nb13_id", suffixes = c(".ezh2", ".h3k27me3"))
head(addf)
dim(addf)
save.image()
### Add in next raw score

rddf <- read.table(file = "../../Scoring/FOXA1/NB13_FOXA1_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## "
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.foxa1", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.foxa1", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)

### NEW TARGET - Copy e.g. ER processing down to end of list below and modify appropriately.

### ER

rddf <- read.table(file = "../../Scoring/ESR1_ER/NB13_ESR1_ER_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## "
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.esr1", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.esr1", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)

### BCL2

rddf <- read.table(file = "../../Scoring/BCL2/NB13_BCL2_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## " quote to control ESS
head(rddf)
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.bcl2", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.bcl2", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)

### MKI67_Ki67

rddf <- read.table(file = "../../Scoring/MKI67_Ki67/NB13_MKI67_Ki67_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## " quote to control ESS
head(rddf)
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.mki67", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.mki67", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)

### PGR_PR

rddf <- read.table(file = "../../Scoring/PGR_PR/NB13_PGR_PR_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## " quote to control ESS
head(rddf)
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.pgr", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.pgr", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)

### TP53_p53

rddf <- read.table(file = "../../Scoring/TP53_p53/NB13_TP53_p53_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## " quote to control ESS
head(rddf)
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.tp53", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.tp53", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)

### KRT5_CK5

rddf <- read.table(file = "../../Scoring/KRT5_CK5/NB13_KRT5_CK5_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## " quote to control ESS
head(rddf)
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.ck5", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.ck5", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)


### CDH1_Ecad

rddf <- read.table(file = "../../Scoring/CDH1_Ecad/NB13_CDH1_Ecad_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## " quote to control ESS
head(rddf)
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.ecad", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.ecad", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)

### EGFR

rddf <- read.table(file = "../../Scoring/EGFR/NB13_EGFR_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## " quote to control ESS
head(rddf)
tail(rddf)
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.egfr", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.egfr", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)


### HER2
### c = cytoplasmic staining in addition to membranous staining.

rddf <- read.table(file = "../../Scoring/ERBB2_HER2/NB13_ERBB2_HER2_v01.csv", sep = ",", quote = '"',
                   stringsAsFactors = FALSE, header = TRUE) ## "
names(rddf) <- gsub("NB13ID", "nb13_id", names(rddf))
names(rddf) <- gsub("PctPosCells", "PctPosCells.her2", names(rddf))
names(rddf) <- gsub("Intensity", "Intensity.her2", names(rddf))

addfbackup <- addf
addf <- merge(addfbackup, rddf[, -match(c("SurgNo", "CoreNo", "NumCores"), names(rddf))], by = "nb13_id")
head(addf)
dim(addf)



save.image()



### Select out the randomly selected side for patients with bilateral samples.
sddf <- addf[addf$isUnpaired == 1, ]
sddf$age_at_diagnosis <- sddf$Age

table(sddf$PctPosCells.egfr)
table(sddf$Intensity.egfr)
table(sddf$PctPosCells.ezh2)
table(sddf$PctPosCells.h3k27me3)
table(sddf$PctPosCells.foxa1)
table(sddf$PctPosCells.esr1)
table(sddf$PctPosCells.bcl2)
table(sddf$PctPosCells.mki67)
table(sddf$PctPosCells.pgr)
table(sddf$PctPosCells.tp53)
table(sddf$PctPosCells.ck5)
table(sddf$PctPosCells.ecad)
table(sddf$PctPosCells.her2)

### s denotes small amount of tissue, with 50-200 evaluable cells.
### No gsub() on "s" will have the effect of setting such values to NA.
### Tomo Osako asked to see results with "s" cores removed.


sddf$int.egfr <- as.numeric(sddf$Intensity.egfr)
with(sddf, table(Intensity.egfr, int.egfr, useNA = "always"))

###sddf$pp.ezh2 <- as.numeric(gsub("s","",sddf$PctPosCells.ezh2))
sddf$pp.ezh2 <- as.numeric(sddf$PctPosCells.ezh2)
with(sddf, table(PctPosCells.ezh2, pp.ezh2, useNA = "always"))
###sddf$pp.h3k27me3 <- as.numeric(gsub("s","",sddf$PctPosCells.h3k27me3))
sddf$pp.h3k27me3 <- as.numeric(sddf$PctPosCells.h3k27me3)
with(sddf, table(PctPosCells.h3k27me3, pp.h3k27me3, useNA = "always"))

sddf$pp.foxa1 <- as.numeric(sddf$PctPosCells.foxa1)
with(sddf, table(PctPosCells.foxa1, pp.foxa1, useNA = "always"))

sddf$pp.esr1 <- as.numeric(sddf$PctPosCells.esr1)
with(sddf, table(PctPosCells.esr1, pp.esr1, useNA = "always"))

sddf$pp.bcl2 <- as.numeric(sddf$PctPosCells.bcl2)
with(sddf, table(PctPosCells.bcl2, pp.bcl2, useNA = "always"))

sddf$pp.mki67 <- as.numeric(gsub("<1","0.5",sddf$PctPosCells.mki67))
with(sddf, table(PctPosCells.mki67, pp.mki67, useNA = "always"))

sddf$pp.pgr <- as.numeric(gsub("<1","0.5",sddf$PctPosCells.pgr))
with(sddf, table(PctPosCells.pgr, pp.pgr, useNA = "always"))

sddf$pp.tp53 <- as.numeric(gsub("<1","0.5",sddf$PctPosCells.tp53))
with(sddf, table(PctPosCells.tp53, pp.tp53, useNA = "always"))

sddf$pp.ck5 <- as.numeric(gsub("<1","0.5",sddf$PctPosCells.ck5))
with(sddf, table(PctPosCells.ck5, pp.ck5, useNA = "always"))

sddf$pp.ecad <- as.numeric(sddf$PctPosCells.ecad)
with(sddf, table(PctPosCells.ecad, pp.ecad, useNA = "always"))

sddf$pp.her2 <- as.numeric(gsub("c","",sddf$PctPosCells.her2))
##sddf$pp.her2 <- as.numeric(sddf$PctPosCells.her2)
with(sddf, table(PctPosCells.her2, pp.her2, useNA = "always"))



### No "s" entries:
### > median(sddf$pp.h3k27me3, na.rm = TRUE)
### [1] 85  # Same as when including "s" scores.
### 
### > median(sddf$pp.ezh2, na.rm = TRUE)
### [1] 3


plot( sddf$Age, sddf$pp.h3k27me3)
lines(supsmu( sddf$Age, sddf$pp.h3k27me3, bass = 3))

plot( sddf$Age, sddf$pp.h3k27me3)
lines(supsmu( sddf$Age, sddf$pp.h3k27me3, bass = 10))

plot( sddf$Age, sddf$pp.ezh2)
lines(supsmu( sddf$Age, sddf$pp.h3k27me3, bass = 10))


sddf$h3k27me3_bgtm_v0n <- (1.0 * (sddf$pp.h3k27me3 > median(sddf$pp.h3k27me3, na.rm = TRUE)))
sddf$h3k27me3_bgtm_v0f <- factor(as.character((1.0 * (sddf$pp.h3k27me3 > median(sddf$pp.h3k27me3, na.rm = TRUE)))),
                                 levels = c("0", "1"),
                                 labels = c("H3K27me3 <= median", "H3K27me3 > median")
                                 )
set.seed(3117); sddf$jpp.h3k27me3 <- jitter(1.0 * (sddf$pp.h3k27me3 > median(sddf$pp.h3k27me3, na.rm = TRUE)), factor = 0.1)

sddf$ezh2_bgtm_v1n <- (1.0 * (sddf$pp.ezh2 > median(sddf$pp.ezh2, na.rm = TRUE)))
sddf$ezh2_bgtm_v1f <- factor(as.character((1.0 * (sddf$pp.ezh2 > median(sddf$pp.ezh2, na.rm = TRUE)))),
                                 levels = c("0", "1"),
                                 labels = c("EZH2 <= median", "EZH2 > median")
                                 )
set.seed(7743); sddf$jpp.ezh2 <- jitter(1.0 * (sddf$pp.ezh2 > median(sddf$pp.ezh2, na.rm = TRUE)), factor = 0.1)

sddf$foxa1_bgtm_v1n <- (1.0 * (sddf$pp.foxa1 > median(sddf$pp.foxa1, na.rm = TRUE)))
sddf$foxa1_bgtm_v1f <- factor(as.character((1.0 * (sddf$pp.foxa1 > median(sddf$pp.foxa1, na.rm = TRUE)))),
                                 levels = c("0", "1"),
                                 labels = c("FOXA1 <= median", "FOXA1 > median")
                                 )
set.seed(7745); sddf$jpp.foxa1 <- jitter(1.0 * (sddf$pp.foxa1 > median(sddf$pp.foxa1, na.rm = TRUE)), factor = 0.1)

sddf$esr1_bgtm_v1n <- (1.0 * (sddf$pp.esr1 > median(sddf$pp.esr1, na.rm = TRUE)))
sddf$esr1_bgtm_v1f <- factor(as.character((1.0 * (sddf$pp.esr1 > median(sddf$pp.esr1, na.rm = TRUE)))),
                                 levels = c("0", "1"),
                                 labels = c("ESR1 <= median", "ESR1 > median")
                                 )
set.seed(7747); sddf$jpp.esr1 <- jitter(1.0 * (sddf$pp.esr1 > median(sddf$pp.esr1, na.rm = TRUE)), factor = 0.1)

plot( sddf$Age, sddf$jpp.h3k27me3 )
lines(supsmu( sddf$Age, sddf$jpp.h3k27me3, bass = 10))

### Contingency table analysis:
fscvars <- c("ezh2_bgtm_v1f", "h3k27me3_bgtm_v0f", "foxa1_bgtm_v1f", "esr1_bgtm_v1f")
bscvars <- c("ezh2_bgtm_v1n", "h3k27me3_bgtm_v0n", "foxa1_bgtm_v1n", "esr1_bgtm_v1n")

cscvars <- c(
  "pp.bcl2",
  "pp.ck5",
  "pp.ecad",
  "int.egfr",
  "pp.esr1",
  "pp.ezh2",
  "pp.foxa1",
  "pp.h3k27me3",
  "pp.her2",
  "pp.mki67",
  "pp.pgr",
  "pp.tp53"
  )


cscvars4 <- c(
  "pp.esr1",
  "pp.ezh2",
  "pp.foxa1",
  "pp.h3k27me3"
  )



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
    feout <- if (min(ctout$chisq$expected) < 3) {
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

pdf("./Plots/BiomarkerRatesPiecewiseByAge_NB13_No_s_v03.pdf", width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 0, 0), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binregfit, binchisqlist = bsctableoutadj )
dev.off()




### Randomization p-value:
### If no trend, ages exchangeable.  Randomize ages.  Calculate regression slope.
### How often does observed slope exceed randomization slope?



### Age histograms - BigSeries, NormalBreast13

win.metafile(file = "BigSeries_NB13_Age_densityplot_v01.wmf", width = 6, height = 5)

plot(c(10, 90), c(0, .035), type = "n", xlab = "Patient Age (Years)", ylab = "Frequency",
     main = "TMA age distributions")
lines(density(wmadf$age_at_diagnosis), lty = 1, lwd = 3)
lines(density(sddf$Age, na.rm = TRUE), lty = 2, lwd = 5)
legend("topleft", legend = c("Big Series TMA", "Normal Breast TMA"), lty = c(1, 2),
       lwd = c(3, 5), seg.len = 10, text.width = 30)
dev.off()
dev.print(pdf, file = "BigSeries_NB13_Age_densityplot_v01.pdf", width = 6, height = 6)

### Continuous plots
simsmufit <- function(x, xdf = sddf, ydf = NULL, xvar = "age_at_diagnosis",
                      binchisqlist = NULL, ybinchisqlist = NULL,
                      tlbl = "", nsim = 2000, alpha = 0.05) {
  rn <- nsim
  alphalev <- alpha/2.0
  UCLn <- trunc(rn * (1 - alphalev))
  LCLn <- trunc(rn * alphalev)
  agesv <- sort(unique(xdf[, xvar]), na.last = NA)
  simmat <- matrix(NA_real_, nrow = rn, ncol = length(agesv))
  ylims <- range(xdf[, x], na.rm = TRUE)
  ylimd <- abs(diff(ylims))
  ylimlo <- ylims[1] - (0.2 * ylimd)
  ylimhi <- ylims[2] + (0.2 * ylimd)
  valdx <- !(is.na(xdf[, x]) | is.na(xdf[, xvar]))
###   rngy <- range(xdf[valdx, x])
###   if (rngy[2] < 5) {
###     ylims <- c(-1, 5)
###   } else {
###     ylims <- c(-5, 110)
###   }
###   plot( xdf[valdx, xvar], xdf[valdx, x], xlab = xvar, ylab = x, xlim = c(15, 95),
###        ylim = ylims, yaxt = "n", pch = ".", col = "grey60" ) 
  xlbc <- "Dx age (years)"
  ylbc <- toupper(strsplit(x, split = "\\.")[[1]][2])
  plot( xdf[valdx, xvar], xdf[valdx, x], xlab = xlbc, ylab = ylbc,
       xlim = c(15, 95), ylim = c(ylimlo, ylimhi), pch = ".", col = "grey60" )
###   if (rngy[2] < 5) {
###     axis(side = 2, at = c(0, 1, 2, 3)) ## Intensity
###   } else {
###     axis(side = 2, at = c(0, 20, 40, 60, 80, 100)) ## Pct cells stained
###   }
  for (i in seq(rn) ) {
    ssout <- supsmu( sample(xdf[valdx, xvar]), xdf[valdx, x], bass = 10)
    spout <- predict(smooth.spline(ssout$x, ssout$y), x = agesv)
    simmat[i, ] <- spout$y
###     lines(ssout, lty = 3, lwd = 1, col = "lightgrey")
  }
  ##points( xdf[, xvar], xdf[, x])
  srtmat <- apply(simmat, 2, sort)
###   lines(agesv, srtmat[LCLn, ], lty = 3, lwd = 3)
###   lines(agesv, srtmat[UCLn, ], lty = 3, lwd = 3)
  polygon(c(agesv, rev(agesv)), c(srtmat[LCLn, ], rev(srtmat[UCLn, ])), col = "wheat3", border = "wheat3")
###   lines(supsmu( xdf[valdx, xvar], xdf[valdx, x], bass = 10), lwd = 3, type = "l")
  ssout <- supsmu( xdf[valdx, xvar], xdf[valdx, x], bass = 10)
  lines(ssout, lwd = 2, type = "l")
  srtmatcols <- match(ssout$x, agesv)
  if ( length(ssout$x) <= length(agesv) ) {
    if ( any( ssout$y < srtmat[LCLn, srtmatcols] ) |
         any( ssout$y > srtmat[UCLn, srtmatcols] ) ) {
      idxs <- c( which( ssout$y < srtmat[LCLn, srtmatcols] ),
                 which( ssout$y > srtmat[UCLn, srtmatcols] ) )
      points(ssout$x[idxs], ssout$y[idxs], pch = 20)
      text(x = 20, y = ylimhi, labels = paste("Age-dependent p<", alpha, sep = ""), adj = 0)
    } else {
      text(x = 20, y = ylimhi, labels = paste("Not age-dependent p>", alpha, sep = ""), adj = 0)
    }
  } else {
    text(x = 20, y = ylimhi, labels = "FIX ME missing ages", adj = 0)
  }
###   legend("bottomright", legend = paste(trunc(100*(1-alpha)), "% Confidence Bounds", sep = ""),
###          lty = 3, lwd = 3, seg.len = 3, text.width = 60)
}

pdf(file = "NormalBreast13_IHC_Trends_99p9pct_40kSimCIs_shaded_v02.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, nsim = 40000, alpha = 0.001)
dev.off()

save.image()

pdf(file = "NormalBreast13_IHC_Trends_95pct_SimCIs_4Tgts_shaded_v02.pdf", width = 8, height = 8)
par(mfrow = c(2, 2), oma = c(1,1,3,1))
sapply(cscvars4, simsmufit)
mtext("NormalBreast13", outer=TRUE)
dev.off()

dev.print(pdf, file = "NormalBreast13_IHC_Trends_95pct_SimCIs_v01.pdf", width = 7, height = 7)
dev.print(win.metafile, file = "NormalBreast13_IHC_Trends_95pct_SimCIs_v03.wmf", width = 7, height = 7)

### Patient characteristics
quantile(sddf$Age, na.rm = TRUE)
sddf$Agef <- cut(sddf$Age, breaks = c(0, 30, 40, 50, 60, 70, 80, 100), right = FALSE)
table(sddf$Agef)
table(sddf$Agef)/sum(table(sddf$Agef))*100
