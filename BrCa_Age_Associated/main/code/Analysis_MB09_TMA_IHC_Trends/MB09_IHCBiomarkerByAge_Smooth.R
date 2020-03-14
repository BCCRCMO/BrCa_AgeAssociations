### er_bgt0_v1, ki67_bgt14_v1, [H3K27me3 not yet available 20130621] Smooth fit against age.
### bcl2 ccnd1 ecad foxa1 h3k27me3

require("splines")
require("party")

if ( !file.exists("Plots") ) { dir.create("Plots") }

tmadf <- read.table("../BuildData/MB09_wholedata.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "")

tmadf$BrCa4f <- factor(tmadf$BrCa4,
                       levels = c("Luminalp", "Luminal/HER2+", "HER2+/ER-PR-", "TNP", "Unassigned"))

## BCL2_v1.pp
tmadf$bcl2_pp_v1n <- as.numeric(tmadf$BCL2_v1.pp)

## birc5_v1
tmadf$birc5_pp_v1n <- as.numeric(tmadf$BIRC5_c_v1)
tmadf$birc5_bgt3_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$birc5_bgt3_v1n[as.numeric(tmadf$BIRC5_c_v1) <= 3] <- 0
tmadf$birc5_bgt3_v1n[as.numeric(tmadf$BIRC5_c_v1) > 3] <- 1

## cd163_v1
tmadf$cd163_c_v1n <- as.numeric(tmadf$cd163_c_v1)

## ck56_v1
tmadf$ck56_c2v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ck56_c2v1n[tmadf$ck56_c2v1 == "Negative"] <- 0
tmadf$ck56_c2v1n[tmadf$ck56_c2v1 == "Any weak/moderate staining"] <- 1
tmadf$ck56_c2v1n[tmadf$ck56_c2v1 == "Strong staining in >20% of tumour cells"] <- 2


## ck5_v1
tmadf$ck5_c2v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ck5_c2v1n[tmadf$ck5_c2v1 == "Negative"] <- 0
tmadf$ck5_c2v1n[tmadf$ck5_c2v1 == "Any weak/moderate staining"] <- 1
tmadf$ck5_c2v1n[tmadf$ck5_c2v1 == "Strong staining in >20% of tumour cells"] <- 2


## ckit_v1
ckit_c_v1
tmadf$ckit_c_v1n <- as.numeric(tmadf$ckit_c_v1)

## cryab_v1
tmadf$cryab_c_v1n <- as.numeric(tmadf$cryab_c_v1)

## cyclind1 CCND1
tmadf$cyclind1_pp_v1n <- as.numeric(tmadf$CyclinD1_v1.pp)

## ecad_v1
tmadf$ecad_pp_v1n <- as.numeric(tmadf$ecad_v1.pp)

## egfr_v1
tmadf$egfr_tsipp_v1n <- as.numeric(tmadf$egfr_v1.c1)

tmadf$egfr_bgt2_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$egfr_bgt2_v1n[tmadf$egfr_b0v12_v1 == "Negative staining"] <- 0
tmadf$egfr_bgt2_v1n[tmadf$egfr_b0v12_v1 == "Any positive staining"] <- 1

## er_v1
tmadf$er_pp_v1n <- as.numeric(tmadf$er_v1.c1)

tmadf$er_bgt0_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$er_bgt0_v1n[tmadf$er_bgt0_v1 == "Negative staining"] <- 0
tmadf$er_bgt0_v1n[tmadf$er_bgt0_v1 == "Any positive staining"] <- 1

## er by dextran coated charcoal er_dcc_v1
tmadf$er_dcc_v1n <- as.numeric(substring(tmadf$er_result, 2, 4))

## ezh2_v1
tmadf$ezh2_pp_v1n <- as.numeric(tmadf$EZH2_v1.c1)
tmadf$EZH2_ble05vgt05_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$EZH2_ble05vgt05_v1n[tmadf$EZH2_ble05vgt05_v1.pp == "0 to 5%"] <- 0
tmadf$EZH2_ble05vgt05_v1n[tmadf$EZH2_ble05vgt05_v1.pp == "10 to 100%"] <- 1

## fascin_v1
tmadf$fascin_c_v1n <- as.numeric(tmadf$fascin_c_v1)


## foxa1_v1
tmadf$foxa1_pp_v1n <- as.numeric(tmadf$foxa1_v1.pp)

## h3k27me3_v1
tmadf$h3k27me3_pp_v1n <- as.numeric(tmadf$h3k27me3_v1.pp)

## her2_v1
tmadf$her2_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$her2_v1n[tmadf$her2_v1 == "Negative: no staining"] <- 0
tmadf$her2_v1n[tmadf$her2_v1 == "Negative: weak staining"] <- 0
tmadf$her2_v1n[tmadf$her2_v1 == "Equivocal: moderate staining"] <- 1
tmadf$her2_v1n[tmadf$her2_v1 == "Equivocal: 3+ staining in <30% cells"] <- 1
tmadf$her2_v1n[tmadf$her2_v1 == "Positive: Strong positive staining in 30% or more of cells"] <- 2


## inpp4b_v1

tmadf$inpp4b_c_v1.ppn <- as.numeric(tmadf$inpp4b_c_v1.pp)

## ki67_v1    ki67_c1v3_v1.6n ki67_bgt14_v1n  ki67_c1v3_v1.6 ki67_bgt14_v1
tmadf$ki67_pp_v1n <- as.numeric(tmadf$ki67_v1.c1)

tmadf$ki67_bgt14_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ki67_bgt14_v1n[tmadf$ki67_bgt14_v1 == "Negative: 0-14%"] <- 0
tmadf$ki67_bgt14_v1n[tmadf$ki67_bgt14_v1 == "Positive: 15% or higher"] <- 1

## nestin_v1  
tmadf$nestin_tsipp_v1n <- as.numeric(tmadf$nestin_v1.c1)

tmadf$nestin_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$nestin_b0v12_v1n[tmadf$nestin_b0v12_v1 == "Negative staining"] <- 0
tmadf$nestin_b0v12_v1n[tmadf$nestin_b0v12_v1 == "Any positive staining"] <- 1

## p16_v1  
tmadf$p16_qsipp_v1n <- as.numeric(tmadf$p16_c4v1.c1)

tmadf$p16_c2v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$p16_c2v1n[tmadf$p16_c2v1 == "0%"] <- 0
tmadf$p16_c2v1n[tmadf$p16_c2v1 == "1-50%"] <- 0
tmadf$p16_c2v1n[tmadf$p16_c2v1 == ">50%"] <- 1

## pr_v1
tmadf$pr_c_v1n <- as.numeric(tmadf$pr_c_v1)

## tp53_v1 
tmadf$tp53_pp_v1n <- as.numeric(tmadf$TP53_c_v1)

tmadf$TP53_bgt0_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$TP53_bgt0_v1n[tmadf$TP53_bgt0_v1 == "Negative staining"] <- 0
tmadf$TP53_bgt0_v1n[tmadf$TP53_bgt0_v1 == "Any positive staining"] <- 1

### Segmented by ER+/HER2- ER+/HER2+ ER-/HER2+ ER-/HER2-   Her2.Expr her2_v1n
tmadf$her2_v2n <- rep(1, nrow(tmadf))
tmadf$her2_v2n[tmadf$her2_v1n == 2] <- 2
tmadf$her2_v2n[tmadf$her2_v1n == 0] <- 0

tmadf$BrCaEH <- rep(NA_character_, nrow(tmadf))
tmadf[tmadf$erposneg == "Positive" & tmadf$her2_v2n == 2, "BrCaEH"] <- "ER+/HER2+"
tmadf[tmadf$erposneg == "Positive" & tmadf$her2_v2n == 0, "BrCaEH"] <- "ER+/HER2-"
tmadf[tmadf$erposneg == "Negative" & tmadf$her2_v2n == 2, "BrCaEH"] <- "ER-/HER2+"
tmadf[tmadf$erposneg == "Negative" & tmadf$her2_v2n == 0, "BrCaEH"] <- "ER-/HER2-"
tmadf$BrCaEHf <- factor(tmadf$BrCaEH, levels =  c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))

write.csv(tmadf, file = "TableS4_METABRIC_TMA_MB09_raw_data.csv")



### Look at age relationship with other biomarkers.

## Binary scatterplot smooths



bscvars <- c(
  "birc5_bgt3_v1n",
  "egfr_bgt2_v1n",
  "er_bgt0_v1n",
  "EZH2_ble05vgt05_v1n",
  "ki67_bgt14_v1n",      ##ki67_bgt14_v1  ##ki67
  "nestin_b0v12_v1n",
  "p16_c2v1n",
  "TP53_bgt0_v1n"  
  )

source("smCrossTable.R")
binchisq <- function(x,
                     xdf = tmadf,
                     agecat = c(20, 40, 55, 70, 100) ) {
  ##browser()
  xdf$agecat <- cut(xdf$age_at_diagnosis, breaks=agecat,
                    labels = c(paste("<=", agecat[2]),
                      paste(agecat[2], "-", agecat[3]),
                      paste(agecat[3], "-", agecat[4]),
                      paste(">", agecat[4])
                      ), right = FALSE,
                    include.lowest = TRUE)
  idxp <- !is.na(xdf[, x])
  
  cat(paste("\n\nMETABRIC TMA set:", x, "\n"))
  ## print(chiout <- chisq.test(xdf[idxp, x], xdf[idxp, "agecat"]))
  ## mtext(text = paste("Proportion", x), side = 2, outer = TRUE)
  ## print(chiout$observed)
  ## print(chiout$expected)
  ctout <- smCrossTable(xdf[idxp, "agecat"], xdf[idxp, x], expected = TRUE, prop.r = TRUE,
                        prop.c = TRUE, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE,
                        format = "SPSS", digits = 1,
                        dnn = c("Age",
                          paste(toupper(strsplit(x, split = "_")[[1]][1]), "(", x, ")")))
  feout <- if (min(ctout$chisq$expected) < 2) {
    fisher.test(xdf[idxp, "agecat"], xdf[idxp, x])
  } else { NULL }
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
           bsctableout[[x]]
         })
names(bsctableoutadj) <- names(bsctableout)




## Scatterplots with smooths
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

binsupsmu <- function(x, xdf = tmadf, ydf = NULL,
                      binchisqlist = NULL, ybinchisqlist = NULL,
                      tlbl = "") {

  idxp <- !is.na(xdf[, x])
  if ( length( unique(xdf[idxp, x]) ) == 2 ) {
    plot(xdf[idxp, "age_at_diagnosis"], jitter(xdf[idxp, x], factor = 0.1),
         ylab = paste("Proportion", x), xlab = "Age (years)", xlim = xlims, ylim = ylims)
    lines(supsmu(xdf[idxp, "age_at_diagnosis"], xdf[idxp, x], bass = 3), lwd = 3)
    abline(h = 0.5, v = 50, lty = 2)
    xdfttl <- "METABRIC TMA set"
    if (tlbl != "") xdfttl <- paste(xdfttl, "(", tlbl, ")" )
    title(main = xdfttl)
    if ( !is.null(binchisqlist) && !is.null(binchisqlist[[x]]) ) {
      ypos <- ifelse(binchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
      xpos <- 75
      text(xpos, ypos, paste(toupper(strsplit(x, split = "_")[[1]][1]),
                             "by AgeCat: p =",
                             format(binchisqlist[[x]]$ctout$chisq$p.value, digits = 2)))
      if ( !is.null( binchisqlist[[x]]$ctout$chisq$p.adjust ) ) {
        ypos <- ypos - 0.06
        xpos <- xpos + 7
        text(xpos, ypos, paste("p.adjust =",
                               format(binchisqlist[[x]]$ctout$chisq$p.adjust, digits = 2)))
      }
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
      lines(supsmu(ydf[idxp, "age_at_diagnosis"], ydf[idxp, x], bass = 3), lwd = 3)
      abline(h = 0.5, v = 50, lty = 2)
      ydfttl <- "Validation set"
      if (tlbl != "") ydfttl <- paste(ydfttl, "(", tlbl, ")" )
      title(main = ydfttl)
      if ( !is.null(ybinchisqlist) && !is.null(ybinchisqlist[[x]]) ) {
        ypos <- ifelse(ybinchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
        xpos <- 75
        text(xpos, ypos, paste(toupper(strsplit(x, split = "_")[[1]][1]),
                               "by AgeCat: p =",
                               format(ybinchisqlist[[x]]$ctout$chisq$p.value, digits = 2)))
        if ( !is.null( ybinchisqlist[[x]]$ctout$chisq$p.adjust ) ) {
          ypos <- ypos - 0.06
          xpos <- xpos + 7
          text(xpos, ypos, paste("p.adjust =",
                                 format(ybinchisqlist[[x]]$ctout$chisq$p.adjust, digits = 2)))
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

binsupsmu(bscvars[1], binchisqlist = bsctableout)

pdf("./Plots/MB09_BiomarkerRatesByAgev02.pdf", width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 0, 0), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu, binchisqlist = bsctableout )
dev.off()



### Luminalp
BrCa4Subtype <- "Luminalp"
lump_t_tblout <- sapply(bscvars, binchisq, xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)


pdf(paste("./Plots/MB09_", BrCa4Subtype, "_BiomarkerRatesByAgev01.pdf", sep = ""), width = 8, height = 10)
par(mfrow = c(2, 2), oma = c(2, 2, 2, 2), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu,
       xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
       binchisqlist = lump_t_tblout, 
       tlbl = BrCa4Subtype )
dev.off()



### TNP
BrCa4Subtype <- "TNP"
lump_t_tblout <- sapply(bscvars, binchisq, xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)


pdf(paste("./Plots/MB09_", BrCa4Subtype, "_BiomarkerRatesByAgev01.pdf", sep = ""), width = 8, height = 10)
par(mfrow = c(2, 2), oma = c(2, 2, 2, 2), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu,
       xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
       binchisqlist = lump_t_tblout,
       tlbl = BrCa4Subtype )
dev.off()


### HER2+       "Luminal/HER2+" "HER2+/ER-PR-"
BrCa4Subtype <- "HER2+"
lump_t_tblout <- sapply(bscvars, binchisq,
                        xdf = tmadf[((tmadf$BrCa4 == "Luminal/HER2+") |
                                     (tmadf$BrCa4 == "HER2+/ER-PR-") ), ],
                        simplify = FALSE, USE.NAMES = TRUE)

lump_v_tblout <- sapply(bscvars, binchisq,
                        xdf = vmadf[((vmadf$BrCa4 == "Luminal/HER2+") |
                                     (vmadf$BrCa4 == "HER2+/ER-PR-") ), ],
                        simplify = FALSE, USE.NAMES = TRUE)

pdf(paste("./Plots/tv_", BrCa4Subtype, "_BiomarkerRatesByAgev01.pdf", sep = ""), width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 2, 2), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu,
       xdf = tmadf[((tmadf$BrCa4 == "Luminal/HER2+") |
                    (tmadf$BrCa4 == "HER2+/ER-PR-") ), ],
       ydf = vmadf[((vmadf$BrCa4 == "Luminal/HER2+") |
                    (vmadf$BrCa4 == "HER2+/ER-PR-") ), ],
       binchisqlist = lump_t_tblout, ybinchisqlist = lump_v_tblout,
       tlbl = BrCa4Subtype )
dev.off()

### Need to adjust the p-values for the tests of association of biomarker
### with age.  40 biomarkers, training and validation, so adjust the 80 p-values for multiple comparisons.


### Need beanplots boxplots of age by subtype

require("beanplot")
win.graph()
beanplot(age_at_diagnosis ~ BrCa4f, data = tmadf,
         col = list("purple", "cyan", "orange", "red", "green"))
binchisq("BrCa4f")
### $ctout$prop.col
###          y
### x           Luminalp Luminal/HER2+ HER2+/ER-PR-        TNP Unassigned
###   <= 40   0.05299539    0.08264463   0.13333333 0.15723270 0.04225352
###   40 - 55 0.27956989    0.34710744   0.35000000 0.36792453 0.23943662
###   55 - 70 0.37403994    0.38842975   0.37500000 0.30188679 0.40845070
###   > 70    0.29339478    0.18181818   0.14166667 0.17295597 0.30985915

binchisq("BrCa4f", xdf = vmadf)
### $ctout$prop.col
###          y
### x           Luminalp Luminal/HER2+ HER2+/ER-PR-        TNP Unassigned
###   <= 40   0.04887510    0.06796117   0.12307692 0.13782051 0.09090909
###   40 - 55 0.30178433    0.35922330   0.28461538 0.37820513 0.27272727
###   55 - 70 0.38789760    0.35922330   0.44615385 0.31410256 0.38961039
###   > 70    0.26144298    0.21359223   0.14615385 0.16987179 0.24675325

cscvars <- matrix(c(
  "bcl2_pp_v1n", "Non binding",
  "birc5_pp_v1n", "Non binding",
  "cd163_c_v1n", "Non binding",
  "ck56_c2v1n", "Non binding",
  "ck5_c2v1n", "Non binding",
  "ckit_c_v1n", "Non binding",
  "cyclind1_pp_v1n", "ER binding",
  "ecad_pp_v1n", "Non binding",
  "egfr_tsipp_v1n", "Non binding",
  "er_pp_v1n", "ER binding",
  "er_dcc_v1n", "ER binding",
  "ezh2_pp_v1n", "Non binding",
  "foxa1_pp_v1n", "ER binding",
  "h3k27me3_pp_v1n", "Epigenetic",
  "her2_v1n", "Non binding",
  "inpp4b_c_v1.ppn", "Non binding",
  "ki67_pp_v1n", "Non binding",
  "nestin_tsipp_v1n", "Non binding",
  "p16_qsipp_v1n", "Non binding",
  "pr_c_v1n", "ER binding",
  "tp53_pp_v1n", "Non binding"
  ), ncol = 2, byrow = TRUE)

cscvars4 <- c(
  "er_pp_v1n",
  "ezh2_pp_v1n",
  "foxa1_pp_v1n",
  "h3k27me3_pp_v1n"
  )

cscvarsEH <- matrix(c(
  "bcl2_pp_v1n", "Non binding",
  "birc5_pp_v1n", "Non binding",
  "cd163_c_v1n", "Non binding",
  "ck56_c2v1n", "Non binding",
  "ck5_c2v1n", "Non binding",
  "ckit_c_v1n", "Non binding",
  "cyclind1_pp_v1n", "ER binding",
  "ecad_pp_v1n", "Non binding",
  "egfr_tsipp_v1n", "Non binding",
  "ezh2_pp_v1n", "Non binding",
  "foxa1_pp_v1n", "ER binding",
  "h3k27me3_pp_v1n", "Epigenetic",
  "inpp4b_c_v1.ppn", "Non binding",
  "ki67_pp_v1n", "Non binding",
  "nestin_tsipp_v1n", "Non binding",
  "p16_qsipp_v1n", "Non binding",
  "pr_c_v1n", "ER binding",
  "tp53_pp_v1n", "Non binding"
  ), ncol = 2, byrow = TRUE)

### Continuous plots
### simsmufit <- function(x, xdf = tmadf, ydf = NULL, xvar = "age_at_diagnosis",
###                       binchisqlist = NULL, ybinchisqlist = NULL,
###                       tlbl = "", nsim = 2000, alpha = 0.05) {
###   rn <- nsim
###   alphalev <- alpha/2.0
###   UCLn <- trunc(rn * (1 - alphalev))
###   LCLn <- trunc(rn * alphalev)
###   agesv <- sort(unique(xdf[, xvar]), na.last = NA)
###   simmat <- matrix(NA_real_, nrow = rn, ncol = length(agesv))
###   ylims <- range(xdf[, x], na.rm = TRUE)
###   ylimd <- abs(diff(ylims))
###   ylimlo <- ylims[1] - (0.2 * ylimd)
###   ylimhi <- ylims[2] + (0.2 * ylimd)
###   valdx <- !(is.na(xdf[, x]) | is.na(xdf[, xvar]))
###   plot( xdf[valdx, xvar], xdf[valdx, x], xlab = xvar, ylab = x, xlim = c(15, 95),
###         ylim = c(ylimlo, ylimhi), pch = ".", col = "grey60" )
### ###  axis(side = 2, at = c(0, 20, 40, 60, 80, 100))
###   for (i in seq(rn) ) {
###     ssout <- supsmu( sample(xdf[valdx, xvar]), xdf[valdx, x], bass = 10)
###     spout <- predict(smooth.spline(ssout$x, ssout$y), x = agesv)
###     simmat[i, ] <- spout$y
### ###     lines(ssout, lty = 3, lwd = 1, col = "lightgrey")
###   }
###   ##points( xdf[, xvar], xdf[, x])
###   srtmat <- apply(simmat, 2, sort)
### ###   lines(agesv, srtmat[LCLn, ], lty = 3, lwd = 3)
### ###   lines(agesv, srtmat[UCLn, ], lty = 3, lwd = 3)
###   polygon(c(agesv, rev(agesv)), c(srtmat[LCLn, ], rev(srtmat[UCLn, ])), col = "wheat3", border = "wheat3")
###   lines(supsmu( xdf[valdx, xvar], xdf[valdx, x], bass = 10), lwd = 3)
### ###   legend("bottomright", legend = paste(trunc(100*(1-alpha)), "% Confidence Bounds", sep = ""),
### ###          lty = 3, lwd = 3, seg.len = 3, text.width = 60)
### }

###
simsmufit <- function(x, xdf = tmadf, ydf = NULL, xvar = "age_at_diagnosis",
                      tlbl = "", nsim = 2000, alpha = 0.05,
                      pdfOutFilep = FALSE, pdfFileNamePrefix = NULL, plotDirName = NULL,
                      ERbindingstatus = NULL, subgroup = NULL
                      ) {
  rn <- nsim
  alphalev <- alpha/2.0
  UCLn <- trunc(rn * (1 - alphalev))
  LCLn <- trunc(rn * alphalev)
  agesv <- sort(unique(xdf[, xvar]), na.last = NA)
  simmat <- matrix(NA_real_, nrow = rn, ncol = length(agesv))
  ylims <- range(xdf[, x], na.rm = TRUE)
  ylimd <- abs(diff(ylims))
###   ylimlo <- ylims[1] - (0.5 * ylimd)
  ylimlo <- ylims[1] - (0.1 * ylimd)
  ylimhi <- ylims[2] + (0.1 * ylimd)
  valdx <- !is.na(xdf[, x])
###   xlb <- strsplit(xvar, split = "_")[[1]][1]
###   xlbc <- paste(toupper(substring(xlb, 1, 1)), tolower(substring(xlb, 2)), sep = "")
  xlbc <- "Dx age (years)"
  ylbc <- paste(toupper(strsplit(x, split = "_")[[1]][1]), "(", ERbindingstatus, ")")
  ylbnc <- paste(subgroup, ": ", length(xdf[valdx, xvar]), "cases")
  if ( pdfOutFilep ) {
      fnm <- gsub(" ", "_", paste(pdfFileNamePrefix, x, ".pdf", sep = ""))
      pdf(file=paste(plotDirName, fnm, sep = "/"), width = 6, height = 6, useDingbats = FALSE)
      on.exit(dev.off())
  }
  
  plot( xdf[valdx, xvar], xdf[valdx, x], xlab = xlbc, main = "", ylab = "",
       xlim = c(15, 95), ylim = c(ylimlo, ylimhi), type = "n" )
  points( xdf[valdx, xvar], xdf[valdx, x], pch = 1, cex = 0.3,
         col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=60, max=255)) ## was "grey60"
  title(main = ylbnc, line = 1)
  title(main = ylbc, line = 2)
###  axis(side = 2, at = c(0, 20, 40, 60, 80, 100))
  for (i in seq(rn) ) {
###     ssout <- supsmu( sample(xdf[valdx, xvar]), xdf[valdx, x], bass = 10)
###     spout <- predict(smooth.spline(ssout$x, ssout$y), x = agesv)
###     simmat[i, ] <- spout$y
      simmat[i, ] <- predict(smooth.spline(supsmu( sample(xdf[valdx, xvar]), xdf[valdx, x], bass = 10)), x = agesv)$y

###     lines(ssout, lty = 3, lwd = 1, col = "lightgrey")
  }
  ##points( xdf[, xvar], xdf[, x])
  srtmat <- apply(simmat, 2, sort)
###   lines(supsmu(agesv, srtmat[LCLn, ], bass = 10), lty = 1, lwd = 3, type = "l", col = "grey80")
###   lines(supsmu(agesv, srtmat[UCLn, ], bass = 10), lty = 1, lwd = 3, type = "l", col = "grey80")
  ## Hex colour #DBD1BA "#CDBA9680" instead of wheat3
  polygon(c(agesv, rev(agesv)), c(srtmat[LCLn, ], rev(srtmat[UCLn, ])),
          col = rgb(t(col2rgb("wheat3", alpha = FALSE)), alpha=128, max=255),
          border = rgb(t(col2rgb("wheat3", alpha = FALSE)), alpha=128, max=255)) ## Hex colour #CDBA9680
  ssout <- supsmu( xdf[valdx, xvar], xdf[valdx, x], bass = 10)
  lines(ssout, lwd = 2, type = "l", col = "#AA1010")  ## #9A1520
  srtmatcols <- match(ssout$x, agesv)
  if ( length(ssout$x) <= length(agesv) ) {
    if ( any( ssout$y < srtmat[LCLn, srtmatcols] ) |
         any( ssout$y > srtmat[UCLn, srtmatcols] ) ) {
      idxs <- c( which( ssout$y < srtmat[LCLn, srtmatcols] ),
                 which( ssout$y > srtmat[UCLn, srtmatcols] ) )
      points(ssout$x[idxs], ssout$y[idxs], pch = 20, col = "#AA1010")
      text(x = 20, y = ylimhi, labels = paste("Age-dependent p<", alpha, sep = ""), adj = 0)
    } else {
      text(x = 20, y = ylimhi, labels = paste("Not age-dependent p>", alpha, sep = ""), adj = 0)
    }
  } else {
    text(x = 20, y = ylimhi, labels = "FIX ME missing ages", adj = 0)
  }
###   legend("bottomright", legend = paste(trunc(100*(1-alpha)), "% Confidence Bounds", sep = ""),
###          lty = 1, lwd = 4, seg.len = 3, text.width = 60, col = "wheat3")
  if ( !is.null(ydf) ) {
    agesv <- sort(unique(ydf[, xvar]), na.last = NA)
    simmat <- matrix(NA_real_, nrow = rn, ncol = length(agesv))
    ylims <- range(ydf[, x], na.rm = TRUE)
    ylimd <- abs(diff(ylims))
    ylimlo <- ylims[1] - (0.2 * ylimd)
    ylimhi <- ylims[2] + (0.2 * ylimd)
    valdy <- !is.na(ydf[, x])
    plot( ydf[valdy, xvar], ydf[valdy, x], xlab = xvar, ylab = x,
         xlim = c(15, 95), ylim = c(ylimlo, ylimhi), pch = ".", col = "grey60" )
###  axis(side = 2, at = c(0, 20, 40, 60, 80, 100))
    for (i in seq(rn) ) {
      ssout <- supsmu( sample(ydf[valdy, xvar]), ydf[valdy, x], bass = 10)
      spout <- predict(smooth.spline(ssout$x, ssout$y), x = agesv)
      simmat[i, ] <- spout$y
###      lines(ssout, lty = 3, lwd = 1, col = "lightgrey")
    }
    ##points( ydf[, xvar], ydf[, x])
    srtmat <- apply(simmat, 2, sort)
    polygon(c(agesv, rev(agesv)), c(srtmat[LCLn, ], rev(srtmat[UCLn, ])), col = "wheat3", border = "wheat3")
###     lines(agesv, srtmat[LCLn, ], lty = 1, lwd = 3, type = "l", col = "grey80")
###     lines(agesv, srtmat[UCLn, ], lty = 1, lwd = 3, type = "l", col = "grey80")
    lines(supsmu( ydf[valdy, xvar], ydf[valdy, x], bass = 10), lwd = 3, type = "l", col = "#AA1010")
###     legend("bottomright", legend = paste(trunc(100*(1-alpha)), "% Confidence Bounds", sep = ""),
###            lty = 1, lwd = 4, seg.len = 3, text.width = 60, col = "wheat3")
  }
}





pdf(file = "MB09_IHC_Trends_95pct_SimCIs_shaded_v03.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
sapply(cscvars, simsmufit)
dev.off()
##BigSeries_TNP_WholeSeries_NoMB09_NoMBex_IHC_Trends_99p99pct_400kSimCIs_shaded_v01.pdf

date()
pdf(file = "MB09_IHC_Trends_95pct_2kSimCIs_shaded_v01.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = tmadf, nsim = 2000, alpha = 0.05)
dev.off()
date()

date()
plDrNm <- "./Plots/MB09_IHC_Trends_Individual"
if ( !file.exists(plDrNm) ) { dir.create(plDrNm) }
set.seed(743); sapply(cscvars, simsmufit, xdf = tmadf, nsim = 2000, alpha = 0.05,
                      pdfOutFilep = TRUE, plotDirName = plDrNm,
                      pdfFileNamePrefix = "MB09_IHC_Trends_95pct_2kSimCIs_shaded_")
date()



date()
pdf(file = "MB09_IHC_Trends_99pct_10kSimCIs_shaded_v01.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = tmadf, nsim = 10000, alpha = 0.01)
dev.off()
date()

date()
plDrNm <- "./Plots/MB09_IHC_Trends_Individual"
if ( !file.exists(plDrNm) ) { dir.create(plDrNm) }
icinm <- "All cases"
set.seed(743); sapply(seq(nrow(cscvars)),
                      function(idx) {
                          simsmufit(cscvars[idx, 1], xdf = tmadf,
                                    nsim = 10000, alpha = 0.01,
                                    pdfOutFilep = TRUE, plotDirName = plDrNm,
                                    pdfFileNamePrefix =
                                      paste("MB09_IHC_Trends_99pct_10kSimCIs_shaded_",
                                            gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                            "_", sep = ""),
                                    ERbindingstatus = cscvars[idx, 2],
                                    subgroup = icinm)
                      } )
date()





date()
pdf(file = "MB09_IHC_Trends_99p9pct_40kSimCIs_shaded_v01.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = tmadf, nsim = 40000, alpha = 0.001)
dev.off()
date()

date()
plDrNm <- "./Plots/MB09_IHC_Trends_Individual"
if ( !file.exists(plDrNm) ) { dir.create(plDrNm) }
set.seed(743); sapply(cscvars, simsmufit, xdf = tmadf, nsim = 40000, alpha = 0.001,
                      pdfOutFilep = TRUE, plotDirName = plDrNm,
                      pdfFileNamePrefix = "MB09_IHC_Trends_99p9pct_40kSimCIs_shaded_")
date()





date()
pdf(file = "MB09_IHC_Trends_99p99pct_400kSimCIs_shaded_v01.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = tmadf, nsim = 400000, alpha = 0.0001)
dev.off()
date()

date()
plDrNm <- "./Plots/MB09_IHC_Trends_Individual"
if ( !file.exists(plDrNm) ) { dir.create(plDrNm) }
set.seed(743); sapply(cscvars, simsmufit, xdf = tmadf, nsim = 400000, alpha = 0.0001,
                      pdfOutFilep = TRUE, plotDirName = plDrNm,
                      pdfFileNamePrefix = "MB09_IHC_Trends_99p99pct_400kSimCIs_shaded_")
date()



date()
### TODO: 20170322  Add to title:  EH subgroup, Sample size, ERbinding status of target
for ( ici in seq( along = levels(tmadf$BrCaEHf) ) ) { ## ici <- 1

    icinm <- levels(tmadf$BrCaEHf)[ici]
    cat("\n\n", icinm, "\n\n")
                
    plDrNm <- "./Plots/MB09_IHC_Trends_Individual"
    if ( !file.exists(plDrNm) ) { dir.create(plDrNm) }
    set.seed(963); sapply(seq(nrow(cscvarsEH)),
                                   function(idx) {
                                       simsmufit(cscvarsEH[idx, 1], xdf = tmadf[which(tmadf$BrCaEHf == icinm), ],
                                                 nsim = 10000, alpha = 0.01,
                                                 pdfOutFilep = TRUE, plotDirName = plDrNm,
                                                 pdfFileNamePrefix =
                                                   paste("MB09_IHC_Trends_99pct_10kSimCIs_shaded_",
                                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                                         "_", sep = ""),
                                                 ERbindingstatus = cscvarsEH[idx, 2],
                                                 subgroup = icinm)
                                       } )
}
date()





pdf(file = "MB09_IHC_Trends_95pct_SimCIs_4Tgts_shaded_v01.pdf", width = 8, height = 8)
par(mfrow = c(2, 2), oma = c(1,1,3,1))
sapply(cscvars4, simsmufit)
mtext("MB09", outer = TRUE)
dev.off()




dev.print(pdf, file = "MB09_IHC_Trends_95pct_SimCIs_v01.pdf", width = 7, height = 7)



### Patient characteristics
quantile(tmadf$age_at_diagnosis, na.rm = TRUE)
tmadf$Agef <- cut(tmadf$age_at_diagnosis, breaks = c(0, 30, 40, 50, 60, 70, 80, 100), right = FALSE)
table(tmadf$Agef)
table(tmadf$Agef)/sum(table(tmadf$Agef))*100

table(tmadf$BrCa4)
table(tmadf$BrCa4)/sum(table(tmadf$BrCa4))*100




quantile(tmadf$age_at_diagnosis, na.rm = TRUE)
tmadf$Agef <- cut(tmadf$age_at_diagnosis, breaks = c(0, 30, 40, 50, 60, 70, 80, 100), right = FALSE)
table(tmadf$Agef)
table(tmadf$Agef)/sum(table(tmadf$Agef))*100

table(tmadf$BrCa4f, tmadf$Agef, useNA = "always")


## meno_status
table(tmadf$meno_status)
table(tmadf$meno_status)/sum(table(tmadf$meno_status))*100


### BrCa subtype

smCrossTable(x = tmadf$pam50_subtype, y = tmadf$BrCa4)
smCrossTable(y = tmadf$pam50_subtype, x = tmadf$BrCa4)

table(tmadf$BrCa4)
table(tmadf$BrCa4)/sum(table(tmadf$BrCa4))*100

## histcat
table(tmadf$histcat)
table(tmadf$histcat)/sum(table(tmadf$histcat))*100

## grade
table(tmadf$grade)
table(tmadf$grade)/sum(table(tmadf$grade))*100

## lvnnew
table(tmadf$lvnnew)
table(tmadf$lvnnew)/sum(table(tmadf$lvnnew))*100

## erposneg
table(tmadf$erposneg)
table(tmadf$erposneg)/sum(table(tmadf$erposneg))*100

## pgr pr_c4v1
table(tmadf$pr_c4v1)
table(tmadf$pr_c4v1)/sum(table(tmadf$pr_c4v1))*100



### HER2 status  her2_v1

table(tmadf$her2_fish_result)
table(tmadf$her2_v1)
table(tmadf$her2_v1)/sum(table(tmadf$her2_v1))*100


### Clinical T N M

table(tmadf$tnm_clin_t_4cat, useNA="always")
table(tmadf$tnm_clin_t_4cat, useNA="always")/sum(table(tmadf$tnm_clin_t_4cat, useNA="always"))*100

table(tmadf$tnm_clin_n, useNA="always")
table(tmadf$tnm_clin_n, useNA="always")/sum(table(tmadf$tnm_clin_n, useNA="always"))*100

table(tmadf$tnm_clin_m, useNA="always")
table(tmadf$tnm_clin_m, useNA="always")/sum(table(tmadf$tnm_clin_m, useNA="always"))*100

### Pathological T N M

table(tmadf$tnm_surg_t_4cat, useNA="always")
table(tmadf$tnm_surg_t_4cat, useNA="always")/sum(table(tmadf$tnm_surg_t_4cat, useNA="always"))*100

table(tmadf$tnm_surg_n, useNA="always")
table(tmadf$tnm_surg_n, useNA="always")/sum(table(tmadf$tnm_surg_n, useNA="always"))*100

table(tmadf$tnm_surg_m, useNA="always")
table(tmadf$tnm_surg_m, useNA="always")/sum(table(tmadf$tnm_surg_m, useNA="always"))*100

## Tx

[143] "ad"                                    
[144] "partial"                               
[145] "complete"                              
[146] "finsurg"                               
[147] "have_complete"                         
[148] "localtx"                               
[149] "rtintent"                              
[150] "brchwrt"                               
[151] "nodalrt"                               
[152] "boost"                                 
[153] "finrt"                                 
[154] "rt"                                    
[155] "behavior"                              
[156] "chemflag"                              
[157] "any_chem"                              
[158] "chemtype"                              
[159] "chemtype_cat2"                         
[160] "init_chemo"                            
[161] "hormflag"                              
[162] "hormtype"                              
[163] "any_horm"                              
[164] "systemic"      
table(tmadf$finsurg, useNA="always")
table(tmadf$localtx, useNA="always")
table(tmadf$localtx, useNA="always")/sum(table(tmadf$localtx, useNA="always"))*100

table(tmadf$systemic, useNA="always")
table(tmadf$systemic, useNA="always")/sum(table(tmadf$systemic, useNA="always"))*100
