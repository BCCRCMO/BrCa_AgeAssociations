# This script was adapted from TCGAStudy.R and used to regenerate TCGA_BrCa
# loess scatterplots with updated figure and file name conventions


# Load --------------------------------------------------------------------
library(tidyverse)
library(ageassn)

# Load data
data_root <- "main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends"
tcga_load(organ = "Br")

# Constants
FDRalphalevels <- c(0.01, 0.05)
FCthresholds <- c(1.25, 2, 4)
root <- "main/code/Analysis_TCGA_BrCa_RNASeq_Expression_Trends/"


# Process -----------------------------------------------------------------

# Process data
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, icdf, gender = "female")

# Add IntClust variables
sedf$iClust <- paste0("iClust", str_pad(sedf$IntClust, width = 2, pad = 0))
sedf$iClustf <- factor(sedf$iClust)

tcgaarnmsAll <- tcgaarnms
tcgaarnms <- c("ELF5")

# All Cases and ER HER2  ---------------------------------------------------

for (FDRalphai in seq_along(FDRalphalevels)) {
  FDRalpha <- FDRalphalevels[FDRalphai]
  
  for (FCthreshi in seq_along(FCthresholds)) {
    FCthresh <- FCthresholds[FCthreshi]
    
    for (nBrCaEHi in c(1, length(levels(sedf$BrCaEHf)))) {
      
      for (ici in seq(nBrCaEHi)) {
        if (nBrCaEHi == 1) {
          icinm <- "AllCases"
        } else {
          icinm <- levels(sedf$BrCaEHf)[ici]
        }
        if (icinm == "AllCases") {
          seBrCaEHidf <- sedf
        } else {
          seBrCaEHidf <- sedf[sedf$BrCaEHf == icinm,]
          seBrCaEHidf <- seBrCaEHidf[!is.na(seBrCaEHidf$BrCaEHf),]
        }
        clKeepVars <- c("tcgaid", "AGE", "ER_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS", "Age", "BrCaEH", "BrCaEHf")
        eKeepVars <- c("Hugo_Symbol", "Entrez_Gene_Id", "Gene_symbol", "SNP", "CHR", "BP", "refseqID")
        arBrCaEHioutdf <- exandf[, eKeepVars]
        arBrCaEHioutmatcolnames <- c("LE60_slope", "GT60_slope", "AllAges_slope",
                                     "LE60_pval", "GT60_pval", "AllAges_pval",
                                     "AgeDependentp")
        arBrCaEHioutmat <- matrix(NA_real_, nrow = nrow(arBrCaEHioutdf), ncol = length(arBrCaEHioutmatcolnames))
        LE60_slope_col <- match("LE60_slope", arBrCaEHioutmatcolnames)
        GT60_slope_col <- match("GT60_slope", arBrCaEHioutmatcolnames)
        AllAges_slope_col <- match("AllAges_slope", arBrCaEHioutmatcolnames)
        LE60_pval_col <- match("LE60_pval", arBrCaEHioutmatcolnames)
        GT60_pval_col <- match("GT60_pval", arBrCaEHioutmatcolnames)
        AllAges_pval_col <- match("AllAges_pval", arBrCaEHioutmatcolnames)
        AgeDependentp_col <- match("AgeDependentp", arBrCaEHioutmatcolnames)
        arBrCaEHioutdf$LE60_slope <- rep(NA_real_ , nrow(arBrCaEHioutdf))
        arBrCaEHioutdf$GT60_slope <- rep(NA_real_ , nrow(arBrCaEHioutdf))
        arBrCaEHioutdf$AllAges_slope <- rep(NA_real_ , nrow(arBrCaEHioutdf))
        arBrCaEHioutdf$LE60_pval <- rep(NA_real_ , nrow(arBrCaEHioutdf))
        arBrCaEHioutdf$GT60_pval <- rep(NA_real_ , nrow(arBrCaEHioutdf))
        arBrCaEHioutdf$AllAges_pval <- rep(NA_real_ , nrow(arBrCaEHioutdf))
        arBrCaEHioutdf$AgeDependentp <- rep(NA , nrow(arBrCaEHioutdf))
        biosigslopeLE60 <-  log2(FCthresh)/(60 - range(seBrCaEHidf$Age)[1] + 1)
        biosigslopeGT60 <-  log2(FCthresh)/(range(seBrCaEHidf$Age)[2] - 60)
        biosigslopeAllAges <-  log2(FCthresh)/(diff(range(seBrCaEHidf$Age)) + 1)
        multcompPval <- FDRalpha/(2 * nrow(arBrCaEHioutdf))
        arlmdf <- data.frame(Age = seBrCaEHidf$Age,
                             probesetni = seBrCaEHidf[, "ESR1"])
        for (ni in seq_along(exandf$Hugo_Symbol)) {
          psi <- arBrCaEHioutdf$Hugo_Symbol[ni]
          drci <- match(psi, names(seBrCaEHidf))
          arlmdf$probesetni <- seBrCaEHidf[, drci]
          if (ni %% 100 == 0) {
            cat(ni, ", ")
          }
          lmifitageLE60 <- lm(probesetni ~ Age, data = arlmdf, na.action = na.exclude, subset = Age <= 60)
          lmifitageGT60 <- lm(probesetni ~ Age, data = arlmdf, na.action = na.exclude, subset = Age > 60)
          lmifitallages <- lm(probesetni ~ Age, data = arlmdf, na.action = na.exclude)
          smrylmifitageLE60 <- summary(lmifitageLE60)
          smrylmifitageGT60 <- summary(lmifitageGT60)
          smrylmifitallages <- summary(lmifitallages)
          
          arBrCaEHioutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["Age", "Estimate"]
          arBrCaEHioutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["Age", "Estimate"]
          arBrCaEHioutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["Age", "Estimate"]
          arBrCaEHioutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["Age", "Pr(>|t|)"]
          arBrCaEHioutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["Age", "Pr(>|t|)"]
          arBrCaEHioutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["Age", "Pr(>|t|)"]
          arBrCaEHioutmat[ni, AgeDependentp_col] <-
            ( ( ( abs(arBrCaEHioutmat[ni, LE60_slope_col])    > biosigslopeLE60 ) ||
                  ( abs(arBrCaEHioutmat[ni, GT60_slope_col])    > biosigslopeGT60 ) ||
                  ( abs(arBrCaEHioutmat[ni, AllAges_slope_col]) > biosigslopeAllAges ) ) &&
                ( ( arBrCaEHioutmat[ni, LE60_pval_col]    < multcompPval ) ||
                    ( arBrCaEHioutmat[ni, GT60_pval_col]    < multcompPval ) ||
                    ( arBrCaEHioutmat[ni, AllAges_pval_col] < multcompPval ) ) )
        }
        dimnames(arBrCaEHioutmat)[[2]] <- arBrCaEHioutmatcolnames
        dimnames(arBrCaEHioutmat)[[1]] <- arBrCaEHioutdf$Hugo_Symbol
        arBrCaEHioutdf[, arBrCaEHioutmatcolnames] <- arBrCaEHioutmat
        cat("\n\n\n")
        ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
        arBrCaEHioutdf$LE60_BHadj_pval <- p.adjust(arBrCaEHioutdf$LE60_pval, method="BH")
        arBrCaEHioutdf$GT60_BHadj_pval <- p.adjust(arBrCaEHioutdf$GT60_pval, method="BH")
        arBrCaEHioutdf$AllAges_BHadj_pval <- p.adjust(arBrCaEHioutdf$AllAges_pval, method="BH")
        
        arBrCaEHioutdf$BHadj_and_AgeDependentp <-
          ( ( ( abs(arBrCaEHioutdf$LE60_slope) > biosigslopeLE60 )  &  ( arBrCaEHioutdf$LE60_BHadj_pval < FDRalpha ) ) |
              ( ( abs(arBrCaEHioutdf$GT60_slope) > biosigslopeGT60 )  &  ( arBrCaEHioutdf$GT60_BHadj_pval < FDRalpha ) ) |
              ( ( abs(arBrCaEHioutdf$AllAges_slope) > biosigslopeAllAges ) &
                  ( arBrCaEHioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
        
        arBrCaEHioutdf$BHadj_signifp <- ( ( ( arBrCaEHioutdf$LE60_BHadj_pval < FDRalpha ) |
                                              ( arBrCaEHioutdf$GT60_BHadj_pval < FDRalpha ) |
                                              ( arBrCaEHioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
        
        arBrCaEHioutdf$LE60_log2FC <- arBrCaEHioutdf$LE60_slope * (60 - range(seBrCaEHidf$Age)[1] + 1)
        arBrCaEHioutdf$GT60_log2FC <- arBrCaEHioutdf$GT60_slope * (range(seBrCaEHidf$Age)[2] - 60)
        arBrCaEHioutdf$AllAges_log2FC <- arBrCaEHioutdf$AllAges_slope * (diff(range(seBrCaEHidf$Age)) + 1)
        
        arBrCaEHioutdf$Best_log2FC <- arBrCaEHioutdf$AllAges_log2FC
        arBrCaEHioutdf$Best_BHadj_pval <- arBrCaEHioutdf$AllAges_BHadj_pval
        
        ## Find the largest significant fold change:
        
        LE60_Bestp <- ( ( abs(arBrCaEHioutdf$LE60_slope) > biosigslopeLE60 ) &
                          ( arBrCaEHioutdf$LE60_BHadj_pval < FDRalpha ) )
        arBrCaEHioutdf[LE60_Bestp, ]$Best_log2FC <- arBrCaEHioutdf[LE60_Bestp, ]$LE60_log2FC
        arBrCaEHioutdf[LE60_Bestp, ]$Best_BHadj_pval <- arBrCaEHioutdf[LE60_Bestp, ]$LE60_BHadj_pval
        
        GT60_Bestp <- ( ( abs(arBrCaEHioutdf$GT60_slope) > biosigslopeGT60 ) &
                          ( arBrCaEHioutdf$GT60_BHadj_pval < FDRalpha ) &
                          ( abs(arBrCaEHioutdf$GT60_slope) > abs(arBrCaEHioutdf$LE60_slope) ) )
        arBrCaEHioutdf[GT60_Bestp, ]$Best_log2FC <- arBrCaEHioutdf[GT60_Bestp, ]$GT60_log2FC
        arBrCaEHioutdf[GT60_Bestp, ]$Best_BHadj_pval <- arBrCaEHioutdf[GT60_Bestp, ]$GT60_BHadj_pval
        
        AllAges_Bestp <- ( ( abs(arBrCaEHioutdf$AllAges_slope) > biosigslopeAllAges ) &
                             ( arBrCaEHioutdf$AllAges_BHadj_pval < FDRalpha ) &
                             ( ( abs(arBrCaEHioutdf$AllAges_log2FC) > abs(arBrCaEHioutdf$LE60_log2FC) ) |
                                 ( abs(arBrCaEHioutdf$AllAges_log2FC) > abs(arBrCaEHioutdf$GT60_log2FC) ) ) )
        arBrCaEHioutdf[AllAges_Bestp, ]$Best_log2FC <- arBrCaEHioutdf[AllAges_Bestp, ]$AllAges_log2FC
        arBrCaEHioutdf[AllAges_Bestp, ]$Best_BHadj_pval <- arBrCaEHioutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
        arBrCaEHioutdf$Abs_Best_log2FC <- abs( arBrCaEHioutdf$Best_log2FC )
        arBrCaEHioutdf$Abs_FoldChange <- 2^arBrCaEHioutdf$Abs_Best_log2FC
        arBrCaEHioutdf$FoldChange_Direction <- ifelse(arBrCaEHioutdf$Best_log2FC > 0, "Up", "Down")
        
        ## Use this Best_BHadj_pval for manhattan plots as well.
        
        arBrCaEHioutdf$ERbinding <- arBrCaEHioutdf$Hugo_Symbol %in% erbnms
        arBrCaEHioutdf$arGeneSet_BHadj_and_AgeDependentp <- ((arBrCaEHioutdf$Gene_symbol %in% arnms) &
                                                               arBrCaEHioutdf$BHadj_and_AgeDependentp)
        
        ## Get single P value for volcano and manhattan plots
        arBrCaEHioutdf$P <- apply(arBrCaEHioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
        arBrCaEHioutdf$PBHadj <- apply(arBrCaEHioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
        
        ## Plot expression data by age for all probes in the Age-Related EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.
        AgeDependent_and_BHadj_FC1.25_ProbeIds <-
          arBrCaEHioutdf[arBrCaEHioutdf$BHadj_and_AgeDependentp, "Hugo_Symbol"]
        tcgaarnms <- sort(tcgaarnms)
        
        ## Whole cohort, all cases
        for (singlePlotsp in c(TRUE, FALSE)) {
          if (icinm == "AllCases") {
            subdir <- "AllCases"
            case <- paste0("All", nrow(seBrCaEHidf), "Cases")
          } else {
            subdir <- "ER_HER2"
            case <- paste0(str_sign(icinm, type = "long"), "_", nrow(seBrCaEHidf), "Cases")
          }
          
          if (!singlePlotsp) {
            pdf(file = paste0(root, "Plots/ProbeLevel/", subdir, "/TCGA_BrCa_AgeRelated_",
                              case, "_FC_", FCthresh, "_FDR_", str_decimal(format(FDRalpha, digits = 5)),
                              "_loess.pdf"),
                width = 8, height = 10, useDingbats = FALSE)
          }
          
          xlims <- c(24, 92)
          ylims <- c(4, 17)
          if (!singlePlotsp) {
            par(mfrow = c(2, 2))
          }
          
          for (ni in seq_along(tcgaarnms)) {
            arnmsi <- tcgaarnms[ni]
            
            if (ni == 1) {
              cat("\n\n### --- ", arnmsi)
            } else {
              cat(" - ", arnmsi)
            }
            if (singlePlotsp) {
              pdf(file = paste0(root, "Plots/ProbeLevel/", subdir, "/SinglePlots/TCGA_BrCa_AgeRelated_",
                                case, "_", arnmsi, "_FC_", FCthresh, "_FDR_", str_decimal(format(FDRalpha, digits = 5)),
                                "_loess.pdf"),
                  width = 5, height = 6, useDingbats = FALSE)
            }
            regression_plot(
              data = seBrCaEHidf,
              x = "Age",
              y = arnmsi,
              age = 60,
              fc = FCthresh,
              fdr = FDRalpha,
              ylab = "log2(Raw expression)",
              xlim = xlims,
              ylim = ylims,
              xlim.pval = rep(23, 3),
              ylim.pval = rev(c(15.1, 15.8, 16.5)),
              xlim.trend = 92,
              ylim.trend = c(16.5, 15.5),
              conf.level = 1 - FDRalpha,
              show.avg = TRUE,
              annotate = AgeDependent_and_BHadj_FC1.25_ProbeIds,
              erbtxt = ifelse(arnmsi %in% erbnms, "(ER binding)", "(non-ER binding)")
            )
            if (singlePlotsp) {
              dev.off()
            }
          }
          if (!singlePlotsp) {
            dev.off()
          }
          cat("\n\n")
        }
      }
    }
  }
}


# IntClust ----------------------------------------------------------------
for (FDRalphai in seq_along(FDRalphalevels)) {
  FDRalpha <- FDRalphalevels[FDRalphai]
  
  for (FCthreshi in seq_along(FCthresholds)) {
    FCthresh <- FCthresholds[FCthreshi]
    
    for (niClusti in length(levels(sedf$iClustf))) {
      
      for (ici in seq(niClusti)) {
        icinm <- levels(sedf$iClustf)[ici]
        cat("\n\n", icinm, "\n\n")
        seiClustidf <- sedf[sedf$iClustf == icinm, ]
        seiClustidf <- seiClustidf[!is.na(seiClustidf$iClustf), ]
        
        clKeepVars <- c("tcgaid", "AGE", "ER_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS", "Age", "iClust", "iClustf", "IntClust")
        eKeepVars <- c("Hugo_Symbol", "Entrez_Gene_Id", "Gene_symbol", "SNP", "CHR", "BP", "refseqID")
        ariClustioutdf <- exandf[, eKeepVars]
        ariClustioutmatcolnames <- c("LE60_slope", "GT60_slope", "AllAges_slope",
                                     "LE60_pval", "GT60_pval", "AllAges_pval",
                                     "AgeDependentp")
        ariClustioutmat <- matrix(NA_real_, nrow = nrow(ariClustioutdf), ncol = length(ariClustioutmatcolnames))
        LE60_slope_col <- match("LE60_slope", ariClustioutmatcolnames)
        GT60_slope_col <- match("GT60_slope", ariClustioutmatcolnames)
        AllAges_slope_col <- match("AllAges_slope", ariClustioutmatcolnames)
        LE60_pval_col <- match("LE60_pval", ariClustioutmatcolnames)
        GT60_pval_col <- match("GT60_pval", ariClustioutmatcolnames)
        AllAges_pval_col <- match("AllAges_pval", ariClustioutmatcolnames)
        AgeDependentp_col <- match("AgeDependentp", ariClustioutmatcolnames)
        ariClustioutdf$LE60_slope <- rep(NA_real_ , nrow(ariClustioutdf))
        ariClustioutdf$GT60_slope <- rep(NA_real_ , nrow(ariClustioutdf))
        ariClustioutdf$AllAges_slope <- rep(NA_real_ , nrow(ariClustioutdf))
        ariClustioutdf$LE60_pval <- rep(NA_real_ , nrow(ariClustioutdf))
        ariClustioutdf$GT60_pval <- rep(NA_real_ , nrow(ariClustioutdf))
        ariClustioutdf$AllAges_pval <- rep(NA_real_ , nrow(ariClustioutdf))
        ariClustioutdf$AgeDependentp <- rep(NA , nrow(ariClustioutdf))
        biosigslopeLE60 <-  log2(FCthresh)/(60 - range(seiClustidf$Age)[1] + 1)
        biosigslopeGT60 <-  log2(FCthresh)/(range(seiClustidf$Age)[2] - 60)
        biosigslopeAllAges <-  log2(FCthresh)/(diff(range(seiClustidf$Age)) + 1)
        multcompPval <- FDRalpha/(2 * nrow(ariClustioutdf))
        arlmdf <- data.frame(Age = seiClustidf$Age,
                             probesetni = seiClustidf[, "ESR1"])
        
        for (ni in seq_along(exandf$Hugo_Symbol)) {
          psi <- ariClustioutdf$Hugo_Symbol[ni]
          drci <- match(psi, names(seiClustidf))
          arlmdf$probesetni <- seiClustidf[, drci]
          if (ni %% 100 == 0) {
            cat(ni, ", ")
          }
          lmifitageLE60 <- lm(probesetni ~ Age, data = arlmdf, na.action = na.exclude, subset = Age <= 60)
          lmifitageGT60 <- lm(probesetni ~ Age, data = arlmdf, na.action = na.exclude, subset = Age > 60)
          lmifitallages <- lm(probesetni ~ Age, data = arlmdf, na.action = na.exclude)
          smrylmifitageLE60 <- summary(lmifitageLE60)
          smrylmifitageGT60 <- summary(lmifitageGT60)
          smrylmifitallages <- summary(lmifitallages)
          
          ariClustioutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["Age", "Estimate"]
          ariClustioutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["Age", "Estimate"]
          ariClustioutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["Age", "Estimate"]
          ariClustioutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["Age", "Pr(>|t|)"]
          ariClustioutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["Age", "Pr(>|t|)"]
          ariClustioutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["Age", "Pr(>|t|)"]
          ariClustioutmat[ni, AgeDependentp_col] <-
            ( ( ( abs(ariClustioutmat[ni, LE60_slope_col])    > biosigslopeLE60 ) ||
                  ( abs(ariClustioutmat[ni, GT60_slope_col])    > biosigslopeGT60 ) ||
                  ( abs(ariClustioutmat[ni, AllAges_slope_col]) > biosigslopeAllAges ) ) &&
                ( ( ariClustioutmat[ni, LE60_pval_col]    < multcompPval ) ||
                    ( ariClustioutmat[ni, GT60_pval_col]    < multcompPval ) ||
                    ( ariClustioutmat[ni, AllAges_pval_col] < multcompPval ) ) )
        }
        dimnames(ariClustioutmat)[[2]] <- ariClustioutmatcolnames
        dimnames(ariClustioutmat)[[1]] <- ariClustioutdf$Hugo_Symbol
        ariClustioutdf[, ariClustioutmatcolnames] <- ariClustioutmat
        cat("\n\n\n")
        ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
        ariClustioutdf$LE60_BHadj_pval <- p.adjust(ariClustioutdf$LE60_pval, method = "BH")
        ariClustioutdf$GT60_BHadj_pval <- p.adjust(ariClustioutdf$GT60_pval, method = "BH")
        ariClustioutdf$AllAges_BHadj_pval <- p.adjust(ariClustioutdf$AllAges_pval, method = "BH")
        
        ariClustioutdf$BHadj_and_AgeDependentp <-
          ( ( ( abs(ariClustioutdf$LE60_slope) > biosigslopeLE60 )  &  ( ariClustioutdf$LE60_BHadj_pval < FDRalpha ) ) |
              ( ( abs(ariClustioutdf$GT60_slope) > biosigslopeGT60 )  &  ( ariClustioutdf$GT60_BHadj_pval < FDRalpha ) ) |
              ( ( abs(ariClustioutdf$AllAges_slope) > biosigslopeAllAges ) &
                  ( ariClustioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
        
        ariClustioutdf$BHadj_signifp <- ( ( ( ariClustioutdf$LE60_BHadj_pval < FDRalpha ) |
                                              ( ariClustioutdf$GT60_BHadj_pval < FDRalpha ) |
                                              ( ariClustioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
        
        ariClustioutdf$LE60_log2FC <- ariClustioutdf$LE60_slope * (60 - range(seiClustidf$Age)[1] + 1)
        ariClustioutdf$GT60_log2FC <- ariClustioutdf$GT60_slope * (range(seiClustidf$Age)[2] - 60)
        ariClustioutdf$AllAges_log2FC <- ariClustioutdf$AllAges_slope * (diff(range(seiClustidf$Age)) + 1)
        
        ariClustioutdf$Best_log2FC <- ariClustioutdf$AllAges_log2FC
        ariClustioutdf$Best_BHadj_pval <- ariClustioutdf$AllAges_BHadj_pval
        
        ## Find the largest significant fold change:
        
        LE60_Bestp <- ( ( abs(ariClustioutdf$LE60_slope) > biosigslopeLE60 ) &
                          ( ariClustioutdf$LE60_BHadj_pval < FDRalpha ) )
        ariClustioutdf[LE60_Bestp, ]$Best_log2FC <- ariClustioutdf[LE60_Bestp, ]$LE60_log2FC
        ariClustioutdf[LE60_Bestp, ]$Best_BHadj_pval <- ariClustioutdf[LE60_Bestp, ]$LE60_BHadj_pval
        
        GT60_Bestp <- ( ( abs(ariClustioutdf$GT60_slope) > biosigslopeGT60 ) &
                          ( ariClustioutdf$GT60_BHadj_pval < FDRalpha ) &
                          ( abs(ariClustioutdf$GT60_slope) > abs(ariClustioutdf$LE60_slope) ) )
        ariClustioutdf[GT60_Bestp, ]$Best_log2FC <- ariClustioutdf[GT60_Bestp, ]$GT60_log2FC
        ariClustioutdf[GT60_Bestp, ]$Best_BHadj_pval <- ariClustioutdf[GT60_Bestp, ]$GT60_BHadj_pval
        
        AllAges_Bestp <- ( ( abs(ariClustioutdf$AllAges_slope) > biosigslopeAllAges ) &
                             ( ariClustioutdf$AllAges_BHadj_pval < FDRalpha ) &
                             ( ( abs(ariClustioutdf$AllAges_log2FC) > abs(ariClustioutdf$LE60_log2FC) ) |
                                 ( abs(ariClustioutdf$AllAges_log2FC) > abs(ariClustioutdf$GT60_log2FC) ) ) )
        ariClustioutdf[AllAges_Bestp, ]$Best_log2FC <- ariClustioutdf[AllAges_Bestp, ]$AllAges_log2FC
        ariClustioutdf[AllAges_Bestp, ]$Best_BHadj_pval <- ariClustioutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
        ariClustioutdf$Abs_Best_log2FC <- abs( ariClustioutdf$Best_log2FC )
        ariClustioutdf$Abs_FoldChange <- 2^ariClustioutdf$Abs_Best_log2FC
        ariClustioutdf$FoldChange_Direction <- ifelse(ariClustioutdf$Best_log2FC > 0, "Up", "Down")
        
        ## Use this Best_BHadj_pval for manhattan plots as well.
        
        ariClustioutdf$ERbinding <- ariClustioutdf$Hugo_Symbol %in% erbnms
        ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp <- ((ariClustioutdf$Gene_symbol %in% arnms) &
                                                               ariClustioutdf$BHadj_and_AgeDependentp)
        ## Get single P value for volcano and manhattan plots
        ariClustioutdf$P <- apply(ariClustioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
        ariClustioutdf$PBHadj <- apply(ariClustioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
        
        ## Plot expression data by age for all probes in the Age-Related EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.
        
        AgeDependent_and_BHadj_FC1.25_ProbeIds <-
          ariClustioutdf[ariClustioutdf$BHadj_and_AgeDependentp, "Hugo_Symbol"]
        tcgaarnms <- sort(tcgaarnms)
        
        for (singlePlotsp in c(TRUE, FALSE)) {
          subdir <- "intClust"
          case <- paste0(icinm, "_", nrow(seiClustidf), "Cases")
          
          if (!singlePlotsp) {
            pdf(file = paste0(root, "Plots/ProbeLevel/intClust/TCGA_BrCa_AgeRelated_",
                              case, "_FC_", FCthresh, "_FDR_", str_decimal(format(FDRalpha, digits = 5)),
                              "_loess.pdf"),
                width = 8, height = 10, useDingbats = FALSE)
          }
          
          xlims <- c(24, 92)
          ylims <- c(4, 17)
          if (!singlePlotsp) {
            par(mfrow = c(2, 2))
          }
          
          for (ni in seq_along(tcgaarnms)) {
            arnmsi <- tcgaarnms[ni]
            
            if (ni == 1) {
              cat("\n\n### --- ", arnmsi)
            } else {
              cat(" - ", arnmsi)
            }
            if (singlePlotsp) {
              pdf(file = paste0(root, "Plots/ProbeLevel/intClust/SinglePlots/TCGA_BrCa_AgeRelated_",
                                case, "_", arnmsi, "_FC_", FCthresh, "_FDR_", str_decimal(format(FDRalpha, digits = 5)),
                                "_loess.pdf"),
                  width = 5, height = 6, useDingbats = FALSE)
            }
            regression_plot(
              data = seiClustidf,
              x = "Age",
              y = arnmsi,
              age = 60,
              fc = FCthresh,
              fdr = FDRalpha,
              ylab = "log2(Raw expression)",
              xlim = xlims,
              ylim = ylims,
              xlim.pval = rep(23, 3),
              ylim.pval = rev(c(15.1, 15.8, 16.5)),
              xlim.trend = 92,
              ylim.trend = c(16.5, 15.5),
              conf.level = 1 - FDRalpha,
              show.avg = TRUE,
              annotate = AgeDependent_and_BHadj_FC1.25_ProbeIds,
              erbtxt = ifelse(arnmsi %in% erbnms, "(ER binding)", "(non-ER binding)")
            )
            if (singlePlotsp) {
              dev.off()
            }
          }
          if (!singlePlotsp) {
            dev.off()
          }
          cat("\n\n")
        }
      }
    }
  }
}


# Remove old files --------------------------------------------------------

all_files <-
  list.files(
    path = paste0(root, "Plots/ProbeLevel/AllCases/"),
    pattern = "lm",
    full.names = TRUE,
    recursive = TRUE
  )
file.remove(all_files)

er_her2_files <-
  list.files(
    path = paste0(root, "Plots/ProbeLevel/ER_HER2/"),
    pattern = "lm",
    full.names = TRUE,
    recursive = TRUE
  )
file.remove(er_her2_files)

intclust_files <-
  list.files(
    path = paste0(root, "Plots/ProbeLevel/intClust/"),
    pattern = "lm",
    full.names = TRUE,
    recursive = TRUE
  )
file.remove(intclust_files)
