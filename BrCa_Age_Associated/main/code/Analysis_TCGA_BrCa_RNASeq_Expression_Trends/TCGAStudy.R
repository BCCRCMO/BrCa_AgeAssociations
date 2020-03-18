### TCGA BrCa study set
library(tidyverse)
library(here)
library(fs)
plot_dirs <- path(
  here("main/code/Analysis_TCGA_BrCa_RNASeq_Expression_Trends/Plots"),
  c("",
    "ProbeLevel/AllCases/SinglePlots",
    "ProbeLevel/ER_HER2/SinglePlots",
    "ProbeLevel/intClust/SinglePlots"
  )
)
dir_create(plot_dirs)

## a useful function: rev() for strings
strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")

leftPad0 <- function(x, length = max(nchar(x))) {
  ### Note:  for numeric data, formatC(numbervec, width = 3, flag = "0")
  ###        will left-pad numerics with zeroes.
  if ( is.numeric(x) ) {
    return( formatC(x, width = length, flag = "0") )
  } else {
    zeroes <- paste(rep("0", length), collapse = "")
    return( strReverse(substring(strReverse(paste0(zeroes, x)), 1, length)) )
  }
}

sm_manhattan <- 
  function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", gsym = "Gene_symbol",
###             col = c("gray10",  "gray60"),
            col = c("black",  "red3", "green4", "blue", "cyan3", "magenta3", "orange", "gray40"),
            chrlabs = NULL, suggestiveline = -log10(1e-05), 
            genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, plotpointidxp = rep(TRUE, nrow(x)),
            plotpointhiliteidxp = rep(FALSE, nrow(x)),  hilitelbls = NULL,  ...) {
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x))) 
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!(gsym %in% names(x))) 
        warning(paste("No gene names column found. OK unless trying to label."))
    if (!is.numeric(x[[chr]])) 
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    d$ppch <- rep(NA_integer_, nrow(d))
    d$ppch[plotpointidxp] <- 20
    if ( any ( plotpointhiliteidxp, na.rm = TRUE ) ) {
        d$ppch[which(plotpointidxp & plotpointhiliteidxp)] <- 8 }

    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])
    if (!is.null(x[[gsym]])) 
        d = transform(d, GSYM = x[[gsym]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        options(scipen = 999)
        d$pos = d$BP/1e+06
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index == 
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                  lastbase
            }
            ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR == 
                i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
            ceiling(1.1*max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logp, col = col[icol], pch = ppch, ...))
            with(d[d$index == unique(d$index)[i], ], abline(v = max(pos), 
                lty = "dashed", col = "grey80"))
            icol = icol + 1
        }
    }
    if (suggestiveline) 
        abline(h = suggestiveline, col = "blue")
    if (genomewideline) 
        abline(h = genomewideline, col = "red")
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
###         with(d.highlight, points(pos, logp, col = "green3", pch = 20, ...))
###         with(d.highlight, text(pos, 1.07*logp, labels=GSYM, cex = 0.5, ...))
###         with(d.highlight, thigmophobe.labels(pos, logp, labels=GSYM, cex = 0.8, ...))
        require("maptools"); with(d.highlight, pointLabel(pos, logp, labels=GSYM, cex = 0.8, ...))
        if (any ( plotpointhiliteidxp, na.rm = TRUE ) && !is.null(hilitelbls) ) {
            legend("topleft", legend = hilitelbls, pch = c(8, 20), lwd = 2, lty = 0)
        }
    }
}




### Data inputs

## Expression data
ealldf <- read.table(file = "../Data/data_RNA_Seq_v2_expression_median.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

### Many entries are 0.  Jitter the medians so regression model fits are not singular.  Add a small amount of noise.

annoColNames <- c("Hugo_Symbol", "Entrez_Gene_Id")
annoColidxs <- match(annoColNames, names(ealldf))
dataColidxs <- setdiff(seq(ncol(ealldf)), annoColidxs)
neDataRows <- nrow(ealldf)
neDataCols <- length(dataColidxs)
set.seed(719); jitterMat <- matrix(runif(neDataRows * neDataCols, min = 1e-06, max = 1e-02), nrow = neDataRows, ncol = neDataCols)
ealldf[, dataColidxs] <- ealldf[, dataColidxs] + jitterMat
### Drop probesets with no variation
ealldfSDs <- apply(log2(16.0 + ealldf[, -c(1, 2)]), 1, function(x) { sqrt(var(x, na.rm = TRUE) ) } )

dim(ealldf)
sum(ealldfSDs > 1e-6)
VarProbesidxs <- which(ealldfSDs > 1e-6)
edf <- ealldf[VarProbesidxs, ]
dim(edf)

### This annotation data does not have Entrez gene IDs so was not useful
### hgndf <- read.table(file = "../Data/HUGO_Biomart_GeneIDs_results.txt",
###                     header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = '"', comment = "")

hgandf <- read.table(file = "../Data/gene_with_protein_product.txt", sep = "\t", stringsAsFactors = FALSE,
                     comment = "", quote = '"', header = TRUE)

annodf <- read.delim(file = paste0("../Data/",
                       "Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09.txt"),
                     header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = '"', comment = "")

## IntClust grouping data from Cambridge
icdf <- read.table(file = "../Data/TCGAIntClust1100.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

###> table(substring(icdf$ID, 14, 16))
### 
###  01A  01B  06A 
### 1080   13    7 
### 
### Drop the 06 entries
icdf <- icdf[which(!(substring(icdf$ID, 14, 15) == "06") ), ]

icdf$tcgaid <- paste(substring(icdf$ID, 1, 12), "01", sep = "-")

eadf <- edf[, c("Hugo_Symbol", "Entrez_Gene_Id")]
table(eadf$Entrez_Gene_Id %in% annodf$Entrez)
table(eadf$Entrez_Gene_Id %in% annodf$Original_Entrez)
table(eadf$Entrez_Gene_Id %in% annodf$Entrez_Gene_ID_0)  ## Match to this for maximal number of matches
annodf$Entrez_Gene_Id <- annodf$Entrez_Gene_ID_0
annodf$Chrcr <- sapply(strsplit(annodf$Genomic_location, split = ":"), function(x) unlist(x)[1])
annodf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "", gsub("_h2_hap1", "", gsub("_random", "", annodf$Chrcr))))
annodf$Chrn <- as.numeric(gsub("chr", "", annodf$Chrc))
annodf$Chrn[annodf$Chrc == "chrX"] <- 23
annodf$Chrn[annodf$Chrc == "chrY"] <- 24
with(annodf, table(Chrcr, Chrn, useNA = "always"))
with(annodf, table(Chrc, Chrn, useNA = "always"))
annodf$CHR <- annodf$Chrn
annodf$Startcr <- sapply(strsplit(annodf$Genomic_location, split = ":"), function(x) unlist(x)[2])
annodf$Stopcr <- sapply(strsplit(annodf$Genomic_location, split = ":"), function(x) unlist(x)[3])
annodf$Startn <- as.numeric(annodf$Startcr)
annodf$Stopn <- as.numeric(annodf$Stopcr)
annodf$BP <- trunc((annodf$Startn + annodf$Stopn)/2)
annodf$SNP <- annodf$Entrez_Gene_Id

###eandf <- merge(eadf, annodf[, c("Entrez_Gene_Id",  "RefSeq_ID_0", "Accession_0", "Gene_symbol", "SNP", "CHR", "BP")], all.x = TRUE)
eadfidxp <- eadf$Entrez_Gene_Id %in% annodf$Entrez_Gene_Id
annodfidxs <- match(eadf$Entrez_Gene_Id, annodf$Entrez_Gene_Id)
eandf <- cbind(eadf[eadfidxp, ],
               annodf[annodfidxs[eadfidxp], c("Entrez_Gene_ID_0",  "RefSeq_ID_0", "Accession_0", "Gene_symbol", "SNP", "CHR", "BP")])

eandf$TCGAsymbolEQMETABRICsymbol <- (eandf$Hugo_Symbol == eandf$Gene_symbol)
head(eandf)
## pep.m$Gene <- sub(".*?GN=(.*?)( .*|$)", "\\1", pep.m$Description)

eandf$refseqID <- sub("^(.*?)(\\..*|$)",  "\\1", eandf$RefSeq_ID_0)
### table(ardf$refseqID %in% eandf$refseqID)  ## Data below
### > table(ardf$refseqID %in% eandf$refseqID)
### 
### FALSE  TRUE 
###   255   341
### Quite low match rate with refseqIDs.  Use Hugo gene name matches.

table(eandf$refseqID == "")
table(is.na(eandf$refseqID))
head(eandf[is.na(eandf$refseqID), ])
dim(eandf)  ## [1] 19038    11
### Drop unknown chromosome entries
eandf <- eandf[!is.na(eandf$CHR), ]
### Drop chromosome 24 entries
eandf <- eandf[eandf$CHR != 24, ]
table(eandf$CHR, useNA = "always")
dim(eandf)  # [1] 18950    11

## Clinical data with header info
clhdf <- read.table(file = "../Data/data_bcr_clinical_data_patient.txt",
                   header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = '"',
                   comment.char = "")
## Clinical data
cldf <- read.table(file = "../Data/data_bcr_clinical_data_patient.txt",
                   header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = '"', skip = 5,
                   comment.char = "")
names(cldf) <- clhdf[4, ]
### Keep only female cases
cldf <- cldf[cldf$GENDER == "FEMALE", ]

## No annotation data available from TCGA for this BrCa data collection.
## Get list of gene names for age related EZH2 associated genes: AgeRelated_Accession_Numbersv02.csv

ardf <- read.table(file = "../Data/AgeRelated_Accession_Numbersv02.csv",
                   header = TRUE, sep = ",", stringsAsFactors = FALSE)
arnms <- unique(unlist(strsplit(c(ardf$GeneName, ardf$IHCName), "/" )))

sum( toupper(edf$Hugo_Symbol) %in% toupper(arnms) )

tcgaarnms <- edf$Hugo_Symbol[which( toupper(edf$Hugo_Symbol) %in% toupper(arnms) )]

## EZH2 related gene list
ezdf <- read.table(file = "../Data/EZH2relatedGenesv08.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
eznms <- unique(unlist(strsplit(c(ezdf$GeneName, ezdf$IHCName), "/" )))

sum( toupper(edf$Hugo_Symbol) %in% toupper(eznms) )

tcgaeznms <- edf$Hugo_Symbol[which( toupper(edf$Hugo_Symbol) %in% toupper(eznms) )]

all(tcgaeznms  %in% tcgaarnms)
table(tcgaeznms  %in% tcgaarnms)
tcgaarnms <- union(tcgaarnms, tcgaeznms)

### ER binding genes from ChIP-Seq
ebdf <- read.table(file = "../Data/ER_binding_gene_hg19_threshold_1e-5.txt",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = '"',
                   comment.char = "")
erbnms <- sort(unique(ebdf$hgnc_symbol[!(is.na(ebdf$hgnc_symbol) | ebdf$hgnc_symbol == "")]))
length(erbnms)
## Match clinical ID to RNASeq IDs
## ER, PR HER2 status are in data fields, availability TBD.

cldf$tcgaid <- paste(cldf$PATIENT_ID, "01", sep = "-")
head(cldf$PATIENT_ID)
head(cldf$tcgaid)
head(names(edf)[-c(1, 2)])

sum( names(edf)[-c(1, 2)] %in% cldf$tcgaid )
sum( cldf$tcgaid %in% names(edf)[-c(1, 2)] )
length( names(edf)[-c(1, 2)] %in% cldf$tcgaid )
length( cldf$tcgaid %in% names(edf)[-c(1, 2)] )
length(cldf$tcgaid)

clKeepVars <- c("tcgaid", "AGE", "GENDER", "MENOPAUSE_STATUS", "ER_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS",
                "HISTORY_OTHER_MALIGNANCY", "HISTORY_NEOADJUVANT_TRTYN", "RADIATION_TREATMENT_ADJUVANT",
                "PHARMACEUTICAL_TX_ADJUVANT", "SURGICAL_PROCEDURE_FIRST", "LYMPH_NODES_EXAMINED",
                "LYMPH_NODE_EXAMINED_COUNT", "LYMPH_NODES_EXAMINED_HE_COUNT", "AJCC_TUMOR_PATHOLOGIC_PT",
                "AJCC_NODES_PATHOLOGIC_PN", "AJCC_METASTASIS_PATHOLOGIC_PM", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                "METASTATIC_SITE", "HISTOLOGICAL_DIAGNOSIS", "OS_STATUS", "OS_MONTHS" )
sdf <- cldf[, clKeepVars]
sdf$Age <- as.numeric(sdf$AGE)
sdf <- sdf[!is.na(sdf$Age), ]
dim(sdf)
### 1083 female cases to work with
sidf <- merge(sdf, icdf, all.x = TRUE, all.y = FALSE)
sematch <- match(sidf$tcgaid, names(edf))
sum(is.na(sematch))
## 4 cases in clinical data with no sample - remove
which(is.na(sematch))
sidf <- sidf[!is.na(sematch), ]
sematch <- match(sidf$tcgaid, names(edf))
sum(is.na(sematch))
dim(sidf) ## 1079 cases
## Restrict expression data to 1079 cases
edf <- edf[, c(match(c("Hugo_Symbol", "Entrez_Gene_Id"), names(edf)), match(sidf$tcgaid, names(edf)))]
dim(edf)
## Restrict expression data to annotated data in eandf
edfidxp <- edf$Entrez_Gene_Id %in% eandf$Entrez_Gene_Id
eandfidxs <- match(edf$Entrez_Gene_Id, eandf$Entrez_Gene_Id)
edfnms <- names(edf)
eandfnms <- names(eandf)
edfkeepnms <- setdiff(edfnms, eandfnms)
exandf <- cbind(eandf[eandfidxs[edfidxp], ], edf[edfidxp, edfkeepnms])
dim(exandf)  ## [1] 18950  1090


all.equal(seq(1079) + 2, match(sidf$tcgaid, names(edf))) ## TRUE  ## 1079 cases plus hugo and entrez cols
all.equal(seq(1079) + 11, match(sidf$tcgaid, names(exandf))) ## TRUE  ## 1079 cases plus 11 annotation cols

### ### Match sedf to annotation data in hgandf

sedf <- cbind(sidf, log2(16.0 + t( as.matrix(exandf[, match(sidf$tcgaid, names(exandf))]) ) ) )
names(sedf) <- c(names(sidf), exandf$Hugo_Symbol)
dim(sedf) ## [1]  1079 18976 ## Down from 1091 with 12 men out
length(c(names(sidf), exandf$Hugo_Symbol))

####################################################################
###
### ER/HER2 subsets
###

###
### BrCaEH
###

### ER+/HER2-  ER+/HER2+  ER-/HER2+  ER-/HER2-  EpHn  EpHp  EnHp  EnHn
### Cull subtypes, run regressions, . . .

sedf$BrCaEH <- rep(NA_character_, nrow(sedf))

sedf[((sedf$ER_STATUS_BY_IHC == "Positive") &
      ((sedf$IHC_HER2 == "Negative") |
       ((sedf$IHC_HER2 == "Equivocal") & (sedf$HER2_FISH_STATUS == "Negative") ) ) ),
     "BrCaEH"]  <- "ER+/HER2-"

sedf[((sedf$ER_STATUS_BY_IHC == "Positive") &
      ((sedf$IHC_HER2 == "Positive") |
       ((sedf$IHC_HER2 == "Equivocal") & (sedf$HER2_FISH_STATUS == "Positive") ) ) ),
     "BrCaEH"]  <- "ER+/HER2+"

sedf[((sedf$ER_STATUS_BY_IHC == "Negative") &
      ((sedf$IHC_HER2 == "Positive") |
       ((sedf$IHC_HER2 == "Equivocal") & (sedf$HER2_FISH_STATUS == "Positive") ) ) ),
     "BrCaEH"]  <- "ER-/HER2+"

sedf[((sedf$ER_STATUS_BY_IHC == "Negative") &
      ((sedf$IHC_HER2 == "Negative") |
       ((sedf$IHC_HER2 == "Equivocal") & (sedf$HER2_FISH_STATUS == "Negative") ) ) ),
     "BrCaEH"]  <- "ER-/HER2-"

sedf$BrCaEHf <- factor(sedf$BrCaEH, levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))

### > table(sedf$BrCaEHf, useNA = "always")
### 
### ER+/HER2- ER+/HER2+ ER-/HER2+ ER-/HER2-      <NA> 
###       538       134        43       155       209 
### 

#####################################################################
### ER/HER2 subsets
FDRalphalevels <- c(0.01, 0.05)
FCthresholds <- c(4, 2, 1.25)
Verbosep <- FALSE ## TRUE
DoMemoizep <- TRUE
SingleOutputFilep <- TRUE
arScatterPlotsp <- TRUE
EqAxesScales <- c(FALSE, TRUE)  ## Must be in order F, T so appropriate scale range can be calculated across conditions
nNotAgeDepToPlot <- 1500
propLowDens <- 0.66
propHiDens <- (1.0 - propLowDens)

for ( FDRalphai in seq(along = FDRalphalevels) ) {
    FDRalpha <- FDRalphalevels[FDRalphai]

    for ( FCthreshi in seq(along = FCthresholds) ) {
        FCthresh <- FCthresholds[FCthreshi]


        for ( EqAxesScalesi in seq(along = EqAxesScales) ) {
            EqAxesScalesp <- EqAxesScales[EqAxesScalesi]

            for ( nBrCaEHi in c( 1, length(levels(sedf$BrCaEHf)) ) ) {
                if (nBrCaEHi == 1) {
                    BrCaEHiAgeAssociatedProbesetsGenesdf <-
                      data.frame(BrCaEHi = "All cases",
                                 N_WholeSeries = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_GeneNames = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nBrCaEHi)
                                 )
                    
                } else {
                    BrCaEHiAgeAssociatedProbesetsGenesdf <-
                      data.frame(BrCaEHi = levels(sedf$BrCaEHf),
                                 N_WholeSeries = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_GeneNames = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nBrCaEHi),
                                 N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nBrCaEHi)
                                 )
                }
                
                for ( ici in seq( nBrCaEHi ) ) { ## ici <- 1

                    if ( nBrCaEHi == 1 ) {
                        icinm <- "AllCases"
                    } else {
                        icinm <- levels(sedf$BrCaEHf)[ici]
                    }

                    cat("\n\n", icinm, "\n\n")
                    
### Biomarker shows age-dependent trend if
### ( ( abs(Slope for age < 60) > log2(1.5)/35 or
###     abs(Slope for age < 60) > log2(1.5)/35 or
###     abs(Slope for age) > log2(1.5)/70 )
###   AND (p-val < 0.05/(2 * num biomarkers) ) )
### Save
###  3 slopes:,
###  3 pvalues:,
###  ILMN probeset ID: Probe_id, Probe_sequence,
###  gene name: Gene_symbol, Gene_synonyms, Synonyms_0, ILMN_Gene_0, Chromosome_0, Cytoband_0, Original_genomic_annotation
###  refseq nm_ code:  RefSeq_transcripts, RefSeq_ID_0

                    ## BrCaEH Subtypes - extract subtype dataframes
                    if ( icinm == "AllCases" ) {
                        seBrCaEHidf <- sedf
                    } else {
                        seBrCaEHidf <- sedf[sedf$BrCaEHf == icinm, ]
                        seBrCaEHidf <- seBrCaEHidf[!is.na(seBrCaEHidf$BrCaEHf), ]
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

                    cat("\n\n\n")
                    arBrCaEHioutdf_memoizedfilename <-
                      paste0("TCGA_BrCa_",
                             gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                             "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                             "_outdf.csv")

                    arBrCaEHioutdf_memoizedfilenamep <- FALSE

                    if ( (!DoMemoizep) && file.exists(arBrCaEHioutdf_memoizedfilename) ) {
                        arBrCaEHioutdf <- read.table(file = arBrCaEHioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                     sep = ",", header = TRUE, comment.char = "")
                        arBrCaEHioutdf_memoizedfilenamep <- TRUE
                    } else {
                        for ( ni in seq(along = exandf$Hugo_Symbol ) ) {
                            psi <- arBrCaEHioutdf$Hugo_Symbol[ni]
                            drci <- match(psi, names(seBrCaEHidf))
                            arlmdf$probesetni <- seBrCaEHidf[, drci]
                            if (ni %% 100 == 0 ) { cat(ni, ", ") }
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

                    }

                    ## Get single P value for volcano and manhattan plots
                    arBrCaEHioutdf$P <- apply(arBrCaEHioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                    arBrCaEHioutdf$PBHadj <- apply(arBrCaEHioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                    ## Cache data for reuse to avoid running 20000 model fits when possible
                    if ( DoMemoizep || !arBrCaEHioutdf_memoizedfilenamep ) {
                        write.csv(arBrCaEHioutdf,
                                  file = arBrCaEHioutdf_memoizedfilename )
                        write.csv(arBrCaEHioutdf[arBrCaEHioutdf$AgeDependentp == 1, ],
                                  file = paste0("TCGA_BrCa_AgeDependent_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                               "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_outdf.csv") )
                        
                        write.csv(arBrCaEHioutdf[arBrCaEHioutdf$BHadj_and_AgeDependentp, ],
                                  file = paste0("TCGA_BrCa_AgeDependent_BHadj_and_FC", FCthresh, "_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_outdf.csv") )
                    }


                    arBrCaEHioutdfarGeneSetManhattanidxp <- arBrCaEHioutdf$arGeneSet_BHadj_and_AgeDependentp & arBrCaEHioutdf$CHR <= 23
                    arBrCaEHioutdfarGeneSetManhattanidxp[is.na(arBrCaEHioutdfarGeneSetManhattanidxp)] <- FALSE
                    arBrCaEHioutdfManhattanidxp <- arBrCaEHioutdf$BHadj_and_AgeDependentp & arBrCaEHioutdf$CHR <= 23
                    arBrCaEHioutdfManhattanidxp[is.na(arBrCaEHioutdfManhattanidxp)] <- FALSE
                    ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                    arBrCaEHioutdfManhattanVarsidxp <-
                      ( apply(arBrCaEHioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                              function(x) !any(is.na(x)))  & ( arBrCaEHioutdf$CHR <= 23 ) )
                    arBrCaEHioutdfManhattanVarsidxp[!arBrCaEHioutdfManhattanVarsidxp] <- FALSE

                    ## Plot expression data by age for all probes in the Age-Related EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.

                    AgeDependent_and_BHadj_FC1.25_ProbeIds <-
                      arBrCaEHioutdf[arBrCaEHioutdf$BHadj_and_AgeDependentp, "Hugo_Symbol"]

                    ## How many AgeDependent_and_BHadj_FC1.25_ProbeIds are ER binding?

                    ERbinding_AgeDependent_and_BHadj_FC1.25_ProbeIds <-
                      intersect(erbnms, AgeDependent_and_BHadj_FC1.25_ProbeIds)

                    tcgaarnms <- sort(tcgaarnms)
                    ## Whole cohort, all cases
                    for ( singlePlotsp in c(TRUE, FALSE) ) {
                        if (!singlePlotsp) {
                            pdf(file = paste0("./Plots/ProbeLevel/ER_HER2/TCGA_BrCa_AgeRelated_",
                                             gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                             nrow(seBrCaEHidf),
                                             "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                             "_lm.pdf"), width = 8, height = 10, useDingbats = FALSE)
                        }
                        
                        xlims <- c(24, 92)
                        ylims <- c(4, 23)
                        if (!singlePlotsp) {
                            par(mfrow = c(2, 2))
                        }

                        for ( ni in seq( along = tcgaarnms ) ) {
                            arnmsi <- tcgaarnms[ni]
                            if ( ni == 1 ) { cat("\n\n### --- ", arnmsi) } else { cat(" - ", arnmsi)}
                            if (singlePlotsp) {
                              pdf(file = paste0("./Plots/ProbeLevel/ER_HER2/SinglePlots/TCGA_BrCa_AgeRelated_",
                                                gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                                nrow(seBrCaEHidf),
                                                "_", arnmsi,
                                                "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                                "_lm.pdf"),
                                    width = 5, height = 6, useDingbats = FALSE)
                            }
                            
                            plot(seBrCaEHidf[, "Age"], seBrCaEHidf[, arnmsi], type = "n",
                                 ylab = "log2(Raw expression)", xlab = "Age at diagnosis",
                                 main = "", xlim = xlims, ylim = ylims)
                            title(main = paste(arnmsi, "(",
                                               ifelse(arnmsi %in% erbnms, "ER binding", "non-ER binding"), ")"), line = 2)
                            points(seBrCaEHidf[, "Age"], seBrCaEHidf[, arnmsi],
                                   pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha = 30, max = 255) )
                            lines(supsmu(seBrCaEHidf[, "Age"], seBrCaEHidf[, arnmsi], span = 0.4, bass = 10),
                                  lwd = 3, col = "#AA1010")
                            abline(h = 0, v = 60, lty = 2)
                            text(x = 23, y = 21, adj = 0, cex = 0.85,
                                 labels = paste0("[<=60] p = ",
                                                 format(arBrCaEHioutdf[arBrCaEHioutdf$Hugo_Symbol == arnmsi, "LE60_pval"], digits = 5)))
                            text(x = 23, y = 20, adj = 0, cex = 0.85,
                                 labels = paste0("[>60] p = ",
                                                 format(arBrCaEHioutdf[arBrCaEHioutdf$Hugo_Symbol == arnmsi, "GT60_pval"], digits = 5)))
                            text(x = 23, y = 19, adj = 0, cex = 0.85,
                                 labels = paste0("[All] p = ",
                                                 format(arBrCaEHioutdf[arBrCaEHioutdf$Hugo_Symbol == arnmsi, "AllAges_pval"], digits = 5)))
                            if ( arnmsi %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
                                text(x = 92, y = 21, labels = "Age-dependent trend:", adj = 1, cex = 0.85)
                                text(x = 92, y = 19,
                                     labels = paste0("|FC|>", FCthresh, " and adjPval<",
                                                     format(FDRalpha, digits = (-log10(FDRalpha)) + 1)),
                                     adj = 1, cex = 0.85)
                            } else {
                                text(x = 92, y = 21, labels = "No detectable trend:", adj = 1, cex = 0.85)
                                text(x = 92, y = 19,
                                     labels = paste0("|FC|<", FCthresh, " or adjPval>",
                                                     format(FDRalpha, digits = (-log10(FDRalpha)) + 1)),
                                     adj = 1, cex = 0.85)
                            }
                            if (singlePlotsp) {
                                dev.off()
                            }
                        }
                        if (!singlePlotsp) {
                            dev.off()
                        }
                        cat("\n\n")
                    }
                    
                    ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                    arBrCaEHioutdfManhattanVarsidxp <-
                      ( apply(arBrCaEHioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                              function(x) !any(is.na(x)))  & ( arBrCaEHioutdf$CHR <= 23 ) )

                    ## Label the top 20 largest fold change with age association and adjusted p-val significance
                    ohdf <- arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, ][
                        order(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, "Abs_FoldChange"], decreasing = TRUE), ]
                    ohdf <- ohdf[ohdf$BHadj_and_AgeDependentp, ]

                    if ( nrow(ohdf) ) {
                        ObjsToHighlight <-
                          unique(c(ohdf[ohdf$FoldChange_Direction == "Up", ][1:10, "SNP"],
                                   ohdf[ohdf$FoldChange_Direction == "Down", ][1:10, "SNP"],
                                   ohdf[1:20, "SNP"]
                                   ))
                        ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]
                    } else {
                        ObjsToHighlight <- NULL
                    }
                    
                    pdf(file = paste0("./Plots/Manhattan/ER_HER2/TCGA_BrCa_AgeRelated_Manhattan_",
                                      gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                      "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                      ".pdf"),
                        width = 12, height = 6, useDingbats = FALSE)
                    mainTitle <- paste("TCGA BrCa", icinm)

                    if ( length( ObjsToHighlight ) ) {
                        sm_manhattan(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p="PBHadj", chrlabs = c(1:22, "X"),
                                     suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                     highlight = ObjsToHighlight,
                                     plotpointidxp = arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp,
                                                                    c( "BHadj_and_AgeDependentp")],
                                     plotpointhiliteidxp =
                                       if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                           arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(Benjamini-Hochberg adjusted P-values)",
                                     main = mainTitle )
                    } else {
                        sm_manhattan(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p="PBHadj", chrlabs = c(1:22, "X"),
                                     suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                     plotpointidxp = arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp,
                                                                    c( "BHadj_and_AgeDependentp")],
                                     plotpointhiliteidxp =
                                       if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                           arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(Benjamini-Hochberg adjusted P-values)",
                                     main = mainTitle )
                    }
                    dev.off()

### Volcano plot using same labels as METABRIC - METABRIC data in
### ../Data/AgeRelated_Volcano_2_40_3_25_4_10_All1992Cases.csv

### Modify METABRIC code for TCGA:

                    alphacrit <- FDRalpha
                    ## FCthresh <- FCthresh  ## 1.25
                    
                    MBobjsToHighlightfilen <-
                      which( grepl(paste0(gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_"),
                                   MBobjsToHighlightlistfiles) &
                             grepl(paste0("_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5))),
                                   MBobjsToHighlightlistfiles) &
                             grepl("Volcano", MBobjsToHighlightlistfiles) &
                             grepl("All1992Cases", MBobjsToHighlightlistfiles) 
                            )
                    if ( length(MBobjsToHighlightfilen) > 1 ) {
                        cat("\n\nMultiple files located\n\n")
                        cat(MBobjsToHighlightlistfiles[MBobjsToHighlightfilen])
                        cat("\n\nNumber of files located = ", length(MBobjsToHighlightfilen), "\n\n" )
                        stop("More than one highlight files found\n\n")
                    }
                    MBobjsToHighlightfile <- MBobjsToHighlightlistfiles[MBobjsToHighlightfilen]
                    
                    MBobjsToHighlight <-
                      read.csv(file = paste0("../Data/", MBobjsToHighlightfile), stringsAsFactors = FALSE)

                    pdf(file = paste0("./Plots/Volcano/ER_HER2/TCGA_BrCa_AgeRelated_Volcano_",
                                      if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                      gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                      "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                      ".pdf"),
                        width = 6, height = 6, useDingbats = FALSE)

                    mainTitle <- paste("TCGA BrCa", icinm)

                    ERbindingAgeDep <- which( ( arBrCaEHioutdf$ERbinding &
                                                (abs(arBrCaEHioutdf$Best_log2FC ) > log2(FCthresh)) &
                                                (arBrCaEHioutdf$Best_BHadj_pval < alphacrit) ) )
                    ERnonbindingAgeDep <- which( ( !arBrCaEHioutdf$ERbinding &
                                                   (abs(arBrCaEHioutdf$Best_log2FC ) > log2(FCthresh)) &
                                                   (arBrCaEHioutdf$Best_BHadj_pval < alphacrit) ) )
                    NotAgeDepidxs <- setdiff(seq(nrow(arBrCaEHioutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                    mainTitleN <- paste("Cases:", nrow(seBrCaEHidf),
                                        "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                    if ( ! EqAxesScalesp ) {
                        if  (ici == 1) {
                            xlims_vpEVS <- max(abs(range(arBrCaEHioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                            ylims_vpEVS <- c( 0, 1.1*max(-log10(arBrCaEHioutdf$Best_BHadj_pval), na.rm = TRUE) )
                        } else {
###                         xlimsi <- max(abs(range(arBrCaEHioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
###                         xlims_vpEVS[1] <- min(xlims_vpEVS[1], xlimsi[1], na.rm = TRUE)
###                         xlims_vpEVS[2] <- max(xlims_vpEVS[2], xlimsi[2], na.rm = TRUE)
                            xlims_vpEVS <- c(-log2(48), log2(48))
                            ylimsi <- c( 0, 1.1*max(-log10(arBrCaEHioutdf$Best_BHadj_pval), na.rm = TRUE) )
                            ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                            ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                        }
                    }

                    ## Get density information to thin out volcano plot
                    densityEst <- hclust(dist(cbind(arBrCaEHioutdf$Best_log2FC[NotAgeDepidxs],
                                                    jitter(-log10(arBrCaEHioutdf$Best_BHadj_pval[NotAgeDepidxs] ) ) ) ),
                                         method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * min(nNotAgeDepToPlot, nNotAgeDep))
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot, nNotAgeDep))
                    
                    lowDensityidxs <- (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityidxs <- sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityidxs), size = nHiDens)
                    ## Need also to ensure that points to be labeled from METABRIC are included
                    evcpmatchidxs <- match(MBobjsToHighlight$Entrez, arBrCaEHioutdf$Entrez_Gene_Id)
                    MBlblidxs <- evcpmatchidxs[evcpmatchidxs %in% NotAgeDepidxs]

                    plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs, MBlblidxs)
                    
                    plot(arBrCaEHioutdf$Best_log2FC, -log10(arBrCaEHioutdf$Best_BHadj_pval), type = "n", main = "", 
                         xlim =  if ( EqAxesScalesp ) {
                                     xlims_vpEVS
                                 } else {
                                     max(abs(range(arBrCaEHioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                                 },
                         xlab = "(Down)        FC        (Up)   ",
                         ylab = "-log10(BH adjusted P-values)", xaxt = "n",
                         ylim = if ( EqAxesScalesp ) {
                                    ylims_vpEVS
                                } else {
                                    c( 0, 1.1*max(-log10(arBrCaEHioutdf$Best_BHadj_pval), na.rm = TRUE) )
                                }
                         )

                    title(main = mainTitle, line = 2)
                    title(main = mainTitleN, line = 1)
                    axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), labels = c(32, 16, 8, 4, 2, 1, 2, 4, 8, 16, 32) ) 
                    
                    points(arBrCaEHioutdf$Best_log2FC[NotAgeDepidxs][plotNotAgeDepidxs],
                           -log10(arBrCaEHioutdf$Best_BHadj_pval[NotAgeDepidxs][plotNotAgeDepidxs]),
                           pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))

                    points(arBrCaEHioutdf$Best_log2FC[ERbindingAgeDep], -log10(arBrCaEHioutdf$Best_BHadj_pval[ERbindingAgeDep]),
                           pch = 20, col = rgb(t(col2rgb("red", alpha = FALSE)), alpha=30, max=255))
                    points(arBrCaEHioutdf$Best_log2FC[ERnonbindingAgeDep], -log10(arBrCaEHioutdf$Best_BHadj_pval[ERnonbindingAgeDep]),
                           pch = 20, col = rgb(t(col2rgb("blue", alpha = FALSE)), alpha=30, max=255))
                    legend("topright", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                    
                    arBrCaEHioutdf$MB_Gene_symbol <- rep("", nrow(arBrCaEHioutdf))
                    MBobjsToHighlight$TCGAmatches <- match(MBobjsToHighlight$Entrez, arBrCaEHioutdf$Entrez_Gene_Id)
                    MBinTCGAidxp <- !is.na(MBobjsToHighlight$TCGAmatches)
                    arBrCaEHioutdf[MBobjsToHighlight[MBinTCGAidxp, ]$TCGAmatches, "MB_Gene_symbol"] <-
                      MBobjsToHighlight[MBinTCGAidxp, ]$Gene_symbol

                    evcp <- arBrCaEHioutdf$Entrez_Gene_Id %in% MBobjsToHighlight$Entrez
                    ## Cull out duplicate gene names
                    evcidx <- unlist(tapply(which(evcp), arBrCaEHioutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                    if ( length(evcidx) > 0 ) {
                        evcERbidx <- evcidx[which(arBrCaEHioutdf[evcidx, "ERbinding"])]
                        evcERnbidx <- setdiff(evcidx, evcERbidx)
                        lblcol <- rep("red", length(evcidx))
                        lblcol[evcidx %in% evcERnbidx] <- "blue"
                        require("maptools"); ## for pointLabel
                        pointLabel(arBrCaEHioutdf[evcidx, ]$Best_log2FC, -log10(arBrCaEHioutdf[evcidx, ]$Best_BHadj_pval),
                                   labels = arBrCaEHioutdf[evcidx, ]$MB_Gene_symbol, cex = 0.5, col = lblcol)
                    }
                    dev.off()
                }
            }
        }
    }
}

### Fisher's exact test for enrichment
### Read memoized results, cull data for Fisher's test

FDRalphalevels <- c(0.05, 0.01)

for ( FDRalphai in seq(along = FDRalphalevels) ) { ## FDRalphai <- 1

    FDRalpha <- FDRalphalevels[FDRalphai]

    for ( nBrCaEHi in c( 1, length(levels(sedf$BrCaEHf)) ) ) {  ## nBrCaEHi  <- 1
        for ( ici in seq( nBrCaEHi ) ) { ## ici <- 1
            
            if ( nBrCaEHi == 1 ) {
                icinm <- "All cases"
                NCases <- nrow(sedf)
            } else {
                icinm <- levels(sedf$BrCaEHf)[ici]
                NCases <- nrow(sedf[sedf$BrCaEHf == icinm, ])
            }

            cat("\n\n", icinm, " : Number cases = ", NCases, " : FDR = ", FDRalpha, "\n\n")

            arBrCaEHioutdf_memoizedfilename <-
              paste0("TCGA_BrCa_",
                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                     "_outdf.csv")

            arBrCaEHioutdf <- read.table(file = arBrCaEHioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                         sep = ",", header = TRUE, comment.char = "")
            ## Get Age associated and ER binding number
            ERbindingTotal <- sum(arBrCaEHioutdf$ERbinding, na.rm = TRUE)
            AgeAssociatedTotal <- sum(arBrCaEHioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
            AgeAssociatedAndERbinding <- sum(arBrCaEHioutdf$ERbinding & arBrCaEHioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
            CaseTotal <- sum(!is.na(arBrCaEHioutdf$BHadj_and_AgeDependentp))
            FETmat <- matrix(c(AgeAssociatedAndERbinding,
                               AgeAssociatedTotal - AgeAssociatedAndERbinding,
                               ERbindingTotal - AgeAssociatedAndERbinding,
                               CaseTotal - AgeAssociatedTotal - ERbindingTotal + AgeAssociatedAndERbinding
                               ), nrow = 2, ncol = 2,
                             dimnames = list("ER_binding" = c("ER binding", "Not ER binding"),
                                             "Age association" = c("Age associated", "Not age associated")))
            require("gmodels")
            print(CrossTable(FETmat, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE, fisher = TRUE, format = "SPSS"))
            cat("\n\n", "-------------------------------------------------------------------------------", "\n\n")        
        }
    }
}

#### ---- ####

### NEED iClust memoized files of gene names identified in Volcano plots for METABRIC for intClust
###
### intClust
###

### iClust_1 . . . iClust_10
### Cull subtypes, run regressions, . . .

sedf$iClust <- paste0("iClust_", sedf$IntClust)
table(sedf$iClust, useNA = "always")




sedf$iClustf <-
  factor(sedf$iClust,
         levels = c("iClust_1", "iClust_2", "iClust_3", "iClust_4", "iClust_5",
                    "iClust_6", "iClust_7", "iClust_8", "iClust_9", "iClust_10"))

### > table(sedf$iClustf, useNA = "always")
### 
###  iClust_1  iClust_2  iClust_3  iClust_4  iClust_5  iClust_6  iClust_7  iClust_8 
###        75        38       181       165        84        60       100       145 
###  iClust_9 iClust_10      <NA> 
###        74       157         0 

######################################################################
### intClust subsets
FDRalphalevels <- c(0.01, 0.05)
FCthresholds <- c(4, 2, 1.25)
Verbosep <- FALSE ## TRUE
DoMemoizep <- TRUE
SingleOutputFilep <- TRUE
arScatterPlotsp <- TRUE
EqAxesScales <- c(FALSE, TRUE)  ## Must be in order F, T so appropriate scale range can be calculated across conditions
nNotAgeDepToPlot <- 1500
propLowDens <- 0.66
propHiDens <- (1.0 - propLowDens)

for ( FDRalphai in seq(along = FDRalphalevels) ) {
    FDRalpha <- FDRalphalevels[FDRalphai]
    
    for ( FCthreshi in seq(along = FCthresholds) ) {
        FCthresh <- FCthresholds[FCthreshi]

        for ( EqAxesScalesi in seq(along = EqAxesScales) ) {
            EqAxesScalesp <- EqAxesScales[EqAxesScalesi]

            for ( niClusti in c( 1, length(levels(sedf$iClustf)) ) ) {
                if (niClusti == 1) {
                    iClustiAgeAssociatedProbesetsGenesdf <-
                      data.frame(iClusti = "All cases",
                                 N_WholeSeries = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti)
                                 )
                    
                } else {
                    iClustiAgeAssociatedProbesetsGenesdf <-
                      data.frame(iClusti = levels(sedf$iClustf),
                                 N_WholeSeries = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                                 N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti)
                                 )
                }
                
                for ( ici in seq( niClusti ) ) { ## ici <- 1

                    if ( niClusti == 1 ) {
                        icinm <- "AllCases"
                    } else {
                        icinm <- levels(sedf$iClustf)[ici]
                    }

                    cat("\n\n", icinm, "\n\n")
                    
### Biomarker shows age-dependent trend if
### ( ( abs(Slope for age < 60) > log2(1.5)/35 or
###     abs(Slope for age < 60) > log2(1.5)/35 or
###     abs(Slope for age) > log2(1.5)/70 )
###   AND (p-val < 0.05/(2 * num biomarkers) ) )
### Save
###  3 slopes:,
###  3 pvalues:,
###  ILMN probeset ID: Probe_id, Probe_sequence,
###  gene name: Gene_symbol, Gene_synonyms, Synonyms_0, ILMN_Gene_0, Chromosome_0, Cytoband_0, Original_genomic_annotation
###  refseq nm_ code:  RefSeq_transcripts, RefSeq_ID_0

                    ## iClust Subtypes - extract subtype dataframes
                    if ( icinm == "AllCases" ) {
                        seiClustidf <- sedf
                    } else {
                        seiClustidf <- sedf[sedf$iClustf == icinm, ]
                        seiClustidf <- seiClustidf[!is.na(seiClustidf$iClustf), ]
                    }
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

                    cat("\n\n\n")
                    ariClustioutdf_memoizedfilename <-
                      paste0("TCGA_BrCa_",
                             gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                             "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                             "_outdf.csv")

                    ariClustioutdf_memoizedfilenamep <- FALSE

                    if ( (!DoMemoizep) && file.exists(ariClustioutdf_memoizedfilename) ) {
                        ariClustioutdf <- read.table(file = ariClustioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                     sep = ",", header = TRUE, comment.char = "")
                        ariClustioutdf_memoizedfilenamep <- TRUE
                    } else {
                        for ( ni in seq(along = exandf$Hugo_Symbol ) ) {
                            psi <- ariClustioutdf$Hugo_Symbol[ni]
                            drci <- match(psi, names(seiClustidf))
                            arlmdf$probesetni <- seiClustidf[, drci]
                            if (ni %% 100 == 0 ) { cat(ni, ", ") }
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
                        ariClustioutdf$LE60_BHadj_pval <- p.adjust(ariClustioutdf$LE60_pval, method="BH")
                        ariClustioutdf$GT60_BHadj_pval <- p.adjust(ariClustioutdf$GT60_pval, method="BH")
                        ariClustioutdf$AllAges_BHadj_pval <- p.adjust(ariClustioutdf$AllAges_pval, method="BH")

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

                    }

                    ## Get single P value for volcano and manhattan plots
                    ariClustioutdf$P <- apply(ariClustioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                    ariClustioutdf$PBHadj <- apply(ariClustioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                    ## Cache data for reuse to avoid running 20000 model fits when possible
                    if ( DoMemoizep || !ariClustioutdf_memoizedfilenamep ) {
                        write.csv(ariClustioutdf,
                                  file = ariClustioutdf_memoizedfilename )
                        write.csv(ariClustioutdf[ariClustioutdf$AgeDependentp == 1, ],
                                  file = paste0("TCGA_BrCa_AgeDependent_",
                                                gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                                "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                                "_outdf.csv") )
                        
                        write.csv(ariClustioutdf[ariClustioutdf$BHadj_and_AgeDependentp, ],
                                  file = paste0("TCGA_BrCa_AgeDependent_BHadj_and_FC", FCthresh, "_",
                                                gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                                "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                                "_outdf.csv") )
                    }


                    ariClustioutdfarGeneSetManhattanidxp <- ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp & ariClustioutdf$CHR <= 23
                    ariClustioutdfarGeneSetManhattanidxp[is.na(ariClustioutdfarGeneSetManhattanidxp)] <- FALSE
                    ariClustioutdfManhattanidxp <- ariClustioutdf$BHadj_and_AgeDependentp & ariClustioutdf$CHR <= 23
                    ariClustioutdfManhattanidxp[is.na(ariClustioutdfManhattanidxp)] <- FALSE
                    ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                    ariClustioutdfManhattanVarsidxp <-
                      ( apply(ariClustioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                              function(x) !any(is.na(x)))  & ( ariClustioutdf$CHR <= 23 ) )
                    ariClustioutdfManhattanVarsidxp[!ariClustioutdfManhattanVarsidxp] <- FALSE

                    ## Plot expression data by age for all probes in the Age-Related EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.

                    AgeDependent_and_BHadj_FC1.25_ProbeIds <-
                      ariClustioutdf[ariClustioutdf$BHadj_and_AgeDependentp, "Hugo_Symbol"]

                    ## How many AgeDependent_and_BHadj_FC1.25_ProbeIds are ER binding?

                    ERbinding_AgeDependent_and_BHadj_FC1.25_ProbeIds <-
                      intersect(erbnms, AgeDependent_and_BHadj_FC1.25_ProbeIds)

                    tcgaarnms <- sort(tcgaarnms)
                    ## Whole cohort, all cases
                    for ( singlePlotsp in c(TRUE, FALSE) ) {
                        if (!singlePlotsp) {
                          pdf(file = paste0("./Plots/ProbeLevel/intClust/TCGA_BrCa_AgeRelated_",
                                            gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                            nrow(seiClustidf),
                                            "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                            "_lm.pdf"), width = 8, height = 10, useDingbats = FALSE)
                        }
                        
                        xlims <- c(24, 92)
                        ylims <- c(4, 23)
                        if (!singlePlotsp) {
                            par(mfrow = c(2, 2))
                        }

                        for ( ni in seq( along = tcgaarnms ) ) {
                            arnmsi <- tcgaarnms[ni]
                            if ( ni == 1 ) { cat("\n\n### --- ", arnmsi) } else { cat(" - ", arnmsi)}
                            if (singlePlotsp) {
                              pdf(file = paste0("./Plots/ProbeLevel/intClust/SinglePlots/TCGA_BrCa_AgeRelated_",
                                                gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                                nrow(seiClustidf),
                                                "_", arnmsi,
                                                "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                                "_lm.pdf"),
                                    width = 5, height = 6, useDingbats = FALSE)
                            }
                            
                            plot(seiClustidf[, "Age"], seiClustidf[, arnmsi], type = "n",
                                 ylab = "log2(Raw expression)", xlab = "Age at diagnosis",
                                 main = "", xlim = xlims, ylim = ylims)
                            title(main = paste(arnmsi, "(",
                                               ifelse(arnmsi %in% erbnms, "ER binding", "non-ER binding"), ")"), line = 2)
                            points(seiClustidf[, "Age"], seiClustidf[, arnmsi],
                                   pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha = 30, max = 255) )
                            lines(supsmu(seiClustidf[, "Age"], seiClustidf[, arnmsi], span = 0.4, bass = 10),
                                  lwd = 3, col = "#AA1010")
                            abline(h = 0, v = 60, lty = 2)
                            text(x = 23, y = 21, adj = 0, cex = 0.85,
                                 labels = paste0("[<=60] p = ",
                                                 format(ariClustioutdf[ariClustioutdf$Hugo_Symbol == arnmsi, "LE60_pval"], digits = 5)))
                            text(x = 23, y = 20, adj = 0, cex = 0.85,
                                 labels = paste0("[>60] p = ",
                                                 format(ariClustioutdf[ariClustioutdf$Hugo_Symbol == arnmsi, "GT60_pval"], digits = 5)))
                            text(x = 23, y = 19, adj = 0, cex = 0.85,
                                 labels = paste0("[All] p = ",
                                                 format(ariClustioutdf[ariClustioutdf$Hugo_Symbol == arnmsi, "AllAges_pval"], digits = 5)))
                            if ( arnmsi %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
                                text(x = 92, y = 21, labels = "Age-dependent trend:", adj = 1, cex = 0.85)
                                text(x = 92, y = 19,
                                     labels = paste0("|FC|>", FCthresh, " and adjPval<",
                                                     format(FDRalpha, digits = (-log10(FDRalpha)) + 1)),
                                     adj = 1, cex = 0.85)
                            } else {
                                text(x = 92, y = 21, labels = "No detectable trend:", adj = 1, cex = 0.85)
                                text(x = 92, y = 19,
                                     labels = paste0("|FC|<", FCthresh, " or adjPval>",
                                                     format(FDRalpha, digits = (-log10(FDRalpha)) + 1)),
                                     adj = 1, cex = 0.85)
                            }
                            if (singlePlotsp) {
                                dev.off()
                            }
                        }
                        if (!singlePlotsp) {
                            dev.off()
                        }
                        cat("\n\n")
                    }
                    
                    ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                    ariClustioutdfManhattanVarsidxp <-
                      ( apply(ariClustioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                              function(x) !any(is.na(x)))  & ( ariClustioutdf$CHR <= 23 ) )

                    ## Label the top 20 largest fold change with age association and adjusted p-val significance
                    ohdf <- ariClustioutdf[ariClustioutdfManhattanVarsidxp, ][
                        order(ariClustioutdf[ariClustioutdfManhattanVarsidxp, "Abs_FoldChange"], decreasing = TRUE), ]
                    ohdf <- ohdf[ohdf$BHadj_and_AgeDependentp, ]

                    if ( nrow(ohdf) ) {
                        ObjsToHighlight <-
                          unique(c(ohdf[ohdf$FoldChange_Direction == "Up", ][1:10, "SNP"],
                                   ohdf[ohdf$FoldChange_Direction == "Down", ][1:10, "SNP"],
                                   ohdf[1:20, "SNP"]
                                   ))
                        ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]
                    } else {
                        ObjsToHighlight <- NULL
                    }

                    pdf(file = paste0("./Plots/Manhattan/intClust/TCGA_BrCa_AgeRelated_Manhattan_",
                                      gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                      "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                      ".pdf"),
                        width = 12, height = 6, useDingbats = FALSE)
                    mainTitle <- paste("TCGA BrCa", icinm)

                    if ( length( ObjsToHighlight ) ) {
                        sm_manhattan(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p = "PBHadj", chrlabs = c(1:22, "X"),
                                     suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                     highlight = ObjsToHighlight,
                                     plotpointidxp = ariClustioutdf[ariClustioutdfManhattanVarsidxp,
                                                                    c( "BHadj_and_AgeDependentp")],
                                     plotpointhiliteidxp =
                                       if (sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                           ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if (sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(Benjamini-Hochberg adjusted P-values)",
                                     main = mainTitle )
                    } else {
                        sm_manhattan(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p = "PBHadj", chrlabs = c(1:22, "X"),
                                     suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                     plotpointidxp = ariClustioutdf[ariClustioutdfManhattanVarsidxp,
                                                                    c( "BHadj_and_AgeDependentp")],
                                     plotpointhiliteidxp =
                                       if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                           ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(Benjamini-Hochberg adjusted P-values)",
                                     main = mainTitle )
                    }
                    dev.off()

### Volcano plot using same labels as METABRIC - METABRIC data in
### ../Data/AgeRelated_Volcano_2_40_3_25_4_10_All1992Cases.csv

### Modify METABRIC code for TCGA:

                    alphacrit <- FDRalpha
                    ## FCthresh <- FCthresh  ## 1.25
                    
                    MBobjsToHighlightfilen <-
                      which( grepl(paste0(gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_"),
                                   MBobjsToHighlightlistfiles) &
                             grepl(paste0("_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5))),
                                   MBobjsToHighlightlistfiles) &
                             grepl("Volcano", MBobjsToHighlightlistfiles) &
                             grepl("All1992Cases", MBobjsToHighlightlistfiles) 
                            )
                    if ( length(MBobjsToHighlightfilen) > 1 ) {
                        cat("\n\nMultiple files located\n\n")
                        cat(MBobjsToHighlightlistfiles[MBobjsToHighlightfilen])
                        cat("\n\nNumber of files located = ", length(MBobjsToHighlightfilen), "\n\n" )
                        stop("More than one highlight files found\n\n")
                    }
                    MBobjsToHighlightfile <- MBobjsToHighlightlistfiles[MBobjsToHighlightfilen]
                    
                    MBobjsToHighlight <-
                      read.csv(file = paste0("../Data/", MBobjsToHighlightfile), stringsAsFactors = FALSE)

                    pdf(file = paste0("./Plots/Volcano/intClust/TCGA_BrCa_AgeRelated_Volcano_",
                                      if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                      gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                      "_FC_", FCthresh, "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                      ".pdf"),
                        width = 6, height = 6, useDingbats = FALSE)
                    
                    mainTitle <- paste("TCGA BrCa", icinm)
                    
                    ERbindingAgeDep <- which( ( ariClustioutdf$ERbinding &
                                                (abs(ariClustioutdf$Best_log2FC ) > log2(FCthresh)) &
                                                (ariClustioutdf$Best_BHadj_pval < alphacrit) ) )
                    ERnonbindingAgeDep <- which( ( !ariClustioutdf$ERbinding &
                                                   (abs(ariClustioutdf$Best_log2FC ) > log2(FCthresh)) &
                                                   (ariClustioutdf$Best_BHadj_pval < alphacrit) ) )
                    NotAgeDepidxs <- setdiff(seq(nrow(ariClustioutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                    mainTitleN <- paste("Cases:", nrow(seiClustidf),
                                        "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                    if ( ! EqAxesScalesp ) {
                        if  (ici == 1) {
                            xlims_vpEVS <- max(abs(range(ariClustioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                            ylims_vpEVS <- c( 0, 1.1*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
                        } else {
###                     xlimsi <- max(abs(range(ariClustioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
###                     xlims_vpEVS[1] <- min(xlims_vpEVS[1], xlimsi[1], na.rm = TRUE)
###                     xlims_vpEVS[2] <- max(xlims_vpEVS[2], xlimsi[2], na.rm = TRUE)
                            xlims_vpEVS <- c(-log2(48), log2(48))
                            ylimsi <- c( 0, 1.1*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
                            ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                            ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                        }
                    }

                    ## Get density information to thin out volcano plot
                    densityEst <- hclust(dist(cbind(ariClustioutdf$Best_log2FC[NotAgeDepidxs],
                                                    jitter(-log10(ariClustioutdf$Best_BHadj_pval[NotAgeDepidxs] ) ) ) ),
                                         method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * min(nNotAgeDepToPlot, nNotAgeDep))
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot, nNotAgeDep))
                    
                    lowDensityidxs <- (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityidxs <- sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityidxs), size = nHiDens)
                    ## Need also to ensure that points to be labeled from METABRIC are included
                    evcpmatchidxs <- match(MBobjsToHighlight$Entrez, ariClustioutdf$Entrez_Gene_Id)
                    MBlblidxs <- evcpmatchidxs[evcpmatchidxs %in% NotAgeDepidxs]

                    plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs, MBlblidxs)

                    plot(ariClustioutdf$Best_log2FC, -log10(ariClustioutdf$Best_BHadj_pval), type = "n", main = "", 
                         xlim =  if ( EqAxesScalesp ) {
                                     xlims_vpEVS
                                 } else {
                                     max(abs(range(ariClustioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                                 },
                         xlab = "(Down)        FC        (Up)   ",
                         ylab = "-log10(BH adjusted P-values)", xaxt = "n",
                         ylim = if ( EqAxesScalesp ) {
                                    ylims_vpEVS
                                } else {
                                    c( 0, 1.1*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
                                }
                         )

                    title(main = mainTitle, line = 2)
                    title(main = mainTitleN, line = 1)
                    axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), labels = c(32, 16, 8, 4, 2, 1, 2, 4, 8, 16, 32) ) 
                    
                    points(ariClustioutdf$Best_log2FC[NotAgeDepidxs][plotNotAgeDepidxs],
                           -log10(ariClustioutdf$Best_BHadj_pval[NotAgeDepidxs][plotNotAgeDepidxs]),
                           pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                    
                    points(ariClustioutdf$Best_log2FC[ERbindingAgeDep], -log10(ariClustioutdf$Best_BHadj_pval[ERbindingAgeDep]),
                           pch = 20, col = rgb(t(col2rgb("red", alpha = FALSE)), alpha=30, max=255))
                    points(ariClustioutdf$Best_log2FC[ERnonbindingAgeDep], -log10(ariClustioutdf$Best_BHadj_pval[ERnonbindingAgeDep]),
                           pch = 20, col = rgb(t(col2rgb("blue", alpha = FALSE)), alpha=30, max=255))
                    legend("topright", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                    
                    ariClustioutdf$MB_Gene_symbol <- rep("", nrow(ariClustioutdf))
                    MBobjsToHighlight$TCGAmatches <- match(MBobjsToHighlight$Entrez, ariClustioutdf$Entrez_Gene_Id)
                    MBinTCGAidxp <- !is.na(MBobjsToHighlight$TCGAmatches)
                    ariClustioutdf[MBobjsToHighlight[MBinTCGAidxp, ]$TCGAmatches, "MB_Gene_symbol"] <-
                      MBobjsToHighlight[MBinTCGAidxp, ]$Gene_symbol

                    evcp <- ariClustioutdf$Entrez_Gene_Id %in% MBobjsToHighlight$Entrez
                    ## Cull out duplicate gene names
                    evcidx <- unlist(tapply(which(evcp), ariClustioutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                    if ( length(evcidx) > 0 ) {
                        evcERbidx <- evcidx[which(ariClustioutdf[evcidx, "ERbinding"])]
                        evcERnbidx <- setdiff(evcidx, evcERbidx)
                        lblcol <- rep("red", length(evcidx))
                        lblcol[evcidx %in% evcERnbidx] <- "blue"
                        require("maptools"); ## for pointLabel
                        pointLabel(ariClustioutdf[evcidx, ]$Best_log2FC, -log10(ariClustioutdf[evcidx, ]$Best_BHadj_pval),
                                   labels = ariClustioutdf[evcidx, ]$MB_Gene_symbol, cex = 0.5, col = lblcol)
                    }
                    dev.off()
                }
            }
        }
    }
}
    
### Fisher's exact test for enrichment
### Read memoized results, cull data for Fisher's test

FDRalphalevels <- c(0.05, 0.01)

for ( FDRalphai in seq(along = FDRalphalevels) ) { ## FDRalphai <- 1

    FDRalpha <- FDRalphalevels[FDRalphai]

    for ( niClusti in c( 1, length(levels(sedf$iClustf)) ) ) {  ## niClusti  <- 1
        for ( ici in seq( niClusti ) ) { ## ici <- 1
            
            if ( niClusti == 1 ) {
                icinm <- "All cases"
                NCases <- nrow(sedf)
            } else {
                icinm <- levels(sedf$iClustf)[ici]
                NCases <- nrow(sedf[sedf$iClustf == icinm, ])
            }

            cat("\n\n", icinm, " : Number cases = ", NCases, " : FDR = ", FDRalpha, "\n\n")

            ariClustioutdf_memoizedfilename <-
              paste0("TCGA_BrCa_",
                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                     "_outdf.csv")

                ariClustioutdf <- read.table(file = ariClustioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                             sep = ",", header = TRUE, comment.char = "")
            ## Get Age associated and ER binding number
            ERbindingTotal <- sum(ariClustioutdf$ERbinding, na.rm = TRUE)
            AgeAssociatedTotal <- sum(ariClustioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
            AgeAssociatedAndERbinding <- sum(ariClustioutdf$ERbinding & ariClustioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
            CaseTotal <- sum(!is.na(ariClustioutdf$BHadj_and_AgeDependentp))
            FETmat <- matrix(c(AgeAssociatedAndERbinding,
                               AgeAssociatedTotal - AgeAssociatedAndERbinding,
                               ERbindingTotal - AgeAssociatedAndERbinding,
                               CaseTotal - AgeAssociatedTotal - ERbindingTotal + AgeAssociatedAndERbinding
                               ), nrow = 2, ncol = 2,
                             dimnames = list("ER_binding" = c("ER binding", "Not ER binding"),
                                             "Age association" = c("Age associated", "Not age associated")))
            require("gmodels")
            print(CrossTable(FETmat, prop.r = TRUE, prop.t = FALSE, prop.chisq = FALSE, fisher = TRUE, format = "SPSS"))
            cat("\n\n", "-------------------------------------------------------------------------------", "\n\n")        
        }
    }
}

### Patient characteristics ##########-----------------------------------------------------

DropCol1c <- "HMGB1P1"
DropColNc <- "AKR1C6P"
DropCol1n <- match(DropCol1c, names(sedf))
DropColNn <- match(DropColNc, names(sedf))
secldf <- sedf[, -c(DropCol1n:DropColNn)]
head(secldf)

lapply(secldf, table, useNA="always")

table(cut(as.numeric(secldf$AGE), breaks = c(0, (2:8)*10, 100)), useNA = "always")
median(as.numeric(secldf$AGE), na.rm = TRUE)
range(as.numeric(secldf$AGE), na.rm = TRUE)

with(secldf, table(SURGICAL_PROCEDURE_FIRST, RADIATION_TREATMENT_ADJUVANT, useNA = "always"))
