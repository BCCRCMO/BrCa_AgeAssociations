### METABRIC study Expression data analysis for Age Associated targets.

### https://ghr.nlm.nih.gov/primer/basics/gene
### NATURE |VOL 431 | 21 OCTOBER 2004 |www.nature.com/nature
### The Human Genome Project has estimated that humans have between 20,000 and 25,000 genes

### Human Molecular Genetics, 2014 1–13
### doi:10.1093/hmg/ddu309
### Multiple evidence strands suggest that there
### may be as few as 19 000 human protein-coding genes



library(here)
library(fs)
library(gmodels)  ## for e.g. CrossTable()


### Make directory to hold plots
plot_dirs <- path(here("main/code/Analysis_METABRIC_Expression_Trends/Plots"),
                  c("", "ProbeLevel"))
dir_create(plot_dirs)


### Constants and other variables for use in analysis below.

PAM50plotcols <- c("purple", "cyan", "orange", "red", "green")
iClustplotcols <- c("#FF5500", "#00EE76", "#CD3278", "#00C5CD", "#8B0000",
                    "#FFFF40", "#0000CD", "#FFAA00", "#EE82EE", "#7D26CD")

### Gene signature inputs
### e.g. from http://www.genecards.org/cgi-bin/carddisp.pl?gene=CDKN1A
### REFSEQ mRNAs for CDKN1A gene (4 alternative transcripts): 
###     NM_000389.4  NM_001220777.1  NM_001220778.1  NM_078467.2  

prdf <- read.table(file = "../Data/PAT_Accession_Numbers_AllAlternatives.csv",
                   header = TRUE, sep = ",", stringsAsFactors = FALSE)

### Standard known list of breast cancer signature gene targets
### ER, PR, HER2, Ki67, EGFR, CK5/6, AURKA, Survivin
brcadf <- read.csv(file = "../Data/BrCaAndNegSignatureGeneTargets.csv", stringsAsFactors = FALSE)

### Cases not in the MB09 TMA
### When studying METABRIC expression (MBex) alongside the Big Series and the MB09 TMA,
### individual patients that show up in all three data sets must be handled.
### Age-associated study uses Big Series, MB09 TMA, and MBex.
### MB09 TMA keeps all cases.
### MBex loses MB09 TMA cases.
### Big Series loses MB09 TMA and MBex cases.
### MBex_MBid_NoOverlap_MB09_BigSeries.csv
### Age study unique df  asudf$MBid are the 1161 cases not in MB09 TMA.
asudf <- read.csv(file = "../Data/MBex_MBid_NoOverlap_MB09_BigSeries.csv", stringsAsFactors = FALSE)

### Transcript - Protein correlation targets
tpcdf <- read.csv(file = "../Data/AgeRelatedTranscriptProteinAssocVars.csv", stringsAsFactors = FALSE)
tpcdf <- tpcdf[!is.na(tpcdf$Illumina_ProbeID), ]

### ER_binding_gene_hg19_threshold_1e-5.txt
### ER binding genes from ChIP-Seq
ebdf <- read.table(file = "../Data/ER_binding_gene_hg19_threshold_1e-5.txt",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = '"',
                   comment.char = "")
erbnms <- sort(unique(ebdf$hgnc_symbol[!(is.na(ebdf$hgnc_symbol) | ebdf$hgnc_symbol == "")]))
length(erbnms)
### ensemble gene id is in illumina annotation




### Illumina annotation data
### NOTE annodf not in same order as metabric data mbexdf
annodf <- read.delim(file = paste("../Data/",
                       "Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09.txt",
                       sep = ""),
                     header = TRUE, stringsAsFactors = FALSE, sep = "\t")
### Set up commonly named ProbeId variable
annodf$ProbeId <- annodf$Probe_id
### > dim(annodf)
### [1] 48803   110

### Additional annotation example file derived during METABRIC study
###
### ANNOTATION <- read.csv(paste("../Data/",
###                              "ANNOTATION.csv", sep = ""), header = TRUE)


### metabric ID coding schema
### 
### 
### Addenbrooke's
### MB-AD-0000	MB-0000
### 
### Nottingham
### MB-NT-2506	MB-2506
### MB-NT-3001	MB-3001
### 
### Vancouver
### MB-VC-4000	MB-4000
### MB-VC-5000	MB-5000
### 
### Manitoba (Add 6000)
### MB-MT-0001	MB-6001
### MT-MB-0001	MB-6001
### 
### Guy's (Add 6000)
### MB-GU-1000	MB-7000

### 
### Need a function to take in any ID coding and spit out new coding.
###
setMBid <- function(x) {
  ## setMBid:  Calculate MBid.new based on given ID.
  ##
  ## First, substitute a period for minus and underscore.
  ## Switch MT.MB to be MB.MT
  ## Drop "MB." prefix
  xu <- toupper(x)
  xr <- gsub("MB\\.", "", gsub("MT\\.MB", "MB\\.MT", gsub( "-", ".", gsub( "_", ".", xu ) ) ))
  ## Extract final digits
  posdot <- unlist(regexec("\\.", xr))
  xrrs <- xr
  if ( length( idxs <- which(posdot > 0) ) ) { xrrs[idxs] <- substring(xr[idxs], posdot[idxs] + 1) }
  xn <- as.numeric(xrrs)
  ## Add 6000 for Manitoba and Guy's
  if ( length( sidxs <- c( grep(toupper("MT"), xu), grep(toupper("GU"), xu) ) ) ) { xn[sidxs] <- xn[sidxs] + 6000 }
  return(paste("MB.", formatC(xn, width = 4, flag = "0"), sep = ""))
}

### Need a function to take in any ID coding and spit out treatment centre.
### 
MBid2Ctr <- function(x) {
  if ( any( is.na( as.numeric( gsub("MB\\.", "", x) ) ) ) ) {
    stop("x must be MBid")
  }
  cut(as.numeric(gsub(toupper("MB\\."), "", toupper(x))),
      c(0, 2000, 4000, 6000, 7000, 10000) - 1,
      labels = c("AD", "NT", "VC", "MT", "GU")) 
}

### Fields with refseq ids
### RefSeq_transcripts, Proportion_RefSeq_transcripts, Proportion_UCSC_transcripts, GenBank_transcripts, 
### Original_transcriptomic_annotation, Lumi_transcriptomic_annotation, Transcriptomic_2nd_matches, Gene_symbol,
### Gene_description, Source_Reference_ID_0, RefSeq_ID_0, Accession_0, Definition_0
### 
### Choose fields as appropriate for study.
### 
### NOTE annodf not in same order as metabric data mbexdf

refSeqRowMatches_annodf <- function(x, mdf = annodf) {
  refSeqCols <- c("RefSeq_transcripts",
                  "Proportion_RefSeq_transcripts",
                  "Proportion_UCSC_transcripts",
                  "GenBank_transcripts",
                  "Original_transcriptomic_annotation",
                  "Lumi_transcriptomic_annotation",
                  "Transcriptomic_2nd_matches",
                  "Transcriptomic_match_overall_._similarity",
                  "Gene_symbol",
                  "Gene_description",
                  "Source_Reference_ID_0",
                  "RefSeq_ID_0",
                  "Accession_0",
                  "Definition_0"
                  )
  rsidxs <- unique(unlist(lapply(refSeqCols,
                       function(x, pattrn, adf) {
                         ## Strip decimal from end of RefSeq ID
                         pattrn <- strsplit(pattrn, split = "\\.")[[1]][1]
                         c(
                           grep(pattern = paste("^", pattrn, "$", sep = ""), x = adf[, x]),
                           grep(pattern = paste("^", pattrn, "\\.", sep = ""), x = adf[, x])
                           )
                       },
                       x, mdf)))
}

geneNameRowMatches_annodf <- function(x, mdf = annodf) {
    ## x:  Character vector of gene symbol names
  geneNameCols <-
    c(
      "Proportion_RefSeq_transcripts",
      "Proportion_UCSC_transcripts",
      "Gene_symbol",
      "Gene_synonyms",
      "Gene_description",
      "Original_genic_annotation",
      "Other_gene_ids",
      "X2nd_match_gene_ids",
      "Symbol_0",
      "Definition_0",
      "Synonyms_0",
      "Obsolete_Probe_Id_0"
      )
  rsidxs <- unique(unlist(lapply(geneNameCols,
                       function(x, pattrn, adf) {
                         ## Strip decimal from end of RefSeq ID
                         pattrn <- strsplit(pattrn, split = "\\.")[[1]][1]
                         #browser()
                         c(
                           grep(pattern = paste("\\(", pattrn, "\\)", sep = ""), x = adf[, x]),
                           grep(pattern = paste("^", pattrn, "$", sep = ""), x = adf[, x]),
                           grep(pattern = paste("\\s", pattrn, "$", sep = ""), x = adf[, x]),
                           grep(pattern = paste("^", pattrn, "\\s", sep = ""), x = adf[, x]),
                           grep(pattern = paste("^", pattrn, "\\;", sep = ""), x = adf[, x]),
                           grep(pattern = paste("\\s", pattrn, "\\;", sep = ""), x = adf[, x])
                           )
                       },
                       x, mdf)))
}


### Utility functions for use in heatmaps
wclust <- function(x) hclust(d = x, method = "ward")
sclust <- function(x) hclust(d = x, method = "single")
aclust <- function(x) hclust(d = x, method = "average")
cclust <- function(x) hclust(d = x, method = "complete")

bestprobeset <- function(eMat, gMat, kdf, kdfg, plotp = FALSE) {
### bestprobeset:  Function to identify probeset associating best with grouping under evaluation,
###                when a single probeset is to be used in a heatmap.
### 
### eMat     Matrix of heatmap values.  Rows:  cases    Cols: probeset
### gMat     Matrix of index labels:  3 columns  1: input gene name 2: probeset  3: annotation gene name
###            obtained from call to refSeqRowMatches_annodf processed into matrix.
### kdf      data frame of METABRIC cases with columns "MBid" and kdfg
### kdfg     name of data frame column with factor variable defining groups to appear in the heatmap
### 
  tapply(seq(nrow(gMat)), gMat[, 1],
         function(x, eM = eMat, gM = gMat, kd = kdf, kg = kdfg) {
           if ( length(x) ) {
             kwout <- rep(0, length(x))
             for ( xi in seq(along = x) ) {
               xind <- x[xi]
               kwdf <- data.frame(ev = eM[ , dimnames(eM)[[2]] == gM[xind, 2]],
                                  gv = kdf[match(kdf$MBid, dimnames(eM)[[1]]), kg])
               kwout[xi] <- kruskal.test(ev ~ gv, data = kwdf)$statistic
               if (plotp) {
                 require("beanplot")
                 beanplot(ev ~ gv, data = kwdf, main = gM[xind, 2])
               }
             }
             return(gM[x, 2][which.max(kwout)])
           } else {
             return(NA_character_)
           }
         } )
}

### Copy number data METABRIC tumour cases
date()
mbcnaddf <- read.delim(file = paste("../",
                        "EGA_Data/EGAS00000000083/EGAD00010000213/",
                        "discovery_CNA_CBS.txt", sep = ""),
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")
date()
mbcnavdf <- read.delim(file = paste("../",
                        "EGA_Data/EGAS00000000083/EGAD00010000215/",
                        "validation_CNA_CBS.txt", sep = ""),
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")
date()
mbcnvddf <- read.delim(file = paste("../",
                        "EGA_Data/EGAS00000000083/EGAD00010000214/",
                        "discovery_CNV_CBS.txt", sep = ""),
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")
date()
mbcnvvdf <- read.delim(file = paste("../",
                        "EGA_Data/EGAS00000000083/EGAD00010000216/",
                        "validation_CNV_CBS.txt", sep = ""),
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")
date()
all.equal(names(mbcnaddf), names(mbcnavdf)) ## TRUE
all.equal(names(mbcnvddf), names(mbcnvvdf)) ## TRUE
mbcnaddf$MBid <- setMBid(mbcnaddf$METABRIC_ID)
mbcnaddf$isDiscovery <- TRUE
mbcnavdf$MBid <- setMBid(mbcnavdf$METABRIC_ID)
mbcnavdf$isDiscovery <- FALSE
mbcnvddf$MBid <- setMBid(mbcnvddf$METABRIC_ID)
mbcnvddf$isDiscovery <- TRUE
mbcnvvdf$MBid <- setMBid(mbcnvvdf$METABRIC_ID)
mbcnvvdf$isDiscovery <- FALSE
nmsvec <- names(mbcnaddf)
nmsvec <- c(nmsvec[grepl("MBid", nmsvec)], nmsvec[!grepl("MBid", nmsvec)])

mbcnadf <- rbind(mbcnaddf[, nmsvec], mbcnavdf[, nmsvec])
mbcnvdf <- rbind(mbcnvddf[, nmsvec], mbcnvvdf[, nmsvec])

NULL
### Expression data METABRIC tumour cases
### Data from European Genome-Phenome Archive from the METABRIC study
### Curtis et al., The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature 2012.
### Note:  About 15 minutes for R to read in 800MB expression data text file.
date()  ##  [1] "Fri Jun 28 16:21:52 2013"
mbexddf <- read.delim(file = paste("../",
                       "EGA_Data/EGAS00000000083/EGAD00010000210/",
                       "discovery_ExpressionMatrix.txt", sep = ""),
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")
date()  ##  [1] "Fri Jun 28 16:37:31 2013"
mbexddf <- mbexddf[order(row.names(mbexddf)), ]
### save(mbexddf, file = paste("../",
###                        "EGA_Data/EGAS00000000083/EGAD00010000210/",
###                        "discovery_ExpressionMatrix.RData", sep = ""))
date()  ##  [1] "Fri Jun 28 16:37:32 2013"

mbexvdf <- read.delim(file = paste("../",
                       "EGA_Data/EGAS00000000083/EGAD00010000211/",
                       "validation_ExpressionMatrix.txt", sep = ""),
                      header = TRUE, stringsAsFactors = FALSE, sep = " ")
date()  ##  [1] "Fri Jun 28 16:40:44 2013"  ## But only 3 mins to read validation.
mbexvdf <- mbexvdf[order(row.names(mbexvdf)), ]
### save(mbexvdf, file = paste("../",
###                        "EGA_Data/EGAS00000000083/EGAD00010000211/",
###                        "validation_ExpressionMatrix.RData", sep = ""))
date()  ##  [1] "Fri Jun 28 16:40:46 2013"

all.equal(row.names(mbexddf), row.names(mbexvdf)) ## [1] TRUE
intersect(names(mbexddf), names(mbexvdf)) ## character(0)

### "Normals"
###
date()  ## [1] "Mon Aug  5 18:22:15 2013"
mbexndf <- read.delim(file = paste("../",
                       "EGA_Data/EGAS00000000083/EGAD00010000212/",
                       "normals_ExpressionMatrix.txt", sep = ""),
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")
date()  ##  [1] "Mon Aug  5 18:23:42 2013"  1.5 mins to read in data
mbexndf <- mbexndf[order(row.names(mbexndf)), ]
###  save(mbexndf, file = paste("../",
###                         "EGA_Data/EGAS00000000083/EGAD00010000212/",
###                         "normals_ExpressionMatrix.RData", sep = ""))

### Outcomes data in cldf below
names(mbexndf) <- setMBid(names(mbexndf))
any(names(mbexndf) %in% names(mbexddf))
any(names(mbexndf) %in% names(mbexvdf))
### any(names(mbexndf) %in% cldf$MBid)
### > any(names(mbexndf) %in% cldf$MBid)
### [1] FALSE
### any(names(mbexndf) %in% setMBid(cldf$match_normal_metabric_id))  ## [1] TRUE
### all(names(mbexndf) %in% setMBid(cldf$match_normal_metabric_id))  ## [1] FALSE
### > sum(names(mbexndf) %in% setMBid(cldf$match_normal_metabric_id))
### [1] 127
### > names(mbexndf)[!(names(mbexndf) %in% setMBid(cldf$match_normal_metabric_id))]
### 
###  [1] "MB.1022" "MB.1063" "MB.0752" "MB.0718" "MB.0938" "MB.0726" "MB.1027"
###  [8] "MB.0706" "MB.1085" "MB.0707" "MB.1018" "MB.0682" "MB.0717" "MB.0782"
### [15] "MB.0669" "MB.0766" "MB.0822"
### 
### So Normals are tissue samples of non-cancerous tissue from breast cancer cases.



### Extract expression data into matrix for faster R processing.
### Discovery set:
mbexdmat <- t(as.matrix(mbexddf))
dtrtcentre <- MBid2Ctr(dimnames(mbexdmat)[[1]])
### Validation set:
mbexvmat <- t(as.matrix(mbexvdf))
vtrtcentre <- MBid2Ctr(dimnames(mbexvmat)[[1]])


sum(is.na(mbexdmat))  ##  [1] 0
sum(is.na(mbexvmat))  ##  [1] 22
vmis <- is.na(mbexvmat)
## Rows with missing data
vmisr <- apply(vmis, 1, function(x) any(x) )
which(vmisr)
length(vmisr)
## Cols with missing data
vmisc <- apply(vmis, 2, function(x) any(x) )
which(vmisc)
length(vmisc)


mbexvmat[vmisr, vmisc]
row.names(mbexvmat[vmisr, vmisc])


### Redo  Zscores from scratch
rzscore.strat.transform <- function(datamat, stratavec)
{
  ## do robust zscore by levels of stratavec
  ##browser()
  ## numGenes <-  ncol(datamat)
  if ( any ( is.na(stratavec) ) ) stop("Cannot handle missing values in strata vec")
  datamatout <- matrix(NA_real_ , nrow = nrow(datamat), ncol = ncol(datamat))
  dimnames(datamatout)[[2]] <- dimnames(datamat)[[2]]
  strataVals <- unique(stratavec)
  numStrata <- length(strataVals)
  for ( si in seq(along = strataVals ) ) {
    idxp <- stratavec == strataVals[si]
    ## Use scale() instead?
    ## datamatout[idxp, ] <-
    ##   sweep(sweep(datamat[idxp, ], 2,
    ##               apply(datamat[idxp, ], 2, function(x) {median(x, na.rm = TRUE)}), "-"), 2,
    ##         apply(datamat[idxp, ], 2, function(x) {diff(quantile(x, probs = c(0.25, 0.75), na.rm = TRUE))}), "/")
    ## A handful of probesets show a high count of a single value, in particular ILMN_1343291 in the validation
    ## set, yielding zero interquartile spread (or underinflated variability).  Get spread of unique values.  
    datamatout[idxp, ] <-
      sweep(sweep(datamat[idxp, ], 2,
                  apply(datamat[idxp, ], 2, function(x) {median(x, na.rm = TRUE)}), "-"), 2,
            apply(datamat[idxp, ], 2, function(x) {diff(quantile(unique(x), probs = c(0.25, 0.75), na.rm = TRUE))}), "/")
    dimnames(datamatout)[[1]][idxp] <- dimnames(datamat)[[1]][idxp]
    
  }
  return(datamatout)    
}

### Create standardized values using robust location/scale
DatasetI.rzs <- rzscore.strat.transform(mbexdmat, dtrtcentre)
DatasetI.rzs <- DatasetI.rzs[, order(dimnames(DatasetI.rzs)[[2]])]
DatasetII.rzs <- rzscore.strat.transform(mbexvmat, vtrtcentre)
DatasetII.rzs <- DatasetII.rzs[, order(dimnames(DatasetII.rzs)[[2]])]
all.equal(dimnames(DatasetI.rzs)[[2]], dimnames(DatasetII.rzs)[[2]]) ## [1] TRUE

Dataset.rzs <- rbind(DatasetI.rzs, DatasetII.rzs)
if ( FALSE ) { rm(DatasetI.rzs, DatasetII.rzs) }
### > all.equal(dimnames(mbexdmat)[[2]], dimnames(mbexvmat)[[2]])   ##  [1] TRUE
mbexdvmat <- rbind(mbexdmat, mbexvmat)
### > all.equal(dimnames(mbexdvmat)[[1]], row.names(sdf))  ##  [1] TRUE
if ( FALSE ) { rm(mbexdmat, mbexvmat) }


Dataset.r <- rbind(mbexdmat, mbexvmat)  ## Note: Same as mbexdvmat


NULL
### > all.equal(dimnames(mbexdmat)[[2]], dimnames(mbexvmat)[[2]])
### [1] TRUE
### > Dataset.r <- rbind(mbexdmat, mbexvmat)
### > all.equal(dimnames(Dataset.rzs), dimnames(Dataset.r))
### [1] TRUE
###

### Match ER binding genes to Illumina IDs

table(ebdf$ensembl_gene_id  %in% annodf$Ensembl_gene_id, useNA = "always")


### Plotting Section -----------------------------------------------------------

### Find genes in annotation data
table(ebdf$ensembl_gene_id  %in% annodf$Ensembl_gene_id, useNA = "always")
table(((ebdf$entrezgene  %in% annodf$Entrez) |
       (ebdf$entrezgene  %in% annodf$Original_Entrez) |
       (ebdf$entrezgene  %in% annodf$Entrez_Gene_ID_0) ),
      useNA = "always")
## No Illumina gene match
ebNoIgdf <- ebdf[!(ebdf$ensembl_gene_id  %in% annodf$Ensembl_gene_id), ]

ebdf$I_Ensembl_gene_id <- match(ebdf$ensembl_gene_id, annodf$Ensembl_gene_id)
is.na(ebdf$I_Ensembl_gene_id[is.na(ebdf$ensembl_gene_id)]) <- TRUE
ebdf$I_Ensembl_gene_id_Gene_symbol <- annodf$Gene_symbol[ebdf$I_Ensembl_gene_id]

ebdf$I_Entrez <- match(ebdf$entrezgene, annodf$Entrez)
is.na(ebdf$I_Entrez[is.na(ebdf$entrezgene)]) <- TRUE
ebdf$I_Entrez_Gene_symbol <- annodf$Gene_symbol[ebdf$I_Entrez]

ebdf$I_Original_Entrez <- match(ebdf$entrezgene, annodf$Original_Entrez)
is.na(ebdf$I_Original_Entrez[is.na(ebdf$entrezgene)]) <- TRUE
ebdf$I_Original_Entrez_Gene_symbol <- annodf$Gene_symbol[ebdf$I_Original_Entrez]

ebdf$I_Entrez_Gene_ID_0 <- match(ebdf$entrezgene, annodf$Entrez_Gene_ID_0)
is.na(ebdf$I_Entrez_Gene_ID_0[is.na(ebdf$entrezgene)]) <- TRUE
ebdf$I_Entrez_Gene_ID_0_Gene_symbol <- annodf$Gene_symbol[ebdf$I_Entrez_Gene_ID_0]

all.equal(ebdf$I_Ensembl_gene_id, ebdf$I_Entrez)
ebdf$NumMatch <-
    sapply(seq(nrow(ebdf)), 
           function(x) {
               allegid <- unlist(ebdf[x, c("I_Ensembl_gene_id", "I_Entrez", "I_Original_Entrez", "I_Entrez_Gene_ID_0")])
               if ( all(is.na(allegid)) ) { 0 } else {
               egid <- allegid[!is.na(allegid)]
               length(unique(egid)) }
               } )
table(ebdf$NumMatch, useNA = "always")
head(ebdf[ebdf$NumMatch == 3, ])

ebdf$NumGenesymbolMatch <-
    sapply(seq(nrow(ebdf)), 
           function(x) {
               allegid <- unlist(ebdf[x, c("I_Ensembl_gene_id_Gene_symbol", "I_Entrez_Gene_symbol",
                                           "I_Original_Entrez_Gene_symbol", "I_Entrez_Gene_ID_0_Gene_symbol")])
               if ( all(is.na(allegid)) ) { 0 } else {
               egid <- allegid[!is.na(allegid)]
               length(unique(egid)) }
               } )
table(ebdf$NumGenesymbolMatch, useNA = "always")

ebdf$Genesymbols <-
    sapply(seq(nrow(ebdf)), 
           function(x) {
               allegid <- unlist(ebdf[x, c("I_Ensembl_gene_id_Gene_symbol", "I_Entrez_Gene_symbol",
                                           "I_Original_Entrez_Gene_symbol", "I_Entrez_Gene_ID_0_Gene_symbol")])
               if ( all(is.na(allegid)) ) { "" } else {
               egid <- allegid[(!(is.na(allegid) | allegid == ""))]
               paste(unique(egid), collapse = " | ") }
               } )

### Mark all annodf rows that match ER binding ChIP-Seq findings:
annodf$ERbinding <- rep(FALSE, nrow(annodf))

### "Ensembl_gene_id"
annodf$ERbinding[which(annodf$Ensembl_gene_id %in% ebdf$ensembl_gene_id[!is.na(ebdf$ensembl_gene_id)])] <- TRUE
### "Entrez" "Original_Entrez" "Entrez_Gene_ID_0"
annodf$ERbinding[which(annodf$Entrez %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)])] <- TRUE
annodf$ERbinding[which(annodf$Original_Entrez %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)])] <- TRUE
annodf$ERbinding[which(annodf$Entrez_Gene_ID_0 %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)])] <- TRUE

### > length(unique(annodf$Ensembl_gene_id))
### [1] 19410
### > length(unique(annodf$Entrez))
### [1] 19267
### > length(unique(annodf$OLMN_Gene_0))
### [1] 0
### > length(unique(annodf$ILMN_Gene_0))
### [1] 37848
### > length(unique(annodf$Symbol_0))
### [1] 25205
### > length(unique(annodf[annodf$ERbinding,]$Ensembl_gene_id))
### [1] 1770
### > length(unique(annodf[!annodf$ERbinding,]$Ensembl_gene_id))
### [1] 17658
### > 17658 + 1770
### [1] 19428
### > 1770/19428*100
### [1] 9.110562
### > # About 9% of all genes show ER binding affinity
### > 

### annodf$ERbinding[ebdf$I_Ensembl_gene_id[!is.na(ebdf$I_Ensembl_gene_id)]] <- TRUE
### annodf$ERbinding[ebdf$I_Entrez[!is.na(ebdf$I_Entrez)]] <- TRUE
### annodf$ERbinding[ebdf$I_Original_Entrez[!is.na(ebdf$I_Original_Entrez)]] <- TRUE
### annodf$ERbinding[ebdf$I_Entrez_Gene_ID_0[!is.na(ebdf$I_Entrez_Gene_ID_0)]] <- TRUE
table(annodf$ERbinding, useNA = "always")
### Now use annodf$ERbinding to annotate graphs, calculate enrichment.

date()
iidxslist <- lapply(prdf$refseqID, function(x) refSeqRowMatches_annodf(x) )
date()
names(iidxslist) <- prdf$GeneName

iidxs <- unlist(iidxslist)
names(iidxs) <- rep(names(iidxslist), lapply(iidxslist, length))
iidxs <- iidxs[!duplicated(iidxs)] ## remove duplicated probeset entries.
iidxlabs <- cbind(names(iidxs), annodf$ProbeId[iidxs], annodf$Gene_symbol[iidxs])
annodf$ProbeId[iidxs]
annodf$Gene_symbol[iidxs]

NULL
### zDHHC23 shows as AK127260
### > annodf$ProbeId[iidxs]
###  [1] "ILMN_1693853" "ILMN_2337058" "ILMN_1761112" "ILMN_2094313" "ILMN_1694514"
###  [6] "ILMN_1774020" "ILMN_1693411" "ILMN_1703370" "ILMN_1684663" "ILMN_1785831"
### [11] "ILMN_2243697" "ILMN_2309926" "ILMN_1737320" "ILMN_2362457" "ILMN_1763568"
### [16] "ILMN_1664946" "ILMN_1697153" "ILMN_1668270" "ILMN_1766896" "ILMN_1769783"
### [21] "ILMN_1654141" "ILMN_1715526" "ILMN_1654606" "ILMN_1736901" "ILMN_1687626"
### [26] "ILMN_2201347" "ILMN_1683544" "ILMN_1679358" "ILMN_1739659" "ILMN_2046003"
### [31] "ILMN_1730568" "ILMN_1789492" "ILMN_2254574" "ILMN_1803824" "ILMN_1748803"
### > annodf$Gene_symbol[iidxs]
###  [1] "HHAT"     "PORCN"    "PORCN"    "ZDHHC1"   "ZDHHC11"  "ZDHHC11" 
###  [7] "ZDHHC11"  "ZDHHC12"  "ZDHHC13"  "ZDHHC13"  "ZDHHC13"  "ZDHHC14" 
### [13] "ZDHHC15"  "ZDHHC16"  "ZDHHC16"  "ZDHHC16"  "ZDHHC17"  "ZDHHC18" 
### [19] "ZDHHC19"  "ZDHHC2"   "ZDHHC20"  "ZDHHC21"  "ZDHHC22"  "AK127260"
### [25] "ZDHHC24"  "ZDHHC3"   "ZDHHC4"   "ZDHHC5"   "ZDHHC6"   "ZDHHC6"  
### [31] "ZDHHC7"   "ZDHHC8"   "ZDHHC9"   "ZDHHC9"   "ZDHHC9"  
### >
### > annodf[iidxs, c("Gene_symbol", "ProbeId", "Gene_synonyms")]
###       Gene_symbol      ProbeId     Gene_synonyms
### 40631        HHAT ILMN_1693853                  
### 20996       PORCN ILMN_2337058 DQ895869 DQ892632
### 38381       PORCN ILMN_1761112                  
### 15526      ZDHHC1 ILMN_2094313                  
### 25554     ZDHHC11 ILMN_1694514                  
### 25553     ZDHHC11 ILMN_1774020                  
### 25555     ZDHHC11 ILMN_1693411                  
### 34650     ZDHHC12 ILMN_1703370          AK092843
### 2626      ZDHHC13 ILMN_1684663                  
### 2627      ZDHHC13 ILMN_1785831                  
### 29809     ZDHHC13 ILMN_2243697                  
### 33600     ZDHHC14 ILMN_2309926                  
### 27174     ZDHHC15 ILMN_1737320                  
### 33020     ZDHHC16 ILMN_2362457                  
### 26609     ZDHHC16 ILMN_1763568                  
### 17026     ZDHHC16 ILMN_1664946                  
### 31748     ZDHHC17 ILMN_1697153                  
### 28027     ZDHHC18 ILMN_1668270                  
### 32211     ZDHHC19 ILMN_1766896                  
### 18565      ZDHHC2 ILMN_1769783                  
### 16107     ZDHHC20 ILMN_1654141                  
### 18019     ZDHHC21 ILMN_1715526                  
### 14910     ZDHHC22 ILMN_1654606                  
### 37693    AK127260 ILMN_1736901           ZDHHC23
### 42213     ZDHHC24 ILMN_1687626          BC002598
### 22115      ZDHHC3 ILMN_2201347                  
### 19852      ZDHHC4 ILMN_1683544                  
### 45647      ZDHHC5 ILMN_1679358          KIAA1748
### 43049      ZDHHC6 ILMN_1739659                  
### 43050      ZDHHC6 ILMN_2046003                  
### 24223      ZDHHC7 ILMN_1730568                  
### 4314       ZDHHC8 ILMN_1789492                  
### 5194       ZDHHC9 ILMN_2254574                  
### 44090      ZDHHC9 ILMN_1803824                  
### 18725      ZDHHC9 ILMN_1748803


### Find corresponding indexes in expression data matrix
izidxs <- rep(NA_integer_, length(iidxs))
for ( i in seq(along = izidxs)) {
  izidxs[i] <- match(annodf$ProbeId[iidxs[i]], dimnames(Dataset.rzs)[[2]])
}
names(izidxs) <- names(iidxs)


bidxslist <- lapply(brcadf$refseqID, function(x) refSeqRowMatches_annodf(x) )
names(bidxslist) <- brcadf$GeneName
bidxslist

bidxs <- unlist(bidxslist)
names(bidxs) <- rep(names(bidxslist), lapply(bidxslist, length))
bidxs <- bidxs[!duplicated(bidxs)] ## Remove duplicate entries
bidxs
bidxlabs <- cbind(names(bidxs), annodf$ProbeId[bidxs], annodf$Gene_symbol[bidxs])
annodf$ProbeId[bidxs]
annodf$Gene_symbol[bidxs]

bzidxs <- rep(NA_integer_, length(bidxs))
for ( i in seq(along = bzidxs)) {
  bzidxs[i] <- match(annodf$ProbeId[bidxs[i]], dimnames(Dataset.rzs)[[2]])
}
bzidxs

require("heatmap.plus")
require("plotrix")
require("RColorBrewer")
require("gplots")

### Clinical data.  T ends up in cldvdf - document how.  It is in these two data files.
###   -rw-r--r--  1 stevenmckinney  staff     123663 Jun 28 17:56 clinical_datasetII_public_FINAL_032112.txt
###   -rw-r--r--  1 stevenmckinney  staff     130604 Jun 28 17:56 clinical_datasetI_public_FINAL_032112.txt

### METABRIC data
### PAM50
### BrCa classification
### cp /Volumes/KilroyHD/kilroy/Projects/MOlab/GavinHa/METABRIC/MaskedRefExpt/DatasetI997/DatasetI997_PAM50classification_ordered.txt .

### brcacldf <- read.delim(file = paste("/Volumes/KilroyHD/kilroy/Projects/MOlab/GavinHa/METABRIC/",
###                          "MaskedRefExpt/DatasetI997/",
###                          "DatasetI997_PAM50classification_ordered.txt", sep = ""),
###                        header = FALSE, stringsAsFactors = FALSE, sep = "\t")
### dim(brcacldf) ## [1] 997   2
### class(brcacldf)
### head(brcacldf)
### names(brcacldf) <- c("mb_id", "BrCaSubtype")
### brcacldf$mbidname <- gsub("-", ".", brcacldf$mb_id)
### table(brcacldf$BrCaSubtype, useNA = "always")
### ### 
### ###  Basal   Her2   LumA   LumB Normal   <NA> 
### ###    118     87    466    268     58      0 
### ###

### brcacldf$BrCaf <- factor(brcacldf$BrCaSubtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal" ))

### cp /Volumes/KilroyHD/kilroy/Projects/MOlab/METABRIC/Data_Clinical/MergeFiles/clinical_datasetI_public_FINAL_032112.txt .

brcaIclalldf <- read.delim(file = paste("../",
                           "Data/",
                           "clinical_datasetI_public_FINAL_032112.txt" , sep = ""),
                       header = TRUE, stringsAsFactors = FALSE, sep = "\t")

brcaIIclalldf <- read.delim(file = paste("../",
                           "Data/",
                           "clinical_datasetII_public_FINAL_032112.txt" , sep = ""),
                       header = TRUE, stringsAsFactors = FALSE, sep = "\t")
clvars <- intersect(names(brcaIclalldf), names(brcaIIclalldf))
cldf <- rbind(brcaIclalldf[, clvars], brcaIIclalldf[, clvars])
cldf$MBid <- setMBid(cldf$METABRIC_ID)




### /Volumes/KilroyHD/kilroy/Projects/MOlab/METABRIC/Data_Clinical/MergeFiles/clinic_all_2010-09-01.txt

clinalldf <- read.delim(file = paste("../",
                           "Data/",
                           "clinic_all_2010-09-01.txt" , sep = ""),
                       header = TRUE, stringsAsFactors = FALSE, sep = "\t")
clinalldf$MBid <- setMBid(clinalldf$metabric_cl_id)

## Get dates from clinalldf
##
## age_at_diagnosis, dob, date_of_diagnosis, dod, last_follow_up_date, last_follow_up_status, first_relapse_type,
## time_until_first_local_regional_relapse, time_until_first_distant_relapse
## histological_type, distant_metastasis_sites_byCarlos, 
##
mvarsclinall <- c("MBid", "age_at_diagnosis", "dob", "date_of_diagnosis", "dod", "last_follow_up_date",
                  "last_follow_up_status", "first_relapse_type",
                  "time_until_first_local_regional_relapse", "time_until_first_distant_relapse",
                  "histological_type", "distant_metastasis_sites_byCarlos")

cldf <- merge(cldf, clinalldf[, mvarsclinall], by = "MBid", suffixes = c("", ".y")  )

head(cldf)


### Cluster analysis identification
clusterIcldf <- read.delim(file = paste("../Data/",
                            "cambridge_clusters_ordered.txt", sep = ""),
                       header = FALSE, stringsAsFactors = FALSE, sep = "\t")
dim(clusterIcldf)
class(clusterIcldf)
head(clusterIcldf)
names(clusterIcldf) <- c("mb_id", "iClusterGroup")
### clusterIcldf$mbidname <- gsub("-", ".", clusterIcldf$mb_id)
table(clusterIcldf$iClusterGroup, useNA = "always")
clusterIcldf$MBid <- setMBid(clusterIcldf$mb_id)
clusterIcldf$Ctr <- MBid2Ctr(clusterIcldf$MBid)

head(clusterIcldf)

### iClusterDatasetIIPredictions.txt
### has intClust calls for validation set.
### 
clusterIIcldf <- read.delim(file = paste("../Data/",
                              "iClusterDatasetIIPredictions.txt", sep = ""),
                            header = TRUE, stringsAsFactors = FALSE, sep = "\t")
names(clusterIIcldf) <- c("mb_id", "iClusterGroup")

clusterIIcldf$MBid <- setMBid(clusterIIcldf$mb_id)
clusterIIcldf$Ctr <- MBid2Ctr(clusterIIcldf$MBid)
head(clusterIIcldf)

clustercldf <- rbind(clusterIcldf, clusterIIcldf)
head(clustercldf)
dim(clustercldf)
clustercldf$iClustf <- factor(clustercldf$iClusterGroup, levels = as.character(1:10))



### Merge clustercldf and brcacldf
### brcacldf <- rbind(brcaIcldf, brcaIIcldf)
mbcldf <- merge(clustercldf, brcacldf)
mbcldf$isDiscovery <- mbcldf$MBid %in% clusterIcldf$MBid
dim(mbcldf)
head(mbcldf)
with(mbcldf, table(iClusterGroup, BrCaf, useNA = "always"))

### Merge BrCa subtype groupings with clinical data
### cldvdf <- merge(mbcldf, cldf, by = "METABRIC_ID", suffixes = c("", ".y"))
cldvdf <- merge(clustercldf, cldf, by = "MBid", suffixes = c("", ".y"))
cldvdf$isDiscovery <- cldvdf$MBid %in% clusterIcldf$MBid
cldvdf$BrCaf <- factor(cldvdf$Pam50Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal" ))

with(cldvdf, table(iClusterGroup, BrCaf, isDiscovery, useNA = "always"))

dim(cldvdf)
names(cldvdf)
### with(cldvdf, all.equal(Pam50Subtype, Pam50Subtype.y))   ## [1] TRUE
cldvdf$survstat <- 1.0*(!(cldvdf$last_follow_up_status == "a"))
cldvdf$survyrs <- cldvdf$T/365.25
with(cldvdf, table(last_follow_up_status, survstat, useNA = "always"))
### cldvdf$iClustf <- factor(cldvdf$iClusterGroup, levels = 1:10, labels = paste("iClust", seq(10)))

cldvdf$gradeb <- 1.0 * (cldvdf$grade == "3")
is.na(cldvdf$gradeb) <- cldvdf$grade == "null"
cldvdf$gradebf <- factor(cldvdf$gradeb, levels = c(0, 1), labels = c("grade 1+2", "grade 3"))

require("survival")
plot(survfit(Surv(survyrs, survstat) ~ BrCaf, data = cldvdf), col = PAM50plotcols)
legend("bottomleft", legend = levels(cldvdf$BrCaf), col = PAM50plotcols, lwd = 3)

plot(survfit(Surv(survyrs, survstat) ~ iClustf, data = cldvdf), col = iClustplotcols)
legend("bottomleft", legend = levels(cldvdf$iClustf), col = iClustplotcols, lwd = 3)

survdiff(Surv(survyrs, survstat) ~ BrCaf, data = cldvdf, rho = 1)  # higher weight to earlier diffs - better for this data
                                        # because of crossing of survival curves at and beyond 7 years.
survdiff(Surv(survyrs, survstat) ~ BrCaf, data = cldvdf, rho = 0)  ## logrank test  rho = 0
survdiff(Surv(survyrs, survstat) ~ iClustf, data = cldvdf, rho = 1)  ## "Breslow" test  rho = 1
### ## "Breslow" test  rho = 1:  Breslow test is close to G-rho test with rho = 1.
survdiff(Surv(survyrs, survstat) ~ iClustf, data = cldvdf, rho = 0)  ## logrank test  rho = 0


cmf <- coxph(Surv(survyrs, survstat) ~ gradebf + iClustf, data = cldvdf,
             subset = !( is.na(gradebf) | is.na(iClustf) | is.na(T) | is.na(survstat) ) )
cmr <- coxph(Surv(survyrs, survstat) ~ iClustf, data = cldvdf,
             subset = !( is.na(gradebf) | is.na(iClustf) | is.na(T) | is.na(survstat) ) )
anova(cmr, cmf)
cmsf <- coxph(Surv(survyrs, survstat) ~ gradebf + strata(iClustf), data = cldvdf,
             subset = !( is.na(gradebf) | is.na(iClustf) | is.na(T) | is.na(survstat) ) )
cmsr <- coxph(Surv(survyrs, survstat) ~ strata(iClustf), data = cldvdf,
             subset = !( is.na(gradebf) | is.na(iClustf) | is.na(T) | is.na(survstat) ) )

anova(cmsf)

### Make directory to hold plots
if ( !file.exists("./Plots") ) dir.create("./Plots")

pdf(file = "./Plots/Surv_iClustf_gradeb12v3_bystrata_v01.pdf", width = 8, height = 10)
par(mfrow = c(4, 3))
nst <- length(levels(cldvdf$iClustf))
for ( si in seq(nst) ) {
  plot(survfit(Surv(survyrs, survstat) ~ gradebf, data = cldvdf,
             subset = (!( is.na(gradebf) | is.na(iClustf) | is.na(T) | is.na(survstat) ) & as.numeric(iClustf) == si) ),
       col = c("green", "red"), main = paste("intClust ", si), xlim = c(0, 8000) )
  ssd <- survdiff(Surv(survyrs, survstat) ~ gradebf, data = cldvdf,
             subset = (!( is.na(gradebf) | is.na(iClustf) | is.na(T) | is.na(survstat) ) & as.numeric(iClustf) == si), rho = 1 )
  ssddf <- (sum(1 * (ssd$exp > 0))) - 1
  ssdpval <- format(signif(1 - pchisq(ssd$chisq, ssddf), 3))
  legend("bottomleft", legend = levels(cldvdf$gradebf), col = c("green", "red"), lwd = 3)
  text(x = 5000, y = 0.95, labels = paste("Breslow p=", ssdpval))
}
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n", main = "", sub = "", xlab = "", ylab = "")
soasd <- survdiff(Surv(survyrs, survstat) ~ gradebf + strata(iClustf), data = cldvdf, rho = 1)
soasddf <- (sum(1 * (soasd$n > 0))) - 1
soasdpval <- format(signif(1 - pchisq(soasd$chisq, soasddf), 3))
text(x = 0.5, y = 0.5, labels = paste("Whole cohort: Breslow p=", soasdpval))
dev.off()

### Datasets I and II
### intClust groups

mbexiClustordidxs <- match(clustercldf$MBid[order(clustercldf$iClusterGroup)], dimnames(Dataset.rzs)[[1]])
mbexBrCaordidxs <- match(brcacldf$MBid[order(brcacldf$BrCaf)], dimnames(Dataset.rzs)[[1]])

iidxslist <- lapply(prdf$refseqID, function(x) refSeqRowMatches_annodf(x) )
names(iidxslist) <- prdf$GeneName
iidxs <- unlist(iidxslist)
names(iidxs) <- rep(names(iidxslist), lapply(iidxslist, length))
iidxlabs <- cbind(names(iidxs), annodf$ProbeId[iidxs], annodf$Gene_symbol[iidxs])

izidxs <- rep(NA_integer_, length(iidxs))
for ( i in seq(along = izidxs)) {
  izidxs[i] <- match(annodf$ProbeId[iidxs[i]], dimnames(Dataset.rzs)[[2]])
}
bidxsbpv <- bestprobeset( eMat = Dataset.rzs, gMat = iidxlabs, kdf = clustercldf, kdfg = "iClustf" )
bestizidxs <- match( bidxsbpv, dimnames(Dataset.rzs)[[2]] )
bestiidxlabs <- iidxlabs[match( bidxsbpv, iidxlabs[, 2] ), ]

require("heatmap.plus")
require("plotrix")
require("RColorBrewer")
require("gplots")

dmat <- t(Dataset.rzs[mbexiClustordidxs, rev(izidxs)])
dmat[dmat < (-3)] <- (-3)
dmat[dmat > 3] <- 3
plotcols <- c("#FF5500", "#00EE76", "#CD3278", "#00C5CD", "#8B0000",
                  "#FFFF40", "#0000CD", "#FFAA00", "#EE82EE", "#7D26CD")
iClustplotcols <- plotcols
pdf(file = "./Plots/dvScreenHitsiClustHeatmap_v01.pdf", width = 8, height = 10)
heatmap(x = dmat, Rowv=NA, Colv=NA, xaxt = "n",
          add.expr, symm = FALSE, revC = FALSE,
          scale = "none", na.rm = TRUE,
          margins = c(5, 5),
          ColSideColors = plotcols[sort(clustercldf$iClusterGroup)],
          labRow = rev(iidxlabs[, 1]),
          labCol = rep("", length(mbexiClustordidxs)), main = NULL,
          xlab = NULL, ylab = NULL,
          keep.dendro = FALSE, verbose = getOption("verbose"),
          col = rev(redgreen(129))  ##rg.colors(16) ##, nc = length(BrCaordidxs)
          )
title(main = "1    2            3                  4               5        6         7             8             9           10                  ", line = -1.5)
color.legend(xl = 0.1, xr = 0.4, yb = 0, yt = 0.02, rect.col = rev(redgreen(129)),
               legend = c("                                          Down   Fold Change    Up"), gradient = "x", align = "rb")

dev.off()

### Use only one probeset per gene target
dbestmat <- t(Dataset.rzs[mbexiClustordidxs, rev(bestizidxs)])
dbestmat[dbestmat < (-3)] <- (-3)
dbestmat[dbestmat > 3] <- 3
iClustdbestmat <- dbestmat

intClustplotcols <- c("#FF5500", "#00EE76", "#CD3278", "#00C5CD", "#8B0000",
                      "#FFFF40", "#0000CD", "#FFAA00", "#EE82EE", "#7D26CD")

pdf(file = "./Plots/dvScreenHitsiClustmaxkwHeatmap_v01.pdf", width = 8, height = 10)
heatmap(x = dbestmat, Rowv=NA, Colv=NA, xaxt = "n",
          add.expr, symm = FALSE, revC = FALSE,
          scale = "none", na.rm = TRUE,
          margins = c(5, 8),
          ColSideColors = intClustplotcols[sort(clustercldf$iClusterGroup)],
          labRow = rev(bestiidxlabs[, 1]),
          labCol = rep("", length(mbexiClustordidxs)), main = NULL,
          xlab = NULL, ylab = NULL,
          keep.dendro = FALSE, verbose = getOption("verbose"),
          col = rev(redgreen(129))  ##rg.colors(16) ##, nc = length(BrCaordidxs)
          )
title(main = "1    2            3                  4             5       6       7             8             9           10                        ", line = -1.5)
require(plotrix)
color.legend(xl = 0.1, xr = 0.4, yb = 0, yt = 0.02, rect.col = rev(redgreen(129)),
               legend = c("                                          Down   Fold Change    Up"), gradient = "x", align = "rb")

dev.off()

### PAM50 breast cancer subtypes

### Datasets I and II
### intClust groups

### mbexiClustordidxs <- match(clustercldf$MBid[order(clustercldf$iClusterGroup)], dimnames(Dataset.rzs)[[1]])
mbexBrCaordidxs <- match(brcacldf$MBid[order(brcacldf$BrCaf)], dimnames(Dataset.rzs)[[1]])

### iidxslist <- lapply(prdf$refseqID, function(x) refSeqRowMatches_annodf(x) )
### names(iidxslist) <- prdf$GeneName
### iidxs <- unlist(iidxslist)
### names(iidxs) <- rep(names(iidxslist), lapply(iidxslist, length))
### iidxlabs <- cbind(names(iidxs), annodf$ProbeId[iidxs], annodf$Gene_symbol[iidxs])
### 
### izidxs <- rep(NA_integer_, length(iidxs))
### for ( i in seq(along = izidxs)) {
###   izidxs[i] <- match(annodf$ProbeId[iidxs[i]], dimnames(Dataset.rzs)[[2]])
### }
p50bidxsbpv <- bestprobeset( eMat = Dataset.rzs, gMat = iidxlabs, kdf = brcacldf, kdfg = "BrCaf" )
p50bestizidxs <- match( p50bidxsbpv, dimnames(Dataset.rzs)[[2]] )
p50bestiidxlabs <- iidxlabs[match( p50bidxsbpv, iidxlabs[, 2] ), ]

dpmat <- t(Dataset.rzs[mbexBrCaordidxs, rev(izidxs)])
dpmat[dpmat < (-3)] <- (-3)
dpmat[dpmat > 3] <- 3
PAM50plotcols <- c("purple", "cyan", "orange", "red", "green")

pdf(file = "./Plots/dvScreenHitsPAM50Heatmap_v01.pdf", width = 8, height = 10)
heatmap(x = dpmat[, 1:1986], Rowv=NA, Colv=NA, xaxt = "n",
          add.expr, symm = FALSE, revC = FALSE,
          scale = "none", na.rm = TRUE,
          margins = c(5, 8),
          ColSideColors = PAM50plotcols[sort(mbcldf[mbcldf$Pam50Subtype != "NC", ]$BrCaf)],
          labRow = rev(iidxlabs[, 1]),
          labCol = rep("", ncol(dpmat[, 1:1986])), main = NULL,
          xlab = NULL, ylab = NULL,
          keep.dendro = FALSE, verbose = getOption("verbose"),
          col = blrdcols
          )
title(main = "   LuminalA                       LuminalB        Her2        Basal          Normal              ",
      line = -1.5, cex = 0.7)
color.legend(xl = 0.1, xr = 0.4, yb = 0, yt = 0.02, rect.col = blrdcols,
               legend = c("                                          Down   Fold Change    Up"), gradient = "x", align = "rb")

dev.off()


### With clustering of genes
### One probeset per gene:

### iidxslist <- lapply(prdf$refseqID, function(x) refSeqRowMatches_annodf(x) )
### names(iidxslist) <- prdf$GeneName
### iidxs <- unlist(iidxslist)
### names(iidxs) <- rep(names(iidxslist), lapply(iidxslist, length))
### iidxlabs <- cbind(names(iidxs), annodf$ProbeId[iidxs], annodf$Gene_symbol[iidxs])
### 
### izidxs <- rep(NA_integer_, length(iidxs))
### for ( i in seq(along = izidxs)) {
###   izidxs[i] <- match(annodf$ProbeId[iidxs[i]], dimnames(Dataset.rzs)[[2]])
### }
mbexBrCaordidxs <- match(brcacldf$MBid[order(brcacldf$BrCaf)], dimnames(Dataset.rzs)[[1]])
p50bidxsbpv <- bestprobeset( eMat = Dataset.rzs, gMat = iidxlabs, kdf = brcacldf, kdfg = "BrCaf" )
p50bestizidxs <- match( p50bidxsbpv, dimnames(Dataset.rzs)[[2]] )
p50bestiidxlabs <- iidxlabs[match( p50bidxsbpv, iidxlabs[, 2] ), ]


dbestmat <- t(Dataset.rzs[mbexBrCaordidxs, rev(p50bestizidxs)])
dbestmat[dbestmat < (-3)] <- (-3)
dbestmat[dbestmat > 3] <- 3
PAM50dbestmat <- dbestmat[, 1:1986]

PAM50plotcols <- c("purple", "cyan", "orange", "red", "green")
blrdcols <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(64))
pdf(file = "./Plots/dvScreenHitsPAM50_maxkw_cclustrow_Heatmap_v01.pdf", width = 8, height = 10)
heatmap(x = dbestmat[, 1:1986], Colv=NA, xaxt = "n", hclustfun = cclust, ## Rowv=NA, 
          symm = FALSE, revC = FALSE,
          scale = "none", na.rm = TRUE,
          margins = c(5, 8),
          ColSideColors = PAM50plotcols[sort(brcacldf[brcacldf$Pam50Subtype != "NC", ]$BrCaf)],
          labRow = rev(p50bestiidxlabs[, 1]),
          labCol = rep("", ncol(dbestmat[, 1:1986])), main = NULL,
          xlab = NULL, ylab = NULL,
          keep.dendro = FALSE, verbose = getOption("verbose"),
          col = blrdcols  
          )
title(main = "                  LuminalA                       LuminalB     Her2     Basal    Normal           ",
      line = -5.2, cex = 0.7)
require(plotrix)
color.legend(xl = 0.2, xr = 0.6, yb = 0.1, yt = 0.12, rect.col = blrdcols, ## rev(redgreen(129)),
               legend = c("                                                      Down        Fold Change         Up"),
             gradient = "x", align = "rb")

dev.off()


### Age by iClust group
cldvdf$ageb40 <- 1.0 * (cldvdf$age_at_diagnosis > 40.0)
cldvdf$ageb40f <- factor(cldvdf$ageb40, labels = c("Dx Age <= 40", "Dx Age > 40"))
with(cldvdf, table(ageb40, ageb40f))

with(cldvdf, table(iClustf, ageb40f))

### Table of intClust by Age at Dx
require("gmodels")
CrossTable(cldvdf$iClustf, cldvdf$ageb40f, digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=FALSE, mcnemar=FALSE,
           resid=FALSE, sresid=TRUE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS")
### Women under 40 are overrepresented in iClust10 and iClust5, and underrepresented in iClust8 and iClust7 and iClust3

date()
save.image()
date()



### Age related

## EZH2 NM_001203247.1  NM_001203248.1  NM_001203249.1  NM_004456.4  NM_152998.2
## MKI67 NM_001145966.1  NM_002417.4
## SKP2 NM_001243120.1  NM_005983.3  NM_032637.3
## YBX1 NM_004559.3  

## List of EZH2 related genes Tomo Osako wanted to focus on.
## This list ends up matching 467 probesets on the Illumina chip.
## Age-Related EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.
## Nov 13 2015:  Naming not optimal - EZH2 associated gene set to assess in age associated investigation is what
##               this gene set is about.  "arGeneSet" file tag denotes analysis focused on this EZH2 H3K27me3
##               relevant gene set compiled by Tomo Osako.
## Some investigations look over all 48803 probesets for age association.
## Other investigations focus on this EZH2-related gene subset.
## Since the paper is about the age-associated nature of EZH2 and H3K27me3 Tomo wanted focus on this gene set.
## Other aspects of this study now include assessing age-associated genes from the whole Illumina probeset list.

ardf <- read.table(file = "../Data/AgeRelated_Accession_Numbers.csv",
                   header = TRUE, sep = ",", stringsAsFactors = FALSE)
ardf <- read.table(file = "../Data/AgeRelated_Accession_Numbersv02.csv",
                   header = TRUE, sep = ",", stringsAsFactors = FALSE)
## ardf: 592 unique RefSeq IDs corresponding to 245 gene / transcript names for EZH2 H3K27me3 relevant gene set.

arplord <- order(toupper(ardf$IHCName))
ardf[arplord, ]
date()
aridxslist <- lapply(ardf$refseqID, function(x) refSeqRowMatches_annodf(x) )
date()
names(aridxslist) <- ardf$GeneName

aridxs <- unlist(aridxslist)
names(aridxs) <- rep(names(aridxslist), lapply(aridxslist, length))
aridxs <- aridxs[!duplicated(aridxs)] ## remove duplicated probeset entries.
aridxlabs <- cbind(names(aridxs), annodf$ProbeId[aridxs], annodf$Gene_symbol[aridxs],
                   ardf[match(names(aridxs), ardf$GeneName), ]$IHCName)
annodf$ProbeId[aridxs]
annodf$Gene_symbol[aridxs]
### > annodf$ProbeId[aridxs]
###  [1] "ILMN_1652913" "ILMN_2364529" "ILMN_1708105" "ILMN_1804654" "ILMN_1734827"
###  [6] "ILMN_1801391" "ILMN_1791002" "ILMN_1665538" "ILMN_1669424" "ILMN_2124769"
### > annodf$Gene_symbol[aridxs]
###  [1] "EZH2"  "EZH2"  "EZH2"  "MKI67" "MKI67" "SKP2"  "SKP2"  "SKP2"  "YBX1" 
### [10] "YBX1" 

### Find corresponding indexes in expression data matrix
arzidxs <- rep(NA_integer_, length(aridxs))
for ( i in seq(along = arzidxs)) {
  arzidxs[i] <- match(annodf$ProbeId[aridxs[i]], dimnames(Dataset.rzs)[[2]])
}
names(arzidxs) <- names(aridxs)

ar.rzs <- Dataset.rzs[, arzidxs]
ar.r <- Dataset.r[, arzidxs]
dim(ar.rzs)
dim(ar.r)
ar.rzs[ar.rzs > 4] <- 4
ar.rzs[ar.rzs < (-4)] <- (-4)

### Calculate AICs at each cut point for each biomarker
### Match to age study unique IDs (no overlap with MB09 TMA and Big Series subset)
### ar.rzs <- ar.rzs[match(asudf$MBid, dimnames(ar.rzs)[[1]]), ]
### ar.r <- ar.r[match(asudf$MBid, dimnames(ar.r)[[1]]), ]
## sardf:  z-score expression data merged with outcomes data
sardf <- cbind(ar.rzs, cldvdf[match(dimnames(ar.rzs)[[1]], cldvdf$MBid), ])
names(sardf)[seq(length(names(arzidxs)))] <-
  paste(ardf[match(names(aridxs), ardf$GeneName), ]$IHCName, names(sardf)[seq(length(names(arzidxs)))], sep = "|")
## sarrdf:  raw expression data merged with outcomes data
sarrdf <- cbind(ar.r, cldvdf[match(dimnames(ar.r)[[1]], cldvdf$MBid), ])
names(sarrdf)[seq(length(names(arzidxs)))] <-
  paste(ardf[match(names(aridxs), ardf$GeneName), ]$IHCName, names(sarrdf)[seq(length(names(arzidxs)))], sep = "|")
###   paste(names(arzidxs), names(sardf)[seq(length(names(arzidxs)))], sep = "|")

## No overlap with MB09 TMA or Big Series
sarnovdf <- sardf[sardf$MBid %in% asudf$MBid, ]
sarrnovdf <- sarrdf[sardf$MBid %in% asudf$MBid, ]




### PAM50 Subtypes Luminal A and Basal
lumaarrdf   <- sarrdf[sarrdf$Pam50Subtype == "LumA", ]
lumbarrdf   <- sarrdf[sarrdf$Pam50Subtype == "LumB", ]
her2arrdf   <- sarrdf[sarrdf$Pam50Subtype == "Her2", ]
basalarrdf  <- sarrdf[sarrdf$Pam50Subtype == "Basal", ]
normalarrdf <- sarrdf[sarrdf$Pam50Subtype == "Normal", ]

lumaarrnovdf   <- sarrnovdf[sarrnovdf$Pam50Subtype == "LumA", ]
lumbarrnovdf   <- sarrnovdf[sarrnovdf$Pam50Subtype == "LumB", ]
her2arrnovdf   <- sarrnovdf[sarrnovdf$Pam50Subtype == "Her2", ]
basalarrnovdf  <- sarrnovdf[sarrnovdf$Pam50Subtype == "Basal", ]
normalarrnovdf <- sarrnovdf[sarrnovdf$Pam50Subtype == "Normal", ]


### End subtypes


### --- Power calcs for ER+ in <=40yrs and 41-49yrs

sdf$yageb40n <- rep(NA_integer_, nrow(sdf))
sdf$yageb40n[sdf$age_at_diagnosis <= 40.9999] <- 1 
sdf$yageb40n[((sdf$age_at_diagnosis > 40.9999) & (sdf$age_at_diagnosis <= 49.9999))] <- 2
sdf$yageb40f <- factor(sdf$yageb40n, labels = c("Age <= 40", "Age 41-49"))
with(sdf, table(yageb40n, yageb40f))
with(sdf, table(yageb40n, yageb40f, ER.Expr))

sdf$BrCaLumn <- rep(NA_integer_, nrow(sdf))
sdf$BrCaLumn[sdf$BrCaf == "LumA" | sdf$BrCaf == "LumB"] <- 1
sdf$BrCaLumn[!(sdf$BrCaf == "LumA" | sdf$BrCaf == "LumB")] <- 2
sdf$BrCaLumf <- factor(sdf$BrCaLumn, labels = c("Luminal", "Other"))
with(sdf, table(yageb40n, yageb40f, BrCaLumf))

### Frequency of the integrative clusters within the Luminal breast cancers vs Other breast cancers
with(sdf, table(iClustf, BrCaLumf))
require("gmodels")
CrossTable(sdf$iClustf, sdf$BrCaLumf, digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=FALSE, mcnemar=FALSE,
           resid=FALSE, sresid=TRUE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS")

### Frequency of the integrative clusters within the Luminal breast cancers
with(sdf[sdf$BrCaLumf == "Luminal", ], table(yageb40f, useNA = "always"))
with(sdf[sdf$BrCaLumf == "Luminal", ], table(age_at_diagnosis, useNA = "always"))
with(sdf[sdf$BrCaLumf == "Luminal", ],
     CrossTable(iClustf, yageb40f, digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=TRUE, mcnemar=FALSE,
           resid=FALSE, sresid=TRUE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS"))
with(sdf[sdf$BrCaLumf == "Luminal", ],
     fisher.test(iClustf, yageb40f, workspace = 20000000))
prop.test(5, 182)
binom.test(5, 182)

### MA12 and MA21 had cases under 60 years old.
with(sdf[sdf$BrCaLumf == "Luminal" & sdf$age_at_diagnosis < 60, ],
     CrossTable(iClustf, yageb40f, digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=TRUE, mcnemar=FALSE,
           resid=FALSE, sresid=TRUE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS"))
with(sdf[sdf$BrCaLumf == "Luminal" & sdf$age_at_diagnosis < 60, ],
     fisher.test(iClustf, yageb40f, workspace = 20000000))

o1vec <- c(rep(0, 1000 - 50), rep(1, 50))
o2vec <- c(rep(0, 1000 - 120), rep(1, 120))

smat <- function(n1=243, n2=283) {
  x1 <- sum(sample(o1vec, n1, replace = TRUE));
  x2 <- sum(sample(o2vec, n2, replace = TRUE));
  matrix(c(x1, n1 - x1, x2, n2 - x2), ncol = 2)
}
set.seed(7537); mean(sapply(seq(10000), function(x) fisher.test(smat(), workspace = 20000000)$p.value < 0.05))
### [1] 0.8136

require("Hmisc")
### intClust2 arm shows about 75% DFS survival at 75 months and 50% survival at 120 months
### intClust5 arm shows about 60% DFS survival at 50 months and 45% survival at  80 months
i2sc <- Weibull2(c(75,120), c(.75,.5))
plot(i2sc(times = 1:150))
i5sc <- Weibull2(c(50,80), c(.6,.45))
plot(i5sc(times = 1:150))

f2 <- Quantile2(i2sc, 
      hratio=function(x)(.1),
      dropin=function(x)(.05*((x^2)/(100 + x^2))),
      dropout=function(x)(.05*((x^2)/(100 + x^2))) )

rcens <- function(n) {6 + (150-6) * (runif(n) ^ .5)}
rcens <- function(n) {60 + (150-60) * (runif(n) ^ .5)}
# To check this, type hist(rcens(10000), nclass=50)
hist(rcens(10000), nclass=50)
### MA12 paper (Annals of Oncology) do not describe censoring patterns.

r2control <- function(n) f2(n, 'control')
r2interv  <- function(n) f2(n, 'intervention')



set.seed(211)
nsim <- 200 ## 100 ## 1000
spowerTable <- rbind(
c(50,
spower(r2control, r2interv, rcens, nc=trunc(.50*50), ni=trunc(.50*50),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*50), ni=trunc(.60*50),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*50), ni=trunc(.70*50),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.25*50), ni=trunc(.75*50),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*50), ni=trunc(.90*50),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*50), ni=trunc(.95*50),test=logrank, nsim=nsim)
),
c(75,
spower(r2control, r2interv, rcens, nc=trunc(.50*75), ni=trunc(.50*75),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*75), ni=trunc(.60*75),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*75), ni=trunc(.70*75),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.20*75), ni=trunc(.80*75),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*75), ni=trunc(.90*75),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*75), ni=trunc(.95*75),test=logrank, nsim=nsim)
),
c(100,
spower(r2control, r2interv, rcens, nc=trunc(.50*100), ni=trunc(.50*100),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*100), ni=trunc(.60*100),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*100), ni=trunc(.70*100),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.20*100), ni=trunc(.80*100),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*100), ni=trunc(.90*100),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*100), ni=trunc(.95*100),test=logrank, nsim=nsim)
),
c(125,
spower(r2control, r2interv, rcens, nc=trunc(.50*125), ni=trunc(.50*125),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*125), ni=trunc(.60*125),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*125), ni=trunc(.70*125),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.20*125), ni=trunc(.80*125),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*125), ni=trunc(.90*125),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*125), ni=trunc(.95*125),test=logrank, nsim=nsim)
),
c(150,
spower(r2control, r2interv, rcens, nc=trunc(.50*150), ni=trunc(.50*150),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*150), ni=trunc(.60*150),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*150), ni=trunc(.70*150),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.20*150), ni=trunc(.80*150),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*150), ni=trunc(.90*150),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*150), ni=trunc(.95*150),test=logrank, nsim=nsim)
),
c(170,
spower(r2control, r2interv, rcens, nc=trunc(.50*170), ni=trunc(.50*170),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*170), ni=trunc(.60*170),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*170), ni=trunc(.70*170),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.20*170), ni=trunc(.80*170),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*170), ni=trunc(.90*170),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*170), ni=trunc(.95*170),test=logrank, nsim=nsim)
),
c(300,
spower(r2control, r2interv, rcens, nc=trunc(.50*300), ni=trunc(.50*300),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*300), ni=trunc(.60*300),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*300), ni=trunc(.70*300),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.20*300), ni=trunc(.80*300),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*300), ni=trunc(.90*300),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*300), ni=trunc(.95*300),test=logrank, nsim=nsim)
),
c(550,
spower(r2control, r2interv, rcens, nc=trunc(.50*550), ni=trunc(.50*550),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*550), ni=trunc(.60*550),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*550), ni=trunc(.70*550),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.25*550), ni=trunc(.75*550),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*550), ni=trunc(.90*550),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*550), ni=trunc(.95*550),test=logrank, nsim=nsim)
),
c(800,
spower(r2control, r2interv, rcens, nc=trunc(.50*800), ni=trunc(.50*800),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.40*800), ni=trunc(.60*800),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.30*800), ni=trunc(.70*800),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.20*800), ni=trunc(.80*800),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.10*800), ni=trunc(.90*800),test=logrank, nsim=nsim),
spower(r2control, r2interv, rcens, nc=trunc(.05*800), ni=trunc(.95*800),test=logrank, nsim=nsim)
)
 )                    


dimnames(spowerTable) <- list(c(), c("Total sample size", "50%/50%", "60%/40%", "70%/30%","80%/20%", "90%/10%","95%/5%"))

write.csv(spowerTable, file = "spowerTable_ywgrant_v01.csv")

### -------------------------

midf <- read.csv(file = "IntClustCentroids_SMmodsv02.csv", stringsAsFactors = FALSE)

mmidf <- merge(midf, annodf, all.x = TRUE, all.y = FALSE, sort = FALSE)
mmidf$GeneTarget <- mmidf$Gene_symbol

write.csv(mmidf[order(mmidf$FeatureNo), ], file = "IntClustCentroidsAnnotated.csv")

###
sdf$age_4cat <- rep(NA_character_, nrow(sdf))
sdf$age_4cat[sdf$age_at_diagnosis <= 40.999] <- "Age <= 40"
sdf$age_4cat[sdf$age_at_diagnosis > 40.999 & sdf$age_at_diagnosis <= 49.999 ] <- "Age 40-49"
sdf$age_4cat[sdf$age_at_diagnosis > 49.999 & sdf$age_at_diagnosis <= 59.999 ] <- "Age 50-59"
sdf$age_4cat[sdf$age_at_diagnosis > 59.999] <- "Age >= 60"
sdf$age_4catf <- factor(sdf$age_4cat, levels = c("Age <= 40", "Age 40-49", "Age 50-59", "Age >= 60"))

CrossTable(sdf$iClustf[sdf$isDiscovery == TRUE], sdf$age_4catf[sdf$isDiscovery == TRUE],
           digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=FALSE, mcnemar=FALSE,
           resid=FALSE, sresid=FALSE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS")

sardf$age_decf <- cut(sardf$age_at_diagnosis, breaks = c(0, 30, 40, 50, 60, 70, 80, 100), right = FALSE)
table(sardf$BrCaf, sardf$age_decf, useNA = "always")


CrossTable(sdf$iClustf[sdf$isDiscovery == FALSE], sdf$age_4catf[sdf$isDiscovery == FALSE],
           digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=FALSE, mcnemar=FALSE,
           resid=FALSE, sresid=FALSE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS")


CrossTable(sdf$iClustf, sdf$age_4catf,
           digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=FALSE, mcnemar=FALSE,
           resid=FALSE, sresid=FALSE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS")


### Correlate transcript and protein - Age related paper 2014
### Transcript - mRNA MBex expression data
### Protein - IHC measures


### MB09 METABRIC Variable Name - protein quantified.  Associate with transcript mRNA expression.
### bcl2_pp_v1n
### birc5_pp_v1n
### cd163_c_v1n
### ck5_c2v1n
### ck56_c2v1n
### ckit_c_v1n
### cyclind1_pp_v1n
### ecad_pp_v1n
### egfr_tsipp_v1n
### er_pp_v1n
### ezh2_pp_v1n
### foxa1_pp_v1n
### h3k27me3_pp_v1n
### her2_v1n
### inpp4b_c_v1.ppn
### ki67_pp_v1n
### nestin_tsipp_v1n
### p16_qsipp_v1n
### tp53_pp_v1n
### pr_c_v1n

### Read in MB09 protein IHC data.   MB09_wholedata.CSV

mtdf <- read.table(file = "../Data/MB09_wholedata.CSV", sep = ",", quote='"', header = TRUE, stringsAsFactors = FALSE)


### MB09_VarName,Illumina ProbeID
### bcl2_pp_v1n,ILMN_1801119
### birc5_pp_v1n,ILMN_1803124
### cd163_c_v1n,ILMN_1733270
### ck5_c2v1n,ILMN_1733270
### ck56_c2v1n,NA
### ckit_c_v1n,ILMN_2229379
### cyclind1_pp_v1n,ILMN_1688480
### ecad_pp_v1n,ILMN_1770940
### egfr_tsipp_v1n,ILMN_1798975
### er_pp_v1n,ILMN_1678535
### ezh2_pp_v1n,ILMN_2364529
### foxa1_pp_v1n,ILMN_1766650
### h3k27me3_pp_v1n,NA
### her2_v1n,ILMN_2352131
### inpp4b_c_v1.ppn,ILMN_2198878
### ki67_pp_v1n,ILMN_1734827
### nestin_tsipp_v1n,ILMN_1738147
### p16_qsipp_v1n,ILMN_1717714
### tp53_pp_v1n,ILMN_1779356
### pr_c_v1n,ILMN_1811014
### 
### 
###

ehdf <- read.table(file = "../Data/EZH2relatedGenesv08.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

date()
ehidxslist <- lapply(ehdf$refseqID, function(x) refSeqRowMatches_annodf(x) )
date()
names(ehidxslist) <- ehdf$GeneName

ehidxs <- unlist(ehidxslist)
names(ehidxs) <- rep(names(ehidxslist), lapply(ehidxslist, length))
ehidxs <- ehidxs[!duplicated(ehidxs)] ## remove duplicated probeset entries.
ehidxlabs <- cbind(names(ehidxs), annodf$ProbeId[ehidxs], annodf$Gene_symbol[ehidxs],
                   ehdf[match(names(ehidxs), ehdf$GeneName), ]$IHCName)
annodf$ProbeId[ehidxs]
annodf$Gene_symbol[ehidxs]


### Find corresponding indexes in expression data matrix
ehzidxs <- rep(NA_integer_, length(ehidxs))
for ( i in seq(along = ehzidxs)) {
  ehzidxs[i] <- match(annodf$ProbeId[ehidxs[i]], dimnames(Dataset.rzs)[[2]])
}
names(ehzidxs) <- names(ehidxs)


eh.rzs <- Dataset.rzs[, ehzidxs]
eh.r <- Dataset.r[, ehzidxs]
dim(eh.rzs)
dim(eh.r)
eh.rzs[eh.rzs > 4] <- 4
eh.rzs[eh.rzs < (-4)] <- (-4)

### Split cases at 60 years age at diagnosis.

### Trend plots first.


### Whole cohort

sehdf <- cbind(eh.rzs, cldvdf[match(dimnames(eh.rzs)[[1]], cldvdf$MBid), ])
sehrdf <- cbind(eh.r, cldvdf[match(dimnames(eh.r)[[1]], cldvdf$MBid), ])
names(sehdf)[seq(length(names(ehzidxs)))] <-
  paste(ehdf[match(names(ehidxs), ehdf$GeneName), ]$IHCName, names(sehdf)[seq(length(names(ehzidxs)))], sep = "|")
names(sehrdf)[seq(length(names(ehzidxs)))] <-
  paste(ehdf[match(names(ehidxs), ehdf$GeneName), ]$IHCName, names(sehrdf)[seq(length(names(ehzidxs)))], sep = "|r|")

### Regression for age <= 60, age > 60   i <- 10; j <- 4
ehwclfitageLE60 <- vector("list", length(names(ehzidxs)))
ehwclfit <- vector("list", length(names(ehzidxs)))
stj <- "All cases"
pdf(file = "./Plots/EZH2_All_AgeRelated_lm_v01.pdf", width = 8, height = 10)
xlims <- c(20, 100)
ylims <- c(-5, 6)
par(mfrow = c(2, 2))
ageLE60idxp <- (sehdf$age_at_diagnosis <= 60)
### for ( i in seq(along = names(ehzidxs)) ) {
for ( ni in seq(along = names(ehzidxs)) ) {
    i <- order(toupper(names(sehdf)[seq(length(names(ehzidxs)))]))[ni]
    cat("\n\n### --- ", names(sehdf)[i])
    lmifitageLE60 <- lm(sehdf[ageLE60idxp, i] ~ sehdf[ageLE60idxp, "age_at_diagnosis"])
    smrylmifitageLE60 <- summary(lmifitageLE60)
    ehwclfitageLE60[[i]] <- lmifitageLE60
    print(smrylmifitageLE60)
    lmifit <- lm(sehdf[, i] ~ sehdf[, "age_at_diagnosis"])
    smrylmifit <- summary(lmifit)
    ehwclfit[[i]] <- lmifit
    print(smrylmifit)
    plot(sehdf[, "age_at_diagnosis"], sehdf[, i],
         ylab = "Normalized Expression (-4, 4)", xlab = "Age at diagnosis",
         main = paste(stj, names(sehdf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
    lines(supsmu(sehdf[, "age_at_diagnosis"], sehdf[, i], span = 0.4, bass = 10), lwd = 3)
    abline(h = 0, v = 60, lty = 2)
    text(x = 40, y = 5.6, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""))
    text(x = 70, y = 5.0, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""))
  }
dev.off()

### Raw data
### Regression for age <= 60, age > 60   i <- 10; j <- 4
ehrwclfitageLE60 <- vector("list", length(names(ehzidxs)))
ehrwclfit <- vector("list", length(names(ehzidxs)))
stj <- "All cases"
pdf(file = "./Plots/EZH2_All_AgeRelated_lm_raw_v01.pdf", width = 8, height = 10)
xlims <- c(20, 100)
ylims <- c(4, 17)
par(mfrow = c(2, 2))
ageLE60idxp <- (sehrdf$age_at_diagnosis <= 60)
### for ( i in seq(along = names(ehzidxs)) ) {
for ( ni in seq(along = names(ehzidxs)) ) {
    i <- order(toupper(names(sehrdf)[seq(length(names(ehzidxs)))]))[ni]
    cat("\n\n### --- ", names(sehrdf)[i])
    lmifitageLE60 <- lm(sehrdf[ageLE60idxp, i] ~ sehrdf[ageLE60idxp, "age_at_diagnosis"])
    smrylmifitageLE60 <- summary(lmifitageLE60)
    ehrwclfitageLE60[[i]] <- lmifitageLE60
    print(smrylmifitageLE60)
    lmifit <- lm(sehrdf[, i] ~ sehrdf[, "age_at_diagnosis"])
    smrylmifit <- summary(lmifit)
    ehrwclfit[[i]] <- lmifit
    print(smrylmifit)
    plot(sehrdf[, "age_at_diagnosis"], sehrdf[, i],
         ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
         main = paste(stj, names(sehrdf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
    lines(supsmu(sehrdf[, "age_at_diagnosis"], sehrdf[, i], span = 0.4, bass = 10), lwd = 3)
    abline(h = 0, v = 60, lty = 2)
    text(x = 40, y = 16.3, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""))
    text(x = 70, y = 15.5, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""))
  }
dev.off()
### Slope too small  -0.0037063 -0.001621
### Borderline 0.008203
### Slope ok   0.008711  0.027910 0.008782 0.014513 -0.012963 0.009837
### p-value a function of slope and data spread.
### > log2(1.5)/70
### [1] 0.008356607  So a slope of 0.0084 or larger corresponds to a fold change of 1.5 or greater.
### > log2(1.5)/35
### [1] 0.01671321  So a slope of 0.02 or larger for <60 or >60 corresponds to fold change of 1.5 or greater.
### AR raw = 0.008711 0.027910  z = 0.008095  0.013545
### NOTCH1 raw =  -0.0072142 -0.011460  z = -0.008844 -0.010421

### Biomarker shows age-dependent trend if
### ( ( abs(Slope for age < 60) > log2(1.25)/35 or abs(Slope for age > 60) > log2(1.25)/35 or abs(Slope for age) > log2(1.25)/70 ) AND
###   (p-val < 0.05/(2 * num biomarkers) ) )
### Save
###  3 slopes:,
###  3 pvalues:,
###  ILMN probeset ID: Probe_id, Probe_sequence,
###  gene name: Gene_symbol, Gene_synonyms, Synonyms_0, ILMN_Gene_0, Chromosome_0, Cytoband_0, Original_genomic_annotation
###  refseq nm_ code:  RefSeq_transcripts, RefSeq_ID_0

aroutdf <- annodf
aroutmatcolnames <- c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
aroutmat <- matrix(NA_real_, nrow = nrow(aroutdf), ncol = length(aroutmatcolnames))
LE60_slope_col <- match("LE60_slope", aroutmatcolnames)
GT60_slope_col <- match("GT60_slope", aroutmatcolnames)
AllAges_slope_col <- match("AllAges_slope", aroutmatcolnames)
LE60_pval_col <- match("LE60_pval", aroutmatcolnames)
GT60_pval_col <- match("GT60_pval", aroutmatcolnames)
AllAges_pval_col <- match("AllAges_pval", aroutmatcolnames)
AgeDependentp_col <- match("AgeDependentp", aroutmatcolnames)
aroutdf$LE60_slope <- rep(NA_real_ , nrow(aroutdf))
aroutdf$GT60_slope <- rep(NA_real_ , nrow(aroutdf))
aroutdf$AllAges_slope <- rep(NA_real_ , nrow(aroutdf))
aroutdf$LE60_pval <- rep(NA_real_ , nrow(aroutdf))
aroutdf$GT60_pval <- rep(NA_real_ , nrow(aroutdf))
aroutdf$AllAges_pval <- rep(NA_real_ , nrow(aroutdf))
aroutdf$AgeDependentp <- rep(NA , nrow(aroutdf))
biosigslopeLE60 <-  log2(1.25)/35
biosigslopeGT60 <-  log2(1.25)/35
biosigslopeAllAges <-  log2(1.25)/70
multcompPval <- 0.05/(2 * nrow(aroutdf))
arlmdf <- data.frame(age_at_diagnosis = sehrdf$age_at_diagnosis,
                       probesetni = Dataset.r[, 1])

cat("\n\n\n")
for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
  psi <- aroutdf$ProbeId[ni]
  drci <- match(psi, dimnames(Dataset.r)[[2]])
  arlmdf$probesetni <- Dataset.r[, drci]
  if (ni %% 100 == 0 ) { cat(ni, ", ") }
  lmifitageLE60 <- lm(probesetni ~ age_at_diagnosis, data = arlmdf, na.action = na.exclude, subset = age_at_diagnosis <= 60)
  lmifitageGT60 <- lm(probesetni ~ age_at_diagnosis, data = arlmdf, na.action = na.exclude, subset = age_at_diagnosis > 60)
  lmifitallages <- lm(probesetni ~ age_at_diagnosis, data = arlmdf, na.action = na.exclude)
  smrylmifitageLE60 <- summary(lmifitageLE60)
  smrylmifitageGT60 <- summary(lmifitageGT60)
  smrylmifitallages <- summary(lmifitallages)

  aroutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
  aroutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
  aroutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
  aroutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
  aroutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
  aroutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
  aroutmat[ni, AgeDependentp_col] <- ( ( ( abs(aroutmat[ni, LE60_slope_col])    > biosigslopeLE60 ) ||
                                   ( abs(aroutmat[ni, GT60_slope_col])    > biosigslopeGT60 ) ||
                                   ( abs(aroutmat[ni, AllAges_slope_col]) > biosigslopeAllAges ) ) &&
                                 ( ( aroutmat[ni, LE60_pval_col]    < multcompPval ) ||
                                   ( aroutmat[ni, GT60_pval_col]    < multcompPval ) ||
                                   ( aroutmat[ni, AllAges_pval_col] < multcompPval ) ) )
}
dimnames(aroutmat)[[2]] <- aroutmatcolnames
dimnames(aroutmat)[[1]] <- aroutdf$ProbeId
aroutdf[, aroutmatcolnames] <- aroutmat
cat("\n\n\n")
### Adjust all p-values and redo age-dependent selection on adjusted p-values < 0.05
aroutdf$LE60_BHadj_pval <- p.adjust(aroutdf$LE60_pval, method="BH")
aroutdf$GT60_BHadj_pval <- p.adjust(aroutdf$GT60_pval, method="BH")
aroutdf$AllAges_BHadj_pval <- p.adjust(aroutdf$AllAges_pval, method="BH")
aroutdf$BHadj_AgeDependentp <- ( ( ( abs(aroutdf$LE60_slope) > biosigslopeLE60 ) |
                                   ( abs(aroutdf$GT60_slope) > biosigslopeGT60 ) |
                                   ( abs(aroutdf$AllAges_slope) > biosigslopeAllAges ) ) &
                                 ( ( aroutdf$LE60_BHadj_pval < 0.05 ) |
                                   ( aroutdf$GT60_BHadj_pval < 0.05 ) |
                                   ( aroutdf$AllAges_BHadj_pval < 0.05 ) ) )
### Alternative definition to consider:  Use this version.
aroutdf$BHadj_and_AgeDependentp <-
  ( ( ( abs(aroutdf$LE60_slope) > biosigslopeLE60 )  &  ( aroutdf$LE60_BHadj_pval < 0.05 ) ) |
    ( ( abs(aroutdf$GT60_slope) > biosigslopeGT60 )  &  ( aroutdf$GT60_BHadj_pval < 0.05 ) ) |
    ( ( abs(aroutdf$AllAges_slope) > biosigslopeAllAges ) &  ( aroutdf$AllAges_BHadj_pval < 0.05 ) ) )

aroutdf$BHadj_signifp <- ( ( ( aroutdf$LE60_BHadj_pval < 0.05 ) |
                             ( aroutdf$GT60_BHadj_pval < 0.05 ) |
                             ( aroutdf$AllAges_BHadj_pval < 0.05 ) ) )


aroutdf$LE60_log2FC <- aroutdf$LE60_slope * 35
aroutdf$GT60_log2FC <- aroutdf$GT60_slope * 35
aroutdf$AllAges_log2FC <- aroutdf$AllAges_slope * 70

aroutdf$Best_log2FC <- aroutdf$AllAges_log2FC
aroutdf$Best_BHadj_pval <- aroutdf$AllAges_BHadj_pval

## Find the largest significant fold change:

LE60_Bestp <- ( ( abs(aroutdf$LE60_slope) > biosigslopeLE60 )          & ( aroutdf$LE60_BHadj_pval < 0.05 ) )
aroutdf[LE60_Bestp, ]$Best_log2FC <- aroutdf[LE60_Bestp, ]$LE60_log2FC
aroutdf[LE60_Bestp, ]$Best_BHadj_pval <- aroutdf[LE60_Bestp, ]$LE60_BHadj_pval

GT60_Bestp <- ( ( abs(aroutdf$GT60_slope) > biosigslopeGT60 )          & ( aroutdf$GT60_BHadj_pval < 0.05 ) &
                ( abs(aroutdf$GT60_slope) > abs(aroutdf$LE60_slope) ) )
aroutdf[GT60_Bestp, ]$Best_log2FC <- aroutdf[GT60_Bestp, ]$GT60_log2FC
aroutdf[GT60_Bestp, ]$Best_BHadj_pval <- aroutdf[GT60_Bestp, ]$GT60_BHadj_pval

AllAges_Bestp <- ( ( abs(aroutdf$AllAges_slope) > biosigslopeAllAges ) & ( aroutdf$AllAges_BHadj_pval < 0.05 ) &
                   ( ( abs(aroutdf$AllAges_log2FC) > abs(aroutdf$LE60_log2FC) ) |
                     ( abs(aroutdf$AllAges_log2FC) > abs(aroutdf$GT60_log2FC) ) ) )
aroutdf[AllAges_Bestp, ]$Best_log2FC <- aroutdf[AllAges_Bestp, ]$AllAges_log2FC
aroutdf[AllAges_Bestp, ]$Best_BHadj_pval <- aroutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
NULL
### Use this Best_BHadj_pval for manhattan plots as well.

### Add ER binding info to aroutdf
aroutdf$ERbinding <- rep(FALSE, nrow(aroutdf))
### "Ensembl_gene_id"
aroutdf$ERbinding[which(aroutdf$Ensembl_gene_id %in% ebdf$ensembl_gene_id[!is.na(ebdf$ensembl_gene_id)])] <- TRUE
### "Entrez" "Original_Entrez" "Entrez_Gene_ID_0"
aroutdf$ERbinding[which(aroutdf$Entrez %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)])] <- TRUE
aroutdf$ERbinding[which(aroutdf$Original_Entrez %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)])] <- TRUE
aroutdf$ERbinding[which(aroutdf$Entrez_Gene_ID_0 %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)])] <- TRUE

length(unique(aroutdf[aroutdf$ERbinding,]$Ensembl_gene_id))
length(unique(aroutdf[aroutdf$ERbinding & aroutdf$BHadj_and_AgeDependentp,]$Ensembl_gene_id)) ## 537
length(unique(aroutdf[aroutdf$BHadj_and_AgeDependentp,]$Ensembl_gene_id)) ## 3897
fisher.test(matrix(c(537, 1770 - 537, 3897 - 537, 19428 - 3897 - (1770 - 537)), nrow = 2, byrow = TRUE))

### > length(unique(aroutdf[aroutdf$ERbinding & aroutdf$BHadj_and_AgeDependentp,]$Ensembl_gene_id))
### [1] 537
### > length(unique(aroutdf[aroutdf$BHadj_and_AgeDependentp,]$Ensembl_gene_id)) ## 
### [1] 3897
### > 537/3897
### [1] 0.1377983
### > 
### > require("gmodels")
### Loading required package: gmodels
### > matrix(c(537, 1770 - 537, 3897 - 537, 19428 - 3897 - (1770 - 537)), nrow = 2, byrow = TRUE)
###      [,1]  [,2]
### [1,]  537  1233
### [2,] 3360 14298
### > sum(c(537, 1770 - 537, 3897 - 537, 19428 - 3897 - (1770 - 537)))
### [1] 19428
### > fisher.test(matrix(c(537, 1770 - 537, 3897 - 537, 19428 - 3897 - (1770 - 537)), nrow = 2, byrow = TRUE))
### 
### 	Fisher's Exact Test for Count Data
### 
### data:  
### p-value < 2.2e-16
### alternative hypothesis: true odds ratio is not equal to 1
### 95 percent confidence interval:
###  1.660123 2.066952
### sample estimates:
### odds ratio 
###   1.853342 
### 
### So about 13% of age associated genes are ER binding so ER affinity is enriched in the age associated set.
require("gmodels")
CrossTable(matrix(c(537, 1770 - 537, 3897 - 537, 19428 - 3897 - (1770 - 537)), nrow = 2, byrow = TRUE,
                   dimnames = list("ER_binding" = c("ER binding", "Not ER binding"),
                                   "Age association" = c("Age associated", "Not age associated"))),
           prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE, fisher = TRUE, format = "SPSS")


with(aroutdf, table(BHadj_AgeDependentp, BHadj_and_AgeDependentp, useNA="ifany"))

aroutdf[aroutdf$BHadj_AgeDependentp & !aroutdf$BHadj_and_AgeDependentp, -c(105:107)]

### Work up a cluster implementation of above process.  About 22 hours to do the whole loop. 10 mins with a matrix instead of df.
write.csv(aroutdf, file = "aroutdf.csv")
write.csv(aroutdf[aroutdf$AgeDependentp == 1, ], file = "AgeDependent_aroutdf.csv")

### AgeDependent_BHadj_FC1.25_aroutdf_ProbeIds <- arnovoutdf[aroutdf$BHadj_AgeDependentp, "Probe_id"] same as
AgeDependent_BHadj_FC1.25_aroutdf_ProbeIds <- aroutdf[aroutdf$BHadj_AgeDependentp, "Probe_id"]

write.csv(aroutdf[aroutdf$BHadj_AgeDependentp, ], file = "AgeDependent_BHadj_aroutdf.csv")
write.csv(aroutdf[aroutdf$BHadj_and_AgeDependentp, ], file = "AgeDependent_BHadj_and_FC1.25_aroutdf.csv")

### Plot expression data by age for all probes in the Age-Related EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.

### Whole cohort 1992 cases

### Regression for age <= 60, age > 60   i <- 10; j <- 4
wclfitageLE60 <- vector("list", length(names(arzidxs)))
wclfit <- vector("list", length(names(arzidxs)))
stj <- paste(dim(sardf)[1], "cases")

pdf(file = "./Plots/AgeRelated_All1992Cases_lm_v10.pdf", width = 8, height = 10)
AgeDependent_and_BHadj_FC1.25_ProbeIds <- aroutdf[aroutdf$BHadj_and_AgeDependentp, "Probe_id"]

xlims <- c(20, 100)
ylims <- c(-5, 6)
par(mfrow = c(2, 2))
ageLE60idxp <- (sardf$age_at_diagnosis <= 60)
for ( ni in seq( along = names(arzidxs) ) ) {
  i <- order(toupper(names(sardf)[seq(length(names(arzidxs)))]))[ni]
  cat("\n\n### --- ", names(sardf)[i])
  lmifitageLE60 <- lm(sardf[ageLE60idxp, i] ~ sardf[ageLE60idxp, "age_at_diagnosis"])
  smrylmifitageLE60 <- summary(lmifitageLE60)
  wclfitageLE60[[i]] <- lmifitageLE60
  print(smrylmifitageLE60)
  lmifit <- lm(sardf[, i] ~ sardf[, "age_at_diagnosis"])
  smrylmifit <- summary(lmifit)
  wclfit[[i]] <- lmifit
  print(smrylmifit)
  plot(sardf[, "age_at_diagnosis"], sardf[, i],
       ylab = "Normalized Expression (-4, 4)", xlab = "Age at diagnosis",
       main = paste(stj, names(sardf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
  lines(supsmu(sardf[, "age_at_diagnosis"], sardf[, i], span = 0.4, bass = 10), lwd = 3)
  abline(h = 0, v = 50, lty = 2)
  text(x = 20, y = 5.6, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  text(x = 22, y = 5.0, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  if ( strsplit(names(sardf)[i], split = "\\|")[[1]][2] %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
    text(x = 98, y = 5.6, labels="Age-dependent trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 5.0, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
  } else {
    text(x = 98, y = 5.6, labels="No detectable trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 5.0, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
  }
}
dev.off()

### Whole cohort - no overlap with MB09 or Big Series TMA

### Regression for age <= 60, age > 60   i <- 10; j <- 4
wclfitageLE60 <- vector("list", length(names(arzidxs)))
wclfit <- vector("list", length(names(arzidxs)))
stj <- paste(dim(sarnovdf)[1], "cases")

pdf(file = "./Plots/AgeRelated_NoOverlapMB09BigSeries_lm_v10.pdf", width = 8, height = 10)
AgeDependent_and_BHadj_FC1.25_ProbeIds <- arnovoutdf[arnovoutdf$BHadj_and_AgeDependentp, "Probe_id"]

xlims <- c(20, 100)
ylims <- c(-5, 6)
par(mfrow = c(2, 2))
ageLE60idxp <- (sarnovdf$age_at_diagnosis <= 60)
for ( ni in seq( along = names(arzidxs) ) ) {
  i <- order(toupper(names(sarnovdf)[seq(length(names(arzidxs)))]))[ni]
  cat("\n\n### --- ", names(sarnovdf)[i])
  lmifitageLE60 <- lm(sarnovdf[ageLE60idxp, i] ~ sarnovdf[ageLE60idxp, "age_at_diagnosis"])
  smrylmifitageLE60 <- summary(lmifitageLE60)
  wclfitageLE60[[i]] <- lmifitageLE60
  print(smrylmifitageLE60)
  lmifit <- lm(sarnovdf[, i] ~ sarnovdf[, "age_at_diagnosis"])
  smrylmifit <- summary(lmifit)
  wclfit[[i]] <- lmifit
  print(smrylmifit)
  plot(sarnovdf[, "age_at_diagnosis"], sarnovdf[, i],
       ylab = "Normalized Expression (-4, 4)", xlab = "Age at diagnosis",
       main = paste(stj, names(sarnovdf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
  lines(supsmu(sarnovdf[, "age_at_diagnosis"], sarnovdf[, i], span = 0.4, bass = 10), lwd = 3)
  abline(h = 0, v = 50, lty = 2)
  text(x = 20, y = 5.6, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  text(x = 22, y = 5.0, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  if ( strsplit(names(sarnovdf)[i], split = "\\|")[[1]][2] %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
    text(x = 98, y = 5.6, labels="Age-dependent trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 5.0, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
  } else {
    text(x = 98, y = 5.6, labels="No detectable trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 5.0, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
  }
}
dev.off()

### Raw data  Whole cohort 1992 cases

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
aroutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(aroutdf))
aroutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(aroutdf))
aroutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(aroutdf))
aroutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(aroutdf))

aroutdf$arGeneSetidxp <- ( aroutdf$ProbeId %in% annodf$ProbeId[aridxs] )
  
aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
  p.adjust(aroutdf[aroutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
  p.adjust(aroutdf[aroutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
  p.adjust(aroutdf[aroutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")

aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
  ( ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
     ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < 0.05 ) ) |
    ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
     ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < 0.05 ) ) |
    ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
     ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < 0.05 ) ) )



### Regression for age <= 60, age > 60   i <- 10; j <- 4
wclfitageLE60 <- vector("list", length(names(arzidxs)))
wclfit <- vector("list", length(names(arzidxs)))
stj <- paste(dim(sarrdf)[1], "cases")

pdf(file = "./Plots/AgeRelated_raw_All1992Cases_lm_v11.pdf", width = 8, height = 10)
AgeDependent_and_BHadj_FC1.25_ProbeIds <- aroutdf[which(aroutdf$arGeneSet_BHadj_and_AgeDependentp), "Probe_id"]

xlims <- c(20, 100)
ylims <- c(4, 17)
par(mfrow = c(2, 2))
ageLE60idxp <- (sarrdf$age_at_diagnosis <= 60)

for ( ni in seq( along = names(arzidxs) ) ) {
  i <- order(toupper(names(sarrdf)[seq(length(names(arzidxs)))]))[ni]
  cat("\n\n### --- ", names(sarrdf)[i])
  lmifitageLE60 <- lm(sarrdf[ageLE60idxp, i] ~ sarrdf[ageLE60idxp, "age_at_diagnosis"])
  smrylmifitageLE60 <- summary(lmifitageLE60)
  wclfitageLE60[[i]] <- lmifitageLE60
  print(smrylmifitageLE60)
  lmifit <- lm(sarrdf[, i] ~ sarrdf[, "age_at_diagnosis"])
  smrylmifit <- summary(lmifit)
  wclfit[[i]] <- lmifit
  print(smrylmifit)
  plot(sarrdf[, "age_at_diagnosis"], sarrdf[, i],
       ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
       main = paste(stj, names(sarrdf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
  lines(supsmu(sarrdf[, "age_at_diagnosis"], sarrdf[, i], span = 0.4, bass = 10), lwd = 3)
  abline(h = 0, v = 50, lty = 2)
  text(x = 20, y = 16.3, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  text(x = 22, y = 15.5, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  if ( strsplit(names(sarrdf)[i], split = "\\|")[[1]][2] %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
    text(x = 98, y = 16.3, labels="Age-dependent trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 15.5, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
  } else {
    text(x = 98, y = 16.3, labels="No detectable trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 15.5, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
  }
}
dev.off()

write.csv(aroutdf[which(aroutdf$arGeneSetidxp), -c(105:107)],
          file = "AgeRelated_arGeneSet_All1992Cases_lm_v11.csv")

### Raw data   - no overlap with MB09 or Big Series TMA


### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
arnovoutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(arnovoutdf))
arnovoutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(arnovoutdf))
arnovoutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(arnovoutdf))
arnovoutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(arnovoutdf))

arnovoutdf$arGeneSetidxp <- ( arnovoutdf$ProbeId %in% annodf$ProbeId[aridxs] )
  
arnovoutdf[arnovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
  p.adjust(arnovoutdf[arnovoutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
arnovoutdf[arnovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
  p.adjust(arnovoutdf[arnovoutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
arnovoutdf[arnovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
  p.adjust(arnovoutdf[arnovoutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")

arnovoutdf[arnovoutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
  ( ( ( abs(arnovoutdf[arnovoutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
     ( arnovoutdf[arnovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < 0.05 ) ) |
    ( ( abs(arnovoutdf[arnovoutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
     ( arnovoutdf[arnovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < 0.05 ) ) |
    ( ( abs(arnovoutdf[arnovoutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
     ( arnovoutdf[arnovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < 0.05 ) ) )


### Regression for age <= 60, age > 60   i <- 10; j <- 4
wclfitageLE60 <- vector("list", length(names(arzidxs)))
wclfit <- vector("list", length(names(arzidxs)))
stj <- paste(dim(sarrnovdf)[1], "cases")

pdf(file = "./Plots/AgeRelated_raw_NoOverlapMB09BigSeries_lm_v11.pdf", width = 8, height = 10)
AgeDependent_and_BHadj_FC1.25_ProbeIds <- arnovoutdf[which(arnovoutdf$arGeneSet_BHadj_and_AgeDependentp), "Probe_id"]

xlims <- c(20, 100)
ylims <- c(4, 17)
par(mfrow = c(2, 2))
ageLE60idxp <- (sarrnovdf$age_at_diagnosis <= 60)
## AgeDependent_BHadj_FC1.25_arnovoutdf_ProbeIds <- arnovoutdf[arnovoutdf$BHadj_AgeDependentp, "Probe_id"] below
### for ( i in seq(along = names(arzidxs)) ) {
for ( ni in seq( along = names(arzidxs) ) ) {
  i <- order(toupper(names(sarrnovdf)[seq(length(names(arzidxs)))]))[ni]
  cat("\n\n### --- ", names(sarrnovdf)[i])
  lmifitageLE60 <- lm(sarrnovdf[ageLE60idxp, i] ~ sarrnovdf[ageLE60idxp, "age_at_diagnosis"])
  smrylmifitageLE60 <- summary(lmifitageLE60)
  wclfitageLE60[[i]] <- lmifitageLE60
  print(smrylmifitageLE60)
  lmifit <- lm(sarrnovdf[, i] ~ sarrnovdf[, "age_at_diagnosis"])
  smrylmifit <- summary(lmifit)
  wclfit[[i]] <- lmifit
  print(smrylmifit)
  plot(sarrnovdf[, "age_at_diagnosis"], sarrnovdf[, i],
       ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
       main = paste(stj, names(sarrnovdf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
  lines(supsmu(sarrnovdf[, "age_at_diagnosis"], sarrnovdf[, i], span = 0.4, bass = 10), lwd = 3)
  abline(h = 0, v = 50, lty = 2)
  text(x = 20, y = 16.3, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  text(x = 22, y = 15.5, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
       adj = 0, cex = 0.85)
  if ( strsplit(names(sarrnovdf)[i], split = "\\|")[[1]][2] %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
    text(x = 98, y = 16.3, labels="Age-dependent trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 15.5, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
  } else {
    text(x = 98, y = 16.3, labels="No detectable trend:", adj = 1, cex = 0.85)
    text(x = 98, y = 15.5, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
  }
}
dev.off()

write.csv(arnovoutdf[which(arnovoutdf$arGeneSetidxp), -c(105:107)],
          file = "AgeRelated_arGeneSet_NoOverlapMB09BigSeries_lm_v11.csv")



### Random set of 4000 genes
length(unique(aroutdf$ILMN_Gene_0))  ## [1] 37848
set.seed(4133); write.table(unique(aroutdf$ILMN_Gene_0)[sample(length(unique(aroutdf$ILMN_Gene_0)))[1:4000]],
                          file = "Random4000_GeneNames.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### Pick 250 cases at random, look for age trends.


### Random250


### Biomarker shows age-dependent trend if
### ( ( abs(Slope for age < 60) > log2(1.5)/35 or abs(Slope for age < 60) > log2(1.5)/35 or abs(Slope for age) > log2(1.5)/70 ) AND
###   (p-val < 0.05/(2 * num biomarkers) ) )
### Save
###  3 slopes:,
###  3 pvalues:,
###  ILMN probeset ID: Probe_id, Probe_sequence,
###  gene name: Gene_symbol, Gene_synonyms, Synonyms_0, ILMN_Gene_0, Chromosome_0, Cytoband_0, Original_genomic_annotation
###  refseq nm_ code:  RefSeq_transcripts, RefSeq_ID_0
### LumAnov = Triple negative phenotype, No Overlap
set.seed(7717); sehRandom250df <- sehdf[sample(seq(nrow(sehdf)))[1:250], ]
arRandom250outdf <- annodf
arRandom250outmatcolnames <- c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
arRandom250outmat <- matrix(NA_real_, nrow = nrow(arRandom250outdf), ncol = length(arRandom250outmatcolnames))
LE60_slope_col <- match("LE60_slope", arRandom250outmatcolnames)
GT60_slope_col <- match("GT60_slope", arRandom250outmatcolnames)
AllAges_slope_col <- match("AllAges_slope", arRandom250outmatcolnames)
LE60_pval_col <- match("LE60_pval", arRandom250outmatcolnames)
GT60_pval_col <- match("GT60_pval", arRandom250outmatcolnames)
AllAges_pval_col <- match("AllAges_pval", arRandom250outmatcolnames)
AgeDependentp_col <- match("AgeDependentp", arRandom250outmatcolnames)
arRandom250outdf$LE60_slope <- rep(NA_real_ , nrow(arRandom250outdf))
arRandom250outdf$GT60_slope <- rep(NA_real_ , nrow(arRandom250outdf))
arRandom250outdf$AllAges_slope <- rep(NA_real_ , nrow(arRandom250outdf))
arRandom250outdf$LE60_pval <- rep(NA_real_ , nrow(arRandom250outdf))
arRandom250outdf$GT60_pval <- rep(NA_real_ , nrow(arRandom250outdf))
arRandom250outdf$AllAges_pval <- rep(NA_real_ , nrow(arRandom250outdf))
arRandom250outdf$AgeDependentp <- rep(NA , nrow(arRandom250outdf))
biosigslopeLE60 <-  log2(1.25)/35
biosigslopeGT60 <-  log2(1.25)/35
biosigslopeAllAges <-  log2(1.25)/70
multcompPval <- 0.05/(2 * nrow(arRandom250outdf))
Random250idxs <- match(sehRandom250df$MBid, dimnames(Dataset.r)[[1]])
arRandom250lmdf <- data.frame(age_at_diagnosis = sehRandom250df$age_at_diagnosis,
                       probesetni = Dataset.r[Random250idxs, 1])


cat("\n\n\n")
for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
  psi <- arRandom250outdf$ProbeId[ni]
  drci <- match(psi, dimnames(Dataset.r)[[2]])
  arRandom250lmdf$probesetni <- Dataset.r[Random250idxs, drci]
  if (ni %% 100 == 0 ) { cat(ni, ", ") }
  lmifitageLE60 <- lm(probesetni ~ age_at_diagnosis, data = arRandom250lmdf, na.action = na.exclude, subset = age_at_diagnosis <= 60)
  lmifitageGT60 <- lm(probesetni ~ age_at_diagnosis, data = arRandom250lmdf, na.action = na.exclude, subset = age_at_diagnosis > 60)
  lmifitallages <- lm(probesetni ~ age_at_diagnosis, data = arRandom250lmdf, na.action = na.exclude)
  smrylmifitageLE60 <- summary(lmifitageLE60)
  smrylmifitageGT60 <- summary(lmifitageGT60)
  smrylmifitallages <- summary(lmifitallages)

  arRandom250outmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
  arRandom250outmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
  arRandom250outmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
  arRandom250outmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
  arRandom250outmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
  arRandom250outmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
}

dimnames(arRandom250outmat)[[2]] <- arRandom250outmatcolnames
dimnames(arRandom250outmat)[[1]] <- arRandom250outdf$ProbeId
arRandom250outdf[, arRandom250outmatcolnames] <- arRandom250outmat
cat("\n\n\n")
### Adjust all p-values and redo age-dependent selection on adjusted p-values < 0.05
arRandom250outdf$LE60_BHadj_pval <- p.adjust(arRandom250outdf$LE60_pval, method="BH")
arRandom250outdf$GT60_BHadj_pval <- p.adjust(arRandom250outdf$GT60_pval, method="BH")
arRandom250outdf$AllAges_BHadj_pval <- p.adjust(arRandom250outdf$AllAges_pval, method="BH")
arRandom250outdf$BHadj_and_AgeDependentp <-
  ( ( ( abs(arRandom250outdf$LE60_slope) > biosigslopeLE60 )       & ( arRandom250outdf$LE60_BHadj_pval < 0.05 ) ) |
    ( ( abs(arRandom250outdf$GT60_slope) > biosigslopeGT60 )       & ( arRandom250outdf$GT60_BHadj_pval < 0.05 ) ) |
    ( ( abs(arRandom250outdf$AllAges_slope) > biosigslopeAllAges ) & ( arRandom250outdf$AllAges_BHadj_pval < 0.05 ) ) )
   
arRandom250outdf$BHadj_signifp <- ( ( ( arRandom250outdf$LE60_BHadj_pval < 0.05 ) |
                                 ( arRandom250outdf$GT60_BHadj_pval < 0.05 ) |
                                 ( arRandom250outdf$AllAges_BHadj_pval < 0.05 ) ) )

### Plot age-associated p-values by chromosome position - Manhattan plot

### qqman package has manhattan plot

require("qqman")

### Data frame with CHR, BP, P (and SNP to avoid warnings)
### Genomic_location
### chr2:206352192:206352241:+
arRandom250outdf$Chrcr <- sapply(strsplit(arRandom250outdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
arRandom250outdf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "", gsub("_h2_hap1", "", gsub("_random", "", arRandom250outdf$Chrcr))))
arRandom250outdf$Chrn <- as.numeric(gsub("chr", "", arRandom250outdf$Chrc))
arRandom250outdf$Chrn[arRandom250outdf$Chrc == "chrX"] <- 23
arRandom250outdf$Chrn[arRandom250outdf$Chrc == "chrY"] <- 24
with(arRandom250outdf, table(Chrcr, Chrn, useNA = "always"))
with(arRandom250outdf, table(Chrc, Chrn, useNA = "always"))
arRandom250outdf$CHR <- arRandom250outdf$Chrn
arRandom250outdf$Startcr <- sapply(strsplit(arRandom250outdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
arRandom250outdf$Stopcr <- sapply(strsplit(arRandom250outdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
arRandom250outdf$Startn <- as.numeric(arRandom250outdf$Startcr)
arRandom250outdf$Stopn <- as.numeric(arRandom250outdf$Stopcr)
arRandom250outdf$BP <- trunc((arRandom250outdf$Startn + arRandom250outdf$Stopn)/2)
arRandom250outdf$SNP <- arRandom250outdf$Probe_id
arRandom250outdf$P <- apply(arRandom250outdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
arRandom250outdf$PBHadj <- apply(arRandom250outdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
arRandom250outdfManhattanidxp <- arRandom250outdf$BHadj_and_AgeDependentp & arRandom250outdf$Chrn <= 23
arRandom250outdfManhattanidxp[is.na(arRandom250outdfManhattanidxp)] <- FALSE
### Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
arRandom250outdfManhattanVarsidxp <-
  apply(arRandom250outdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1, function(x) !any(is.na(x)))  & arRandom250outdf$Chrn <= 23

quartz()
manhattan(arRandom250outdf[arRandom250outdfManhattanidxp, c( "SNP", "CHR", "BP", "P")])
quartz()
### Use AgeRelated_Manhattan_AgeDependent_and_BHadjSignificantv01.pdf
### Discard AgeRelated_Manhattan_BHadjSignificantv01.pdf
pdf(file = "./Plots/AgeRelated_Manhattan_Random250_AgeDependent_and_BHadjSignificantv01.pdf", width = 12, height = 6)
sm_manhattan(arRandom250outdf[arRandom250outdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj")], p="PBHadj", chrlabs = c(1:22, "X"),
             suggestiveline = -log10(0.05), genomewideline = FALSE,
             plotpointidxp = arRandom250outdf[arRandom250outdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],
             ylab = "-log10(Benjamini-Hochberg adjusted P-values)",
             main = "Random250: Gene expression showing significant age-related trend: FC > 1.25")
dev.off()
quartz()
nomissidxp <- !( is.na(arRandom250outdf$CHR) | is.na(arRandom250outdf$BP) | is.na(arRandom250outdf$P) )
nomissidxp[is.na(nomissidxp)] <- FALSE
manhattan(arRandom250outdf[nomissidxp,  c( "SNP", "CHR", "BP", "P")])

write.csv(arRandom250outdf[arRandom250outdfManhattanidxp, -c(105:107)][
            -log10(arRandom250outdf[arRandom250outdfManhattanidxp, c( "PBHadj")]) > 1, ]$Gene_symbol,
          file = "arRandom250outdf_GeneSymbol_Top100AgeRelated.csv")
write.csv(arRandom250outdf[arRandom250outdfManhattanidxp, -c(105:107)][
            -log10(arRandom250outdf[arRandom250outdfManhattanidxp, c( "PBHadj")]) > 1, ],
          file = "arRandom250outdf_Top100AgeRelated.csv")
write.csv(arRandom250outdf[arRandom250outdfManhattanidxp, -c(105:107)],
          file = "AgeRelated_Random250_AgeDependent_and_BHadjSignificantv01.csv")

### > arLumAoutdf[arLumAoutdfManhattanidxp, -c(105:107)][
###          -log10(arLumAoutdf[arLumAoutdfManhattanidxp, c( "PBHadj")]) > 20, ]$Gene_symbol
###  [1] "KRT14"   "FOXD2"   "TP63"    "KRT17"   "DIO2"    "ESR1"    "PDGFRA" 
###  [8] "PDGFRA"  "FMO1"    "COL17A1" "COL17A1"

### Examine distribution of age-related.
table(arRandom250outdf$CHR) ### Expected
Ei <- table(arRandom250outdf[arRandom250outdf$Chrn <= 23, ]$CHR)/sum(table(arRandom250outdf[arRandom250outdf$Chrn <= 23, ]$CHR))*100
Oi <- table(arRandom250outdf[arRandom250outdfManhattanidxp, ]$CHR)/sum(table(arRandom250outdf[arRandom250outdfManhattanidxp, ]$CHR))*100
Oi - Ei
Oi/Ei
((Oi - Ei)^2)/Ei
sum(((Oi - Ei)^2)/Ei)
Eni <- table(arRandom250outdf[arRandom250outdf$Chrn <= 23, ]$CHR)
Oni <- table(arRandom250outdf[arRandom250outdfManhattanidxp, ]$CHR)
CrossTable(arRandom250outdf[arRandom250outdf$Chrn <= 23, ]$CHR, arRandom250outdf[arRandom250outdf$Chrn <= 23, ]$BHadj_AgeDependentp,
           digits = 1, expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=FALSE, prop.chisq=FALSE, chisq = TRUE, fisher=FALSE, mcnemar=FALSE,
           resid=FALSE, sresid=TRUE, asresid=FALSE,
           missing.include=FALSE,
           format="SPSS") ## By probeset

## Use tapply to reduce to genes
## Problem - gene names not consistent.  Unclear how to gather probesets into genes.
##   Table of gene names yields over 30,000 symbols for 21,000 genes if indeed 21,000 is even remotely close
##   to the "right" number of genes in the human genome.
### tapply( aroutdf$Gene_symbol, list(aroutdf$Chrn, aroutdf$BHadj_AgeDependentp),





       
### IRanges method:
require("MOutils")
require("IRanges")

sir <- IRanges(start = scdf$start, end = scdf$stop, names = leftPad0(as.character(scdf$chromosome)))

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col. = "red", sep = 0.5, maxf = NULL, yaxlabp = TRUE, ...)
{

### Adapted from plotRanges function in IRangesOverview.pdf (page 11) from Bioconductor, package IRanges
###   An Introduction to IRanges
### Patrick Aboyoun, Michael Lawrence, Herve Pages
### June 30, 2012
  xall <- x
  xlimall <- xlim
  xallnm <- unique(names(xall))
  xallnum <- as.numeric(xallnm)
  xlimallnm <- gsub("chr", "", names(xlimall))
  xlimallnum <- as.numeric(gsub("Y", "24", gsub("X", "23", xlimallnm)))
  nnm <- length(xlimallnm)
  par(mfcol = c(1, nnm), mar = c(0, 0, 0, 0), oma = c(10, 3, 10, 3))
  plot.new()

  for (nmi in seq(along = xlimallnm)) {
    if (nmi > 1) frame()
    height <- 1
    xtlab <- xlimallnm[nmi]
    xallidx <- which(xallnum %in% nmi)
    if ( length(xallidx) ) {
      namei <- xallnm[xallidx]
      x <- xall[names(xall) == namei]
      bins <- disjointBins(IRanges(start(x), end(x) + 1))
    } else {
      x <- NULL
      bins <- NULL
    }

    if (is(xlim, "Ranges")) {
      xlim <- c(min(start(xlim)), max(end(xlim)))
    } else {
      xlim <- c(1, xlimall[nmi])
    }
    if ( is.null(maxf) ) {
      maxbins <- max(bins)
    } else {
      maxbins <- maxf
    }
    plot.window(xlim, c(0, maxbins * (height + sep)))
    if ( !is.null(x) ) {
      ybottom <- bins * (sep + height) - height
      rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col., border = col., ...)
    }
    axis(1, at = mean(xlim), labels = xtlab, tick = TRUE, col = "black")
    abline(h = 0)
    abline(v = max(xlim), col = "gray80", lty = 2)
    if (nmi == 1) axis(2, labels = yaxlabp)
  }
  #title(main)
  #axis(1)
}
### Need a Ranges object showing chromosome min and max vals
library("BSgenome.Hsapiens.UCSC.hg19")
cdat <- seqlengths(Hsapiens)[1:24]/1000000
cdatprop <- cdat/sum(cdat)

pdf(file = "CNVconcordantv01.pdf", width = 8, height = 3)
plotRanges(sir, xlim = cdat[1:23], col. = "red", sep = 0, maxf = 20)
dev.off()

bir <- IRanges(start = bcdf$start, end = bcdf$stop, names = leftPad0(as.character(bcdf$chromosome)))

pdf(file = "CNVdiscordantv01.pdf", width = 8, height = 3)
plotRanges(bir, xlim = cdat[1:23], col. = "blue", sep = 0, maxf = 5)
dev.off()

ugnnvir <- IRanges(start = ugnnvdf$start_position, end = ugnnvdf$end_position, names = leftPad0(as.character(ugnnvdf$chromosome_name)))
pdf(file = "RelAmpOnlyDiffGeneNamesInCNVregionsCNVdiscordantv03.pdf", width = 8, height = 3)
plotRanges(ugnnvir, xlim = cdat[1:23], col. = "blue", sep = 0, maxf = 1, yaxlabp = FALSE)
dev.off()


### EZH2 CNAs by age
 ## CNAs copy number aberrations - somatic copy number events.
 ## Identify CNAs for EZH2.
 ## Cis-EZH2 CNAs
 ## EZH2:  Chr 7 Start: 148,504,464 bp from pter End: 148,581,441 bp from pter
 (148504464 + 148581441)/2 + c(-1, 1)*1500000  ## [1] 147042952 150042952

qchrnum <- 7
qchrstart <- 147042952
qchrend <- 150042952


ezh2.wgd.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ((mbcnadf$call2 == "HETD") | (mbcnadf$call2 == "HOMD")))
ezh2.pgdl.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ((mbcnadf$call2 == "HETD") | (mbcnadf$call2 == "HOMD")))
ezh2.pgdm.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ((mbcnadf$call2 == "HETD") | (mbcnadf$call2 == "HOMD")))
ezh2.pgdr.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ((mbcnadf$call2 == "HETD") | (mbcnadf$call2 == "HOMD")))
ezh2.cnad.posns <- ezh2.wgd.cna.posns | ezh2.pgdl.cna.posns | ezh2.pgdm.cna.posns | ezh2.pgdr.cna.posns

ezh2.cisd.MBids <- sort( unique( mbcnadf[ezh2.cnad.posns, "MBid"] ) )
length(ezh2.cisd.MBids)  ## 188 cases with a HETD or HOMD aberration within a 3Mb window of EZH2


ezh2.wga.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ((mbcnadf$call2 == "GAIN") | (mbcnadf$call2 == "AMP")))
ezh2.pgal.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ((mbcnadf$call2 == "GAIN") | (mbcnadf$call2 == "AMP")))
ezh2.pgam.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ((mbcnadf$call2 == "GAIN") | (mbcnadf$call2 == "AMP")))
ezh2.pgar.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ((mbcnadf$call2 == "GAIN") | (mbcnadf$call2 == "AMP")))
ezh2.cnaa.posns <- ezh2.wga.cna.posns | ezh2.pgal.cna.posns | ezh2.pgam.cna.posns | ezh2.pgar.cna.posns

ezh2.cisa.MBids <- sort( unique( mbcnadf[ezh2.cnaa.posns, "MBid"] ) )
### > length(ezh2.cisa.MBids)  ## [1] 217  cases with EZH2 somatic gain


sehdf$ezh2_cisCNAdel_v1n <- rep(0, nrow(sehdf))
sehdf[sehdf$MBid %in% ezh2.cisd.MBids, ]$ezh2_cisCNAdel_v1n <- 1

sehdf$ezh2_cisCNAgain_v1n <- rep(0, nrow(sehdf))
sehdf[sehdf$MBid %in% ezh2.cisa.MBids, ]$ezh2_cisCNAgain_v1n <- 1

sehrdf$ezh2_cisCNAdel_v1n <- rep(0, nrow(sehrdf))
sehrdf[sehrdf$MBid %in% ezh2.cisd.MBids, ]$ezh2_cisCNAdel_v1n <- 1

sehrdf$ezh2_cisCNAgain_v1n <- rep(0, nrow(sehrdf))
sehrdf[sehrdf$MBid %in% ezh2.cisa.MBids, ]$ezh2_cisCNAgain_v1n <- 1

all.equal(sort(ezh2.cisd.MBids), sort(sehdf[sehdf$ezh2_cisCNAdel_v1n == 1, ]$MBid))
all.equal(sort(ezh2.cisa.MBids), sort(sehdf[sehdf$ezh2_cisCNAgain_v1n == 1, ]$MBid))

simsmufit <- function(x, xdf = sehdf, ydf = NULL, xvar = "age_at_diagnosis", xhlbl = NULL, xlbl = NULL, ylbl = NULL,
                      tlbl = "", nsim = 2000, alpha = 0.05) {
  rn <- nsim
  alphalev <- alpha/2.0
  UCLn <- trunc(rn * (1 - alphalev))
  LCLn <- trunc(rn * alphalev)
  agesv <- sort(unique(xdf[, xvar]), na.last = NA)
  simmat <- matrix(NA_real_, nrow = rn, ncol = length(agesv))
  ylims <- range(xdf[, x], na.rm = TRUE)
  ylimd <- abs(diff(ylims))
###   ylimlo <- ylims[1] - (0.5 * ylimd)
  ylimlo <- ylims[1] - (0.2 * ylimd)
  ylimhi <- ylims[2] + (0.2 * ylimd)
  valdx <- !is.na(xdf[, x])
###   xlb <- strsplit(xvar, split = "_")[[1]][1]
###   xlbc <- paste(toupper(substring(xlb, 1, 1)), tolower(substring(xlb, 2)), sep = "")
  if(missing(xlbl)) {
    xlbc <- "Dx age (years)"
    xhlab <- "Age-dependent"
    xhlclab <- tolower(xhlab)
    xlims <- c(15, 95)
  } else {
   xlbc <- xlbl
   xhlab <- xhlbl
   xhlclab <- xhlab
   xlims <- range(xdf[valdx, xvar])
  }
  if(missing(ylbl)) {
    ylbc <- toupper(strsplit(x, split = "_")[[1]][1])
  } else {
    ylbc <- ylbl
  }
  xlabposn <- xlims[1] + 0.05*diff(xlims)

  plot( xdf[valdx, xvar], xdf[valdx, x], xlab = xlbc, ylab = ylbc,
       xlim = xlims, ylim = c(ylimlo, ylimhi), pch = ".", col = "grey60" )
###  axis(side = 2, at = c(0, 20, 40, 60, 80, 100))
  for (i in seq(rn) ) {
    ssout <- supsmu( sample(xdf[valdx, xvar]), xdf[valdx, x], bass = 10)
    spout <- predict(smooth.spline(ssout$x, ssout$y), x = agesv)
    simmat[i, ] <- spout$y
###     lines(ssout, lty = 3, lwd = 1, col = "lightgrey")
  }
  ##points( xdf[, xvar], xdf[, x])
  srtmat <- apply(simmat, 2, sort)
###   lines(supsmu(agesv, srtmat[LCLn, ], bass = 10), lty = 1, lwd = 3, type = "l", col = "grey80")
###   lines(supsmu(agesv, srtmat[UCLn, ], bass = 10), lty = 1, lwd = 3, type = "l", col = "grey80")
  polygon(c(agesv, rev(agesv)), c(srtmat[LCLn, ], rev(srtmat[UCLn, ])), col = "wheat3", border = "wheat3")
  ssout <- supsmu( xdf[valdx, xvar], xdf[valdx, x], bass = 10)
  lines(ssout, lwd = 2, type = "l")
  srtmatcols <- match(ssout$x, agesv)
  if ( length(ssout$x) <= length(agesv) ) {
    if ( any( ssout$y < srtmat[LCLn, srtmatcols] ) |
         any( ssout$y > srtmat[UCLn, srtmatcols] ) ) {
      idxs <- c( which( ssout$y < srtmat[LCLn, srtmatcols] ),
                 which( ssout$y > srtmat[UCLn, srtmatcols] ) )
      points(ssout$x[idxs], ssout$y[idxs], pch = 20)
      text(x = xlabposn, y = ylimhi, labels = paste(xhlab, " p<", alpha, sep = ""), adj = 0)
    } else {
      text(x = xlabposn, y = ylimhi, labels = paste("Not ", xhlclab, " p>", alpha, sep = ""), adj = 0)
    }
  } else {
    text(x = xlabposn, y = ylimhi, labels = "FIX ME missing ages", adj = 0)
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
         xlim = xlims, ylim = c(ylimlo, ylimhi), pch = ".", col = "grey60" )
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
    lines(supsmu( ydf[valdy, xvar], ydf[valdy, x], bass = 10), lwd = 3, type = "l")
###     legend("bottomright", legend = paste(trunc(100*(1-alpha)), "% Confidence Bounds", sep = ""),
###            lty = 1, lwd = 4, seg.len = 3, text.width = 60, col = "wheat3")
  }
}

cscvars <- c("ezh2_cisCNAdel_v1n")
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2del_Trends_95pct_2kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehdf, ylbl = "cis-EZH2 CNA: Hom or Het Del", nsim = 2000, alpha = 0.05)
dev.off()
date()

cscvars <- c("ezh2_cisCNAgain_v1n")
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Trends_95pct_2kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehdf, ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 2000, alpha = 0.05)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Trends_99pct_10kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehdf, ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 10000, alpha = 0.01)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Trends_99p9pct_40kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehdf, ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 40000, alpha = 0.001)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Trends_99p99pct_400kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehdf, ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 400000, alpha = 0.0001)
dev.off()
date()

### Check expression by CNA.
cscvars <- c("ezh2_cisCNAdel_v1n")
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2del_Expression_Trends_95pct_2kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf,
       xvar = "EZH2|r|ILMN_2364529", xlbl = "EZH2 expression: EZH2|r|ILMN_2364529", xhlbl = "Associated",
       ylbl = "cis-EZH2 CNA: Hom or Het Del", nsim = 2000, alpha = 0.05)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2del_Expression_Trends_99pct_10kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf,
       xvar = "EZH2|r|ILMN_2364529", xlbl = "EZH2 expression: EZH2|r|ILMN_2364529", xhlbl = "Associated",
       ylbl = "cis-EZH2 CNA: Hom or Het Del", nsim = 10000, alpha = 0.01)
dev.off()
date()
cscvars <- c("ezh2_cisCNAdel_v1n")
pdf(file = "./Plots/MBEX_CNA_cisEZH2del_Expression_Trends_99p9pct_40kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf,
       xvar = "EZH2|r|ILMN_2364529", xlbl = "EZH2 expression: EZH2|r|ILMN_2364529", xhlbl = "Associated",
       ylbl = "cis-EZH2 CNA: Hom or Het Del", nsim = 40000, alpha = 0.001)
dev.off()
date()


cscvars <- c("ezh2_cisCNAgain_v1n")
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Expression_Trends_95pct_2kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf,
       xvar = "EZH2|r|ILMN_2364529", xlbl = "EZH2 expression: EZH2|r|ILMN_2364529", xhlbl = "Associated",
       ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 2000, alpha = 0.05)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Expression_Trends_99pct_10kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf,
       xvar = "EZH2|r|ILMN_2364529", xlbl = "EZH2 expression: EZH2|r|ILMN_2364529", xhlbl = "Associated",
       ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 10000, alpha = 0.01)
dev.off()
date()
cscvars <- c("ezh2_cisCNAgain_v1n")
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Expression_Trends_99p9pct_40kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf,
       xvar = "EZH2|r|ILMN_2364529", xlbl = "EZH2 expression: EZH2|r|ILMN_2364529", xhlbl = "Associated",
       ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 40000, alpha = 0.001)
dev.off()  ## 99p9 sims about 10 mins
date()
cscvars <- c("ezh2_cisCNAgain_v1n")
pdf(file = "./Plots/MBEX_CNA_cisEZH2gain_Expression_Trends_99p99pct_400kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf,
       xvar = "EZH2|r|ILMN_2364529", xlbl = "EZH2 expression: EZH2|r|ILMN_2364529", xhlbl = "Associated",
       ylbl = "cis-EZH2 CNA: Gain or High level Amp", nsim = 400000, alpha = 0.0001)
dev.off()  ## 99p99 sims about 130 mins
date()

### EZH2 gene only.
 ## CNAs copy number aberrations - somatic copy number events.
 ## Identify CNAs for EZH2.
 ## Cis-EZH2 CNAs
 ## EZH2:  Chr 7 Start: 148,504,464 bp from pter End: 148,581,441 bp from pter
## (148504464 + 148581441)/2 + c(-1, 1)*1500000  ## [1] 147042952 150042952

qchrnum <- 7
qchrstart <- 148504464
qchrend <- 148581441


ezh2_only.wgv.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ((mbcnadf$call2 == "CNV") ))
ezh2_only.pgvl.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ((mbcnadf$call2 == "CNV") ))
ezh2_only.pgvm.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ((mbcnadf$call2 == "CNV") ))
ezh2_only.pgvr.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ((mbcnadf$call2 == "CNV") ))
ezh2_only.cnav.posns <- ezh2_only.wgv.cna.posns | ezh2_only.pgvl.cna.posns | ezh2_only.pgvm.cna.posns | ezh2_only.pgvr.cna.posns

ezh2_only.cisv.MBids <- sort( unique( mbcnadf[ezh2_only.cnav.posns, "MBid"] ) )
length(ezh2_only.cisv.MBids)  ## 0 cases with a CNV in EZH2_ONLY


ezh2_only.wgn.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ((mbcnadf$call2 == "NEUT") ))
ezh2_only.pgnl.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ((mbcnadf$call2 == "NEUT") ))
ezh2_only.pgnm.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ((mbcnadf$call2 == "NEUT") ))
ezh2_only.pgnr.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ((mbcnadf$call2 == "NEUT") ))
ezh2_only.cnan.posns <- ezh2_only.wgn.cna.posns | ezh2_only.pgnl.cna.posns | ezh2_only.pgnm.cna.posns | ezh2_only.pgnr.cna.posns

ezh2_only.cisn.MBids <- sort( unique( mbcnadf[ezh2_only.cnan.posns, "MBid"] ) )
length(ezh2_only.cisn.MBids)  ## 1810 cases with a NEUT in EZH2_ONLY


ezh2_only.wgdd.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ((mbcnadf$call2 == "HOMD") ))
ezh2_only.pgddl.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ((mbcnadf$call2 == "HOMD") ))
ezh2_only.pgddm.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ((mbcnadf$call2 == "HOMD") ))
ezh2_only.pgddr.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ((mbcnadf$call2 == "HOMD") ))
ezh2_only.cnadd.posns <- ezh2_only.wgdd.cna.posns | ezh2_only.pgddl.cna.posns | ezh2_only.pgddm.cna.posns | ezh2_only.pgddr.cna.posns

ezh2_only.cisdd.MBids <- sort( unique( mbcnadf[ezh2_only.cnadd.posns, "MBid"] ) )
length(ezh2_only.cisdd.MBids)  ## 3 cases with a HOMD in EZH2_ONLY


ezh2_only.wgd.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ((mbcnadf$call2 == "HETD") ))
ezh2_only.pgdl.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ((mbcnadf$call2 == "HETD") ))
ezh2_only.pgdm.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ((mbcnadf$call2 == "HETD") ))
ezh2_only.pgdr.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ((mbcnadf$call2 == "HETD") ))
ezh2_only.cnad.posns <- ezh2_only.wgd.cna.posns | ezh2_only.pgdl.cna.posns | ezh2_only.pgdm.cna.posns | ezh2_only.pgdr.cna.posns

ezh2_only.cisd.MBids <- sort( unique( mbcnadf[ezh2_only.cnad.posns, "MBid"] ) )
length(ezh2_only.cisd.MBids)  ## 59 cases with a HETD in EZH2_ONLY


ezh2_only.wga.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ((mbcnadf$call2 == "GAIN") ))
ezh2_only.pgal.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ((mbcnadf$call2 == "GAIN") ))
ezh2_only.pgam.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ((mbcnadf$call2 == "GAIN") ))
ezh2_only.pgar.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ((mbcnadf$call2 == "GAIN") ))
ezh2_only.cnaa.posns <- ezh2_only.wga.cna.posns | ezh2_only.pgal.cna.posns | ezh2_only.pgam.cna.posns | ezh2_only.pgar.cna.posns

ezh2_only.cisa.MBids <- sort( unique( mbcnadf[ezh2_only.cnaa.posns, "MBid"] ) )
### > length(ezh2_only.cisa.MBids)  ## [1] 140  cases with EZH2_ONLY somatic gain

ezh2_only.wgaa.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrend) &
                          ( (mbcnadf$call2 == "AMP")))
ezh2_only.pgaal.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrstart) & (mbcnadf$loc.end >= qchrstart) &
                           ( (mbcnadf$call2 == "AMP")))
ezh2_only.pgaam.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start >= qchrstart) & (mbcnadf$loc.end <= qchrend) &
                           ( (mbcnadf$call2 == "AMP")))
ezh2_only.pgaar.cna.posns <- ((mbcnadf$chrom == qchrnum) & (mbcnadf$loc.start <= qchrend) & (mbcnadf$loc.end >= qchrend) &
                           ( (mbcnadf$call2 == "AMP")))
ezh2_only.cnaaa.posns <- ezh2_only.wgaa.cna.posns | ezh2_only.pgaal.cna.posns | ezh2_only.pgaam.cna.posns | ezh2_only.pgaar.cna.posns

ezh2_only.cisaa.MBids <- sort( unique( mbcnadf[ezh2_only.cnaaa.posns, "MBid"] ) )
### > length(ezh2_only.cisaa.MBids)  ## [1] 6  cases with EZH2_ONLY somatic gain
intersect(ezh2_only.cisdd.MBids, ezh2_only.cisd.MBids)
intersect(c(ezh2_only.cisdd.MBids, ezh2_only.cisd.MBids), ezh2_only.cisa.MBids)
intersect(c(ezh2_only.cisdd.MBids, ezh2_only.cisd.MBids, ezh2_only.cisa.MBids), ezh2_only.cisaa.MBids)
intersect(ezh2_only.cisa.MBids, ezh2_only.cisaa.MBids)
## [1] "MB.6163"
## Build ezh2_only_cna state var and do against age

sehdf$ezh2_onlyCNA_v1n <- rep(NA_real_, nrow(sehdf))
sehdf[sehdf$MBid %in% ezh2_only.cisn.MBids, ]$ezh2_onlyCNA_v1n <- 2
sehdf[sehdf$MBid %in% ezh2_only.cisd.MBids, ]$ezh2_onlyCNA_v1n <- 1
sehdf[sehdf$MBid %in% ezh2_only.cisa.MBids, ]$ezh2_onlyCNA_v1n <- 4
sehdf[sehdf$MBid %in% ezh2_only.cisdd.MBids, ]$ezh2_onlyCNA_v1n <- 0
sehdf[sehdf$MBid %in% ezh2_only.cisaa.MBids, ]$ezh2_onlyCNA_v1n <- 8

sehrdf$ezh2_onlyCNA_v1n <- rep(NA_real_, nrow(sehrdf))
sehrdf[sehrdf$MBid %in% ezh2_only.cisn.MBids, ]$ezh2_onlyCNA_v1n <- 2
sehrdf[sehrdf$MBid %in% ezh2_only.cisd.MBids, ]$ezh2_onlyCNA_v1n <- 1
sehrdf[sehrdf$MBid %in% ezh2_only.cisa.MBids, ]$ezh2_onlyCNA_v1n <- 4
sehrdf[sehrdf$MBid %in% ezh2_only.cisdd.MBids, ]$ezh2_onlyCNA_v1n <- 0
sehrdf[sehrdf$MBid %in% ezh2_only.cisaa.MBids, ]$ezh2_onlyCNA_v1n <- 8

write.csv(sehrdf, file = "AgeDependent_sehrdf.csv")
table(sehrdf$MBid %in% asudf$MBid, useNA = "always")
sehrNoMB09df <- sehrdf[sehrdf$MBid %in% asudf$MBid, ]
all.equal(sort(sehrNoMB09df$MBid), sort(asudf$MBid)) ## TRUE


cscvars <- c("ezh2_onlyCNA_v1n")
date()
pdf(file = "./Plots/MBEX_CNA_EZH2_only_Age_Trends_sehrdf_95pct_2kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf, ylbl = "EZH2 only CNA: Hom, Het, Neut, Gain, Amp", nsim = 2000, alpha = 0.05)
dev.off()
date()

cscvars <- c("ezh2_onlyCNA_v1n")
date()
pdf(file = "./Plots/MBEX_CNA_EZH2_only_Age_Trends_sehrdf_99pct_10kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xdf = sehrdf, ylbl = "EZH2 only CNA: Hom, Het, Neut, Gain, Amp", nsim = 10000, alpha = 0.01)
dev.off()
date()

### Horizontal axis expression:
cscvars <- c("ezh2_onlyCNA_v1n")
date()
pdf(file = "./Plots/MBEX_CNA_EZH2_only_Expression_Trends_sehrdf_95pct_2kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xvar = "EZH2|r|ILMN_1708105", xlbl = "EZH2 expression",
       xhlbl = "Expression associated", 
       xdf = sehrdf, ylbl = "EZH2 only CNA: Hom, Het, Neut, Gain, Amp", nsim = 2000, alpha = 0.05)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_EZH2_only_Expression_Trends_sehrdf_99pct_10kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xvar = "EZH2|r|ILMN_1708105", xlbl = "EZH2 expression",
       xhlbl = "Expression associated", 
       xdf = sehrdf, ylbl = "EZH2 only CNA: Hom, Het, Neut, Gain, Amp", nsim = 10000, alpha = 0.01)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_EZH2_only_Expression_Trends_sehrdf_99p9pct_40kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xvar = "EZH2|r|ILMN_1708105", xlbl = "EZH2 expression",
       xhlbl = "Expression associated", 
       xdf = sehrdf, ylbl = "EZH2 only CNA: Hom, Het, Neut, Gain, Amp", nsim = 40000, alpha = 0.001)
dev.off()
date()
pdf(file = "./Plots/MBEX_CNA_EZH2_only_Expression_Trends_sehrdf_99p99pct_400kSimCIs_shaded_v01.pdf", width = 6, height = 6)
### par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, xvar = "EZH2|r|ILMN_1708105", xlbl = "EZH2 expression",
       xhlbl = "Expression associated", 
       xdf = sehrdf, ylbl = "EZH2 only CNA: Hom, Het, Neut, Gain, Amp", nsim = 400000, alpha = 0.0001)
dev.off()
date()


NULL

###  TODO: 20150501
###  EZH2_associated_DNAinfo_v03.txt has Gene, Chromosome, Start and Stop for Tomo's gene list.
###  Use this and code above to automate building CNA vars like EZH2 one above.

###
### Venn diagram 20150511
###
table(aroutdf$BHadj_and_AgeDependentp)
table(arLumAoutdf$BHadj_and_AgeDependentp)
table(aroutdf$BHadj_and_AgeDependentp & arLumAoutdf$BHadj_and_AgeDependentp)

arLumAoutdf$Probe_id[which(arLumAoutdf$BHadj_and_AgeDependentp)]
arLumAoutdf$Symbol_0[which(arLumAoutdf$BHadj_and_AgeDependentp)]


AgeRelated_AllMETABRIC_AgeDependent_and_BHadjSignificant_GeneNames_v05.txt
AgeRelated_AllMETABRIC_LumA_AgeDependent_and_BHadjSignificant_GeneSymbols_v01.txt
AgeRelated_AllMETABRIC_LumB_AgeDependent_and_BHadjSignificant_GeneNames_v01.txt

table(duplicated(aroutdf[aroutdfManhattanidxp, ][!duplicated(aroutdf[aroutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0))
table(duplicated(arLumAoutdf[arLumAoutdfManhattanidxp, ][!duplicated(arLumAoutdf[arLumAoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0))
table(duplicated(arLumBoutdf[arLumBoutdfManhattanidxp, ][!duplicated(arLumBoutdf[arLumBoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0))
table(duplicated(arHer2outdf[arHer2outdfManhattanidxp, ][!duplicated(arHer2outdf[arHer2outdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0))
table(duplicated(arBasaloutdf[arBasaloutdfManhattanidxp, ][!duplicated(arBasaloutdf[arBasaloutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0))
table(duplicated(arNormaloutdf[arNormaloutdfManhattanidxp, ][!duplicated(arNormaloutdf[arNormaloutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0))

length(aroutdf[which(!duplicated(aroutdf[aroutdfManhattanidxp, ][!duplicated(aroutdf[aroutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)), ]$ILMN_Gene_0)
WholeCohortGN <- unique(aroutdf[aroutdfManhattanidxp, ][
                   !duplicated(aroutdf[aroutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)
LumAGN <- unique(arLumAoutdf[arLumAoutdfManhattanidxp, ][
                   !duplicated(arLumAoutdf[arLumAoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)
LumBGN <- unique(arLumBoutdf[arLumBoutdfManhattanidxp, ][
                   !duplicated(arLumBoutdf[arLumBoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)
Her2GN <- unique(arHer2outdf[arHer2outdfManhattanidxp, ][
                   !duplicated(arHer2outdf[arHer2outdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)
BasalGN <- unique(arBasaloutdf[arBasaloutdfManhattanidxp, ][
                   !duplicated(arBasaloutdf[arBasaloutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)
NormalGN <- unique(arNormaloutdf[arNormaloutdfManhattanidxp, ][
                   !duplicated(arNormaloutdf[arNormaloutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)


length(intersect(WholeCohortGN, LumAGN))

x <- list(WholeCohortGN = WholeCohortGN,
          LumAGN = LumAGN,
          LumBGN = LumBGN,
          Her2GN = Her2GN,
          BasalGN = BasalGN,
          NormalGN = NormalGN)

WABHBx <- list(WholeCohortGN = WholeCohortGN,
          LumAGN = LumAGN,
          LumBGN = LumBGN,
          Her2GN = Her2GN,
          BasalGN = BasalGN)
require("VennDiagram")
WABHB <- venn.diagram(WABHBx, filename = "./Plots/VennDiagram_WABHB_v01.png", imagetype="png")

WABHNx <- list(WholeCohortGN = WholeCohortGN,
          LumAGN = LumAGN,
          LumBGN = LumBGN,
          Her2GN = Her2GN,
          NormalGN = NormalGN)
require("VennDiagram")
WABHN <- venn.diagram(WABHNx, filename = "./Plots/VennDiagram_WABHN_v01.png", imagetype="png")

WABBNx <- list(WholeCohortGN = WholeCohortGN,
          LumAGN = LumAGN,
          LumBGN = LumBGN,
          BasalGN = BasalGN,
          NormalGN = NormalGN)
require("VennDiagram")
WABBN <- venn.diagram(WABBNx, filename = "./Plots/VennDiagram_WABBN_v01.png", imagetype="png")

ABHBNx <- list(
          LumAGN = LumAGN,
          LumBGN = LumBGN,
          Her2GN = Her2GN,
          BasalGN = BasalGN,
          NormalGN = NormalGN)
require("VennDiagram")
ABHBN <- venn.diagram(ABHBNx, filename = "./Plots/VennDiagram_ABHBN_v01.png", imagetype="png")





NULL
###--------------------------------------------------------------------------------------------------
###
###  Other gene sets to average over
###
###--------------------------------------------------------------------------------------------------

###
EZH2associatedGenesetNamesv <-
  c(
###     "01_Yu2007_Cancer_Cell_Table_S1",
###     "02_Yu2007_Cancer_Cell_Table_S2",
    "03_Yu2007_Cancer_Res_Table_S2",
    "04_Squazzo2006_S6",
    "05_Tan2007",
    "06_Ku2008_Human",
    "07_Bracken2006_S2",
    "08_Xu2012",
    "09_Nuytten2008_CAUTION_EZH2_activated",
    "10_Nuytten2008_CAUTION_EZH2_repressed",
    "11_Lee2006_S9_Common",
    "12_Lee2011_EZH2_activated",
    "13_Lee2011_EZH2_repressed",
    "14_Cheng2011_Common",
    "15_Gupta2010_S7")

for ( oGLi in seq(along = EZH2associatedGenesetNamesv) ) {
  ogaGeneList <- EZH2associatedGenesetNamesv[oGLi] ## "01_Yu2007_Cancer_Cell_Table_S1"
  ogaGeneListFilename <- paste(ogaGeneList, ".csv", sep = "")
  ogadf <- read.table(file = paste("../Data/EZH2associated/", ogaGeneListFilename, sep = ""),
                      sep = ",", header = TRUE, stringsAsFactors = FALSE)

  date()
  ogaidxslist <- lapply(ogadf$GeneName, function(x) geneNameRowMatches_annodf(x) )
  date()
  names(ogaidxslist) <- ogadf$GeneName

  ogaidxs <- unlist(ogaidxslist)
  names(ogaidxs) <- rep(names(ogaidxslist), lapply(ogaidxslist, length))
  ogaidxs <- ogaidxs[!duplicated(ogaidxs)] ## remove duplicated probeset entries.
  ogaidxlabs <- cbind(names(ogaidxs), annodf$ProbeId[ogaidxs], annodf$Gene_symbol[ogaidxs])
  annodf$ProbeId[ogaidxs]
  annodf$Gene_symbol[ogaidxs]


### Find corresponding indexes in expression data matrix
  ogazidxs <- rep(NA_integer_, length(ogaidxs))
  for ( i in seq(along = ogazidxs)) {
    ogazidxs[i] <- match(annodf$ProbeId[ogaidxs[i]], dimnames(Dataset.rzs)[[2]])
  }
  names(ogazidxs) <- names(ogaidxs)


  oga.rzs <- Dataset.rzs[, ogazidxs]
  oga.r <- Dataset.r[, ogazidxs]
  dim(oga.rzs)
  dim(oga.r)
  oga.rzs[oga.rzs > 4] <- 4
  oga.rzs[oga.rzs < (-4)] <- (-4)

  ## Compute average and graph
  ## Data is in sogardf, from oga.r matrix.
  Geneset.r.mean <- rowMeans(oga.r, na.rm = TRUE)

  NULL
### Split cases at 60 years age at diagnosis.
### Trend plots first.
### Whole cohort
### Raw data  Whole cohort 1992 cases
### Have to calculate statistical significance and biological relevance
### across the current gene set (nnn probes).
### Using results from the genome-wide set here is inappropriate.
  aroutdf[, paste("EZH2associated", ogaGeneList, "LE60_BHadj_pval", sep = "_")] <- rep(NA_real_, nrow(aroutdf))
  aroutdf[, paste("EZH2associated", ogaGeneList, "GT60_BHadj_pval", sep = "_")] <- rep(NA_real_, nrow(aroutdf))
  aroutdf[, paste("EZH2associated", ogaGeneList, "AllAges_BHadj_pval", sep = "_")] <- rep(NA_real_, nrow(aroutdf))
  aroutdf[, paste("EZH2associated", ogaGeneList, "BHadj_and_AgeDependentp", sep = "_")] <- rep(NA, nrow(aroutdf))

  aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")] <- ( aroutdf$ProbeId %in% annodf$ProbeId[ogaidxs] )
  
### aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
###   p.adjust(aroutdf[aroutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
  aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("EZH2associated", ogaGeneList, "LE60_BHadj_pval", sep = "_")] <-
            p.adjust(aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$LE60_pval, method="BH")
  aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("EZH2associated", ogaGeneList, "GT60_BHadj_pval", sep = "_")] <-
            p.adjust(aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$GT60_pval, method="BH")
  aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("EZH2associated", ogaGeneList, "AllAges_BHadj_pval", sep = "_")] <-
            p.adjust(aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$AllAges_pval, method="BH")

### aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
###   ( ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
###      ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < 0.05 ) ) |
###     ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
###      ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < 0.05 ) ) |
###     ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
###      ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < 0.05 ) ) )
  aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("EZH2associated", ogaGeneList, "BHadj_and_AgeDependentp", sep = "_")] <-
            ( ( ( abs(aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$LE60_slope) > biosigslopeLE60 )       &
               ( aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")],
                         paste("EZH2associated", ogaGeneList, "LE60_BHadj_pval", sep = "_")] < 0.05 ) ) |
             ( ( abs(aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$GT60_slope) > biosigslopeGT60 )       &
              ( aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")],
                        paste("EZH2associated", ogaGeneList, "GT60_BHadj_pval", sep = "_")] < 0.05 ) ) |
             ( ( abs(aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$AllAges_slope) > biosigslopeAllAges ) &
              ( aroutdf[aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")],
                        paste("EZH2associated", ogaGeneList, "AllAges_BHadj_pval", sep = "_")] < 0.05 ) ) )

  sogadf <- cbind(oga.rzs, cldvdf[match(dimnames(oga.rzs)[[1]], cldvdf$MBid), ])
  sogardf <- cbind(Geneset.r.mean, oga.r, cldvdf[match(dimnames(oga.r)[[1]], cldvdf$MBid), ])
  names(sogadf)[seq(length(names(ogazidxs)))] <-
    paste(ogadf[match(names(ogaidxs), ogadf$GeneName), ], names(sogadf)[seq(length(names(ogazidxs)))], sep = "|")
  names(sogardf)[seq(length(names(ogazidxs))) + 1] <-
    paste(ogadf[match(names(ogaidxs), ogadf$GeneName), ], names(sogardf)[seq(length(names(ogazidxs))) + 1], sep = "|r|")


### Regression for age <= 60, age > 60   i <- 10; j <- 4
  wclfitageLE60 <- vector("list", length(names(ogazidxs)))
  wclfit <- vector("list", length(names(ogazidxs)))
  stj <- paste(dim(sogardf)[1], "cases")

  if ( !file.exists("./Plots/EZH2associated") ) dir.create("./Plots/EZH2associated", recursive = TRUE)

  ogaGeneListOutFilename <- paste("AgeRelated_raw_All1992Cases_", ogaGeneList, "_lm_v04.pdf", sep = "")
  pdf(file = paste("./Plots/EZH2associated/", ogaGeneListOutFilename, sep = ""), width = 8, height = 10)
### AgeDependent_and_BHadj_FC1.25_ProbeIds <-
###   aroutdf[which(aroutdf$arGeneSet_BHadj_and_AgeDependentp), "Probe_id"]
  AgeDependent_and_BHadj_FC1.25_ProbeIds <-
    aroutdf[which(aroutdf[, paste("EZH2associated", ogaGeneList, "BHadj_and_AgeDependentp", sep = "_")]), "Probe_id"]

  xlims <- c(20, 100)
  ylims <- c(4, 17)
  par(mfrow = c(2, 2))
  ageLE60idxp <- (sogardf$age_at_diagnosis <= 60)

  for ( ni in seq( along = c("Geneset.r.mean", names(ogazidxs) ) ) ) {
    i <- c(1, 1 + order( toupper(names(sogardf)[1 + seq(length(names(ogazidxs)))]) ) )[ni]
    cat("\n\n ### --- ", names(sogardf)[i])
    lmifitageLE60 <- lm(sogardf[ageLE60idxp, i] ~ sogardf[ageLE60idxp, "age_at_diagnosis"])
    lmifitageGT60 <- lm(sogardf[!ageLE60idxp, i] ~ sogardf[!ageLE60idxp, "age_at_diagnosis"])
    smrylmifitageLE60 <- summary(lmifitageLE60)
    smrylmifitageGT60 <- summary(lmifitageGT60)
    wclfitageLE60[[i]] <- lmifitageLE60
    print(smrylmifitageLE60)
    print(smrylmifitageGT60)
    lmifit <- lm(sogardf[, i] ~ sogardf[, "age_at_diagnosis"])
    smrylmifit <- summary(lmifit)
    wclfit[[i]] <- lmifit
    print(smrylmifit)
    plot(sogardf[, "age_at_diagnosis"], sogardf[, i],
         ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
         main = paste(stj, names(sogardf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
    lines(supsmu(sogardf[, "age_at_diagnosis"], sogardf[, i], span = 0.4, bass = 10), lwd = 3)
    abline(h = 0, v = 50, lty = 2)
    text(x = 20, y = 16.3, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
         adj = 0, cex = 0.85)
    text(x = 22, y = 15.5, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
         adj = 0, cex = 0.85)
    if ( i == 1 ) {
      ##     biosigslopeLE60    biosigslopeGT60     biosigslopeAllAges
      nOtherGenesetAssociations <- 25
      if ( ( ( abs(smrylmifitageLE60$coefficients[2, 1]) > biosigslopeLE60 ) &
            ( smrylmifitageLE60$coefficients[2, 4] < 0.05/nOtherGenesetAssociations ) ) |
          ( ( abs(smrylmifitageGT60$coefficients[2, 1]) > biosigslopeGT60 ) &
           ( smrylmifitageGT60$coefficients[2, 4] < 0.05/nOtherGenesetAssociations ) ) |
          ( ( abs(smrylmifit$coefficients[2, 1]) > biosigslopeAllAges ) &
           ( smrylmifit$coefficients[2, 4] < 0.05/nOtherGenesetAssociations ) ) ) {
        text(x = 98, y = 16.3, labels="Age-dependent trend:", adj = 1, cex = 0.85)
        text(x = 98, y = 15.5, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
      } else {
        text(x = 98, y = 16.3, labels="No detectable trend:", adj = 1, cex = 0.85)
        text(x = 98, y = 15.5, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
      }
    }
    if ( strsplit(names(sogardf)[i], split = "\\|")[[1]][3] %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
      text(x = 98, y = 16.3, labels="Age-dependent trend:", adj = 1, cex = 0.85)
      text(x = 98, y = 15.5, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
    } else {
      text(x = 98, y = 16.3, labels="No detectable trend:", adj = 1, cex = 0.85)
      text(x = 98, y = 15.5, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
    }
  }
  dev.off()
### [105] "Ontology_Component_0" [106] "Ontology_Process_0"  [107] "Ontology_Function_0"
  outcolidxs <- c(1:104, 108:145, which(grepl(ogaGeneList, names(aroutdf))) )
  write.csv(aroutdf[which(aroutdf[, paste("EZH2associated", ogaGeneList, "GeneSetidxp", sep = "_")]), outcolidxs],
            file = paste("./Plots/EZH2associated/AgeRelated_raw_All1992Cases_", ogaGeneList, "_lm_v04.csv", sep = ""))
}


NULL
###-------------------------------------------------------------------------------------------------------------------
### End EZH2associated gene sets to average over

###
ERalphaassociatedGenesetNamesv <-
  c(
    "01_Kwon_PNAS2007",
    "02_Romano_Mol_Cell_Endocrin2010",
    "03_Lin_Cancer_Res_2007",
    "04_Laganiere_PNAS_2005",
    "05_Cheng_Mol_Cell_2006",
    "06_Bourdeau_Mol_Endocrinol2004",
    "07_Lin_Genome_Biol_2004",
    "08_Lin_PLoS_Genet_2007",
    "09_Bourdeau_Nuc_Acid_Res_2008",
    "10_Yamaga_Horm_Cancer_2013",
    "11_Jin_Nucleic_Acids_Res_2004"
    )

for ( oGLi in seq(along = ERalphaassociatedGenesetNamesv) ) {
  ogaGeneList <- ERalphaassociatedGenesetNamesv[oGLi] ## "01_Yu2007_Cancer_Cell_Table_S1"
  ogaGeneListFilename <- paste(ogaGeneList, ".csv", sep = "")
  ogadf <- read.table(file = paste("../Data/ERalphaassociated/", ogaGeneListFilename, sep = ""),
                      sep = ",", header = TRUE, stringsAsFactors = FALSE)

  date()
  ogaidxslist <- lapply(ogadf$GeneName, function(x) geneNameRowMatches_annodf(x) )
  date()
  names(ogaidxslist) <- ogadf$GeneName

  ogaidxs <- unlist(ogaidxslist)
  names(ogaidxs) <- rep(names(ogaidxslist), lapply(ogaidxslist, length))
  ogaidxs <- ogaidxs[!duplicated(ogaidxs)] ## remove duplicated probeset entries.
  ogaidxlabs <- cbind(names(ogaidxs), annodf$ProbeId[ogaidxs], annodf$Gene_symbol[ogaidxs])
  annodf$ProbeId[ogaidxs]
  annodf$Gene_symbol[ogaidxs]


### Find corresponding indexes in expression data matrix
  ogazidxs <- rep(NA_integer_, length(ogaidxs))
  for ( i in seq(along = ogazidxs)) {
    ogazidxs[i] <- match(annodf$ProbeId[ogaidxs[i]], dimnames(Dataset.rzs)[[2]])
  }
  names(ogazidxs) <- names(ogaidxs)


  oga.rzs <- Dataset.rzs[, ogazidxs]
  oga.r <- Dataset.r[, ogazidxs]
  dim(oga.rzs)
  dim(oga.r)
  oga.rzs[oga.rzs > 4] <- 4
  oga.rzs[oga.rzs < (-4)] <- (-4)

  ## Compute average and graph
  ## Data is in sogardf, from oga.r matrix.
  Geneset.r.mean <- rowMeans(oga.r, na.rm = TRUE)

  NULL
### Split cases at 60 years age at diagnosis.
### Trend plots first.
### Whole cohort
### Raw data  Whole cohort 1992 cases
### Have to calculate statistical significance and biological relevance
### across the current gene set (nnn probes).
### Using results from the genome-wide set here is inappropriate.
  aroutdf[, paste("ERalphaassociated", ogaGeneList, "LE60_BHadj_pval", sep = "_")] <- rep(NA_real_, nrow(aroutdf))
  aroutdf[, paste("ERalphaassociated", ogaGeneList, "GT60_BHadj_pval", sep = "_")] <- rep(NA_real_, nrow(aroutdf))
  aroutdf[, paste("ERalphaassociated", ogaGeneList, "AllAges_BHadj_pval", sep = "_")] <- rep(NA_real_, nrow(aroutdf))
  aroutdf[, paste("ERalphaassociated", ogaGeneList, "BHadj_and_AgeDependentp", sep = "_")] <- rep(NA, nrow(aroutdf))

  aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")] <- ( aroutdf$ProbeId %in% annodf$ProbeId[ogaidxs] )
  
### aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
###   p.adjust(aroutdf[aroutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
  aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("ERalphaassociated", ogaGeneList, "LE60_BHadj_pval", sep = "_")] <-
            p.adjust(aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$LE60_pval, method="BH")
  aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("ERalphaassociated", ogaGeneList, "GT60_BHadj_pval", sep = "_")] <-
            p.adjust(aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$GT60_pval, method="BH")
  aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("ERalphaassociated", ogaGeneList, "AllAges_BHadj_pval", sep = "_")] <-
            p.adjust(aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$AllAges_pval, method="BH")

### aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
###   ( ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
###      ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < 0.05 ) ) |
###     ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
###      ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < 0.05 ) ) |
###     ( ( abs(aroutdf[aroutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
###      ( aroutdf[aroutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < 0.05 ) ) )
  aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")],
          paste("ERalphaassociated", ogaGeneList, "BHadj_and_AgeDependentp", sep = "_")] <-
            ( ( ( abs(aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$LE60_slope) > biosigslopeLE60 )       &
               ( aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")],
                         paste("ERalphaassociated", ogaGeneList, "LE60_BHadj_pval", sep = "_")] < 0.05 ) ) |
             ( ( abs(aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$GT60_slope) > biosigslopeGT60 )       &
              ( aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")],
                        paste("ERalphaassociated", ogaGeneList, "GT60_BHadj_pval", sep = "_")] < 0.05 ) ) |
             ( ( abs(aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")], ]$AllAges_slope) > biosigslopeAllAges ) &
              ( aroutdf[aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")],
                        paste("ERalphaassociated", ogaGeneList, "AllAges_BHadj_pval", sep = "_")] < 0.05 ) ) )

  sogadf <- cbind(oga.rzs, cldvdf[match(dimnames(oga.rzs)[[1]], cldvdf$MBid), ])
  sogardf <- cbind(Geneset.r.mean, oga.r, cldvdf[match(dimnames(oga.r)[[1]], cldvdf$MBid), ])
  names(sogadf)[seq(length(names(ogazidxs)))] <-
    paste(ogadf[match(names(ogaidxs), ogadf$GeneName), ], names(sogadf)[seq(length(names(ogazidxs)))], sep = "|")
  names(sogardf)[seq(length(names(ogazidxs))) + 1] <-
    paste(ogadf[match(names(ogaidxs), ogadf$GeneName), ], names(sogardf)[seq(length(names(ogazidxs))) + 1], sep = "|r|")


### Regression for age <= 60, age > 60   i <- 10; j <- 4
  wclfitageLE60 <- vector("list", length(names(ogazidxs)))
  wclfit <- vector("list", length(names(ogazidxs)))
  stj <- paste(dim(sogardf)[1], "cases")

  if ( !file.exists("./Plots/ERalphaassociated") ) dir.create("./Plots/ERalphaassociated", recursive = TRUE)

  ogaGeneListOutFilename <- paste("AgeRelated_raw_All1992Cases_", ogaGeneList, "_lm_v04.pdf", sep = "")
  pdf(file = paste("./Plots/ERalphaassociated/", ogaGeneListOutFilename, sep = ""), width = 8, height = 10)
### AgeDependent_and_BHadj_FC1.25_ProbeIds <-
###   aroutdf[which(aroutdf$arGeneSet_BHadj_and_AgeDependentp), "Probe_id"]
  AgeDependent_and_BHadj_FC1.25_ProbeIds <-
    aroutdf[which(aroutdf[, paste("ERalphaassociated", ogaGeneList, "BHadj_and_AgeDependentp", sep = "_")]), "Probe_id"]

  xlims <- c(20, 100)
  ylims <- c(4, 17)
  par(mfrow = c(2, 2))
  ageLE60idxp <- (sogardf$age_at_diagnosis <= 60)

  for ( ni in seq( along = c("Geneset.r.mean", names(ogazidxs) ) ) ) {
    i <- c(1, 1 + order( toupper(names(sogardf)[1 + seq(length(names(ogazidxs)))]) ) )[ni]
    cat("\n\n ### --- ", names(sogardf)[i])
    lmifitageLE60 <- lm(sogardf[ageLE60idxp, i] ~ sogardf[ageLE60idxp, "age_at_diagnosis"])
    lmifitageGT60 <- lm(sogardf[!ageLE60idxp, i] ~ sogardf[!ageLE60idxp, "age_at_diagnosis"])
    smrylmifitageLE60 <- summary(lmifitageLE60)
    smrylmifitageGT60 <- summary(lmifitageGT60)
    wclfitageLE60[[i]] <- lmifitageLE60
    print(smrylmifitageLE60)
    print(smrylmifitageGT60)
    lmifit <- lm(sogardf[, i] ~ sogardf[, "age_at_diagnosis"])
    smrylmifit <- summary(lmifit)
    wclfit[[i]] <- lmifit
    print(smrylmifit)
    plot(sogardf[, "age_at_diagnosis"], sogardf[, i],
         ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
         main = paste(stj, names(sogardf)[i], sep = ": "), xlim = xlims, ylim = ylims, pch = ".")
    lines(supsmu(sogardf[, "age_at_diagnosis"], sogardf[, i], span = 0.4, bass = 10), lwd = 3)
    abline(h = 0, v = 50, lty = 2)
    text(x = 20, y = 16.3, labels = paste("[<=60] p = ", format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
         adj = 0, cex = 0.85)
    text(x = 22, y = 15.5, labels = paste("[All] p = ", format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
         adj = 0, cex = 0.85)
    if ( i == 1 ) {
      ##     biosigslopeLE60    biosigslopeGT60     biosigslopeAllAges
      nOtherGenesetAssociations <- 25
      if ( ( ( abs(smrylmifitageLE60$coefficients[2, 1]) > biosigslopeLE60 ) &
            ( smrylmifitageLE60$coefficients[2, 4] < 0.05/nOtherGenesetAssociations ) ) |
          ( ( abs(smrylmifitageGT60$coefficients[2, 1]) > biosigslopeGT60 ) &
           ( smrylmifitageGT60$coefficients[2, 4] < 0.05/nOtherGenesetAssociations ) ) |
          ( ( abs(smrylmifit$coefficients[2, 1]) > biosigslopeAllAges ) &
           ( smrylmifit$coefficients[2, 4] < 0.05/nOtherGenesetAssociations ) ) ) {
        text(x = 98, y = 16.3, labels="Age-dependent trend:", adj = 1, cex = 0.85)
        text(x = 98, y = 15.5, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
      } else {
        text(x = 98, y = 16.3, labels="No detectable trend:", adj = 1, cex = 0.85)
        text(x = 98, y = 15.5, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
      }
    }
    if ( strsplit(names(sogardf)[i], split = "\\|")[[1]][3] %in% AgeDependent_and_BHadj_FC1.25_ProbeIds ) {
      text(x = 98, y = 16.3, labels="Age-dependent trend:", adj = 1, cex = 0.85)
      text(x = 98, y = 15.5, labels="|FC|>1.25 and adjPval<0.05", adj = 1, cex = 0.85)
    } else {
      text(x = 98, y = 16.3, labels="No detectable trend:", adj = 1, cex = 0.85)
      text(x = 98, y = 15.5, labels="|FC|<1.25 or adjPval>0.05", adj = 1, cex = 0.85)
    }
  }
  dev.off()
### [105] "Ontology_Component_0" [106] "Ontology_Process_0"  [107] "Ontology_Function_0"
  outcolidxs <- c(1:104, 108:145, which(grepl(ogaGeneList, names(aroutdf))) )
  write.csv(aroutdf[which(aroutdf[, paste("ERalphaassociated", ogaGeneList, "GeneSetidxp", sep = "_")]), outcolidxs],
            file = paste("./Plots/ERalphaassociated/AgeRelated_raw_All1992Cases_", ogaGeneList, "_lm_v04.csv", sep = ""))
}


NULL
###-------------------------------------------------------------------------------------------------------------------
### End ERalphaassociated gene sets to average over



NULL
###-------------------------------------------------------------------------------------------------------------------
### End Other gene sets to average over


### TODO: 20151102
### IntClust analyses to parallel PAM50 findings
### Set up loop, or convert PAM50 types to IntClust types below.

### Do for each PAM50 subtype:
sm_manhattan <- 
  function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", gsym = "Gene_symbol",
###             col = c("gray10",  "gray60"),
            col = c("black",  "red3", "green4", "blue", "cyan3", "magenta3", "orange", "gray40"),
            chrlabs = NULL, suggestiveline = -log10(1e-05), 
            genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, plotpointidxp = rep(TRUE, nrow(x)),
            plotpointhiliteidxp = rep(FALSE, nrow(x)),  hilitelbls = NULL,  ...) 
{
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

################################################################################ xxxxx zzzzz
###
### BrCaEH
###

### ER+/HER2-  ER+/HER2+  ER-/HER2+  ER-/HER2-  EpHn  EpHp  EnHp  EnHn
### Redo the PAM50 and iClust exercise
### Cull subtypes, run regressions, . . .

sehdf$BrCaEH <- rep(NA_character_, nrow(sehdf))
sehdf[sehdf$ER.Expr == "+" & sehdf$Her2.Expr == "+", "BrCaEH"] <- "ER+/HER2+"
sehdf[sehdf$ER.Expr == "+" & sehdf$Her2.Expr == "-", "BrCaEH"] <- "ER+/HER2-"
sehdf[sehdf$ER.Expr == "-" & sehdf$Her2.Expr == "+", "BrCaEH"] <- "ER-/HER2+"
sehdf[sehdf$ER.Expr == "-" & sehdf$Her2.Expr == "-", "BrCaEH"] <- "ER-/HER2-"
with(sehdf, table(Pam50Subtype, BrCaEH, useNA = "always")) ## 1% difference in classification
sehdf$BrCaEHf <- factor(sehdf$BrCaEH, levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))

sehnovdf$BrCaEH <- rep(NA_character_, nrow(sehnovdf))
sehnovdf[sehnovdf$ER.Expr == "+" & sehnovdf$Her2.Expr == "+", "BrCaEH"] <- "ER+/HER2+"
sehnovdf[sehnovdf$ER.Expr == "+" & sehnovdf$Her2.Expr == "-", "BrCaEH"] <- "ER+/HER2-"
sehnovdf[sehnovdf$ER.Expr == "-" & sehnovdf$Her2.Expr == "+", "BrCaEH"] <- "ER-/HER2+"
sehnovdf[sehnovdf$ER.Expr == "-" & sehnovdf$Her2.Expr == "-", "BrCaEH"] <- "ER-/HER2-"
with(sehnovdf, table(Pam50Subtype, BrCaEH, useNA = "always")) ## 1% difference in classification
sehnovdf$BrCaEHf <- factor(sehnovdf$BrCaEH, levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))

sehrdf$BrCaEH <- rep(NA_character_, nrow(sehdf))
sehrdf[sehrdf$ER.Expr == "+" & sehrdf$Her2.Expr == "+", "BrCaEH"] <- "ER+/HER2+"
sehrdf[sehrdf$ER.Expr == "+" & sehrdf$Her2.Expr == "-", "BrCaEH"] <- "ER+/HER2-"
sehrdf[sehrdf$ER.Expr == "-" & sehrdf$Her2.Expr == "+", "BrCaEH"] <- "ER-/HER2+"
sehrdf[sehrdf$ER.Expr == "-" & sehrdf$Her2.Expr == "-", "BrCaEH"] <- "ER-/HER2-"
with(sehrdf, table(BrCa4c, BrCaEH, useNA = "always")) ## 1% difference in classification
sehrdf$BrCaEHf <- factor(sehrdf$BrCaEH, levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))

###                    BrCaEH
### BrCa4c              ER-/HER2- ER-/HER2+ ER+/HER2- ER+/HER2+ <NA>
###   ER-PR-HER2-             320         0         0         0    0
###   ER-PR-HER2+               0       134         0         0    0
###   ER+ or PR+, HER2-        15         0      1408         0    0
###   ER+ or PR+, HER2+         0         5         0       110    0
###   <NA>                      0         0         0         0    0
### 

with(sehrdf, table(Pam50Subtype, BrCaEHf, useNA = "always"))
with(sehrdf, table(Pam50Subtype, BrCaEH, useNA = "always"))

###             BrCaEH
### Pam50Subtype ER-/HER2- ER-/HER2+ ER+/HER2- ER+/HER2+ <NA>
###       Basal        263        27        38         3    0
###       Her2          40       102        65        33    0
###       LumA           3         1       695        22    0
###       LumB           0         0       447        45    0
###       NC             0         0         6         0    0
###       Normal        29         9       157         7    0
###       <NA>           0         0         0         0    0

## approx Line 6840 ### ar set:  LumA from all 1992 cases

sarrdf$BrCaEH <- rep(NA_character_, nrow(sehdf))
sarrdf[sarrdf$ER.Expr == "+" & sarrdf$Her2.Expr == "+", "BrCaEH"] <- "ER+/HER2+"
sarrdf[sarrdf$ER.Expr == "+" & sarrdf$Her2.Expr == "-", "BrCaEH"] <- "ER+/HER2-"
sarrdf[sarrdf$ER.Expr == "-" & sarrdf$Her2.Expr == "+", "BrCaEH"] <- "ER-/HER2+"
sarrdf[sarrdf$ER.Expr == "-" & sarrdf$Her2.Expr == "-", "BrCaEH"] <- "ER-/HER2-"
with(sarrdf, table(BrCa4c, BrCaEH, useNA = "always")) ## 'BrCa4c' not found   See 1% difference in classification above
with(sarrdf, table(Pam50Subtype, BrCaEH, useNA = "always"))  ### Same table as above.
sarrdf$BrCaEHf <- factor(sarrdf$BrCaEH, levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))

sarrnovdf$BrCaEH <- rep(NA_character_, nrow(sehdf))
sarrnovdf[sarrnovdf$ER.Expr == "+" & sarrnovdf$Her2.Expr == "+", "BrCaEH"] <- "ER+/HER2+"
sarrnovdf[sarrnovdf$ER.Expr == "+" & sarrnovdf$Her2.Expr == "-", "BrCaEH"] <- "ER+/HER2-"
sarrnovdf[sarrnovdf$ER.Expr == "-" & sarrnovdf$Her2.Expr == "+", "BrCaEH"] <- "ER-/HER2+"
sarrnovdf[sarrnovdf$ER.Expr == "-" & sarrnovdf$Her2.Expr == "-", "BrCaEH"] <- "ER-/HER2-"
with(sarrnovdf, table(BrCa4c, BrCaEH, useNA = "always")) ## 'BrCa4c' not found   See 1% difference in classification above
with(sarrnovdf, table(Pam50Subtype, BrCaEH, useNA = "always"))  ### Same table as above.
sarrnovdf$BrCaEHf <- factor(sarrnovdf$BrCaEH, levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))


################################################################################
### ER/HER2 groups

### BrCaEH loop:

### arGeneSet tag denotes 467 Illumina probesets corresponding to the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.
### Otherwise analysis is for all probesets.
### Set up output for age-associated transcript expression trends
### - across whole genome
### - across arGeneSet of 244 genes (467 probesets) for the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.

### BrCaEH i across whole genome Loop Begin   ici <- 0
FDRalphalevels <- c(0.05, 0.01)
Verbosep <- FALSE ## TRUE
DoMemoizep <- FALSE ## TRUE
SingleOutputFilep <- TRUE
arScatterPlotsp <- FALSE ## TRUE
EqAxesScales <- c(FALSE, TRUE)  ## Must be in order F, T so appropriate scale range can be calculated across conditions
nNotAgeDepToPlot <- 1500
propLowDens <- 0.66
propHiDens <- (1.0 - propLowDens)
maxClusterN <- 15000


for ( FDRalphai in seq(along = FDRalphalevels) ) {

    for ( EqAxesScalesi in seq(along = EqAxesScales) ) {

        FDRalpha <- FDRalphalevels[FDRalphai]
        EqAxesScalesp <- EqAxesScales[EqAxesScalesi]

        for ( nBrCaEHi in c( 1, length(levels(sehdf$BrCaEHf)) ) ) {
            if (nBrCaEHi == 1) {
                BrCaEHiAgeAssociatedProbesetsGenesdf <-
                  data.frame(BrCaEHi = "All cases",
                             N_WholeSeries = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_Probesets = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_GeneNames = rep(NA_integer_, nBrCaEHi),
                             N_NoOverlapBigSeries = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_Probesets = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_GeneNames = rep(NA_integer_, nBrCaEHi),

                             N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nBrCaEHi),

                             N_NOv_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nBrCaEHi)
                             )

            } else {
                BrCaEHiAgeAssociatedProbesetsGenesdf <-
                  data.frame(BrCaEHi = levels(sehdf$BrCaEHf),
                             N_WholeSeries = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_Probesets = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_GeneNames = rep(NA_integer_, nBrCaEHi),
                             N_NoOverlapBigSeries = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_Probesets = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_GeneNames = rep(NA_integer_, nBrCaEHi),

                             N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nBrCaEHi),

                             N_NOv_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nBrCaEHi),
                             N_NOv_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nBrCaEHi)
                             )
            }
            
            for ( ici in seq( nBrCaEHi ) ) { ## ici <- 1

                if ( nBrCaEHi == 1 ) {
                    icinm <- "All cases"
                } else {
                    icinm <- levels(sehdf$BrCaEHf)[ici]
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
### BrCaEHinov = BrCaEH group i phenotype, No Overlap with Big Series subset

                ## BrCaEH Subtypes - extract subtype dataframes
                if ( icinm == "All cases" ) {
                    BrCaEHiarrdf   <- sarrdf
                    BrCaEHiarrnovdf   <- sarrnovdf
                    sehBrCaEHidf <- sehdf
                } else {
                    BrCaEHiarrdf   <- sarrdf[sarrdf$BrCaEHf == icinm, ]
                    BrCaEHiarrnovdf   <- sarrnovdf[sarrnovdf$BrCaEHf == icinm, ]
                    sehBrCaEHidf <- sehdf[sehdf$BrCaEHf == icinm, ]
                }
                annodfNames <- names(annodf)
                annodfNamesFirst <- c("Probe_id", "Search_key", "Gene_symbol" )
                annodfNamesLast <- setdiff(annodfNames, annodfNamesFirst)
                arBrCaEHioutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]

                ## Set up matrix to hold regression results for fits to all probes
                arBrCaEHioutmatcolnames <-
                  c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
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
                ## Biologically significant change in expression:  25% increase or decrease in fold change over time
                biosigslopeLE60 <-  log2(1.25)/35
                biosigslopeGT60 <-  log2(1.25)/35
                biosigslopeAllAges <-  log2(1.25)/70
                ## Case ids for this subset
                BrCaEHiidxs <- match(sehBrCaEHidf$MBid, dimnames(Dataset.r)[[1]])
                ## Data frame for regressions
                arBrCaEHilmdf <- data.frame(age_at_diagnosis = sehBrCaEHidf$age_at_diagnosis,
                                            probesetni = Dataset.r[BrCaEHiidxs, 1])
                ## Loop across all probes and fit regression lines - All cases
                cat("\n\n\n") 
                cat(icinm)  
                cat("\n\n\n") 
                ## Begin all 48k
                arBrCaEHioutdf_memoizedfilename <-
                  paste("AgeRelated_AllProbes_",
                        gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                        "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                        "_All1992Cases_lm_v01.csv", sep = "")

                arBrCaEHioutdf_memoizedfilenamep <- FALSE

                if ( (!DoMemoizep) && file.exists(arBrCaEHioutdf_memoizedfilename) ) {
                    arBrCaEHioutdf <- read.table(file = arBrCaEHioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                 sep = ",", header = TRUE, comment.char = "")
                    arBrCaEHioutdf_memoizedfilenamep <- TRUE
                } else {
                    for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
                        psi <- arBrCaEHioutdf$ProbeId[ni]  ## Probeset i
                        drci <- match(psi, dimnames(Dataset.r)[[2]])
                        arBrCaEHilmdf$probesetni <- Dataset.r[BrCaEHiidxs, drci]
                        if (ni %% 100 == 0 ) { cat(ni, ", ") }
                        lmifitageLE60 <- lm(probesetni ~ age_at_diagnosis, data = arBrCaEHilmdf,
                                            na.action = na.exclude, subset = age_at_diagnosis <= 60)
                        lmifitageGT60 <- lm(probesetni ~ age_at_diagnosis, data = arBrCaEHilmdf,
                                            na.action = na.exclude, subset = age_at_diagnosis > 60)
                        lmifitallages <- lm(probesetni ~ age_at_diagnosis, data = arBrCaEHilmdf,
                                            na.action = na.exclude)
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        smrylmifitallages <- summary(lmifitallages)
                        
                        arBrCaEHioutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                        arBrCaEHioutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                        arBrCaEHioutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                        arBrCaEHioutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arBrCaEHioutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arBrCaEHioutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                    }
                    
                    dimnames(arBrCaEHioutmat)[[2]] <- arBrCaEHioutmatcolnames
                    dimnames(arBrCaEHioutmat)[[1]] <- arBrCaEHioutdf$ProbeId
### End all 48k
                    arBrCaEHioutdf[, arBrCaEHioutmatcolnames] <- arBrCaEHioutmat
                    cat("\n\n\n")
                    ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                    arBrCaEHioutdf$LE60_BHadj_pval <- p.adjust(arBrCaEHioutdf$LE60_pval, method="BH")
                    arBrCaEHioutdf$GT60_BHadj_pval <- p.adjust(arBrCaEHioutdf$GT60_pval, method="BH")
                    arBrCaEHioutdf$AllAges_BHadj_pval <- p.adjust(arBrCaEHioutdf$AllAges_pval, method="BH")
                    arBrCaEHioutdf$BHadj_and_AgeDependentp <-
                      ( ( ( abs(arBrCaEHioutdf$LE60_slope) > biosigslopeLE60 )       & ( arBrCaEHioutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arBrCaEHioutdf$GT60_slope) > biosigslopeGT60 )       & ( arBrCaEHioutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arBrCaEHioutdf$AllAges_slope) > biosigslopeAllAges ) & ( arBrCaEHioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arBrCaEHioutdf$BHadj_signifp <- ( ( ( arBrCaEHioutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                        ( arBrCaEHioutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                        ( arBrCaEHioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arBrCaEHioutdf$LE60_log2FC <- arBrCaEHioutdf$LE60_slope * 35
                    arBrCaEHioutdf$GT60_log2FC <- arBrCaEHioutdf$GT60_slope * 35
                    arBrCaEHioutdf$AllAges_log2FC <- arBrCaEHioutdf$AllAges_slope * 70
                    
                    arBrCaEHioutdf$Best_log2FC <- arBrCaEHioutdf$AllAges_log2FC
                    arBrCaEHioutdf$Best_BHadj_pval <- arBrCaEHioutdf$AllAges_BHadj_pval
                    
                    ## Find the largest significant fold change:

                    LE60_Bestp <- ( ( abs(arBrCaEHioutdf$LE60_slope) > biosigslopeLE60 )          &
                                    ( arBrCaEHioutdf$LE60_BHadj_pval < FDRalpha ) )
                    arBrCaEHioutdf[LE60_Bestp, ]$Best_log2FC <- arBrCaEHioutdf[LE60_Bestp, ]$LE60_log2FC
                    arBrCaEHioutdf[LE60_Bestp, ]$Best_BHadj_pval <- arBrCaEHioutdf[LE60_Bestp, ]$LE60_BHadj_pval
                    
                    GT60_Bestp <- ( ( abs(arBrCaEHioutdf$GT60_slope) > biosigslopeGT60 )          &
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
                }

                ## No overlap 1161 cases:  No Overlap with Big Series and MB09 TMA cases
                ## BrCaEHinov = Pam50 BrCaEHi phenotype, from No Overlap with Big Series and MB09 TMA cases
                if ( icinm == "All cases" ) {
                    sehBrCaEHinovdf <- sehnovdf
                } else {
                    sehBrCaEHinovdf <- sehnovdf[sehnovdf$BrCaEHf == icinm, ]
                }
                arBrCaEHinovoutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]
                arBrCaEHinovoutmatcolnames <-
                  c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
                arBrCaEHinovoutmat <- matrix(NA_real_, nrow = nrow(arBrCaEHinovoutdf), ncol = length(arBrCaEHinovoutmatcolnames))
                LE60_slope_col <- match("LE60_slope", arBrCaEHinovoutmatcolnames)
                GT60_slope_col <- match("GT60_slope", arBrCaEHinovoutmatcolnames)
                AllAges_slope_col <- match("AllAges_slope", arBrCaEHinovoutmatcolnames)
                LE60_pval_col <- match("LE60_pval", arBrCaEHinovoutmatcolnames)
                GT60_pval_col <- match("GT60_pval", arBrCaEHinovoutmatcolnames)
                AllAges_pval_col <- match("AllAges_pval", arBrCaEHinovoutmatcolnames)
                AgeDependentp_col <- match("AgeDependentp", arBrCaEHinovoutmatcolnames)
                arBrCaEHinovoutdf$LE60_slope <- rep(NA_real_ , nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$GT60_slope <- rep(NA_real_ , nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$AllAges_slope <- rep(NA_real_ , nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$LE60_pval <- rep(NA_real_ , nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$GT60_pval <- rep(NA_real_ , nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$AllAges_pval <- rep(NA_real_ , nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$AgeDependentp <- rep(NA , nrow(arBrCaEHinovoutdf))
                biosigslopeLE60 <-  log2(1.25)/35
                biosigslopeGT60 <-  log2(1.25)/35
                biosigslopeAllAges <-  log2(1.25)/70
                BrCaEHinovidxs <- match(sehBrCaEHinovdf$MBid, dimnames(Dataset.r)[[1]])
                arBrCaEHinovlmdf <- data.frame(age_at_diagnosis = sehBrCaEHinovdf$age_at_diagnosis,
                                               probesetni = Dataset.r[BrCaEHinovidxs, 1])

                cat("\n\n\n")
                cat(icinm)
                cat("\n\n\n")
                arBrCaEHinovoutdf_memoizedfilename <-
                  paste("AgeRelated_AllProbes_",
                        gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                        "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                        "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "")
                arBrCaEHinovoutdf_memoizedfilenamep <- FALSE

                if ( (!DoMemoizep) && file.exists(arBrCaEHinovoutdf_memoizedfilename) ) {
                    arBrCaEHinovoutdf <- read.table(file = arBrCaEHinovoutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                    sep = ",", header = TRUE, comment.char = "")
                    arBrCaEHinovoutdf_memoizedfilenamep <- TRUE
                } else {
                    for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
                        psi <- arBrCaEHinovoutdf$ProbeId[ni]
                        drci <- match(psi, dimnames(Dataset.r)[[2]])
                        arBrCaEHinovlmdf$probesetni <- Dataset.r[BrCaEHinovidxs, drci]
                        if (ni %% 100 == 0 ) { cat(ni, ", ") }
                        lmifitageLE60 <-
                          lm(probesetni ~ age_at_diagnosis, data = arBrCaEHinovlmdf, na.action = na.exclude, subset = age_at_diagnosis <= 60)
                        lmifitageGT60 <-
                          lm(probesetni ~ age_at_diagnosis, data = arBrCaEHinovlmdf, na.action = na.exclude, subset = age_at_diagnosis > 60)
                        lmifitallages <-
                          lm(probesetni ~ age_at_diagnosis, data = arBrCaEHinovlmdf, na.action = na.exclude)
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        smrylmifitallages <- summary(lmifitallages)
                        
                        arBrCaEHinovoutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                        arBrCaEHinovoutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                        arBrCaEHinovoutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                        arBrCaEHinovoutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arBrCaEHinovoutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arBrCaEHinovoutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                    }

                    dimnames(arBrCaEHinovoutmat)[[2]] <- arBrCaEHinovoutmatcolnames
                    dimnames(arBrCaEHinovoutmat)[[1]] <- arBrCaEHinovoutdf$ProbeId
                    arBrCaEHinovoutdf[, arBrCaEHinovoutmatcolnames] <- arBrCaEHinovoutmat
                    cat("\n\n\n")
                    ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                    arBrCaEHinovoutdf$LE60_BHadj_pval <- p.adjust(arBrCaEHinovoutdf$LE60_pval, method="BH")
                    arBrCaEHinovoutdf$GT60_BHadj_pval <- p.adjust(arBrCaEHinovoutdf$GT60_pval, method="BH")
                    arBrCaEHinovoutdf$AllAges_BHadj_pval <- p.adjust(arBrCaEHinovoutdf$AllAges_pval, method="BH")
                    arBrCaEHinovoutdf$BHadj_and_AgeDependentp <-
                      ( ( ( abs(arBrCaEHinovoutdf$LE60_slope) > biosigslopeLE60 )       &
                          ( arBrCaEHinovoutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arBrCaEHinovoutdf$GT60_slope) > biosigslopeGT60 )       &
                          ( arBrCaEHinovoutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arBrCaEHinovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                          ( arBrCaEHinovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arBrCaEHinovoutdf$BHadj_signifp <- ( ( ( arBrCaEHinovoutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                           ( arBrCaEHinovoutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                           ( arBrCaEHinovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arBrCaEHinovoutdf$LE60_log2FC <- arBrCaEHinovoutdf$LE60_slope * 35
                    arBrCaEHinovoutdf$GT60_log2FC <- arBrCaEHinovoutdf$GT60_slope * 35
                    arBrCaEHinovoutdf$AllAges_log2FC <- arBrCaEHinovoutdf$AllAges_slope * 70

                    arBrCaEHinovoutdf$Best_log2FC <- arBrCaEHinovoutdf$AllAges_log2FC
                    arBrCaEHinovoutdf$Best_BHadj_pval <- arBrCaEHinovoutdf$AllAges_BHadj_pval

                    ## Find the largest significant fold change:

                    LE60_Bestp <- ( ( abs(arBrCaEHinovoutdf$LE60_slope) > biosigslopeLE60 ) &
                                    ( arBrCaEHinovoutdf$LE60_BHadj_pval < FDRalpha ) )
                    arBrCaEHinovoutdf[LE60_Bestp, ]$Best_log2FC <- arBrCaEHinovoutdf[LE60_Bestp, ]$LE60_log2FC
                    arBrCaEHinovoutdf[LE60_Bestp, ]$Best_BHadj_pval <- arBrCaEHinovoutdf[LE60_Bestp, ]$LE60_BHadj_pval

                    GT60_Bestp <- ( ( abs(arBrCaEHinovoutdf$GT60_slope) > biosigslopeGT60 )          &
                                    ( arBrCaEHinovoutdf$GT60_BHadj_pval < FDRalpha )                     &
                                    ( abs(arBrCaEHinovoutdf$GT60_slope) > abs(arBrCaEHinovoutdf$LE60_slope) ) )
                    arBrCaEHinovoutdf[GT60_Bestp, ]$Best_log2FC <- arBrCaEHinovoutdf[GT60_Bestp, ]$GT60_log2FC
                    arBrCaEHinovoutdf[GT60_Bestp, ]$Best_BHadj_pval <- arBrCaEHinovoutdf[GT60_Bestp, ]$GT60_BHadj_pval

                    AllAges_Bestp <- ( ( abs(arBrCaEHinovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                                       ( arBrCaEHinovoutdf$AllAges_BHadj_pval < FDRalpha )               &
                                       ( ( abs(arBrCaEHinovoutdf$AllAges_log2FC) > abs(arBrCaEHinovoutdf$LE60_log2FC) ) |
                                         ( abs(arBrCaEHinovoutdf$AllAges_log2FC) > abs(arBrCaEHinovoutdf$GT60_log2FC) ) ) )
                    arBrCaEHinovoutdf[AllAges_Bestp, ]$Best_log2FC <- arBrCaEHinovoutdf[AllAges_Bestp, ]$AllAges_log2FC
                    arBrCaEHinovoutdf[AllAges_Bestp, ]$Best_BHadj_pval <- arBrCaEHinovoutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
                    arBrCaEHinovoutdf$Abs_Best_log2FC <- abs( arBrCaEHinovoutdf$Best_log2FC )
                    arBrCaEHinovoutdf$Abs_FoldChange <- 2^arBrCaEHinovoutdf$Abs_Best_log2FC
                    arBrCaEHinovoutdf$FoldChange_Direction <- ifelse(arBrCaEHinovoutdf$Best_log2FC > 0, "Up", "Down")
                }

### Use this Best_BHadj_pval for manhattan plots as well.

### Plot scatterplots of gene sets:
### ar = age related study general gene set of interest for Tomo
### p  = polycomb independent
### prc2

### ar set:  BrCaEHi from all 1992 cases
### BrCaEHiarrdf <- sarrdf[sarrdf$BrCaEHf == icinm, ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                arBrCaEHioutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(arBrCaEHioutdf))
                arBrCaEHioutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(arBrCaEHioutdf))
                arBrCaEHioutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(arBrCaEHioutdf))
                arBrCaEHioutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(arBrCaEHioutdf))
                arBrCaEHioutdf$arGeneSet_PBHadj <- rep(NA, nrow(arBrCaEHioutdf))
                
                arBrCaEHioutdf$arGeneSetidxp <- ( arBrCaEHioutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                
                arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                  p.adjust(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                  p.adjust(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                  p.adjust(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                  apply(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp,
                                       c("arGeneSet_LE60_BHadj_pval", "arGeneSet_GT60_BHadj_pval", "arGeneSet_AllAges_BHadj_pval")], 1, min)


                arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                  ( ( ( abs(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
                      ( arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < FDRalpha ) ) |
                    ( ( abs(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
                      ( arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < FDRalpha ) ) |
                    ( ( abs(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                      ( arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                


### BrCaEHi All 1992 cases
### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in arBrCaEHioutdf
                BrCaEHiclfitageLE60 <- vector("list", length(names(arzidxs)))
                BrCaEHiclfitageGT60 <- vector("list", length(names(arzidxs)))
                BrCaEHiclfit <- vector("list", length(names(arzidxs)))
                stj <- paste("N =", dim(BrCaEHiarrdf)[1])

                ## Scatterplots for EZH2 H3K27me3 pathway genes arGeneSet
                if ( SingleOutputFilep && arScatterPlotsp ) {
                    pdf(file = paste("./Plots/MBEX_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                        width = 8, height = 10, useDingbats = FALSE)
                    par(mfrow = c(2, 2))
                }

                xlims <- c(20, 100)
                ylims <- c(4, 17)

                ageLE60idxp <- (BrCaEHiarrdf$age_at_diagnosis <= 60)
                ageGT60idxp <- (!ageLE60idxp)
                for ( ni in seq( along = names(arzidxs) ) ) {

                    i <- order(toupper(names(BrCaEHiarrdf)[seq(length(names(arzidxs)))]))[ni]
                    cat("\n\n### --- ", names(BrCaEHiarrdf)[i])
                    if ( !SingleOutputFilep && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/ProbeLevel/MBEX_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_", gsub("/", "_", names(BrCaEHiarrdf)[i]),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                            width = 6, height = 7, useDingbats = FALSE)
                        par(mfrow = c(1, 1))
                    }
                    lmifitageLE60 <- lm(BrCaEHiarrdf[ageLE60idxp, i] ~ BrCaEHiarrdf[ageLE60idxp, "age_at_diagnosis"])
                    smrylmifitageLE60 <- summary(lmifitageLE60)
                    BrCaEHiclfitageLE60[[i]] <- lmifitageLE60
                    if ( Verbosep ) print(smrylmifitageLE60)
                    lmifitageGT60 <- lm(BrCaEHiarrdf[ageGT60idxp, i] ~ BrCaEHiarrdf[ageGT60idxp, "age_at_diagnosis"])
                    smrylmifitageGT60 <- summary(lmifitageGT60)
                    BrCaEHiclfitageGT60[[i]] <- lmifitageGT60
                    if ( Verbosep ) print(smrylmifitageGT60)
                    lmifit <- lm(BrCaEHiarrdf[, i] ~ BrCaEHiarrdf[, "age_at_diagnosis"])
                    smrylmifit <- summary(lmifit)
                    BrCaEHiclfit[[i]] <- lmifit
                    if ( Verbosep ) print(smrylmifit)
                    if ( arScatterPlotsp ) {
                        plot(BrCaEHiarrdf[, "age_at_diagnosis"], BrCaEHiarrdf[, i],
                             ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                             main = "", xlim = xlims, ylim = ylims, type = "n")
                        title( main = paste(icinm, stj, sep = ": "), line = 2)
                        title( main =
                                 paste(names(BrCaEHiarrdf)[i], "(",
                                       ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(BrCaEHiarrdf)[i],
                                                                                   split = "\\|")[[1]][2], "ERbinding"],
                                              "ER binding", "non-ER binding"), ")"), line = 1)
                        points(BrCaEHiarrdf[, "age_at_diagnosis"], BrCaEHiarrdf[, i],
                               pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                        lines(supsmu(BrCaEHiarrdf[, "age_at_diagnosis"], BrCaEHiarrdf[, i], span = 0.4, bass = 10),
                              lwd = 3, col = "#AA1010")
                        abline(h = 0, v = 50, lty = 2)
                        text(x = 20, y = 16.3,
                             labels = paste("[<=60] p = ",
                                            format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 20, y = 15.5,
                             labels = paste("[>60] p = ",
                                            format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 98, y = 16.3,
                             labels = paste("[All] p = ",
                                            format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 1, cex = 0.85)
                        ## Need BrCaEHi subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                        if ( strsplit(names(BrCaEHiarrdf)[i], split = "\\|")[[1]][2] %in%
                             arBrCaEHioutdf[arBrCaEHioutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                            text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|>1.25 and adjPval<",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        } else {
                            text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels=paste("|FC|<1.25 or adjPval>",
                                              format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        }
                    }
                    if (!SingleOutputFilep && arScatterPlotsp ) { dev.off() }
                }
                if (SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arBrCaEHioutdf) %in% skipVarNms )
                if ( DoMemoizep || (! arBrCaEHioutdf_memoizedfilenamep) ) {
                    write.csv(arBrCaEHioutdf[, -c(skipVarCols)], file = arBrCaEHioutdf_memoizedfilename )
                    write.csv(arBrCaEHioutdf[arBrCaEHioutdf$arGeneSetidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_arGeneSet_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_All1992Cases_lm_v01.csv", sep = "") )
                }


### ar set:  BrCaEHi from 1161 no overlap cases (shown in asudf$MBid)
### BrCaEHiarrnovdf <- sarrdf[rownames(sarrdf) %in% asudf$MBid & sarrdf$Pam50Subtype == "BrCaEHi", ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                arBrCaEHinovoutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(arBrCaEHinovoutdf))
                arBrCaEHinovoutdf$arGeneSet_PBHadj <- rep(NA, nrow(arBrCaEHinovoutdf))
                
                arBrCaEHinovoutdf$arGeneSetidxp <- ( arBrCaEHinovoutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                
                arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                  p.adjust(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                  p.adjust(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                  p.adjust(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                  apply(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp,
                                          c("arGeneSet_LE60_BHadj_pval",
                                            "arGeneSet_GT60_BHadj_pval",
                                            "arGeneSet_AllAges_BHadj_pval")], 1, min)

                
                arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                  ( ( ( abs(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$LE60_slope)    > biosigslopeLE60 )    &
                      ( arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval    < FDRalpha ) ) |
                    ( ( abs(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$GT60_slope)    > biosigslopeGT60 )    &
                      ( arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval    < FDRalpha ) ) |
                    ( ( abs(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                      ( arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                

### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in arBrCaEHinovoutdf
                BrCaEHiclfitageLE60 <- vector("list", length(names(arzidxs)))
                BrCaEHiclfitageGT60 <- vector("list", length(names(arzidxs)))
                BrCaEHiclfit <- vector("list", length(names(arzidxs)))
                stj <- paste("N =", dim(BrCaEHiarrnovdf)[1], "(No overlap)")
                
                if ( SingleOutputFilep && arScatterPlotsp ) {
                    pdf(file = paste("./Plots/MBEX_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                        width = 8, height = 10, useDingbats = FALSE)
                    par(mfrow = c(2, 2))
                }
                xlims <- c(20, 100)
                ylims <- c(4, 17)

                ageLE60idxp <- (BrCaEHiarrnovdf$age_at_diagnosis <= 60)
                ageGT60idxp <- (!ageLE60idxp)
                for ( ni in seq( along = names(arzidxs) ) ) {
                    i <- order(toupper(names(BrCaEHiarrnovdf)[seq(length(names(arzidxs)))]))[ni]
                    cat("\n\n### --- ", names(BrCaEHiarrnovdf)[i])
                    
                    if ( (!SingleOutputFilep) && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/ProbeLevel/MBEX_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), 
                                         "_", gsub("/", "_", names(BrCaEHiarrnovdf)[i]), 
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                            width = 6, height = 7, useDingbats = FALSE)
                        par(mfrow = c(1, 1))
                    }

                    lmifitageLE60 <- lm(BrCaEHiarrnovdf[ageLE60idxp, i] ~ BrCaEHiarrnovdf[ageLE60idxp, "age_at_diagnosis"])
                    smrylmifitageLE60 <- summary(lmifitageLE60)
                    BrCaEHiclfitageLE60[[i]] <- lmifitageLE60
                    if ( Verbosep ) print(smrylmifitageLE60)
                    lmifitageGT60 <- lm(BrCaEHiarrnovdf[ageGT60idxp, i] ~ BrCaEHiarrnovdf[ageGT60idxp, "age_at_diagnosis"])
                    smrylmifitageGT60 <- summary(lmifitageGT60)
                    BrCaEHiclfitageGT60[[i]] <- lmifitageGT60
                    if ( Verbosep ) print(smrylmifitageGT60)
                    lmifit <- lm(BrCaEHiarrnovdf[, i] ~ BrCaEHiarrnovdf[, "age_at_diagnosis"])
                    smrylmifit <- summary(lmifit)
                    BrCaEHiclfit[[i]] <- lmifit
                    if ( Verbosep ) print(smrylmifit)
                    if ( arScatterPlotsp ) {
                        plot(BrCaEHiarrnovdf[, "age_at_diagnosis"], BrCaEHiarrnovdf[, i],
                             ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                             main = "", xlim = xlims, ylim = ylims, type = "n")
                        title( main = paste(icinm, stj, sep = ": "), line = 2)
                        title( main =
                                 paste(names(BrCaEHiarrnovdf)[i], "(",
                                       ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(BrCaEHiarrnovdf)[i],
                                                                                   split = "\\|")[[1]][2], "ERbinding"],
                                              "ER binding", "non-ER binding"), ")"), line = 1)
                        points(BrCaEHiarrnovdf[, "age_at_diagnosis"], BrCaEHiarrnovdf[, i],
                               pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                        lines(supsmu(BrCaEHiarrnovdf[, "age_at_diagnosis"], BrCaEHiarrnovdf[, i], span = 0.4, bass = 10),
                              lwd = 3, col = "#AA1010")
                        abline(h = 0, v = 50, lty = 2)
                        text(x = 20, y = 16.3,
                             labels = paste("[<=60] p = ",
                                            format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 20, y = 15.5,
                             labels = paste("[>60] p = ",
                                            format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 98, y = 16.3,
                             labels = paste("[All] p = ",
                                            format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 1, cex = 0.85)
                        ## Need BrCaEHi subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                        if ( strsplit(names(BrCaEHiarrnovdf)[i], split = "\\|")[[1]][2] %in%
                             arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                            text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|>1.25 and adjPval<",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        } else {
                            text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|<1.25 or adjPval>",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        }
                    }
                    if ( (!SingleOutputFilep) && arScatterPlotsp ) { dev.off() }
                }
                if ( SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arBrCaEHinovoutdf) %in% skipVarNms )
                if ( DoMemoizep || (! arBrCaEHinovoutdf_memoizedfilenamep ) ) {
                    write.csv(arBrCaEHinovoutdf[, -c(skipVarCols)],
                              file = arBrCaEHinovoutdf_memoizedfilename )
                    write.csv(arBrCaEHinovoutdf[arBrCaEHinovoutdf$arGeneSetidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_arGeneSet_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
                }



### Plot age-associated p-values by chromosome position - Manhattan plot

### qqman package has manhattan plot

                require("qqman")

### Data frame with CHR, BP, P (and SNP to avoid warnings)
### Genomic_location
### chr2:206352192:206352241:+
                arBrCaEHioutdf$Chrcr <- sapply(strsplit(arBrCaEHioutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                arBrCaEHioutdf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "",
                                                                  gsub("_h2_hap1", "", gsub("_random", "", arBrCaEHioutdf$Chrcr))))
                arBrCaEHioutdf$Chrn <- as.numeric(gsub("chr", "", arBrCaEHioutdf$Chrc))
                arBrCaEHioutdf$Chrn[arBrCaEHioutdf$Chrc == "chrX"] <- 23
                arBrCaEHioutdf$Chrn[arBrCaEHioutdf$Chrc == "chrY"] <- 24
                with(arBrCaEHioutdf, table(Chrcr, Chrn, useNA = "always"))
                with(arBrCaEHioutdf, table(Chrc, Chrn, useNA = "always"))
                arBrCaEHioutdf$CHR <- arBrCaEHioutdf$Chrn
                arBrCaEHioutdf$Startcr <- sapply(strsplit(arBrCaEHioutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                arBrCaEHioutdf$Stopcr <- sapply(strsplit(arBrCaEHioutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                arBrCaEHioutdf$Startn <- as.numeric(arBrCaEHioutdf$Startcr)
                arBrCaEHioutdf$Stopn <- as.numeric(arBrCaEHioutdf$Stopcr)
                arBrCaEHioutdf$BP <- trunc((arBrCaEHioutdf$Startn + arBrCaEHioutdf$Stopn)/2)
                arBrCaEHioutdf$SNP <- arBrCaEHioutdf$Probe_id
                arBrCaEHioutdf$P <- apply(arBrCaEHioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                arBrCaEHioutdf$PBHadj <- apply(arBrCaEHioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                arBrCaEHioutdfarGeneSetManhattanidxp <- arBrCaEHioutdf$arGeneSet_BHadj_and_AgeDependentp & arBrCaEHioutdf$Chrn <= 23
                arBrCaEHioutdfarGeneSetManhattanidxp[is.na(arBrCaEHioutdfarGeneSetManhattanidxp)] <- FALSE
                arBrCaEHioutdfManhattanidxp <- arBrCaEHioutdf$BHadj_and_AgeDependentp & arBrCaEHioutdf$Chrn <= 23
                arBrCaEHioutdfManhattanidxp[is.na(arBrCaEHioutdfManhattanidxp)] <- FALSE
                ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                arBrCaEHioutdfManhattanVarsidxp <-
                  ( apply(arBrCaEHioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                          function(x) !any(is.na(x)))  & ( arBrCaEHioutdf$Chrn <= 23 ) )

                arBrCaEHioutdf[arBrCaEHioutdf$Gene_symbol == "", "Gene_symbol"] <-
                  arBrCaEHioutdf[arBrCaEHioutdf$Gene_symbol == "", "ILMN_Gene_0"]
                ## Number of probesets showing evidence of association with age
                NpsBrCaEHiAssocAge <- sum(arBrCaEHioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets"] <- NpsBrCaEHiAssocAge
                ## Number of probesets showing evidence of association with age with genomic location
                NpsBrCaEHiAssocAgecGL <- sum(arBrCaEHioutdfManhattanidxp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc"] <- NpsBrCaEHiAssocAgecGL
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgBrCaEHiAssocAgecGL <- length(unique(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, "ILMN_Gene_0"]))
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames"] <- NgBrCaEHiAssocAgecGL
                ## arGeneSet subset
                ## Number of probesets showing evidence of association with age
                NpsBrCaEHiAssocAge_arGeneSet <- sum(arBrCaEHioutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_arGeneSet"] <- NpsBrCaEHiAssocAge_arGeneSet
                ## Number of probesets showing evidence of association with age with genomic location
                NpsBrCaEHiAssocAgecGL_arGeneSet <- sum(arBrCaEHioutdfarGeneSetManhattanidxp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <- NpsBrCaEHiAssocAgecGL_arGeneSet
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgBrCaEHiAssocAgecGL_arGeneSet <- length(unique(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames_arGeneSet"] <- NgBrCaEHiAssocAgecGL_arGeneSet

                ## Non-overlapping with Big Series
                arBrCaEHinovoutdf$Chrcr <- sapply(strsplit(arBrCaEHinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                arBrCaEHinovoutdf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "",
                                                                     gsub("_h2_hap1", "", gsub("_random", "", arBrCaEHinovoutdf$Chrcr))))
                arBrCaEHinovoutdf$Chrn <- as.numeric(gsub("chr", "", arBrCaEHinovoutdf$Chrc))
                arBrCaEHinovoutdf$Chrn[arBrCaEHinovoutdf$Chrc == "chrX"] <- 23
                arBrCaEHinovoutdf$Chrn[arBrCaEHinovoutdf$Chrc == "chrY"] <- 24
                with(arBrCaEHinovoutdf, table(Chrcr, Chrn, useNA = "always"))
                with(arBrCaEHinovoutdf, table(Chrc, Chrn, useNA = "always"))
                arBrCaEHinovoutdf$CHR <- arBrCaEHinovoutdf$Chrn
                arBrCaEHinovoutdf$Startcr <- sapply(strsplit(arBrCaEHinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                arBrCaEHinovoutdf$Stopcr <- sapply(strsplit(arBrCaEHinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                arBrCaEHinovoutdf$Startn <- as.numeric(arBrCaEHinovoutdf$Startcr)
                arBrCaEHinovoutdf$Stopn <- as.numeric(arBrCaEHinovoutdf$Stopcr)
                arBrCaEHinovoutdf$BP <- trunc((arBrCaEHinovoutdf$Startn + arBrCaEHinovoutdf$Stopn)/2)
                arBrCaEHinovoutdf$SNP <- arBrCaEHinovoutdf$Probe_id
                arBrCaEHinovoutdf$P <- apply(arBrCaEHinovoutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                arBrCaEHinovoutdf$PBHadj <- apply(arBrCaEHinovoutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                arBrCaEHinovoutdfarGeneSetManhattanidxp <-
                  arBrCaEHinovoutdf$arGeneSet_BHadj_and_AgeDependentp & arBrCaEHinovoutdf$Chrn <= 23
                arBrCaEHinovoutdfarGeneSetManhattanidxp[is.na(arBrCaEHinovoutdfarGeneSetManhattanidxp)] <- FALSE
                arBrCaEHinovoutdfManhattanidxp <- arBrCaEHinovoutdf$BHadj_and_AgeDependentp & arBrCaEHinovoutdf$Chrn <= 23
                arBrCaEHinovoutdfManhattanidxp[is.na(arBrCaEHinovoutdfManhattanidxp)] <- FALSE
### Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                arBrCaEHinovoutdfManhattanVarsidxp <-
                  ( apply(arBrCaEHinovoutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                          function(x) !any(is.na(x)))  & ( arBrCaEHinovoutdf$Chrn <= 23 ) )

                arBrCaEHinovoutdf[arBrCaEHinovoutdf$Gene_symbol == "", "Gene_symbol"] <-
                  arBrCaEHinovoutdf[arBrCaEHinovoutdf$Gene_symbol == "", "ILMN_Gene_0"]
                ## Number of probesets showing evidence of association with age
                NpsBrCaEHinovAssocAge <- sum(arBrCaEHinovoutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets"] <- NpsBrCaEHinovAssocAge
                ## Number of probesets showing evidence of association with age with genomic location
                NpsBrCaEHinovAssocAgecGL <- sum(arBrCaEHinovoutdfManhattanidxp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc"] <- NpsBrCaEHinovAssocAgecGL
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgBrCaEHinovAssocAgecGL <- length(unique(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, "ILMN_Gene_0"]))
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames"] <- NgBrCaEHinovAssocAgecGL
                ## arGeneSet subset
                ## Number of probesets showing evidence of association with age
                NpsBrCaEHinovAssocAge_arGeneSet <- sum(arBrCaEHinovoutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_arGeneSet"] <- NpsBrCaEHinovAssocAge_arGeneSet
                ## Number of probesets showing evidence of association with age with genomic location
                NpsBrCaEHinovAssocAgecGL_arGeneSet <- sum(arBrCaEHinovoutdfarGeneSetManhattanidxp, na.rm = TRUE)
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <-
                  NpsBrCaEHinovAssocAgecGL_arGeneSet
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgBrCaEHinovAssocAgecGL_arGeneSet <-
                  length(unique(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                BrCaEHiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames_arGeneSet"] <-
                  NgBrCaEHinovAssocAgecGL_arGeneSet
                
                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arBrCaEHioutdf) %in% skipVarNms )

                if ( sum(arBrCaEHioutdfManhattanidxp) ) {
                    write.csv(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, -c(skipVarCols)][
                        order(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, "Abs_FoldChange"], decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.csv(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.table(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, ][
                        !duplicated(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, ][
                            !duplicated(arBrCaEHioutdf[arBrCaEHioutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arBrCaEHinovoutdf) %in% skipVarNms )

                if ( sum(arBrCaEHinovoutdfManhattanidxp) ) {
                    write.csv(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, -c(skipVarCols)][
                        order(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.csv(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.table(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, ][
                        !duplicated(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, ][
                            !duplicated(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_outdf_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                ## arGeneSet subset
                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arBrCaEHioutdf) %in% skipVarNms )
                if ( sum ( arBrCaEHioutdfarGeneSetManhattanidxp ) ) {
                    write.csv(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                        order(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.csv(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_arGeneSet_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.table(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, ][
                        !duplicated(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, ][
                            !duplicated(arBrCaEHioutdf[arBrCaEHioutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arBrCaEHinovoutdf) %in% skipVarNms )
                if ( sum ( arBrCaEHinovoutdfarGeneSetManhattanidxp ) ) {
                    write.csv(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                        order(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv",
                                     sep = "") )
                    write.csv(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_arGeneSet_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.table(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, ][
                        !duplicated(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, ][
                               !duplicated(arBrCaEHinovoutdf[arBrCaEHinovoutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0),
                               ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                
                ## Do stacked Manhattan and volcano plot.  Plot log(FC) on X axis, log10(pval) on Y axis

                if ( !SingleOutputFilep ) {
                    pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 8, height = 8, useDingbats = FALSE)
                    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                } else {
                    pdf(file = paste("./Plots/AgeRelated_Manhattan_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                }

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
                
                if ( length( ObjsToHighlight ) ) {
                    sm_manhattan(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj", chrlabs = c(1:22, "X"),
                                 suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                 highlight = ObjsToHighlight,
                                 plotpointidxp = arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp,
                                                                c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp = 
                                   if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                } else {
                    sm_manhattan(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj", chrlabs = c(1:22, "X"),
                                 suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                 plotpointidxp = arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],
                                 plotpointhiliteidxp =
                                   if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arBrCaEHioutdf[arBrCaEHioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                }
                title(paste("METABRIC:", icinm), line = 2)
                title("Expression showing age-related trend, FC > 1.25", line = 1)
                
                write.csv(arBrCaEHioutdf[which(arBrCaEHioutdf$SNP %in% ObjsToHighlight), ],
                          file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_All1992Cases.csv", sep = "") )

                if ( SingleOutputFilep ) {
                    dev.off()
                    pdf(file = paste("./Plots/AgeRelated_Volcano_",
                                 if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   

                }

                ## Volcano

                stairstepvals <-
                  switch(as.character(FDRalpha),
                         "0.05" = switch(as.character(icinm),
                                         "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                                         "ER+/HER2-" = c(log2(3), 10, log2(2.5), 20, log2(1.25), 25), 
                                         "ER+/HER2+" = c(log2(4), 10, log2(3), 25, log2(2), 40), ## Nothing to label
                                         "ER-/HER2+" = c(log2(4), 10, log2(3), 25, log2(2), 40), ## Nothing to label
                                         "ER-/HER2-" = c(log2(2.25), 1.5, log2(1.75), 1.75, log2(1.25), 2)
                                         ),
                         "0.01" = switch(as.character(icinm),
                                         "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                                         "ER+/HER2-" = c(log2(3), 10, log2(2.5), 20, log2(1.25), 25), 
                                         "ER+/HER2+" = c(log2(4), 10, log2(3), 25, log2(2), 40), ## Nothing to label
                                         "ER-/HER2+" = c(log2(4), 10, log2(3), 25, log2(2), 40), ## Nothing to label
                                         "ER-/HER2-" = c(log2(2.25),  -log10(FDRalpha),
                                                         log2(1.75),  -log10(FDRalpha), log2(1.25),  -log10(FDRalpha))
                                         )
                         )

                volcanoaddlidxp <- ( ( ( abs(arBrCaEHioutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                       ( -log10(arBrCaEHioutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                     ( ( abs(arBrCaEHioutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                       ( -log10(arBrCaEHioutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                     ( ( abs(arBrCaEHioutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                       ( -log10(arBrCaEHioutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                ObjsToHighlight <- arBrCaEHioutdf[volcanoaddlidxp, "SNP"]
                
                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                mainTitle <- paste("METABRIC BrCa", icinm)

                ERbindingAgeDep <- which( ( arBrCaEHioutdf$ERbinding &
                                            (abs(arBrCaEHioutdf$Best_log2FC ) > log2(FCthresh)) &
                                            (arBrCaEHioutdf$Best_BHadj_pval < alphacrit) ) )
                ERnonbindingAgeDep <- which( ( !arBrCaEHioutdf$ERbinding &
                                               (abs(arBrCaEHioutdf$Best_log2FC ) > log2(FCthresh)) &
                                               (arBrCaEHioutdf$Best_BHadj_pval < alphacrit) ) )
                NotAgeDepidxs <- setdiff(seq(nrow(arBrCaEHioutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                mainTitleN <- paste("Cases:", nrow(sehBrCaEHidf),
                                    "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                if ( ! EqAxesScalesp ) {
                    if  (nBrCaEHi == 1) {
                        xlims_vpEVS <- max(abs(range(arBrCaEHioutdf$Best_log2FC, na.rm = TRUE))) * c(-1, 1) * 1.1
                        ylims_vpEVS <- c( 0, 1.1 * max(-log10(arBrCaEHioutdf$Best_BHadj_pval), na.rm = TRUE) )
                    } else {
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylimsi <- c( 0, 1.2 * max(-log10(arBrCaEHioutdf$Best_BHadj_pval), na.rm = TRUE) )
                        ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                        ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                    }
                }
                ## Will need to cull data points at root of volcano - 48000 points too much.
                ## Get density information to thin out volcano plot
                ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                NADn <- length(NotAgeDepidxs)
                ClustRepsnRem <- NADn %% maxClusterN
                ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered
                

                for ( csi in seq(ClustRepsn) ) {
                    cat("   Cluster rep: ", csi, "  ")
                    if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                    csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                
                    densityEst <-
                      hclust(dist(cbind(arBrCaEHioutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                        jitter(-log10(arBrCaEHioutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                     method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                    lowDensityiidxs <-
                      (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityiidxs <-
                      sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                    if ( csi == 1 ) {
                        lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                        hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                    } else {
                        lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                        hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                    }
                }
                    
                plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                
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
                                c( 0, 1.2 * max(-log10(arBrCaEHioutdf$Best_BHadj_pval), na.rm = TRUE) )
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
                legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                VolcanoAddlObjsToHighlight <- c()
                evcp <- arBrCaEHioutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                ## Cull out duplicate gene names
                evcidx <- unlist(tapply(which(evcp), arBrCaEHioutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                evcERbidx <- evcidx[which(arBrCaEHioutdf[evcidx, "ERbinding"])]
                evcERnbidx <- setdiff(evcidx, evcERbidx)
                lblcol <- rep("red", length(evcidx))
                lblcol[evcidx %in% evcERnbidx] <- "blue"
                if ( length(evcidx)  ) {
                    require("maptools"); ## for pointLabel
                    pointLabel(arBrCaEHioutdf[evcidx, ]$Best_log2FC, -log10(arBrCaEHioutdf[evcidx, ]$Best_BHadj_pval),
                               labels = arBrCaEHioutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                }
                ## Save ojects labeled for labeling in TCGA plots
                write.csv(arBrCaEHioutdf[evcidx, ],
                          file = paste("../Data/AgeRelated_Volcano_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       gsub(" ", "",
                                            paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                           2^stairstepvals[3], stairstepvals[4],
                                                           2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_All1992Cases.csv", sep = "") )
                dev.off()

                ## End All1992Cases

                ## No overlap with Big Series or MB09
                ## Do stacked Manhattan and volcano plot.  Plot log(FC) on X axis, log10(pval) on Y axis
                if ( !SingleOutputFilep ) {
                    pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 8, height = 8, useDingbats = FALSE)
                    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                } else {
                    pdf(file = paste("./Plots/AgeRelated_Manhattan_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                }

                ## Label the top 20 largest fold change with age association and adjusted p-val significance
                ohdf <- arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, ][
                    order(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, "Abs_FoldChange"], decreasing = TRUE), ]
                ohdf <- ohdf[ohdf$BHadj_and_AgeDependentp, ]
                
                ObjsToHighlight <-
                  unique(c(ohdf[ohdf$FoldChange_Direction == "Up", ][1:10, "SNP"],
                           ohdf[ohdf$FoldChange_Direction == "Down", ][1:10, "SNP"],
                           ohdf[1:20, "SNP"]
                           ))

                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]
                
                if ( length( ObjsToHighlight ) ) {
                    sm_manhattan(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj",
                                 chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha), genomewideline = FALSE, highlight = ObjsToHighlight,
                                 plotpointidxp = arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp =
                                   if(sum(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                       arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                } else {
                    sm_manhattan(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj",
                                 chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha), genomewideline = FALSE, 
                                 plotpointidxp = arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp =
                                   if(sum(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                       arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arBrCaEHinovoutdf[arBrCaEHinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                }
                title(paste("METABRIC:", icinm), line = 2)
                title("Expression showing age-related trend, FC > 1.25", line = 1)
                
                write.csv(arBrCaEHinovoutdf[which(arBrCaEHinovoutdf$SNP %in% ObjsToHighlight), ],
                          file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_NoOverlapMB09BigSeries.csv", sep = "") )

                if ( SingleOutputFilep ) {
                    dev.off()
                    pdf(file = paste("./Plots/AgeRelated_Volcano_",
                                 if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1)) 
                }
                
                ## Volcano
                
                stairstepvals <-
                  switch(as.character(FDRalpha),
                         "0.05" =  switch(as.character(icinm),
                                          "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                          "ER+/HER2-" = c(log2(3), 5, log2(2.5), 7, log2(1.25), 9), 
                                          "ER+/HER2+" = c(log2(4),  -log10(FDRalpha), log2(3), -log10(FDRalpha), log2(2), -log10(FDRalpha)), 
                                          "ER-/HER2+" = c(log2(4), 10, log2(3), 25, log2(2), 40), ## Nothing to label
                                          "ER-/HER2-" = c(log2(2.25), 1.5, log2(1.75), 1.75, log2(1.25), 2) ## Nothing to label
                                          ),
                         "0.01" =  switch(as.character(icinm),
                                          "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                          "ER+/HER2-" = c(log2(3), 5, log2(2.5), 7, log2(1.25), 10), 
                                          "ER+/HER2+" = c(log2(4),  -log10(FDRalpha), log2(3), -log10(FDRalpha), log2(2), -log10(FDRalpha)), 
                                          "ER-/HER2+" = c(log2(4), 10, log2(3), 25, log2(2), 40), ## Nothing to label
                                          "ER-/HER2-" = c(log2(2.25), 2, log2(1.75), 2, log2(1.25), 2) ## Nothing to label
                                          )
                         )
                
                volcanoaddlidxp <- ( ( ( abs(arBrCaEHinovoutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                       ( -log10(arBrCaEHinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                     ( ( abs(arBrCaEHinovoutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                       ( -log10(arBrCaEHinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                     ( ( abs(arBrCaEHinovoutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                       ( -log10(arBrCaEHinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                ObjsToHighlight <- arBrCaEHinovoutdf[volcanoaddlidxp, "SNP"]

                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                mainTitle <- paste("METABRIC BrCa", icinm)

                ERbindingAgeDep <- which( ( arBrCaEHinovoutdf$ERbinding &
                                            (abs(arBrCaEHinovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                            (arBrCaEHinovoutdf$Best_BHadj_pval < alphacrit) ) )
                ERnonbindingAgeDep <- which( ( !arBrCaEHinovoutdf$ERbinding &
                                               (abs(arBrCaEHinovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                               (arBrCaEHinovoutdf$Best_BHadj_pval < alphacrit) ) )
                NotAgeDepidxs <- setdiff(seq(nrow(arBrCaEHinovoutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                mainTitleN <- paste("Cases:", nrow(sehBrCaEHinovdf),
                                    "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                if ( ! EqAxesScalesp ) {
                    if  (nBrCaEHi == 1) {
                        xlims_vpEVS <- max(abs(range(arBrCaEHinovoutdf$Best_log2FC, na.rm = TRUE))) * c(-1, 1) * 1.1
                        ylims_vpEVS <- c( 0, 1.1 * max(-log10(arBrCaEHinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                    } else {
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylimsi <- c( 0, 1.2 * max(-log10(arBrCaEHinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                        ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                        ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                    }
                }
                ## Will need to cull data points at root of volcano - 48000 points too much.
                ## Get density information to thin out volcano plot
                ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                NADn <- length(NotAgeDepidxs)
                ClustRepsnRem <- NADn %% maxClusterN
                ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered
                

                for ( csi in seq(ClustRepsn) ) {
                    cat("   Cluster rep: ", csi, "  ")
                    if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                    csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                
                    densityEst <-
                      hclust(dist(cbind(arBrCaEHinovoutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                        jitter(-log10(arBrCaEHinovoutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                     method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                    lowDensityiidxs <-
                      (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityiidxs <-
                      sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                    if ( csi == 1 ) {
                        lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                        hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                    } else {
                        lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                        hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                    }
                }
                    
                plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                
                plot(arBrCaEHinovoutdf$Best_log2FC, -log10(arBrCaEHinovoutdf$Best_BHadj_pval), type = "n", main = "", 
                     xlim =  if ( EqAxesScalesp ) {
                                 xlims_vpEVS
                             } else {
                                 max(abs(range(arBrCaEHinovoutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                             },
                     xlab = "(Down)        FC        (Up)   ",
                     ylab = "-log10(BH adjusted P-values)", xaxt = "n",
                     ylim = if ( EqAxesScalesp ) {
                                ylims_vpEVS
                            } else {
                                c( 0, 1.2 * max(-log10(arBrCaEHinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                            }
                     )

                title(main = mainTitle, line = 2)
                title(main = mainTitleN, line = 1)
                axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), labels = c(32, 16, 8, 4, 2, 1, 2, 4, 8, 16, 32) )
                
                points(arBrCaEHinovoutdf$Best_log2FC[NotAgeDepidxs][plotNotAgeDepidxs],
                       -log10(arBrCaEHinovoutdf$Best_BHadj_pval[NotAgeDepidxs][plotNotAgeDepidxs]),
                       pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))

                points(arBrCaEHinovoutdf$Best_log2FC[ERbindingAgeDep], -log10(arBrCaEHinovoutdf$Best_BHadj_pval[ERbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("red", alpha = FALSE)), alpha=30, max=255))
                points(arBrCaEHinovoutdf$Best_log2FC[ERnonbindingAgeDep], -log10(arBrCaEHinovoutdf$Best_BHadj_pval[ERnonbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("blue", alpha = FALSE)), alpha=30, max=255))
                legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                VolcanoAddlObjsToHighlight <- c()
                evcp <- arBrCaEHinovoutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                ## Cull out duplicate gene names
                evcidx <- unlist(tapply(which(evcp), arBrCaEHinovoutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                evcERbidx <- evcidx[which(arBrCaEHinovoutdf[evcidx, "ERbinding"])]
                evcERnbidx <- setdiff(evcidx, evcERbidx)
                lblcol <- rep("red", length(evcidx))
                lblcol[evcidx %in% evcERnbidx] <- "blue"
                if ( length(evcidx)  ) {
                    require("maptools"); ## for pointLabel
                    pointLabel(arBrCaEHinovoutdf[evcidx, ]$Best_log2FC, -log10(arBrCaEHinovoutdf[evcidx, ]$Best_BHadj_pval),
                               labels = arBrCaEHinovoutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                }
                ## Save ojects labeled for labeling in TCGA plots
                write.csv(arBrCaEHinovoutdf[evcidx, ],
                          file = paste("../Data/AgeRelated_Volcano_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       gsub(" ", "",
                                            paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                           2^stairstepvals[3], stairstepvals[4],
                                                           2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_NoOverlapMB09BigSeries.csv", sep = "") )
                dev.off()
                ## End No overlap with Big Series or MB09
            }   ## End BrCaEH i loop
        }   ## End  nBrCaEHi loop

        BrCaEHiAgeAssociatedProbesetsGenesdf$N_WholeSeries <- as.vector(table(sehdf$BrCaEHf))
        BrCaEHiAgeAssociatedProbesetsGenesdf$N_NoOverlapBigSeries <- as.vector(table(sehnovdf$BrCaEHf))
        
        ## Write out data on counts of probesets/genes showing age association
        write.csv(BrCaEHiAgeAssociatedProbesetsGenesdf,
                  file = paste("AgeAssociatedProbesetGenes_BrCaEH_All",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               ".csv", sep = "") )
        
### AgeRelated_arGeneSet_ERpHER2m_All1992Cases_lm_v01.csv
### AgeRelated_arGeneSet_ERpHER2m_NoOverlapMB09BigSeries_lm_v01.csv
### Read in all these files, add BrCaEH col, save as one file.
### ici <- 1

        for ( ici in seq(levels(sehdf$BrCaEHf)) ) {
            icinm <- levels(sehdf$BrCaEHf)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$BrCaEHf <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_BrCaEH_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               "_All1992Cases_lm_v01.csv", sep = "") )

        
        for ( ici in seq(levels(sehnovdf$BrCaEHf)) ) {
            icinm <- levels(sehnovdf$BrCaEHf)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$BrCaEHf <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_BrCaEH_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
    }
}

################################################################################
### intClust groups

### iClust loop:

### arGeneSet tag denotes 467 Illumina probesets corresponding to the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.
### Otherwise analysis is for all probesets.
### Set up output for age-associated transcript expression trends
### - across whole genome
### - across arGeneSet of 244 genes (467 probesets) for the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.

### iClust i across whole genome Loop Begin   ici <- 0
FDRalphalevels <- c(0.05, 0.01)
Verbosep <- FALSE ## TRUE
DoMemoizep <- FALSE ## TRUE
SingleOutputFilep <- TRUE
arScatterPlotsp <- FALSE ## TRUE
EqAxesScales <- c(FALSE, TRUE)  ## Must be in order F, T so appropriate scale range can be calculated across conditions
nNotAgeDepToPlot <- 1500
propLowDens <- 0.66
propHiDens <- (1.0 - propLowDens)
maxClusterN <- 15000


for ( FDRalphai in seq(along = FDRalphalevels) ) {

    for ( EqAxesScalesi in seq(along = EqAxesScales) ) {

        FDRalpha <- FDRalphalevels[FDRalphai]
        EqAxesScalesp <- EqAxesScales[EqAxesScalesi]

        for ( niClusti in c( 1, length(levels(sehdf$iClustf)) ) ) { ## niClusti <- 1
            if (niClusti == 1) {
                iClustiAgeAssociatedProbesetsGenesdf <-
                  data.frame(iClusti = "All cases",
                             N_WholeSeries = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),
                             N_NoOverlapBigSeries = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),

                             N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti),

                             N_NOv_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti)
                             )

            } else {
                iClustiAgeAssociatedProbesetsGenesdf <-
                  data.frame(iClusti = levels(sehdf$iClustf),
                             N_WholeSeries = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),
                             N_NoOverlapBigSeries = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),

                             N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                             N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti),

                             N_NOv_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                             N_NOv_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti)
                             )
            }
            
            for ( ici in seq( niClusti ) ) { ## ici <- 1

                if ( niClusti == 1 ) {
                    icinm <- "All cases"
                } else {
                    icinm <- levels(sehdf$iClustf)[ici]
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
### iClustinov = iClust group i phenotype, No Overlap with Big Series subset

                ## iClust Subtypes - extract subtype dataframes
                if ( icinm == "All cases" ) {
                    iClustiarrdf   <- sarrdf
                    iClustiarrnovdf   <- sarrnovdf
                    sehiClustidf <- sehdf
                } else {
                    iClustiarrdf   <- sarrdf[sarrdf$iClustf == icinm, ]
                    iClustiarrnovdf   <- sarrnovdf[sarrnovdf$iClustf == icinm, ]
                    sehiClustidf <- sehdf[sehdf$iClustf == icinm, ]
                }
                annodfNames <- names(annodf)
                annodfNamesFirst <- c("Probe_id", "Search_key", "Gene_symbol" )
                annodfNamesLast <- setdiff(annodfNames, annodfNamesFirst)
                ariClustioutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]

                ## Set up matrix to hold regression results for fits to all probes
                ariClustioutmatcolnames <-
                  c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
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
                ## Biologically significant change in expression:  25% increase or decrease in fold change over time
                biosigslopeLE60 <-  log2(1.25)/35
                biosigslopeGT60 <-  log2(1.25)/35
                biosigslopeAllAges <-  log2(1.25)/70
                ## Case ids for this subset
                iClustiidxs <- match(sehiClustidf$MBid, dimnames(Dataset.r)[[1]])
                ## Data frame for regressions
                ariClustilmdf <- data.frame(age_at_diagnosis = sehiClustidf$age_at_diagnosis,
                                            probesetni = Dataset.r[iClustiidxs, 1])
                ## Loop across all probes and fit regression lines - All cases
                cat("\n\n\n") 
                cat(icinm)  
                cat("\n\n\n") 
                ## Begin all 48k
                ariClustioutdf_memoizedfilename <-
                  paste("AgeRelated_AllProbes_",
                        gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                        "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                        "_All1992Cases_lm_v01.csv", sep = "")

                ariClustioutdf_memoizedfilenamep <- FALSE

                if ( (!DoMemoizep) && file.exists(ariClustioutdf_memoizedfilename) ) {
                    ariClustioutdf <- read.table(file = ariClustioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                 sep = ",", header = TRUE, comment.char = "")
                    ariClustioutdf_memoizedfilenamep <- TRUE
                } else {
                    for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
                        psi <- ariClustioutdf$ProbeId[ni]  ## Probeset i
                        drci <- match(psi, dimnames(Dataset.r)[[2]])
                        ariClustilmdf$probesetni <- Dataset.r[iClustiidxs, drci]
                        if (ni %% 100 == 0 ) { cat(ni, ", ") }
                        lmifitageLE60 <- lm(probesetni ~ age_at_diagnosis, data = ariClustilmdf,
                                            na.action = na.exclude, subset = age_at_diagnosis <= 60)
                        lmifitageGT60 <- lm(probesetni ~ age_at_diagnosis, data = ariClustilmdf,
                                            na.action = na.exclude, subset = age_at_diagnosis > 60)
                        lmifitallages <- lm(probesetni ~ age_at_diagnosis, data = ariClustilmdf,
                                            na.action = na.exclude)
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        smrylmifitallages <- summary(lmifitallages)
                        
                        ariClustioutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                        ariClustioutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                        ariClustioutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                        ariClustioutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        ariClustioutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        ariClustioutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                    }
                    
                    dimnames(ariClustioutmat)[[2]] <- ariClustioutmatcolnames
                    dimnames(ariClustioutmat)[[1]] <- ariClustioutdf$ProbeId
### End all 48k
                    ariClustioutdf[, ariClustioutmatcolnames] <- ariClustioutmat
                    cat("\n\n\n")
                    ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                    ariClustioutdf$LE60_BHadj_pval <- p.adjust(ariClustioutdf$LE60_pval, method="BH")
                    ariClustioutdf$GT60_BHadj_pval <- p.adjust(ariClustioutdf$GT60_pval, method="BH")
                    ariClustioutdf$AllAges_BHadj_pval <- p.adjust(ariClustioutdf$AllAges_pval, method="BH")
                    ariClustioutdf$BHadj_and_AgeDependentp <-
                      ( ( ( abs(ariClustioutdf$LE60_slope) > biosigslopeLE60 )       & ( ariClustioutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(ariClustioutdf$GT60_slope) > biosigslopeGT60 )       & ( ariClustioutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(ariClustioutdf$AllAges_slope) > biosigslopeAllAges ) & ( ariClustioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    ariClustioutdf$BHadj_signifp <- ( ( ( ariClustioutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                        ( ariClustioutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                        ( ariClustioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    ariClustioutdf$LE60_log2FC <- ariClustioutdf$LE60_slope * 35
                    ariClustioutdf$GT60_log2FC <- ariClustioutdf$GT60_slope * 35
                    ariClustioutdf$AllAges_log2FC <- ariClustioutdf$AllAges_slope * 70
                    
                    ariClustioutdf$Best_log2FC <- ariClustioutdf$AllAges_log2FC
                    ariClustioutdf$Best_BHadj_pval <- ariClustioutdf$AllAges_BHadj_pval
                    
                    ## Find the largest significant fold change:

                    LE60_Bestp <- ( ( abs(ariClustioutdf$LE60_slope) > biosigslopeLE60 )          &
                                    ( ariClustioutdf$LE60_BHadj_pval < FDRalpha ) )
                    ariClustioutdf[LE60_Bestp, ]$Best_log2FC <- ariClustioutdf[LE60_Bestp, ]$LE60_log2FC
                    ariClustioutdf[LE60_Bestp, ]$Best_BHadj_pval <- ariClustioutdf[LE60_Bestp, ]$LE60_BHadj_pval
                    
                    GT60_Bestp <- ( ( abs(ariClustioutdf$GT60_slope) > biosigslopeGT60 )          &
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
                }

                ## No overlap 1161 cases:  No Overlap with Big Series and MB09 TMA cases
                ## iClustinov = Pam50 iClusti phenotype, from No Overlap with Big Series and MB09 TMA cases
                if ( icinm == "All cases" ) {
                    sehiClustinovdf <- sehnovdf
                } else {
                    sehiClustinovdf <- sehnovdf[sehnovdf$iClustf == icinm, ]
                }
                ariClustinovoutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]
                ariClustinovoutmatcolnames <-
                  c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
                ariClustinovoutmat <- matrix(NA_real_, nrow = nrow(ariClustinovoutdf), ncol = length(ariClustinovoutmatcolnames))
                LE60_slope_col <- match("LE60_slope", ariClustinovoutmatcolnames)
                GT60_slope_col <- match("GT60_slope", ariClustinovoutmatcolnames)
                AllAges_slope_col <- match("AllAges_slope", ariClustinovoutmatcolnames)
                LE60_pval_col <- match("LE60_pval", ariClustinovoutmatcolnames)
                GT60_pval_col <- match("GT60_pval", ariClustinovoutmatcolnames)
                AllAges_pval_col <- match("AllAges_pval", ariClustinovoutmatcolnames)
                AgeDependentp_col <- match("AgeDependentp", ariClustinovoutmatcolnames)
                ariClustinovoutdf$LE60_slope <- rep(NA_real_ , nrow(ariClustinovoutdf))
                ariClustinovoutdf$GT60_slope <- rep(NA_real_ , nrow(ariClustinovoutdf))
                ariClustinovoutdf$AllAges_slope <- rep(NA_real_ , nrow(ariClustinovoutdf))
                ariClustinovoutdf$LE60_pval <- rep(NA_real_ , nrow(ariClustinovoutdf))
                ariClustinovoutdf$GT60_pval <- rep(NA_real_ , nrow(ariClustinovoutdf))
                ariClustinovoutdf$AllAges_pval <- rep(NA_real_ , nrow(ariClustinovoutdf))
                ariClustinovoutdf$AgeDependentp <- rep(NA , nrow(ariClustinovoutdf))
                biosigslopeLE60 <-  log2(1.25)/35
                biosigslopeGT60 <-  log2(1.25)/35
                biosigslopeAllAges <-  log2(1.25)/70
                iClustinovidxs <- match(sehiClustinovdf$MBid, dimnames(Dataset.r)[[1]])
                ariClustinovlmdf <- data.frame(age_at_diagnosis = sehiClustinovdf$age_at_diagnosis,
                                               probesetni = Dataset.r[iClustinovidxs, 1])

                cat("\n\n\n")
                cat(icinm)
                cat("\n\n\n")
                ariClustinovoutdf_memoizedfilename <-
                  paste("AgeRelated_AllProbes_",
                        gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                        "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                        "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "")
                ariClustinovoutdf_memoizedfilenamep <- FALSE

                if ( (!DoMemoizep) && file.exists(ariClustinovoutdf_memoizedfilename) ) {
                    ariClustinovoutdf <- read.table(file = ariClustinovoutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                    sep = ",", header = TRUE, comment.char = "")
                    ariClustinovoutdf_memoizedfilenamep <- TRUE
                } else {
                    for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
                        psi <- ariClustinovoutdf$ProbeId[ni]
                        drci <- match(psi, dimnames(Dataset.r)[[2]])
                        ariClustinovlmdf$probesetni <- Dataset.r[iClustinovidxs, drci]
                        if (ni %% 100 == 0 ) { cat(ni, ", ") }
                        lmifitageLE60 <-
                          lm(probesetni ~ age_at_diagnosis, data = ariClustinovlmdf, na.action = na.exclude,
                             subset = age_at_diagnosis <= 60)
                        lmifitageGT60 <-
                          lm(probesetni ~ age_at_diagnosis, data = ariClustinovlmdf, na.action = na.exclude,
                             subset = age_at_diagnosis > 60)
                        lmifitallages <-
                          lm(probesetni ~ age_at_diagnosis, data = ariClustinovlmdf, na.action = na.exclude)
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        smrylmifitallages <- summary(lmifitallages)
                        
                        ariClustinovoutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                        ariClustinovoutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                        ariClustinovoutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                        ariClustinovoutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        ariClustinovoutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        ariClustinovoutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                    }

                    dimnames(ariClustinovoutmat)[[2]] <- ariClustinovoutmatcolnames
                    dimnames(ariClustinovoutmat)[[1]] <- ariClustinovoutdf$ProbeId
                    ariClustinovoutdf[, ariClustinovoutmatcolnames] <- ariClustinovoutmat
                    cat("\n\n\n")
                    ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                    ariClustinovoutdf$LE60_BHadj_pval <- p.adjust(ariClustinovoutdf$LE60_pval, method="BH")
                    ariClustinovoutdf$GT60_BHadj_pval <- p.adjust(ariClustinovoutdf$GT60_pval, method="BH")
                    ariClustinovoutdf$AllAges_BHadj_pval <- p.adjust(ariClustinovoutdf$AllAges_pval, method="BH")
                    ariClustinovoutdf$BHadj_and_AgeDependentp <-
                      ( ( ( abs(ariClustinovoutdf$LE60_slope) > biosigslopeLE60 )       &
                          ( ariClustinovoutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(ariClustinovoutdf$GT60_slope) > biosigslopeGT60 )       &
                          ( ariClustinovoutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(ariClustinovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                          ( ariClustinovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    ariClustinovoutdf$BHadj_signifp <- ( ( ( ariClustinovoutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                           ( ariClustinovoutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                           ( ariClustinovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    ariClustinovoutdf$LE60_log2FC <- ariClustinovoutdf$LE60_slope * 35
                    ariClustinovoutdf$GT60_log2FC <- ariClustinovoutdf$GT60_slope * 35
                    ariClustinovoutdf$AllAges_log2FC <- ariClustinovoutdf$AllAges_slope * 70

                    ariClustinovoutdf$Best_log2FC <- ariClustinovoutdf$AllAges_log2FC
                    ariClustinovoutdf$Best_BHadj_pval <- ariClustinovoutdf$AllAges_BHadj_pval

                    ## Find the largest significant fold change:

                    LE60_Bestp <- ( ( abs(ariClustinovoutdf$LE60_slope) > biosigslopeLE60 ) &
                                    ( ariClustinovoutdf$LE60_BHadj_pval < FDRalpha ) )
                    ariClustinovoutdf[LE60_Bestp, ]$Best_log2FC <- ariClustinovoutdf[LE60_Bestp, ]$LE60_log2FC
                    ariClustinovoutdf[LE60_Bestp, ]$Best_BHadj_pval <- ariClustinovoutdf[LE60_Bestp, ]$LE60_BHadj_pval

                    GT60_Bestp <- ( ( abs(ariClustinovoutdf$GT60_slope) > biosigslopeGT60 )          &
                                    ( ariClustinovoutdf$GT60_BHadj_pval < FDRalpha )                     &
                                    ( abs(ariClustinovoutdf$GT60_slope) > abs(ariClustinovoutdf$LE60_slope) ) )
                    ariClustinovoutdf[GT60_Bestp, ]$Best_log2FC <- ariClustinovoutdf[GT60_Bestp, ]$GT60_log2FC
                    ariClustinovoutdf[GT60_Bestp, ]$Best_BHadj_pval <- ariClustinovoutdf[GT60_Bestp, ]$GT60_BHadj_pval

                    AllAges_Bestp <- ( ( abs(ariClustinovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                                       ( ariClustinovoutdf$AllAges_BHadj_pval < FDRalpha )               &
                                       ( ( abs(ariClustinovoutdf$AllAges_log2FC) > abs(ariClustinovoutdf$LE60_log2FC) ) |
                                         ( abs(ariClustinovoutdf$AllAges_log2FC) > abs(ariClustinovoutdf$GT60_log2FC) ) ) )
                    ariClustinovoutdf[AllAges_Bestp, ]$Best_log2FC <- ariClustinovoutdf[AllAges_Bestp, ]$AllAges_log2FC
                    ariClustinovoutdf[AllAges_Bestp, ]$Best_BHadj_pval <- ariClustinovoutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
                    ariClustinovoutdf$Abs_Best_log2FC <- abs( ariClustinovoutdf$Best_log2FC )
                    ariClustinovoutdf$Abs_FoldChange <- 2^ariClustinovoutdf$Abs_Best_log2FC
                    ariClustinovoutdf$FoldChange_Direction <- ifelse(ariClustinovoutdf$Best_log2FC > 0, "Up", "Down")
                }

### Use this Best_BHadj_pval for manhattan plots as well.

### Plot scatterplots of gene sets:
### ar = age related study general gene set of interest for Tomo
### p  = polycomb independent
### prc2

### ar set:  iClusti from all 1992 cases
### iClustiarrdf <- sarrdf[sarrdf$iClustf == icinm, ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                ariClustioutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(ariClustioutdf))
                ariClustioutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(ariClustioutdf))
                ariClustioutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(ariClustioutdf))
                ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(ariClustioutdf))
                ariClustioutdf$arGeneSet_PBHadj <- rep(NA, nrow(ariClustioutdf))
                
                ariClustioutdf$arGeneSetidxp <- ( ariClustioutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                
                ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                  p.adjust(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                  p.adjust(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                  p.adjust(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                  apply(ariClustioutdf[ariClustioutdf$arGeneSetidxp,
                                       c("arGeneSet_LE60_BHadj_pval", "arGeneSet_GT60_BHadj_pval", "arGeneSet_AllAges_BHadj_pval")], 1, min)


                ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                  ( ( ( abs(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
                      ( ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < FDRalpha ) ) |
                    ( ( abs(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
                      ( ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < FDRalpha ) ) |
                    ( ( abs(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                      ( ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                


### iClusti All 1992 cases
### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in ariClustioutdf
                iClusticlfitageLE60 <- vector("list", length(names(arzidxs)))
                iClusticlfitageGT60 <- vector("list", length(names(arzidxs)))
                iClusticlfit <- vector("list", length(names(arzidxs)))
                stj <- paste("N =", dim(iClustiarrdf)[1])

                ## Scatterplots for EZH2 H3K27me3 pathway genes arGeneSet
                if ( SingleOutputFilep && arScatterPlotsp ) {
                    pdf(file = paste("./Plots/MBEX_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                        width = 8, height = 10, useDingbats = FALSE)
                    par(mfrow = c(2, 2))
                }

                xlims <- c(20, 100)
                ylims <- c(4, 17)

                ageLE60idxp <- (iClustiarrdf$age_at_diagnosis <= 60)
                ageGT60idxp <- (!ageLE60idxp)
                for ( ni in seq( along = names(arzidxs) ) ) {

                    i <- order(toupper(names(iClustiarrdf)[seq(length(names(arzidxs)))]))[ni]
                    cat("\n\n### --- ", names(iClustiarrdf)[i])
                    if ( !SingleOutputFilep && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/ProbeLevel/MBEX_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_", gsub("/", "_", names(iClustiarrdf)[i]),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                            width = 6, height = 7, useDingbats = FALSE)
                        par(mfrow = c(1, 1))
                    }
                    lmifitageLE60 <- lm(iClustiarrdf[ageLE60idxp, i] ~ iClustiarrdf[ageLE60idxp, "age_at_diagnosis"])
                    smrylmifitageLE60 <- summary(lmifitageLE60)
                    iClusticlfitageLE60[[i]] <- lmifitageLE60
                    if ( Verbosep ) print(smrylmifitageLE60)
                    lmifitageGT60 <- lm(iClustiarrdf[ageGT60idxp, i] ~ iClustiarrdf[ageGT60idxp, "age_at_diagnosis"])
                    smrylmifitageGT60 <- summary(lmifitageGT60)
                    iClusticlfitageGT60[[i]] <- lmifitageGT60
                    if ( Verbosep ) print(smrylmifitageGT60)
                    lmifit <- lm(iClustiarrdf[, i] ~ iClustiarrdf[, "age_at_diagnosis"])
                    smrylmifit <- summary(lmifit)
                    iClusticlfit[[i]] <- lmifit
                    if ( Verbosep ) print(smrylmifit)
                    if ( arScatterPlotsp ) {
                        plot(iClustiarrdf[, "age_at_diagnosis"], iClustiarrdf[, i],
                             ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                             main = "", xlim = xlims, ylim = ylims, type = "n")
                        title( main = paste(icinm, stj, sep = ": "), line = 2)
                        title( main =
                                 paste(names(iClustiarrdf)[i], "(",
                                       ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(iClustiarrdf)[i],
                                                                                   split = "\\|")[[1]][2], "ERbinding"],
                                              "ER binding", "non-ER binding"), ")"), line = 1)
                        points(iClustiarrdf[, "age_at_diagnosis"], iClustiarrdf[, i],
                               pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                        lines(supsmu(iClustiarrdf[, "age_at_diagnosis"], iClustiarrdf[, i], span = 0.4, bass = 10),
                              lwd = 3, col = "#AA1010")
                        abline(h = 0, v = 50, lty = 2)
                        text(x = 20, y = 16.3,
                             labels = paste("[<=60] p = ",
                                            format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 20, y = 15.5,
                             labels = paste("[>60] p = ",
                                            format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 98, y = 16.3,
                             labels = paste("[All] p = ",
                                            format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 1, cex = 0.85)
                        ## Need iClusti subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                        if ( strsplit(names(iClustiarrdf)[i], split = "\\|")[[1]][2] %in%
                             ariClustioutdf[ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                            text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|>1.25 and adjPval<",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        } else {
                            text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels=paste("|FC|<1.25 or adjPval>",
                                              format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        }
                    }
                    if (!SingleOutputFilep && arScatterPlotsp ) { dev.off() }
                }
                if (SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(ariClustioutdf) %in% skipVarNms )
                if ( DoMemoizep || (! ariClustioutdf_memoizedfilenamep ) ) {
                    write.csv(ariClustioutdf[, -c(skipVarCols)], file = ariClustioutdf_memoizedfilename )
                    write.csv(ariClustioutdf[ariClustioutdf$arGeneSetidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_arGeneSet_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_All1992Cases_lm_v01.csv", sep = "") )
                }


### ar set:  iClusti from 1161 no overlap cases (shown in asudf$MBid)
### iClustiarrnovdf <- sarrdf[rownames(sarrdf) %in% asudf$MBid & sarrdf$Pam50Subtype == "iClusti", ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                ariClustinovoutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(ariClustinovoutdf))
                ariClustinovoutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(ariClustinovoutdf))
                ariClustinovoutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(ariClustinovoutdf))
                ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(ariClustinovoutdf))
                ariClustinovoutdf$arGeneSet_PBHadj <- rep(NA, nrow(ariClustinovoutdf))
                
                ariClustinovoutdf$arGeneSetidxp <- ( ariClustinovoutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                
                ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                  p.adjust(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                  p.adjust(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                  p.adjust(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                  apply(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp,
                                          c("arGeneSet_LE60_BHadj_pval",
                                            "arGeneSet_GT60_BHadj_pval",
                                            "arGeneSet_AllAges_BHadj_pval")], 1, min)

                
                ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                  ( ( ( abs(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$LE60_slope)    > biosigslopeLE60 )    &
                      ( ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval    < FDRalpha ) ) |
                    ( ( abs(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$GT60_slope)    > biosigslopeGT60 )    &
                      ( ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval    < FDRalpha ) ) |
                    ( ( abs(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                      ( ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                

### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in ariClustinovoutdf
                iClusticlfitageLE60 <- vector("list", length(names(arzidxs)))
                iClusticlfitageGT60 <- vector("list", length(names(arzidxs)))
                iClusticlfit <- vector("list", length(names(arzidxs)))
                stj <- paste("N =", dim(iClustiarrnovdf)[1], "(No overlap)")
                
                if ( SingleOutputFilep && arScatterPlotsp ) {
                    pdf(file = paste("./Plots/MBEX_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                        width = 8, height = 10, useDingbats = FALSE)
                    par(mfrow = c(2, 2))
                }
                xlims <- c(20, 100)
                ylims <- c(4, 17)

                ageLE60idxp <- (iClustiarrnovdf$age_at_diagnosis <= 60)
                ageGT60idxp <- (!ageLE60idxp)
                for ( ni in seq( along = names(arzidxs) ) ) {
                    i <- order(toupper(names(iClustiarrnovdf)[seq(length(names(arzidxs)))]))[ni]
                    cat("\n\n### --- ", names(iClustiarrnovdf)[i])
                    
                    if ( (!SingleOutputFilep) && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/ProbeLevel/MBEX_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), 
                                         "_", gsub("/", "_", names(iClustiarrnovdf)[i]), 
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                            width = 6, height = 7, useDingbats = FALSE)
                        par(mfrow = c(1, 1))
                    }

                    lmifitageLE60 <- lm(iClustiarrnovdf[ageLE60idxp, i] ~ iClustiarrnovdf[ageLE60idxp, "age_at_diagnosis"])
                    smrylmifitageLE60 <- summary(lmifitageLE60)
                    iClusticlfitageLE60[[i]] <- lmifitageLE60
                    if ( Verbosep ) print(smrylmifitageLE60)
                    lmifitageGT60 <- lm(iClustiarrnovdf[ageGT60idxp, i] ~ iClustiarrnovdf[ageGT60idxp, "age_at_diagnosis"])
                    smrylmifitageGT60 <- summary(lmifitageGT60)
                    iClusticlfitageGT60[[i]] <- lmifitageGT60
                    if ( Verbosep ) print(smrylmifitageGT60)
                    lmifit <- lm(iClustiarrnovdf[, i] ~ iClustiarrnovdf[, "age_at_diagnosis"])
                    smrylmifit <- summary(lmifit)
                    iClusticlfit[[i]] <- lmifit
                    if ( Verbosep ) print(smrylmifit)
                    if ( arScatterPlotsp ) {
                        plot(iClustiarrnovdf[, "age_at_diagnosis"], iClustiarrnovdf[, i],
                             ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                             main = "", xlim = xlims, ylim = ylims, type = "n")
                        title( main = paste(icinm, stj, sep = ": "), line = 2)
                        title( main =
                                 paste(names(iClustiarrnovdf)[i], "(",
                                       ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(iClustiarrnovdf)[i],
                                                                                   split = "\\|")[[1]][2], "ERbinding"],
                                              "ER binding", "non-ER binding"), ")"), line = 1)
                        points(iClustiarrnovdf[, "age_at_diagnosis"], iClustiarrnovdf[, i],
                               pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                        lines(supsmu(iClustiarrnovdf[, "age_at_diagnosis"], iClustiarrnovdf[, i], span = 0.4, bass = 10),
                              lwd = 3, col = "#AA1010")
                        abline(h = 0, v = 50, lty = 2)
                        text(x = 20, y = 16.3,
                             labels = paste("[<=60] p = ",
                                            format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 20, y = 15.5,
                             labels = paste("[>60] p = ",
                                            format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 98, y = 16.3,
                             labels = paste("[All] p = ",
                                            format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 1, cex = 0.85)
                        ## Need iClusti subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                        if ( strsplit(names(iClustiarrnovdf)[i], split = "\\|")[[1]][2] %in%
                             ariClustinovoutdf[ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                            text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|>1.25 and adjPval<",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        } else {
                            text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|<1.25 or adjPval>",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        }
                    }
                    if ( (!SingleOutputFilep) && arScatterPlotsp ) { dev.off() }
                }
                if ( SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(ariClustinovoutdf) %in% skipVarNms )
                if ( DoMemoizep || (! ariClustinovoutdf_memoizedfilenamep ) ) {
                    write.csv(ariClustinovoutdf[, -c(skipVarCols)],
                              file = ariClustinovoutdf_memoizedfilename )
                    write.csv(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_arGeneSet_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
                }



### Plot age-associated p-values by chromosome position - Manhattan plot

### qqman package has manhattan plot

                require("qqman")

### Data frame with CHR, BP, P (and SNP to avoid warnings)
### Genomic_location
### chr2:206352192:206352241:+
                ariClustioutdf$Chrcr <- sapply(strsplit(ariClustioutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                ariClustioutdf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "",
                                                                  gsub("_h2_hap1", "", gsub("_random", "", ariClustioutdf$Chrcr))))
                ariClustioutdf$Chrn <- as.numeric(gsub("chr", "", ariClustioutdf$Chrc))
                ariClustioutdf$Chrn[ariClustioutdf$Chrc == "chrX"] <- 23
                ariClustioutdf$Chrn[ariClustioutdf$Chrc == "chrY"] <- 24
                with(ariClustioutdf, table(Chrcr, Chrn, useNA = "always"))
                with(ariClustioutdf, table(Chrc, Chrn, useNA = "always"))
                ariClustioutdf$CHR <- ariClustioutdf$Chrn
                ariClustioutdf$Startcr <- sapply(strsplit(ariClustioutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                ariClustioutdf$Stopcr <- sapply(strsplit(ariClustioutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                ariClustioutdf$Startn <- as.numeric(ariClustioutdf$Startcr)
                ariClustioutdf$Stopn <- as.numeric(ariClustioutdf$Stopcr)
                ariClustioutdf$BP <- trunc((ariClustioutdf$Startn + ariClustioutdf$Stopn)/2)
                ariClustioutdf$SNP <- ariClustioutdf$Probe_id
                ariClustioutdf$P <- apply(ariClustioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                ariClustioutdf$PBHadj <- apply(ariClustioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                ariClustioutdfarGeneSetManhattanidxp <- ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp & ariClustioutdf$Chrn <= 23
                ariClustioutdfarGeneSetManhattanidxp[is.na(ariClustioutdfarGeneSetManhattanidxp)] <- FALSE
                ariClustioutdfManhattanidxp <- ariClustioutdf$BHadj_and_AgeDependentp & ariClustioutdf$Chrn <= 23
                ariClustioutdfManhattanidxp[is.na(ariClustioutdfManhattanidxp)] <- FALSE
                ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                ariClustioutdfManhattanVarsidxp <-
                  ( apply(ariClustioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                          function(x) !any(is.na(x)))  & ( ariClustioutdf$Chrn <= 23 ) )

                ariClustioutdf[ariClustioutdf$Gene_symbol == "", "Gene_symbol"] <-
                  ariClustioutdf[ariClustioutdf$Gene_symbol == "", "ILMN_Gene_0"]
                ## Number of probesets showing evidence of association with age
                NpsiClustiAssocAge <- sum(ariClustioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets"] <- NpsiClustiAssocAge
                ## Number of probesets showing evidence of association with age with genomic location
                NpsiClustiAssocAgecGL <- sum(ariClustioutdfManhattanidxp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc"] <- NpsiClustiAssocAgecGL
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgiClustiAssocAgecGL <- length(unique(ariClustioutdf[ariClustioutdfManhattanidxp, "ILMN_Gene_0"]))
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames"] <- NgiClustiAssocAgecGL
                ## arGeneSet subset
                ## Number of probesets showing evidence of association with age
                NpsiClustiAssocAge_arGeneSet <- sum(ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_arGeneSet"] <- NpsiClustiAssocAge_arGeneSet
                ## Number of probesets showing evidence of association with age with genomic location
                NpsiClustiAssocAgecGL_arGeneSet <- sum(ariClustioutdfarGeneSetManhattanidxp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <- NpsiClustiAssocAgecGL_arGeneSet
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgiClustiAssocAgecGL_arGeneSet <- length(unique(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames_arGeneSet"] <- NgiClustiAssocAgecGL_arGeneSet

                ## Non-overlapping with Big Series
                ariClustinovoutdf$Chrcr <- sapply(strsplit(ariClustinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                ariClustinovoutdf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "",
                                                                     gsub("_h2_hap1", "", gsub("_random", "", ariClustinovoutdf$Chrcr))))
                ariClustinovoutdf$Chrn <- as.numeric(gsub("chr", "", ariClustinovoutdf$Chrc))
                ariClustinovoutdf$Chrn[ariClustinovoutdf$Chrc == "chrX"] <- 23
                ariClustinovoutdf$Chrn[ariClustinovoutdf$Chrc == "chrY"] <- 24
                with(ariClustinovoutdf, table(Chrcr, Chrn, useNA = "always"))
                with(ariClustinovoutdf, table(Chrc, Chrn, useNA = "always"))
                ariClustinovoutdf$CHR <- ariClustinovoutdf$Chrn
                ariClustinovoutdf$Startcr <- sapply(strsplit(ariClustinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                ariClustinovoutdf$Stopcr <- sapply(strsplit(ariClustinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                ariClustinovoutdf$Startn <- as.numeric(ariClustinovoutdf$Startcr)
                ariClustinovoutdf$Stopn <- as.numeric(ariClustinovoutdf$Stopcr)
                ariClustinovoutdf$BP <- trunc((ariClustinovoutdf$Startn + ariClustinovoutdf$Stopn)/2)
                ariClustinovoutdf$SNP <- ariClustinovoutdf$Probe_id
                ariClustinovoutdf$P <- apply(ariClustinovoutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                ariClustinovoutdf$PBHadj <- apply(ariClustinovoutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                ariClustinovoutdfarGeneSetManhattanidxp <- ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp & ariClustinovoutdf$Chrn <= 23
                ariClustinovoutdfarGeneSetManhattanidxp[is.na(ariClustinovoutdfarGeneSetManhattanidxp)] <- FALSE
                ariClustinovoutdfManhattanidxp <- ariClustinovoutdf$BHadj_and_AgeDependentp & ariClustinovoutdf$Chrn <= 23
                ariClustinovoutdfManhattanidxp[is.na(ariClustinovoutdfManhattanidxp)] <- FALSE
### Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                ariClustinovoutdfManhattanVarsidxp <-
                  ( apply(ariClustinovoutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                          function(x) !any(is.na(x)))  & ( ariClustinovoutdf$Chrn <= 23 ) )

                ariClustinovoutdf[ariClustinovoutdf$Gene_symbol == "", "Gene_symbol"] <-
                  ariClustinovoutdf[ariClustinovoutdf$Gene_symbol == "", "ILMN_Gene_0"]
                ## Number of probesets showing evidence of association with age
                NpsiClustinovAssocAge <- sum(ariClustinovoutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets"] <- NpsiClustinovAssocAge
                ## Number of probesets showing evidence of association with age with genomic location
                NpsiClustinovAssocAgecGL <- sum(ariClustinovoutdfManhattanidxp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc"] <- NpsiClustinovAssocAgecGL
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgiClustinovAssocAgecGL <- length(unique(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, "ILMN_Gene_0"]))
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames"] <- NgiClustinovAssocAgecGL
                ## arGeneSet subset
                ## Number of probesets showing evidence of association with age
                NpsiClustinovAssocAge_arGeneSet <- sum(ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_arGeneSet"] <- NpsiClustinovAssocAge_arGeneSet
                ## Number of probesets showing evidence of association with age with genomic location
                NpsiClustinovAssocAgecGL_arGeneSet <- sum(ariClustinovoutdfarGeneSetManhattanidxp, na.rm = TRUE)
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <-
                  NpsiClustinovAssocAgecGL_arGeneSet
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgiClustinovAssocAgecGL_arGeneSet <-
                  length(unique(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames_arGeneSet"] <-
                  NgiClustinovAssocAgecGL_arGeneSet
                
                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(ariClustioutdf) %in% skipVarNms )

                if ( sum(ariClustioutdfManhattanidxp) ) {
                    write.csv(ariClustioutdf[ariClustioutdfManhattanidxp, -c(skipVarCols)][
                        order(ariClustioutdf[ariClustioutdfManhattanidxp, "Abs_FoldChange"], decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.csv(ariClustioutdf[ariClustioutdfManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.table(ariClustioutdf[ariClustioutdfManhattanidxp, ][
                        !duplicated(ariClustioutdf[ariClustioutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(ariClustioutdf[ariClustioutdfManhattanidxp, ][
                            !duplicated(ariClustioutdf[ariClustioutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(ariClustinovoutdf) %in% skipVarNms )

                if ( sum(ariClustinovoutdfManhattanidxp) ) {
                    write.csv(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, -c(skipVarCols)][
                        order(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.csv(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.table(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, ][
                        !duplicated(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, ][
                            !duplicated(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_outdf_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                ## arGeneSet subset
                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(ariClustioutdf) %in% skipVarNms )
                if ( sum ( ariClustioutdfarGeneSetManhattanidxp ) ) {
                    write.csv(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                        order(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.csv(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_arGeneSet_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.table(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, ][
                        !duplicated(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, ][
                            !duplicated(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(ariClustinovoutdf) %in% skipVarNms )
                if ( sum ( ariClustinovoutdfarGeneSetManhattanidxp ) ) {
                    write.csv(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                        order(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv",
                                     sep = "") )
                    write.csv(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_arGeneSet_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.table(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, ][
                      !duplicated(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                        !duplicated(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, ][
                          !duplicated(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }
                
                ## Do stacked Manhattan and volcano plot for All1992Cases.  Plot log(FC) on X axis, log10(pval) on Y axis

                if ( !SingleOutputFilep ) {
                    pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 8, height = 8, useDingbats = FALSE)
                    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                } else {
                    pdf(file = paste("./Plots/AgeRelated_Manhattan_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                }

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
                
                if ( length( ObjsToHighlight ) ) {
                    sm_manhattan(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj", chrlabs = c(1:22, "X"),
                                 suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                 highlight = ObjsToHighlight,
                                 plotpointidxp = ariClustioutdf[ariClustioutdfManhattanVarsidxp,
                                                                c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp =
                                   if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                } else {
                    sm_manhattan(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj", chrlabs = c(1:22, "X"),
                                 suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                 plotpointidxp = ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],
                                 plotpointhiliteidxp =
                                   if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                }
                title(paste("METABRIC:", icinm), line = 2)
                title("Expression showing age-related trend, FC > 1.25", line = 1)

                write.csv(ariClustioutdf[which(ariClustioutdf$SNP %in% ObjsToHighlight), ],
                          file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_All1992Cases.csv", sep = "") )
                
                if ( SingleOutputFilep ) {
                    dev.off()
                    pdf(file = paste("./Plots/AgeRelated_Volcano_",
                                 if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   

                }
                ## Volcano
                ## Label top values
                ## Determine objects to highlight from volcano plot, label in both manhattan and volcano

                stairstepvals <-
                  switch(as.character(FDRalpha),
                    "0.05" = switch(as.character(icinm),
                               "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                               "iClust 1" = c(log2(4), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                               "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                               "iClust 3" = c(log2(4), 1.8, log2(3), 3, log2(1.25), 5), 
                               "iClust 4" = c(log2(3), -log10(FDRalpha), log2(2.5), 2.5, log2(1.25), 3.5), 
                               "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                               "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                               "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                               "iClust 8" = c(log2(3), -log10(FDRalpha), log2(2.5), 4, log2(1.25), 5), 
                               "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                               "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                               ),
                     "0.01" = switch(as.character(icinm),
                                "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                                "iClust 1" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                "iClust 3" = c(log2(4), 2, log2(3), 3, log2(1.25), 5), 
                                "iClust 4" = c(log2(4), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/5), log2(1.25), 3.5), 
                                "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                "iClust 8" = c(log2(4), -log10(FDRalpha), log2(3), 3, log2(1.25), 5), 
                                "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                                "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                )
                         )
                
                volcanoaddlidxp <- ( ( ( abs(ariClustioutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                       ( -log10(ariClustioutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                     ( ( abs(ariClustioutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                       ( -log10(ariClustioutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                     ( ( abs(ariClustioutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                       ( -log10(ariClustioutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                ObjsToHighlight <- ariClustioutdf[volcanoaddlidxp, "SNP"]
                
                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                mainTitle <- paste("METABRIC BrCa", icinm)

                ERbindingAgeDep <- which( ( ariClustioutdf$ERbinding &
                                            (abs(ariClustioutdf$Best_log2FC ) > log2(FCthresh)) &
                                            (ariClustioutdf$Best_BHadj_pval < alphacrit) ) )
                ERnonbindingAgeDep <- which( ( !ariClustioutdf$ERbinding &
                                               (abs(ariClustioutdf$Best_log2FC ) > log2(FCthresh)) &
                                               (ariClustioutdf$Best_BHadj_pval < alphacrit) ) )
                NotAgeDepidxs <- setdiff(seq(nrow(ariClustioutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                mainTitleN <- paste("Cases:", nrow(sehiClustidf),
                                    "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                
                if ( ! EqAxesScalesp ) {
###                    if  (niClusti == 1) {
                    if  (ici == 1) {
###                         xlims_vpEVS <- max(abs(range(ariClustioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylims_vpEVS <- c( 0, 1.2*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
                    } else {
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylimsi <- c( 0, 1.2*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
                        ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                        ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                    }
                }
                ## Will need to cull data points at root of volcano - 48000 points too much.
                ## Get density information to thin out volcano plot
                ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                NADn <- length(NotAgeDepidxs)
                ClustRepsnRem <- NADn %% maxClusterN
                ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered
                

                for ( csi in seq(ClustRepsn) ) {
                    cat("   Cluster rep: ", csi, "  ")
                    if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                    csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                
                    densityEst <-
                      hclust(dist(cbind(ariClustioutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                        jitter(-log10(ariClustioutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                     method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                    lowDensityiidxs <-
                      (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityiidxs <-
                      sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                    if ( csi == 1 ) {
                        lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                        hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                    } else {
                        lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                        hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                    }
                }
                    
                plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                
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
                                c( 0, 1.2*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
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
                legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                VolcanoAddlObjsToHighlight <- c()
                evcp <- ariClustioutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                ## Cull out duplicate gene names
                evcidx <- unlist(tapply(which(evcp), ariClustioutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                evcERbidx <- evcidx[which(ariClustioutdf[evcidx, "ERbinding"])]
                evcERnbidx <- setdiff(evcidx, evcERbidx)
                lblcol <- rep("red", length(evcidx))
                lblcol[evcidx %in% evcERnbidx] <- "blue"
                if ( length(evcidx)  ) {
                    require("maptools"); ## for pointLabel
                    pointLabel(ariClustioutdf[evcidx, ]$Best_log2FC, -log10(ariClustioutdf[evcidx, ]$Best_BHadj_pval),
                               labels = ariClustioutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                }
                write.csv(ariClustioutdf[evcidx, ],
                          file = paste("../Data/AgeRelated_Volcano_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       gsub(" ", "",
                                            paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                           2^stairstepvals[3], stairstepvals[4],
                                                           2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_All1992Cases.csv", sep = "") )

                dev.off()


                ## End All1992Cases

                ## No overlap with Big Series or MB09
                ## Do stacked Manhattan and volcano plot for NoOverlapMB09BigSeries.  Plot log(FC) on X axis, log10(pval) on Y axis

                if ( !SingleOutputFilep ) {
                    pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 8, height = 8, useDingbats = FALSE)
                    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                } else {
                    pdf(file = paste("./Plots/AgeRelated_Manhattan_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                }

                ## Label the top 20 largest fold change with age association and adjusted p-val significance
                ohdf <- ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, ][
                    order(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, "Abs_FoldChange"], decreasing = TRUE), ]
                ohdf <- ohdf[ohdf$BHadj_and_AgeDependentp, ]
                
                ObjsToHighlight <-
                  unique(c(ohdf[ohdf$FoldChange_Direction == "Up", ][1:10, "SNP"],
                           ohdf[ohdf$FoldChange_Direction == "Down", ][1:10, "SNP"],
                           ohdf[1:20, "SNP"]
                           ))

                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                if ( length( ObjsToHighlight ) ) {
                    sm_manhattan(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj",
                                 chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha),
                                 genomewideline = FALSE, highlight = ObjsToHighlight,
                                 plotpointidxp = ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp =
                                   if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                       ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                } else {
                    sm_manhattan(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj",
                                 chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha), genomewideline = FALSE, 
                                 plotpointidxp = ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp =
                                   if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                       ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                }
                title(paste("METABRIC:", icinm), line = 2)
                title("Expression showing age-related trend, FC > 1.25", line = 1)
                
                write.csv(ariClustinovoutdf[which(ariClustinovoutdf$SNP %in% ObjsToHighlight), ],
                          file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_NoOverlapMB09BigSeries.csv", sep = "") )
                
                if ( SingleOutputFilep ) {
                    dev.off()
                    pdf(file = paste("./Plots/AgeRelated_Volcano_",
                                 if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1)) 
                }
                ## Volcano
                ## Determine objects to highlight from volcano plot, label in both manhattan and volcano

                stairstepvals <-
                  switch(as.character(FDRalpha),
                         "0.05" = switch(as.character(icinm),
                                         "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                         "iClust 1" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 3" = c(log2(4), -log10(FDRalpha), log2(3), -log10(FDRalpha/5), log2(1.25), 2.5), 
                                         "iClust 4" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                         "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                         "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 8" = c(log2(3), -log10(FDRalpha), log2(2.5), 2, log2(1.25), 2.5), 
                                         "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                                         "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                         ),
                         "0.01" = switch(as.character(icinm),
                                         "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                         "iClust 1" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 3" = c(log2(4), -log10(FDRalpha), log2(2.5), 2, log2(1.25), 2.5),
                                         "iClust 4" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "iClust 8" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                         "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                                         "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                         )
                         )
                
                volcanoaddlidxp <- ( ( ( abs(ariClustinovoutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                       ( -log10(ariClustinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                     ( ( abs(ariClustinovoutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                       ( -log10(ariClustinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                     ( ( abs(ariClustinovoutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                       ( -log10(ariClustinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                ObjsToHighlight <- ariClustinovoutdf[volcanoaddlidxp, "SNP"]

                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                mainTitle <- paste("METABRIC BrCa", icinm)

                ERbindingAgeDep <- which( ( ariClustinovoutdf$ERbinding &
                                            (abs(ariClustinovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                            (ariClustinovoutdf$Best_BHadj_pval < alphacrit) ) )
                ERnonbindingAgeDep <- which( ( !ariClustinovoutdf$ERbinding &
                                               (abs(ariClustinovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                               (ariClustinovoutdf$Best_BHadj_pval < alphacrit) ) )
                NotAgeDepidxs <- setdiff(seq(nrow(ariClustinovoutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                mainTitleN <- paste("Cases:", nrow(sehiClustinovdf),
                                    "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                if ( ! EqAxesScalesp ) {
                    if  (ici == 1) {
###                     if  (niClusti == 1) {
###                         xlims_vpEVS <- max(abs(range(ariClustinovoutdf$Best_log2FC, na.rm = TRUE))) * c(-1, 1) * 1.1
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylims_vpEVS <- c( 0, 1.2 * max(-log10(ariClustinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                    } else {
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylimsi <- c( 0, 1.2 * max(-log10(ariClustinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                        ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                        ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                    }
                }
                ## Will need to cull data points at root of volcano - 48000 points too much.
                ## Get density information to thin out volcano plot
                ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                NADn <- length(NotAgeDepidxs)
                ClustRepsnRem <- NADn %% maxClusterN
                ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered
                

                for ( csi in seq(ClustRepsn) ) {
                    cat("   Cluster rep: ", csi, "  ")
                    if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                    csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                
                    densityEst <-
                      hclust(dist(cbind(ariClustinovoutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                        jitter(-log10(ariClustinovoutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                     method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                    lowDensityiidxs <-
                      (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityiidxs <-
                      sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                    if ( csi == 1 ) {
                        lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                        hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                    } else {
                        lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                        hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                    }
                }
                    
                plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                
                plot(ariClustinovoutdf$Best_log2FC, -log10(ariClustinovoutdf$Best_BHadj_pval), type = "n", main = "", 
                     xlim =  if ( EqAxesScalesp ) {
                                 xlims_vpEVS
                             } else {
                                 max(abs(range(ariClustinovoutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                             },
                     xlab = "(Down)        FC        (Up)   ",
                     ylab = "-log10(BH adjusted P-values)", xaxt = "n",
                     ylim = if ( EqAxesScalesp ) {
                                ylims_vpEVS
                            } else {
                                c( 0, 1.2 * max(-log10(ariClustinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                            }
                     )

                title(main = mainTitle, line = 2)
                title(main = mainTitleN, line = 1)
                axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), labels = c(32, 16, 8, 4, 2, 1, 2, 4, 8, 16, 32) )
                
                points(ariClustinovoutdf$Best_log2FC[NotAgeDepidxs][plotNotAgeDepidxs],
                       -log10(ariClustinovoutdf$Best_BHadj_pval[NotAgeDepidxs][plotNotAgeDepidxs]),
                       pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))

                points(ariClustinovoutdf$Best_log2FC[ERbindingAgeDep], -log10(ariClustinovoutdf$Best_BHadj_pval[ERbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("red", alpha = FALSE)), alpha=30, max=255))
                points(ariClustinovoutdf$Best_log2FC[ERnonbindingAgeDep], -log10(ariClustinovoutdf$Best_BHadj_pval[ERnonbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("blue", alpha = FALSE)), alpha=30, max=255))
                legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                VolcanoAddlObjsToHighlight <- c()
                evcp <- ariClustinovoutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                ## Cull out duplicate gene names
                evcidx <- unlist(tapply(which(evcp), ariClustinovoutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                evcERbidx <- evcidx[which(ariClustinovoutdf[evcidx, "ERbinding"])]
                evcERnbidx <- setdiff(evcidx, evcERbidx)
                lblcol <- rep("red", length(evcidx))
                lblcol[evcidx %in% evcERnbidx] <- "blue"
                if ( length(evcidx)  ) {
                    require("maptools"); ## for pointLabel
                    pointLabel(ariClustinovoutdf[evcidx, ]$Best_log2FC, -log10(ariClustinovoutdf[evcidx, ]$Best_BHadj_pval),
                               labels = ariClustinovoutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                }
                ## Save ojects labeled for labeling in TCGA plots
                write.csv(ariClustinovoutdf[evcidx, ],
                          file = paste("../Data/AgeRelated_Volcano_PointsToLabel_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       gsub(" ", "",
                                            paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                           2^stairstepvals[3], stairstepvals[4],
                                                           2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_NoOverlapMB09BigSeries.csv", sep = "") )
                dev.off()
                ## End No overlap with Big Series or MB09
            }   ## End iClust i loop
        }   ## End niClusti loop


        iClustiAgeAssociatedProbesetsGenesdf$N_WholeSeries <- as.vector(table(sehdf$iClustf))
        iClustiAgeAssociatedProbesetsGenesdf$N_NoOverlapBigSeries <- as.vector(table(sehnovdf$iClustf))
        
        ## Write out data on counts of probesets/genes showing age association
        write.csv(iClustiAgeAssociatedProbesetsGenesdf,
                  file = paste("AgeAssociatedProbesetGenes_iClust_All",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               ".csv", sep = "") )
        
### AgeRelated_arGeneSet_xxx_All1992Cases_lm_v01.csv
### AgeRelated_arGeneSet_xxx_NoOverlapMB09BigSeries_lm_v01.csv
### Read in all these files, add iClust col, save as one file.

        for ( ici in seq(levels(sehdf$iClustf)) ) {
            icinm <- levels(sehdf$iClustf)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$iClustf <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_iClust_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               "_All1992Cases_lm_v01.csv", sep = "") )

        
        for ( ici in seq(levels(sehnovdf$iClustf)) ) {
            icinm <- levels(sehnovdf$iClustf)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$iClustf <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_iClust_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
    }
}


################################################################################
### PAM50 groups

sehdf$Pam50f <- factor(sehdf$Pam50Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal") )
sehnovdf$Pam50f <- factor(sehnovdf$Pam50Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal") )
sehrdf$Pam50f <- factor(sehrdf$Pam50Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal") )
sarrdf$Pam50f <- factor(sarrdf$Pam50Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal") )
sarrnovdf$Pam50f <- factor(sarrnovdf$Pam50Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal") )

################################################################################
### PAM50 groups

### Pam50 loop:

### arGeneSet tag denotes 467 Illumina probesets corresponding to the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.
### Otherwise analysis is for all probesets.
### Set up output for age-associated transcript expression trends
### - across whole genome
### - across arGeneSet of 244 genes (467 probesets) for the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.

### Pam50 i across whole genome Loop Begin   ici <- 0
##
## TODO 20170208:  Get asterisk labeling in manhattan plots to identify ERbinding targets - for intClust, BrCaEH, All
##                 as in sm_manhattan calls below for Pam50 - Rerun Pam50 to get proper manhattan plots
### Pam50
FDRalphalevels <- c(0.05, 0.01)
Verbosep <- FALSE ## TRUE
DoMemoizep <-  FALSE ## TRUE
SingleOutputFilep <- TRUE
arScatterPlotsp <- TRUE
EqAxesScales <- c(FALSE, TRUE)  ## Must be in order F, T so appropriate scale range can be calculated across conditions
nNotAgeDepToPlot <- 1500
propLowDens <- 0.66
propHiDens <- (1.0 - propLowDens)
maxClusterN <- 15000

for ( FDRalphai in seq(along = FDRalphalevels) ) {

    for ( EqAxesScalesi in seq(along = EqAxesScales) ) {

        FDRalpha <- FDRalphalevels[FDRalphai]
        EqAxesScalesp <- EqAxesScales[EqAxesScalesi]

        for ( nPam50i in c( 1, length(levels(sehdf$Pam50f)) ) ) {
            if (nPam50i == 1) {
                Pam50iAgeAssociatedProbesetsGenesdf <-
                  data.frame(Pam50i = "All cases",
                             N_WholeSeries = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_Probesets = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_GeneNames = rep(NA_integer_, nPam50i),
                             N_NoOverlapBigSeries = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_Probesets = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_GeneNames = rep(NA_integer_, nPam50i),

                             N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nPam50i),

                             N_NOv_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nPam50i)
                             )

            } else {
                Pam50iAgeAssociatedProbesetsGenesdf <-
                  data.frame(Pam50i = levels(sehdf$Pam50f),
                             N_WholeSeries = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_Probesets = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_GeneNames = rep(NA_integer_, nPam50i),
                             N_NoOverlapBigSeries = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_Probesets = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_GeneNames = rep(NA_integer_, nPam50i),

                             N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nPam50i),
                             N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nPam50i),

                             N_NOv_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, nPam50i),
                             N_NOv_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, nPam50i)
                             )
            }
            
            for ( ici in seq( nPam50i ) ) { ## ici <- 1

                if ( nPam50i == 1 ) {
                    icinm <- "All cases"
                } else {
                    icinm <- levels(sehdf$Pam50f)[ici]
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
### Pam50inov = Pam50 group i phenotype, No Overlap with Big Series subset

                ## Pam50 Subtypes - extract subtype dataframes
                if ( icinm == "All cases" ) {
                    Pam50iarrdf   <- sarrdf
                    Pam50iarrnovdf   <- sarrnovdf
                    sehPam50idf <- sehdf
                } else {
                    Pam50iarrdf   <- sarrdf[sarrdf$Pam50f == icinm, ]
                    Pam50iarrnovdf   <- sarrnovdf[sarrnovdf$Pam50f == icinm, ]
                    sehPam50idf <- sehdf[sehdf$Pam50f == icinm, ]
                }
                annodfNames <- names(annodf)
                annodfNamesFirst <- c("Probe_id", "Search_key", "Gene_symbol" )
                annodfNamesLast <- setdiff(annodfNames, annodfNamesFirst)
                arPam50ioutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]

                ## Set up matrix to hold regression results for fits to all probes
                arPam50ioutmatcolnames <-
                  c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
                arPam50ioutmat <- matrix(NA_real_, nrow = nrow(arPam50ioutdf), ncol = length(arPam50ioutmatcolnames))
                LE60_slope_col <- match("LE60_slope", arPam50ioutmatcolnames)
                GT60_slope_col <- match("GT60_slope", arPam50ioutmatcolnames)
                AllAges_slope_col <- match("AllAges_slope", arPam50ioutmatcolnames)
                LE60_pval_col <- match("LE60_pval", arPam50ioutmatcolnames)
                GT60_pval_col <- match("GT60_pval", arPam50ioutmatcolnames)
                AllAges_pval_col <- match("AllAges_pval", arPam50ioutmatcolnames)
                AgeDependentp_col <- match("AgeDependentp", arPam50ioutmatcolnames)
                arPam50ioutdf$LE60_slope <- rep(NA_real_ , nrow(arPam50ioutdf))
                arPam50ioutdf$GT60_slope <- rep(NA_real_ , nrow(arPam50ioutdf))
                arPam50ioutdf$AllAges_slope <- rep(NA_real_ , nrow(arPam50ioutdf))
                arPam50ioutdf$LE60_pval <- rep(NA_real_ , nrow(arPam50ioutdf))
                arPam50ioutdf$GT60_pval <- rep(NA_real_ , nrow(arPam50ioutdf))
                arPam50ioutdf$AllAges_pval <- rep(NA_real_ , nrow(arPam50ioutdf))
                arPam50ioutdf$AgeDependentp <- rep(NA , nrow(arPam50ioutdf))
                ## Biologically significant change in expression:  25% increase or decrease in fold change over time
                biosigslopeLE60 <-  log2(1.25)/35
                biosigslopeGT60 <-  log2(1.25)/35
                biosigslopeAllAges <-  log2(1.25)/70
                ## Case ids for this subset
                Pam50iidxs <- match(sehPam50idf$MBid, dimnames(Dataset.r)[[1]])
                ## Data frame for regressions
                arPam50ilmdf <- data.frame(age_at_diagnosis = sehPam50idf$age_at_diagnosis,
                                           probesetni = Dataset.r[Pam50iidxs, 1])
                ## Loop across all probes and fit regression lines - All cases
                cat("\n\n\n") 
                cat(icinm)  
                cat("\n\n\n") 
                ## Begin all 48k
                arPam50ioutdf_memoizedfilename <-
                  paste("AgeRelated_AllProbes_Pam50_",
                        gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                        "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                        "_All1992Cases_lm_v01.csv", sep = "")

                arPam50ioutdf_memoizedfilenamep <- FALSE

                if ( (!DoMemoizep) && file.exists(arPam50ioutdf_memoizedfilename) ) {
                    arPam50ioutdf <- read.table(file = arPam50ioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                sep = ",", header = TRUE, comment.char = "")
                    arPam50ioutdf_memoizedfilenamep <- TRUE
                } else {
                    for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
                        psi <- arPam50ioutdf$ProbeId[ni]  ## Probeset i
                        drci <- match(psi, dimnames(Dataset.r)[[2]])
                        arPam50ilmdf$probesetni <- Dataset.r[Pam50iidxs, drci]
                        if (ni %% 100 == 0 ) { cat(ni, ", ") }
                        lmifitageLE60 <- lm(probesetni ~ age_at_diagnosis, data = arPam50ilmdf,
                                            na.action = na.exclude, subset = age_at_diagnosis <= 60)
                        lmifitageGT60 <- lm(probesetni ~ age_at_diagnosis, data = arPam50ilmdf,
                                            na.action = na.exclude, subset = age_at_diagnosis > 60)
                        lmifitallages <- lm(probesetni ~ age_at_diagnosis, data = arPam50ilmdf,
                                            na.action = na.exclude)
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        smrylmifitallages <- summary(lmifitallages)
                        
                        arPam50ioutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                        arPam50ioutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                        arPam50ioutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                        arPam50ioutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arPam50ioutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arPam50ioutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                    }
                    
                    dimnames(arPam50ioutmat)[[2]] <- arPam50ioutmatcolnames
                    dimnames(arPam50ioutmat)[[1]] <- arPam50ioutdf$ProbeId
### End all 48k
                    arPam50ioutdf[, arPam50ioutmatcolnames] <- arPam50ioutmat
                    cat("\n\n\n")
                    ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                    arPam50ioutdf$LE60_BHadj_pval <- p.adjust(arPam50ioutdf$LE60_pval, method="BH")
                    arPam50ioutdf$GT60_BHadj_pval <- p.adjust(arPam50ioutdf$GT60_pval, method="BH")
                    arPam50ioutdf$AllAges_BHadj_pval <- p.adjust(arPam50ioutdf$AllAges_pval, method="BH")
                    arPam50ioutdf$BHadj_and_AgeDependentp <-
                      ( ( ( abs(arPam50ioutdf$LE60_slope) > biosigslopeLE60 )       & ( arPam50ioutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arPam50ioutdf$GT60_slope) > biosigslopeGT60 )       & ( arPam50ioutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arPam50ioutdf$AllAges_slope) > biosigslopeAllAges ) & ( arPam50ioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arPam50ioutdf$BHadj_signifp <- ( ( ( arPam50ioutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                       ( arPam50ioutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                       ( arPam50ioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arPam50ioutdf$LE60_log2FC <- arPam50ioutdf$LE60_slope * 35
                    arPam50ioutdf$GT60_log2FC <- arPam50ioutdf$GT60_slope * 35
                    arPam50ioutdf$AllAges_log2FC <- arPam50ioutdf$AllAges_slope * 70
                    
                    arPam50ioutdf$Best_log2FC <- arPam50ioutdf$AllAges_log2FC
                    arPam50ioutdf$Best_BHadj_pval <- arPam50ioutdf$AllAges_BHadj_pval
                    
                    ## Find the largest significant fold change:

                    LE60_Bestp <- ( ( abs(arPam50ioutdf$LE60_slope) > biosigslopeLE60 )          &
                                    ( arPam50ioutdf$LE60_BHadj_pval < FDRalpha ) )
                    arPam50ioutdf[LE60_Bestp, ]$Best_log2FC <- arPam50ioutdf[LE60_Bestp, ]$LE60_log2FC
                    arPam50ioutdf[LE60_Bestp, ]$Best_BHadj_pval <- arPam50ioutdf[LE60_Bestp, ]$LE60_BHadj_pval
                    
                    GT60_Bestp <- ( ( abs(arPam50ioutdf$GT60_slope) > biosigslopeGT60 )          &
                                    ( arPam50ioutdf$GT60_BHadj_pval < FDRalpha ) &
                                    ( abs(arPam50ioutdf$GT60_slope) > abs(arPam50ioutdf$LE60_slope) ) )
                    arPam50ioutdf[GT60_Bestp, ]$Best_log2FC <- arPam50ioutdf[GT60_Bestp, ]$GT60_log2FC
                    arPam50ioutdf[GT60_Bestp, ]$Best_BHadj_pval <- arPam50ioutdf[GT60_Bestp, ]$GT60_BHadj_pval
                    
                    AllAges_Bestp <- ( ( abs(arPam50ioutdf$AllAges_slope) > biosigslopeAllAges ) &
                                       ( arPam50ioutdf$AllAges_BHadj_pval < FDRalpha ) &
                                       ( ( abs(arPam50ioutdf$AllAges_log2FC) > abs(arPam50ioutdf$LE60_log2FC) ) |
                                         ( abs(arPam50ioutdf$AllAges_log2FC) > abs(arPam50ioutdf$GT60_log2FC) ) ) )
                    arPam50ioutdf[AllAges_Bestp, ]$Best_log2FC <- arPam50ioutdf[AllAges_Bestp, ]$AllAges_log2FC
                    arPam50ioutdf[AllAges_Bestp, ]$Best_BHadj_pval <- arPam50ioutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
                    arPam50ioutdf$Abs_Best_log2FC <- abs( arPam50ioutdf$Best_log2FC )
                    arPam50ioutdf$Abs_FoldChange <- 2^arPam50ioutdf$Abs_Best_log2FC
                    arPam50ioutdf$FoldChange_Direction <- ifelse(arPam50ioutdf$Best_log2FC > 0, "Up", "Down")
                    
                    ## Use this Best_BHadj_pval for manhattan plots as well.
                }

                ## No overlap 1161 cases:  No Overlap with Big Series and MB09 TMA cases
                ## Pam50inov = Pam50 Pam50i phenotype, from No Overlap with Big Series and MB09 TMA cases
                if ( icinm == "All cases" ) {
                    sehPam50inovdf <- sehnovdf
                } else {
                    sehPam50inovdf <- sehnovdf[sehnovdf$Pam50f == icinm, ]
                }
                arPam50inovoutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]
                arPam50inovoutmatcolnames <-
                  c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
                arPam50inovoutmat <- matrix(NA_real_, nrow = nrow(arPam50inovoutdf), ncol = length(arPam50inovoutmatcolnames))
                LE60_slope_col <- match("LE60_slope", arPam50inovoutmatcolnames)
                GT60_slope_col <- match("GT60_slope", arPam50inovoutmatcolnames)
                AllAges_slope_col <- match("AllAges_slope", arPam50inovoutmatcolnames)
                LE60_pval_col <- match("LE60_pval", arPam50inovoutmatcolnames)
                GT60_pval_col <- match("GT60_pval", arPam50inovoutmatcolnames)
                AllAges_pval_col <- match("AllAges_pval", arPam50inovoutmatcolnames)
                AgeDependentp_col <- match("AgeDependentp", arPam50inovoutmatcolnames)
                arPam50inovoutdf$LE60_slope <- rep(NA_real_ , nrow(arPam50inovoutdf))
                arPam50inovoutdf$GT60_slope <- rep(NA_real_ , nrow(arPam50inovoutdf))
                arPam50inovoutdf$AllAges_slope <- rep(NA_real_ , nrow(arPam50inovoutdf))
                arPam50inovoutdf$LE60_pval <- rep(NA_real_ , nrow(arPam50inovoutdf))
                arPam50inovoutdf$GT60_pval <- rep(NA_real_ , nrow(arPam50inovoutdf))
                arPam50inovoutdf$AllAges_pval <- rep(NA_real_ , nrow(arPam50inovoutdf))
                arPam50inovoutdf$AgeDependentp <- rep(NA , nrow(arPam50inovoutdf))
                biosigslopeLE60 <-  log2(1.25)/35
                biosigslopeGT60 <-  log2(1.25)/35
                biosigslopeAllAges <-  log2(1.25)/70
                Pam50inovidxs <- match(sehPam50inovdf$MBid, dimnames(Dataset.r)[[1]])
                arPam50inovlmdf <- data.frame(age_at_diagnosis = sehPam50inovdf$age_at_diagnosis,
                                              probesetni = Dataset.r[Pam50inovidxs, 1])

                cat("\n\n\n")
                cat(icinm)
                cat("\n\n\n")
                arPam50inovoutdf_memoizedfilename <-
                  paste("AgeRelated_AllProbes_Pam50_",
                        gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                        "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                        "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "")
                arPam50inovoutdf_memoizedfilenamep <- FALSE

                if ( (!DoMemoizep) && file.exists(arPam50inovoutdf_memoizedfilename) ) {
                    arPam50inovoutdf <- read.table(file = arPam50inovoutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                   sep = ",", header = TRUE, comment.char = "")
                    arPam50inovoutdf_memoizedfilenamep <- TRUE
                } else {
                    for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
                        psi <- arPam50inovoutdf$ProbeId[ni]
                        drci <- match(psi, dimnames(Dataset.r)[[2]])
                        arPam50inovlmdf$probesetni <- Dataset.r[Pam50inovidxs, drci]
                        if (ni %% 100 == 0 ) { cat(ni, ", ") }
                        lmifitageLE60 <-
                          lm(probesetni ~ age_at_diagnosis, data = arPam50inovlmdf, na.action = na.exclude, subset = age_at_diagnosis <= 60)
                        lmifitageGT60 <-
                          lm(probesetni ~ age_at_diagnosis, data = arPam50inovlmdf, na.action = na.exclude, subset = age_at_diagnosis > 60)
                        lmifitallages <-
                          lm(probesetni ~ age_at_diagnosis, data = arPam50inovlmdf, na.action = na.exclude)
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        smrylmifitallages <- summary(lmifitallages)
                        
                        arPam50inovoutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                        arPam50inovoutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                        arPam50inovoutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                        arPam50inovoutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arPam50inovoutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        arPam50inovoutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                    }

                    dimnames(arPam50inovoutmat)[[2]] <- arPam50inovoutmatcolnames
                    dimnames(arPam50inovoutmat)[[1]] <- arPam50inovoutdf$ProbeId
                    arPam50inovoutdf[, arPam50inovoutmatcolnames] <- arPam50inovoutmat
                    cat("\n\n\n")
                    ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                    arPam50inovoutdf$LE60_BHadj_pval <- p.adjust(arPam50inovoutdf$LE60_pval, method="BH")
                    arPam50inovoutdf$GT60_BHadj_pval <- p.adjust(arPam50inovoutdf$GT60_pval, method="BH")
                    arPam50inovoutdf$AllAges_BHadj_pval <- p.adjust(arPam50inovoutdf$AllAges_pval, method="BH")
                    arPam50inovoutdf$BHadj_and_AgeDependentp <-
                      ( ( ( abs(arPam50inovoutdf$LE60_slope) > biosigslopeLE60 )       &
                          ( arPam50inovoutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arPam50inovoutdf$GT60_slope) > biosigslopeGT60 )       &
                          ( arPam50inovoutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(arPam50inovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                          ( arPam50inovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arPam50inovoutdf$BHadj_signifp <- ( ( ( arPam50inovoutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                          ( arPam50inovoutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                          ( arPam50inovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                    
                    arPam50inovoutdf$LE60_log2FC <- arPam50inovoutdf$LE60_slope * 35
                    arPam50inovoutdf$GT60_log2FC <- arPam50inovoutdf$GT60_slope * 35
                    arPam50inovoutdf$AllAges_log2FC <- arPam50inovoutdf$AllAges_slope * 70

                    arPam50inovoutdf$Best_log2FC <- arPam50inovoutdf$AllAges_log2FC
                    arPam50inovoutdf$Best_BHadj_pval <- arPam50inovoutdf$AllAges_BHadj_pval

                    ## Find the largest significant fold change:

                    LE60_Bestp <- ( ( abs(arPam50inovoutdf$LE60_slope) > biosigslopeLE60 ) &
                                    ( arPam50inovoutdf$LE60_BHadj_pval < FDRalpha ) )
                    arPam50inovoutdf[LE60_Bestp, ]$Best_log2FC <- arPam50inovoutdf[LE60_Bestp, ]$LE60_log2FC
                    arPam50inovoutdf[LE60_Bestp, ]$Best_BHadj_pval <- arPam50inovoutdf[LE60_Bestp, ]$LE60_BHadj_pval

                    GT60_Bestp <- ( ( abs(arPam50inovoutdf$GT60_slope) > biosigslopeGT60 )          &
                                    ( arPam50inovoutdf$GT60_BHadj_pval < FDRalpha )                     &
                                    ( abs(arPam50inovoutdf$GT60_slope) > abs(arPam50inovoutdf$LE60_slope) ) )
                    arPam50inovoutdf[GT60_Bestp, ]$Best_log2FC <- arPam50inovoutdf[GT60_Bestp, ]$GT60_log2FC
                    arPam50inovoutdf[GT60_Bestp, ]$Best_BHadj_pval <- arPam50inovoutdf[GT60_Bestp, ]$GT60_BHadj_pval

                    AllAges_Bestp <- ( ( abs(arPam50inovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                                       ( arPam50inovoutdf$AllAges_BHadj_pval < FDRalpha )               &
                                       ( ( abs(arPam50inovoutdf$AllAges_log2FC) > abs(arPam50inovoutdf$LE60_log2FC) ) |
                                         ( abs(arPam50inovoutdf$AllAges_log2FC) > abs(arPam50inovoutdf$GT60_log2FC) ) ) )
                    arPam50inovoutdf[AllAges_Bestp, ]$Best_log2FC <- arPam50inovoutdf[AllAges_Bestp, ]$AllAges_log2FC
                    arPam50inovoutdf[AllAges_Bestp, ]$Best_BHadj_pval <- arPam50inovoutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
                    arPam50inovoutdf$Abs_Best_log2FC <- abs( arPam50inovoutdf$Best_log2FC )
                    arPam50inovoutdf$Abs_FoldChange <- 2^arPam50inovoutdf$Abs_Best_log2FC
                    arPam50inovoutdf$FoldChange_Direction <- ifelse(arPam50inovoutdf$Best_log2FC > 0, "Up", "Down")
                }

### Use this Best_BHadj_pval for manhattan plots as well.

### Plot scatterplots of gene sets:
### ar = age related study general gene set of interest for Tomo
### p  = polycomb independent
### prc2

### ar set:  Pam50i from all 1992 cases
### Pam50iarrdf <- sarrdf[sarrdf$Pam50f == icinm, ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                arPam50ioutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(arPam50ioutdf))
                arPam50ioutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(arPam50ioutdf))
                arPam50ioutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(arPam50ioutdf))
                arPam50ioutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(arPam50ioutdf))
                arPam50ioutdf$arGeneSet_PBHadj <- rep(NA, nrow(arPam50ioutdf))
                
                arPam50ioutdf$arGeneSetidxp <- ( arPam50ioutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                
                arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                  p.adjust(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                  p.adjust(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                  p.adjust(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                  apply(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp,
                                      c("arGeneSet_LE60_BHadj_pval", "arGeneSet_GT60_BHadj_pval", "arGeneSet_AllAges_BHadj_pval")], 1, min)

                arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                  ( ( ( abs(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
                      ( arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < FDRalpha ) ) |
                    ( ( abs(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
                      ( arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < FDRalpha ) ) |
                    ( ( abs(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                      ( arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                


### Pam50i All 1992 cases
### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in arPam50ioutdf
                Pam50iclfitageLE60 <- vector("list", length(names(arzidxs)))
                Pam50iclfitageGT60 <- vector("list", length(names(arzidxs)))
                Pam50iclfit <- vector("list", length(names(arzidxs)))
                stj <- paste("N =", dim(Pam50iarrdf)[1])

                ## Scatterplots for EZH2 H3K27me3 pathway genes arGeneSet
                if ( SingleOutputFilep && arScatterPlotsp ) {
                    pdf(file = paste("./Plots/MBEX_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                        width = 8, height = 10, useDingbats = FALSE)
                    par(mfrow = c(2, 2))
                }

                xlims <- c(20, 100)
                ylims <- c(4, 17)

                ageLE60idxp <- (Pam50iarrdf$age_at_diagnosis <= 60)
                ageGT60idxp <- (!ageLE60idxp)
                for ( ni in seq( along = names(arzidxs) ) ) {

                    i <- order(toupper(names(Pam50iarrdf)[seq(length(names(arzidxs)))]))[ni]
                    cat("\n\n### --- ", names(Pam50iarrdf)[i])
                    if ( !SingleOutputFilep && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/ProbeLevel/PAM50/MBEX_Pam50_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_", gsub("/", "_", names(Pam50iarrdf)[i]),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                            width = 6, height = 7, useDingbats = FALSE)
                        par(mfrow = c(1, 1))
                    }
                    lmifitageLE60 <- lm(Pam50iarrdf[ageLE60idxp, i] ~ Pam50iarrdf[ageLE60idxp, "age_at_diagnosis"])
                    smrylmifitageLE60 <- summary(lmifitageLE60)
                    Pam50iclfitageLE60[[i]] <- lmifitageLE60
                    if ( Verbosep ) print(smrylmifitageLE60)
                    lmifitageGT60 <- lm(Pam50iarrdf[ageGT60idxp, i] ~ Pam50iarrdf[ageGT60idxp, "age_at_diagnosis"])
                    smrylmifitageGT60 <- summary(lmifitageGT60)
                    Pam50iclfitageGT60[[i]] <- lmifitageGT60
                    if ( Verbosep ) print(smrylmifitageGT60)
                    lmifit <- lm(Pam50iarrdf[, i] ~ Pam50iarrdf[, "age_at_diagnosis"])
                    smrylmifit <- summary(lmifit)
                    Pam50iclfit[[i]] <- lmifit
                    if ( Verbosep ) print(smrylmifit)
                    if ( arScatterPlotsp ) {
                        plot(Pam50iarrdf[, "age_at_diagnosis"], Pam50iarrdf[, i],
                             ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                             main = "", xlim = xlims, ylim = ylims, type = "n")
                        title( main = paste(icinm, stj, sep = ": "), line = 2)
                        title( main =
                                 paste(names(Pam50iarrdf)[i], "(",
                                       ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(Pam50iarrdf)[i],
                                                                                   split = "\\|")[[1]][2], "ERbinding"],
                                              "ER binding", "non-ER binding"), ")"), line = 1)
                        points(Pam50iarrdf[, "age_at_diagnosis"], Pam50iarrdf[, i],
                               pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                        lines(supsmu(Pam50iarrdf[, "age_at_diagnosis"], Pam50iarrdf[, i], span = 0.4, bass = 10),
                              lwd = 3, col = "#AA1010")
                        abline(h = 0, v = 50, lty = 2)
                        text(x = 20, y = 16.3,
                             labels = paste("[<=60] p = ",
                                            format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 20, y = 15.5,
                             labels = paste("[>60] p = ",
                                            format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 98, y = 16.3,
                             labels = paste("[All] p = ",
                                            format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 1, cex = 0.85)
                        ## Need Pam50i subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                        if ( strsplit(names(Pam50iarrdf)[i], split = "\\|")[[1]][2] %in%
                             arPam50ioutdf[arPam50ioutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                            text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|>1.25 and adjPval<",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        } else {
                            text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels=paste("|FC|<1.25 or adjPval>",
                                              format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        }
                    }
                    if (!SingleOutputFilep && arScatterPlotsp ) { dev.off() }
                }
                if (SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arPam50ioutdf) %in% skipVarNms )
                if ( DoMemoizep || (! arPam50ioutdf_memoizedfilenamep ) ) {
                    write.csv(arPam50ioutdf[, -c(skipVarCols)], file = arPam50ioutdf_memoizedfilename )
                    write.csv(arPam50ioutdf[arPam50ioutdf$arGeneSetidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_arGeneSet_Pam50_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_All1992Cases_lm_v01.csv", sep = "") )
                }


### ar set:  Pam50i from 1161 no overlap cases (shown in asudf$MBid)
### Pam50iarrnovdf <- sarrdf[rownames(sarrdf) %in% asudf$MBid & sarrdf$Pam50Subtype == "Pam50i", ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                arPam50inovoutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(arPam50inovoutdf))
                arPam50inovoutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(arPam50inovoutdf))
                arPam50inovoutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(arPam50inovoutdf))
                arPam50inovoutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(arPam50inovoutdf))
                arPam50inovoutdf$arGeneSet_PBHadj <- rep(NA, nrow(arPam50inovoutdf))
                
                arPam50inovoutdf$arGeneSetidxp <- ( arPam50inovoutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                
                arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                  p.adjust(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                  p.adjust(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                  p.adjust(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                  apply(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp,
                                         c("arGeneSet_LE60_BHadj_pval",
                                           "arGeneSet_GT60_BHadj_pval",
                                           "arGeneSet_AllAges_BHadj_pval")], 1, min)

                
                arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                  ( ( ( abs(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$LE60_slope)    > biosigslopeLE60 )    &
                      ( arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval    < FDRalpha ) ) |
                    ( ( abs(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$GT60_slope)    > biosigslopeGT60 )    &
                      ( arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval    < FDRalpha ) ) |
                    ( ( abs(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                      ( arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                

### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in arPam50inovoutdf
                Pam50iclfitageLE60 <- vector("list", length(names(arzidxs)))
                Pam50iclfitageGT60 <- vector("list", length(names(arzidxs)))
                Pam50iclfit <- vector("list", length(names(arzidxs)))
                stj <- paste("N =", dim(Pam50iarrnovdf)[1], "(No overlap)")
                
                if ( SingleOutputFilep && arScatterPlotsp ) {
                    pdf(file = paste("./Plots/MBEX_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                        width = 8, height = 10, useDingbats = FALSE)
                    par(mfrow = c(2, 2))
                }
                xlims <- c(20, 100)
                ylims <- c(4, 17)

                ageLE60idxp <- (Pam50iarrnovdf$age_at_diagnosis <= 60)
                ageGT60idxp <- (!ageLE60idxp)
                for ( ni in seq( along = names(arzidxs) ) ) {
                    i <- order(toupper(names(Pam50iarrnovdf)[seq(length(names(arzidxs)))]))[ni]
                    cat("\n\n### --- ", names(Pam50iarrnovdf)[i])
                    
                    if ( (!SingleOutputFilep) && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/ProbeLevel/PAM50/MBEX_Pam50_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), 
                                         "_", gsub("/", "_", names(Pam50iarrnovdf)[i]), 
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                            width = 6, height = 7, useDingbats = FALSE)
                        par(mfrow = c(1, 1))
                    }

                    lmifitageLE60 <- lm(Pam50iarrnovdf[ageLE60idxp, i] ~ Pam50iarrnovdf[ageLE60idxp, "age_at_diagnosis"])
                    smrylmifitageLE60 <- summary(lmifitageLE60)
                    Pam50iclfitageLE60[[i]] <- lmifitageLE60
                    if ( Verbosep ) print(smrylmifitageLE60)
                    lmifitageGT60 <- lm(Pam50iarrnovdf[ageGT60idxp, i] ~ Pam50iarrnovdf[ageGT60idxp, "age_at_diagnosis"])
                    smrylmifitageGT60 <- summary(lmifitageGT60)
                    Pam50iclfitageGT60[[i]] <- lmifitageGT60
                    if ( Verbosep ) print(smrylmifitageGT60)
                    lmifit <- lm(Pam50iarrnovdf[, i] ~ Pam50iarrnovdf[, "age_at_diagnosis"])
                    smrylmifit <- summary(lmifit)
                    Pam50iclfit[[i]] <- lmifit
                    if ( Verbosep ) print(smrylmifit)
                    if ( arScatterPlotsp ) {
                        plot(Pam50iarrnovdf[, "age_at_diagnosis"], Pam50iarrnovdf[, i],
                             ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                             main = "", xlim = xlims, ylim = ylims, type = "n")
                        title( main = paste(icinm, stj, sep = ": "), line = 2)
                        title( main =
                                 paste(names(Pam50iarrnovdf)[i], "(",
                                       ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(Pam50iarrnovdf)[i],
                                                                                   split = "\\|")[[1]][2], "ERbinding"],
                                              "ER binding", "non-ER binding"), ")"), line = 1)
                        points(Pam50iarrnovdf[, "age_at_diagnosis"], Pam50iarrnovdf[, i],
                               pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                        lines(supsmu(Pam50iarrnovdf[, "age_at_diagnosis"], Pam50iarrnovdf[, i], span = 0.4, bass = 10),
                              lwd = 3, col = "#AA1010")
                        abline(h = 0, v = 50, lty = 2)
                        text(x = 20, y = 16.3,
                             labels = paste("[<=60] p = ",
                                            format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 20, y = 15.5,
                             labels = paste("[>60] p = ",
                                            format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 0, cex = 0.85)
                        text(x = 98, y = 16.3,
                             labels = paste("[All] p = ",
                                            format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                             adj = 1, cex = 0.85)
                        ## Need Pam50i subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                        if ( strsplit(names(Pam50iarrnovdf)[i], split = "\\|")[[1]][2] %in%
                             arPam50inovoutdf[arPam50inovoutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                            text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|>1.25 and adjPval<",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        } else {
                            text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                            text(x = 98, y = 15.0,
                                 labels = paste("|FC|<1.25 or adjPval>",
                                                format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                 adj = 1, cex = 0.85)
                        }
                    }
                    if ( (!SingleOutputFilep) && arScatterPlotsp ) { dev.off() }
                }
                if ( SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arPam50inovoutdf) %in% skipVarNms )
                if ( DoMemoizep || (! arPam50inovoutdf_memoizedfilenamep ) ) {
                    write.csv(arPam50inovoutdf[, -c(skipVarCols)],
                              file = arPam50inovoutdf_memoizedfilename )
                    write.csv(arPam50inovoutdf[arPam50inovoutdf$arGeneSetidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_arGeneSet_Pam50_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
                }



### Plot age-associated p-values by chromosome position - Manhattan plot

### qqman package has manhattan plot

                require("qqman")

### Data frame with CHR, BP, P (and SNP to avoid warnings)
### Genomic_location
### chr2:206352192:206352241:+
                arPam50ioutdf$Chrcr <- sapply(strsplit(arPam50ioutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                arPam50ioutdf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "",
                                                                 gsub("_h2_hap1", "", gsub("_random", "", arPam50ioutdf$Chrcr))))
                arPam50ioutdf$Chrn <- as.numeric(gsub("chr", "", arPam50ioutdf$Chrc))
                arPam50ioutdf$Chrn[arPam50ioutdf$Chrc == "chrX"] <- 23
                arPam50ioutdf$Chrn[arPam50ioutdf$Chrc == "chrY"] <- 24
                with(arPam50ioutdf, table(Chrcr, Chrn, useNA = "always"))
                with(arPam50ioutdf, table(Chrc, Chrn, useNA = "always"))
                arPam50ioutdf$CHR <- arPam50ioutdf$Chrn
                arPam50ioutdf$Startcr <- sapply(strsplit(arPam50ioutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                arPam50ioutdf$Stopcr <- sapply(strsplit(arPam50ioutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                arPam50ioutdf$Startn <- as.numeric(arPam50ioutdf$Startcr)
                arPam50ioutdf$Stopn <- as.numeric(arPam50ioutdf$Stopcr)
                arPam50ioutdf$BP <- trunc((arPam50ioutdf$Startn + arPam50ioutdf$Stopn)/2)
                arPam50ioutdf$SNP <- arPam50ioutdf$Probe_id
                arPam50ioutdf$P <- apply(arPam50ioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                arPam50ioutdf$PBHadj <- apply(arPam50ioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                arPam50ioutdfarGeneSetManhattanidxp <- arPam50ioutdf$arGeneSet_BHadj_and_AgeDependentp & arPam50ioutdf$Chrn <= 23
                arPam50ioutdfarGeneSetManhattanidxp[is.na(arPam50ioutdfarGeneSetManhattanidxp)] <- FALSE
                arPam50ioutdfManhattanidxp <- arPam50ioutdf$BHadj_and_AgeDependentp & arPam50ioutdf$Chrn <= 23
                arPam50ioutdfManhattanidxp[is.na(arPam50ioutdfManhattanidxp)] <- FALSE
                ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                arPam50ioutdfManhattanVarsidxp <-
                  ( apply(arPam50ioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                          function(x) !any(is.na(x)))  & ( arPam50ioutdf$Chrn <= 23 ) )

                arPam50ioutdf[arPam50ioutdf$Gene_symbol == "", "Gene_symbol"] <-
                  arPam50ioutdf[arPam50ioutdf$Gene_symbol == "", "ILMN_Gene_0"]
                ## Number of probesets showing evidence of association with age
                NpsPam50iAssocAge <- sum(arPam50ioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets"] <- NpsPam50iAssocAge
                ## Number of probesets showing evidence of association with age with genomic location
                NpsPam50iAssocAgecGL <- sum(arPam50ioutdfManhattanidxp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc"] <- NpsPam50iAssocAgecGL
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgPam50iAssocAgecGL <- length(unique(arPam50ioutdf[arPam50ioutdfManhattanidxp, "ILMN_Gene_0"]))
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames"] <- NgPam50iAssocAgecGL
                ## arGeneSet subset
                ## Number of probesets showing evidence of association with age
                NpsPam50iAssocAge_arGeneSet <- sum(arPam50ioutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_arGeneSet"] <- NpsPam50iAssocAge_arGeneSet
                ## Number of probesets showing evidence of association with age with genomic location
                NpsPam50iAssocAgecGL_arGeneSet <- sum(arPam50ioutdfarGeneSetManhattanidxp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <- NpsPam50iAssocAgecGL_arGeneSet
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgPam50iAssocAgecGL_arGeneSet <- length(unique(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames_arGeneSet"] <- NgPam50iAssocAgecGL_arGeneSet

                ## Non-overlapping with Big Series
                arPam50inovoutdf$Chrcr <- sapply(strsplit(arPam50inovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                arPam50inovoutdf$Chrc <- gsub("_qbl_hap2", "", gsub("_cox_hap1", "",
                                                                    gsub("_h2_hap1", "", gsub("_random", "", arPam50inovoutdf$Chrcr))))
                arPam50inovoutdf$Chrn <- as.numeric(gsub("chr", "", arPam50inovoutdf$Chrc))
                arPam50inovoutdf$Chrn[arPam50inovoutdf$Chrc == "chrX"] <- 23
                arPam50inovoutdf$Chrn[arPam50inovoutdf$Chrc == "chrY"] <- 24
                with(arPam50inovoutdf, table(Chrcr, Chrn, useNA = "always"))
                with(arPam50inovoutdf, table(Chrc, Chrn, useNA = "always"))
                arPam50inovoutdf$CHR <- arPam50inovoutdf$Chrn
                arPam50inovoutdf$Startcr <- sapply(strsplit(arPam50inovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                arPam50inovoutdf$Stopcr <- sapply(strsplit(arPam50inovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                arPam50inovoutdf$Startn <- as.numeric(arPam50inovoutdf$Startcr)
                arPam50inovoutdf$Stopn <- as.numeric(arPam50inovoutdf$Stopcr)
                arPam50inovoutdf$BP <- trunc((arPam50inovoutdf$Startn + arPam50inovoutdf$Stopn)/2)
                arPam50inovoutdf$SNP <- arPam50inovoutdf$Probe_id
                arPam50inovoutdf$P <- apply(arPam50inovoutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                arPam50inovoutdf$PBHadj <- apply(arPam50inovoutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                arPam50inovoutdfarGeneSetManhattanidxp <- arPam50inovoutdf$arGeneSet_BHadj_and_AgeDependentp & arPam50inovoutdf$Chrn <= 23
                arPam50inovoutdfarGeneSetManhattanidxp[is.na(arPam50inovoutdfarGeneSetManhattanidxp)] <- FALSE
                arPam50inovoutdfManhattanidxp <- arPam50inovoutdf$BHadj_and_AgeDependentp & arPam50inovoutdf$Chrn <= 23
                arPam50inovoutdfManhattanidxp[is.na(arPam50inovoutdfManhattanidxp)] <- FALSE
                ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                arPam50inovoutdfManhattanVarsidxp <-
                  ( apply(arPam50inovoutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                          function(x) !any(is.na(x)))  & ( arPam50inovoutdf$Chrn <= 23 ) )

                arPam50inovoutdf[arPam50inovoutdf$Gene_symbol == "", "Gene_symbol"] <-
                  arPam50inovoutdf[arPam50inovoutdf$Gene_symbol == "", "ILMN_Gene_0"]
                ## Number of probesets showing evidence of association with age
                NpsPam50inovAssocAge <- sum(arPam50inovoutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets"] <- NpsPam50inovAssocAge
                ## Number of probesets showing evidence of association with age with genomic location
                NpsPam50inovAssocAgecGL <- sum(arPam50inovoutdfManhattanidxp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc"] <- NpsPam50inovAssocAgecGL
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgPam50inovAssocAgecGL <- length(unique(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, "ILMN_Gene_0"]))
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames"] <- NgPam50inovAssocAgecGL
                ## arGeneSet subset
                ## Number of probesets showing evidence of association with age
                NpsPam50inovAssocAge_arGeneSet <- sum(arPam50inovoutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_arGeneSet"] <- NpsPam50inovAssocAge_arGeneSet
                ## Number of probesets showing evidence of association with age with genomic location
                NpsPam50inovAssocAgecGL_arGeneSet <- sum(arPam50inovoutdfarGeneSetManhattanidxp, na.rm = TRUE)
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <-
                  NpsPam50inovAssocAgecGL_arGeneSet
                ## Number of unique gene names corresponding to probesets assoc with age
                ## Note:  This is arbitrary as gene name variables are poorly maintained
                ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                ##        for consistency
                NgPam50inovAssocAgecGL_arGeneSet <-
                  length(unique(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                Pam50iAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames_arGeneSet"] <-
                  NgPam50inovAssocAgecGL_arGeneSet
                
                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arPam50ioutdf) %in% skipVarNms )

                if ( sum(arPam50ioutdfManhattanidxp) ) {
                    write.csv(arPam50ioutdf[arPam50ioutdfManhattanidxp, -c(skipVarCols)][
                        order(arPam50ioutdf[arPam50ioutdfManhattanidxp, "Abs_FoldChange"], decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.csv(arPam50ioutdf[arPam50ioutdfManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_Pam50_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.table(arPam50ioutdf[arPam50ioutdfManhattanidxp, ][
                        !duplicated(arPam50ioutdf[arPam50ioutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arPam50ioutdf[arPam50ioutdfManhattanidxp, ][
                            !duplicated(arPam50ioutdf[arPam50ioutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arPam50inovoutdf) %in% skipVarNms )

                if ( sum(arPam50inovoutdfManhattanidxp) ) {
                    write.csv(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, -c(skipVarCols)][
                        order(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.csv(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_Pam50_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.table(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, ][
                        !duplicated(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, ][
                            !duplicated(arPam50inovoutdf[arPam50inovoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_outdf_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                ## arGeneSet subset
                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arPam50ioutdf) %in% skipVarNms )
                if ( sum ( arPam50ioutdfarGeneSetManhattanidxp ) ) {
                    write.csv(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                        order(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.csv(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_Pam50_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_arGeneSet_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                    write.table(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, ][
                        !duplicated(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, ][
                            !duplicated(arPam50ioutdf[arPam50ioutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }

                skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                skipVarCols <- which( names(arPam50inovoutdf) %in% skipVarNms )
                if ( sum ( arPam50inovoutdfarGeneSetManhattanidxp ) ) {
                    write.csv(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                        order(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                              decreasing = TRUE), ][1:100, "Gene_symbol"],
                        file = paste("AgeRelated_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv",
                                     sep = "") )
                    write.csv(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                              file = paste("AgeRelated_Pam50_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_arGeneSet_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                    write.table(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, ][
                        !duplicated(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0[
                          !duplicated(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, ][
                            !duplicated(arPam50inovoutdf[arPam50inovoutdfarGeneSetManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                        file = paste("AgeRelated_AllMETABRIC_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                        row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                }
                
                ## Do stacked Manhattan and volcano plot for All1992Cases.  Plot log(FC) on X axis, log10(pval) on Y axis

                if ( !SingleOutputFilep ) {
                    pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 8, height = 8, useDingbats = FALSE)
                    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                } else {
                    pdf(file = paste("./Plots/AgeRelated_Manhattan_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                }

                ## Label the top 20 largest fold change with age association and adjusted p-val significance
                ohdf <- arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, ][
                    order(arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, "Abs_FoldChange"], decreasing = TRUE), ]
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
                
                if ( length( ObjsToHighlight ) ) {
                    sm_manhattan(arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj", chrlabs = c(1:22, "X"),
                                 suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                 highlight = ObjsToHighlight,
                                 plotpointidxp = arPam50ioutdf[arPam50ioutdfManhattanVarsidxp,
                                                               c( "BHadj_and_AgeDependentp")],
                                 plotpointhiliteidxp =
                                   if(sum(arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )

                } else {
                    sm_manhattan(arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj", chrlabs = c(1:22, "X"),
                                 suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                 plotpointidxp = arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],
                                 plotpointhiliteidxp =
                                   if(sum(arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arPam50ioutdf[arPam50ioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                }
                title(paste("METABRIC:", icinm), line = 2)
                title("Expression showing age-related trend, FC > 1.25", line = 1)

                write.csv(arPam50ioutdf[which(arPam50ioutdf$SNP %in% ObjsToHighlight), ],
                          file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_Pam50_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_All1992Cases.csv", sep = "") )
                
                if ( SingleOutputFilep ) {
                    dev.off()
                    pdf(file = paste("./Plots/AgeRelated_Volcano_Pam50_",
                                 if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   

                }
                ## Volcano
                ## Label top values
                ## Determine objects to highlight from volcano plot, label in both manhattan and volcano

                stairstepvals <-
                  switch(as.character(FDRalpha),
                         "0.05" = switch(as.character(icinm),
                                         "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                                         "LumA" = c(log2(3.5), -log10(FDRalpha), log2(2.5), 10, log2(1.25), 15), 
                                         "LumB" = c(log2(3), -log10(FDRalpha), log2(2), 4, log2(1.25), 5), 
                                         "Her2" = c(log2(4), -log10(FDRalpha), log2(3), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "Basal" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "Normal" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                         ),
                         "0.01" = switch(as.character(icinm),
                                         "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                                         "LumA" = c(log2(3.5), -log10(FDRalpha), log2(2.5), 10, log2(1.25), 13),
                                         "LumB" = c(log2(3), -log10(FDRalpha), log2(2), 4, log2(1.25), 5),
                                         "Her2" = c(log2(4), -log10(FDRalpha), log2(3), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                         "Basal" = c(log2(4), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/5), log2(1.25), 3.5), 
                                         "Normal" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha))
                                         )
                         )
                
                volcanoaddlidxp <- ( ( ( abs(arPam50ioutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                       ( -log10(arPam50ioutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                     ( ( abs(arPam50ioutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                       ( -log10(arPam50ioutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                     ( ( abs(arPam50ioutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                       ( -log10(arPam50ioutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                ObjsToHighlight <- arPam50ioutdf[volcanoaddlidxp, "SNP"]
                
                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                mainTitle <- paste("METABRIC BrCa", icinm)

                ERbindingAgeDep <- which( ( arPam50ioutdf$ERbinding &
                                            (abs(arPam50ioutdf$Best_log2FC ) > log2(FCthresh)) &
                                            (arPam50ioutdf$Best_BHadj_pval < alphacrit) ) )
                ERnonbindingAgeDep <- which( ( !arPam50ioutdf$ERbinding &
                                               (abs(arPam50ioutdf$Best_log2FC ) > log2(FCthresh)) &
                                               (arPam50ioutdf$Best_BHadj_pval < alphacrit) ) )
                NotAgeDepidxs <- setdiff(seq(nrow(arPam50ioutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                mainTitleN <- paste("Cases:", nrow(sehPam50idf),
                                    "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                
                if ( ! EqAxesScalesp ) {
                    if  (ici == 1) {
                        xlims_vpEVS <- max(abs(range(arPam50ioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                        ylims_vpEVS <- c( 0, 1.2*max(-log10(arPam50ioutdf$Best_BHadj_pval), na.rm = TRUE) )
                    } else {
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylimsi <- c( 0, 1.2*max(-log10(arPam50ioutdf$Best_BHadj_pval), na.rm = TRUE) )
                        ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                        ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                    }
                }
                ## Will need to cull data points at root of volcano - 48000 points too much.
                ## Get density information to thin out volcano plot
                ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                NADn <- length(NotAgeDepidxs)
                ClustRepsnRem <- NADn %% maxClusterN
                ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered
                

                for ( csi in seq(ClustRepsn) ) {
                    cat("   Cluster rep: ", csi, "  ")
                    if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                    csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                
                    densityEst <-
                      hclust(dist(cbind(arPam50ioutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                        jitter(-log10(arPam50ioutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                     method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                    lowDensityiidxs <-
                      (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityiidxs <-
                      sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                    if ( csi == 1 ) {
                        lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                        hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                    } else {
                        lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                        hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                    }
                }
                    
                plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                
                plot(arPam50ioutdf$Best_log2FC, -log10(arPam50ioutdf$Best_BHadj_pval), type = "n", main = "", 
                     xlim =  if ( EqAxesScalesp ) {
                                 xlims_vpEVS
                             } else {
                                 max(abs(range(arPam50ioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                             },
                     xlab = "(Down)        FC        (Up)   ",
                     ylab = "-log10(BH adjusted P-values)", xaxt = "n",
                     ylim = if ( EqAxesScalesp ) {
                                ylims_vpEVS
                            } else {
                                c( 0, 1.2 * max(-log10(arPam50ioutdf$Best_BHadj_pval), na.rm = TRUE) )
                            }
                     )

                title(main = mainTitle, line = 2)
                title(main = mainTitleN, line = 1)
                axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), labels = c(32, 16, 8, 4, 2, 1, 2, 4, 8, 16, 32) ) 
                

                points(arPam50ioutdf$Best_log2FC[NotAgeDepidxs][plotNotAgeDepidxs],
                       -log10(arPam50ioutdf$Best_BHadj_pval[NotAgeDepidxs][plotNotAgeDepidxs]),
                       pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))

                points(arPam50ioutdf$Best_log2FC[ERbindingAgeDep], -log10(arPam50ioutdf$Best_BHadj_pval[ERbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("red", alpha = FALSE)), alpha=30, max=255))
                points(arPam50ioutdf$Best_log2FC[ERnonbindingAgeDep], -log10(arPam50ioutdf$Best_BHadj_pval[ERnonbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("blue", alpha = FALSE)), alpha=30, max=255))
                legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                VolcanoAddlObjsToHighlight <- c()
                evcp <- arPam50ioutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                ## Cull out duplicate gene names
                evcidx <- unlist(tapply(which(evcp), arPam50ioutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                evcERbidx <- evcidx[which(arPam50ioutdf[evcidx, "ERbinding"])]
                evcERnbidx <- setdiff(evcidx, evcERbidx)
                lblcol <- rep("red", length(evcidx))
                lblcol[evcidx %in% evcERnbidx] <- "blue"
                if ( length(evcidx)  ) {
                    require("maptools"); ## for pointLabel
                    pointLabel(arPam50ioutdf[evcidx, ]$Best_log2FC, -log10(arPam50ioutdf[evcidx, ]$Best_BHadj_pval),
                               labels = arPam50ioutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                }
                ## Save ojects labeled for labeling in TCGA plots
                write.csv(arPam50ioutdf[evcidx, ],
                          file = paste("../Data/AgeRelated_Volcano_PointsToLabel_Pam50_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       gsub(" ", "",
                                            paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                           2^stairstepvals[3], stairstepvals[4],
                                                           2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_All1992Cases.csv", sep = "") )

                dev.off()


                ## Do stacked Manhattan and volcano plot for NoOverlapMB09BigSeries.  Plot log(FC) on X axis, log10(pval) on Y axis

                if ( !SingleOutputFilep ) {
                    pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 8, height = 8, useDingbats = FALSE)
                    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                } else {
                    pdf(file = paste("./Plots/AgeRelated_Manhattan_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                }

                ## Label the top 20 largest fold change with age association and adjusted p-val significance
                ohdf <- arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, ][
                    order(arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, "Abs_FoldChange"], decreasing = TRUE), ]
                ohdf <- ohdf[ohdf$BHadj_and_AgeDependentp, ]
                
                ObjsToHighlight <-
                  unique(c(ohdf[ohdf$FoldChange_Direction == "Up", ][1:10, "SNP"],
                           ohdf[ohdf$FoldChange_Direction == "Down", ][1:10, "SNP"],
                           ohdf[1:20, "SNP"]
                           ))

                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                if ( length( ObjsToHighlight ) ) {
                    sm_manhattan(arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj",
                                 chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha), genomewideline = FALSE, highlight = ObjsToHighlight,
                                 plotpointidxp = arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp =
                                   if(sum(arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                } else {
                    sm_manhattan(arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                 p="PBHadj",
                                 chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha), genomewideline = FALSE, 
                                 plotpointidxp = arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                 plotpointhiliteidxp =
                                   if(sum(arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                       arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                 hilitelbls =
                                   if(sum(arPam50inovoutdf[arPam50inovoutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                       c("ER binding", "Non binding") } else { NULL },
                                 ylab = "-log10(BH adjusted P-values)",
                                 main = "" )
                }
                title(paste("METABRIC:", icinm), line = 2)
                title("Expression showing age-related trend, FC > 1.25", line = 1)
                
                write.csv(arPam50inovoutdf[which(arPam50inovoutdf$SNP %in% ObjsToHighlight), ],
                          file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_Pam50_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), 
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_NoOverlapMB09BigSeries.csv", sep = "") )
                
                if ( SingleOutputFilep ) {
                    dev.off()
                    pdf(file = paste("./Plots/AgeRelated_Volcano_Pam50_",
                                 if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                        width = 6, height = 6, useDingbats = FALSE)
                    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1)) 
                }
                ## Volcano
                ## Determine objects to highlight from volcano plot, label in both manhattan and volcano

                stairstepvals <-
                  switch(as.character(FDRalpha),
                         "0.05" =
                           switch(as.character(icinm),
                                  "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                  "LumA" = c(log2(3), -log10(FDRalpha), log2(2), 5, log2(1.25), 7), 
                                  "LumB" = c(log2(3), -log10(FDRalpha), log2(2), 2, log2(1.25), 3), 
                                  "Her2" = c(log2(4), -log10(FDRalpha), log2(3), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                  "Basal" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                  "Normal" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5))
                                         ),
                         "0.01" =
                           switch(as.character(icinm),
                                  "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                  "LumA" = c(log2(3), -log10(FDRalpha), log2(2), 5, log2(1.25), 7.5), 
                                  "LumB" = c(log2(2.5), -log10(FDRalpha), log2(1.75), 2.5, log2(1.25), 3), 
                                  "Her2" = c(log2(4), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),
                                  "Basal" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                  "Normal" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha))
                                         )
                         )
                
                volcanoaddlidxp <- ( ( ( abs(arPam50inovoutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                       ( -log10(arPam50inovoutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                     ( ( abs(arPam50inovoutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                       ( -log10(arPam50inovoutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                     ( ( abs(arPam50inovoutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                       ( -log10(arPam50inovoutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                ObjsToHighlight <- arPam50inovoutdf[volcanoaddlidxp, "SNP"]

                ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                mainTitle <- paste("METABRIC BrCa", icinm)

                ERbindingAgeDep <- which( ( arPam50inovoutdf$ERbinding &
                                            (abs(arPam50inovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                            (arPam50inovoutdf$Best_BHadj_pval < alphacrit) ) )
                ERnonbindingAgeDep <- which( ( !arPam50inovoutdf$ERbinding &
                                               (abs(arPam50inovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                               (arPam50inovoutdf$Best_BHadj_pval < alphacrit) ) )
                NotAgeDepidxs <- setdiff(seq(nrow(arPam50inovoutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                mainTitleN <- paste("Cases:", nrow(sehPam50inovdf),
                                    "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                if ( ! EqAxesScalesp ) {
                    if  (nPam50i == 1) {
                        xlims_vpEVS <- max(abs(range(arPam50inovoutdf$Best_log2FC, na.rm = TRUE))) * c(-1, 1) * 1.1
                        ylims_vpEVS <- c( 0, 1.1 * max(-log10(arPam50inovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                    } else {
                        xlims_vpEVS <- c(-log2(48), log2(48))
                        ylimsi <- c( 0, 1.2 * max(-log10(arPam50inovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                        ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                        ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                    }
                }
                ## Will need to cull data points at root of volcano - 48000 points too much.
                ## Get density information to thin out volcano plot
                ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                NADn <- length(NotAgeDepidxs)
                ClustRepsnRem <- NADn %% maxClusterN
                ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered
                

                for ( csi in seq(ClustRepsn) ) {
                    cat("   Cluster rep: ", csi, "  ")
                    if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                    csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                
                    densityEst <-
                      hclust(dist(cbind(arPam50inovoutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                        jitter(-log10(arPam50inovoutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                     method = "single")
                    nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                    nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                    nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                    lowDensityiidxs <-
                      (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                    hiDensityiidxs <-
                      sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                    if ( csi == 1 ) {
                        lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                        hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                    } else {
                        lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                        hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                    }
                }
                    
                plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                
                plot(arPam50inovoutdf$Best_log2FC, -log10(arPam50inovoutdf$Best_BHadj_pval), type = "n", main = "", 
                     xlim =  if ( EqAxesScalesp ) {
                                 xlims_vpEVS
                             } else {
                                 max(abs(range(arPam50inovoutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                             },
                     xlab = "(Down)        FC        (Up)   ",
                     ylab = "-log10(BH adjusted P-values)", xaxt = "n",
                     ylim = if ( EqAxesScalesp ) {
                                ylims_vpEVS
                            } else {
                                c( 0, 1.2 * max(-log10(arPam50inovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                            }
                     )

                title(main = mainTitle, line = 2)
                title(main = mainTitleN, line = 1)
                axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), labels = c(32, 16, 8, 4, 2, 1, 2, 4, 8, 16, 32) )
                
                points(arPam50inovoutdf$Best_log2FC[NotAgeDepidxs][plotNotAgeDepidxs],
                       -log10(arPam50inovoutdf$Best_BHadj_pval[NotAgeDepidxs][plotNotAgeDepidxs]),
                       pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))

                points(arPam50inovoutdf$Best_log2FC[ERbindingAgeDep], -log10(arPam50inovoutdf$Best_BHadj_pval[ERbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("red", alpha = FALSE)), alpha=30, max=255))
                points(arPam50inovoutdf$Best_log2FC[ERnonbindingAgeDep], -log10(arPam50inovoutdf$Best_BHadj_pval[ERnonbindingAgeDep]),
                       pch = 20, col = rgb(t(col2rgb("blue", alpha = FALSE)), alpha=30, max=255))
                legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                VolcanoAddlObjsToHighlight <- c()
                evcp <- arPam50inovoutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                ## Cull out duplicate gene names
                evcidx <- unlist(tapply(which(evcp), arPam50inovoutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                evcERbidx <- evcidx[which(arPam50inovoutdf[evcidx, "ERbinding"])]
                evcERnbidx <- setdiff(evcidx, evcERbidx)
                lblcol <- rep("red", length(evcidx))
                lblcol[evcidx %in% evcERnbidx] <- "blue"
                if ( length(evcidx)  ) {
                    require("maptools"); ## for pointLabel
                    pointLabel(arPam50inovoutdf[evcidx, ]$Best_log2FC, -log10(arPam50inovoutdf[evcidx, ]$Best_BHadj_pval),
                               labels = arPam50inovoutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                }
                ## Save ojects labeled for labeling in TCGA plots
                write.csv(arPam50inovoutdf[evcidx, ],
                          file = paste("../Data/AgeRelated_Volcano_PointsToLabel_Pam50_",
                                       gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                       gsub(" ", "",
                                            paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                           2^stairstepvals[3], stairstepvals[4],
                                                           2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                       "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                       "_NoOverlapMB09BigSeries.csv", sep = "") )
                dev.off()
                ## End No overlap with Big Series or MB09
            }   ## End Pam50 i loop
        }   ## End  nPam50i loop

        Pam50iAgeAssociatedProbesetsGenesdf$N_WholeSeries <- as.vector(table(sehdf$Pam50f))
        Pam50iAgeAssociatedProbesetsGenesdf$N_NoOverlapBigSeries <- as.vector(table(sehnovdf$Pam50f))
        
        ## Write out data on counts of probesets/genes showing age association
        write.csv(Pam50iAgeAssociatedProbesetsGenesdf,
                  file = paste("AgeAssociatedProbesetGenes_Pam50_All",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               ".csv", sep = "") )
        
### AgeRelated_arGeneSet_xxx_All1992Cases_lm_v01.csv
### AgeRelated_arGeneSet_xxx_NoOverlapMB09BigSeries_lm_v01.csv
### Read in all these files, add Pam50 col, save as one file.

        for ( ici in seq(levels(sehdf$Pam50f)) ) {
            icinm <- levels(sehdf$Pam50f)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_All1992Cases_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$Pam50f <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_Pam50_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               "_All1992Cases_lm_v01.csv", sep = "") )

        
        for ( ici in seq(levels(sehnovdf$Pam50f)) ) {
            icinm <- levels(sehnovdf$Pam50f)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_Pam50_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     "_NoOverlapMB09BigSeries_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$Pam50f <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_Pam50_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
    }
}

#########################################################################
### Bootstrap bbbbbbbbbb
#########################################################################

################################################################################
### intClust groups

### iClust loop:  Bootstrap

### arGeneSet tag denotes 467 Illumina probesets corresponding to the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.
### Otherwise analysis is for all probesets.
### Set up output for age-associated transcript expression trends
### - across whole genome
### - across arGeneSet of 244 genes (467 probesets) for the EZH2 H3K27me3 relevant gene set compiled by Tomo Osako.

### iClust i across whole genome Bootstrap Loop Begin   ici <- 0
FDRalphalevels <- c(0.05, 0.01)
Verbosep <- FALSE ## TRUE
DoMemoizep <- TRUE
SingleOutputFilep <- TRUE
arScatterPlotsp <- TRUE
EqAxesScales <- c(FALSE, TRUE)  ## Must be in order F, T so appropriate scale range can be calculated across conditions
nNotAgeDepToPlot <- 1500
propLowDens <- 0.66
propHiDens <- (1.0 - propLowDens)
maxClusterN <- 15000
set.seed(3511)
Bboot <- 20
NiBoot <- c(300, 200, 150, 100, 75, 50)

for ( FDRalphai in seq(along = FDRalphalevels) ) {  ## FDRalphai <- 1

    for ( EqAxesScalesi in seq(along = EqAxesScales) ) {  ## EqAxesScalesi <- 1

        FDRalpha <- FDRalphalevels[FDRalphai]
        EqAxesScalesp <- EqAxesScales[EqAxesScalesi]

        niClusti <- length(levels(sehdf$iClustf))
        iClustiAgeAssociatedProbesetsGenesdf <-
          data.frame(iClusti = levels(sehdf$iClustf),
                     N_WholeSeries = rep(NA_integer_, niClusti),
                     N_WS_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                     N_WS_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                     N_WS_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),
                     N_NoOverlapBigSeries = rep(NA_integer_, niClusti),
                     N_NOv_AgeAssoc_Probesets = rep(NA_integer_, niClusti),
                     N_NOv_AgeAssoc_Probesets_cGenomicLoc = rep(NA_integer_, niClusti),
                     N_NOv_AgeAssoc_GeneNames = rep(NA_integer_, niClusti),
                     
                     N_WS_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                     N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                     N_WS_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti),
                     
                     N_NOv_AgeAssoc_Probesets_arGeneSet = rep(NA_integer_, niClusti),
                     N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet = rep(NA_integer_, niClusti),
                     N_NOv_AgeAssoc_GeneNames_arGeneSet = rep(NA_integer_, niClusti)
                     )
        
        for ( ici in seq( niClusti ) ) { ## ici <- 2
            
            for ( NiBi in seq(along = NiBoot) ) {  ## NiBi <- 1    Bootstrap sample sizes
                NiBooti <- NiBoot[NiBi]
                
                for ( Bbi in seq( Bboot ) ) {  ## Bbi <- 1  Bootstrap replicates
                    
                    icinm <- levels(sehdf$iClustf)[ici]
                    
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
### iClustinov = iClust group i phenotype, No Overlap with Big Series subset

                    ## iClust Subtypes - extract subtype dataframes
                    iClustiarrdf   <- sarrdf[sarrdf$iClustf == icinm, ]
                    iClustiarrnovdf   <- sarrnovdf[sarrnovdf$iClustf == icinm, ]
                    sehiClustidf <- sehdf[sehdf$iClustf == icinm, ]
                    ## No overlap 1161 cases:  No Overlap with Big Series and MB09 TMA cases
                    sehiClustinovdf <- sehnovdf[sehnovdf$iClustf == icinm, ]

                    ## Set up bootstrap data
                    bootorigdataidxs <- seq(nrow(iClustiarrdf))
                    bootisampleidxs <- sample(bootorigdataidxs, size = NiBooti, replace = TRUE)
                    bootorigdatanovidxs <- seq(nrow(iClustiarrnovdf))
                    bootisamplenovidxs <- sample(bootorigdatanovidxs, size = NiBooti, replace = TRUE)
                    iClustiarrdfbackup       <- iClustiarrdf
                    iClustiarrnovdfbackup    <- iClustiarrnovdf
                    sehiClustidfbackup       <- sehiClustidf
                    sehiClustinovdfbackup       <- sehiClustinovdf
                    iClustiarrdf   <- iClustiarrdfbackup[bootisampleidxs, ]
                    iClustiarrnovdf   <- iClustiarrnovdfbackup[bootisamplenovidxs, ]
                    sehiClustidf    <- sehiClustidfbackup[bootisampleidxs, ]
                    sehiClustinovdf <- sehiClustinovdfbackup[bootisamplenovidxs, ]
                    
                    annodfNames <- names(annodf)
                    annodfNamesFirst <- c("Probe_id", "Search_key", "Gene_symbol" )
                    annodfNamesLast <- setdiff(annodfNames, annodfNamesFirst)
                    ariClustioutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]

                    ## Set up matrix to hold regression results for fits to all probes
                    ariClustioutmatcolnames <-
                      c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
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
                    ## Biologically significant change in expression:  25% increase or decrease in fold change over time
                    biosigslopeLE60 <-  log2(1.25)/35
                    biosigslopeGT60 <-  log2(1.25)/35
                    biosigslopeAllAges <-  log2(1.25)/70
                    ## Case ids for this subset
                    iClustiidxs <- match(sehiClustidf$MBid, dimnames(Dataset.r)[[1]])
                    ## Data frame for regressions
                    ariClustilmdf <- data.frame(age_at_diagnosis = sehiClustidf$age_at_diagnosis,
                                                probesetni = Dataset.r[iClustiidxs, 1])
                    ## Loop across all probes and fit regression lines
                    cat("\n\n\n") 
                    cat(icinm)  
                    cat("\n\n\n") 
                    ## Begin all 48k
                    ariClustioutdf_memoizedfilename <-
                      paste("AgeRelated_AllProbes_",
                            gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                            "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                            "_Boot_", Bbi, 
                            "_All1992Cases_lm_v01.csv", sep = "")

                    ariClustioutdf_memoizedfilenamep <- FALSE

                    if ( (!DoMemoizep) && file.exists(ariClustioutdf_memoizedfilename) ) {
                        ariClustioutdf <- read.table(file = ariClustioutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                     sep = ",", header = TRUE, comment.char = "")
                        ariClustioutdf_memoizedfilenamep <- TRUE
                    } else {
                        for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) { ## ni <- 1
                            psi <- ariClustioutdf$ProbeId[ni]  ## Probeset i
                            drci <- match(psi, dimnames(Dataset.r)[[2]])
                            ariClustilmdf$probesetni <- Dataset.r[iClustiidxs, drci]
                            if (ni %% 100 == 0 ) { cat(ni, ", ") }
                            lmifitageLE60 <- lm(probesetni ~ age_at_diagnosis, data = ariClustilmdf,
                                                na.action = na.exclude, subset = age_at_diagnosis <= 60)
                            lmifitageGT60 <- lm(probesetni ~ age_at_diagnosis, data = ariClustilmdf,
                                                na.action = na.exclude, subset = age_at_diagnosis > 60)
                            lmifitallages <- lm(probesetni ~ age_at_diagnosis, data = ariClustilmdf,
                                                na.action = na.exclude)
                            smrylmifitageLE60 <- summary(lmifitageLE60)
                            smrylmifitageGT60 <- summary(lmifitageGT60)
                            smrylmifitallages <- summary(lmifitallages)
                            
                            ariClustioutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                            ariClustioutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                            ariClustioutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                            ariClustioutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                            ariClustioutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                            ariClustioutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        }
                        
                        dimnames(ariClustioutmat)[[2]] <- ariClustioutmatcolnames
                        dimnames(ariClustioutmat)[[1]] <- ariClustioutdf$ProbeId
### End all 48k
                        ariClustioutdf[, ariClustioutmatcolnames] <- ariClustioutmat
                        cat("\n\n\n")
                        ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                        ariClustioutdf$LE60_BHadj_pval <- p.adjust(ariClustioutdf$LE60_pval, method="BH")
                        ariClustioutdf$GT60_BHadj_pval <- p.adjust(ariClustioutdf$GT60_pval, method="BH")
                        ariClustioutdf$AllAges_BHadj_pval <- p.adjust(ariClustioutdf$AllAges_pval, method="BH")
                        ariClustioutdf$BHadj_and_AgeDependentp <-
                          ( ( ( abs(ariClustioutdf$LE60_slope) > biosigslopeLE60 )       &
                              ( ariClustioutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                            ( ( abs(ariClustioutdf$GT60_slope) > biosigslopeGT60 )       &
                              ( ariClustioutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                            ( ( abs(ariClustioutdf$AllAges_slope) > biosigslopeAllAges ) &
                              ( ariClustioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                        
                        ariClustioutdf$BHadj_signifp <- ( ( ( ariClustioutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                            ( ariClustioutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                            ( ariClustioutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                        
                        ariClustioutdf$LE60_log2FC <- ariClustioutdf$LE60_slope * 35
                        ariClustioutdf$GT60_log2FC <- ariClustioutdf$GT60_slope * 35
                        ariClustioutdf$AllAges_log2FC <- ariClustioutdf$AllAges_slope * 70
                        
                        ariClustioutdf$Best_log2FC <- ariClustioutdf$AllAges_log2FC
                        ariClustioutdf$Best_BHadj_pval <- ariClustioutdf$AllAges_BHadj_pval
                        
                        ## Find the largest significant fold change:

                        LE60_Bestp <- ( ( abs(ariClustioutdf$LE60_slope) > biosigslopeLE60 )          &
                                        ( ariClustioutdf$LE60_BHadj_pval < FDRalpha ) )
                        ariClustioutdf[LE60_Bestp, ]$Best_log2FC <- ariClustioutdf[LE60_Bestp, ]$LE60_log2FC
                        ariClustioutdf[LE60_Bestp, ]$Best_BHadj_pval <- ariClustioutdf[LE60_Bestp, ]$LE60_BHadj_pval
                        
                        GT60_Bestp <- ( ( abs(ariClustioutdf$GT60_slope) > biosigslopeGT60 )          &
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
                    }


                    ariClustinovoutdf <- annodf[, c(annodfNamesFirst, annodfNamesLast)]
                    ariClustinovoutmatcolnames <-
                      c("LE60_slope", "GT60_slope", "AllAges_slope", "LE60_pval", "GT60_pval", "AllAges_pval", "AgeDependentp")
                    ariClustinovoutmat <- matrix(NA_real_, nrow = nrow(ariClustinovoutdf), ncol = length(ariClustinovoutmatcolnames))
                    LE60_slope_col <- match("LE60_slope", ariClustinovoutmatcolnames)
                    GT60_slope_col <- match("GT60_slope", ariClustinovoutmatcolnames)
                    AllAges_slope_col <- match("AllAges_slope", ariClustinovoutmatcolnames)
                    LE60_pval_col <- match("LE60_pval", ariClustinovoutmatcolnames)
                    GT60_pval_col <- match("GT60_pval", ariClustinovoutmatcolnames)
                    AllAges_pval_col <- match("AllAges_pval", ariClustinovoutmatcolnames)
                    AgeDependentp_col <- match("AgeDependentp", ariClustinovoutmatcolnames)
                    ariClustinovoutdf$LE60_slope <- rep(NA_real_ , nrow(ariClustinovoutdf))
                    ariClustinovoutdf$GT60_slope <- rep(NA_real_ , nrow(ariClustinovoutdf))
                    ariClustinovoutdf$AllAges_slope <- rep(NA_real_ , nrow(ariClustinovoutdf))
                    ariClustinovoutdf$LE60_pval <- rep(NA_real_ , nrow(ariClustinovoutdf))
                    ariClustinovoutdf$GT60_pval <- rep(NA_real_ , nrow(ariClustinovoutdf))
                    ariClustinovoutdf$AllAges_pval <- rep(NA_real_ , nrow(ariClustinovoutdf))
                    ariClustinovoutdf$AgeDependentp <- rep(NA , nrow(ariClustinovoutdf))
                    biosigslopeLE60 <-  log2(1.25)/35
                    biosigslopeGT60 <-  log2(1.25)/35
                    biosigslopeAllAges <-  log2(1.25)/70
                    iClustinovidxs <- match(sehiClustinovdf$MBid, dimnames(Dataset.r)[[1]])
                    ariClustinovlmdf <- data.frame(age_at_diagnosis = sehiClustinovdf$age_at_diagnosis,
                                                   probesetni = Dataset.r[iClustinovidxs, 1])
                    ## Loop across all probes and fit regression lines
                    cat("\n\n\n")
                    cat(icinm)
                    cat("\n\n\n")
                    ariClustinovoutdf_memoizedfilename <-
                      paste("AgeRelated_AllProbes_",
                            gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                            "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                            "_Boot_", Bbi, 
                            "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "")
                    ariClustinovoutdf_memoizedfilenamep <- FALSE

                    if ( (!DoMemoizep) && file.exists(ariClustinovoutdf_memoizedfilename) ) {
                        ariClustinovoutdf <- read.table(file = ariClustinovoutdf_memoizedfilename, stringsAsFactors = FALSE,
                                                        sep = ",", header = TRUE, comment.char = "")
                        ariClustinovoutdf_memoizedfilenamep <- TRUE
                    } else {
                        for ( ni in seq(along = dimnames(Dataset.r)[[2]] ) ) {
                            psi <- ariClustinovoutdf$ProbeId[ni]
                            drci <- match(psi, dimnames(Dataset.r)[[2]])
                            ariClustinovlmdf$probesetni <- Dataset.r[iClustinovidxs, drci]
                            if (ni %% 100 == 0 ) { cat(ni, ", ") }
                            lmifitageLE60 <-
                              lm(probesetni ~ age_at_diagnosis, data = ariClustinovlmdf, na.action = na.exclude,
                                 subset = age_at_diagnosis <= 60)
                            lmifitageGT60 <-
                              lm(probesetni ~ age_at_diagnosis, data = ariClustinovlmdf, na.action = na.exclude,
                                 subset = age_at_diagnosis > 60)
                            lmifitallages <-
                              lm(probesetni ~ age_at_diagnosis, data = ariClustinovlmdf, na.action = na.exclude)
                            smrylmifitageLE60 <- summary(lmifitageLE60)
                            smrylmifitageGT60 <- summary(lmifitageGT60)
                            smrylmifitallages <- summary(lmifitallages)
                            
                            ariClustinovoutmat[ni, LE60_slope_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Estimate"]
                            ariClustinovoutmat[ni, GT60_slope_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Estimate"]
                            ariClustinovoutmat[ni, AllAges_slope_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Estimate"]
                            ariClustinovoutmat[ni, LE60_pval_col] <- smrylmifitageLE60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                            ariClustinovoutmat[ni, GT60_pval_col] <- smrylmifitageGT60$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                            ariClustinovoutmat[ni, AllAges_pval_col] <- smrylmifitallages$coefficients["age_at_diagnosis", "Pr(>|t|)"]
                        }

                        dimnames(ariClustinovoutmat)[[2]] <- ariClustinovoutmatcolnames
                        dimnames(ariClustinovoutmat)[[1]] <- ariClustinovoutdf$ProbeId
                        ariClustinovoutdf[, ariClustinovoutmatcolnames] <- ariClustinovoutmat
                        cat("\n\n\n")
                        ## Adjust all p-values and redo age-dependent selection on adjusted p-values < FDRalpha
                        ariClustinovoutdf$LE60_BHadj_pval <- p.adjust(ariClustinovoutdf$LE60_pval, method="BH")
                        ariClustinovoutdf$GT60_BHadj_pval <- p.adjust(ariClustinovoutdf$GT60_pval, method="BH")
                        ariClustinovoutdf$AllAges_BHadj_pval <- p.adjust(ariClustinovoutdf$AllAges_pval, method="BH")
                        ariClustinovoutdf$BHadj_and_AgeDependentp <-
                          ( ( ( abs(ariClustinovoutdf$LE60_slope) > biosigslopeLE60 )       &
                              ( ariClustinovoutdf$LE60_BHadj_pval    < FDRalpha ) ) |
                            ( ( abs(ariClustinovoutdf$GT60_slope) > biosigslopeGT60 )       &
                              ( ariClustinovoutdf$GT60_BHadj_pval    < FDRalpha ) ) |
                            ( ( abs(ariClustinovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                              ( ariClustinovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                        
                        ariClustinovoutdf$BHadj_signifp <- ( ( ( ariClustinovoutdf$LE60_BHadj_pval    < FDRalpha ) |
                                                               ( ariClustinovoutdf$GT60_BHadj_pval    < FDRalpha ) |
                                                               ( ariClustinovoutdf$AllAges_BHadj_pval < FDRalpha ) ) )
                        
                        ariClustinovoutdf$LE60_log2FC <- ariClustinovoutdf$LE60_slope * 35
                        ariClustinovoutdf$GT60_log2FC <- ariClustinovoutdf$GT60_slope * 35
                        ariClustinovoutdf$AllAges_log2FC <- ariClustinovoutdf$AllAges_slope * 70

                        ariClustinovoutdf$Best_log2FC <- ariClustinovoutdf$AllAges_log2FC
                        ariClustinovoutdf$Best_BHadj_pval <- ariClustinovoutdf$AllAges_BHadj_pval

                        ## Find the largest significant fold change:

                        LE60_Bestp <- ( ( abs(ariClustinovoutdf$LE60_slope) > biosigslopeLE60 ) &
                                        ( ariClustinovoutdf$LE60_BHadj_pval < FDRalpha ) )
                        ariClustinovoutdf[LE60_Bestp, ]$Best_log2FC <- ariClustinovoutdf[LE60_Bestp, ]$LE60_log2FC
                        ariClustinovoutdf[LE60_Bestp, ]$Best_BHadj_pval <- ariClustinovoutdf[LE60_Bestp, ]$LE60_BHadj_pval

                        GT60_Bestp <- ( ( abs(ariClustinovoutdf$GT60_slope) > biosigslopeGT60 )          &
                                        ( ariClustinovoutdf$GT60_BHadj_pval < FDRalpha )                     &
                                        ( abs(ariClustinovoutdf$GT60_slope) > abs(ariClustinovoutdf$LE60_slope) ) )
                        ariClustinovoutdf[GT60_Bestp, ]$Best_log2FC <- ariClustinovoutdf[GT60_Bestp, ]$GT60_log2FC
                        ariClustinovoutdf[GT60_Bestp, ]$Best_BHadj_pval <- ariClustinovoutdf[GT60_Bestp, ]$GT60_BHadj_pval

                        AllAges_Bestp <- ( ( abs(ariClustinovoutdf$AllAges_slope) > biosigslopeAllAges ) &
                                           ( ariClustinovoutdf$AllAges_BHadj_pval < FDRalpha )               &
                                           ( ( abs(ariClustinovoutdf$AllAges_log2FC) > abs(ariClustinovoutdf$LE60_log2FC) ) |
                                             ( abs(ariClustinovoutdf$AllAges_log2FC) > abs(ariClustinovoutdf$GT60_log2FC) ) ) )
                        ariClustinovoutdf[AllAges_Bestp, ]$Best_log2FC <- ariClustinovoutdf[AllAges_Bestp, ]$AllAges_log2FC
                        ariClustinovoutdf[AllAges_Bestp, ]$Best_BHadj_pval <- ariClustinovoutdf[AllAges_Bestp, ]$AllAges_BHadj_pval
                        ariClustinovoutdf$Abs_Best_log2FC <- abs( ariClustinovoutdf$Best_log2FC )
                        ariClustinovoutdf$Abs_FoldChange <- 2^ariClustinovoutdf$Abs_Best_log2FC
                        ariClustinovoutdf$FoldChange_Direction <- ifelse(ariClustinovoutdf$Best_log2FC > 0, "Up", "Down")
                    }

                    ## Use this Best_BHadj_pval for manhattan plots as well.

### Plot scatterplots of gene sets:
### ar = age related study general gene set of interest for Tomo
### p  = polycomb independent
### prc2

### ar set:  iClusti from all 1992 cases
### iClustiarrdf <- sarrdf[sarrdf$iClustf == icinm, ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                    ariClustioutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(ariClustioutdf))
                    ariClustioutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(ariClustioutdf))
                    ariClustioutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(ariClustioutdf))
                    ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(ariClustioutdf))
                    ariClustioutdf$arGeneSet_PBHadj <- rep(NA, nrow(ariClustioutdf))
                    
                    ariClustioutdf$arGeneSetidxp <- ( ariClustioutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                    
                    ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                      p.adjust(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                    ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                      p.adjust(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                    ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                      p.adjust(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                    ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                      apply(ariClustioutdf[ariClustioutdf$arGeneSetidxp,
                                           c("arGeneSet_LE60_BHadj_pval",
                                             "arGeneSet_GT60_BHadj_pval",
                                             "arGeneSet_AllAges_BHadj_pval")], 1, min)

                    ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                      ( ( ( abs(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$LE60_slope) > biosigslopeLE60 )       &
                          ( ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval < FDRalpha ) ) |
                        ( ( abs(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$GT60_slope) > biosigslopeGT60 )       &
                          ( ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval < FDRalpha ) ) |
                        ( ( abs(ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                          ( ariClustioutdf[ariClustioutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                    
### iClusti All 1992 cases
### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in ariClustioutdf
                    iClusticlfitageLE60 <- vector("list", length(names(arzidxs)))
                    iClusticlfitageGT60 <- vector("list", length(names(arzidxs)))
                    iClusticlfit <- vector("list", length(names(arzidxs)))
                    stj <- paste("N =", dim(iClustiarrdf)[1])

                    ## Scatterplots for EZH2 H3K27me3 pathway genes arGeneSet
                    if ( SingleOutputFilep && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/MBEX_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                            width = 8, height = 10, useDingbats = FALSE)
                        par(mfrow = c(2, 2))
                    }

                    xlims <- c(20, 100)
                    ylims <- c(4, 17)

                    ageLE60idxp <- (iClustiarrdf$age_at_diagnosis <= 60)
                    ageGT60idxp <- (!ageLE60idxp)
                    for ( ni in seq( along = names(arzidxs) ) ) {

                        i <- order(toupper(names(iClustiarrdf)[seq(length(names(arzidxs)))]))[ni]
                        cat("\n\n### --- ", names(iClustiarrdf)[i])
                        if ( !SingleOutputFilep && arScatterPlotsp ) {
                            pdf(file = paste("./Plots/ProbeLevel/MBEX_",
                                             gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                             "_", gsub("/", "_", names(iClustiarrdf)[i]),
                                             "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                             "_Boot_", Bbi, 
                                             "_arGeneSet_raw_All1992Cases_lm_v01.pdf", sep = ""),
                                width = 6, height = 7, useDingbats = FALSE)
                            par(mfrow = c(1, 1))
                        }
                        lmifitageLE60 <- lm(iClustiarrdf[ageLE60idxp, i] ~ iClustiarrdf[ageLE60idxp, "age_at_diagnosis"])
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        iClusticlfitageLE60[[i]] <- lmifitageLE60
                        if ( Verbosep ) print(smrylmifitageLE60)
                        lmifitageGT60 <- lm(iClustiarrdf[ageGT60idxp, i] ~ iClustiarrdf[ageGT60idxp, "age_at_diagnosis"])
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        iClusticlfitageGT60[[i]] <- lmifitageGT60
                        if ( Verbosep ) print(smrylmifitageGT60)
                        lmifit <- lm(iClustiarrdf[, i] ~ iClustiarrdf[, "age_at_diagnosis"])
                        smrylmifit <- summary(lmifit)
                        iClusticlfit[[i]] <- lmifit
                        if ( Verbosep ) print(smrylmifit)
                        if ( arScatterPlotsp ) {
                            plot(iClustiarrdf[, "age_at_diagnosis"], iClustiarrdf[, i],
                                 ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                                 main = "", xlim = xlims, ylim = ylims, type = "n")
                            title( main = paste(icinm, stj, sep = ": "), line = 2)
                            title( main =
                                     paste(names(iClustiarrdf)[i], "(",
                                           ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(iClustiarrdf)[i],
                                                                                       split = "\\|")[[1]][2], "ERbinding"],
                                                  "ER binding", "non-ER binding"), ")"), line = 1)
                            points(iClustiarrdf[, "age_at_diagnosis"], iClustiarrdf[, i],
                                   pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                            lines(supsmu(iClustiarrdf[, "age_at_diagnosis"], iClustiarrdf[, i], span = 0.4, bass = 10),
                                  lwd = 3, col = "#AA1010")
                            abline(h = 0, v = 50, lty = 2)
                            text(x = 20, y = 16.3,
                                 labels = paste("[<=60] p = ",
                                                format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                                 adj = 0, cex = 0.85)
                            text(x = 20, y = 15.5,
                                 labels = paste("[>60] p = ",
                                                format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                                 adj = 0, cex = 0.85)
                            text(x = 98, y = 16.3,
                                 labels = paste("[All] p = ",
                                                format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                                 adj = 1, cex = 0.85)
                            ## Need iClusti subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                            if ( strsplit(names(iClustiarrdf)[i], split = "\\|")[[1]][2] %in%
                                 ariClustioutdf[ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                                text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                                text(x = 98, y = 15.0,
                                     labels = paste("|FC|>1.25 and adjPval<",
                                                    format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                     adj = 1, cex = 0.85)
                            } else {
                                text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                                text(x = 98, y = 15.0,
                                     labels=paste("|FC|<1.25 or adjPval>",
                                                  format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                     adj = 1, cex = 0.85)
                            }
                        }
                        if (!SingleOutputFilep && arScatterPlotsp ) { dev.off() }
                    }
                    if (SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                    skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                    skipVarCols <- which( names(ariClustioutdf) %in% skipVarNms )
                    if ( DoMemoizep || (! ariClustioutdf_memoizedfilenamep ) ) {
                        write.csv(ariClustioutdf[, -c(skipVarCols)], file = ariClustioutdf_memoizedfilename )
                        write.csv(ariClustioutdf[ariClustioutdf$arGeneSetidxp, -c(skipVarCols)],
                                  file = paste("AgeRelated_arGeneSet_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_Boot_", Bbi, 
                                               "_All1992Cases_lm_v01.csv", sep = "") )
                    }


### ar set:  iClusti from 1161 no overlap cases (shown in asudf$MBid)
### iClustiarrnovdf <- sarrdf[rownames(sarrdf) %in% asudf$MBid & sarrdf$Pam50Subtype == "iClusti", ]

### Have to calculate statistical significance and biological relevance across the ar gene set (467 probes).
### Using results from the genome-wide set here is inappropriate.
                    ariClustinovoutdf$arGeneSet_LE60_BHadj_pval <- rep(NA_real_, nrow(ariClustinovoutdf))
                    ariClustinovoutdf$arGeneSet_GT60_BHadj_pval <- rep(NA_real_, nrow(ariClustinovoutdf))
                    ariClustinovoutdf$arGeneSet_AllAges_BHadj_pval <- rep(NA_real_, nrow(ariClustinovoutdf))
                    ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp <- rep(NA, nrow(ariClustinovoutdf))
                    ariClustinovoutdf$arGeneSet_PBHadj <- rep(NA, nrow(ariClustinovoutdf))
                    
                    ariClustinovoutdf$arGeneSetidxp <- ( ariClustinovoutdf$ProbeId %in% annodf$ProbeId[aridxs] )
                    
                    ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval <-
                      p.adjust(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$LE60_pval, method="BH")
                    ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval <-
                      p.adjust(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$GT60_pval, method="BH")
                    ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval <-
                      p.adjust(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$AllAges_pval, method="BH")
                    ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_PBHadj <-
                      apply(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp,
                                              c("arGeneSet_LE60_BHadj_pval",
                                                "arGeneSet_GT60_BHadj_pval",
                                                "arGeneSet_AllAges_BHadj_pval")], 1, min)

                    
                    ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_BHadj_and_AgeDependentp <-
                      ( ( ( abs(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$LE60_slope)    > biosigslopeLE60 )    &
                          ( ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_LE60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$GT60_slope)    > biosigslopeGT60 )    &
                          ( ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_GT60_BHadj_pval    < FDRalpha ) ) |
                        ( ( abs(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$AllAges_slope) > biosigslopeAllAges ) &
                          ( ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, ]$arGeneSet_AllAges_BHadj_pval < FDRalpha ) ) )
                    

### Regression for age <= 60, age > 60   i <- 10; j <- 4  Results in ariClustinovoutdf
                    iClusticlfitageLE60 <- vector("list", length(names(arzidxs)))
                    iClusticlfitageGT60 <- vector("list", length(names(arzidxs)))
                    iClusticlfit <- vector("list", length(names(arzidxs)))
                    stj <- paste("N =", dim(iClustiarrnovdf)[1], "(No overlap)")
                    
                    if ( SingleOutputFilep && arScatterPlotsp ) {
                        pdf(file = paste("./Plots/MBEX_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                            width = 8, height = 10, useDingbats = FALSE)
                        par(mfrow = c(2, 2))
                    }
                    xlims <- c(20, 100)
                    ylims <- c(4, 17)

                    ageLE60idxp <- (iClustiarrnovdf$age_at_diagnosis <= 60)
                    ageGT60idxp <- (!ageLE60idxp)
                    for ( ni in seq( along = names(arzidxs) ) ) {
                        i <- order(toupper(names(iClustiarrnovdf)[seq(length(names(arzidxs)))]))[ni]
                        cat("\n\n### --- ", names(iClustiarrnovdf)[i])
                        
                        if ( (!SingleOutputFilep) && arScatterPlotsp ) {
                            pdf(file = paste("./Plots/ProbeLevel/MBEX_",
                                             gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), 
                                             "_", gsub("/", "_", names(iClustiarrnovdf)[i]), 
                                             "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                             "_Boot_", Bbi, 
                                             "_arGeneSet_raw_NoOverlapMB09BigSeries_lm_v01.pdf", sep = ""),
                                width = 6, height = 7, useDingbats = FALSE)
                            par(mfrow = c(1, 1))
                        }

                        lmifitageLE60 <- lm(iClustiarrnovdf[ageLE60idxp, i] ~ iClustiarrnovdf[ageLE60idxp, "age_at_diagnosis"])
                        smrylmifitageLE60 <- summary(lmifitageLE60)
                        iClusticlfitageLE60[[i]] <- lmifitageLE60
                        if ( Verbosep ) print(smrylmifitageLE60)
                        lmifitageGT60 <- lm(iClustiarrnovdf[ageGT60idxp, i] ~ iClustiarrnovdf[ageGT60idxp, "age_at_diagnosis"])
                        smrylmifitageGT60 <- summary(lmifitageGT60)
                        iClusticlfitageGT60[[i]] <- lmifitageGT60
                        if ( Verbosep ) print(smrylmifitageGT60)
                        lmifit <- lm(iClustiarrnovdf[, i] ~ iClustiarrnovdf[, "age_at_diagnosis"])
                        smrylmifit <- summary(lmifit)
                        iClusticlfit[[i]] <- lmifit
                        if ( Verbosep ) print(smrylmifit)
                        if ( arScatterPlotsp ) {
                            plot(iClustiarrnovdf[, "age_at_diagnosis"], iClustiarrnovdf[, i],
                                 ylab = "Raw Expression (4, 16)", xlab = "Age at diagnosis",
                                 main = "", xlim = xlims, ylim = ylims, type = "n")
                            title( main = paste(icinm, stj, sep = ": "), line = 2)
                            title( main =
                                     paste(names(iClustiarrnovdf)[i], "(",
                                           ifelse(aroutdf[aroutdf$Probe_id == strsplit(names(iClustiarrnovdf)[i],
                                                                                       split = "\\|")[[1]][2], "ERbinding"],
                                                  "ER binding", "non-ER binding"), ")"), line = 1)
                            points(iClustiarrnovdf[, "age_at_diagnosis"], iClustiarrnovdf[, i],
                                   pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))
                            lines(supsmu(iClustiarrnovdf[, "age_at_diagnosis"], iClustiarrnovdf[, i], span = 0.4, bass = 10),
                                  lwd = 3, col = "#AA1010")
                            abline(h = 0, v = 50, lty = 2)
                            text(x = 20, y = 16.3,
                                 labels = paste("[<=60] p = ",
                                                format(smrylmifitageLE60$coefficients[2, 4], digits = 5), sep = ""),
                                 adj = 0, cex = 0.85)
                            text(x = 20, y = 15.5,
                                 labels = paste("[>60] p = ",
                                                format(smrylmifitageGT60$coefficients[2, 4], digits = 5), sep = ""),
                                 adj = 0, cex = 0.85)
                            text(x = 98, y = 16.3,
                                 labels = paste("[All] p = ",
                                                format(smrylmifit$coefficients[2, 4], digits = 5), sep = ""),
                                 adj = 1, cex = 0.85)
                            ## Need iClusti subset results for AgeDependent_BHadj_FC1.25_ProbeIds 
                            if ( strsplit(names(iClustiarrnovdf)[i], split = "\\|")[[1]][2] %in%
                                 ariClustinovoutdf[ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp, "ProbeId"] ) {
                                text(x = 98, y = 15.5, labels="Age-dependent trend:", adj = 1, cex = 0.85)
                                text(x = 98, y = 15.0,
                                     labels = paste("|FC|>1.25 and adjPval<",
                                                    format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                     adj = 1, cex = 0.85)
                            } else {
                                text(x = 98, y = 15.5, labels="No detectable trend:", adj = 1, cex = 0.85)
                                text(x = 98, y = 15.0,
                                     labels = paste("|FC|<1.25 or adjPval>",
                                                    format(FDRalpha, digits = (-log10(FDRalpha)) + 1), sep = ""),
                                     adj = 1, cex = 0.85)
                            }
                        }
                        if ( (!SingleOutputFilep) && arScatterPlotsp ) { dev.off() }
                    }
                    if ( SingleOutputFilep && arScatterPlotsp ) { dev.off() }

                    skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                    skipVarCols <- which( names(ariClustinovoutdf) %in% skipVarNms )
                    if ( DoMemoizep || (! ariClustinovoutdf_memoizedfilenamep ) ) {
                        write.csv(ariClustinovoutdf[, -c(skipVarCols)],
                                  file = ariClustinovoutdf_memoizedfilename )
                        write.csv(ariClustinovoutdf[ariClustinovoutdf$arGeneSetidxp, -c(skipVarCols)],
                                  file = paste("AgeRelated_arGeneSet_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_Boot_", Bbi, 
                                               "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
                    }

### Plot age-associated p-values by chromosome position - Manhattan plot
### qqman package has manhattan plot

                    require("qqman")

### Data frame with CHR, BP, P (and SNP to avoid warnings)
### Genomic_location
### chr2:206352192:206352241:+
                    ariClustioutdf$Chrcr <- sapply(strsplit(ariClustioutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                    ariClustioutdf$Chrc <- gsub("_qbl_hap2", "",
                                                gsub("_cox_hap1", "",
                                                     gsub("_h2_hap1", "",
                                                          gsub("_random", "", ariClustioutdf$Chrcr))))
                    ariClustioutdf$Chrn <- as.numeric(gsub("chr", "", ariClustioutdf$Chrc))
                    ariClustioutdf$Chrn[ariClustioutdf$Chrc == "chrX"] <- 23
                    ariClustioutdf$Chrn[ariClustioutdf$Chrc == "chrY"] <- 24
                    with(ariClustioutdf, table(Chrcr, Chrn, useNA = "always"))
                    with(ariClustioutdf, table(Chrc, Chrn, useNA = "always"))
                    ariClustioutdf$CHR <- ariClustioutdf$Chrn
                    ariClustioutdf$Startcr <- sapply(strsplit(ariClustioutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                    ariClustioutdf$Stopcr <- sapply(strsplit(ariClustioutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                    ariClustioutdf$Startn <- as.numeric(ariClustioutdf$Startcr)
                    ariClustioutdf$Stopn <- as.numeric(ariClustioutdf$Stopcr)
                    ariClustioutdf$BP <- trunc((ariClustioutdf$Startn + ariClustioutdf$Stopn)/2)
                    ariClustioutdf$SNP <- ariClustioutdf$Probe_id
                    ariClustioutdf$P <- apply(ariClustioutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                    ariClustioutdf$PBHadj <- apply(ariClustioutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                    ariClustioutdfarGeneSetManhattanidxp <- ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp & ariClustioutdf$Chrn <= 23
                    ariClustioutdfarGeneSetManhattanidxp[is.na(ariClustioutdfarGeneSetManhattanidxp)] <- FALSE
                    ariClustioutdfManhattanidxp <- ariClustioutdf$BHadj_and_AgeDependentp & ariClustioutdf$Chrn <= 23
                    ariClustioutdfManhattanidxp[is.na(ariClustioutdfManhattanidxp)] <- FALSE
                    ## Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                    ariClustioutdfManhattanVarsidxp <-
                      ( apply(ariClustioutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                              function(x) !any(is.na(x)))  & ( ariClustioutdf$Chrn <= 23 ) )

                    ariClustioutdf[ariClustioutdf$Gene_symbol == "", "Gene_symbol"] <-
                      ariClustioutdf[ariClustioutdf$Gene_symbol == "", "ILMN_Gene_0"]
                    ## Number of probesets showing evidence of association with age
                    NpsiClustiAssocAge <- sum(ariClustioutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets"] <- NpsiClustiAssocAge
                    ## Number of probesets showing evidence of association with age with genomic location
                    NpsiClustiAssocAgecGL <- sum(ariClustioutdfManhattanidxp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc"] <- NpsiClustiAssocAgecGL
                    ## Number of unique gene names corresponding to probesets assoc with age
                    ## Note:  This is arbitrary as gene name variables are poorly maintained
                    ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                    ##        for consistency
                    NgiClustiAssocAgecGL <- length(unique(ariClustioutdf[ariClustioutdfManhattanidxp, "ILMN_Gene_0"]))
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames"] <- NgiClustiAssocAgecGL
                    ## arGeneSet subset
                    ## Number of probesets showing evidence of association with age
                    NpsiClustiAssocAge_arGeneSet <- sum(ariClustioutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_arGeneSet"] <- NpsiClustiAssocAge_arGeneSet
                    ## Number of probesets showing evidence of association with age with genomic location
                    NpsiClustiAssocAgecGL_arGeneSet <- sum(ariClustioutdfarGeneSetManhattanidxp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <- NpsiClustiAssocAgecGL_arGeneSet
                    ## Number of unique gene names corresponding to probesets assoc with age
                    ## Note:  This is arbitrary as gene name variables are poorly maintained
                    ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                    ##        for consistency
                    NgiClustiAssocAgecGL_arGeneSet <- length(unique(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_WS_AgeAssoc_GeneNames_arGeneSet"] <- NgiClustiAssocAgecGL_arGeneSet

                    ## Non-overlapping with Big Series
                    ariClustinovoutdf$Chrcr <-
                      sapply(strsplit(ariClustinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[1])
                    ariClustinovoutdf$Chrc <-
                      gsub("_qbl_hap2", "",
                           gsub("_cox_hap1", "",
                                gsub("_h2_hap1", "",
                                     gsub("_random", "", ariClustinovoutdf$Chrcr))))
                    ariClustinovoutdf$Chrn <- as.numeric(gsub("chr", "", ariClustinovoutdf$Chrc))
                    ariClustinovoutdf$Chrn[ariClustinovoutdf$Chrc == "chrX"] <- 23
                    ariClustinovoutdf$Chrn[ariClustinovoutdf$Chrc == "chrY"] <- 24
                    with(ariClustinovoutdf, table(Chrcr, Chrn, useNA = "always"))
                    with(ariClustinovoutdf, table(Chrc, Chrn, useNA = "always"))
                    ariClustinovoutdf$CHR <- ariClustinovoutdf$Chrn
                    ariClustinovoutdf$Startcr <-
                      sapply(strsplit(ariClustinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[2])
                    ariClustinovoutdf$Stopcr <-
                      sapply(strsplit(ariClustinovoutdf$Genomic_location, split = ":"), function(x) unlist(x)[3])
                    ariClustinovoutdf$Startn <- as.numeric(ariClustinovoutdf$Startcr)
                    ariClustinovoutdf$Stopn <- as.numeric(ariClustinovoutdf$Stopcr)
                    ariClustinovoutdf$BP <- trunc((ariClustinovoutdf$Startn + ariClustinovoutdf$Stopn)/2)
                    ariClustinovoutdf$SNP <- ariClustinovoutdf$Probe_id
                    ariClustinovoutdf$P <-
                      apply(ariClustinovoutdf[, c("LE60_pval", "GT60_pval", "AllAges_pval")], 1, min)
                    ariClustinovoutdf$PBHadj <-
                      apply(ariClustinovoutdf[, c("LE60_BHadj_pval", "GT60_BHadj_pval", "AllAges_BHadj_pval")], 1, min)
                    ariClustinovoutdfarGeneSetManhattanidxp <-
                      ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp & ariClustinovoutdf$Chrn <= 23
                    ariClustinovoutdfarGeneSetManhattanidxp[is.na(ariClustinovoutdfarGeneSetManhattanidxp)] <- FALSE
                    ariClustinovoutdfManhattanidxp <-
                      ariClustinovoutdf$BHadj_and_AgeDependentp & ariClustinovoutdf$Chrn <= 23
                    ariClustinovoutdfManhattanidxp[is.na(ariClustinovoutdfManhattanidxp)] <- FALSE
### Get indexes with no missing data for the Manhattan variables to avoid problems with manhattan plot
                    ariClustinovoutdfManhattanVarsidxp <-
                      ( apply(ariClustinovoutdf[, c( "SNP", "CHR", "BP", "PBHadj")], 1,
                              function(x) !any(is.na(x)))  & ( ariClustinovoutdf$Chrn <= 23 ) )

                    ariClustinovoutdf[ariClustinovoutdf$Gene_symbol == "", "Gene_symbol"] <-
                      ariClustinovoutdf[ariClustinovoutdf$Gene_symbol == "", "ILMN_Gene_0"]
                    ## Number of probesets showing evidence of association with age
                    NpsiClustinovAssocAge <- sum(ariClustinovoutdf$BHadj_and_AgeDependentp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets"] <- NpsiClustinovAssocAge
                    ## Number of probesets showing evidence of association with age with genomic location
                    NpsiClustinovAssocAgecGL <- sum(ariClustinovoutdfManhattanidxp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc"] <-
                      NpsiClustinovAssocAgecGL
                    ## Number of unique gene names corresponding to probesets assoc with age
                    ## Note:  This is arbitrary as gene name variables are poorly maintained
                    ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                    ##        for consistency
                    NgiClustinovAssocAgecGL <- length(unique(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, "ILMN_Gene_0"]))
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames"] <- NgiClustinovAssocAgecGL
                    ## arGeneSet subset
                    ## Number of probesets showing evidence of association with age
                    NpsiClustinovAssocAge_arGeneSet <- sum(ariClustinovoutdf$arGeneSet_BHadj_and_AgeDependentp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_arGeneSet"] <- NpsiClustinovAssocAge_arGeneSet
                    ## Number of probesets showing evidence of association with age with genomic location
                    NpsiClustinovAssocAgecGL_arGeneSet <- sum(ariClustinovoutdfarGeneSetManhattanidxp, na.rm = TRUE)
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_Probesets_cGenomicLoc_arGeneSet"] <-
                      NpsiClustinovAssocAgecGL_arGeneSet
                    ## Number of unique gene names corresponding to probesets assoc with age
                    ## Note:  This is arbitrary as gene name variables are poorly maintained
                    ##        but since this is what we started using months ago for this number, use ILMN_Gene_0
                    ##        for consistency
                    NgiClustinovAssocAgecGL_arGeneSet <-
                      length(unique(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, "ILMN_Gene_0"]))
                    iClustiAgeAssociatedProbesetsGenesdf[ici, "N_NOv_AgeAssoc_GeneNames_arGeneSet"] <-
                      NgiClustinovAssocAgecGL_arGeneSet
                    
                    skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                    skipVarCols <- which( names(ariClustioutdf) %in% skipVarNms )

                    if ( sum(ariClustioutdfManhattanidxp) ) {
                        write.csv(ariClustioutdf[ariClustioutdfManhattanidxp, -c(skipVarCols)][
                            order(ariClustioutdf[ariClustioutdfManhattanidxp,
                                                 "Abs_FoldChange"], decreasing = TRUE), ][1:100, "Gene_symbol"],
                            file = paste("AgeRelated_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                        write.csv(ariClustioutdf[ariClustioutdfManhattanidxp, -c(skipVarCols)],
                                  file = paste("AgeRelated_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_Boot_", Bbi, 
                                               "_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                        write.table(ariClustioutdf[ariClustioutdfManhattanidxp,
                                                   ][!duplicated(ariClustioutdf[
                                                          ariClustioutdfManhattanidxp,
                                                          ]$Source_Reference_ID_0),
                                                     ]$ILMN_Gene_0[
                                                           !duplicated(ariClustioutdf[
                                                                ariClustioutdfManhattanidxp,
                                                                ][
                                                                !duplicated(ariClustioutdf[
                                                                     ariClustioutdfManhattanidxp,
                                                                     ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                                    file = paste("AgeRelated_AllMETABRIC_",
                                                 gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                                 "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                                 "_Boot_", Bbi, 
                                                 "_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                                    row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                    }

                    skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                    skipVarCols <- which( names(ariClustinovoutdf) %in% skipVarNms )

                    if ( sum(ariClustinovoutdfManhattanidxp) ) {
                        write.csv(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, -c(skipVarCols)][
                            order(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, "Abs_FoldChange"],
                                  decreasing = TRUE), ][1:100, "Gene_symbol"],
                            file = paste("AgeRelated_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                        write.csv(ariClustinovoutdf[ariClustinovoutdfManhattanidxp, -c(skipVarCols)],
                                  file = paste("AgeRelated_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_Boot_", Bbi, 
                                               "_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                        write.table(ariClustinovoutdf[
                            ariClustinovoutdfManhattanidxp, ][
                            !duplicated(ariClustinovoutdf[
                                 ariClustinovoutdfManhattanidxp,
                                 ]$Source_Reference_ID_0),
                            ]$ILMN_Gene_0[
                                  !duplicated(ariClustinovoutdf[
                                       ariClustinovoutdfManhattanidxp,
                                       ][
                                       !duplicated(ariClustinovoutdf[
                                            ariClustinovoutdfManhattanidxp, ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                            file = paste("AgeRelated_AllMETABRIC_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_outdf_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                            row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                    }

                    ## arGeneSet subset
                    skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                    skipVarCols <- which( names(ariClustioutdf) %in% skipVarNms )
                    if ( sum ( ariClustioutdfarGeneSetManhattanidxp ) ) {
                        write.csv(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                            order(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                                  decreasing = TRUE), ][1:100, "Gene_symbol"],
                            file = paste("AgeRelated_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                        write.csv(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                                  file = paste("AgeRelated_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_Boot_", Bbi, 
                                               "_arGeneSet_AgeDependent_and_BHadjSignificant_All1992Cases.csv", sep = "") )
                        write.table(ariClustioutdf[
                            ariClustioutdfarGeneSetManhattanidxp,
                            ][!duplicated(ariClustioutdf[
                                   ariClustioutdfarGeneSetManhattanidxp,
                                   ]$Source_Reference_ID_0),
                              ]$ILMN_Gene_0[
                                    !duplicated(ariClustioutdf[
                                         ariClustioutdfarGeneSetManhattanidxp, ][
                                         !duplicated(ariClustioutdf[ariClustioutdfarGeneSetManhattanidxp,
                                                                    ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                            file = paste("AgeRelated_AllMETABRIC_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_All1992Cases.txt", sep = ""),
                            row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                    }

                    skipVarNms <- c("Ontology_Component_0", "Ontology_Process_0",   "Ontology_Function_0" )
                    skipVarCols <- which( names(ariClustinovoutdf) %in% skipVarNms )
                    if ( sum ( ariClustinovoutdfarGeneSetManhattanidxp ) ) {
                        write.csv(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, -c(skipVarCols)][
                            order(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, "Abs_FoldChange"],
                                  decreasing = TRUE), ][1:100, "Gene_symbol"],
                            file = paste("AgeRelated_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_arGeneSet_GeneSymbol_Top100_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv",
                                         sep = "") )
                        write.csv(ariClustinovoutdf[ariClustinovoutdfarGeneSetManhattanidxp, -c(skipVarCols)],
                                  file = paste("AgeRelated_",
                                               gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                               "_Boot_", Bbi, 
                                               "_arGeneSet_AgeDependent_and_BHadjSignificant_NoOverlapMB09BigSeries.csv", sep = "") )
                        write.table(ariClustinovoutdf[
                            ariClustinovoutdfarGeneSetManhattanidxp,
                            ][!duplicated(ariClustinovoutdf[
                                   ariClustinovoutdfarGeneSetManhattanidxp,
                                   ]$Source_Reference_ID_0),
                              ]$ILMN_Gene_0[
                                    !duplicated(ariClustinovoutdf[
                                         ariClustinovoutdfarGeneSetManhattanidxp,
                                         ][
                                         !duplicated(ariClustinovoutdf[
                                              ariClustinovoutdfarGeneSetManhattanidxp,
                                              ]$Source_Reference_ID_0), ]$ILMN_Gene_0)],
                            file = paste("AgeRelated_AllMETABRIC_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_arGeneSet_AgeDependent_and_BHadjSignificant_ILMNGeneNames_NoOverlapMB09BigSeries.txt", sep = ""),
                            row.names = FALSE, col.names = FALSE, quote = FALSE) ##  age related gene symbols for Cytoscape
                    }
                    
                    ## Do stacked Manhattan and volcano plot for All1992Cases.  Plot log(FC) on X axis, log10(pval) on Y axis

                    if ( !SingleOutputFilep ) {
                        pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_All1992Cases_v01.pdf", sep = ""),
                            width = 8, height = 8, useDingbats = FALSE)
                        par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                    } else {
                        pdf(file = paste("./Plots/AgeRelated_Manhattan_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_All1992Cases_v01.pdf", sep = ""),
                            width = 6, height = 6, useDingbats = FALSE)
                        par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                    }

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
                    
                    if ( length( ObjsToHighlight ) ) {
                        sm_manhattan(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p="PBHadj", chrlabs = c(1:22, "X"),
                                     suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                     highlight = ObjsToHighlight,
                                     plotpointidxp = ariClustioutdf[ariClustioutdfManhattanVarsidxp,
                                                                    c( "BHadj_and_AgeDependentp")],  
                                     plotpointhiliteidxp =
                                       if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                           ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(BH adjusted P-values)",
                                     main = "" )
                    } else {
                        sm_manhattan(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p="PBHadj", chrlabs = c(1:22, "X"),
                                     suggestiveline = -log10(FDRalpha), genomewideline = FALSE,
                                     plotpointidxp = ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],
                                     plotpointhiliteidxp =
                                       if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0) {
                                           ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if(sum(ariClustioutdf[ariClustioutdfManhattanVarsidxp, c( "ERbinding" )], na.rm = TRUE) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(BH adjusted P-values)",
                                     main = "" )
                    }
                    title(paste("METABRIC:", icinm), line = 2)
                    title("Expression showing age-related trend, FC > 1.25", line = 1)

                    write.csv(ariClustioutdf[which(ariClustioutdf$SNP %in% ObjsToHighlight), ],
                              file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_Boot_", Bbi, 
                                           "_All1992Cases.csv", sep = "") )
                    
                    if ( SingleOutputFilep ) {
                        dev.off()
                        pdf(file = paste("./Plots/AgeRelated_Volcano_",
                                         if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_All1992Cases_v01.pdf", sep = ""),
                            width = 6, height = 6, useDingbats = FALSE)
                        par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   

                    }
                    ## Volcano
                    ## Label top values
                    ## Determine objects to highlight from volcano plot, label in both manhattan and volcano

                    stairstepvals <-
                      switch(as.character(FDRalpha),
                             "0.05" = switch(as.character(icinm),
                                             "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                                             "iClust 1" = c(log2(4), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                             "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 3" = c(log2(4), 1.8, log2(3), 3, log2(1.25), 5), 
                                             "iClust 4" = c(log2(3), -log10(FDRalpha), log2(2.5), 2.5, log2(1.25), 3.5), 
                                             "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                             "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                             "iClust 8" = c(log2(3), -log10(FDRalpha), log2(2.5), 4, log2(1.25), 5), 
                                             "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                                             "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                             ),
                             "0.01" = switch(as.character(icinm),
                                             "All cases" = c(log2(4), 10, log2(3), 25, log2(2), 40),
                                             "iClust 1" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 3" = c(log2(4), 2, log2(3), 3, log2(1.25), 5), 
                                             "iClust 4" = c(log2(4), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/5), log2(1.25), 3.5), 
                                             "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 8" = c(log2(4), -log10(FDRalpha), log2(3), 3, log2(1.25), 5), 
                                             "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                                             "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                             )
                             )
                    
                    volcanoaddlidxp <- ( ( ( abs(ariClustioutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                           ( -log10(ariClustioutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                         ( ( abs(ariClustioutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                           ( -log10(ariClustioutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                         ( ( abs(ariClustioutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                           ( -log10(ariClustioutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                    ObjsToHighlight <- ariClustioutdf[volcanoaddlidxp, "SNP"]
                    
                    ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                    mainTitle <- paste("METABRIC BrCa", icinm)

                    ERbindingAgeDep <- which( ( ariClustioutdf$ERbinding &
                                                (abs(ariClustioutdf$Best_log2FC ) > log2(FCthresh)) &
                                                (ariClustioutdf$Best_BHadj_pval < alphacrit) ) )
                    ERnonbindingAgeDep <- which( ( !ariClustioutdf$ERbinding &
                                                   (abs(ariClustioutdf$Best_log2FC ) > log2(FCthresh)) &
                                                   (ariClustioutdf$Best_BHadj_pval < alphacrit) ) )
                    NotAgeDepidxs <- setdiff(seq(nrow(ariClustioutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                    mainTitleN <- paste("Cases:", nrow(sehiClustidf),
                                        "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                    
                    if ( ! EqAxesScalesp ) {
                        if  (ici == 1) {
###                         xlims_vpEVS <- max(abs(range(ariClustioutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                            xlims_vpEVS <- c(-log2(48), log2(48))
                            ylims_vpEVS <- c( 0, 1.2*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
                        } else {
                            xlims_vpEVS <- c(-log2(48), log2(48))
                            ylimsi <- c( 0, 1.2*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
                            ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                            ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                        }
                    }
                    ## Will need to cull data points at root of volcano - 48000 points too much.
                    ## Get density information to thin out volcano plot
                    ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                    NADn <- length(NotAgeDepidxs)
                    ClustRepsnRem <- NADn %% maxClusterN
                    ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                    Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                    Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                    ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered

                    for ( csi in seq(ClustRepsn) ) {
                        cat("   Cluster rep: ", csi, "  ")
                        if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                        csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                        
                        densityEst <-
                          hclust(dist(cbind(ariClustioutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                            jitter(-log10(ariClustioutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                 method = "single")
                        nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                        nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                        nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                        lowDensityiidxs <-
                          (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                        hiDensityiidxs <-
                          sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                        if ( csi == 1 ) {
                            lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                            hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                        } else {
                            lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                            hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                        }
                    }
                    
                    plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                    cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                    
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
                                    c( 0, 1.2*max(-log10(ariClustioutdf$Best_BHadj_pval), na.rm = TRUE) )
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
                    legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                    VolcanoAddlObjsToHighlight <- c()
                    evcp <- ariClustioutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                    ## Cull out duplicate gene names
                    evcidx <- unlist(tapply(which(evcp), ariClustioutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                    evcERbidx <- evcidx[which(ariClustioutdf[evcidx, "ERbinding"])]
                    evcERnbidx <- setdiff(evcidx, evcERbidx)
                    lblcol <- rep("red", length(evcidx))
                    lblcol[evcidx %in% evcERnbidx] <- "blue"
                    if ( length(evcidx)  ) {
                        require("maptools"); ## for pointLabel
                        pointLabel(ariClustioutdf[evcidx, ]$Best_log2FC, -log10(ariClustioutdf[evcidx, ]$Best_BHadj_pval),
                                   labels = ariClustioutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                    }
                    write.csv(ariClustioutdf[evcidx, ],
                              file = paste("../Data/AgeRelated_Volcano_PointsToLabel_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                           gsub(" ", "",
                                                paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                               2^stairstepvals[3], stairstepvals[4],
                                                               2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_Boot_", Bbi, 
                                           "_All1992Cases.csv", sep = "") )

                    dev.off()


                    ## End All1992Cases

                    ## No overlap with Big Series or MB09
                    ## Do stacked Manhattan and volcano plot for NoOverlapMB09BigSeries.  Plot log(FC) on X axis, log10(pval) on Y axis

                    if ( !SingleOutputFilep ) {
                        pdf(file = paste("./Plots/AgeRelated_ManhattanAndVolcano_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                            width = 8, height = 8, useDingbats = FALSE)
                        par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
                    } else {
                        pdf(file = paste("./Plots/AgeRelated_Manhattan_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                            width = 6, height = 6, useDingbats = FALSE)
                        par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))   
                    }

                    ## Label the top 20 largest fold change with age association and adjusted p-val significance
                    ohdf <- ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, ][
                        order(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, "Abs_FoldChange"], decreasing = TRUE), ]
                    ohdf <- ohdf[ohdf$BHadj_and_AgeDependentp, ]
                    
                    ObjsToHighlight <-
                      unique(c(ohdf[ohdf$FoldChange_Direction == "Up", ][1:10, "SNP"],
                               ohdf[ohdf$FoldChange_Direction == "Down", ][1:10, "SNP"],
                               ohdf[1:20, "SNP"]
                               ))

                    ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                    if ( length( ObjsToHighlight ) ) {
                        sm_manhattan(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p="PBHadj",
                                     chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha),
                                     genomewideline = FALSE, highlight = ObjsToHighlight,
                                     plotpointidxp = ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                     plotpointhiliteidxp =
                                       if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                           ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(BH adjusted P-values)",
                                     main = "" )
                    } else {
                        sm_manhattan(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "SNP", "CHR", "BP", "PBHadj", "Gene_symbol")],
                                     p="PBHadj",
                                     chrlabs = c(1:22, "X"), suggestiveline = -log10(FDRalpha), genomewideline = FALSE, 
                                     plotpointidxp = ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "BHadj_and_AgeDependentp")],  
                                     plotpointhiliteidxp =
                                       if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0) {
                                           ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )] } else { NULL },
                                     hilitelbls =
                                       if(sum(ariClustinovoutdf[ariClustinovoutdfManhattanVarsidxp, c( "ERbinding" )]) > 0 ) {
                                           c("ER binding", "Non binding") } else { NULL },
                                     ylab = "-log10(BH adjusted P-values)",
                                     main = "" )
                    }
                    title(paste("METABRIC:", icinm), line = 2)
                    title("Expression showing age-related trend, FC > 1.25", line = 1)
                    
                    write.csv(ariClustinovoutdf[which(ariClustinovoutdf$SNP %in% ObjsToHighlight), ],
                              file = paste("../Data/AgeRelated_Manhattan_PointsToLabel_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_Boot_", Bbi, 
                                           "_NoOverlapMB09BigSeries.csv", sep = "") )
                    
                    if ( SingleOutputFilep ) {
                        dev.off()
                        pdf(file = paste("./Plots/AgeRelated_Volcano_",
                                         if ( EqAxesScalesp ) { "EqAxes_" } else { "" },
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_Boot_", Bbi, 
                                         "_NoOverlapMB09BigSeries_v01.pdf", sep = ""),
                            width = 6, height = 6, useDingbats = FALSE)
                        par(mfrow = c(1, 1), mar = c(4, 4, 3, 1)) 
                    }
                    ## Volcano
                    ## Determine objects to highlight from volcano plot, label in both manhattan and volcano

                    stairstepvals <-
                      switch(as.character(FDRalpha),
                             "0.05" = switch(as.character(icinm),
                                             "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                             "iClust 1" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 3" = c(log2(4), -log10(FDRalpha), log2(3), -log10(FDRalpha/5), log2(1.25), 2.5), 
                                             "iClust 4" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                             "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                             "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 8" = c(log2(3), -log10(FDRalpha), log2(2.5), 2, log2(1.25), 2.5), 
                                             "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                                             "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                             ),
                             "0.01" = switch(as.character(icinm),
                                             "All cases" = c(log2(4), 5, log2(3), 9, log2(2), 13),
                                             "iClust 1" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 2" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 3" = c(log2(4), -log10(FDRalpha), log2(2.5), 2, log2(1.25), 2.5),
                                             "iClust 4" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 5" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 6" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 7" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)), 
                                             "iClust 8" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha/2), log2(1.25), -log10(FDRalpha/5)), 
                                             "iClust 9" = c(log2(3), -log10(FDRalpha), log2(2.5), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)),   
                                             "iClust 10" = c(log2(2.25), -log10(FDRalpha), log2(1.75), -log10(FDRalpha), log2(1.25), -log10(FDRalpha)) 
                                             )
                             )
                    
                    volcanoaddlidxp <- ( ( ( abs(ariClustinovoutdf[, "Best_log2FC"]) > stairstepvals[1] ) &
                                           ( -log10(ariClustinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[2] ) ) |
                                         ( ( abs(ariClustinovoutdf[, "Best_log2FC"]) > stairstepvals[3] ) &
                                           ( -log10(ariClustinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[4] ) ) |
                                         ( ( abs(ariClustinovoutdf[, "Best_log2FC"]) > stairstepvals[5] ) &
                                           ( -log10(ariClustinovoutdf[, "Best_BHadj_pval"]) > stairstepvals[6] ) ) )
                    ObjsToHighlight <- ariClustinovoutdf[volcanoaddlidxp, "SNP"]

                    ObjsToHighlight <- ObjsToHighlight[!is.na(ObjsToHighlight)]

                    mainTitle <- paste("METABRIC BrCa", icinm)

                    ERbindingAgeDep <- which( ( ariClustinovoutdf$ERbinding &
                                                (abs(ariClustinovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                                (ariClustinovoutdf$Best_BHadj_pval < alphacrit) ) )
                    ERnonbindingAgeDep <- which( ( !ariClustinovoutdf$ERbinding &
                                                   (abs(ariClustinovoutdf$Best_log2FC ) > log2(FCthresh)) &
                                                   (ariClustinovoutdf$Best_BHadj_pval < alphacrit) ) )
                    NotAgeDepidxs <- setdiff(seq(nrow(ariClustinovoutdf)), c(ERbindingAgeDep, ERnonbindingAgeDep))

                    mainTitleN <- paste("Cases:", nrow(sehiClustinovdf),
                                        "  Age-dependent genes:", length(ERbindingAgeDep) + length(ERnonbindingAgeDep))

                    if ( ! EqAxesScalesp ) {
                        if  (ici == 1) {
###                     if  (niClusti == 1) {
###                         xlims_vpEVS <- max(abs(range(ariClustinovoutdf$Best_log2FC, na.rm = TRUE))) * c(-1, 1) * 1.1
                            xlims_vpEVS <- c(-log2(48), log2(48))
                            ylims_vpEVS <- c( 0, 1.2 * max(-log10(ariClustinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                        } else {
                            xlims_vpEVS <- c(-log2(48), log2(48))
                            ylimsi <- c( 0, 1.2 * max(-log10(ariClustinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                            ylims_vpEVS[1] <- min(ylims_vpEVS[1], ylimsi[1], na.rm = TRUE)
                            ylims_vpEVS[2] <- max(ylims_vpEVS[2], ylimsi[2], na.rm = TRUE)
                        }
                    }
                    ## Will need to cull data points at root of volcano - 48000 points too much.
                    ## Get density information to thin out volcano plot
                    ## Could split data into say 3 equal subsets, cluster each of those and combine results.
                    NADn <- length(NotAgeDepidxs)
                    ClustRepsnRem <- NADn %% maxClusterN
                    ClustRepsn <- if (ClustRepsnRem == 0) { NADn %/% maxClusterN } else { (NADn %/% maxClusterN) + 1 }
                    Clustni <- trunc( NADn / ClustRepsn ) ## Size for first few clusterings
                    Clustnlast <- NADn - ((ClustRepsn - 1) * Clustni) ## Size for final clustering
                    ClustOrdidxs <- sample(NADn) ## Random ordering of data to be clustered
                    

                    for ( csi in seq(ClustRepsn) ) {
                        cat("   Cluster rep: ", csi, "  ")
                        if ( csi == ClustRepsn ) { csin <- Clustnlast } else { csin <- Clustni }
                        csiidxs <- ClustOrdidxs[ (csi - 1) * Clustni + seq(csin) ]
                        
                        densityEst <-
                          hclust(dist(cbind(ariClustinovoutdf$Best_log2FC[NotAgeDepidxs][csiidxs],
                                            jitter(-log10(ariClustinovoutdf$Best_BHadj_pval[NotAgeDepidxs][csiidxs] ) ) ) ),
                                 method = "single")
                        nNotAgeDep <- sum(densityEst$merge < 0, na.rm = TRUE)
                        nLowDens <- trunc(propLowDens * ( min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep) ) )
                        nHiDens <- trunc(propHiDens * min(nNotAgeDepToPlot/ClustRepsn, nNotAgeDep))

                        lowDensityiidxs <-
                          (-1 * rev(densityEst$merge[densityEst$merge < 0])[seq(nLowDens)] )
                        hiDensityiidxs <-
                          sample( setdiff( (-1 * densityEst$merge[densityEst$merge < 0] ), lowDensityiidxs), size = nHiDens)
                        if ( csi == 1 ) {
                            lowDensityidxs <- NotAgeDepidxs[csiidxs][lowDensityiidxs]
                            hiDensityidxs <- NotAgeDepidxs[csiidxs][hiDensityiidxs]
                        } else {
                            lowDensityidxs <- c(lowDensityidxs,  NotAgeDepidxs[csiidxs][lowDensityiidxs])
                            hiDensityidxs <- c(hiDensityidxs,  NotAgeDepidxs[csiidxs][hiDensityiidxs])
                        }
                    }
                    
                    plotNotAgeDepidxs <- c(lowDensityidxs, hiDensityidxs)

                    cat("\nlength(unique(plotNotAgeDepidxs)): ",length(unique(plotNotAgeDepidxs)), "\n")
                    
                    plot(ariClustinovoutdf$Best_log2FC, -log10(ariClustinovoutdf$Best_BHadj_pval), type = "n", main = "", 
                         xlim =  if ( EqAxesScalesp ) {
                                     xlims_vpEVS
                                 } else {
                                     max(abs(range(ariClustinovoutdf$Best_log2FC, na.rm = TRUE)))*c(-1, 1)*1.1
                                 },
                         xlab = "(Down)        FC        (Up)   ",
                         ylab = "-log10(BH adjusted P-values)", xaxt = "n",
                         ylim = if ( EqAxesScalesp ) {
                                    ylims_vpEVS
                                } else {
                                    c( 0, 1.2 * max(-log10(ariClustinovoutdf$Best_BHadj_pval), na.rm = TRUE) )
                                }
                         )

                    title(main = mainTitle, line = 2)
                    title(main = mainTitleN, line = 1)
                    axis(side = 1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), labels = c(32, 16, 8, 4, 2, 1, 2, 4, 8, 16, 32) )
                    
                    points(ariClustinovoutdf$Best_log2FC[NotAgeDepidxs][plotNotAgeDepidxs],
                           -log10(ariClustinovoutdf$Best_BHadj_pval[NotAgeDepidxs][plotNotAgeDepidxs]),
                           pch = 20, col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255))

                    points(ariClustinovoutdf$Best_log2FC[ERbindingAgeDep], -log10(ariClustinovoutdf$Best_BHadj_pval[ERbindingAgeDep]),
                           pch = 20, col = rgb(t(col2rgb("red", alpha = FALSE)), alpha=30, max=255))
                    points(ariClustinovoutdf$Best_log2FC[ERnonbindingAgeDep], -log10(ariClustinovoutdf$Best_BHadj_pval[ERnonbindingAgeDep]),
                           pch = 20, col = rgb(t(col2rgb("blue", alpha = FALSE)), alpha=30, max=255))
                    legend("topleft", legend = c("ER binding", "Non binding"), col = c("red", "blue"), lwd = 3, pch = 1 )

                    VolcanoAddlObjsToHighlight <- c()
                    evcp <- ariClustinovoutdf$SNP  %in% c(ObjsToHighlight, VolcanoAddlObjsToHighlight)
                    ## Cull out duplicate gene names
                    evcidx <- unlist(tapply(which(evcp), ariClustinovoutdf[which(evcp), "Gene_symbol"], function(x) x[1]))
                    evcERbidx <- evcidx[which(ariClustinovoutdf[evcidx, "ERbinding"])]
                    evcERnbidx <- setdiff(evcidx, evcERbidx)
                    lblcol <- rep("red", length(evcidx))
                    lblcol[evcidx %in% evcERnbidx] <- "blue"
                    if ( length(evcidx)  ) {
                        require("maptools"); ## for pointLabel
                        pointLabel(ariClustinovoutdf[evcidx, ]$Best_log2FC, -log10(ariClustinovoutdf[evcidx, ]$Best_BHadj_pval),
                                   labels = ariClustinovoutdf[evcidx, ]$Gene_symbol, cex = 0.5, col = lblcol)
                    }
                    ## Save ojects labeled for labeling in TCGA plots
                    write.csv(ariClustinovoutdf[evcidx, ],
                              file = paste("../Data/AgeRelated_Volcano_PointsToLabel_",
                                           gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))), "_",
                                           gsub(" ", "",
                                                paste(format(c(2^stairstepvals[1], stairstepvals[2],
                                                               2^stairstepvals[3], stairstepvals[4],
                                                               2^stairstepvals[5], stairstepvals[6]), digits = 3), collapse = "_")),
                                           "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                           "_Boot_", Bbi, 
                                           "_NoOverlapMB09BigSeries.csv", sep = "") )
                    dev.off()
                    ## End No overlap with Big Series or MB09
                }
            }
            
        }   ## End iClust i loop

        iClustiAgeAssociatedProbesetsGenesdf$N_WholeSeries <- as.vector(table(sehdf$iClustf))
        iClustiAgeAssociatedProbesetsGenesdf$N_NoOverlapBigSeries <- as.vector(table(sehnovdf$iClustf))
        
        ## Write out data on counts of probesets/genes showing age association
        write.csv(iClustiAgeAssociatedProbesetsGenesdf,
                  file = paste("AgeAssociatedProbesetGenes_iClust_All",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               ## "_Boot_", Bbi, 
                               ".csv", sep = "") )
        
### AgeRelated_arGeneSet_xxx_All1992Cases_lm_v01.csv
### AgeRelated_arGeneSet_xxx_NoOverlapMB09BigSeries_lm_v01.csv
### Read in all these files, add iClust col, save as one file.

        for ( ici in seq(levels(sehdf$iClustf)) ) {
            icinm <- levels(sehdf$iClustf)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     ## "_Boot_", Bbi, 
                                     "_All1992Cases_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$iClustf <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_iClust_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               ## "_Boot_", Bbi, 
                               "_All1992Cases_lm_v01.csv", sep = "") )

        
        for ( ici in seq(levels(sehnovdf$iClustf)) ) {
            icinm <- levels(sehnovdf$iClustf)[ici]
            argsicidf <-
              read.csv( file = paste("AgeRelated_arGeneSet_",
                                     gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
                                     "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                     ## "_Boot_", Bbi, 
                                     "_NoOverlapMB09BigSeries_lm_v01.csv", sep = ""),
                       stringsAsFactors = FALSE)
            argsicidf$iClustf <- icinm
            if ( ici == 1 ) {
                argsicalldf <- argsicidf
            } else {
                argsicalldf <- rbind(argsicalldf, argsicidf)
            }
        }
        write.csv(argsicalldf,
                  file = paste("AgeRelated_arGeneSet_iClust_All_",
                               "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                               ## "_Boot_", Bbi, 
                               "_NoOverlapMB09BigSeries_lm_v01.csv", sep = "") )
    } ## End EqAxesScalesi loop
} ## End FDRalphai  loop


################################################################################

### Process Derek's bootstrap results
### /Users/stevenmckinney/git/github/BrCa_Age_Associated/main/code/Analysis_METABRIC_Expression_Trends
### AgeDependent_BootstrapStats_NB300_Rep01.rds   Reps 01 to 20

AgeDependent_BootstrapStats_NB300_Rep01 <-
  readRDS(file = paste("/Users/stevenmckinney/git/github",
                       "/BrCa_Age_Associated/main/code/Analysis_METABRIC_Expression_Trends",
                       "/", "AgeDependent_BootstrapStats_NB300_Rep01.rds", sep = "") )

### Compare top genes - large abs fold change and age associated trend
### How stable are gene lists?
### Get top 100 genes from each bootstrap sample
### intersect across bootstrap samples
### Intersect each bootstrap sample with observed sample and record proportion overlap
### Memoized:  AgeRelated_AllProbes_iClust_2_FDR_0p01_All1992Cases_lm_v01.csv

Ntop <- 100
FDRalpha <- 0.01
Bboot <- 20
NiBoot <- 300

niClusti <- length(levels(sehdf$iClustf))
for ( ici in seq( niClusti ) ) { ## ici <- 1
    icinm <- levels(sehdf$iClustf)[ici]
    cat("\n\n", icinm, "\n\n")

    ariClustioutdf_memoizedfilename <-
      paste("AgeRelated_AllProbes_",
            gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "n", gsub("\\+", "p", icinm)))),
            "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
            "_All1992Cases_lm_v01.csv", sep = "")

    for ( Bbi in seq( Bboot ) ) {
    require("MOutils")
        AgeDependent_BootstrapStats_NBNiBoot_RepBbi <-
          readRDS(file = paste("/Users/stevenmckinney/git/github",
                               "/BrCa_Age_Associated/main/code/Analysis_METABRIC_Expression_Trends",
                               "/", "AgeDependent_BootstrapStats_NB", NiBoot, "_Rep",
                               leftPad0(Bbi, 2), ".rds", sep = "") )

        }

}
### -----------------------------------------------------------------------------
### 
### Create ESR1 graphs with 99% CIs to illustrate intClust group variability
### ER−alpha | ILMN_1678535   sehrdf
###
### Bootstrap and generate 99% CI from supsmu fits.  If CI contains horizontal line thru raw expression mean
### then no evidence of association.

###  ER-alpha|ILMN_1678535

### For each iClust group, fit one regression line across all ages and plot with 99%CI

FDRalpha <- 0.01

sehrdf$ESR1 <- sehrdf[, "ER-alpha|r|ILMN_1678535"]
niClusti <- length(levels(sehrdf$iClustf))
for ( ici in seq( niClusti ) ) { ## ici <- 2
    icinm <- levels(sehrdf$iClustf)[ici]
    iClustiarrdf <- sehrdf[which(sehrdf$iClustf == icinm), ]
    cat("\n\n", icinm, "\n\n")
###     lmifitallages <- lm(ESR1 ~ age_at_diagnosis, data = iClustiarrdf,
###                         na.action = na.exclude)
    library("ggplot2")
    pdf(file = paste("./Plots/MBEX_lm_ESR1_",
                                         gsub(" ", "_", gsub( "\\/", "", gsub( "\\-", "neg", gsub("\\+", "pos", icinm)))),
                                         "_FDR_", gsub("\\.", "p", format(FDRalpha, digits = 5)),
                                         "_arGeneSet_raw_All1992Cases.pdf", sep = ""),
        width = 8, height = 10, useDingbats = FALSE)
    tsize <- 4.0
    print(ggplot(iClustiarrdf, aes(x=age_at_diagnosis, y=ESR1)) + ylim(4, 14) + xlim(25, 95) +
          geom_point(col = rgb(t(col2rgb("black", alpha = FALSE)), alpha=30, max=255)) +
          geom_smooth(method=lm, se=TRUE, level = (1.0 - FDRalpha/niClusti), col = "#AA1010") +
          geom_hline(aes(yintercept=mean(iClustiarrdf$ESR1)), linetype = 2) +
          annotate("text", x = 32, y = 13.5, size = tsize, label = icinm) +
          annotate("text", x = 32, y = 13.2, size = tsize, label = paste("N =", nrow(iClustiarrdf))) +
          annotate("text", x = 32, y = 12.9, size = tsize, label = paste("alpha =", format(FDRalpha/niClusti, digits = 5))) +
          annotate("text", x = 32, y = 12.6, size = tsize, label = paste("FDR =", format(FDRalpha, digits = 5)))
          )
    dev.off()
}


### ------------------------------------------------------
###
### Merge hg19 info to annodf
###
###
###
### ------------------------------------------------------

head -n 500000 /Volumes/KilroyHD/kilroy/Projects/GenomicSignatures/GRCh38/GRCh38_latest_genomic.gff | tail

h38df <- read.table(file = "/Volumes/KilroyHD/kilroy/Projects/GenomicSignatures/GRCh38/GRCh38_latest_genomic.gff",
                    sep = "\t", header = FALSE, skip = 9, quote = "\"", comment.char = "",
                    fill = TRUE, stringsAsFactors = FALSE, colClasses = "character")
V9list <- lapply(h38df$V9, function(x) strsplit(x, split = ";"))
efn <- function(x, mstr = "gene=") {
    xv <- unlist(strsplit(x, split = ";"))
    gsub(mstr, "", xv[charmatch(mstr, xv)])
}
V9gn <- sapply(h38df$V9, egn, USE.NAMES = FALSE)
V9ID <- sapply(h38df$V9, efn, mstr = "ID=", USE.NAMES = FALSE)
V9IDc <- gsub("[0-9]", "", V9ID)
mgv <- match("gene", V9IDc)
length(mgv)  ## [1] 60178
length(unique(V9gn[mgv])) ## [1] 54267
length(V9gn)
nrow(h38df)
length(unique(V9gn))  ## [1] 54268
tail(V9ID)  ## [1] "id1871770" "gene60177" "rna165305" "id1871771" "id1871772" NA 
## "gene60177"
## So there's nothing close to 20000 for gene count in h38 genomic annotation.
## There's more codes here than probesets in Illumina mRNA chip data.
## Identifying genes appears nowhere near any exact science here.
### > length(unique(annodf[, "Ensembl_gene_id"]))
### [1] 19410
### This appears to be the only count near 20000

ens38df <- 
  read.table(file = "/Volumes/KilroyHD/kilroy/Projects/GenomicSignatures/GRCh38/Ensembl/Homo_sapiens.GRCh38.89.gff3",
             sep = "\t", header = FALSE, skip = 215, quote = "\"", comment.char = "",
             fill = TRUE, stringsAsFactors = FALSE, colClasses = "character")

ens38dfV9ID <- sapply(ens38df$V9, efn, mstr = "ID=gene:", USE.NAMES = FALSE)
ens38df$V9IDgene <- sapply(ens38df$V9, efn, mstr = "ID=gene:", USE.NAMES = FALSE)
ens38df$V9IDName <- sapply(ens38df$V9, efn, mstr = "Name=", USE.NAMES = FALSE)
ens38df$V9IDbiotype <- sapply(ens38df$V9, efn, mstr = "biotype=", USE.NAMES = FALSE)
ens38df$V9IDlogic_name <- sapply(ens38df$V9, efn, mstr = "logic_name=", USE.NAMES = FALSE)
length(unique(ens38df$V9IDbiotype))
length(unique(ens38df$V9IDlogic_name))
length(unique(ens38df$V3))

length(unique( ens38df[which(ens38df$V9IDlogic_name == "ensembl_havana_gene"), "V9IDName"]))

head(ens38df[which(ens38df$V9IDlogic_name == "ensembl"), ])

table(  ens38df[which(ens38df$V9IDlogic_name == "ensembl_havana_gene"), "V9IDgene"]   %in% annodf$Ensembl_gene_id)
table( annodf$Ensembl_gene_id %in%  ens38df$V9IDgene )

### > table(table(annodf$Ensembl_gene_id))
### 
###     1     2     3     4     5     6     7     8     9    10    11    13    22 
### 12690  4133  1776   549   150    64    24    10     9     1     1     1     1 
### 18804 
### 1

## 18804 Illumina probes have no Ensembl gene id so we can't use annodf$Ensembl_gene_id to group into gene sets.

table(annodf$Ensembl_gene_id)[which(table(annodf$Ensembl_gene_id) == 18804)]
table(annodf$Ensembl_gene_id)[which(table(annodf$Ensembl_gene_id) == 22)]
annodf[which(annodf$Ensembl_gene_id == "ENSG00000157423"), -c(105:107)]
annodf[which(annodf$Ensembl_gene_id == "ENSG00000189064"), -c(105:107)]
head(annodf[which(annodf$Ensembl_gene_id == ""), -c(105:107)])
annodf[which(annodf$Ensembl_gene_id == ""), -c(105:107)][sample(seq(length(which(annodf$Ensembl_gene_id == ""))), size=10),]

hg19df <- read.table(file = paste("/Volumes/KilroyHD/kilroy/Projects/",
                                  "MOlab/TomoOsako/AgeAssociatedPaper2014/TCGAStudy/Data/",
                                  "HUGO_Biomart_GeneIDs_results.txt", sep = ""),
                     sep = "\t", header = TRUE, quote = "\"", comment.char = "",
                     fill = TRUE, stringsAsFactors = FALSE, colClasses = "character")
## Filter out withdrawn entries
dim(hg19df)
hg19df <- hg19df[hg19df$Status == "Approved", ]
dim(hg19df)
hg19df <- hg19df[!(hg19df$RefSeq.accession == "" & hg19df$Ensembl.gene.ID == ""), ]
dim(hg19df)


any(hg19df$Approved.symbol == "")  ## [1] FALSE
any(annodf$Gene_symbol == "")
head(which(annodf$Gene_symbol == ""))
annodf$hg19RefSeqGeneSymbol <- rep(NA_character_, nrow(annodf))
annodf$hg19reGeneSymbol <- rep(NA_character_, nrow(annodf))

rsmg <-
  cbind(
      annodf$Gene_symbol,
      hg19df$Approved.symbol[match(annodf$Ensembl_gene_id, hg19df$Ensembl.gene.ID)],
      hg19df$Approved.symbol[match(annodf$RefSeq_ID_0, hg19df$RefSeq.accession)],
      hg19df$Approved.symbol[match(annodf$RefSeq_transcripts, hg19df$RefSeq.accession)],
      hg19df$Approved.symbol[match(annodf$Original_transcriptomic_annotation, hg19df$RefSeq.accession)],
      hg19df$Approved.symbol[match(annodf$Lumi_transcriptomic_annotation, hg19df$RefSeq.accession)]
  )
is.na(rsmg[, 1]) <- annodf$Gene_symbol == ""
is.na(rsmg[, 2]) <- annodf$Ensembl_gene_id == ""
is.na(rsmg[, 3]) <- annodf$RefSeq_ID_0 == ""
is.na(rsmg[, 4]) <- annodf$RefSeq_transcripts == ""
is.na(rsmg[, 5]) <- annodf$Original_transcriptomic_annotation == ""
is.na(rsmg[, 6]) <- annodf$Lumi_transcriptomic_annotation == ""


annodf$hg19reGeneSymbol <-
  apply(rsmg, 1, function(x) paste(unique(x[!is.na(x)]), collapse = "|"))

annodf[sample(seq(nrow(annodf)), size = 10), c("Gene_symbol", "hg19RefSeqGeneSymbol", "hg19reGeneSymbol")]

### ------------------------------------------


tcgaandf <- read.csv(file = paste("/Volumes/KilroyHD/kilroy/Projects/MOlab/",
                                  "TomoOsako/AgeAssociatedPaper2014/TCGAStudy/",
                                  "Data/",
                                  "TCGA_BRCA_mRNAex_annotation.csv",
                                  sep = ""))
### dim(tcgaandf)
### dim(annodf)
### names(annodf)[which(grepl("ntrez", names(annodf)))]
### annodf[1:10, which(grepl("ntrez", names(annodf)))]
### 
### table(tcgaandf$Entrez_Gene_Id %in% annodf$Entrez)
### table(tcgaandf$Entrez_Gene_Id %in% annodf$Original_Entrez)
### table(tcgaandf$Entrez_Gene_Id %in% annodf$Entrez_Gene_ID_0)
### 
### table(annodf$Entrez %in%   tcgaandf$Entrez_Gene_Id)
### table(annodf$Original_Entrez  %in%  tcgaandf$Entrez_Gene_Id )
### table(annodf$Entrez_Gene_ID_0 %in%   tcgaandf$Entrez_Gene_Id  )

all.equal(sort(which(annodf$Original_Entrez  %in%  tcgaandf$Entrez_Gene_Id )),
          sort(which(annodf$Entrez_Gene_ID_0 %in%   tcgaandf$Entrez_Gene_Id))) ## [1] TRUE

### Need an Entrez_Gene_ID_matchp to match 18950 Illumina probes to RNAseq TCGA data.
### Match the probe with matching Entrez gene ID and the largest average value and the strongest trend if present.
### 3 loops.  1. match on entrez   2. For multiple entrez, match on strongest age associated.
###           3. For multiple entrez, match on largest average intensity.
exmedians <- apply(Dataset.r, 2, function(x) median(x, na.rm=TRUE))
### table(names(exmedians) %in% annodf$Probe_id)
### names(exmedians)[match(annodf$Probe_id, names(exmedians))][1:5]
### annodf$Probe_id[1:5]
annodf$raw_expression_medians <- exmedians[match(annodf$Probe_id, names(exmedians))]
aroutdf$raw_expression_medians <- exmedians[match(aroutdf$Probe_id, names(exmedians))]
### annodf$raw_expression_medians[1:5]
### median(Dataset.r[, "ILMN_1725881"])
### dimnames(Dataset.r)[[2]][1]

tcgaandf$MatchingIlluminaProbeId <- rep(NA_character_, nrow(tcgaandf))
aroutdf$MatchingRNASeqEntrezGeneIDp <- rep(FALSE, nrow(aroutdf))

for ( ei in seq(along = tcgaandf$Entrez_Gene_Id) ) {
    if (ei %% 100 == 0 ) { cat(ei, ", ") }
    
    entrezi <- tcgaandf$Entrez_Gene_Id[ei]
    entreziidxs <- which(aroutdf$Entrez_Gene_ID_0 %in% entrezi)
    
    if ( length(entreziidxs) == 1 ) {
        tcgaandf$MatchingIlluminaProbeId[ei] <- aroutdf$Probe_id[entreziidxs]
        aroutdf$MatchingRNASeqEntrezGeneIDp[entreziidxs] <- TRUE
    } else {
        if ( length(entreziidxs) > 1 ) {
            
            maxabsslopeidx <- which.max(abs(aroutdf$Best_log2FC[entreziidxs]))
            bestpvalidx <- which.min(abs(aroutdf$Best_BHadj_pval[entreziidxs]))
            maxexpressionidx <- which.max(aroutdf$raw_expression_medians[entreziidxs])
            
            if ( ((maxabsslopeidx == bestpvalidx) && 
                  (aroutdf$raw_expression_medians[entreziidxs][bestpvalidx] > 6.0)) ) {
                
                tcgaandf$MatchingIlluminaProbeId[ei] <- aroutdf$Probe_id[entreziidxs][maxabsslopeidx]
                aroutdf$MatchingRNASeqEntrezGeneIDp[entreziidxs][maxabsslopeidx] <- TRUE
            } else {
                if (aroutdf$raw_expression_medians[entreziidxs][bestpvalidx] > 6.0) {
                    
                    tcgaandf$MatchingIlluminaProbeId[ei] <- aroutdf$Probe_id[entreziidxs][bestpvalidx]
                    aroutdf$MatchingRNASeqEntrezGeneIDp[entreziidxs][bestpvalidx] <- TRUE
                } else {
                    
                    tcgaandf$MatchingIlluminaProbeId[ei] <- aroutdf$Probe_id[entreziidxs][maxexpressionidx]
                    aroutdf$MatchingRNASeqEntrezGeneIDp[entreziidxs][maxexpressionidx] <- TRUE
                }
            }
            
        } else {
            cat(paste("\nNo match for ei =", ei, ": Entrez_Gene_Id =", tcgaandf$Entrez_Gene_Id[ei], "\n"))
        }
    }
}

write.csv(tcgaandf,
          file = paste("/Volumes/KilroyHD/kilroy/Projects/MOlab/",
                       "TomoOsako/AgeAssociatedPaper2014/TCGAStudy/",
                       "Data/",
                       "TCGA_BRCA_mRNAex_annotation.csv",
                       sep = ""))
