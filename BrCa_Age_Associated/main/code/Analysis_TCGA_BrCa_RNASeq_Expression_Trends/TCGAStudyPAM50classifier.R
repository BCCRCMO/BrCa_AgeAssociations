### TCGA BrCa study set


require("genefu")

if( ! file.exists("Plots") ) { dir.create("Plots") }
if( ! file.exists("Plots/SinglePlots") ) { dir.create("Plots/SinglePlots") }


## a useful function: rev() for strings
strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

leftPad0 <- function(x,
                     length = max(nchar(x)))
{
  ### Note:  for numeric data, formatC(numbervec, width = 3, flag = "0")
  ###        will left-pad numerics with zeroes.
  if ( is.numeric(x) ) {
    return( formatC(x, width = length, flag = "0") )
  } else {
    zeroes <- paste(rep("0", length), collapse = "")
    return( strReverse(substring(strReverse(paste(zeroes, x, sep = "")), 1, length)) )
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

annodf <- read.delim(file = paste("../Data/",
                       "Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09.txt",
                       sep = ""),
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

icdf$tcgaid <- paste(substring(gsub("\\-", "\\.", icdf$ID), 1, 12), "01", sep = ".")

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

cldf$tcgaid <- paste(gsub("-", ".", cldf$PATIENT_ID), "01", sep = ".")
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


all.equal(seq(1079)+2, match(sidf$tcgaid, names(edf))) ## TRUE  ## 1079 cases plus hugo and entrez cols
all.equal(seq(1079)+11, match(sidf$tcgaid, names(exandf))) ## TRUE  ## 1079 cases plus 11 annotation cols

### ### Match sedf to annotation data in hgandf

sedf <- cbind(sidf, log2(16.0 + t( as.matrix(exandf[, match(sidf$tcgaid, names(exandf))]) ) ) )
names(sedf) <- c(names(sidf), exandf$Hugo_Symbol)
dim(sedf) ## [1]  1079 18976 ## Down from 1091 with 12 men out
length(c(names(sidf), exandf$Hugo_Symbol))


### PAM50

### data  Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.

exan50df <- exandf[which(exandf$Entrez_Gene_Id %in% pam50$centroids.map$EntrezGene.ID), ]
emat <- log2(as.matrix(t(exan50df[, 12:1090])))
dimnames(emat)[[2]] <- exan50df$Gene_symbol
ermat <- apply(emat, 2, rescale)
amat <- exan50df[, c("Entrez_Gene_Id", "Gene_symbol", "Gene_symbol")]
names(amat) <- c("EntrezGene.ID", "Gene.Symbol", "probe")
PAM50Preds <- molecular.subtyping(sbt.model = "pam50", data = emat, annot = amat, do.mapping = FALSE)
ucscdf <- read.csv(file = "../Data/UCSC_TCGA_PAM50.csv")
ucscdf$tcgaid <- gsub("-", ".", ucscdf$sampleid)
table(sedf$tcgaid %in% ucscdf$tcgaid)

table(ucscdf$PAM50[match(sedf$tcgaid, ucscdf$tcgaid)])

all.equal(sedf$tcgaid, names(PAM50Preds$subtype)) ## TRUE
sedf$PAM50 <- PAM50Preds$subtype
sedf$survstat <- 1.0 * (sedf$OS_STATUS == "DECEASED")
sedf$survtime <- as.numeric(sedf$OS_MONTHS)
### fit a Kaplan-Meier and plot it
### One neg surv time, several 0 survtimes
is.na(sedf$survtime) <- sedf$survtime <= 0
fit <- survfit(Surv(survtime, survstat) ~ PAM50, data = sedf)
survreg(Surv(survtime, survstat) ~ PAM50, data = sedf)
survdiff(Surv(survtime, survstat) ~ PAM50, data = sedf)
pdf(file = "./Plots/TCGA_PAM50_genefu_initial.pdf", width = 7, height = 7, useDingbats = FALSE)
plot(fit, col =c("#006d2c", "#8856a7","#a50f15", "#08519c", "#000000"),lty = 1,lwd = 3,
     xlab = "Time (months)",ylab = "Probability of Survival")
legend("topright",
       fill = c("#006d2c", "#8856a7","#a50f15", "#08519c", "#000000"),
       legend = levels(sedf$PAM50))
dev.off()
