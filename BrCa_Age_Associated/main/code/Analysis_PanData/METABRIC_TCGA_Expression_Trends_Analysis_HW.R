## Version to merge METABRIC and TCGA and Tier 1/2
## Then can get sets of overlapping genes and do volcano plots
## Merge TCGA Kidney, Lung, Prostate, Thyroid as well so can look at common genes between
##  breast and others.

## Version to address genes of interest to Hakwoo Lee
## Added legend position for volcano plots

# Setup -------------------------------------------------------------------

# Packages
library(tidyverse)
library(magrittr)
library(glue)
library(survival)
library(descr)
library(fs)
library(ageassn)

##

arlmdf <- readRDS(file = paste("/Users/stevenmckinney/git/github/",
                  "BrCa_Age_Associated2/main/code/",
                  "Analysis_METABRIC_Expression_Trends/Memoized/",
                  "AgeRelated_AllProbes_All_cases_FDR_0p01_All1992Cases_lm_v01.rds", sep = ""))
arlmdf <- data.frame(arlmdf)

tlmdf <- read.csv(file = paste("/Users/stevenmckinney/git/github/",
                               "BrCa_Age_Associated2/main/code/",
                               "Analysis_TCGA_BrCa_RNASeq_Expression_Trends/Memoized/",
                               "TCGA_BrCa_All_cases_FDR_0p01_outdf.csv", sep = "")
                  )
tadf <- read.csv(file = paste("/Users/stevenmckinney/git/github/",
                              "BrCa_Age_Associated2/main/data/",
                              "Data_TCGA_BrCa_RNASeq_Expression_Trends/",
                              "TCGA_BRCA_mRNAex_annotation.csv", sep = ""))
talmdf <- merge(tadf, tlmdf, by = "Entrez_Gene_Id", suffixes = c("", ".tcga"))

ebdf_t1 <- as.data.frame(read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_1_EntrezID_and_Symbol_ER_at_promotors.txt"))
ebdf_t2 <- as.data.frame(read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_2_EntrezID_and_Symbol_total_ER_binding_via_20kb.txt"))

names(ebdf_t1) <- c("Entrez_Gene_Id", "Gene_symbol")
names(ebdf_t2) <- c("Entrez_Gene_Id", "Gene_symbol")
ebdf_t12 <- talmdf[, c("Entrez_Gene_Id", "Gene_symbol")]
ebdf_t12$ERbt12 <- rep("Non binding", nrow(talmdf))
ebdf_t12$ERbt12[which(ebdf_t12$Entrez_Gene_Id %in% ebdf_t2$Entrez_Gene_Id)] <- "Tier 2 Only"
ebdf_t12$ERbt12[which(ebdf_t12$Entrez_Gene_Id %in% ebdf_t1$Entrez_Gene_Id)] <- "Tier 1"
# Entrez_Gene_Id  Gene_symbol

tbalmdf <- merge(talmdf, ebdf_t12, by = "Entrez_Gene_Id", suffixes = c("", ".nki"))

## METABRIC and TCGA breast data merged
mtdf <- merge(arlmdf, tbalmdf, 
              by.x = "ProbeId", 
              by.y = "MatchingIlluminaProbeId", 
              suffixes = c(".MB", ".TCGA_Br"))

## Merge also Kidney, Lung Prostate, Thyroid
## Kidney
tkadf <- read.csv(file = paste("/Users/stevenmckinney/git/github/",
                               "BrCa_Age_Associated2/main/code/",
                               "Analysis_TCGA_KidneyCa_RNASeq_Expression_Trends/AllCases/tables/",
                               "aroutdf_FDR_0p01_FC_1p25.csv", sep = ""))
tkadf$Abs_FoldChange <- 2^abs(tkadf$Best_log2FC)
tkadf$FoldChange_Direction <- ifelse(tkadf$Best_log2FC > 0, "Up", "Down")

mtkdf <- merge(mtdf, tkadf, by="Entrez_Gene_Id", suffixes = c("", ".TCGA_Kidney"))

## Lung
tladf <- read.csv(file = paste("/Users/stevenmckinney/git/github/",
                               "BrCa_Age_Associated2/main/code/",
                               "Analysis_TCGA_LungCa_RNASeq_Expression_Trends/AllCases/tables/",
                               "aroutdf_FDR_0p01_FC_1p25.csv", sep = ""))
tladf$Abs_FoldChange <- 2^abs(tladf$Best_log2FC)
tladf$FoldChange_Direction <- ifelse(tladf$Best_log2FC > 0, "Up", "Down")

mtkldf <- merge(mtkdf, tladf, by="Entrez_Gene_Id", suffixes = c("", ".TCGA_Lung"))

## Prostate
tpadf <- read.csv(file = paste("/Users/stevenmckinney/git/github/",
                               "BrCa_Age_Associated2/main/code/",
                               "Analysis_TCGA_ProstateCa_RNASeq_Expression_Trends/tables/",
                               "aroutdf_FDR_0p01_FC_1p25.csv", sep = ""))
tpadf$Abs_FoldChange <- 2^abs(tpadf$Best_log2FC)
tpadf$FoldChange_Direction <- ifelse(tpadf$Best_log2FC > 0, "Up", "Down")

mtklpdf <- merge(mtkldf, tpadf, by="Entrez_Gene_Id", suffixes = c("", ".TCGA_Prostate"))


## Thyroid
ttadf <- read.csv(file = paste("/Users/stevenmckinney/git/github/",
                               "BrCa_Age_Associated2/main/code/",
                               "Analysis_TCGA_ThyroidCa_RNASeq_Expression_Trends/AllCases/tables/",
                               "aroutdf_FDR_0p01_FC_1p25.csv", sep = ""))
ttadf$Abs_FoldChange <- 2^abs(ttadf$Best_log2FC)
ttadf$FoldChange_Direction <- ifelse(ttadf$Best_log2FC > 0, "Up", "Down")

mtklptdf <- merge(mtklpdf, ttadf, by="Entrez_Gene_Id", suffixes = c("", ".TCGA_Thyroid"))

## UpSetR Venn Sets
venndf <- read.table(file = paste("./main/code/Analysis_PanData/",
                                            "UpSetR_MB_TCGA_Gene_VennSet.csv", 
                                            sep = ""),
                      sep = ",", header = TRUE, quote = '"',
                     stringsAsFactors = FALSE
)
#head(venndf)
#table(venndf$Hugo_Symbol %in% mtklptdf$Hugo_Symbol, useNA = "always")
mtklptvdf <- merge(mtklptdf, venndf, by = "Hugo_Symbol", all = TRUE)
mtklptvdf$AgeAssociatedVennSet[which(is.na(mtklptvdf$AgeAssociatedVennSet))] <- "Not age associated"
#dim(mtklptvdf)
#head(mtklptvdf, 1)
#mtklptvdf[which(mtklptvdf$Hugo_Symbol == "TRIM5"), ]
#table(mtklptvdf$AgeAssociatedVennSet, useNA = "always")
#names(mtklptvdf)
mtklptvdf$X1 <- mtklptvdf$X <- mtklptvdf$X.1 <- mtklptvdf$X.tcga <- NULL

mtklptvdf$ERbinding.TCGA_Thyroid <- mtklptvdf$ERbinding <- mtklptvdf$ERbinding.TCGA_Lung <- NULL
mtklptvdf$ERbinding.TCGA_Br <- mtklptvdf$ERbinding.TCGA_Prostate <- NULL
mtklptvdf$ERbinding_NKI.TCGA_Thyroid <- mtklptvdf$ERbinding_NKI.TCGA_Lung <- NULL
mtklptvdf$ERbinding_NKI.TCGA_Prostate <- mtklptvdf$ERbt12 <- NULL
mtklptvdf$arGeneSet_BHadj_and_AgeDependentp.TCGA_Br <- mtklptvdf$Hugo_Symbol.tcga <- NULL
mtklptvdf$Entrez_Gene_ID_0.TCGA_Br <- mtklptvdf$RefSeq_ID_0.TCGA_Br <- NULL
mtklptvdf$Accession_0.TCGA_Br <- mtklptvdf$Gene_symbol.TCGA_Br <- mtklptvdf$ProbeId <- NULL
mtklptvdf$AgeDependentp.MB <- mtklptvdf$arGeneSet_LE60_BHadj_pval <- mtklptvdf$arGeneSet_GT60_BHadj_pval <- NULL
mtklptvdf$arGeneSet_AllAges_BHadj_pval <- mtklptvdf$arGeneSet_BHadj_and_AgeDependentp.MB <- NULL
mtklptvdf$arGeneSet_PBHadj  <- mtklptvdf$arGeneSetidxp <- mtklptvdf$TCGAsymbolEQMETABRICsymbol <- NULL
mtklptvdf$AgeDependentp.TCGA_Br <- mtklptvdf$BHadj_signifp.TCGA_Br <- NULL
mtklptvdf$arGeneSet_BHadj_and_AgeDependentp.TCGA_Br <- mtklptvdf$BHadj_signifp <- NULL
mtklptvdf$BHadj_signifp.TCGA_Prostate <- NULL
mtklptvdf$ERbinding.MB  <- mtklptvdf$BHadj_signifp.MB <- mtklptvdf$SNP.tcga  <- mtklptvdf$CHR.tcga <- NULL
mtklptvdf$BP.tcga <- mtklptvdf$refseqID.tcga <- mtklptvdf$BHadj_signifp.TCGA_Lung <- NULL
mtklptvdf$AgeDependentp <- mtklptvdf$AgeDependentp.TCGA_Lung <- NULL
mtklptvdf$P <- mtklptvdf$PBHadj <- NULL
mtklptvdf$AgeDependentp.TCGA_Thyroid  <- mtklptvdf$AgeDependentp.TCGA_Prostate <- NULL
mtklptvdf$Hugo_Symbol.TCGA_Thyroid <- NULL
mtklptvdf$Hugo_Symbol.TCGA_Prostate  <- mtklptvdf$Hugo_Symbol.TCGA_Lung  <- NULL
mtklptvdf$Gene_symbol.nki	 <- mtklptvdf$Hugo_Symbol.TCGA_Kidney <- mtklptvdf$Gene_symbol.tcga <- NULL
mtklptvdf$BHadj_signifp.TCGA_Thyroid <- NULL

mtklptvdf$Probe_sequence <- mtklptvdf$Probe_length <- mtklptvdf$Probe_type <- NULL
mtklptvdf$Lumi_id <- NULL
mtklptvdf$Genomic_location <- mtklptvdf$Original_genomic_annotation <- NULL
mtklptvdf$Original_genomic_match <- NULL
mtklptvdf$UCSC_genomic_match <- mtklptvdf$Repeat_Mask <- NULL
mtklptvdf$Nr_nuc_RMasked <- mtklptvdf$CpG_islands	<- mtklptvdf$SNPs <- NULL
mtklptvdf$Proportion_RefSeq_transcripts <- mtklptvdf$UCSC_transcripts <- NULL
mtklptvdf$Proportion_UCSC_transcripts <- mtklptvdf$GenBank_transcripts <- NULL
mtklptvdf$Proportion_GenBank_transcripts <- mtklptvdf$Exons <- NULL
mtklptvdf$Proportion_Ensembl_transcripts <- mtklptvdf$Original_transcriptomic_annotation <- NULL
mtklptvdf$Original_Unigene <- mtklptvdf$Original_Unigene_match <- NULL
mtklptvdf$Original_Entrez <- mtklptvdf$Original_Entrez_match <- NULL
mtklptvdf$Original_genic_annotation <- mtklptvdf$Original_genic_match <- NULL 
mtklptvdf$Other_gene_ids <- mtklptvdf$X2nd_match_gene_ids <- mtklptvdf$Unique_transcriptomic <- NULL
mtklptvdf$Coding_zone <- mtklptvdf$Perfect_match <- mtklptvdf$GC_contents_. <- mtklptvdf$SNP_number <- NULL
mtklptvdf$Quality_score <- mtklptvdf$Quality_comments <- NULL
mtklptvdf$Species_0 <- mtklptvdf$Source_0	 <- mtklptvdf$Search_Key_0 <- NULL
mtklptvdf$Transcript_0 <- mtklptvdf$ILMN_Gene_0 <- NULL
mtklptvdf$Source_Reference_ID_0	 <- mtklptvdf$RefSeq_ID_0.MB <- NULL
mtklptvdf$Unigene_ID_0 <- mtklptvdf$Entrez_Gene_ID_0.MB <- mtklptvdf$GI_0	 <- NULL
mtklptvdf$Accession_0.MB <- mtklptvdf$Symbol_0 <- mtklptvdf$Protein_Product_0 <- NULL
mtklptvdf$Probe_Id_0 <- mtklptvdf$Array_Address_Id_0 <- NULL
mtklptvdf$Probe_Type_0 <- mtklptvdf$Probe_Start_0 <- mtklptvdf$Probe_Sequence_0 <- NULL
mtklptvdf$Chromosome_0 <- mtklptvdf$Probe_Chr_Orientation_0 <- NULL
mtklptvdf$Probe_Coordinates_0 <- mtklptvdf$Cytoband_0 <- mtklptvdf$Definition_0 <- NULL
mtklptvdf$Synonyms_0 <- mtklptvdf$Obsolete_Probe_Id_0 <- NULL
mtklptvdf$SNP <- mtklptvdf$CHR <- mtklptvdf$BP <- mtklptvdf$refseqID <- NULL
mtklptvdf$Search_key <- mtklptvdf$Gene_symbol.MB <- NULL
mtklptvdf$miRNAs <- mtklptvdf$Other_genomic_matches <- mtklptvdf$Genomic_match_scheme <- NULL
mtklptvdf$Genomic_match_length <- mtklptvdf$Genomic_match_._similarity <- NULL
mtklptvdf$Genomic_match_overall_._similarity <- mtklptvdf$Genomic_2nd_matches <- NULL
mtklptvdf$Genomic_2nd_match_scheme <- mtklptvdf$Genomic_2nd_match_length <- NULL
mtklptvdf$Genomic_2nd_match_._similarity <- mtklptvdf$Genomic_2nd_match_overall_._similarity <- NULL
mtklptvdf$RefSeq_transcripts <- mtklptvdf$Ensembl_transcripts <-NULL
mtklptvdf$Original_transcriptomic_match <- NULL
mtklptvdf$GIs <- mtklptvdf$Original_GI <- mtklptvdf$Original_GI_match <- NULL
mtklptvdf$CCDSs <- mtklptvdf$Proteins <- NULL	
mtklptvdf$Original_protein <- mtklptvdf$Original_protein_match	 <- mtklptvdf$SwissProts <- NULL
mtklptvdf$Ensembl_transcriptomic_annotation <- mtklptvdf$Ensembl_transcriptomic_match <- NULL
mtklptvdf$Lumi_transcriptomic_annotation <- mtklptvdf$Lumi_transcriptomic_match	<- NULL
mtklptvdf$Splice_junction <- NULL
mtklptvdf$Transcriptomic_match_scheme <- mtklptvdf$Transcriptomic_match_length <- NULL
mtklptvdf$Transcriptomic_match_._similarity <- mtklptvdf$Transcriptomic_match_overall_._similarity <- NULL	
mtklptvdf$Transcriptomic_2nd_matches <- mtklptvdf$Transcriptomic_2nd_match_scheme <- NULL
mtklptvdf$Transcriptomic_2nd_match_length <- mtklptvdf$Transcriptomic_2nd_match_._similarity <- NULL	
mtklptvdf$Transcriptomic_2nd_match_overall_._similarity <- mtklptvdf$Gene_synonyms <- NULL




## Rename with .TCGA_Kidney  "Hugo_Symbol.TCGA_Kidney"                      
krvars <- c("LE60_slope", "GT60_slope", "AllAges_slope",                                 
            "LE60_pval", "GT60_pval", "AllAges_pval",                                  
            "LE60_BHadj_pval", "GT60_BHadj_pval",                          
            "AllAges_BHadj_pval", "BHadj_and_AgeDependentp", "LE60_log2FC",                                  
            "GT60_log2FC", "AllAges_log2FC", "Best_log2FC",                                  
            "Best_BHadj_pval", "Abs_FoldChange", "FoldChange_Direction")

if( !any(is.na(match(krvars, names(mtklptvdf)))) ) {
  names(mtklptvdf)[match(krvars, names(mtklptvdf))] <- paste(krvars, ".TCGA_Kidney", sep = "")
}

# move ERbinding_NKI to end
mtklptvdf$ERbinding_NKI_end <- mtklptvdf$ERbinding_NKI
mtklptvdf$ERbinding_NKI <- NULL
mtklptvdf$ERbinding_NKI <- mtklptvdf$ERbinding_NKI_end
mtklptvdf$ERbinding_NKI_end <- NULL

head(mtklptvdf, 2)

write.table(mtklptvdf, 
          file = paste("/Users/stevenmckinney/git/github/",
                       "BrCa_Age_Associated2/main/code/Analysis_PanData/",
                       "Supplemental_AgeAssociatedFoldChanges_AllCases_MB_TCGA.csv", 
                       sep = "" ),
          sep = ",", quote = TRUE, row.names = FALSE )

##

# Raw and intermediate data inputs
annodf <- read_rds("main/data/Data_METABRIC_Expression_Trends/Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09.rds")
ebdf <- read.delim("main/data/Data_JasonCarroll/ER_binding_gene_hg19_threshold_1e-5.txt")  # Raw file has row names, can't use read_tsv
ebdf_t1 <- read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_1_EntrezID_and_Symbol_ER_at_promotors.txt")
ebdf_t2 <- read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_2_EntrezID_and_Symbol_total_ER_binding_via_20kb.txt")
#ehdf <- read_csv("main/data/Data_TomoOsako_BrCa_AgeRelated_EZH2Related/EZH2relatedGenesv08.csv")
### Age study unique df  asudf$MBid are the 1161 cases not in MB09 TMA.
asudf <- read_csv("main/data/Data_BrCaSubtypeDistributionBySeries/MBex_MBid_NoOverlap_MB09_BigSeries.csv")
#sehdf <- read_rds("main/data/Data_METABRIC_Expression_Trends/sehdf.rds")
sehrdf <- read_rds("main/data/Data_METABRIC_Expression_Trends/sehrdf.rds")
ehzidxs <- read_rds("main/data/Data_METABRIC_Expression_Trends/ehzidxs.rds")
#sardf <- read_rds("main/data/Data_METABRIC_Expression_Trends/sardf.rds")
sarrdf <- read_rds("main/data/Data_METABRIC_Expression_Trends/sarrdf.rds")
# aridxs maps to rows of annotation df annodf
aridxs <- read_rds("main/data/Data_METABRIC_Expression_Trends/aridxs.rds")
# arzidxs maps to columns of expression matrix (probesets are columns)
arzidxs <- read_rds("main/data/Data_METABRIC_Expression_Trends/arzidxs.rds")
Dataset_r <- read_rds("main/data/Data_METABRIC_Expression_Trends/Dataset_r.rds")
## Make sehrdf, sarrdf, ehzidxs, aridxs, arzidxs for Hakwoo's list and ?biphasic?

# Need probeids, rows in annodf (aridxs), cols in Dataset_r (arzidxs)

gndf <- data.frame(GeneName = c("ELF5", "ELF5"),
                   IHCName = c("ELF5", "ELF5"),
                   refseqID = c("NM_001422","NM_198381"),
                   ProbeIDs = c("ILMN_1684699",  "ILMN_1813270"),
                   stringsAsFactors = FALSE)
hw.rdf <- 
  (Dataset_r[, which(dimnames(Dataset_r)[[2]]  %in% 
                       gndf$ProbeIDs)] %>% 
     data.frame()
  )
names(hw.rdf) <- paste(gndf$GeneName[1], names(hw.rdf), sep = "|")

## Add Hakwoo's genes into the set
if ( all.equal(row.names(hw.rdf), sarrdf$MBid) ) {
  sarrdf <- bind_cols(hw.rdf, sarrdf %>% select(., MBid:rev(names(sarrdf))[1])) } else {
    stop("Row order mismatch for sarrdf construction")
  }
aridxs <- match ( gndf$ProbeIDs , annodf$Probe_id)
arzidxs <- match ( gndf$ProbeIDs , dimnames(Dataset_r)[[2]])

# aroutdf objects
root <- "main/code/Analysis_METABRIC_Expression_Trends/aroutdf"

# All Cases
aroutdf_all <- 
  dir_ls(path(root, "AllCases")) %>% 
  purrr::set_names(gsub(".*FC(.*)_.*", "FC_\\1", .)) %>%
  map(read_csv)

# IntClust
aroutdf_intclust <- 
  dir_ls(path(root, "intClust")) %>%
  purrr::set_names(gsub(".*FC(.*)_.*", "FC_\\1", .)) %>%
  map(read_csv)

# PAM50
aroutdf_pam50 <- 
  dir_ls(path(root, "Pam50")) %>% 
  purrr::set_names(gsub(".*FC(.*)_.*", "FC_\\1", .)) %>%
  map(read_csv)

# ER/HER2
aroutdf_erher2 <-
  dir_ls(path(root, "ER_HER2")) %>% 
  purrr::set_names(gsub(".*FC(.*)_.*", "FC_\\1", .)) %>%
  map(read_csv)

# ER+/ER-
aroutdf_er <- 
  dir_ls(path(root, "ER")) %>% 
  purrr::set_names(gsub(".*FC(.*)_.*", "FC_\\1", .)) %>%
  map(read_csv)


# Prep --------------------------------------------------------------------

# Add variables to ebdf, annodf, and sarrdf
ebdf <- ebdf %>% 
  mutate(I_Ensembl_gene_id = match(ensembl_gene_id, annodf$Ensembl_gene_id) %>% 
           ifelse(is.na(ensembl_gene_id), NA, .),
         I_Ensembl_gene_id_Gene_symbol = annodf$Gene_symbol[I_Ensembl_gene_id],
         I_Entrez = match(entrezgene, annodf$Entrez) %>% 
           ifelse(is.na(entrezgene), NA, .),
         I_Entrez_Gene_symbol = annodf$Gene_symbol[I_Entrez],
         I_Original_Entrez = match(entrezgene, annodf$Original_Entrez) %>% 
           ifelse(is.na(entrezgene), NA, .),
         I_Original_Entrez_Gene_symbol = annodf$Gene_symbol[I_Original_Entrez],
         I_Entrez_Gene_ID_0 = match(entrezgene, annodf$Entrez_Gene_ID_0) %>% 
           ifelse(is.na(entrezgene), NA, .),
         I_Entrez_Gene_ID_0_Gene_symbol = annodf$Gene_symbol[I_Entrez_Gene_ID_0],
         NumMatch = strsplit(paste(I_Ensembl_gene_id, I_Entrez,
                                   I_Original_Entrez, I_Entrez_Gene_ID_0), " ") %>% 
           vapply(function(x) {
             xn <- x[x != "NA"]
             return(ifelse(length(xn) == 0, 0, length(unique(xn))))
           }, numeric(1)),
         NumGenesymbolMatch = strsplit(paste(
           I_Ensembl_gene_id_Gene_symbol, I_Entrez_Gene_symbol,
           I_Original_Entrez_Gene_symbol, I_Entrez_Gene_ID_0_Gene_symbol), " ") %>% 
           vapply(function(x) {
             xn <- x[x != "NA"]
             return(ifelse(length(xn) == 0, 0, length(unique(xn))))
           }, numeric(1)),
         Genesymbols = strsplit(paste(
           I_Ensembl_gene_id_Gene_symbol, I_Entrez_Gene_symbol,
           I_Original_Entrez_Gene_symbol, I_Entrez_Gene_ID_0_Gene_symbol), " ") %>% 
           vapply(function(x) {
             xn <- x[x != "NA"]
             return(ifelse(length(xn) == 0, "", paste(unique(x[!(x %in% c("NA", ""))]),
                                                      collapse = " | ")))
           }, character(1))
  ) %>% 
  set_rownames(as.character(row.names(.)))

annodf <- annodf %>% 
  mutate(ProbeId = Probe_id,
         ERbinding = ifelse(Ensembl_gene_id %in% ebdf$ensembl_gene_id[!is.na(ebdf$ensembl_gene_id)] |
                              Entrez %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)] |
                              Original_Entrez %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)] |
                              Entrez_Gene_ID_0 %in% ebdf$entrezgene[!is.na(ebdf$entrezgene)], TRUE, FALSE))

sarrdf <- sarrdf %>% 
  mutate(
    BrCaEHf = case_when(
      ER.Expr == "+" & Her2.Expr == "+" ~ "ER+/HER2+",
      ER.Expr == "+" & Her2.Expr == "-" ~ "ER+/HER2-",
      ER.Expr == "-" & Her2.Expr == "+" ~ "ER-/HER2+",
      ER.Expr == "-" & Her2.Expr == "-" ~ "ER-/HER2-",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-")),
    BrCaEf = case_when(
      ER.Expr == "+" ~ "ER+",
      ER.Expr == "-" ~ "ER-",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("ER+", "ER-"))
  )

# BrCa Datasets
brcaIclalldf <- read_tsv("main/data/Data_METABRIC_Expression_Trends/clinical_datasetI_public_FINAL_032112.txt")
brcaIIclalldf <- read_tsv("main/data/Data_METABRIC_Expression_Trends/clinical_datasetII_public_FINAL_032112.txt")
clinalldf <- read_tsv("main/data/Data_METABRIC_Expression_Trends/clinic_all_2010-09-01.txt")

clvars <- intersect(names(brcaIclalldf), names(brcaIIclalldf))
cldf <- rbind(brcaIclalldf[, clvars], brcaIIclalldf[, clvars]) %>% 
  mutate(MBid = setMBid(METABRIC_ID))
clinalldf <- clinalldf %>% 
  mutate(MBid = setMBid(metabric_cl_id))

mvarsclinall <- c("MBid", "age_at_diagnosis", "dob", "date_of_diagnosis", "dod", "last_follow_up_date",
                  "last_follow_up_status", "first_relapse_type",
                  "time_until_first_local_regional_relapse", "time_until_first_distant_relapse",
                  "histological_type", "distant_metastasis_sites_byCarlos")

cldf <- merge(cldf, clinalldf[, mvarsclinall], by = "MBid", suffixes = c("", ".y"))

# Cluster analysis identification
clusterIcldf <- read_tsv("main/data/Data_METABRIC_Expression_Trends/cambridge_clusters_ordered.txt",
                         col_names = FALSE)
clusterIcldf <- clusterIcldf %>% 
  set_names(c("mb_id", "iClusterGroup")) %>% 
  mutate(MBid = setMBid(mb_id),
         Ctr = MBid2Ctr(MBid))

# intClust calls for validation set
clusterIIcldf <- read_tsv("main/data/Data_METABRIC_Expression_Trends/iClusterDatasetIIPredictions.txt")
clusterIIcldf <- clusterIIcldf %>% 
  set_names(c("mb_id", "iClusterGroup")) %>% 
  mutate(MBid = setMBid(mb_id),
         Ctr = MBid2Ctr(MBid))

# Combine cluster dfs
clustercldf <- rbind(clusterIcldf, clusterIIcldf) %>% 
  mutate(iClustf = factor(iClusterGroup, levels = as.character(1:10)))

# PAM50 BrCa Classification
brcacldf <- read_tsv("main/data/Data_METABRIC_Expression_Trends/DatasetI997_PAM50classification_ordered.txt",
                     col_names = FALSE)
brcacldf <- brcacldf %>% 
  set_names(c("mb_id", "BrCaSubtype")) %>% 
  mutate(mbidname = gsub("-", ".", mb_id),
         BrCaf = factor(BrCaSubtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal")))

# Merge clustercldf and brcacldf
mbcldf <- merge(clustercldf, brcacldf) %>% 
  mutate(isDiscovery = MBid %in% clusterIcldf$MBid)

# Merge BrCa subtype groupings with clinical data
cldvdf <- merge(clustercldf, cldf, by = "MBid", suffixes = c("", ".y")) %>% 
  mutate(isDiscovery = MBid %in% clusterIcldf$MBid,
         BrCaf = factor(Pam50Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal")),
         survstat = 1L * !(last_follow_up_status == "a"),
         survyrs = `T` / 365.25,
         gradeb = ifelse(grade == "null", NA, 1L * (grade == "3")),
         gradebf = factor(gradeb, levels = 0:1, labels = c("grade 1+2", "grade 3")))

## No overlap with MB09 TMA or Big Series
#sarnovdf <- sardf[sardf$MBid %in% asudf$MBid, ]
sarrnovdf <- sarrdf[sarrdf$MBid %in% asudf$MBid, ]

# Calculate statistical significance and biological relevance across the ar gene set (467 probes).
aroutdf_all <- aroutdf_all %>% 
  map(metabric_ar, aridxs = aridxs, annodf = annodf, fdr = 0.01)
aroutdf_intclust <- aroutdf_intclust %>% 
  map(metabric_ar, aridxs = aridxs, annodf = annodf, fdr = 0.01)
aroutdf_pam50 <- aroutdf_pam50 %>% 
  map(metabric_ar, aridxs = aridxs, annodf = annodf, fdr = 0.01)
aroutdf_erher2 <- aroutdf_erher2 %>% 
  map(metabric_ar, aridxs = aridxs, annodf = annodf, fdr = 0.01)
aroutdf_er <- aroutdf_er %>% 
  map(metabric_ar, aridxs = aridxs, annodf = annodf, fdr = 0.01)

# Tidy gene lists
erbnms <- tidy_genes(ebdf, "hgnc_symbol")
erbnms_t1 <- tidy_genes(ebdf_t1, "Symbol")
erbnms_t2 <- tidy_genes(ebdf_t2, "Symbol")


# Survival ----------------------------------------------------------------

# Kaplan-Meier plots
PAM50plotcols <- c("purple", "cyan", "orange", "red", "green")
plot(survfit(Surv(survyrs, survstat) ~ BrCaf, data = cldvdf), col = PAM50plotcols)
legend("bottomleft", legend = levels(cldvdf$BrCaf), col = PAM50plotcols, lwd = 3)

iClustplotcols <- c("#FF5500", "#00EE76", "#CD3278", "#00C5CD", "#8B0000",
                    "#FFFF40", "#0000CD", "#FFAA00", "#EE82EE", "#7D26CD")
plot(survfit(Surv(survyrs, survstat) ~ iClustf, data = cldvdf), col = iClustplotcols)
legend("bottomleft", legend = levels(cldvdf$iClustf), col = iClustplotcols, lwd = 3)

# Higher weight to earlier diffs - better for this data
# because of crossing of survival curves at and beyond 7 years.
# "Breslow" test  rho = 1:  Breslow test is close to G-rho test with rho = 1.
survdiff(Surv(survyrs, survstat) ~ BrCaf, data = cldvdf, rho = 1)
survdiff(Surv(survyrs, survstat) ~ BrCaf, data = cldvdf, rho = 0)
survdiff(Surv(survyrs, survstat) ~ iClustf, data = cldvdf, rho = 1)
survdiff(Surv(survyrs, survstat) ~ iClustf, data = cldvdf, rho = 0)

# Cox models
cmf <- coxph(Surv(survyrs, survstat) ~ gradebf + iClustf, data = cldvdf,
             subset = !(is.na(gradebf) | is.na(iClustf) | is.na(`T`) | is.na(survstat)))
cmr <- coxph(Surv(survyrs, survstat) ~ iClustf, data = cldvdf,
             subset = !(is.na(gradebf) | is.na(iClustf) | is.na(`T`) | is.na(survstat)))
anova(cmr, cmf)

# Stratified by iClustf
cmsf <- coxph(Surv(survyrs, survstat) ~ gradebf + strata(iClustf), data = cldvdf,
              subset = !(is.na(gradebf) | is.na(iClustf) | is.na(`T`) | is.na(survstat)))
cmsr <- coxph(Surv(survyrs, survstat) ~ strata(iClustf), data = cldvdf,
              subset = !(is.na(gradebf) | is.na(iClustf) | is.na(`T`) | is.na(survstat)))

anova(cmsf)


# Plots -------------------------------------------------------------------

# Set plot path and common arguments
plot_path <- "main/code/Analysis_METABRIC_Expression_Trends/Plots/"
pargs <- list(
  x = "age_at_diagnosis",
  age = 60,
  conf.level = 0.99,
  show.avg = TRUE
)
# 
# responses <- sort(names(sehdf)[seq_along(ehzidxs)])
# pdf(file = glue("{plot_path}EZH2_All_AgeRelated_lm_v01_DC.pdf"),
#     width = 8, height = 10)
# par(mfrow = c(2, 2))
# responses %>%
#   walk(~ invoke(
#     regression_plot,
#     pargs,
#     y = .,
#     data = sehdf,
#     main = "All cases:"
#   ))
# dev.off()
# 
# responses <- sort(names(sehrdf)[seq_along(ehzidxs)])
# pdf(file = glue("{plot_path}EZH2_All_AgeRelated_lm_raw_v01_DC.pdf"),
#     width = 8, height = 10)
# par(mfrow = c(2, 2))
# responses %>%
#   walk(~ invoke(
#     regression_plot,
#     pargs,
#     y = .,
#     data = sehrdf,
#     main = "All cases:"
#   ))
# dev.off()

# responses <- sort(names(sardf)[seq_along(arzidxs)])
# pdf(file = glue("{plot_path}AgeRelated_All1992Cases_lm_v10_DC.pdf"),
#     width = 8, height = 10)
# par(mfrow = c(2, 2))
# responses %>%
#   walk(~ invoke(
#     regression_plot,
#     pargs,
#     y = .,
#     data = sardf,
#     main = paste(nrow(sardf), "cases:"),
#     annotate = aroutdf_all$FC_1.25_AllCases$Probe_id
#   ))
# dev.off()
# 
# responses <- sort(names(sarnovdf)[seq_along(arzidxs)])
# pdf(file = glue("{plot_path}AgeRelated_NoOverlapMB09BigSeries_lm_v10_DC.pdf"),
#     width = 8, height = 10)
# par(mfrow = c(2, 2))
# responses %>%
#   walk(~ invoke(
#     regression_plot,
#     pargs,
#     y = .,
#     data = sarnovdf,
#     main = paste(nrow(sarnovdf), "cases:"),
#     annotate = aroutdf_all$FC_1.25_AllCases$Probe_id
#   ))
# dev.off()

# responses <- sort(names(sarrdf)[seq_along(arzidxs)])
# pdf(file = glue("{plot_path}AgeRelated_raw_All1992Cases_lm_v11_DC.pdf"),
#     width = 8, height = 10)
# par(mfrow = c(2, 2))
# responses %>%
#   walk(~ invoke(
#     regression_plot,
#     pargs,
#     y = .,
#     data = sarrdf,
#     main = paste(nrow(sarrdf), "cases:"),
#     annotate = with(
#       aroutdf_all$FC_1.25_AllCases, Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)]
#     )
#   ))
# dev.off()


# Loess Scatterplots for arGeneSet ----------------------------------------

# Original Issues before generating results for entire arGeneSet:
# Issue #13: intClust subgroup graphs of ESR1 | ER−alpha | ILMN_1678535 with 99%
# CIs
# Issue #14: PAM50 and ER status subgroup graphs of ESR1 | ER−alpha |
# ILMN_1678535 with 99% CIs

# response <- "ER-alpha|ILMN_1678535"

# ** Setup ----------------------------------------------------------------

# FDR and Fold Changes
fdr <- 0.01
fcs <- c(1.25, 2, 4)

# File paths & names
root <- "main/code/Analysis_METABRIC_Expression_Trends/Plots/ProbeLevel/"
base_name1 <- "MBEX_{s}_{r}_FC_{fc}_FDR_{fdr}_arGeneSet_raw_{c}_loess_v01"
base_name2 <- "MBEX_{s}_{r}_FC_{fc}_FDR_{fdr}_arGeneSet_raw_{c}_loess_v02"
base_names <- list(base_name1, base_name2)
full_name <- "{root}{sub}{name}.pdf"

# PDF arguments for all-in-one/single plots
pdf_args_all <- list(width = 8, height = 10, useDingbats = FALSE)
pdf_args_single <- list(width = 6, height = 7, useDingbats = FALSE)

# Gene targets (with probe IDs), probe IDs only, ER binding genes
targets <- sort(names(sarrdf)[seq_along(arzidxs)])
probeids <- str_split_fixed(targets, "\\|", 2)[, 2]

# Plotting arguments common to all cases/subgroups
plot_args <- probeids %>% 
  purrr::set_names() %>% 
  map(~ list(
    x = "age_at_diagnosis",
    y = grep(., targets, value = TRUE),
    age = 60,
    conf.level = 1 - fdr,
    show.avg = TRUE,
    erbtxt = ifelse(
      with(annodf, Gene_symbol[Probe_id == .]) %in% erbnms,
      "(ER binding)",
      "(non-ER binding)"
    )
  ))

# NKI ER binding gene annotation text conditions
erbtxt_NKI <- probeids %>% 
  purrr::set_names() %>% 
  map_chr(~ {
    t1_lgl <- with(annodf, Gene_symbol[Probe_id == .]) %in% erbnms_t1
    t2_lgl <- with(annodf, Gene_symbol[Probe_id == .]) %in% erbnms_t2
    if (t1_lgl & t2_lgl) {
      "\nNKI Tier 1: ER binding; Tier 2: ER binding"
    } else if (!t1_lgl & t2_lgl) {
      "\nNKI Tier 1: non-ER binding; Tier 2: ER binding"
    } else {
      "\nNKI Tier 1: non-ER binding; Tier 2: non-ER binding"
    }
  })
setdiff(erbnms_t1, erbnms_t2) # only 3 conditions, erbnms_t1 subset of erbnms_t2

# Modified plotting arguments for NKI ER binding gene lists
plot_args_NKI <- map2(plot_args, erbtxt_NKI, function(x, y)
  modify_at(x, "erbtxt", ~ y))

probes <- list(plot_args, plot_args_NKI)


# ** All cases ------------------------------------------------------------

# Annotations for all cases
annotations <- aroutdf_all %>% 
  map(~ with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)]))

# All Cases
walk(fcs, function(fc) {
  walk2(probes, base_names, ~ {
    walk(.x, function(probe) {
      nm <- glue(.y, s = "All_cases", r = gsub("\\/", "_", probe$y),
                 fc = str_decimal(fc), fdr = str_decimal(fdr), c = "All1992Cases")
      fn <- glue(full_name, sub = "AllCases/", name = nm)
      main <- glue("All cases: N = {nrow(sarrdf)}\n\n")
      invoke(pdf, pdf_args_single, file = fn)
      invoke(regression_plot, probe, data = sarrdf, fc = fc, main = main,
             annotate = annotations[[paste0("FC_", fc, "_AllCases")]])
      dev.off()
    })
  })
})


# ** IntClust -------------------------------------------------------------

# Split by IntClust
data_groups <- sarrdf %>% split(.$iClustf)

# Count group sizes, check mean enclosure by loess CI bands (Issue #13)
# map_int(data_groups, nrow)
# targets %>% 
#   purrr::set_names() %>% 
#   map_df(~ map_lgl(data_groups, is_mean_enclosed, x = "age_at_diagnosis", y = .))

# Store directory, subtitles, and ER/HER2 specific annotations from aroutdf
sub <- "intClust"
main_1 <- "iClust {unique(.$iClustf)}: N = {nrow(.)}\n\n"
annotations <- aroutdf_intclust %>%
  map(~ if (nrow(.) == 0) {
    character(0)
  } else {
    with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)])
  })

# One PDF for all IntClust groups
walk(fcs, function(fc) {
  walk2(probes, base_names, ~ {
    walk(.x, function(probe) {
      nm <- glue(.y, s = "iClust_All", r = gsub("\\/", "_", probe$y),
                 fc = str_decimal(fc), fdr = str_decimal(fdr),
                 c = "All1992Cases")
      fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
      invoke(pdf, pdf_args_all, file = fn)
      par(mfrow = c(2, 2))
      iwalk(data_groups,
            ~ invoke(regression_plot, probe, data = .x, fc = fc,
                     main = glue(main_1, . = .x),
                     annotate = annotations[[paste("FC", fc, paste0("iClust", str_pad(.y, width = 2, pad = 0)), sep = "_")]]))
      dev.off()
    })
  })
})

# One PDF for each IntClust group
walk(fcs, function(fc) {
  walk(plot_args, function(pargs) {
    iwalk(data_groups, ~ {
      nm <- glue(
        base_name1,
        s = glue("iClust_{n}", n = str_pad(.y, width = 2, pad = 0)),
        r = gsub("\\/", "_", pargs$y),
        fc = str_decimal(fc),
        fdr = str_decimal(fdr),
        c = paste0(nrow(.x), "Cases")
      )
      fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
      invoke(pdf, pdf_args_single, file = fn)
      invoke(regression_plot, pargs, data = .x, fc = fc,
             main = glue(main_1, . = .x),
             annotate = annotations[[paste("FC", fc, paste0("iClust", str_pad(.y, width = 2, pad = 0)),
                                           sep = "_")]])
      dev.off()
    })
  })
})


# ** PAM50 ----------------------------------------------------------------

# Split by PAM50
data_groups <- sarrdf %>% split(.$BrCaf)

# Count group sizes, check mean enclosure by loess CI bands (Issue #13)
# map_int(data_groups, nrow)
# targets %>% 
#   purrr::set_names() %>% 
#   map_df(~ map_lgl(data_groups, is_mean_enclosed, x = "age_at_diagnosis", y = .))

# Store directory, subtitles, and ER/HER2 specific annotations from aroutdf
sub <- "Pam50"
main_1 <- "PAM50 {unique(.$BrCaf)}: N = {nrow(.)}\n\n"
annotations <- aroutdf_pam50 %>%
  map(~ if (nrow(.) == 0) {
    character(0)
  } else {
    with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)])
  })

# One PDF for all PAM50 groups
walk(fcs, function(fc) {
  walk2(probes, base_names, ~ {
    walk(.x, function(probe) {
      nm <- glue(.y, s = paste0(sub, "_All"),
                 r = gsub("\\/", "_", probe$y), fc = str_decimal(fc),
                 fdr = str_decimal(fdr), c = "All1992Cases")
      fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
      invoke(pdf, pdf_args_all, file = fn)
      par(mfrow = c(2, 2))
      iwalk(data_groups,
            ~ invoke(regression_plot, probe, data = .x, fc = fc,
                     main = glue(main_1, . = .x),
                     annotate = annotations[[paste("FC", fc, .y, sep = "_")]]))
      dev.off()
    })
  })
})

# One PDF for each PAM50 group
walk(fcs, function(fc) {
  walk(plot_args, function(pargs) {
    iwalk(data_groups, ~ {
      nm <- glue(base_name1,
                 s = glue("{sub}_{g}", g = .y),
                 r = gsub("\\/", "_", pargs$y),
                 fc = str_decimal(fc),
                 fdr = str_decimal(fdr),
                 c = paste0(nrow(.x), "Cases"))
      fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
      invoke(pdf, pdf_args_single, file = fn)
      invoke(regression_plot, pargs, data = .x, fc = fc,
             main = glue(main_1, . = .x),
             annotate = annotations[[paste("FC", fc, .y, sep = "_")]])
      dev.off()
    })
  })
})


# ** ER/HER2 --------------------------------------------------------------

# Split by ER/HER2
data_groups <- sarrdf %>% split(.$BrCaEHf)

# Count group sizes, check mean enclosure by loess CI bands (Issue #13)
# map_int(data_groups, nrow)
# targets %>% 
#   purrr::set_names() %>% 
#   map_df(~ map_lgl(data_groups, is_mean_enclosed, x = "age_at_diagnosis", y = .))

# Store directory, subtitles, and ER/HER2 specific annotations from aroutdf
sub <- "ER_HER2"
main_1 <- "{unique(.$BrCaEHf)}: N = {nrow(.)}\n\n"
annotations <- aroutdf_erher2 %>%
  map(~ if (nrow(.) == 0) {
    character(0)
  } else {
    with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)])
  })
# 
# # One PDF for all ER/HER2 status groups
# walk(fcs, function(fc) {
#   walk2(probes, base_names, ~ {
#     walk(.x, function(probe) {
#       nm <- glue(.y, s = paste0(sub, "_All"),
#                  r = gsub("\\/", "_", probe$y), fc = str_decimal(fc),
#                  fdr = str_decimal(fdr), c = "All1992Cases")
#       fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
#       invoke(pdf, pdf_args_all, file = fn)
#       par(mfrow = c(2, 2))
#       iwalk(data_groups,
#             ~ invoke(regression_plot, probe, data = .x, fc = fc,
#                      main = glue(main_1, . = .x),
#                      annotate = annotations[[paste("FC", fc, str_sign(.y), sep = "_")]]))
#       dev.off()
#     })
#   })
# })

# One PDF for each ER/HER2 status group
walk(fcs, function(fc) {
  walk(plot_args, function(pargs) {
    iwalk(data_groups, ~ {
      nm <- glue(base_name1, s = str_sign(.y, type = "long"),
                 r = gsub("\\/", "_", pargs$y), fc = str_decimal(fc),
                 fdr = str_decimal(fdr), c = paste0(nrow(.x), "Cases"))
      fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
      invoke(pdf, pdf_args_single, file = fn)
      invoke(regression_plot, pargs, data = .x, fc = fc,
             main = glue(main_1, . = .x),
             annotate = annotations[[paste("FC", fc, str_sign(.y), sep = "_")]])
      dev.off()
    })
  })
})


# ** ER+/ER- --------------------------------------------------------------

# Split by ER+/ER-
data_groups <- sarrdf %>% split(.$BrCaEf)

# Count group sizes, check mean enclosure by loess CI bands (Issue #13)
# map_int(data_groups, nrow)
# targets %>% 
#   purrr::set_names() %>% 
#   map_df(~ map_lgl(data_groups, is_mean_enclosed, x = "age_at_diagnosis", y = .))

# Store directory, subtitles, and ER+/ER- specific annotations from aroutdf
sub <- "ER"
main_1 <- "{unique(.$BrCaEf)}: N = {nrow(.)}\n\n"
annotations <- aroutdf_er %>%
  map(~ if (nrow(.) == 0) {
    character(0)
  } else {
    with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)])
  })

# # One PDF for all ER+/ER- status groups
# walk(fcs, function(fc) {
#   walk2(probes, base_names, ~ {
#     walk(.x, function(probe) {
#       nm <- glue(.y, s = paste0(sub, "_All"),
#                  r = gsub("\\/", "_", probe$y), fc = str_decimal(fc),
#                  fdr = str_decimal(fdr), c = "All1992Cases")
#       fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
#       invoke(pdf, pdf_args_all, file = fn)
#       par(mfrow = c(2, 2))
#       iwalk(data_groups,
#             ~ invoke(regression_plot, probe, data = .x, fc = fc,
#                      main = glue(main_1, . = .x),
#                      annotate = annotations[[paste("FC", fc, str_sign(.y), sep = "_")]]))
#       dev.off()
#     })
#   })
# })

# One PDF for each ER+/ER- status group
walk(fcs, function(fc) {
  walk(plot_args, function(pargs) {
    iwalk(data_groups, ~ {
      nm <- glue(base_name1, s = str_sign(.y, type = "long"),
                 r = gsub("\\/", "_", pargs$y), fc = str_decimal(fc),
                 fdr = str_decimal(fdr), c = paste0(nrow(.x), "Cases"))
      fn <- glue(full_name, sub = paste0(sub, "/"), name = nm)
      invoke(pdf, pdf_args_single, file = fn)
      invoke(regression_plot, pargs, data = .x, fc = fc,
             main = glue(main_1, . = .x),
             annotate = annotations[[paste("FC", fc, str_sign(.y), sep = "_")]])
      dev.off()
    })
  })
})


# Volcano Plots -----------------------------------------------------------

memoized_path <- "main/code/Analysis_METABRIC_Expression_Trends/Memoized"
fdr <- 0.01
legend_position <- "topleft"
set.seed(317) ## hclust thinning uses sample()

# ** All cases ------------------------------------------------------------

probes_all <- map(dir_ls(path = memoized_path, regexp = "All_cases.*rds"), readRDS)
probes_all <- probes_all %>% 
  map(~ mutate(
    ., 
    !!"PBHadj_log" := -log10(.data$Best_BHadj_pval),
    !!"ERbinding1" := .data$Gene_symbol %in% erbnms_t1,
    !!"ERbinding2" := .data$Gene_symbol %in% erbnms_t2,
    !!"vcp" := abs(.data$Best_log2FC) > log2(1.25) & .data$Best_BHadj_pval < fdr,
    !!"binding_group" := case_when(
      .data$vcp & .data$ERbinding1  & .data$ERbinding2 ~ "Tier 1",
      .data$vcp & !.data$ERbinding1 & .data$ERbinding2 ~ "Tier 2 Only",
      .data$vcp & !.data$ERbinding1 & !.data$ERbinding2 ~ "Non binding",
      TRUE ~ " "
    )
  )) %>% 
  map(metabric_evcp, c(log2(4), 10, log2(3), 25, log2(2), 40))

root.dir <- "main/code/Analysis_METABRIC_Expression_Trends/Plots/Volcano/AllCases"
plot.suffix <- "AgeRelated_Volcano_All_cases_FDR_0p01_All1992Cases_v02.pdf"
title <- "METABRIC: All cases"
metabric_volcano(probes_all, root.dir, plot.suffix, title, thin = TRUE, legloc = legend_position)


# ** IntClust -------------------------------------------------------------

probes_intclust <- map(dir_ls(path = memoized_path, regexp = "iClust.*rds"), readRDS)

th_intclust <- list(
  iClust01 = c(log2(3), -log10(0.01), log2(2.5), -log10(0.01), log2(1.25), -log10(0.01)), 
  iClust02 = c(log2(3), -log10(0.01), log2(2.5), -log10(0.01), log2(1.25), -log10(0.01)), 
  iClust03 = c(log2(4), 2, log2(3), 3, log2(1.25), 5), 
  iClust04 = c(log2(4), -log10(0.01), log2(2.5), -log10(0.01/5), log2(1.25), 3.5), 
  iClust05 = c(log2(3), -log10(0.01), log2(2.5), -log10(0.01), log2(1.25), -log10(0.01)), 
  iClust06 = c(log2(3), -log10(0.01), log2(2.5), -log10(0.01), log2(1.25), -log10(0.01)), 
  iClust07 = c(log2(3), -log10(0.01), log2(2.5), -log10(0.01), log2(1.25), -log10(0.01)), 
  iClust08 = c(log2(4), -log10(0.01), log2(3), 3, log2(1.25), 5), 
  iClust09 = c(log2(3), -log10(0.01), log2(2.5), -log10(0.01), log2(1.25), -log10(0.01)),   
  iClust10 = c(log2(2.25), -log10(0.01), log2(1.75), -log10(0.01), log2(1.25), -log10(0.01))
)

probes_intclust <- probes_intclust %>% 
  `[`(gtools::mixedsort(names(.))) %>% 
  map(~ mutate(
    ., 
    !!"PBHadj_log" := -log10(.data$Best_BHadj_pval),
    !!"ERbinding1" := .data$Gene_symbol %in% erbnms_t1,
    !!"ERbinding2" := .data$Gene_symbol %in% erbnms_t2,
    !!"vcp" := abs(.data$Best_log2FC) > log2(1.25) & .data$Best_BHadj_pval < fdr,
    !!"binding_group" := case_when(
      .data$vcp & .data$ERbinding1  & .data$ERbinding2 ~ "Tier 1",
      .data$vcp & !.data$ERbinding1 & .data$ERbinding2 ~ "Tier 2 Only",
      .data$vcp & !.data$ERbinding1 & !.data$ERbinding2 ~ "Non binding",
      TRUE ~ " "
    )
  )) %>% 
  map2(th_intclust, metabric_evcp)

root.dir <- "main/code/Analysis_METABRIC_Expression_Trends/Plots/Volcano/intClust"
title <- paste("METABRIC: iClust", seq_along(probes_intclust))

# Scale free
plot.suffix <- paste0("AgeRelated_Volcano_iClust_", seq_along(probes_intclust),
                      "_FDR_0p01_All1992Cases_v02.pdf")
metabric_volcano(probes_intclust, root.dir, plot.suffix, title, thin = TRUE, legloc = legend_position)

# Scale fixed
plot.suffix <- paste0("AgeRelated_Volcano_iClust_", seq_along(probes_intclust),
                      "_FDR_0p01_All1992Cases_v03.pdf")
xlim <- max_xlim(probes_intclust)
ylim <- max_ylim(probes_intclust)
metabric_volcano(probes_intclust, root.dir, plot.suffix, title,
                 xlim = xlim, ylim = ylim, thin = TRUE, legloc = legend_position)


# ** PAM50 ----------------------------------------------------------------

probes_pam50 <- map(dir_ls(path = memoized_path, regexp = "Pam50.*rds"), readRDS)

th_pam50 <- list(
  Basal = c(log2(4), -log10(0.01), log2(2.5), -log10(0.01/5), log2(1.25), 3.5), 
  Her2 = c(log2(4), -log10(0.01), log2(3), -log10(0.01), log2(1.25), -log10(0.01)), 
  LumA = c(log2(3.5), -log10(0.01), log2(2.5), 10, log2(1.25), 13),
  LumB = c(log2(3), -log10(0.01), log2(2), 4, log2(1.25), 5),
  Normal = c(log2(3), -log10(0.01), log2(2.5), -log10(0.01), log2(1.25), -log10(0.01))
)

probes_pam50 <- probes_pam50 %>% 
  map(~ mutate(
    ., 
    !!"PBHadj_log" := -log10(.data$Best_BHadj_pval),
    !!"ERbinding1" := .data$Gene_symbol %in% erbnms_t1,
    !!"ERbinding2" := .data$Gene_symbol %in% erbnms_t2,
    !!"vcp" := abs(.data$Best_log2FC) > log2(1.25) & .data$Best_BHadj_pval < fdr,
    !!"binding_group" := case_when(
      .data$vcp & .data$ERbinding1  & .data$ERbinding2 ~ "Tier 1",
      .data$vcp & !.data$ERbinding1 & .data$ERbinding2 ~ "Tier 2 Only",
      .data$vcp & !.data$ERbinding1 & !.data$ERbinding2 ~ "Non binding",
      TRUE ~ " "
    )
  )) %>% 
  map2(th_pam50, metabric_evcp)

root.dir <- "main/code/Analysis_METABRIC_Expression_Trends/Plots/Volcano/Pam50"
title <- paste("METABRIC:", names(th_pam50))

# Scale free
plot.suffix <- paste0("AgeRelated_Volcano_Pam50_", names(th_pam50),
                      "_FDR_0p01_All1992Cases_v02.pdf")
metabric_volcano(probes_pam50, root.dir, plot.suffix, title, thin = TRUE, legloc = legend_position)

# Scale fixed
plot.suffix <- paste0("AgeRelated_Volcano_Pam50_", names(th_pam50),
                      "_FDR_0p01_All1992Cases_v03.pdf")
xlim <- max_xlim(probes_pam50)
ylim <- max_ylim(probes_pam50)
metabric_volcano(probes_pam50, root.dir, plot.suffix, title,
                 xlim = xlim, ylim = ylim, thin = TRUE, legloc = legend_position)


# ** ER/HER2 --------------------------------------------------------------

probes_erher2 <- map(dir_ls(path = memoized_path, regexp = "HER.*rds"), readRDS)

th_erher2 <- list(
  `ER-/HER2-` = c(log2(2.25), 2, log2(1.75), 2, log2(1.25), 2),
  `ER-/HER2+` = c(log2(4), 10, log2(3), 25, log2(2), 40),
  `ER+/HER2-` = c(log2(3), 10, log2(2.5), 20, log2(1.25), 25), 
  `ER+/HER2+` = c(log2(4), 10, log2(3), 25, log2(2), 40)
)

probes_erher2 <- probes_erher2 %>% 
  map(~ mutate(
    ., 
    !!"PBHadj_log" := -log10(.data$Best_BHadj_pval),
    !!"ERbinding1" := .data$Gene_symbol %in% erbnms_t1,
    !!"ERbinding2" := .data$Gene_symbol %in% erbnms_t2,
    !!"vcp" := abs(.data$Best_log2FC) > log2(1.25) & .data$Best_BHadj_pval < fdr,
    !!"binding_group" := case_when(
      .data$vcp & .data$ERbinding1  & .data$ERbinding2 ~ "Tier 1",
      .data$vcp & !.data$ERbinding1 & .data$ERbinding2 ~ "Tier 2 Only",
      .data$vcp & !.data$ERbinding1 & !.data$ERbinding2 ~ "Non binding",
      TRUE ~ " "
    )
  )) %>% 
  map2(th_erher2, metabric_evcp)

root.dir <- "main/code/Analysis_METABRIC_Expression_Trends/Plots/Volcano/ER_HER2"
title <- paste("METABRIC:", names(th_erher2))

# Scale free
plot.suffix <- paste0("AgeRelated_Volcano_", str_sign(names(th_erher2), type = "long"),
                      "_FDR_0p01_All1992Cases_v02.pdf")
metabric_volcano(probes_erher2, root.dir, plot.suffix, title, thin = TRUE, legloc = legend_position)

# Scale fixed
plot.suffix <- paste0("AgeRelated_Volcano_", str_sign(names(th_erher2), type = "long"),
                      "_FDR_0p01_All1992Cases_v03.pdf")
xlim <- max_xlim(probes_erher2)
ylim <- max_ylim(probes_erher2)
metabric_volcano(probes_erher2, root.dir, plot.suffix, title,
                 xlim = xlim, ylim = ylim, thin = TRUE, legloc = legend_position)


# Loess Scatterplots for volcano labels -----------------------------------

# ** Setup ----------------------------------------------------------------

# FDR and Fold Changes
fdr <- 0.01
fcs <- c(1.25, 2, 4)

# File paths & names
root <- "main/code/Analysis_METABRIC_Expression_Trends/Plots/ProbeLevel/VolcanoLabs/"
dir_create(root)
base_name1 <- "MBEX_{s}_{r}_FC_{fc}_FDR_{fdr}_VolcanoLabs_{c}_loess_v01"
base_name2 <- gsub("v01", "v02", base_name1)
base_names <- list(base_name1, base_name2)
full_name <- "{root}{sub}{name}.pdf"

# PDF arguments for all-in-one/single plots
pdf_args_all <- list(width = 8, height = 10, useDingbats = FALSE)
pdf_args_single <- list(width = 6, height = 7, useDingbats = FALSE)

# Common plotting arguments
common_args <- list(
  x = "age_at_diagnosis",
  age = 60,
  conf.level = 1 - fdr,
  show.avg = TRUE
)

# Reorder clinical data by same order as probes
cldvdf_match <- cldvdf %>% 
  extract(match(rownames(Dataset_r), .$MBid), )


# ** All cases ------------------------------------------------------------

dir_create(path(root, "AllCases"))

# Extract gene and probe names for volcano labels
vlabs <- probes_all %>% 
  map(~ {
    dplyr::filter(., .data$evcp & .data$binding_group != " ") %>%
      dplyr::select(.data$ERbinding, .data$binding_group, .data$Gene_symbol, .data$Probe_id) %>% 
      dplyr::mutate(GP = paste(.data$Gene_symbol, .data$Probe_id, sep = "|"),
                    ERbinding = ifelse(.data$ERbinding, "(ER binding)", "(non-ER binding)"),
                    binding_group = case_when(
                      .data$binding_group == "Tier 1" ~ "\nNKI Tier 1: ER binding; Tier 2: ER binding",
                      .data$binding_group == "Tier 2 Only" ~ "\nNKI Tier 1: non-ER binding; Tier 2: ER binding",
                      .data$binding_group == "Non binding" ~ "\nNKI Tier 1: non-ER binding; Tier 2: non-ER binding"
                    )) %>%
      dplyr::arrange(.data$Gene_symbol)
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::distinct()

# Plotting arguments for vlabs
plot_args2 <- vlabs %>% 
  pmap(~ list(y = ..5, erbtxt = ..1)) %>% 
  map(list_merge, !!!common_args) %>% 
  map(`[`, c("x", "y", "age", "conf.level", "show.avg", "erbtxt")) %>%
  set_names(vlabs$Probe_id)

# Modified plotting arguments for NKI ER binding gene lists
plot_args_NKI2 <- map2(plot_args2, vlabs$binding_group, ~ list_modify(.x, erbtxt = .y))
probes2 <- list(plot_args2, plot_args_NKI2)

# Volcano labels sarrdf
vsarrdf <- Dataset_r[, vlabs$Probe_id] %>%
  `colnames<-`(vlabs$GP) %>%
  cbind(cldvdf_match)

# Annotations for all cases
annotations <- aroutdf_all %>% 
  map(~ with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)]))

# All Cases
walk(fcs, function(fc) {
  walk2(probes2, base_names, ~ {
    walk(.x, function(probe) {
      pfc <- str_decimal(fc)
      pfdr <- str_decimal(fdr)
      nm <- glue(.y, s = "All_cases", r = gsub("\\/", "_", probe$y),
                 fc = pfc, fdr = pfdr, c = "All1992Cases")
      fn <- glue(full_name, sub = "AllCases/", name = nm)
      main <- glue("All cases: N = {nrow(vsarrdf)}\n\n")
      invoke(pdf, pdf_args_single, file = fn)
      invoke(regression_plot, probe, data = vsarrdf, fc = fc, main = main,
             annotate = annotations[[paste0("FC_", fc, "_AllCases")]])
      dev.off()
    })
  })
})


# ** IntClust -------------------------------------------------------------

dir_create(path(root, "intClust"))

# Extract gene and probe names for volcano labels
vlabs <- probes_intclust %>%
  map(~ {
    dplyr::filter(., .data$evcp & .data$binding_group != " ") %>%
      dplyr::select(.data$ERbinding, .data$binding_group, .data$Gene_symbol, .data$Probe_id) %>% 
      dplyr::mutate(GP = paste(.data$Gene_symbol, .data$Probe_id, sep = "|"),
                    ERbinding = ifelse(.data$ERbinding, "(ER binding)", "(non-ER binding)"),
                    binding_group = case_when(
                      .data$binding_group == "Tier 1" ~ "\nNKI Tier 1: ER binding; Tier 2: ER binding",
                      .data$binding_group == "Tier 2 Only" ~ "\nNKI Tier 1: non-ER binding; Tier 2: ER binding",
                      .data$binding_group == "Non binding" ~ "\nNKI Tier 1: non-ER binding; Tier 2: non-ER binding"
                    )) %>%
      dplyr::arrange(.data$Gene_symbol)
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::distinct()

# Plotting arguments for vlabs
plot_args2 <- vlabs %>% 
  pmap(~ list(y = ..5, erbtxt = ..1)) %>% 
  map(list_merge, !!!common_args) %>% 
  map(`[`, c("x", "y", "age", "conf.level", "show.avg", "erbtxt")) %>%
  set_names(vlabs$Probe_id)

# Modified plotting arguments for NKI ER binding gene lists
plot_args_NKI2 <- map2(plot_args2, vlabs$binding_group, ~ list_modify(.x, erbtxt = .y))
probes2 <- list(plot_args2, plot_args_NKI2)

# Volcano labels sarrdf
vsarrdf <- Dataset_r[, vlabs$Probe_id] %>%
  `colnames<-`(vlabs$GP) %>%
  cbind(cldvdf_match)

# Split by IntClust
data_groups <- vsarrdf %>% split(.$iClustf)

# Store directory, subtitles, and ER/HER2 specific annotations from aroutdf
ssub <- "intClust"
main_1 <- "iClust {unique(.$iClustf)}: N = {nrow(.)}\n\n"
annotations <- aroutdf_intclust %>%
  map(~ if (nrow(.) == 0) {
    character(0)
  } else {
    with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)])
  })

# One PDF for all IntClust groups
walk(fcs, function(fc) {
  walk2(probes2, base_names, ~ {
    walk(.x, function(probe) {
      pfc <- str_decimal(fc)
      pfdr <- str_decimal(fdr)
      nm <- glue(.y, s = "iClust_All", r = gsub("\\/", "_", probe$y),
                 fc = pfc, fdr = pfdr, c = "All1992Cases")
      fn <- glue(full_name, sub = paste0(ssub, "/"), name = nm)
      invoke(pdf, pdf_args_all, file = fn)
      par(mfrow = c(2, 2))
      iwalk(data_groups,
            ~ invoke(regression_plot, probe, data = .x, fc = fc,
                     main = glue(main_1, . = .x),
                     annotate = annotations[[paste("FC", fc, paste0("iClust", str_pad(.y, width = 2, pad = 0)), sep = "_")]]))
      dev.off()
    })
  })
})


# ** PAM50 ----------------------------------------------------------------

dir_create(path(root, "Pam50"))

vlabs <- probes_pam50 %>% 
  map(~ {
    dplyr::filter(., .data$evcp & .data$binding_group != " ") %>%
      dplyr::select(.data$ERbinding, .data$binding_group, .data$Gene_symbol, .data$Probe_id) %>% 
      dplyr::mutate(GP = paste(.data$Gene_symbol, .data$Probe_id, sep = "|"),
                    ERbinding = ifelse(.data$ERbinding, "(ER binding)", "(non-ER binding)"),
                    binding_group = case_when(
                      .data$binding_group == "Tier 1" ~ "\nNKI Tier 1: ER binding; Tier 2: ER binding",
                      .data$binding_group == "Tier 2 Only" ~ "\nNKI Tier 1: non-ER binding; Tier 2: ER binding",
                      .data$binding_group == "Non binding" ~ "\nNKI Tier 1: non-ER binding; Tier 2: non-ER binding"
                    )) %>%
      dplyr::arrange(.data$Gene_symbol)
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::distinct()

# Plotting arguments for vlabs
plot_args2 <- vlabs %>% 
  pmap(~ list(y = ..5, erbtxt = ..1)) %>% 
  map(list_merge, !!!common_args) %>% 
  map(`[`, c("x", "y", "age", "conf.level", "show.avg", "erbtxt")) %>%
  set_names(vlabs$Probe_id)

# Modified plotting arguments for NKI ER binding gene lists
plot_args_NKI2 <- map2(plot_args2, vlabs$binding_group, ~ list_modify(.x, erbtxt = .y))
probes2 <- list(plot_args2, plot_args_NKI2)

# Volcano labels sarrdf
vsarrdf <- Dataset_r[, vlabs$Probe_id] %>%
  `colnames<-`(vlabs$GP) %>%
  cbind(cldvdf_match)

# Split by PAM50
data_groups <- vsarrdf %>% split(.$BrCaf)

# Store directory, subtitles, and ER/HER2 specific annotations from aroutdf
ssub <- "Pam50"
main_1 <- "PAM50 {unique(.$BrCaf)}: N = {nrow(.)}\n\n"
annotations <- aroutdf_pam50 %>%
  map(~ if (nrow(.) == 0) {
    character(0)
  } else {
    with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)])
  })

# One PDF for all PAM50 groups
walk(fcs, function(fc) {
  walk2(probes2, base_names, ~ {
    walk(.x, function(probe) {
      pfc <- str_decimal(fc)
      pfdr <- str_decimal(fdr)
      nm <- glue(.y, s = paste0(ssub, "_All"), r = gsub("\\/", "_", probe$y),
                 fc = pfc, fdr = pfdr, c = "All1992Cases")
      fn <- glue(full_name, sub = paste0(ssub, "/"), name = nm)
      invoke(pdf, pdf_args_all, file = fn)
      par(mfrow = c(2, 2))
      iwalk(data_groups,
            ~ invoke(regression_plot, probe, data = .x, fc = fc,
                     main = glue(main_1, . = .x),
                     annotate = annotations[[paste("FC", fc, .y, sep = "_")]]))
      dev.off()
    })
  })
})


# ** ER/HER2 --------------------------------------------------------------

dir_create(path(root, "ER_HER2"))

vlabs <- probes_erher2 %>% 
  map(~ {
    dplyr::filter(., .data$evcp & .data$binding_group != " ") %>%
      dplyr::select(.data$ERbinding, .data$binding_group, .data$Gene_symbol, .data$Probe_id) %>% 
      dplyr::mutate(GP = paste(.data$Gene_symbol, .data$Probe_id, sep = "|"),
                    ERbinding = ifelse(.data$ERbinding, "(ER binding)", "(non-ER binding)"),
                    binding_group = case_when(
                      .data$binding_group == "Tier 1" ~ "\nNKI Tier 1: ER binding; Tier 2: ER binding",
                      .data$binding_group == "Tier 2 Only" ~ "\nNKI Tier 1: non-ER binding; Tier 2: ER binding",
                      .data$binding_group == "Non binding" ~ "\nNKI Tier 1: non-ER binding; Tier 2: non-ER binding"
                    )) %>%
      dplyr::arrange(.data$Gene_symbol)
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::distinct()

# Plotting arguments for vlabs
plot_args2 <- vlabs %>% 
  pmap(~ list(y = ..5, erbtxt = ..1)) %>% 
  map(list_merge, !!!common_args) %>% 
  map(`[`, c("x", "y", "age", "conf.level", "show.avg", "erbtxt")) %>%
  set_names(vlabs$Probe_id)

# Modified plotting arguments for NKI ER binding gene lists
plot_args_NKI2 <- map2(plot_args2, vlabs$binding_group, ~ list_modify(.x, erbtxt = .y))
probes2 <- list(plot_args2, plot_args_NKI2)

# Volcano labels sarrdf
vsarrdf <- Dataset_r[, vlabs$Probe_id] %>%
  `colnames<-`(vlabs$GP) %>%
  cbind(cldvdf_match) %>% 
  mutate(
    BrCaEHf = case_when(
      ER.Expr == "+" & Her2.Expr == "+" ~ "ER+/HER2+",
      ER.Expr == "+" & Her2.Expr == "-" ~ "ER+/HER2-",
      ER.Expr == "-" & Her2.Expr == "+" ~ "ER-/HER2+",
      ER.Expr == "-" & Her2.Expr == "-" ~ "ER-/HER2-",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))
  )

# Split by ER/HER2
data_groups <- vsarrdf %>% split(.$BrCaEHf)

# Store directory, subtitles, and ER/HER2 specific annotations from aroutdf
ssub <- "ER_HER2"
main_1 <- "{unique(.$BrCaEHf)}: N = {nrow(.)}\n\n"
annotations <- aroutdf_erher2 %>%
  map(~ if (nrow(.) == 0) {
    character(0)
  } else {
    with(., Probe_id[which(arGeneSet_BHadj_and_AgeDependentp)])
  })

# One PDF for all ER/HER2 status groups
walk(fcs, function(fc) {
  walk2(probes2, base_names, ~ {
    walk(.x, function(probe) {
      pfc <- str_decimal(fc)
      pfdr <- str_decimal(fdr)
      nm <- glue(.y, s = paste0(ssub, "_All"), r = gsub("\\/", "_", probe$y),
                 fc = pfc, fdr = pfdr, c = "All1992Cases")
      fn <- glue(full_name, sub = paste0(ssub, "/"), name = nm)
      invoke(pdf, pdf_args_all, file = fn)
      par(mfrow = c(2, 2))
      iwalk(data_groups,
            ~ invoke(regression_plot, probe, data = .x, fc = fc,
                     main = glue(main_1, . = .x),
                     annotate = annotations[[paste("FC", fc, str_sign(.y), sep = "_")]]))
      dev.off()
    })
  })
})


# Bootstrap Analysis ------------------------------------------------------

# Set up inputs for bootstrap analysis
data <- sehrdf %>% select(MBid, age_at_diagnosis, iClustf)
x <- "age_at_diagnosis"
y <- Dataset_r
group <- data$iClustf
probes <- annodf$ProbeId
genes <- annodf$Gene_symbol
base_path <- "main/code/Analysis_METABRIC_Expression_Trends/Bootstrap/Counts/"
base_name <- "AgeDependent_Bootstrap_Counts"
base_title <- "Proportion of Age-Associated Gene Targets per IntClust group
on 20 replicates of bootstrapped samples of size"

# First round of bootstrap sizes: nb = 50, 100, 200
for (nb in c(50, 100, 200)) {
  date()
  res <- boot_lm(
    data, x, y, group, probes, genes, nb = nb,
    filename = glue("{base_path}{base_name}_NB{fnb}.rds",
                    fnb = str_pad(nb, width = 3, pad = 0))
  )
  date()
  boot_boxplot(res$counts_probes, title = glue("{base_title} {nb}"))
}

# Run Times
# nb = 50   (~ 25 hours)
# [1] "Thu Apr 13 12:09:14 2017"
# [1] "Fri Apr 14 12:58:24 2017"

# nb = 100   (~ 26 hours)
# [1] "Fri Apr 14 12:58:25 2017"
# [1] "Sat Apr 15 13:54:42 2017"

# nb = 200   (~ 26.5 hours)
# [1] "Sat Apr 15 13:54:43 2017"
# [1] "Sun Apr 16 15:17:31 2017"

# Second round of bootstrap sizes: nb = 75, 150, 300 (better representation of
# actual IntClust group sizes) run on three concurrent RStudio sessions
for (nb in c(75, 150, 300)) {
  date()
  res <- boot_lm(
    data, x, y, group, probes, genes, nb = nb,
    filename = glue("{base_path}{base_name}_NB{fnb}.rds",
                    fnb = str_pad(nb, width = 3, pad = 0))
  )
  date()
  boot_boxplot(res$counts_probes, title = glue("{base_title} {nb}"))
}

# Run Times
# nb <- 75  (~ 30 hours)
# [1] "Tue Apr 18 16:16:02 2017"
# [1] "Wed Apr 19 22:01:02 2017"

# nb <- 150   (~ 30 hours)
# [1] "Tue Apr 18 16:17:25 2017"
# [1] "Wed Apr 19 22:09:23 2017"

# nb <- 300   (~ 30 hours)
# [1] "Tue Apr 18 16:23:57 2017"
# [1] "Wed Apr 19 22:48:55 2017"


# Additional Bootstrap Checks ---------------------------------------------

# Test 1 is a sanity check to ensure that when we use the raw
# IntClust sizes the results are what we expect from previous analyses
# Test 2 extracts extra metadata for the nb = 300 case
base_path <- "main/code/Analysis_METABRIC_Expression_Trends/Bootstrap/Counts/"
base_name <- "AgeDependent_Bootstrap_Counts"

# Test 1: One replicate of group-specific bootstrap sizes (~ 1.4 hours)
date()  # [1] "Tue Apr 25 14:54:26 2017"
res.test1 <- boot_lm(data, x, y, group, probes, genes, reps = 1,
                     filename = glue("{base_path}{base_name}_NBtest1.rds"))
date()  # [1] "Tue Apr 25 16:20:14 2017"
boot_boxplot(res.test1$counts_probes,
             title = "Proportion of Age-Associated Gene Targets per IntClust group
on Group-specific bootstrapped sample sizes")

# Test 2: Two reps of NB = 300, with counts for each group (~ 2.8 hours)
date()  # [1] "Tue Apr 25 14:53:47 2017"
res.test2 <- boot_lm(data, x, y, group, probes, genes, reps = 2, nb = 300,
                     filename = glue("{base_path{base_name}_NBtest2.rds"))
date()  # [1] "Tue Apr 25 17:36:00 2017"
boot_boxplot(res.test2$counts_probes,
             title = "Proportion of Age-Associated Gene Targets per IntClust group
on 2 replicates of bootstrapped samples of size 300")


# Linear Regression Plots -------------------------------------------------

# Constant inputs
seed <- 1
reps <- 1
nb <- 300
subset <- 60
fc <- 1.25
fdr <- 0.01
biosigslopeLE <- log2(fc) / 35
biosigslopeGT <- log2(fc) / 35
biosigslopeAll <- log2(fc) / 70

# Read in data inputs if they don't already exist in current environment
base.dir <- "main/data/Data_METABRIC_Expression_Trends/"
walk(c("sehrdf", "ehzidxs", "Dataset_r"),
     ~ if (!exists(.)) assign(., read_rds(glue("{base.dir}{.}.rds")), 1))

# Object inputs
data <- sehrdf %>% select(MBid, age_at_diagnosis, iClustf)
x <- "age_at_diagnosis"
y <- Dataset_r
group <- data$iClustf
responses <- sort(names(sehrdf)[seq_along(ehzidxs)])
ilmn_id <- gsub(".*\\|", "", responses)

# Construct bootstrap
set.seed(seed)
boot_dat <- levels(group) %>% 
  purrr::set_names() %>% 
  map(~ data[group == .x, ]) %>% 
  map(~ rerun(reps, sample_n(.x, nb, replace = TRUE)))

# Get lm test results
test <- boot_dat %>% 
  transpose() %>% 
  extract2(1) %>% 
  map(function(d) {
    map(ilmn_id, function(p) {
      df <- cbind(d, probeset = y[match(d[["MBid"]], rownames(y)),
                                  match(p, colnames(y))])
      list(df[[x]] <= subset, df[[x]] > subset, NULL) %>%
        map(~ lm(df[["probeset"]] ~ df[[x]], na.action = na.exclude,
                 subset = .x)) %>%
        lapply(broom::tidy) %>%
        bind_rows() %>%
        magrittr::extract(.$term %in% "df[[x]]", c("estimate", "p.value")) %>% 
        flatten_dbl() %>% 
        purrr::set_names(c("LE_slope", "GT_slope", "All_slope",
                           "LE_pval", "GT_pval", "All_pval"))
    }) %>% 
      invoke(rbind, .) %>% 
      as.data.frame() %>% 
      mutate(
        LE_pval_BHadj = p.adjust(LE_pval, method = "BH"),
        GT_pval_BHadj = p.adjust(GT_pval, method = "BH"),
        All_pval_BHadj = p.adjust(All_pval, method = "BH")
      ) %>%
      metabric_signif("BHadj_and_AgeDependentp") %>%
      mutate(PROBE = ilmn_id, GENE = responses)
  })


# Add annotations
df.lm <- test %>% 
  map(~ {
    .x %>% 
      mutate_at(vars(matches("pval")),
                funs(vapply(., format, digits = 5, character(1)))) %>% 
      transmute(
        GP = GENE,
        PROBE = PROBE,
        LE60_pval = glue("[<={subset}] p = {LE_pval}"),
        GT60_pval = glue("[>{subset}] p = {GT_pval}"),
        AllAges_pval = glue("[All] p = {All_pval}"),
        Trend = ifelse(
          BHadj_and_AgeDependentp,
          glue("Age-dependent trend:\n\n|FC|>{fc} and adjPval<{fdr}"),
          glue("No detectable trend:\n\n|FC|<{fc} or adjPval>{fdr}")
        )
      )
  })

# Make into long format for plotting
dat.plot <- boot_dat %>%
  map(1) %>%
  map(~ {
    .x %>%
      rename(Age = age_at_diagnosis) %>%
      cbind(y[match(.[["MBid"]], rownames(y)), ilmn_id]) %>%
      gather(PROBE, Expression, -c(1:3))
  }) %>% 
  map2(df.lm, merge) %>% 
  map(~ split(.x, .x[["GP"]]))

# Create trend plots for each IntClust (01-10) group (~ 20 mins)
dat.plot %>% 
  walk2(., str_pad(names(.), width = 2, pad = 0), ~ {
    ap <- .x %>% imap(~ scatterfit_plot(x = .x, title = .y))
    grobs <- gridExtra::marrangeGrob(ap, nrow = 2, ncol = 2)
    ggsave(
      filename = glue("{plot_path}AgeRelated_lm_bootstrap_IntClust{.y}.pdf"),
      plot = grobs, width = 8, height = 10)
  })


# Bootstrap Result Matrices -----------------------------------------------

# Read in data inputs if they don't already exist in current environment
base.dir <- "main/data/Data_METABRIC_Expression_Trends/"
walk(c("annodf", "sehrdf", "aridxs", "Dataset_r"), ~
       if (!exists(.)) {
         nm <- ifelse(. == "annodf",
                      "Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09",
                      .)
         assign(., read_rds(glue("{base.dir}{nm}.rds")), 1)
       })

# Object inputs
data <- sehrdf %>% select(MBid, age_at_diagnosis, iClustf)
x <- "age_at_diagnosis"
y <- Dataset_r
group <- data$iClustf
probes <- annodf$ProbeId
genes <- annodf$Gene_symbol

# Write all outputs to file (6 NB * 20 Reps = 120 Objects)
reps <- 20
padded_reps <- str_pad(seq_len(reps), width = 2, pad = 0)
base_path <- "main/code/Analysis_METABRIC_Expression_Trends/Bootstrap/Stats/"
base_name <- "AgeDependent_Bootstrap_Stats"

# Generate bootstrap result statistics matrices for six NB sizes
for (nb in c(50, 75, 100, 150, 200, 300)) {
  date()
  boot_data <- boot_generate(data, group, nb = nb, reps = reps)
  vals <- boot_stats(boot_data, x, y, probes, genes, aridxs)
  walk2(vals, padded_reps,
        ~ write_rds(.x,
                    path = glue("{base_path}{base_name}_NB{fnb}_Rep{.y}.rds",
                                fnb = str_pad(nb, width = 3, pad = 0)),
                    compress = "xz")
  )
  date()
}

# Run Times
# nb = 50
# [1] "Fri May 12 09:59:53 2017"
# [1] "Sat May 13 11:03:58 2017"

# nb = 75
# [1] "Sat May 13 11:03:58 2017"
# [1] "Sun May 14 12:38:41 2017"

# nb = 100
# [1] "Sun May 14 12:38:41 2017"
# [1] "Mon May 15 14:05:13 2017"

# nb = 150
# [1] "Tue May 16 09:41:07 2017"
# [1] "Wed May 17 12:18:11 2017"

# nb = 200
# [1] "Wed May 17 12:18:12 2017"
# [1] "Thu May 18 14:32:31 2017"

# nb = 300
# [1] "Wed May  3 13:55:06 2017"
# [1] "Thu May  4 17:29:27 2017"


# Aggregate Probes --------------------------------------------------------

# Aggregate probes across bootstrap replicates for most frequently chosen ones
# Count number of probes chosen in at least 10, 15, 19, 20 times out of 20 reps
group <- pluck(sehrdf, "iClustf")
freqs <- purrr::set_names(c(10, 15, 19, 20))
sizes <- c(50, 75, 100, 150, 200, 300)
base_path <- "main/code/Analysis_METABRIC_Expression_Trends/Bootstrap/TopProbes/"
base_name <- "AgeDependent_Bootstrap_TopProbes"

# Read in bootstrap statistics for each NB size and collect into list of lists
vals.50 <- list.files(pattern = "Bootstrap_Stats_NB050", recursive = TRUE) %>% 
  map(read_rds)
vals.75 <- list.files(pattern = "Bootstrap_Stats_NB075", recursive = TRUE) %>% 
  map(read_rds)
vals.100 <- list.files(pattern = "Bootstrap_Stats_NB100", recursive = TRUE) %>% 
  map(read_rds)
vals.150 <- list.files(pattern = "Bootstrap_Stats_NB150", recursive = TRUE) %>% 
  map(read_rds)
vals.200 <- list.files(pattern = "Bootstrap_Stats_NB200", recursive = TRUE) %>% 
  map(read_rds)
vals.300 <- list.files(pattern = "Bootstrap_Stats_NB300", recursive = TRUE) %>% 
  map(read_rds)

## FC > 1.25 (original bootstrap statistics objects)
fc <- 1.25
for (nb in sizes) {
  agg <- get(glue("vals.{nb}")) %>% 
    boot_top(group, freqs)
  fn <- glue("{base_path}{base_name}_NB{fnb}_FC{ffc}.rds",
             fnb = str_pad(nb, width = 3, pad = 0),
             ffc = str_decimal(fc))
  print(bind_rows(modify_depth(agg, 2, length)))
  write_rds(agg, path = fn)
}

## FC > 2 and FC > 4 (modify significance condition in objects first)
for (fc in c(2, 4)) {
  for (nb in sizes) {
    agg <- get(glue("vals.{nb}")) %>% 
      modify_depth(2, metabric_signif, var = "BHadj_and_AgeDependentp",
                   fc = fc) %>% 
      boot_top(group, freqs)
    fn <- glue("{base_path}{base_name}_NB{fnb}_FC{ffc}.rds",
               fnb = str_pad(nb, width = 3, pad = 0),
               ffc = str_decimal(fc))
    print(bind_rows(modify_depth(agg, 2, length)))
    write_rds(agg, path = fn)
  }
}
