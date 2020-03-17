# Load packages and input data
library(tidyverse)
library(fs)
library(ageassn)

ebdf_t1 <- read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_1_EntrezID_and_Symbol_ER_at_promotors.txt")
ebdf_t2 <- read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_2_EntrezID_and_Symbol_total_ER_binding_via_20kb.txt")
TCGA_BRCA_mRNAex_annotation <- read_csv("main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends/TCGA_BRCA_mRNAex_annotation.csv")

memoized_path <- "main/code/Analysis_METABRIC_Expression_Trends/Memoized"
probes_all <- map(dir_ls(path = memoized_path, regexp = "All_cases.*rds"), readRDS)
probes_intclust <- map(dir_ls(path = memoized_path, regexp = "iClust.*rds"), readRDS)
probes_pam50 <- map(dir_ls(path = memoized_path, regexp = "Pam50.*rds"), readRDS)
probes_erher2 <- map(dir_ls(path = memoized_path, regexp = "HER.*rds"), readRDS)

# Extract unique gene labels
erbnms_t1 <- tidy_genes(ebdf_t1, "Entrez ID")
erbnms_t2 <- tidy_genes(ebdf_t2, "Entrez ID")

# Probes to keep from TCGA BrCa expression (18950)
probes <- TCGA_BRCA_mRNAex_annotation$MatchingIlluminaProbeId

# Generate contingency tables for METABRIC datasets
fet_args <- lst(probes, erbnms_t1, erbnms_t2)

fet_all <- map(probes_all, ~ invoke(metabric_fet, fet_args, data = .))
fet_intclust <- map(probes_intclust, ~ invoke(metabric_fet, fet_args, data = .))
fet_pam50 <- map(probes_pam50, ~ invoke(metabric_fet, fet_args, data = .))
fet_erher2 <- map(probes_erher2, ~ invoke(metabric_fet, fet_args, data = .))

# Create FET tables directory
fet_path <- gsub("Memoized", "FET_tables", memoized_path)
dir_create(fet_path)

# Change names to reflect save directory
patterns <- c("Memoized" = "FET_tables",
              "AllProbes" = "TCGAProbes",
              "lm_v01" = "FET_tables")
fet_all %>% 
  `names<-`(str_replace_all(names(.), patterns)) %>% 
  iwalk(saveRDS)
fet_intclust %>% 
  `names<-`(str_replace_all(names(.), patterns)) %>% 
  iwalk(saveRDS)
fet_pam50 %>% 
  `names<-`(str_replace_all(names(.), patterns)) %>% 
  iwalk(saveRDS)
fet_erher2 %>% 
  `names<-`(str_replace_all(names(.), patterns)) %>% 
  iwalk(saveRDS)
