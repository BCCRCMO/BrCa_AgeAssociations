# Load packages
library(tidyverse)
library(glue)
library(fs)
library(ageassn)

# Load data
fdr <- c(0.05, 0.01)
root.dir <- glue("main/code/Analysis_TCGA_BrCa_RNASeq_Expression_Trends/FET_tables/")
tcga_load("Br")

# Prepare data
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, icdf,
             gender = "female")


# ChIP-Seq ER Binding -----------------------------------------------------

# Create full data
tcga_aroutdf(exandf, sedf, erbnms, root.dir, fdr = fdr, suffix = "All")

# Create EHf subsetted data
EHfs <- levels(sedf$BrCaEHf)
walk(EHfs, ~ {
  sedfg <- filter(sedf, BrCaEHf == .)
  tcga_aroutdf(exandf, sedfg, erbnms, root.dir, fdr = fdr,
               suffix = str_sign(.))
})

# Create iCf subsetted data
iCfs <- sort(unique(sedf$IntClust))
walk(iCfs, ~ {
  sedfg <- filter(sedf, IntClust == .)
  tcga_aroutdf(exandf, sedfg, erbnms, root.dir, fdr = fdr,
               suffix = paste0("iClust", str_pad(., width = 2, pad = 0)))
})


# NKI ER Binding ----------------------------------------------------------

# Read in aroutdf objects and add NKI ER binding labels, then resave
dir_ls(root.dir) %>%
  map(read_csv) %>% 
  map(~ mutate(
    .,
    !!"ERbinding_NKI" := case_when(
      .data$Entrez_Gene_Id %in% erbnms_t1 & .data$Entrez_Gene_Id %in% erbnms_t2 ~ "Tier 1",
      !.data$Entrez_Gene_Id %in% erbnms_t1 & .data$Entrez_Gene_Id %in% erbnms_t2 ~ "Tier 2 Only",
      !.data$Entrez_Gene_Id %in% erbnms_t1 & !.data$Entrez_Gene_Id %in% erbnms_t2 ~ "Not ER binding",
    )
  )) %>% 
  iwalk(write_csv)
