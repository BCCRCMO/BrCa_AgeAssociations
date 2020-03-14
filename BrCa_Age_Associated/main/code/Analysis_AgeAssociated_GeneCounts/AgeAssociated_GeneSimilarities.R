# Load packages
library(tidyverse)
library(glue)
library(ageassn)

# Comparisons on fdr and fc levels
comps <- list(list(0.05, 1.25),
              list(0.01, 1.25),
              list(0.01, 2),
              list(0.01, 4)) %>% 
  set_names(gsub("list\\(|\\)|,", "", .)) %>% 
  map(set_names, c("fdr", "fc"))
  
# All Cases
# Filter METABRIC by matching TCGA probes
TCGA_BRCA_mRNAex_annotation <- read_csv("main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends/TCGA_BRCA_mRNAex_annotation.csv")
METABRIC <- read_rds("main/code/Analysis_METABRIC_Expression_Trends/Memoized/AgeRelated_AllProbes_All_cases_FDR_0p01_All1992Cases_lm_v01.rds")
METABRIC <- METABRIC %>% 
  filter(.data$Probe_id %in% TCGA_BRCA_mRNAex_annotation$MatchingIlluminaProbeId)

# TCGA BrCa
code.dir <- "main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/"
BrCa_files <- paste0(code.dir, "Memoized/TCGA_BrCa_All_cases_FDR_{fdr}_outdf.csv")
TCGA_BrCa <- map(set_names(c(0.05, 0.01)), ~ {
  glue(BrCa_files, organ = "Br", fdr = str_decimal(.))
}) %>% 
  map(read_csv)

# TCGA Other
TCGA_Kidney <- comps %>% 
  map(~ tcga_path(
    root.dir = glue(code.dir, organ = "Kidney"),
    file.name = "BHadj_and_AgeDependent_aroutdf",
    fdr = .[["fdr"]],
    fc = .[["fc"]],
    sub.dir = "AllCases/tables/"
  )) %>% 
  map(read_csv)

TCGA_Lung <- comps %>% 
  map(~ tcga_path(
    root.dir = glue(code.dir, organ = "Lung"),
    file.name = "BHadj_and_AgeDependent_aroutdf",
    fdr = .[["fdr"]],
    fc = .[["fc"]],
    sub.dir = "AllCases/tables/"
  )) %>% 
  map(read_csv)

TCGA_Prostate <- comps %>% 
  map(~ tcga_path(
    root.dir = glue(code.dir, organ = "Prostate"),
    file.name = "BHadj_and_AgeDependent_aroutdf",
    fdr = .[["fdr"]],
    fc = .[["fc"]],
    sub.dir = "tables/"
  )) %>% 
  map(read_csv)

TCGA_Thyroid <- comps %>% 
  map(~ tcga_path(
    root.dir = glue(code.dir, organ = "Thyroid"),
    file.name = "BHadj_and_AgeDependent_aroutdf",
    fdr = .[["fdr"]],
    fc = .[["fc"]],
    sub.dir = "AllCases/tables/"
  )) %>% 
  map(read_csv)

# Store all TCGA data in a list
organs <- c("Br", "Kidney", "Lung", "Prostate", "Thyroid")
TCGA_list <-
  list(TCGA_BrCa, TCGA_Kidney, TCGA_Lung, TCGA_Prostate, TCGA_Thyroid) %>%
  set_names(organs)

# METABRIC vs TCGA
args <- map(comps, ~ splice(., metabric = METABRIC))
mvt_tabs <- organs %>%
  map(~ invoke_map(comp_metabric_tcga, args, TCGA_list[[.]], organ = .))

# Pairwise TCGA
organ_pairs <- combn(organs, 2) %>% array_branch(1)
pair_tabs <- organ_pairs %>%
  pmap(
    ~ invoke_map(
      comp_pair_tcga,
      comps,
      x = TCGA_list[[.x]],
      y = TCGA_list[[.y]],
      x_organ = .x,
      y_organ = .y
    )
  )

# FET p-values and adjusted p-values for multiple comparisons
pvals <- splice(mvt_tabs, pair_tabs) %>% 
  modify_depth(2, ~ .[["fisher.ts"]]$p.value) %>%
  unlist()
adj_pvals <- p.adjust.methods %>% 
  set_names() %>% 
  map_df(p.adjust, p = pvals)
