## This script generates results for the aroutdf object, using FC = 1.25, 2, 4

# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse)
library(glue)
library(fs)
library(ageassn)

# Raw data inputs
base_path <- "main/data/Data_METABRIC_Expression_Trends/"
save_path <- "main/code/Analysis_METABRIC_Expression_Trends/aroutdf/"
filename <- "{save_path}{subdir}/AgeDependent_BHadj_and_FC{i}_{g}_aroutdf.csv"
annodf <- read_rds(glue("{base_path}Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09.rds"))
annodf$ProbeId <- annodf$Probe_id # Add alternate probe ID column name
ebdf <- read.delim("main/data/Data_JasonCarroll/ER_binding_gene_hg19_threshold_1e-5.txt")  # Raw file has row names, can't use read_tsv
sehrdf <- read_rds(glue("{base_path}sehrdf.rds"))
Dataset_r <- read_rds(glue("{base_path}Dataset_r.rds"))

# Add subgroup definitions for ER/HER2 and ER+/ER-, modify iClustf levels
sehrdf <- sehrdf %>% 
  mutate(
    BrCaEHf = case_when(
      ER.Expr == "+" & Her2.Expr == "+" ~ "ER+/HER2+",
      ER.Expr == "+" & Her2.Expr == "-" ~ "ER+/HER2-",
      ER.Expr == "-" & Her2.Expr == "+" ~ "ER-/HER2+",
      ER.Expr == "-" & Her2.Expr == "-" ~ "ER-/HER2-",
      TRUE ~ NA_character_
    ) %>% 
      factor(levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-")) %>% 
      fct_relabel(str_sign),
    BrCaEf = case_when(
      ER.Expr == "+" ~ "ER+",
      ER.Expr == "-" ~ "ER-",
      TRUE ~ NA_character_
    ) %>% 
      factor(levels = c("ER+", "ER-")) %>% 
      fct_relabel(str_sign),
    iClustf = paste0("iClust", str_pad(iClustf, width = 2, pad = "0"))
  )


# Make aroutdf ------------------------------------------------------------

# Create directories for aroutdf objects
subdirs <- set_names(c("AllCases", "ER_HER2", "ER", "intClust", "Pam50"))
dir_create(path(save_path, subdirs))

# Common arguments
args <- lst(filename, annodf, ebdf, Dataset_r)

# All cases
subdir <- subdirs["AllCases"]
invoke(.f = metabric_aroutdf, .x = args, data = sehrdf)

# ER/HER2
group <- "BrCaEHf"
subdir <- subdirs["ER_HER2"]
sehrdf %>%
  split(.[[group]]) %>%
  walk(~ invoke(.f = metabric_aroutdf, .x = args, data = ., group = group))

# ER+/ER-
group <- "BrCaEf"
subdir <- subdirs["ER"]
sehrdf %>%
  split(.[[group]]) %>%
  walk(~ invoke(.f = metabric_aroutdf, .x = args, data = ., group = group))

# PAM50
group <- "BrCaf"
subdir <- subdirs["Pam50"]
sehrdf %>%
  split(.[[group]]) %>%
  walk(~ invoke(.f = metabric_aroutdf, .x = args, data = ., group = group))

# IntClust
group <- "iClustf"
subdir <- subdirs["intClust"]
sehrdf %>%
  split(.[[group]]) %>%
  walk(~ invoke(.f = metabric_aroutdf, .x = args, data = ., group = group))


# Add ER binding gene lists -----------------------------------------------

# Add Tier 1/2 ER binding gene lists from NKI
ebdf_t1 <- read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_1_EntrezID_and_Symbol_ER_at_promotors.txt")
ebdf_t2 <- read_tsv("main/data/Data_StaceyJoostenNKI_ERbinding/Tier_2_EntrezID_and_Symbol_total_ER_binding_via_20kb.txt")

dir_ls(path = save_path, recursive = TRUE, type = "file") %>% 
  map(read_csv) %>% 
  map(metabric_ebdf, ebdf_t1, ebdf_t2) %>% 
  map(select, -matches("X1")) %>% 
  iwalk(write.csv, row.names = FALSE)
