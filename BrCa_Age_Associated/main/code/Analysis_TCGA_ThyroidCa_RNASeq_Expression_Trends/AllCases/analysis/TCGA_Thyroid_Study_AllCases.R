# Load packages
library(tidyverse)
library(glue)
library(fs)
library(ageassn)  # devtools::install_github("BCCRCMO/ageassn")

# Set constants, argument list, load data
biosig_args <-
  list(list(0.01, 1.25), list(0.01, 2), list(0.01, 4), list(0.05, 1.25)) %>%
  map(set_names, c("fdr", "fc"))
organ <- "Thyroid"
sub <- "AllCases"
root.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/{sub}/")
tables.dir <- file.path(root.dir, "tables/")
tcga_load(organ)

# Munge data inputs and obtain analysis outputs
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, gender = "all")

# Age-related output data frame for main significance conditions
ar <- tcga_aroutdf(exandf, sedf, erbnms, tables.dir, fdr = 0.01, fc = 1.25,
                   subset = TRUE, overwrite = TRUE)

# Calculate aroutdfs for modified conditions
tcga_modify %>% 
  invoke_map(biosig_args, ar, sedf, tables.dir, subset = TRUE) %>% 
  invisible()

# Add NKI ER binding gene list variable
dir_ls(tables.dir) %>%
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

# Fisher Exact Tests
fdr <- unique(map_dbl(biosig_args, "fdr"))
# tcga_fet(tables.dir, fdr = fdr, nki = FALSE)
# tcga_fet(tables.dir, fdr = fdr, nki = TRUE)

# Linear modelling of age and expression
file.name <- glue("TCGA_{organ}Ca_AgeRelated_{sub}{nrow(sedf)}")
tcga_lm %>%
  invoke_map(biosig_args, root.dir, file.name, sedf, tcgaarnms, erbnms) %>% 
  invisible()

# Gather list of aroutdf objects and create suffix for volcano/manhattan plots
aroutdf_lst <- tcga_augment %>%
  invoke_map(biosig_args, root.dir, annodf, erbnms, erbnms_t1, erbnms_t2, arnms) %>%
  set_names(map_chr(
    biosig_args,
    ~ paste("aroutdf", str_decimal(.$fdr), str_decimal(.$fc), sep = "_")
  ))
plot.suffix1 <- glue("TCGA_{organ}Ca_{sub}_AgeDependent_and_BHadjSignificant")
plot.suffix2 <- glue("TCGA_{organ}Ca_{sub}_AgeDependent_and_BHadjSignificant_v02")

# Volcano Plots
tcga_volcano(aroutdf_lst, root.dir, plot.suffix1,
             title = glue("TCGA {organ}Ca {sub} Volcano Plot\nAge Dependent and BH significant"),
             thin = TRUE)

# Manhattan Plots
tcga_manhattan(aroutdf_lst, root.dir, plot.suffix1,
               title = glue("TCGA {organ}Ca {sub} Manhattan Plot"))

# Volcano Plots with NKI ER binding colours
tcga_volcano(aroutdf_lst, root.dir, plot.suffix2,
             title = glue("TCGA {organ}Ca {sub} Volcano Plot\nAge Dependent and BH significant"),
             nki = TRUE,
             thin = TRUE)
