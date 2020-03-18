# Load packages
library(tidyverse)
library(gtools)
library(glue)
library(fs)
library(ageassn)

# Constants
fdr <- c(0.01, 0.05)
#legend_position <- "topleft"
legend_position <- "bottomleft"
set.seed(3137) ## hclust thinning uses sample()


# Read in memoized data files
organ <- "Br"
root.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
aroutdf_all <- dir_ls(path(root.dir, "Memoized"), regexp = "All_cases") %>% 
  map(read_csv)
aroutdf_erher2 <- dir_ls(path(root.dir, "Memoized"), regexp = "HER2") %>% 
  map(read_csv)
aroutdf_intclust <- dir_ls(path(root.dir, "Memoized"), regexp = "iClust") %>% 
  map(read_csv) %>% 
  `[`(mixedsort(names(.)))

# Read in points to label files
label.dir <- "main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends/"
labels_all <- dir_ls(label.dir, regexp = "All_cases.*All1992Cases") %>% 
  map(read_csv)
labels_erher2 <- dir_ls(label.dir, regexp = "HER2.*All1992Cases") %>% 
  map(read_csv) %>% 
  `[`(mixedorder(gsub("(ER.*HER2.*_+?).*_(FDR)", "\\1\\2", names(.))))
labels_intclust <- dir_ls(label.dir, regexp = "iClust.*All1992Cases") %>% 
  map(read_csv) %>% 
  `[`(mixedorder(gsub("(iClust_.+?_).*_(FDR)", "\\1\\2", names(.))))

# Load BrCa TCGA data and process
tcga_load(organ)
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, icdf)

# Tidy gene lists
erbnms <- tidy_genes(ebdf, "hgnc_symbol")
erbnms_t1 <- tidy_genes(ebdf_t1, "Symbol")
erbnms_t2 <- tidy_genes(ebdf_t2, "Symbol")

# Function for augmenting the data with some variables for plotting
volcano_augment <- function(x, annodf, labels, fc = 1.25, fdr = 0.01) {
  dplyr::mutate(
    x,
    !!"PBHadj" := .data$Best_BHadj_pval,
    !!"PBHadj_log" := -log10(.data$PBHadj),
    !!"symbol_match" := match(.data$Entrez_Gene_Id, annodf$Original_Entrez),
    !!"symbol" := annodf$Gene_symbol[.data$symbol_match],
    !!"ERbinding" := .data$symbol %in% erbnms,
    !!"ERbinding1" := .data$symbol %in% erbnms_t1,
    !!"ERbinding2" := .data$symbol %in% erbnms_t2,
    !!"vcp" := abs(.data$Best_log2FC) > log2(fc) & .data$PBHadj < fdr,
    !!"binding_group" := dplyr::case_when(
      .data$vcp & .data$ERbinding ~ "ER binding",
      .data$vcp & !.data$ERbinding ~ "Non binding",
      TRUE ~ " "
    ),
    !!"binding_group_NKI" := dplyr::case_when(
      .data$vcp & .data$ERbinding1  & .data$ERbinding2 ~ "Tier 1",
      .data$vcp & !.data$ERbinding1 & .data$ERbinding2 ~ "Tier 2 Only",
      .data$vcp & !.data$ERbinding1 & !.data$ERbinding2 ~ "Non binding",
      TRUE ~ " "
    ),
    !!"evcp" := .data$Entrez_Gene_Id %in% labels$Entrez
  )
}

# Add variables
aroutdf_all <- list(aroutdf_all, labels_all,
                    fdr = rep(fdr, length(aroutdf_all) / 2)) %>% 
  pmap(volcano_augment, annodf = annodf)
aroutdf_erher2 <- list(aroutdf_erher2, labels_erher2,
                       fdr = rep(fdr, length(aroutdf_erher2) / 2)) %>% 
  pmap(volcano_augment, annodf = annodf)
aroutdf_intclust <- list(aroutdf_intclust, labels_intclust,
                         fdr = rep(fdr, length(aroutdf_intclust) / 2)) %>% 
  pmap(volcano_augment, annodf = annodf)

# Common arguments to invoke plotting function over
args <- list(
  x = "Best_log2FC",
  y = "PBHadj_log",
  group = "binding_group_NKI",
  label = "symbol",
  subset = "evcp",
  legend = c("Tier 1", "Tier 2 Only", "Non binding"),
  color = c("red", "orange", "blue"),
  thin = TRUE,
  legloc = legend_position
)

# Invoke plotting function on each of the objects, saving to different files
volcano_patterns <- c("BrCa" = "BrCa_AgeRelated_Volcano", "outdf.csv" = "v02.pdf")

aroutdf_all %>%
  list(data = ., 
       filename = path(root.dir, "Plots/Volcano/AllCases",
                       str_replace_all(basename(names(.)), volcano_patterns)),
       group = glue("TCGA {organ}Ca All cases")) %>% 
  pwalk(~ invoke(volcano_plot, args, data = ..1, filename = ..2, title = ..3))

aroutdf_erher2 %>%
  list(data = .,
       filename = path(root.dir, "Plots/Volcano/ER_HER2",
                       str_replace_all(basename(names(.)), volcano_patterns)),
       group = glue('TCGA {organ}Ca {rep(c("ER-/HER2-", "ER-/HER2+", "ER+/HER2-", "ER+/HER2+"), each = 2)}')) %>% 
  pwalk(~ invoke(volcano_plot, args, data = ..1, filename = ..2, title = ..3))

aroutdf_intclust %>%
  list(data = .,
       filename = path(root.dir, "Plots/Volcano/intClust",
                       str_replace_all(basename(names(.)), volcano_patterns)),
       group = glue('TCGA {organ}Ca {rep(paste("iClust", seq_len(10), sep = "_"), each = 2)}')) %>% 
  pwalk(~ invoke(volcano_plot, args, data = ..1, filename = ..2, title = ..3))

# Equal axes plots
volcano_patterns <- c("BrCa" = "BrCa_AgeRelated_Volcano", "outdf.csv" = "v03.pdf")
xlim <- max_xlim(aroutdf_erher2)
ylim <- max_ylim(aroutdf_erher2)
aroutdf_erher2 %>%
  list(data = .,
       filename = path(root.dir, "Plots/Volcano/ER_HER2",
                       str_replace_all(basename(names(.)), volcano_patterns)),
       group = glue('TCGA {organ}Ca {rep(c("ER-/HER2-", "ER-/HER2+", "ER+/HER2-", "ER+/HER2+"), each = 2)}')) %>% 
  pwalk(~ invoke(volcano_plot, args, data = ..1, filename = ..2, title = ..3,
                 xlim = xlim, ylim = ylim))

xlim <- max_xlim(aroutdf_intclust)
ylim <- max_ylim(aroutdf_intclust)
aroutdf_intclust %>%
  list(data = .,
       filename = path(root.dir, "Plots/Volcano/intClust",
                       str_replace_all(basename(names(.)), volcano_patterns)),
       group = glue('TCGA {organ}Ca {rep(paste("iClust", seq_len(10), sep = "_"), each = 2)}')) %>% 
  pwalk(~ invoke(volcano_plot, args, data = ..1, filename = ..2, title = ..3,
                 xlim = xlim, ylim = ylim))
