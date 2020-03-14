
# Setup -------------------------------------------------------------------

# Load packages and set constants
library(tidyverse)
library(glue)
library(tools)
library(mosaic)
FDR <- 0.01
FC <- c(1.25, 2, 4) %>% set_names(paste0("Abs(FC) > ", .))


# TCGA BrCa ---------------------------------------------------------------

# Initial munging
# Source data
cldf <- read_tsv("main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends/data_bcr_clinical_data_patient.txt", skip = 4,
                 col_types = cols(INITIAL_PATHOLOGIC_DX_YEAR = "c",
                                  AGE = "c",
                                  DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS = "c"))  # Missing data coded as characters
ealldf <- read_rds("main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends/data_RNA_Seq_v2_expression_median.rds")
icdf <- read_tsv("main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends/TCGAIntClust1100.txt")

# Many entries are 0
# Jitter the medians so regression model fits are not singular
# Add a small amount of noise
annoColNames <- c("Hugo_Symbol", "Entrez_Gene_Id")
annoColidxs <- match(annoColNames, names(ealldf))
dataColidxs <- setdiff(seq(ncol(ealldf)), annoColidxs)
neDataRows <- nrow(ealldf)
neDataCols <- length(dataColidxs)
set.seed(719)
jitterMat <- matrix(runif(neDataRows * neDataCols, min = 1e-06, max = 1e-02),
                    nrow = neDataRows, ncol = neDataCols)
ealldf[, dataColidxs] <- ealldf[, dataColidxs] + jitterMat

# Drop probesets with no variation
ealldfSD <- apply(log2(16.0 + ealldf[, -c(1, 2)]), 1, sd, na.rm = TRUE)
edf <- ealldf[ealldfSD > 1e-6, ]

# IntClust
icdf <- icdf %>% 
  filter(!(substring(ID, 14, 15) == "06")) %>% 
  mutate(IntClust = factor(paste("iClust", IntClust), paste("iClust", seq(10))),
         tcgaid = paste(substring(ID, 1, 12), "01", sep = "-"))

# Merge IntClust and keep only rows with matching tcgaid samples in edf
sedf <- cldf %>% 
  mutate(tcgaid = paste(PATIENT_ID, "01", sep = "-"),
         Age = as.numeric(AGE)) %>% 
  filter(!is.na(Age) & !is.na(match(tcgaid, names(edf))) & GENDER == "FEMALE") %>% 
  left_join(icdf, by = "tcgaid", all.x = TRUE, all.y = FALSE) %>% 
  mutate(BrCaEHf = mosaic::derivedFactor(
    `ER- HER2-` = ER_STATUS_BY_IHC == "Negative" & (IHC_HER2 == "Negative" | (IHC_HER2 == "Equivocal" & HER2_FISH_STATUS == "Negative")),
    `ER- HER2+` = ER_STATUS_BY_IHC == "Negative" & (IHC_HER2 == "Positive" | (IHC_HER2 == "Equivocal" & HER2_FISH_STATUS == "Positive")),
    `ER+ HER2-` = ER_STATUS_BY_IHC == "Positive" & (IHC_HER2 == "Negative" | (IHC_HER2 == "Equivocal" & HER2_FISH_STATUS == "Negative")),
    `ER+ HER2+` = ER_STATUS_BY_IHC == "Positive" & (IHC_HER2 == "Positive" | (IHC_HER2 == "Equivocal" & HER2_FISH_STATUS == "Positive")),
    .method = "first", .default = NA_character_))

# Root directory
tcga_dir <- "main/code/Analysis_TCGA_BrCa_RNASeq_Expression_Trends/Memoized"

# Case Counts
case_counts <- with(sedf, c(`All Cases` = length(tcgaid),
                            table(BrCaEHf),
                            table(IntClust))) %>%
  enframe(name = "Subtype", value = "Case Count")

# All Cases
dat <- read_csv("main/code/Analysis_TCGA_BrCa_RNASeq_Expression_Trends/Memoized/TCGA_BrCa_All_cases_FDR_0p01_outdf.csv")
Allcounts <- map_int(FC, ~ sum(((dat$Abs_FoldChange > .) & (dat$PBHadj < FDR) & 
                                 (abs(dat$LE60_log2FC - dat$GT60_log2FC) > log2(.)) ) ) )
Allprops <- map_dbl(FC, ~ with(dat, n_distinct(Entrez_Gene_Id[((Abs_FoldChange > .) & (PBHadj < FDR) & 
                                                 (abs(LE60_log2FC - GT60_log2FC) > log2(.)))]) /
                                 n_distinct(Entrez_Gene_Id)))

Allgnms <- map(FC, ~ with(dat, Hugo_Symbol[((Abs_FoldChange > .) & 
                                                     (PBHadj < FDR) & 
                                                     (abs(LE60_log2FC - GT60_log2FC) > log2(.)))]))


# ER status
ERsubs <- list.files(path = tcga_dir, pattern = glue("ER.*{FDR}"),
                     full.names = TRUE) %>% 
  map(read_csv)
ERcounts <- ERsubs %>% 
  map_df(~ map(FC, function(fc) sum(( (.$Abs_FoldChange > fc) & (.$PBHadj < FDR) & 
                                      (abs(.$LE60_log2FC - .$GT60_log2FC) > log2(fc)) ) ) ) )
ERprops <- ERsubs %>% 
  map_df(~ map(FC, function(fc) with(., n_distinct(Entrez_Gene_Id[((Abs_FoldChange > fc) & (PBHadj < FDR) & 
                                                     (abs(LE60_log2FC - GT60_log2FC) > log2(fc)))]) /
                                       n_distinct(Entrez_Gene_Id)) ) )

# IntClust
ICsubs <- list.files(path = tcga_dir, pattern = glue("iClust.*{FDR}"),
                     full.names = TRUE) %>% 
  gtools::mixedsort() %>% 
  map(read_csv)
ICcounts <- ICsubs %>% 
  map_df(~ map(FC, 
               function(fc) sum( ( (.$Abs_FoldChange > fc) & 
                                     (.$PBHadj < FDR) & 
                                     (abs(.$LE60_log2FC - .$GT60_log2FC) > log2(fc)) ) ) ) )
ICprops <- ICsubs %>% 
  map_df(~ map(FC, 
               function(fc) with(., n_distinct(Entrez_Gene_Id[
                 ((Abs_FoldChange > fc) & 
                    (PBHadj < FDR) & 
                    (abs(LE60_log2FC - GT60_log2FC) > log2(fc)))]) /
                   n_distinct(Entrez_Gene_Id) ) ) )

# Result
TCGA_BrCa <- tibble(Source = "TCGA",
                    Site = "BrCa",
                    Gender = "Female") %>%
  cbind(
    case_counts,
    `Gene Count` = n_distinct(na.omit(dat$Hugo_Symbol)),
    `Probe Count` = nrow(dat),
    FDR = FDR,
    rbind(Allcounts, ERcounts, ICcounts),
    rbind(Allprops, ERprops, ICprops) %>%
      magrittr::set_colnames(paste(names(.), "(Prop)"))
  )


# METABRIC ----------------------------------------------------------------

# Case Counts
sehdf <- read_rds("main/data/Data_METABRIC_Expression_Trends/sehdf.rds")

# All cases
all_cc <- c("All Cases", nrow(sehdf))

# ER/HER2 status
er_her2_cc <- sehdf %>% 
  count(ER.Expr, Her2.Expr) %>% 
  transmute(Subtype = paste0("ER", ER.Expr, " HER2", Her2.Expr),
            `Case Count` = n)

# IntClust
intclust_cc <- sehdf %>% 
  count(iClustf) %>% 
  transmute(Subtype = paste("iClust", iClustf),
            `Case Count` = n)

# Pam50
pam50_cc <- sehdf %>%
  count(Pam50Subtype) %>% 
  transmute(Subtype = ifelse(Pam50Subtype == "Her2", "HER2", Pam50Subtype),
            `Case Count` = n) %>% 
  magrittr::extract(match(c("LumA", "LumB", "HER2", "Basal", "Normal"),
                          .$Subtype), )
# SM Q- Why "missing argument to function call" above in magrittr::extract??

# Combine case counts
case_counts <- rbind(all_cc, er_her2_cc, intclust_cc, pam50_cc)

# Target Counts

# If memoized data still in csv, compress into rds and then remove
metabric_dir <- "main/code/Analysis_METABRIC_Expression_Trends/Memoized"
metabric_files <- list_files_with_exts(metabric_dir, "csv")

# Matching Illumina ProbeIds from METABRIC to TCGA mRNA	
TCGA_BRCA_mRNAex_annotation <- read_csv("main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends/TCGA_BRCA_mRNAex_annotation.csv")

if (length(metabric_files)) {
  # Read in csv files, but set their names to have rds extensions, then compress
  metabric_csv <- metabric_files %>% 
    set_names(gsub("csv$", "rds", .)) %>% 
    map(read_csv)
  iwalk(metabric_csv, write_rds, compress = "xz")
  
  # Read in rds objects
  metabric_rds <- list_files_with_exts(metabric_dir, "rds") %>% map(read_rds)
  
  # Verify that the rds objects are identical to the csv before removing csv
  if (all(map2_lgl(metabric_csv, metabric_rds, identical))) {
    walk(metabric_files, file.remove)
  }
}

# Start tabulations
metabric_files <- list_files_with_exts(metabric_dir, "rds") %>% 
  purrr::set_names() %>% # ensure purrr namespace when no arguments given
  gtools::mixedsort() %>% # Move IntClust10 from between 1 and 2
  magrittr::extract(c(1:15, 18:19, 17:16, 20)) %>% # Specific PAM50 order
  map(read_rds) %>% 
  map(~ add_column(
    ., Hugo_Symbol = dat$Hugo_Symbol[match(.$Entrez_Gene_ID_0, dat$Entrez_Gene_Id)],
    .before = "Gene_symbol"
  ))

metabric_tabs <- metabric_files %>% 
  map_df(~ map(FC, function(fc) sum(.$Abs_FoldChange > fc &
                                      .$Best_BHadj_pval < FDR &
                                      (abs(.$LE60_log2FC - .$GT60_log2FC) > log2(fc)) ) ) )

metabric_props <- metabric_files %>% 
  map_df(~ map(FC, function(fc)
    with(., n_distinct(Entrez[((Abs_FoldChange > fc) & (Best_BHadj_pval < FDR) & 
                         (abs(LE60_log2FC - GT60_log2FC) > log2(fc)) )]) /
           n_distinct(Entrez) ) ) )

metabric_gnms <- metabric_files %>% 
  map(~ map(FC, function(fc)
    with(., Gene_symbol[((Abs_FoldChange > fc) & (Best_BHadj_pval < FDR) & 
                                 (abs(LE60_log2FC - GT60_log2FC) > log2(fc)) )]) ) )


metabric_gnms <- metabric_files %>% 
  map(~ map(FC, function(fc)
    with(., Hugo_Symbol[((Abs_FoldChange > fc) & (Best_BHadj_pval < FDR) & 
                           (abs(LE60_log2FC - GT60_log2FC) > log2(fc)) )]) ) )

# Result
METABRIC <- tibble(Source = "METABRIC",
                   Site = "BrCa",
                   Gender = "Female") %>%
  cbind(
    case_counts,
    `Gene Count` = n_distinct(na.omit(metabric_files[[1]]$Hugo_Symbol)),
    `Probe Count` = nrow(metabric_files[[1]]),
    FDR = FDR,
    metabric_tabs,
    metabric_props %>% 
      magrittr::set_colnames(paste(names(.), "(Prop)"))
  )

# Filtered for TCGA mRNA expression Result  ## SM Q - What is this for?
metabric_files_filtered <- metabric_files %>% 
  map(~ dplyr::filter(., Probe_id %in% TCGA_BRCA_mRNAex_annotation$MatchingIlluminaProbeId))
metabric_tabs_filtered  <- metabric_files_filtered %>% 
  map_df(~ map(FC, function(fc) sum(((.$Abs_FoldChange > fc) &
                                      (.$Best_BHadj_pval < FDR) &
                                       (abs(.$LE60_log2FC - .$GT60_log2FC) > log2(fc)) ) ) ) )

metabric_props_filtered  <- metabric_files_filtered %>% 
  map_df(~ map(FC, function(fc)
    with(., n_distinct(Entrez[((Abs_FoldChange > fc) & (Best_BHadj_pval < FDR) & 
                                 (abs(LE60_log2FC - GT60_log2FC) > log2(fc)))]) /
           n_distinct(Entrez) ) ) )


METABRIC_filtered <- tibble(Source = "METABRIC",
                            Site = "BrCa",
                            Gender = "Female") %>%
  cbind(
    case_counts,
    `Gene Count` = n_distinct(na.omit(metabric_files_filtered[[1]]$Hugo_Symbol)),
    `Probe Count` = nrow(metabric_files_filtered[[1]]),
    FDR = FDR,
    metabric_tabs_filtered,
    metabric_props_filtered %>% 
      magrittr::set_colnames(paste(names(.), "(Prop)"))
  )

# TCGA KidneyCa -----------------------------------------------------------

# Data Inputs
organ <- "Kidney"
code.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
data.dir <- glue("main/data/Data_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
np <- read_csv(glue("{code.dir}AllCases/tables/aroutdf_FDR_0p01_FC_1p25.csv"))
cldf <- read_tsv(glue("{data.dir}data_bcr_clinical_data_patient.txt"), skip = 4)
ealldf <- read_rds(glue("{data.dir}data_RNA_Seq_v2_expression_median.rds"))

# Case Counts
cc <- cldf %>% 
  transmute(tcgaid = paste(PATIENT_ID, "01", sep = "-"),
            Gender = str_to_title(GENDER),
            Subtype = "All Cases",
            Age = as.numeric(AGE)) %>% 
  filter(!is.na(Age) & !is.na(match(tcgaid, names(ealldf)))) %>% 
  count(Gender, Subtype) %>% 
  rbind(tibble(Gender = "All Cases",
               Subtype = "All Cases", 
               n = sum(.$n)), .) %>% 
  rename(`Case Count` = n)

# Target Counts
# paths <- c("AllCases", "Females", "Males") %>%
#   map(~ glue::glue("{code.dir}{.x}/tables")) %>%
#   map(list.files, pattern = "BHadj.*0p01", full.names = TRUE)
paths <- c("AllCases", "Females", "Males") %>%
  map(~ glue::glue("{code.dir}{.x}/tables")) %>%
  map(list.files, pattern = "^(aroutdf).*0p01", full.names = TRUE)
# 
# fnpre <- "main/code/Analysis_TCGA_KidneyCa_RNASeq_Expression_Trends/"
# fnmid <- "/tables/aroutdf_FDR_0p01_FC_"
# fnpost <- ".csv"
# ggrps <- c("AllCases", "Females", "Males")
# ggrpslbls <- c("All Cases", "Females", "Males")
# fclevs <- c(1.25, 2, 4)
# fclevsc <- gsub("\\.", "p", as.character(fclevs))
# 
# TCGA_KidneyCa <- TCGA_BrCa[-seq(nrow(TCGA_BrCa)), ]
# rowi <- 0
# for ( gi in seq(along = ggrps)) {
#   rowi <- rowi + 1
#   gci <- ggrps[gi]
#   gclbli <- ggrpslbls[gi]
#   TCGA_KidneyCa[rowi, "Source"] <- "TCGA"
#   TCGA_KidneyCa[rowi, "Site"] <- "KidneyCa"
#   TCGA_KidneyCa[rowi, "Gender"] <- gclbli
#   TCGA_KidneyCa[rowi, "Subtype"] <- "All Cases"  ## No biomarker subtypes here
#   TCGA_KidneyCa[rowi, "Case Count"] <- as.integer(cc[which(cc$Gender %in% gclbli), "Case Count"])
#   TCGA_KidneyCa[rowi, "FDR"] <- FDR
# 
#   for ( fi in seq(along = fclevsc)) {
#     fcci <- fclevsc[fi]
#     fci <- fclevs[fi]
#     fnm <- paste(fnpre, gci, fnmid, fcci, fnpost, sep = "")
#     cat("\n",fnm,"\n")
#     adfi <- read.csv(file = fnm)
#     TCGA_KidneyCa[rowi, "Gene Count"] <- nrow(adfi)
#     TCGA_KidneyCa[rowi, "Probe Count"] <- nrow(adfi)
# 
# 
#   }
# }

# Hugo_Symbol,Entrez_Gene_Id,LE60_slope,GT60_slope,AllAges_slope,LE60_pval,GT60_pval,AllAges_pval,
# AgeDependentp,LE60_BHadj_pval,GT60_BHadj_pval,AllAges_BHadj_pval,BHadj_and_AgeDependentp,
# BHadj_signifp,LE60_log2FC,GT60_log2FC,AllAges_log2FC,Best_log2FC,Best_BHadj_pval,ERbinding,ERbinding_NKI
# paths %>%
#   map( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .))
#   }) %>% map(names)
# 
# paths %>%
#   map( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .))
#   }) %>% map(names) %>%  unlist
#   map(names, ~ {as.numeric(gsub(., "Abs(FC) > ", ""))})
# # I need the FC from the file name to use in the map
# ardfs <-
#   paths %>%
#      map_df( ~ {
#        .x %>%
#          set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                      gsub("p", ".", .) %>%
#                      paste("Abs(FC) >", .)) %>%
#          map(read_csv) %>% 
#          map_df(~ map(FC, function(fc) sum(( (abs(.$Best_log2FC) > log2(fc)) & (.$Best_BHadj_pval < FDR) & 
#                                                (abs(.$LE60_log2FC - .$GT60_log2FC) > log2(fc)) ) ) ) ) # %>%
#                #map(nrow)
#      }) 
#vignette("programming", "dplyr")
# target counts of biphasic genes
tc <-
  paths %>%
  set_names(c("AllCases", "Females", "Males")) %>% 
  map_dfr(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~ sum(( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                                  (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ))    
  })  %>% 
  set_names( ~ paste("Abs(FC) >", .)) %>%
  cbind((. / nrow(np)) %>%
          magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
  mutate(Gender = c("All Cases", "Female", "Male"))

# Gene names of biphasic genes
tcgnms <-
  paths %>%
  set_names(c("AllCases", "Females", "Males")) %>% 
  map(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~ np$Hugo_Symbol[ ( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                     (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ]  )    
  }) 


# %>% 
#   mutate(Gender = c("All Cases", "Female", "Male"))

# %>% 
  # do.call(rbind, .)  ## To get a row named data frame

# tc <- paths %>%
#   map_df( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .)) %>%
#       map(read_csv) %>%
#       map(nrow)
#   }) %>%
#   cbind((. / nrow(np)) %>%
#           magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
#   mutate(Gender = c("All Cases", "Female", "Male"))

# 
# paths %>% map( ~ { .x %>% print } )
# tc <- paths %>% map_df( ~ { .x %>% set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                                            gsub("p", ".", .) %>%
#                                            paste("Abs(FC) >", .)) } ) %>%
#       map(gsub(".*FC_(.*)\\.csv", "\\1", .)) %>%
#                                               map(gsub("p", ".", .)) } )
#   map_df( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .)) %>%
#       map(read_csv) %>% map(str) } )
#        map(with(., extract((abs(.$LE60_log2FC - .$GT60_log2FC) > log2( as.numeric(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                                                        gsub("p", ".", .) %>%
#                                                        paste("Abs(FC) >", .)))), ) ) ) %>%
#       map(nrow)
#   }) %>%
#   cbind((. / nrow(np)) %>%
#           magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
#   mutate(Gender = c("All Cases", "Female", "Male"))
# 
# paths <- c("AllCases", "Females", "Males") %>%
#   map(~ glue::glue("{code.dir}{.x}/tables")) %>%
#   map(list.files, pattern = "BHadj.*0p01", full.names = TRUE)
# tc <- paths %>%
#   map_df( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .)) %>%
#       map(read_csv) %>%
#       map(filter) %>% ## How to spec the filter ???
#       map(nrow)
#   }) %>%
#   cbind((. / nrow(np)) %>%
#           magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
#   mutate(Gender = c("All Cases", "Female", "Male"))
# 
# paths <- c("AllCases", "Females", "Males") %>%
#   map(~ glue::glue("{code.dir}{.x}/tables")) %>%
#   map(list.files, pattern = "BHadj.*0p01", full.names = TRUE)
# 
# tc <- paths %>%
#   map_df( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .)) %>%
#       map(read_csv) %>%
#       map(nrow)
#   }) %>%
#   cbind((. / nrow(np)) %>%
#           magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
#   mutate(Gender = c("All Cases", "Female", "Male"))

# Result
TCGA_KidneyCa <- tibble(Source = "TCGA",
                        Site = "KidneyCa") %>%
  cbind(cc, `Gene Count` = n_distinct(na.omit(np$Hugo_Symbol)),
        `Probe Count` = nrow(np), FDR = FDR) %>%
  left_join(tc, by = "Gender")


# TCGA LungCa -----------------------------------------------------------

# Data Inputs
organ <- "Lung"
code.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
data.dir <- glue("main/data/Data_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
np <- read_csv(glue("{code.dir}AllCases/tables/aroutdf_FDR_0p01_FC_1p25.csv"))
cldf <- read_tsv(glue("{data.dir}data_bcr_clinical_data_patient.txt"), skip = 4)
ealldf <- read_rds(glue("{data.dir}data_RNA_Seq_v2_expression_median.rds"))

# Case Counts
cc <- cldf %>% 
  transmute(tcgaid = paste(PATIENT_ID, "01", sep = "-"),
            Gender = str_to_title(GENDER),
            Subtype = "All Cases",
            Age = as.numeric(AGE)) %>% 
  filter(!is.na(Age) & !is.na(match(tcgaid, names(ealldf)))) %>% 
  count(Gender, Subtype) %>% 
  rbind(tibble(Gender = "All Cases",
               Subtype = "All Cases", 
               n = sum(.$n)), .) %>% 
  rename(`Case Count` = n)

# Target Counts
# paths <- c("AllCases", "Females", "Males") %>%
#   map(~ glue::glue("{code.dir}{.x}/tables")) %>%
#   map(list.files, pattern = "BHadj.*0p01", full.names = TRUE)
# tc <- paths %>%
#   map_df( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .)) %>%
#       map(read_csv) %>%
#       map(nrow)
#   }) %>%
#   cbind((. / nrow(np)) %>%
#           magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
#   mutate(Gender = c("All Cases", "Female", "Male"))

paths <- c("AllCases", "Females", "Males") %>%
  map(~ glue::glue("{code.dir}{.x}/tables")) %>%
  map(list.files, pattern = "^(aroutdf).*0p01", full.names = TRUE)

tc <-
  paths %>%
  set_names(c("AllCases", "Females", "Males")) %>% 
  map_dfr(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~ sum(( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                     (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ))    
  })  %>% 
  set_names( ~ paste("Abs(FC) >", .)) %>%
  cbind((. / nrow(np)) %>%
          magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
  mutate(Gender = c("All Cases", "Female", "Male"))

# Biphasic gene names
tcgnms <-
  paths %>%
  set_names(c("AllCases", "Females", "Males")) %>% 
  map(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~ .x$Hugo_Symbol[ ( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                     (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ]  )    
  })  %>% 
  set_names( ~ paste("Abs(FC) >", .))


# Result
TCGA_LungCa <- tibble(Source = "TCGA",
                      Site = "LungCa") %>%
  cbind(cc, `Gene Count` = n_distinct(na.omit(np$Hugo_Symbol)),
        `Probe Count` = nrow(np), FDR = FDR) %>%
  left_join(tc, by = "Gender")


# TCGA ProstateCa -----------------------------------------------------------

# Data Inputs
organ <- "Prostate"
code.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
data.dir <- glue("main/data/Data_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
np <- read_csv(glue("{code.dir}tables/aroutdf_FDR_0p01_FC_1p25.csv"))
cldf <- read_tsv(glue("{data.dir}data_bcr_clinical_data_patient.txt"), skip = 4)
ealldf <- read_rds(glue("{data.dir}data_RNA_Seq_v2_expression_median.rds"))

# Case Counts
cc <- cldf %>% 
  transmute(tcgaid = paste(PATIENT_ID, "01", sep = "-"),
            Gender = str_to_title(GENDER),
            Subtype = "All Cases",
            Age = as.numeric(AGE)) %>% 
  filter(!is.na(Age) & !is.na(match(tcgaid, names(ealldf)))) %>% 
  count(Gender, Subtype) %>% 
  rename(`Case Count` = n)

# Target Counts
# paths <- c("") %>%
#   map(~ glue::glue("{code.dir}{.x}/tables")) %>%
#   map(list.files, pattern = "BHadj.*0p01", full.names = TRUE)
# tc <- paths %>%
#   map_df( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .)) %>%
#       map(read_csv) %>%
#       map(nrow)
#   }) %>%
#   cbind((. / nrow(np)) %>%
#           magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
#   mutate(Gender = c("Male"))


paths <- c("") %>%
  map(~ glue::glue("{code.dir}{.x}/tables")) %>%
  map(list.files, pattern = "^(aroutdf).*0p01", full.names = TRUE)

tc <-
  paths %>%
  set_names(c("")) %>% 
  map_dfr(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~ sum(( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                     (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ))    
  })  %>% 
  set_names( ~ paste("Abs(FC) >", .)) %>%
  cbind((. / nrow(np)) %>%
          magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
  mutate(Gender = c("Male"))


# Biphasic gene names
tcgnms <-
  paths %>%
  set_names(c("")) %>% 
  map(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~  .x$Hugo_Symbol[( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                     (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ] )    
  })  %>% 
  set_names( ~ paste("Abs(FC) >", .))


# Result
TCGA_ProstateCa <- tibble(Source = "TCGA",
                          Site = "ProstateCa") %>%
  cbind(cc, `Gene Count` = n_distinct(na.omit(np$Hugo_Symbol)),
        `Probe Count` = nrow(np), FDR = FDR) %>% 
  left_join(tc, by = "Gender")


# TCGA ThyroidCa -----------------------------------------------------------

# Data Inputs
organ <- "Thyroid"
code.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
data.dir <- glue("main/data/Data_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
np <- read_csv(glue("{code.dir}AllCases/tables/aroutdf_FDR_0p01_FC_1p25.csv"))
cldf <- read_tsv(glue("{data.dir}data_bcr_clinical_data_patient.txt"), skip = 4)
ealldf <- read_rds(glue("{data.dir}data_RNA_Seq_v2_expression_median.rds"))

# Case Counts
cc <- cldf %>% 
  transmute(tcgaid = paste(PATIENT_ID, "01", sep = "-"),
            Gender = str_to_title(GENDER),
            Subtype = "All Cases",
            Age = as.numeric(AGE)) %>% 
  filter(!is.na(Age) & !is.na(match(tcgaid, names(ealldf)))) %>% 
  count(Gender, Subtype) %>% 
  rbind(tibble(Gender = "All Cases",
               Subtype = "All Cases", 
               n = sum(.$n)), .) %>% 
  rename(`Case Count` = n)

# Target Counts
# paths <- c("AllCases", "Females", "Males") %>%
#   map(~ glue::glue("{code.dir}{.x}/tables")) %>%
#   map(list.files, pattern = "BHadj.*0p01", full.names = TRUE)
# tc <- paths %>%
#   map_df( ~ {
#     .x %>%
#       set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
#                   gsub("p", ".", .) %>%
#                   paste("Abs(FC) >", .)) %>%
#       map(read_csv) %>%
#       map(nrow)
#   }) %>%
#   cbind((. / nrow(np)) %>%
#           magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
#   mutate(Gender = c("All Cases", "Female", "Male"))


paths <- c("AllCases", "Females", "Males") %>%
  map(~ glue::glue("{code.dir}{.x}/tables")) %>%
  map(list.files, pattern = "^(aroutdf).*0p01", full.names = TRUE)

tc <-
  paths %>%
  set_names(c("AllCases", "Females", "Males")) %>% 
  map_dfr(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~ sum(( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                     (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ))    
  })  %>% 
  set_names( ~ paste("Abs(FC) >", .)) %>%
  cbind((. / nrow(np)) %>%
          magrittr::set_colnames(paste(names(.), "(Prop)"))) %>% 
  mutate(Gender = c("All Cases", "Female", "Male"))


# Biphasic gene names

tcgnms <-
  paths %>%
  set_names(c("AllCases", "Females", "Males")) %>% 
  map(~ {
    .x %>%
      set_names(gsub(".*FC_(.*)\\.csv", "\\1", .) %>%
                  gsub("p", ".", .) %>%
                  paste("Abs(FC) >", .)) %>%
      map(read_csv) %>% 
      set_names(gsub("Abs\\(FC\\) > ", "", names(.))) %>% 
      imap(~  .x$Hugo_Symbol[( (abs(.x$Best_log2FC) > log2(as.numeric(.y))) & (.x$Best_BHadj_pval < FDR) &
                     (abs(.x$LE60_log2FC - .x$GT60_log2FC) > log2(as.numeric(.y))) ) ] )    
  })  %>% 
  set_names( ~ paste("Abs(FC) >", .))



# Result
TCGA_ThyroidCa <- tibble(Source = "TCGA",
                         Site = "ThyroidCa") %>%
  cbind(cc, `Gene Count` = n_distinct(na.omit(np$Hugo_Symbol)),
        `Probe Count` = nrow(np), FDR = FDR) %>% 
  left_join(tc, by = "Gender")


# Final Result ------------------------------------------------------------

# Combine data sources
final <-
  rbind(TCGA_BrCa,
        METABRIC,
        TCGA_KidneyCa,
        TCGA_LungCa,
        TCGA_ProstateCa,
        TCGA_ThyroidCa) %>% tbl_df

# Add percentages
pcts <- final %>%
  select(matches("Prop")) %>%
  mutate_all(funs(case_when(
    . == 0 ~ "0%",
    . > 0 & . < 0.001 ~ "< 0.1%",
    TRUE ~ paste0(format(round(. * 100, digits = 1), nsmall = 1), "%")
  ))) %>%
  magrittr::set_names(gsub("Prop", "Pct", names(.)))

# Format proportions
final <- final %>%
  mutate_at(vars(matches("Prop")),
            funs(case_when(
              . == 0 ~ "0",
              . > 0 & . < 0.001 ~ "< 0.001",
              TRUE ~ format(round(., digits = 3))
            ))) %>%
  cbind(pcts)

# Final version for paper: collapse/remove M/F, remove probe & keep gene column
# Main: overlapping gene loci, Supplemental: full probes
final_main <- rbind(TCGA_BrCa, METABRIC_filtered, TCGA_KidneyCa,
                    TCGA_LungCa, TCGA_ProstateCa, TCGA_ThyroidCa)

# Add percentages
pcts_main <- final_main %>%
  select(matches("Prop")) %>%
  mutate_all(funs(case_when(
    . == 0 ~ "0%",
    . > 0 & . < 0.001 ~ "< 0.1%",
    TRUE ~ paste0(format(round(. * 100, digits = 1), nsmall = 1), "%")
  ))) %>%
  magrittr::set_names(gsub("Prop", "Pct", names(.)))

final_main <- final_main %>% 
  mutate_at(vars(matches("Prop")),
            funs(case_when(
              . == 0 ~ "0",
              . > 0 & . < 0.001 ~ "< 0.001",
              TRUE ~ format(round(., digits = 3))
            ))) %>%
  cbind(pcts_main) %>% 
  filter(Site %in% c("BrCa", "ProstateCa") |
           Site != "BrCa" & Gender == "All Cases"
  ) %>% 
  select(-Gender) %>% 
  select(-contains("Probe Count"))

  #tbl_df %>% 
#  select(-c("Probe Count", "Gender"))

final_supp <- final %>% 
  filter(Site %in% c("BrCa", "ProstateCa") |
           Site != "BrCa" & Gender == "All Cases"
  ) %>% 
  select(-Gender) %>% 
  select(-contains("Gene Count"))

#  select(-c("Gene Count", "Gender"))

write_csv(final, "main/code/Analysis_AgeAssociated_GeneCounts/AgeAssociated_BiphasicGeneCounts.csv")
write_csv(final_main, "main/code/Analysis_AgeAssociated_GeneCounts/AgeAssociated_BiphasicGeneCounts_Main.csv")
write_csv(final_supp, "main/code/Analysis_AgeAssociated_GeneCounts/AgeAssociated_BiphasicGeneCounts_Supp.csv")
