
# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse)
library(here)
library(fs)
library(tmatrends)

# Directory path for all-in-one files
DrNm <- "NormalBreast13/Analysis/PatientCharacteristics"
walk(here(DrNm, c("Plots", "Tables")), create_dir)

# File name
flNm <- "NormalBreast13_IHC_Trends"

# Subdirectory path for individual files
plDrNm <- paste0(DrNm, "/Plots/", flNm, "_Individual")
walk(here(plDrNm), create_dir)

# Load BigSeries whole data, NB13
wmadf <- read.csv(
  file = here("02-008_BigSeries/Analysis/BuildData/02-008_BigSeries_Whole.CSV"),
  stringsAsFactors = FALSE,
  na.strings = " ",
  comment.char = ""
)
pddf <- read.csv(
  file = here("NormalBreast13/NormalBreast13Creation/RandomizedLat4coresBefore2cores_NormalBreastTMA_NB13.csv"),
  stringsAsFactors = FALSE
)

# Load EZH2 and H3K27me3 tsv
score_dir <- "NormalBreast13/Scoring"
eddf <- read_tsv(here(score_dir, "EZH2/NB13_EZH2_v02.tsv"))
hddf <- read_tsv(here(score_dir, "H3K27me3/NB13_H3K27me3_v02.tsv"))

# Repeated columns
repeats <- c("SurgNo", "CoreNo", "NumCores")
hddf <- select(hddf, -one_of(repeats))

# Rest of the target scores in csv
targets <- c("foxa1", "esr1", "bcl2", "mki67", "pgr",
             "tp53", "ck5", "ecad", "egfr", "her2")
rddf <- list.files(here(score_dir), pattern = "^NB13.*csv$", recursive = TRUE,
                   full.names = TRUE) %>%
  `[`(map_int(targets, grep, x = ., ignore.case = TRUE)) %>%
  map(read.csv, stringsAsFactors = FALSE) %>%
  map(select, -one_of(repeats)) %>%
  map2(targets, ~ `names<-`(.x, paste(names(.x), .y, sep = "."))) %>%
  map(~ rename_at(., 3, funs(assign(., "nb13_id"))))

# Merge scores
addf <- c(list(pddf, eddf, hddf), rddf) %>%
  reduce(merge, by = "nb13_id", suffixes = c(".ezh2", ".h3k27me3"))

# Add biomarker definitions
# (Filter for randomly selected side for patients with bilateral samples)
sddf <- define_nb13(addf)


# Trend Analysis ----------------------------------------------------------

# Variables of Interest
cscvars <- matrix(c(
  "pp.bcl2", "Non binding",
  "pp.ck5", "Non binding",
  "pp.ecad", "Non binding",
  "int.egfr", "Non binding",
  "pp.esr1", "ER binding",
  "pp.ezh2", "Non binding",
  "pp.foxa1", "ER binding",
  "pp.h3k27me3", "Epigenetic",
  "pp.her2", "Non binding",
  "pp.mki67", "Non binding",
  "pp.pgr", "ER binding",
  "pp.tp53", "Non binding"
), ncol = 2, byrow = TRUE) %>% 
  array_branch(margin = 2) %>% 
  set_names(c("GPEC_name", "ER_binding"))


# Smoother Fit Correlation ------------------------------------------------

# Correlation measure of proportion of assay range covered on average by cases
# on supsmu() smoother fit
nb13_corr <- list(sddf) %>% 
  purrr::set_names("NormalBreast13_TMA_WholeSeries_Correlation") %>% 
  map(function(data) {
    cscvars$GPEC_name %>% 
      purrr::set_names() %>%
      purrr::map_dbl(supsmu.cor, data = data, x = "age_at_diagnosis")
  }) %>% 
  do.call(cbind.data.frame, .) %>% 
  rownames_to_column("Gene")

write_csv(nb13_corr, path(DrNm, "Tables", "NormalBreast13_TMA_Correlations.csv"))


# 95% CI, 2k sim, all-in-one ----------------------------------------------

ihc_trends_all(sddf, cscvars, nsim = 2000, alpha = 0.05, filename = flNm,
               dirname = DrNm)
# Start: Tue Jan 23 10:09:33 2018
#   End: Tue Jan 23 10:10:12 2018


# 95% CI, 2k sim, individual, All cases -----------------------------------

ihc_trends_all(sddf, cscvars, nsim = 2000, alpha = 0.05, filename = flNm,
               dirname = plDrNm, seed = 743, individual = TRUE)
# Start: Tue Jan 23 10:10:12 2018
#   End: Tue Jan 23 10:10:44 2018


# 99% CI, 10k sim, all-in-one ---------------------------------------------

ihc_trends_all(sddf, cscvars, nsim = 10000, alpha = 0.01, filename = flNm,
               dirname = DrNm)
# Start: Tue Jan 23 10:10:44 2018
#   End: Tue Jan 23 10:13:05 2018


# 99% CI, 10k sim, individual, All cases ----------------------------------

ihc_trends_all(sddf, cscvars, nsim = 10000, alpha = 0.01, filename = flNm,
               dirname = plDrNm, seed = 753, individual = TRUE)
# Start: Tue Jan 23 10:13:05 2018
#   End: Tue Jan 23 10:15:32 2018


# 99.9% CI, 40k sim, all-in-one -------------------------------------------

ihc_trends_all(sddf, cscvars, nsim = 40000, alpha = 0.001, filename = flNm,
               dirname = DrNm)
# Start: Tue Jan 23 10:15:32 2018
#   End: Tue Jan 23 10:25:11 2018


# 99.9% CI, 40k sim, individual, All cases --------------------------------

ihc_trends_all(sddf, cscvars, nsim = 40000, alpha = 0.001, filename = flNm,
               dirname = plDrNm, seed = 763, individual = TRUE)
# Start: Tue Jan 23 10:25:11 2018
#   End: Tue Jan 23 10:34:58 2018


# 99.99% CI, 400k sim, all-in-one -----------------------------------------

ihc_trends_all(sddf, cscvars, nsim = 400000, alpha = 0.0001, filename = flNm,
               dirname = DrNm)
# Start: Tue Jan 23 10:34:59 2018
#   End: Tue Jan 23 12:15:00 2018


# 99.99% CI, 400k sim, individual, All cases ------------------------------

ihc_trends_all(sddf, cscvars, nsim = 400000, alpha = 0.0001, filename = flNm,
               dirname = plDrNm, seed = 773, individual = TRUE)
# Start: Tue Jan 23 12:15:00 2018
#   End: Tue Jan 23 13:55:28 2018

# Significance Summaries --------------------------------------------------

pattern <- "NormalBreast13_IHC" %>% 
  purrr::set_names(paste0(DrNm, "/Tables/", ., "_Significance.csv"))

signif_tables <- purrr::map(pattern, sig_count, path = DrNm,
                            genes = cscvars$GPEC_name)

iwalk(signif_tables, write_csv)


# Copy Files --------------------------------------------------------------

# Copy files to BrCa_Age_Associated GitHub repo
github_path <- "~/Documents/GitHub/BrCa_Age_Associated/main/code/Analysis_NormalBreast13_TMA_IHC_Trends"

# Correlation Table
files <- dir_ls(path = here(DrNm, "Tables"), regexp = "Correlations")
new_path <- path(github_path, "Tables", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# All-in-one plots
files <- dir_ls(path = here(DrNm), regexp = "shaded")
new_path <- path(github_path, "Plots", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# Individual plots
sub_dir <- "NormalBreast13_IHC_Trends_Individual"
files <- dir_ls(path = here(DrNm, "Plots", sub_dir))
new_path <- path(github_path, "Plots", sub_dir, path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# Significance tables
files <- dir_ls(path = here(DrNm, "Tables"), regexp = "Significance")
new_path <- path(github_path, "Tables", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)
