
# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse)
library(here)
library(fs)
library(tmatrends)

# Directory path for all-in-one files
DrNm <- "MB09/Analysis/TMACharacteristics"
walk(here(DrNm, c("Plots", "Tables")), create_dir)

# File name
flNm <- "MB09_IHC_Trends"

# Subdirectory path for individual files
plDrNm <- paste0(DrNm, "/Plots/", flNm, "_Individual")
walk(here(plDrNm), create_dir)

# Load MB09 whole data
tmadf <- read.csv(
  file = here(dirname(DrNm), "BuildData/MB09_wholedata.CSV"),
  stringsAsFactors = FALSE,
  na.strings = " ",
  comment.char = ""
)

# Add biomarker definitions
tmadf <- define_mb09(tmadf)


# Trend Analysis ----------------------------------------------------------

# Variables of interest
cscvars <- matrix(c(
  "bcl2_pp_v1n", "Non binding",
  "birc5_pp_v1n", "Non binding",
  "cd163_c_v1n", "Non binding",
  "ck56_c2v1n", "Non binding",
  "ck5_c2v1n", "Non binding",
  "ckit_c_v1n", "Non binding",
  "cyclind1_pp_v1n", "ER binding",
  "ecad_pp_v1n", "Non binding",
  "egfr_tsipp_v1n", "Non binding",
  "er_pp_v1n", "ER binding",
  "er_dcc_v1n", "ER binding",
  "ezh2_pp_v1n", "Non binding",
  "foxa1_pp_v1n", "ER binding",
  "h3k27me3_pp_v1n", "Epigenetic",
  "her2_v1n", "Non binding",
  "inpp4b_c_v1.ppn", "Non binding",
  "ki67_pp_v1n", "Non binding",
  "nestin_tsipp_v1n", "Non binding",
  "p16_qsipp_v1n", "Non binding",
  "pr_c_v1n", "ER binding",
  "tp53_pp_v1n", "Non binding"
), ncol = 2, byrow = TRUE) %>% 
  array_branch(margin = 2) %>% 
  set_names(c("GPEC_name", "ER_binding"))

cscvars4 <- c(
  "er_pp_v1n",
  "ezh2_pp_v1n",
  "foxa1_pp_v1n",
  "h3k27me3_pp_v1n"
)

cscvarsEH <- matrix(c(
  "bcl2_pp_v1n", "Non binding",
  "birc5_pp_v1n", "Non binding",
  "cd163_c_v1n", "Non binding",
  "ck56_c2v1n", "Non binding",
  "ck5_c2v1n", "Non binding",
  "ckit_c_v1n", "Non binding",
  "cyclind1_pp_v1n", "ER binding",
  "ecad_pp_v1n", "Non binding",
  "egfr_tsipp_v1n", "Non binding",
  "ezh2_pp_v1n", "Non binding",
  "foxa1_pp_v1n", "ER binding",
  "h3k27me3_pp_v1n", "Epigenetic",
  "inpp4b_c_v1.ppn", "Non binding",
  "ki67_pp_v1n", "Non binding",
  "nestin_tsipp_v1n", "Non binding",
  "p16_qsipp_v1n", "Non binding",
  "pr_c_v1n", "ER binding",
  "tp53_pp_v1n", "Non binding"
), ncol = 2, byrow = TRUE) %>% 
  array_branch(margin = 2) %>% 
  set_names(c("GPEC_name", "ER_binding"))

# Var names with and without ER/HER2
var_all <- cscvars[c("GPEC_name", "ER_binding")]
var_sub <- cscvarsEH[c("GPEC_name", "ER_binding")]


# Smoother Fit Correlation ------------------------------------------------

# Correlation measure of proportion of assay range covered on average by cases
# on supsmu() smoother fit
mb09_corr <- list(tmadf) %>% 
  purrr::set_names("MB09_TMA_Training_Correlation") %>% 
  map(function(data) {
    cscvars$GPEC_name %>% 
      purrr::set_names() %>%
      purrr::map_dbl(supsmu.cor, data = data, x = "age_at_diagnosis")
  }) %>% 
  do.call(cbind.data.frame, .) %>% 
  rownames_to_column("Gene")

write_csv(mb09_corr, path(DrNm, "Tables", "MB09_TMA_Correlations.csv"))


cli::cat_line("# Start: ", date())
#dirname <- dirname %||% ""
subgroup <- "BrCaEHf"
purrr::walk(levels(tmadf[[subgroup]]), function(icinm) {
  cli::cat_line(icinm)
  set.seed(473)
  # purrr::pwalk(vars, ~{
  #   simsmufit(x = .x, xdf = data[which(data[[subgroup]] == 
  #                                        icinm), ], nsim = nsim, alpha = alpha, ERbindingstatus = .y, 
  #             subgroup = icinm, pdfOutFilep = TRUE, plotDirName = dirname, 
  #             pdfFileNamePrefix = paste0(fn_append(filename, 
  #                                                  alpha, nsim), "_shaded_", fn_subgroup(icinm)))
  # })
  bs_corr <- list(tmadf) %>% 
    set_names(c("MB09_TMA_WholeSeries_Correlation")) %>% 
    map(function(data) {
      var_sub$GPEC_name %>% 
        set_names() %>%
        map_dbl(supsmu.cor, data = data[which(data[[subgroup]] == 
                                                icinm), ], x = "age_at_diagnosis")
    }) %>% 
    do.call(cbind.data.frame, .) %>% 
    rownames_to_column("Gene")
  
  write_csv(bs_corr, path(DrNm, "Tables", 
                          paste0("MB09_TMA_", 
                                 tmatrends:::fn_subgroup(icinm), 
                                 "_Correlations.csv")))
})
cli::cat_line("#   End: ", date())



# 95% CI, 2k sim, all-in-one ----------------------------------------------

ihc_trends_all(tmadf, cscvars, nsim = 2000, alpha = 0.05, filename = flNm,
               dirname = DrNm)
# Start: Mon Jan 22 14:04:29 2018
#   End: Mon Jan 22 14:06:11 2018


# 95% CI, 2k sim, individual, All cases -----------------------------------
ihc_trends_all(tmadf, cscvars, nsim = 2000, alpha = 0.05, filename = flNm,
               dirname = plDrNm, individual = TRUE, seed = 743)
# Start: Mon Jan 22 14:06:51 2018
#   End: Mon Jan 22 14:08:40 2018


# 95% CI, 2k sim, individual, ER HER2 subgroups ---------------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 2000, alpha = 0.05,
                    filename = flNm, dirname = plDrNm, seed = 963)
# Start: Mon Jan 22 14:09:33 2018
#   End: Mon Jan 22 14:12:30 2018


# 95% CI, 2k sim, individual, ER subgroups --------------------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 2000, alpha = 0.05,
                    filename = flNm, dirname = plDrNm, subgroup = "BrCaEf",
                    seed = 963)
# Start: Tue Jan 30 14:24:33 2018
#   End: Tue Jan 30 14:26:48 2018


# 99% CI, 10k sim, all-in-one ---------------------------------------------

ihc_trends_all(tmadf, cscvars, nsim = 10000, alpha = 0.01, filename = flNm,
               dirname = DrNm)
# Start: Mon Jan 22 14:12:59 2018
#   End: Mon Jan 22 14:21:27 2018


# 99% CI, 10k sim, individual, All cases ----------------------------------

ihc_trends_all(tmadf, cscvars, nsim = 10000, alpha = 0.01, filename = flNm,
               dirname = plDrNm, individual = TRUE, seed = 743)
# Start: Mon Jan 22 14:21:27 2018
#   End: Mon Jan 22 14:30:06 2018


# 99% CI, 10k sim, individual, ER HER2 subgroups --------------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 10000, alpha = 0.01,
                    filename = flNm, dirname = plDrNm, seed = 963)
# Start: Mon Jan 22 14:30:06 2018
#   End: Mon Jan 22 14:44:46 2018


# 99% CI, 10k sim, individual, ER subgroups -------------------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 10000, alpha = 0.01,
                    filename = flNm, dirname = plDrNm, subgroup = "BrCaEf",
                    seed = 963)
# Start: Tue Jan 30 14:26:48 2018
#   End: Tue Jan 30 14:37:38 2018


# 99.9% CI, 40k sim, all-in-one -------------------------------------------

ihc_trends_all(tmadf, cscvars, nsim = 40000, alpha = 0.001, filename = flNm,
               dirname = DrNm)
# Start: Mon Jan 22 14:45:54 2018
#   End: Mon Jan 22 15:19:45 2018


# 99.9% CI, 40k sim, individual, All cases --------------------------------

ihc_trends_all(tmadf, cscvars, nsim = 40000, alpha = 0.001, filename = flNm,
               dirname = plDrNm, individual = TRUE, seed = 743)
# Start: Mon Jan 22 15:19:46 2018
#   End: Mon Jan 22 15:52:59 2018


# 99.9% CI, 40k sim, individual, ER HER2 subgroups ------------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 40000, alpha = 0.001,
                    filename = flNm, dirname = plDrNm, seed = 963)
# Start: Mon Jan 22 15:52:59 2018
#   End: Mon Jan 22 16:48:23 2018


# 99.9% CI, 40k sim, individual, ER subgroups -----------------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 40000, alpha = 0.001,
                    filename = flNm, dirname = plDrNm, subgroup = "BrCaEf",
                    seed = 963)
# Start: Tue Jan 30 14:37:38 2018
#   End: Tue Jan 30 15:24:29 2018


# 99.99% CI, 400k sim, all-in-one -----------------------------------------

ihc_trends_all(tmadf, cscvars, nsim = 400000, alpha = 0.0001, filename = flNm,
               dirname = DrNm)
# Start: Mon Jan 22 16:59:00 2018
#   End: Mon Jan 22 22:23:11 2018


# 99.99% CI, 400k sim, individual, All cases ------------------------------

ihc_trends_all(tmadf, cscvars, nsim = 400000, alpha = 0.0001, filename = flNm,
               dirname = plDrNm, individual = TRUE, seed = 743)
# Start: Mon Jan 22 22:23:11 2018
#   End: Tue Jan 23 03:47:47 2018


# 99.99% CI, 400k sim, individual, ER HER2 subgroups ----------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 400000, alpha = 0.0001,
                    filename = flNm, dirname = plDrNm, seed = 963)
# Start: Tue Jan 23 14:22:09 2018
#   End: Wed Jan 24 00:49:48 2018


# 99.99% CI, 400k sim, individual, ER subgroups ---------------------------
ihc_trends_subgroup(tmadf, cscvarsEH, nsim = 400000, alpha = 0.0001,
                    filename = flNm, dirname = plDrNm, subgroup = "BrCaEf",
                    seed = 963)
# Start: Tue Jan 30 15:24:29 2018
#   End: Tue Jan 30 22:42:22 2018


# Significance Summaries --------------------------------------------------

# Full data
pattern <- "MB09_IHC.*kSimCIs_shaded_v01" %>% 
  purrr::set_names(paste0(
    DrNm,
    "/Tables/",
    stringr::str_split(., "\\.")[[1]][1],
    "_Significance.csv"))

signif_tables <- purrr::map(pattern, sig_count, path = DrNm,
                            genes = cscvars$GPEC_name)
iwalk(signif_tables, write_csv)

# ER/HER2 status
pattern <- paste0(c("ERpHER2n", "ERpHER2p", "ERnHER2p", "ERnHER2n", "ERp", "ERn"), "_")

path <- rep(plDrNm, 6) %>%
  purrr::set_names(paste0(
    DrNm,
    "/Tables/",
    purrr::map2_chr(., pattern, ~ gsub("Individual", .y, basename(.x))),
    "Significance.csv"
  ))

signif_tables2 <- purrr::map2(path, pattern, sig_count,
                              genes = cscvarsEH$GPEC_name)
iwalk(signif_tables2, write_csv)


# Copy Files --------------------------------------------------------------

# Copy files to BrCa_Age_Associated GitHub repo
github_path <- "~/Documents/GitHub/BrCa_Age_Associated/main/code/Analysis_MB09_TMA_IHC_Trends"

# Correlation Table
files <- dir_ls(path = here(DrNm, "Tables"), regexp = "Correlations")
new_path <- path(github_path, "Tables", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# All-in-one plots
files <- dir_ls(path = here(DrNm), regexp = "kSimCIs_shaded_v01")
new_path <- path(github_path, "Plots", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# Individual plots
sub_dir <- "MB09_IHC_Trends_Individual"
files <- dir_ls(path = here(DrNm, "Plots", sub_dir))
new_path <- path(github_path, "Plots", sub_dir, path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# Significance tables
files <- dir_ls(path = here(DrNm, "Tables"), regexp = "Significance")
new_path <- path(github_path, "Tables", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)
