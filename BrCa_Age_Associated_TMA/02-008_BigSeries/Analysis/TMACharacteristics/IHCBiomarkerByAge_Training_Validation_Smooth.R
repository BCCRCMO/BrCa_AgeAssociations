
# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse)
library(here)
library(fs)
library(tmatrends)

# Directory path for all-in-one files
DrNm <- "02-008_BigSeries/Analysis/TMACharacteristics"
walk(here(DrNm, c("Plots", "Tables")), create_dir)

# File names
tflNm <- "BigSeries_Training_IHC_Trends"
vflNm <- "BigSeries_Validation_IHC_Trends"
wflNm <- "BigSeries_WholeSeries_IHC_Trends"

# Subdirectory paths for individual files
tplDrNm <- paste0(DrNm, "/Plots/", tflNm, "_Individual")
vplDrNm <- paste0(DrNm, "/Plots/", vflNm, "_Individual")
wplDrNm <- paste0(DrNm, "/Plots/", wflNm, "_Individual")
walk(here(c(tplDrNm, vplDrNm, wplDrNm)), create_dir)

# Load training and validation BigSeries data
tmadf <- read_csv(
  file = here(dirname(DrNm), "BuildData/02-008_BigSeries_Training.CSV"),
  na = " ",
  trim_ws = FALSE
)
vmadf <- read_csv(
  file = here(dirname(DrNm), "BuildData/02-008_BigSeries_Validation.CSV"),
  na = " ",
  trim_ws = FALSE,
  col_types = list(immuno_stains = col_logical())
)

# Add biomarker definitions
tmadf <- define_bs(tmadf)
vmadf <- define_bs(vmadf)
wmadf <- rbind(tmadf, vmadf)


# Trends Analysis ---------------------------------------------------------

# Variables of interest
cscvars <- matrix(c(
  "bcl2_c_v2.ppnc", "Non binding", "BCL2", "Bcl2",
  "ca9_b0v123_v1n", "Non binding", "CA9", "CA9",
  "ck56_v1nc", "Non binding", "KRT5", "CK5",
  "cldn3_v1nc", "Non binding", "CLDN3", "CLDN3",
  "cryab_b0v12_v1.1n", "Non binding", "CRYAB", "cryab",
  "cycline_v1.1.ppn", "Non binding", "CCNE1", "Cyclin E1",
  "ecad_v1nc", "Non binding", "CDH1", "E-Cadherin",
  "egfr_v1nc", "Non binding", "EGFR", "egfr",
  "ephb4_v1.1nc", "Non binding", "EPHB4", "EPHB4",
  "er_v2nc", "ER binding", "ESR1", "ER-alpha",
  "er_dcc_v1n", "ER binding", "ESR1", "ER-alpha",
  "ezh2_v1.ppnc", "Non binding", "EZH2", "EZH2",
  "foxa1_v1.ppnc", "ER binding", "FOXA1", "FoxA1",
  "gata3_v1nc", "Non binding", "GATA3", "gata3",
  "h3k27me3_v1.1.ppnc", "Epigenetic", "h3k27me3", "h3k27me3",
  "her2_v2nc", "Non binding", "ERBB2", "HER2",
  "her3_v1.1nc", "Non binding", "ERBB3", "HER3",
  "her4_v1nc", "ER binding", "ERBB4", "HER4",
  "hsp27_v1nc", "Non binding", "HSPB1", "HSP27",
  "igf1r_v2.0.1.ppn", "ER binding", "IGF1R", "IGF1R",
  "igfbp2_v1nc", "Non binding", "IGFBP2", "IGFBP2",
  "inpp4b_c_v1.ppn", "Non binding", "INPP4B", "INPP4B",
  "ki67_c_v1.6nc", "Non binding", "MKI67", "Ki67",
  "kit_v1nc", "Non binding", "KIT", "c-Kit",
  "krt5_v1nc", "Non binding", "KRT5", "CK5",
  "ku7080_c_v1nc", "ER binding", "XRCC6", "Ku70/80",
  "mdm2_v1nc", "Non binding", "MDM2", "Mdm2",
  "nestin_v1nc", "Non binding", "NES", "Nestin",
  "p16_c1v1nc", "Non binding", "CDKN2A", "p16INK4A",
  "p27_v1nc", "Non binding", "CDKN1B", "p27KIP1",
  "p53_c_v1.1.1.pnnc", "Non binding", "TP53", "p53",
  "pcad_v1.1.ppnc", "ER binding", "CDH3", "P-Cadherin",
  "pgp_v1nc", "Non binding", "PGP", "PGP",
  "pipki_v1nc", "Non binding", "PIP5K1C", "PIP5K1C",
  "plau_c1v1", "Non binding", "PLAU", "PLAU",
  "podxl_v1nc", "Non binding", "PODXL", "PODXL",
  "pr_v3nc", "ER binding", "PGR", "PgR/pr",
  "psf_c_v1nc", "Non binding", "SFPQ", "PSF",
  "pten_v1nc","Non binding", "PTEN", "PTEN",
  "ret_v1nc", "ER binding", "RET", "Ret",
  "skp2_v1nc", "Non binding", "SKP2", "Skp2",
  "trim29_v1nc", "Non binding", "TRIM29", "TRIM29",
  "yb1_v1.2nc", "Non binding", "YBX1", "YB1"
), ncol = 4, byrow = TRUE) %>% 
  array_branch(margin = 2) %>% 
  set_names(c("GPEC_name", "ER_binding", "Gene_symbol", "IHC_name"))

# Remove ER and HER2
cscvarsEH <- cscvars %>%
  map(`[`, -match(c("er_v2nc", "her2_v2nc"), .$GPEC_name))

# Var names with and without ER/HER2
var_all <- cscvars[c("GPEC_name", "ER_binding")]
var_sub <- cscvarsEH[c("GPEC_name", "ER_binding")]


# Smoother Fit Correlation ------------------------------------------------

# Correlation measure of proportion of assay range covered on average by cases
# on supsmu() smoother fit
bs_corr <- list(tmadf, vmadf, wmadf) %>% 
  set_names(c("BigSeries_TMA_Training_Correlation",
                     "BigSeries_TMA_Validation_Correlation",
                     "BigSeries_TMA_WholeSeries_Correlation")) %>% 
  map(function(data) {
    var_all$GPEC_name %>% 
      set_names() %>%
      map_dbl(supsmu.cor, data = data, x = "age_at_diagnosis")
  }) %>% 
  do.call(cbind.data.frame, .) %>% 
  rownames_to_column("Gene")

write_csv(bs_corr, path(DrNm, "Tables", "BigSeries_TMA_Correlations.csv"))


cli::cat_line("# Start: ", date())
#dirname <- dirname %||% ""
subgroup <- "BrCaEHf"
purrr::walk(levels(wmadf[[subgroup]]), function(icinm) {
  cli::cat_line(icinm)
  set.seed(473)
  # purrr::pwalk(vars, ~{
  #   simsmufit(x = .x, xdf = data[which(data[[subgroup]] == 
  #                                        icinm), ], nsim = nsim, alpha = alpha, ERbindingstatus = .y, 
  #             subgroup = icinm, pdfOutFilep = TRUE, plotDirName = dirname, 
  #             pdfFileNamePrefix = paste0(fn_append(filename, 
  #                                                  alpha, nsim), "_shaded_", fn_subgroup(icinm)))
  # })
  bs_corr <- list(wmadf) %>% 
    set_names(c("BigSeries_TMA_WholeSeries_Correlation")) %>% 
    map(function(data) {
      var_sub$GPEC_name %>% 
        set_names() %>%
        map_dbl(supsmu.cor, data = data[which(data[[subgroup]] == 
                                                icinm), ], x = "age_at_diagnosis")
    }) %>% 
    do.call(cbind.data.frame, .) %>% 
    rownames_to_column("Gene")
  
  write_csv(bs_corr, path(DrNm, "Tables", 
                          paste0("BigSeries_TMA_", 
                                 tmatrends:::fn_subgroup(icinm), 
                                 "_Correlations.csv")))
})
cli::cat_line("#   End: ", date())

# 95% CI, 2k sim, all-in-one ----------------------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 2000, alpha = 0.05, filename = tflNm,
               dirname = DrNm)
# Start: Tue Jan 23 14:38:23 2018
#   End: Tue Jan 23 14:41:29 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 2000, alpha = 0.05, filename = vflNm,
               dirname = DrNm)
# Start: Tue Jan 23 14:41:29 2018
#   End: Tue Jan 23 14:44:44 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 2000, alpha = 0.05, filename = wflNm,
               dirname = DrNm)
# Start: Tue Jan 23 14:44:44 2018
#   End: Tue Jan 23 14:49:47 2018


# 95% CI, 2k sim, individual, All cases -----------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 2000, alpha = 0.05, filename = tflNm,
               dirname = tplDrNm, individual = TRUE)
# Start: Tue Jan 23 14:49:47 2018
#   End: Tue Jan 23 14:52:57 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 2000, alpha = 0.05, filename = vflNm,
               dirname = vplDrNm, individual = TRUE)
# Start: Tue Jan 23 14:52:57 2018
#   End: Tue Jan 23 14:56:19 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 2000, alpha = 0.05, filename = wflNm,
               dirname = wplDrNm, individual = TRUE)
# Start: Tue Jan 23 14:56:19 2018
#   End: Tue Jan 23 15:01:27 2018


# 95% CI, 2k sim, individual, ER HER2 subgroups ---------------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 2000, alpha = 0.05,
                    filename = tflNm, dirname = tplDrNm)
# Start: Tue Jan 23 15:01:27 2018
#   End: Tue Jan 23 15:09:20 2018

# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 2000, alpha = 0.05,
                    filename = vflNm, dirname = vplDrNm)
# Start: Tue Jan 23 15:09:20 2018
#   End: Tue Jan 23 15:16:52 2018

# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 2000, alpha = 0.05,
                    filename = wflNm, dirname = wplDrNm)
# Start: Tue Jan 23 15:16:52 2018
#   End: Tue Jan 23 15:26:15 2018


# 95% CI, 2k sim, individual, ER subgroups --------------------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 2000, alpha = 0.05,
                    filename = tflNm, dirname = tplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 13:57:12 2018
#   End: Mon Jan 29 14:02:35 2018


# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 2000, alpha = 0.05,
                    filename = vflNm, dirname = vplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 14:02:35 2018
#   End: Mon Jan 29 14:07:56 2018


# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 2000, alpha = 0.05,
                    filename = wflNm, dirname = wplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 14:07:56 2018
#   End: Mon Jan 29 14:15:18 2018


# 99% CI, 10k sim, all-in-one ---------------------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 10000, alpha = 0.01, filename = tflNm,
               dirname = DrNm)
# Start: Tue Jan 23 15:26:15 2018
#   End: Tue Jan 23 15:43:07 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 10000, alpha = 0.01, filename = vflNm,
               dirname = DrNm)
# Start: Tue Jan 23 15:43:07 2018
#   End: Tue Jan 23 15:58:36 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 10000, alpha = 0.01, filename = wflNm,
               dirname = DrNm)
# Start: Tue Jan 23 15:58:36 2018
#   End: Tue Jan 23 16:22:54 2018


# 99% CI, 10k sim, individual, All cases ----------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 10000, alpha = 0.01, filename = tflNm,
               dirname = tplDrNm, individual = TRUE)
# Start: Tue Jan 23 16:22:54 2018
#   End: Tue Jan 23 16:38:05 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 10000, alpha = 0.01, filename = vflNm,
               dirname = vplDrNm, individual = TRUE)
# Start: Tue Jan 23 16:38:05 2018
#   End: Tue Jan 23 16:53:59 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 10000, alpha = 0.01, filename = wflNm,
               dirname = wplDrNm, individual = TRUE)
# Start: Tue Jan 23 16:53:59 2018
#   End: Tue Jan 23 17:18:07 2018


# 99% CI, 10k sim, individual, ER HER2 subgroups --------------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 10000, alpha = 0.01,
                    filename = tflNm, dirname = tplDrNm)
# Start: Tue Jan 23 17:18:07 2018
#   End: Tue Jan 23 17:50:23 2018

# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 10000, alpha = 0.01,
                    filename = vflNm, dirname = vplDrNm)
# Start: Tue Jan 23 17:50:23 2018
#   End: Tue Jan 23 18:22:36 2018

# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 10000, alpha = 0.01,
                    filename = wflNm, dirname = wplDrNm)
# Start: Tue Jan 23 18:22:36 2018
#   End: Tue Jan 23 19:02:08 2018


# 99% CI, 10k sim, individual, ER subgroups -------------------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 10000, alpha = 0.01,
                    filename = tflNm, dirname = tplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 14:15:18 2018
#   End: Mon Jan 29 14:40:40 2018


# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 10000, alpha = 0.01,
                    filename = vflNm, dirname = vplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 14:40:40 2018
#   End: Mon Jan 29 15:05:07 2018


# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 10000, alpha = 0.01,
                    filename = wflNm, dirname = wplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 15:05:07 2018
#   End: Mon Jan 29 15:38:16 2018



# 99.9% CI, 40k sim, all-in-one -------------------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 40000, alpha = 0.001, filename = tflNm,
               dirname = DrNm)
# Start: Tue Jan 23 19:02:08 2018
#   End: Tue Jan 23 20:00:46 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 40000, alpha = 0.001, filename = vflNm,
               dirname = DrNm)
# Start: Tue Jan 23 20:00:46 2018
#   End: Tue Jan 23 21:00:55 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 40000, alpha = 0.001, filename = wflNm,
               dirname = DrNm)
# Start: Tue Jan 23 21:00:55 2018
#   End: Tue Jan 23 22:34:51 2018


# 99.9% CI, 40k sim, individual, All cases --------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 40000, alpha = 0.001, filename = tflNm,
               dirname = tplDrNm, individual = TRUE)
# Start: Tue Jan 23 22:34:51 2018
#   End: Tue Jan 23 23:34:03 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 40000, alpha = 0.001, filename = vflNm,
               dirname = vplDrNm, individual = TRUE)
# Start: Tue Jan 23 23:34:03 2018
#   End: Wed Jan 24 00:35:04 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 40000, alpha = 0.001, filename = wflNm,
               dirname = wplDrNm, individual = TRUE)
# Start: Wed Jan 24 00:35:04 2018
#   End: Wed Jan 24 02:05:28 2018


# 99.9% CI, 40k sim, individual, ER HER2 subgroups ------------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 40000, alpha = 0.001,
                    filename = tflNm, dirname = tplDrNm)
# Start: Wed Jan 24 02:05:28 2018
#   End: Wed Jan 24 04:08:02 2018

# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 40000, alpha = 0.001,
                    filename = vflNm, dirname = vplDrNm)
# Start: Wed Jan 24 04:08:02 2018
#   End: Wed Jan 24 06:10:23 2018

# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 40000, alpha = 0.001,
                    filename = wflNm, dirname = wplDrNm)
# Start: Wed Jan 24 06:10:24 2018
#   End: Wed Jan 24 08:40:52 2018


# 99.9% CI, 40k sim, individual, ER subgroups -----------------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 40000, alpha = 0.001,
                    filename = tflNm, dirname = tplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 15:38:16 2018
#   End: Mon Jan 29 17:13:35 2018


# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 40000, alpha = 0.001,
                    filename = vflNm, dirname = vplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 17:13:35 2018
#   End: Mon Jan 29 18:45:27 2018


# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 40000, alpha = 0.001,
                    filename = wflNm, dirname = wplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 18:45:27 2018
#   End: Mon Jan 29 20:53:55 2018


# 99.99% CI, 400k sim, all-in-one -----------------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 400000, alpha = 0.0001, filename = tflNm,
               dirname = DrNm)
# Start: Thu Jan 25 21:59:19 2018
#   End: Fri Jan 26 07:21:32 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 400000, alpha = 0.0001, filename = vflNm,
               dirname = DrNm)
# Start: Thu Jan 25 12:42:34 2018
#   End: Thu Jan 25 22:59:56 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 400000, alpha = 0.0001, filename = wflNm,
               dirname = DrNm)
# Start: Fri Jan 26 09:45:48 2018
#   End: Sat Jan 27 01:55:15 2018


# 99.99% CI, 400k sim, individual, All cases ------------------------------

# Training
ihc_trends_all(tmadf, var_all, nsim = 400000, alpha = 0.0001, filename = tflNm,
               dirname = tplDrNm, individual = TRUE)
# Start: Sat Jan 27 01:55:15 2018
#   End: Sat Jan 27 11:43:59 2018

# Validation
ihc_trends_all(vmadf, var_all, nsim = 400000, alpha = 0.0001, filename = vflNm,
               dirname = vplDrNm, individual = TRUE)
# Start: Sat Jan 27 11:43:59 2018
#   End: Sat Jan 27 21:44:27 2018

# Whole Series
ihc_trends_all(wmadf, var_all, nsim = 400000, alpha = 0.0001, filename = wflNm,
               dirname = wplDrNm, individual = TRUE)
# Start: Sat Jan 27 21:44:27 2018
#   End: Sun Jan 28 13:17:12 2018


# 99.99% CI, 400k sim, individual, ER HER2 subgroups ----------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 400000, alpha = 0.0001,
                    filename = tflNm, dirname = tplDrNm)
# Start: Fri Jan 26 09:46:04 2018
#   End: Sat Jan 27 07:18:10 2018

# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 400000, alpha = 0.0001,
                    filename = vflNm, dirname = vplDrNm)
# Start: Sat Jan 27 07:18:10 2018
#   End: Sun Jan 28 04:14:48 2018

# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 400000, alpha = 0.0001,
                    filename = wflNm, dirname = wplDrNm)
# Start: Sun Jan 28 04:14:48 2018
#   End: Mon Jan 29 05:19:05 2018


# 99.99% CI, 400k sim, individual, ER subgroups ---------------------------

# Training
ihc_trends_subgroup(tmadf, var_sub, nsim = 400000, alpha = 0.0001,
                    filename = tflNm, dirname = tplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 14:01:42 2018
#   End: Tue Jan 30 05:33:46 2018


# Validation
ihc_trends_subgroup(vmadf, var_sub, nsim = 400000, alpha = 0.0001,
                    filename = vflNm, dirname = vplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 14:02:03 2018
#   End: Tue Jan 30 05:27:58 2018


# Whole Series
ihc_trends_subgroup(wmadf, var_sub, nsim = 400000, alpha = 0.0001,
                    filename = wflNm, dirname = wplDrNm, subgroup = "BrCaEf")
# Start: Mon Jan 29 14:03:27 2018
#   End: Tue Jan 30 10:47:51 2018


# TNP and Luminalp Setup --------------------------------------------------

# TNP File names
tntflNm <- "BigSeries_TNP_Training_IHC_Trends"
tnvflNm <- "BigSeries_TNP_Validation_IHC_Trends"
tnwflNm <- "BigSeries_TNP_WholeSeries_IHC_Trends"

# TNP data
tntmadf <- filter(tmadf, BrCa4 == "TNP")
tnvmadf <- filter(vmadf, BrCa4 == "TNP")
tnwmadf <- filter(wmadf, BrCa4 == "TNP")

# Luminalp File names
lptflNm <- "BigSeries_Luminalp_Training_IHC_Trends"
lpvflNm <- "BigSeries_Luminalp_Validation_IHC_Trends"
lpwflNm <- "BigSeries_Luminalp_WholeSeries_IHC_Trends"

# Luminalp data
lptmadf <- filter(tmadf, BrCa4 == "Luminalp")
lpvmadf <- filter(vmadf, BrCa4 == "Luminalp")
lpwmadf <- filter(wmadf, BrCa4 == "Luminalp")


# TNP 95% CI, 2k sim, all-in-one ------------------------------------------

# Training
ihc_trends_all(tntmadf, var_all, nsim = 2000, alpha = 0.05, filename = tntflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:12:12 2018
#   End: Wed Jan 24 10:13:42 2018

# Validation
ihc_trends_all(tnvmadf, var_all, nsim = 2000, alpha = 0.05, filename = tnvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:13:42 2018
#   End: Wed Jan 24 10:15:12 2018

# Whole Series
ihc_trends_all(tnwmadf, var_all, nsim = 2000, alpha = 0.05, filename = tnwflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:15:12 2018
#   End: Wed Jan 24 10:16:57 2018


# TNP 99% CI, 10k sim, all-in-one -----------------------------------------

# Training
ihc_trends_all(tntmadf, var_all, nsim = 10000, alpha = 0.01, filename = tntflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:16:57 2018
#   End: Wed Jan 24 10:24:26 2018

# Validation
ihc_trends_all(tnvmadf, var_all, nsim = 10000, alpha = 0.01, filename = tnvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:24:26 2018
#   End: Wed Jan 24 10:32:04 2018

# Whole Series
ihc_trends_all(tnwmadf, var_all, nsim = 10000, alpha = 0.01, filename = tnwflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:32:04 2018
#   End: Wed Jan 24 10:41:14 2018


# TNP 99.9% CI, 40k sim, all-in-one ---------------------------------------

# Training
ihc_trends_all(tntmadf, var_all, nsim = 40000, alpha = 0.001, filename = tntflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:41:14 2018
#   End: Wed Jan 24 11:14:29 2018

# Validation
ihc_trends_all(tnvmadf, var_all, nsim = 40000, alpha = 0.001, filename = tnvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 11:14:29 2018
#   End: Wed Jan 24 11:46:31 2018

# Whole Series
ihc_trends_all(tnwmadf, var_all, nsim = 40000, alpha = 0.001, filename = tnwflNm,
               dirname = DrNm)
# Start: Wed Jan 24 11:46:31 2018
#   End: Wed Jan 24 12:23:04 2018


# TNP 99.99% CI, 400k sim, all-in-one -------------------------------------

# Training
ihc_trends_all(tntmadf, var_all, nsim = 400000, alpha = 0.0001, filename = tntflNm,
               dirname = DrNm)
# Start: Wed Jan 24 10:29:50 2018
#   End: Wed Jan 24 16:14:55 2018

# Validation
ihc_trends_all(tnvmadf, var_all, nsim = 400000, alpha = 0.0001, filename = tnvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 16:14:56 2018
#   End: Wed Jan 24 21:23:00 2018

# Whole Series
ihc_trends_all(tnwmadf, var_all, nsim = 400000, alpha = 0.0001, filename = tnwflNm,
               dirname = DrNm)
# Start: Wed Jan 24 21:23:00 2018
#   End: Thu Jan 25 03:21:58 2018


# Luminalp 95% CI, 2k sim, all-in-one ------------------------------------------

# Training
ihc_trends_all(lptmadf, var_all, nsim = 2000, alpha = 0.05, filename = lptflNm,
               dirname = DrNm)
# Start: Wed Jan 24 12:23:04 2018
#   End: Wed Jan 24 12:25:53 2018

# Validation
ihc_trends_all(lpvmadf, var_all, nsim = 2000, alpha = 0.05, filename = lpvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 12:25:53 2018
#   End: Wed Jan 24 12:29:01 2018

# Whole Series
ihc_trends_all(lpwmadf, var_all, nsim = 2000, alpha = 0.05, filename = lpwflNm,
               dirname = DrNm)
# Start: Wed Jan 24 12:29:01 2018
#   End: Wed Jan 24 12:33:29 2018


# Luminalp 99% CI, 10k sim, all-in-one -----------------------------------------

# Training
ihc_trends_all(lptmadf, var_all, nsim = 10000, alpha = 0.01, filename = lptflNm,
               dirname = DrNm)
# Start: Wed Jan 24 12:33:29 2018
#   End: Wed Jan 24 12:47:00 2018

# Validation
ihc_trends_all(lpvmadf, var_all, nsim = 10000, alpha = 0.01, filename = lpvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 12:47:00 2018
#   End: Wed Jan 24 13:00:40 2018

# Whole Series
ihc_trends_all(lpwmadf, var_all, nsim = 10000, alpha = 0.01, filename = lpwflNm,
               dirname = DrNm)
# Start: Wed Jan 24 13:00:40 2018
#   End: Wed Jan 24 13:18:59 2018


# Luminalp 99.9% CI, 40k sim, all-in-one ----------------------------------

# Training
ihc_trends_all(lptmadf, var_all, nsim = 40000, alpha = 0.001, filename = lptflNm,
               dirname = DrNm)
# Start: Wed Jan 24 13:18:59 2018
#   End: Wed Jan 24 14:08:19 2018

# Validation
ihc_trends_all(lpvmadf, var_all, nsim = 40000, alpha = 0.001, filename = lpvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 14:08:19 2018
#   End: Wed Jan 24 15:02:45 2018

# Whole Series
ihc_trends_all(lpwmadf, var_all, nsim = 40000, alpha = 0.001, filename = lpwflNm,
               dirname = DrNm)
# Start: Wed Jan 24 15:02:45 2018
#   End: Wed Jan 24 16:18:34 2018


# Luminalp 99.99% CI, 400k sim, all-in-one -------------------------------------

# Training
ihc_trends_all(lptmadf, var_all, nsim = 400000, alpha = 0.0001, filename = lptflNm,
               dirname = DrNm)
# Start: Thu Jan 25 03:21:58 2018
#   End: Thu Jan 25 11:18:18 2018

# Validation
ihc_trends_all(lpvmadf, var_all, nsim = 400000, alpha = 0.0001, filename = lpvflNm,
               dirname = DrNm)
# Start: Wed Jan 24 21:52:54 2018
#   End: Thu Jan 25 06:03:22 2018

# Whole Series
ihc_trends_all(lpwmadf, var_all, nsim = 400000, alpha = 0.0001, filename = lpwflNm,
               dirname = DrNm)
# Start: Thu Jan 25 09:56:38 2018
#   End: Thu Jan 25 21:59:19 2018


# Significance Summaries --------------------------------------------------

# Because of the lengthy computation time, the following summaries are
# constructed by inspection of text annotations in PDF figures

# Full data and TNP/Luminalp subsets
pattern <- c(
  "BigSeries_Training_IHC",
  "BigSeries_Validation_IHC",
  "BigSeries_WholeSeries_IHC",
  "BigSeries_TNP_Training_IHC",
  "BigSeries_TNP_Validation_IHC",
  "BigSeries_TNP_WholeSeries_IHC",
  "BigSeries_Luminalp_Training_IHC",
  "BigSeries_Luminalp_Validation_IHC",
  "BigSeries_Luminalp_WholeSeries_IHC"
) %>%
  purrr::set_names(paste0(DrNm, "/Tables/", ., "_Significance.csv"))

signif_tables <- purrr::map(pattern, sig_count, path = DrNm,
                            genes = var_all$GPEC_name)
iwalk(signif_tables, write_csv)


# ER/HER2 status
pattern <- paste0(rep(c("ERpHER2n", "ERpHER2p", "ERnHER2p", "ERnHER2n",
                             "ERp", "ERn"), 3), "_")

path <- rep(c(tplDrNm, vplDrNm, wplDrNm), each = 6) %>%
  purrr::set_names(paste0(
    DrNm,
    "/Tables/",
    purrr::map2_chr(., pattern, ~ gsub("Individual", .y, basename(.x))),
    "Significance.csv"
  ))

signif_tables2 <- purrr::map2(path, pattern, sig_count,
                              genes = var_sub$GPEC_name)
iwalk(signif_tables2, write_csv)


# Copy Files --------------------------------------------------------------

# Copy files to BrCa_Age_Associated GitHub repo
github_path <- "~/Documents/GitHub/BrCa_Age_Associated/main/code/Analysis_BigSeries_TMA_IHC_Trends"

# Correlation Table
files <- dir_ls(path = here(DrNm, "Tables"), regexp = "Correlations")
new_path <- path(github_path, "Tables", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# All-in-one plots
files <- dir_ls(path = here(DrNm), regexp = "(?<!NoMBex_)IHC.*shaded", perl = TRUE)
new_path <- path(github_path, "Plots", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# Individual plots
sub_dir <- "BigSeries_Training_IHC_Trends_Individual"
files <- dir_ls(path = here(DrNm, "Plots", sub_dir))
new_path <- path(github_path, "Plots", sub_dir, path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

sub_dir <- "BigSeries_Validation_IHC_Trends_Individual"
files <- dir_ls(path = here(DrNm, "Plots", sub_dir))
new_path <- path(github_path, "Plots", sub_dir, path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

sub_dir <- "BigSeries_WholeSeries_IHC_Trends_Individual"
files <- dir_ls(path = here(DrNm, "Plots", sub_dir))
new_path <- path(github_path, "Plots", sub_dir, path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)

# Significance tables
files <- dir_ls(path = here(DrNm, "Tables"), regexp = "Significance")
new_path <- path(github_path, "Tables", path_file(files))
file_copy(path = files, new_path = new_path, overwrite = TRUE)
