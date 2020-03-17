### Age deciles
##aqmat <- floor(readRDS("main/code/Analysis_age_associated_IHC/AgeDeciles.rds"))

aqmat <- floor(readRDS("main/data/Data_PanData/aqmat.rds"))
tmaagedf <- readRDS("main/data/Data_PanData/tmaagedf.rds")


MBdf <- readRDS("main/data/Data_METABRIC_Expression_Trends/sardf.rds")
aqmat <- rbind(aqmat,
               quantile(floor(MBdf$age_at_diagnosis), probs = seq(0, 1, 0.10), na.rm = TRUE))
dimnames(aqmat)[[1]][nrow(aqmat)] <- "METABRIC"

tmaagedf <- rbind(tmaagedf,
                  data.frame(age_at_diagnosis = MBdf$age_at_diagnosis,
                             datasource = rep("METABRIC", length(MBdf$age_at_diagnosis)),
                             stringsAsFactors = FALSE))

ggplot(data=tmaagedf,aes(x=age_at_diagnosis, group=datasource, fill=datasource)) + 
  geom_density(adjust=1.5 , alpha=0.2)

# TCGA Breast

# Load --------------------------------------------------------------------
library(tidyverse)
library(ageassn)
library(ggplot2)
library(ggridges)

# Load data
data_root <- "main/data/Data_TCGA_BrCa_RNASeq_Expression_Trends"
tcga_load(organ = "Br")

# Constants
FDRalphalevels <- c(0.01, 0.05)
FCthresholds <- c(1.25, 2, 4)
root <- "main/code/Analysis_TCGA_BrCa_RNASeq_Expression_Trends/"


# Process -----------------------------------------------------------------

# Process data
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, icdf, gender = "female")

# Add IntClust variables
sedf$iClust <- paste0("iClust", str_pad(sedf$IntClust, width = 2, pad = 0))
sedf$iClustf <- factor(sedf$iClust)

TCGABrCasedf <- sedf
aqmat <- rbind(aqmat,
               quantile(as.numeric(TCGABrCasedf$AGE), probs = seq(0, 1, 0.10), na.rm = TRUE))
dimnames(aqmat)[[1]][nrow(aqmat)] <- "TCGA BrCa"

tmaagedf <- rbind(tmaagedf,
                  data.frame(age_at_diagnosis = as.numeric(TCGABrCasedf$AGE),
                             datasource = rep("TCGA BrCa", length(as.numeric(TCGABrCasedf$AGE))),
                             stringsAsFactors = FALSE))


# TCGA Kidney
# Load packages
library(tidyverse)
library(glue)
library(fs)
library(ageassn)  # devtools::install_github("BCCRCMO/ageassn")

# Set constants, argument list, load data
biosig_args <-
  list(list(0.01, 1.25), list(0.01, 2), list(0.01, 4), list(0.05, 1.25)) %>%
  map(set_names, c("fdr", "fc"))
organ <- "Kidney"
sub <- "AllCases"
root.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/{sub}/")
tables.dir <- file.path(root.dir, "tables/")
tcga_load(organ)

# Munge data inputs and obtain analysis outputs
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, gender = "all")

TCGAKidneyCasedf <- sedf
aqmat <- rbind(aqmat,
               floor(quantile(floor(as.numeric(TCGAKidneyCasedf$AGE)), 
                              probs = seq(0, 1, 0.10), na.rm = TRUE)))
dimnames(aqmat)[[1]][nrow(aqmat)] <- "TCGA KidneyCa"

tmaagedf <- rbind(tmaagedf,
                  data.frame(age_at_diagnosis = as.numeric(TCGAKidneyCasedf$AGE),
                             datasource = rep("TCGA KidneyCa", length(as.numeric(TCGAKidneyCasedf$AGE))),
                             stringsAsFactors = FALSE))


# TCGA Lung
# Load packages
library(tidyverse)
library(glue)
library(fs)
library(ageassn)  # devtools::install_github("BCCRCMO/ageassn")

# Set constants, argument list, load data
biosig_args <-
  list(list(0.01, 1.25), list(0.01, 2), list(0.01, 4), list(0.05, 1.25)) %>%
  map(set_names, c("fdr", "fc"))
organ <- "Lung"
sub <- "AllCases"
root.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/{sub}/")
tables.dir <- file.path(root.dir, "tables/")
tcga_load(organ)

# Munge data inputs and obtain analysis outputs
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, gender = "all")


TCGALungCasedf <- sedf
aqmat <- rbind(aqmat,
               floor(quantile(floor(as.numeric(TCGALungCasedf$AGE)), 
                              probs = seq(0, 1, 0.10), na.rm = TRUE)))
dimnames(aqmat)[[1]][nrow(aqmat)] <- "TCGA LungCa"

tmaagedf <- rbind(tmaagedf,
                  data.frame(age_at_diagnosis = as.numeric(TCGALungCasedf$AGE),
                             datasource = rep("TCGA LungCa", length(as.numeric(TCGALungCasedf$AGE))),
                             stringsAsFactors = FALSE))


# TCGA Prostate
# Load packages
library(tidyverse)
library(glue)
library(fs)
library(ageassn)  # devtools::install_github("BCCRCMO/ageassn")

# Set constants, argument list, load data
biosig_args <-
  list(list(0.01, 1.25), list(0.01, 2), list(0.01, 4), list(0.05, 1.25)) %>%
  map(set_names, c("fdr", "fc"))
organ <- "Prostate"
sub <- "AllCases"
root.dir <- glue("main/code/Analysis_TCGA_{organ}Ca_RNASeq_Expression_Trends/")
tables.dir <- file.path(root.dir, "tables/")
tcga_load(organ)

# Munge data inputs and obtain analysis outputs
tcga_process(cldf, ealldf, annodf, ardf, ezdf, ebdf, ebdf_t1, ebdf_t2, gender = "all")

TCGAProstateCasedf <- sedf
aqmat <- rbind(aqmat,
               floor(quantile(floor(as.numeric(TCGAProstateCasedf$AGE)), 
                              probs = seq(0, 1, 0.10), na.rm = TRUE)))
dimnames(aqmat)[[1]][nrow(aqmat)] <- "TCGA ProstateCa"

tmaagedf <- rbind(tmaagedf,
                  data.frame(age_at_diagnosis = as.numeric(TCGAProstateCasedf$AGE),
                             datasource = rep("TCGA ProstateCa", 
                                              length(as.numeric(TCGAProstateCasedf$AGE))),
                             stringsAsFactors = FALSE))



# TCGA Thyroid
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

TCGAThyroidCasedf <- sedf
aqmat <- rbind(aqmat,
               floor(quantile(floor(as.numeric(TCGAThyroidCasedf$AGE)), 
                              probs = seq(0, 1, 0.10), na.rm = TRUE)))
dimnames(aqmat)[[1]][nrow(aqmat)] <- "TCGA ThyroidCa"

tmaagedf <- rbind(tmaagedf,
                  data.frame(age_at_diagnosis = as.numeric(TCGAThyroidCasedf$AGE),
                             datasource = rep("TCGA ThyroidCa", 
                                              length(as.numeric(TCGAThyroidCasedf$AGE))),
                             stringsAsFactors = FALSE))
#[1] "BigSeries TMA"      "MB09 TMA"           "NormalBreast13 TMA" "METABRIC"           "TCGA BrCa"          "TCGA KidneyCa"      "TCGA LungCa"       
#[8] "TCGA ProstateCa"    "TCGA ThyroidCa"   
horder <- c(8,7,6,5,4,3,9,2,1)
#horder <- c(1:9)
unique(tmaagedf$datasource)[horder]
#> unique(tmaagedf$datasource)[horder]
#[1] "TCGA ProstateCa"    "TCGA LungCa"        "TCGA KidneyCa"      "TCGA BrCa"          "METABRIC"           "NormalBreast13 TMA" "TCGA ThyroidCa"    
#[8] "MB09 TMA"

tmaagedf$datasourcef <- factor(tmaagedf$datasource, levels = unique(tmaagedf$datasource)[horder])

ggplot(data=tmaagedf,aes(x=age_at_diagnosis, group=datasourcef, fill=datasourcef)) + 
  geom_density(adjust=1.5 , alpha=0.7) +
  xlab("Age at diagnosis (years)") + ylab("Density") +
  scale_fill_manual(name = "Data source", 
                     values = c("#FF5500", "#EE82EE",  "#00EE76", "#CD3278", "#8B0000",
                                "#0000CD", "#FFAA00", "#00C5CD", "#FFFF40")[horder])

tmaagedf$datasourcef <- factor(tmaagedf$datasource, levels = sort(unique(tmaagedf$datasource)))
ggplot(data=tmaagedf, aes(x = age_at_diagnosis, y = datasourcef)) + 
  geom_density_ridges() + 
  theme_minimal(base_size = 14) + theme(axis.text.y = element_text(vjust = 0)) +
  xlab("Age at diagnosis (years)") + ylab("Density") +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0))



## iClust pallette:  
## c("#FF5500", "#00EE76", "#CD3278", "#00C5CD", "#8B0000",
##   "#0000CD", "#FFAA00", "#EE82EE", "#FFFF40", "#7D26CD")

