# Load packages and read in clinical data
library(biostatUtil)
library(Hmisc)
library(plyr)
library(dplyr)
dat.expr <- read_rds("main/data/Data_TCGA_ProstateCa_RNASeq_Expression_Trends/data_RNA_Seq_v2_expression_median.rds")
dat <- read_tsv("main/data/Data_TCGA_ProstateCa_RNASeq_Expression_Trends/data_bcr_clinical_data_patient.txt", skip = 4)
dat <- dat %>% 
  mutate(tcgaid = paste(PATIENT_ID, "01", sep = "-"),
         AGE = as.numeric(AGE)) %>% 
  filter(!is.na(AGE) & !is.na(match(tcgaid, names(dat.expr))))

# Clinical Characteristics ------------------------------------------------

# Number of patients
length(dat$PATIENT_ID)

# Age (years): median and distribution
summary(dat$AGE)
dat$AGE %>% 
  cut2(c(min(.), seq(30, 80, 10), max(.))) %>% 
  table(Range = .) %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Gender
dat$GENDER %>% 
  factor(levels = c("MALE", "FEMALE")) %>% 
  table(useNA = "always", dnn = "Gender") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Lymphovascular Invasion
dat$LYMPH_NODES_EXAMINED %>%
  mapvalues(from = "[Discrepancy]", to = "[Not Available]") %>% 
  factor(levels = c("YES", "NO", "[Not Available]")) %>% 
  table(useNA = "always", dnn = "Grade") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))


# TNM Staging -------------------------------------------------------------

# Clinical T classification
dat$CLIN_T_STAGE %>% 
  mapvalues(from = c("T1a", "T1b", "T1c", "T2a", "T2b", "T2c", "T3a", "T3b", "[Not Available]"),
            to = c("T1", "T1", "T1", "T2", "T2", "T2", "T3", "T3", "TX")) %>% 
  table(useNA = "always", dnn = "Tumour") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Clinical N classification
dat$CLIN_N_STAGE %>% 
  table(useNA = "always", dnn = "Nodes")

# Clinical M classification
dat$CLIN_M_STAGE %>% 
  mapvalues(from = c("M1a", "M1b", "M1c", "[Not Available]"),
            to = c("M1", "M1", "M1", "MX")) %>%
  table(useNA = "always", dnn = "Metastasis") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Pathological T classification
dat$PATH_T_STAGE %>% 
  mapvalues(from = c("T2a", "T2b", "T2c", "T3a", "T3b", "[Discrepancy]", "[Not Available]"),
            to = c("T2", "T2", "T2", "T3", "T3", "TX", "TX")) %>%
  table(useNA = "always", dnn = "Tumour") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Pathological N classification
dat$PATH_N_STAGE %>% 
  mapvalues(from = "[Not Available]", to = "NX") %>% 
  table(useNA = "always", dnn = "Nodes") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Pathological M classification
dat$PATH_M_STAGE %>% 
  mapvalues(from = "[Not Applicable]", to = "MX") %>%
  table(useNA = "always", dnn = "Metastasis") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))
