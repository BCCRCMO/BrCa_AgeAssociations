# Load packages and read in clinical data
library(biostatUtil)
library(Hmisc)
library(plyr)
library(dplyr)
dat.expr <- read_rds("main/data/Data_TCGA_ThyroidCa_RNASeq_Expression_Trends/data_RNA_Seq_v2_expression_median.rds")
dat <- read_tsv("main/data/Data_TCGA_ThyroidCa_RNASeq_Expression_Trends/data_bcr_clinical_data_patient.txt", skip = 4)
dat <- dat %>% 
  mutate(tcgaid = paste(PATIENT_ID, "01", sep = "-"),
         AGE = as.numeric(AGE)) %>% 
  filter(!is.na(AGE) & !is.na(match(tcgaid, names(dat.expr))) & GENDER == "MALE")


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
  table(useNA = "always", dnn = "Tumour")

# Clinical N classification
dat$CLIN_N_STAGE %>% 
  table(useNA = "always", dnn = "Nodes")

# Clinical M classification
dat$CLIN_M_STAGE %>% 
  table(useNA = "always", dnn = "Metastasis")

# Pathological T classification
dat$AJCC_TUMOR_PATHOLOGIC_PT %>% 
  mapvalues(from = c("T1a", "T1b", "T4a"),
            to = c("T1", "T1", "T4")) %>%
  table(useNA = "always", dnn = "Tumour") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Pathological N classification
dat$AJCC_NODES_PATHOLOGIC_PN %>% 
  mapvalues(from = c("N1a", "N1b"), to = c("N1", "N1")) %>%
  table(useNA = "always", dnn = "Nodes") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))

# Pathological M classification
dat$AJCC_METASTASIS_PATHOLOGIC_PM %>% 
  table(useNA = "always", dnn = "Metastasis") %>% 
  data.frame() %>% 
  mutate(Percent = colPercent(Freq, pretty.text = TRUE, keep = FALSE, digits = 3))
