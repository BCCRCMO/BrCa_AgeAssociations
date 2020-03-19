### Data for subtype distribution plots/tables

tmadf <- read.table("../BuildData/02-008_BigSeries_Training.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "")
names(tmadf)[1] <- "bseries_id"

vmadf <- read.table("../BuildData/02-008_BigSeries_Validation.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "")
names(vmadf)[1] <- "bseries_id"

### match3df <- read.table(file = "../../Outcomes/BigSeries_bseries_id_NoOverlap_MB09_METABRICexpression.csv",
###                        stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
###                        na.strings = " ", comment.char = "")

wmadf <- rbind(tmadf, vmadf)

VariablesToExport <- c(
    "bseries_id",
    "sex",
    "age_at_diagnosis",
    "meno_status",
    "er_b0v123_v2",
    "her2_b012v23a_v2",
    "pr_b0v123_v3",
    "ck56_b0v12_v1",
    "egfr_b0v12_v1",
    "ki67_c1v3_v1.6",
    "BrCa4",
    "BrCa5",
    "BrCa8"
)


write.csv(wmadf[, VariablesToExport], file = "02-008_BigSeries_wholedata_Age_BrCaSubtype.CSV")

