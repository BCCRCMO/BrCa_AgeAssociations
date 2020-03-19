### Data for subtype distribution plots/tables


tmadf <- read.table("../BuildData/MB09_wholedata.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "")

VariablesToExport <- c(
    "MB09_id",
    "MetabricID",
    "ERNumber_alpha",
    "sex",
    "age_at_diagnosis",
    "meno_status",
    "er_bgt0_v1",
    "her2_b01v23_v1",
    "pr_bgt0_v1",
    "ck56_b01v2_v1",
    "egfr_b0v12_v1",
    "ki67_bgt14_v1",
    "BrCa4",
    "BrCa5",
    "BrCa8"
)

write.csv(tmadf[, VariablesToExport], file = "MB09_wholedata_Age_BrCaSubtype.CSV")


    
    
