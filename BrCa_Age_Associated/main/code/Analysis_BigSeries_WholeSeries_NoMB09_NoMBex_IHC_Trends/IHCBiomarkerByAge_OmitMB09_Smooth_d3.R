### Smooth fits against age.
### Need to omit MB09 cases.  MB09_wholedata_match_to_BigSeries.csv
### c:/kilroy/Projects/TMA/02-008_BigSeries/Analysis/BuildData/MB09_wholedata_match_to_BigSeries.csv
### c:/kilroy/Projects/TMA/02-008_BigSeries/Analysis/BuildData/METABRIC_Expression_wholedata_match_to_BigSeries.csv

require("MOutils")
require("splines")
require("party")

if ( !file.exists("Plots") ) { dir.create("Plots") }

tmadf <- read.table("../BuildData/02-008_BigSeries_Training.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "")
names(tmadf)[1] <- "bseries_id"

vmadf <- read.table("../BuildData/02-008_BigSeries_Validation.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "")
names(vmadf)[1] <- "bseries_id"

match3df <- read.table(file = "../../Outcomes/BigSeries_bseries_id_NoOverlap_MB09_METABRICexpression.csv",
                       stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                       na.strings = " ", comment.char = "")

NULL
### > length(mb09matchdf$bseries_id[which((mb09matchdf$bseries_id %in% c(tmadf$bseries_id, vmadf$bseries_id)))])
### [1] 512    
###
### > length(mb09matchdf$bseries_id[which((mb09matchdf$bseries_id %in% c(tmadf$bseries_id)))])
### [1] 259
### > length(mb09matchdf$bseries_id[which((mb09matchdf$bseries_id %in% c(vmadf$bseries_id)))])
### [1] 253
### >
### > nrow(tmadf)
### [1] 2003
### > 259/2003
### [1] 0.129306
### > nrow(vmadf)
### [1] 1988
### > 253/1988
### [1] 0.1272636
### >   13% of cases in BigSeries appear in MB09.  Remove them and redo smoother runs.


tmadf$BrCa4f <- factor(tmadf$BrCa4,
                       levels = c("Luminalp", "Luminal/HER2+", "HER2+/ER-PR-", "TNP", "Unassigned"))
vmadf$BrCa4f <- factor(vmadf$BrCa4,
                       levels = c("Luminalp", "Luminal/HER2+", "HER2+/ER-PR-", "TNP", "Unassigned"))

##ar
tmadf$ar_c_v1nc  <- as.numeric(tmadf$ar_c_v1)
vmadf$ar_c_v1nc  <- as.numeric(vmadf$ar_c_v1)

tmadf$ar_c1v1_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ar_c1v1_v1n[tmadf$ar_c1v1_v1 == "<75%"] <- 0
tmadf$ar_c1v1_v1n[tmadf$ar_c1v1_v1 == ">=75%"] <- 1

vmadf$ar_c1v1_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$ar_c1v1_v1n[vmadf$ar_c1v1_v1 == "<75%"] <- 0
vmadf$ar_c1v1_v1n[vmadf$ar_c1v1_v1 == ">=75%"] <- 1

tmadf$ar_c_v1nc  <- as.numeric(tmadf$ar_c_v1)
vmadf$ar_c_v1nc  <- as.numeric(vmadf$ar_c_v1)

tmadf$ar_c1v1_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ar_c1v1_v1n[tmadf$ar_c1v1_v1 == "<75%"] <- 0
tmadf$ar_c1v1_v1n[tmadf$ar_c1v1_v1 == ">=75%"] <- 1

vmadf$ar_c1v1_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$ar_c1v1_v1n[vmadf$ar_c1v1_v1 == "<75%"] <- 0
vmadf$ar_c1v1_v1n[vmadf$ar_c1v1_v1 == ">=75%"] <- 1

##bcl2_c1v1_v2
tmadf$bcl2_c_v2.ppnc  <- as.numeric(tmadf$bcl2_c_v2.pp)
vmadf$bcl2_c_v2.ppnc  <- as.numeric(vmadf$bcl2_c_v2.pp)

tmadf$bcl2_bge60_c_v2.ppn <- rep(NA_integer_, nrow(tmadf))
tmadf$bcl2_bge60_c_v2.ppn[tmadf$bcl2_bge60_c_v2.pp == "0% to 50%"] <- 0
tmadf$bcl2_bge60_c_v2.ppn[tmadf$bcl2_bge60_c_v2.pp == ">50% to 100%"] <- 1

vmadf$bcl2_bge60_c_v2.ppn <- rep(NA_integer_, nrow(vmadf))
vmadf$bcl2_bge60_c_v2.ppn[vmadf$bcl2_bge60_c_v2.pp == "0% to 50%"] <- 0
vmadf$bcl2_bge60_c_v2.ppn[vmadf$bcl2_bge60_c_v2.pp == ">50% to 100%"] <- 1

##ca9
tmadf$ca9_b0v123_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ca9_b0v123_v1n[tmadf$ca9_b0v123_v1 == "negative"] <- 0
tmadf$ca9_b0v123_v1n[tmadf$ca9_b0v123_v1 == "any positive"] <- 1

vmadf$ca9_b0v123_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$ca9_b0v123_v1n[vmadf$ca9_b0v123_v1 == "negative"] <- 0
vmadf$ca9_b0v123_v1n[vmadf$ca9_b0v123_v1 == "any positive"] <- 1

## ccnd1_fish_v2  Mostly missing.

## cd163_c_v1
tmadf$cd163_c_v1nc <- as.numeric(tmadf$cd163_c_v1)
vmadf$cd163_c_v1nc <- as.numeric(vmadf$cd163_c_v1)

## ck5/6  krt5
tmadf$ck56_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$ck56_v1nc[tmadf$ck56_v1 == "negative"] <- 0
tmadf$ck56_v1nc[tmadf$ck56_v1 == "weak staining cytoplastmic and/or membranous staining"] <- 1
tmadf$ck56_v1nc[tmadf$ck56_v1 == "strong staining cytoplastmic and/or membranous staining"] <- 2

vmadf$ck56_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$ck56_v1nc[vmadf$ck56_v1 == "negative"] <- 0
vmadf$ck56_v1nc[vmadf$ck56_v1 == "weak staining cytoplastmic and/or membranous staining"] <- 1
vmadf$ck56_v1nc[vmadf$ck56_v1 == "strong staining cytoplastmic and/or membranous staining"] <- 2

tmadf$ck56_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ck56_b0v12_v1n[tmadf$ck56_b0v12_v1 == "negative"] <- 0
tmadf$ck56_b0v12_v1n[tmadf$ck56_b0v12_v1 == "weak or strong cytoplastmic and/or membranous staining"] <- 1

vmadf$ck56_b0v12_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$ck56_b0v12_v1n[vmadf$ck56_b0v12_v1 == "negative"] <- 0
vmadf$ck56_b0v12_v1n[vmadf$ck56_b0v12_v1 == "weak or strong cytoplastmic and/or membranous staining"] <- 1

##cldn3
tmadf$cldn3_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$cldn3_v1nc[tmadf$cldn3_v1 == "negative"] <- 0
tmadf$cldn3_v1nc[tmadf$cldn3_v1 == "any staining"] <- 1
tmadf$cldn3_v1nc[tmadf$cldn3_v1 == "strong stain in >20%"] <- 2

vmadf$cldn3_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$cldn3_v1nc[vmadf$cldn3_v1 == "negative"] <- 0
vmadf$cldn3_v1nc[vmadf$cldn3_v1 == "any staining"] <- 1
vmadf$cldn3_v1nc[vmadf$cldn3_v1 == "strong stain in >20%"] <- 2

tmadf$cldn3_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$cldn3_b0v12_v1n[tmadf$cldn3_v1 == "negative"] <- 0
tmadf$cldn3_b0v12_v1n[tmadf$cldn3_v1 == "any staining"] <- 1
tmadf$cldn3_b0v12_v1n[tmadf$cldn3_v1 == "strong stain in >20%"] <- 1

vmadf$cldn3_b0v12_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$cldn3_b0v12_v1n[vmadf$cldn3_v1 == "negative"] <- 0
vmadf$cldn3_b0v12_v1n[vmadf$cldn3_v1 == "any staining"] <- 1
vmadf$cldn3_b0v12_v1n[vmadf$cldn3_v1 == "strong stain in >20%"] <- 1

tmadf$cldn3_b01v2_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$cldn3_b01v2_v1n[tmadf$cldn3_v1 == "negative"] <- 0
tmadf$cldn3_b01v2_v1n[tmadf$cldn3_v1 == "any staining"] <- 0
tmadf$cldn3_b01v2_v1n[tmadf$cldn3_v1 == "strong stain in >20%"] <- 1

vmadf$cldn3_b01v2_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$cldn3_b01v2_v1n[vmadf$cldn3_v1 == "negative"] <- 0
vmadf$cldn3_b01v2_v1n[vmadf$cldn3_v1 == "any staining"] <- 0
vmadf$cldn3_b01v2_v1n[vmadf$cldn3_v1 == "strong stain in >20%"] <- 1

##cryab alpha-basic crystallin
tmadf$cryab_b0v12_v1.1n <- rep(NA_integer_, nrow(tmadf))
tmadf$cryab_b0v12_v1.1n[tmadf$cryab_b0v12_v1.1 == "negative"] <- 0
tmadf$cryab_b0v12_v1.1n[tmadf$cryab_b0v12_v1.1 == "focal or diffuse positive"] <- 1

vmadf$cryab_b0v12_v1.1n <- rep(NA_integer_, nrow(vmadf))
vmadf$cryab_b0v12_v1.1n[vmadf$cryab_b0v12_v1.1 == "negative"] <- 0
vmadf$cryab_b0v12_v1.1n[vmadf$cryab_b0v12_v1.1 == "focal or diffuse positive"] <- 1

## cyclin E
### table(tmadf$CyclinE_Ariol_Percent_Positive_Nuclein, useNA="always")
tmadf$cycline_v1.1.ppn <- as.numeric(tmadf$CyclinE_Ariol_Percent_Positive_Nuclein)
vmadf$cycline_v1.1.ppn <- as.numeric(vmadf$CyclinE_Ariol_Percent_Positive_Nuclein)


##ecad
tmadf$ecad_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$ecad_v1nc[tmadf$ecad_v1 == "negative"] <- 0
tmadf$ecad_v1nc[tmadf$ecad_v1 == "weak membranous > 50% or strong < 50%"] <- 1
tmadf$ecad_v1nc[tmadf$ecad_v1 == "strong membranous >50%"] <- 2

vmadf$ecad_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$ecad_v1nc[vmadf$ecad_v1 == "negative"] <- 0
vmadf$ecad_v1nc[vmadf$ecad_v1 == "weak membranous > 50% or strong < 50%"] <- 1
vmadf$ecad_v1nc[vmadf$ecad_v1 == "strong membranous >50%"] <- 2

tmadf$ecad_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ecad_b0v12_v1n[tmadf$ecad_b0v12_v1 == "Negative"] <- 0
tmadf$ecad_b0v12_v1n[tmadf$ecad_b0v12_v1 == "Positive"] <- 1

vmadf$ecad_b0v12_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$ecad_b0v12_v1n[vmadf$ecad_b0v12_v1 == "Negative"] <- 0
vmadf$ecad_b0v12_v1n[vmadf$ecad_b0v12_v1 == "Positive"] <- 1

##egfr
tmadf$egfr_v1nc  <- rep(NA_integer_, nrow(tmadf))
tmadf$egfr_v1nc[tmadf$egfr_v1 == "negative"] <- 0
tmadf$egfr_v1nc[tmadf$egfr_v1 == "weak positive (any staining) and cytoplasmic staining"] <- 1
tmadf$egfr_v1nc[tmadf$egfr_v1 == "strong positive (strong staining in >20% of the cells)"] <- 2

vmadf$egfr_v1nc  <- rep(NA_integer_, nrow(vmadf))
vmadf$egfr_v1nc[vmadf$egfr_v1 == "negative"] <- 0
vmadf$egfr_v1nc[vmadf$egfr_v1 == "weak positive (any staining) and cytoplasmic staining"] <- 1
vmadf$egfr_v1nc[vmadf$egfr_v1 == "strong positive (strong staining in >20% of the cells)"] <- 2

tmadf$egfr_b0v12_v1n  <- rep(NA_integer_, nrow(tmadf))
tmadf$egfr_b0v12_v1n[tmadf$egfr_b0v12_v1 == "negative"] <- 0
tmadf$egfr_b0v12_v1n[tmadf$egfr_b0v12_v1 == "weak and strong staining"] <- 1

vmadf$egfr_b0v12_v1n  <- rep(NA_integer_, nrow(vmadf))
vmadf$egfr_b0v12_v1n[vmadf$egfr_b0v12_v1 == "negative"] <- 0
vmadf$egfr_b0v12_v1n[vmadf$egfr_b0v12_v1 == "weak and strong staining"] <- 1


##ephb4
tmadf$ephb4_v1.1nc  <- rep(NA_integer_, nrow(tmadf))
tmadf$ephb4_v1.1nc[tmadf$ephb4_v1.1.int == "Negative"] <- 0
tmadf$ephb4_v1.1nc[tmadf$ephb4_v1.1.int == "Weak stain"] <- 1
tmadf$ephb4_v1.1nc[tmadf$ephb4_v1.1.int == "Intermediate stain"] <- 2
tmadf$ephb4_v1.1nc[tmadf$ephb4_v1.1.int == "Strong stain"] <- 3

vmadf$ephb4_v1.1nc  <- rep(NA_integer_, nrow(vmadf))
vmadf$ephb4_v1.1nc[vmadf$ephb4_v1.1.int == "Negative"] <- 0
vmadf$ephb4_v1.1nc[vmadf$ephb4_v1.1.int == "Weak stain"] <- 1
vmadf$ephb4_v1.1nc[vmadf$ephb4_v1.1.int == "Intermediate stain"] <- 2
vmadf$ephb4_v1.1nc[vmadf$ephb4_v1.1.int == "Strong stain"] <- 3


##er
tmadf$er_v2nc <- rep(NA_integer_, nrow(tmadf))
tmadf$er_v2nc[tmadf$er_v2 == "no nuclei stained or <1%"] <- 0
tmadf$er_v2nc[tmadf$er_v2 == "1-25% of the nuclei stained"] <- 1
tmadf$er_v2nc[tmadf$er_v2 == "25-75% of the nuclei stained"] <- 2
tmadf$er_v2nc[tmadf$er_v2 == ">75% of the nuclei stained"] <- 3

vmadf$er_v2nc <- rep(NA_integer_, nrow(vmadf))
vmadf$er_v2nc[vmadf$er_v2 == "no nuclei stained or <1%"] <- 0
vmadf$er_v2nc[vmadf$er_v2 == "1-25% of the nuclei stained"] <- 1
vmadf$er_v2nc[vmadf$er_v2 == "25-75% of the nuclei stained"] <- 2
vmadf$er_v2nc[vmadf$er_v2 == ">75% of the nuclei stained"] <- 3

tmadf$er_b0v123_v2n <- rep(NA_integer_, nrow(tmadf))
tmadf$er_b0v123_v2n[tmadf$er_b0v123_v2 == "negative"] <- 0
tmadf$er_b0v123_v2n[tmadf$er_b0v123_v2 == "any nuclei staining (score 1-3)"] <- 1

vmadf$er_b0v123_v2n <- rep(NA_integer_, nrow(vmadf))
vmadf$er_b0v123_v2n[vmadf$er_b0v123_v2 == "negative"] <- 0
vmadf$er_b0v123_v2n[vmadf$er_b0v123_v2 == "any nuclei staining (score 1-3)"] <- 1

## er by destran coated charcoal er_dcc

tmadf$er_dcc_v1n <- rep(NA_integer_, nrow(tmadf))
vmadf$er_dcc_v1n <- rep(NA_integer_, nrow(vmadf))

tmadf$er_dcc_v1n <- as.numeric(substring(tmadf$er_result, 2, 4))
vmadf$er_dcc_v1n <- as.numeric(substring(vmadf$er_result, 2, 4))

##ezh2
tmadf$ezh2_v1.ppnc <- as.numeric(tmadf$ezh2_v1.pp)
vmadf$ezh2_v1.ppnc <- as.numeric(vmadf$ezh2_v1.pp)

tmadf$ezh2_blt10vge10_v1.ppn <- rep(NA_integer_, nrow(tmadf))
tmadf$ezh2_blt10vge10_v1.ppn[tmadf$ezh2_blt10vge10_v1.pp == "0 to 5%"] <- 0
tmadf$ezh2_blt10vge10_v1.ppn[tmadf$ezh2_blt10vge10_v1.pp == "10 to 100%"] <- 1
tmadf$ezh2_blt10vge10_v1.ppn[tmadf$ezh2_blt10vge10_v1.pp == "10 to 100 %"] <- 1

vmadf$ezh2_blt10vge10_v1.ppn <- rep(NA_integer_, nrow(vmadf))
vmadf$ezh2_blt10vge10_v1.ppn[vmadf$ezh2_blt10vge10_v1.pp == "0 to 5%"] <- 0
vmadf$ezh2_blt10vge10_v1.ppn[vmadf$ezh2_blt10vge10_v1.pp == "10 to 100%"] <- 1
vmadf$ezh2_blt10vge10_v1.ppn[vmadf$ezh2_blt10vge10_v1.pp == "10 to 100 %"] <- 1

##foxa1
tmadf$foxa1_v1.ppnc <- as.numeric(tmadf$foxa1_v1.percent)
vmadf$foxa1_v1.ppnc <- as.numeric(vmadf$foxa1_v1.percent)

tmadf$foxa1_b0v123_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$foxa1_b0v123_v1n[tmadf$foxa1_v1.intensity == "0"] <- 0
tmadf$foxa1_b0v123_v1n[tmadf$foxa1_v1.intensity %in% c("1", "2", "3")] <- 1

vmadf$foxa1_b0v123_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$foxa1_b0v123_v1n[vmadf$foxa1_v1.intensity == "0"] <- 0
vmadf$foxa1_b0v123_v1n[vmadf$foxa1_v1.intensity %in% c("1", "2", "3")] <- 1

##gata3
tmadf$gata3_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$gata3_v1nc[tmadf$gata3_v1 == "<5% positive nuclei"] <- 0
tmadf$gata3_v1nc[tmadf$gata3_v1 == "5-20% positive nuclei"] <- 1
tmadf$gata3_v1nc[tmadf$gata3_v1 == ">20% positive nuclei"] <- 2

vmadf$gata3_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$gata3_v1nc[vmadf$gata3_v1 == "<5% positive nuclei"] <- 0
vmadf$gata3_v1nc[vmadf$gata3_v1 == "5-20% positive nuclei"] <- 1
vmadf$gata3_v1nc[vmadf$gata3_v1 == ">20% positive nuclei"] <- 2

tmadf$gata3_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$gata3_b0v12_v1n[tmadf$gata3_b0v12_v1 == "<5% positive nuclei"] <- 0
tmadf$gata3_b0v12_v1n[tmadf$gata3_b0v12_v1 == ">=5% positive nuclei"] <- 1

vmadf$gata3_b0v12_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$gata3_b0v12_v1n[vmadf$gata3_b0v12_v1 == "<5% positive nuclei"] <- 0
vmadf$gata3_b0v12_v1n[vmadf$gata3_b0v12_v1 == ">=5% positive nuclei"] <- 1

##h3k27me3
tmadf$h3k27me3_v1.1.ppnc <- as.numeric(tmadf$h3k27me3_v1.1.pp)
vmadf$h3k27me3_v1.1.ppnc <- as.numeric(vmadf$h3k27me3_v1.1.pp)

tmadf$h3k27me3_v1_b012v3n <- rep(NA_integer_, nrow(tmadf))
tmadf$h3k27me3_v1_b012v3n[tmadf$h3k27me3_v1_b012v3 == "Negative to moderate"] <- 0
tmadf$h3k27me3_v1_b012v3n[tmadf$h3k27me3_v1_b012v3 == "Strong"] <- 1

vmadf$h3k27me3_v1_b012v3n <- rep(NA_integer_, nrow(vmadf))
vmadf$h3k27me3_v1_b012v3n[vmadf$h3k27me3_v1_b012v3 == "Negative to moderate"] <- 0
vmadf$h3k27me3_v1_b012v3n[vmadf$h3k27me3_v1_b012v3 == "Strong"] <- 1

tmadf$h3k27me3_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$h3k27me3_v1nc[tmadf$h3k27me3_v1 == "negative (0)"] <- 0
tmadf$h3k27me3_v1nc[tmadf$h3k27me3_v1 == "weak (1+)"] <- 1
tmadf$h3k27me3_v1nc[tmadf$h3k27me3_v1 == "moderate positive (2^)"] <- 2
tmadf$h3k27me3_v1nc[tmadf$h3k27me3_v1 == "moderate positive (2~)"] <- 3
tmadf$h3k27me3_v1nc[tmadf$h3k27me3_v1 == "strong positive (3^)"] <- 4
tmadf$h3k27me3_v1nc[tmadf$h3k27me3_v1 == "strong positive (3~)"] <- 5

vmadf$h3k27me3_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$h3k27me3_v1nc[vmadf$h3k27me3_v1 == "negative (0)"] <- 0
vmadf$h3k27me3_v1nc[vmadf$h3k27me3_v1 == "weak (1+)"] <- 1
vmadf$h3k27me3_v1nc[vmadf$h3k27me3_v1 == "moderate positive (2^)"] <- 2
vmadf$h3k27me3_v1nc[vmadf$h3k27me3_v1 == "moderate positive (2~)"] <- 3
vmadf$h3k27me3_v1nc[vmadf$h3k27me3_v1 == "strong positive (3^)"] <- 4
vmadf$h3k27me3_v1nc[vmadf$h3k27me3_v1 == "strong positive (3~)"] <- 5


##her2
tmadf$her2_v2nc <- rep(NA_integer_, nrow(tmadf))
tmadf$her2_v2nc[tmadf$her2_v2 == "negative"] <- 0
tmadf$her2_v2nc[tmadf$her2_v2 == "weak"] <- 1
tmadf$her2_v2nc[tmadf$her2_v2 == "moderate"] <- 2
tmadf$her2_v2nc[tmadf$her2_v2 == "strong"] <- 3

vmadf$her2_v2nc <- rep(NA_integer_, nrow(vmadf))
vmadf$her2_v2nc[vmadf$her2_v2 == "negative"] <- 0
vmadf$her2_v2nc[vmadf$her2_v2 == "weak"] <- 1
vmadf$her2_v2nc[vmadf$her2_v2 == "moderate"] <- 2
vmadf$her2_v2nc[vmadf$her2_v2 == "strong"] <- 3

tmadf$her2_b012v23a_v2n <- rep(NA_integer_, nrow(tmadf))
tmadf$her2_b012v23a_v2n[tmadf$her2_b012v23a_v2 == "Her2 = {0,1,2 w/ FISH neg}"] <- 0
tmadf$her2_b012v23a_v2n[tmadf$her2_b012v23a_v2 == "Her2 = {2 w/ FISH pos, 3}"] <- 1

vmadf$her2_b012v23a_v2n <- rep(NA_integer_, nrow(vmadf))
vmadf$her2_b012v23a_v2n[vmadf$her2_b012v23a_v2 == "Her2 = {0,1,2 w/ FISH neg}"] <- 0
vmadf$her2_b012v23a_v2n[vmadf$her2_b012v23a_v2 == "Her2 = {2 w/ FISH pos, 3}"] <- 1

##her3
tmadf$her3_v1.1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$her3_v1.1nc[tmadf$her3_v1.1 == "no staining"] <- 0
tmadf$her3_v1.1nc[tmadf$her3_v1.1 == "<20% of tumor cells stained and/or weak staining"] <- 1
tmadf$her3_v1.1nc[tmadf$her3_v1.1 == ">20% of tumor cells stained and/or strong staining"] <- 2

vmadf$her3_v1.1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$her3_v1.1nc[vmadf$her3_v1.1 == "no staining"] <- 0
vmadf$her3_v1.1nc[vmadf$her3_v1.1 == "<20% of tumor cells stained and/or weak staining"] <- 1
vmadf$her3_v1.1nc[vmadf$her3_v1.1 == ">20% of tumor cells stained and/or strong staining"] <- 2


tmadf$her3_b0v12_v1.1n <- rep(NA_integer_, nrow(tmadf))
tmadf$her3_b0v12_v1.1n[tmadf$her3_b0v12_v1.1 == "no staining"] <- 0
tmadf$her3_b0v12_v1.1n[tmadf$her3_b0v12_v1.1 == "any staining"] <- 1

vmadf$her3_b0v12_v1.1n <- rep(NA_integer_, nrow(vmadf))
vmadf$her3_b0v12_v1.1n[vmadf$her3_b0v12_v1.1 == "no staining"] <- 0
vmadf$her3_b0v12_v1.1n[vmadf$her3_b0v12_v1.1 == "any staining"] <- 1

##her4
tmadf$her4_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$her4_v1nc[tmadf$her4_v1 == "no staining"] <- 0
tmadf$her4_v1nc[tmadf$her4_v1 == "weak staining, cytoplasmic and/or membranous"] <- 1
tmadf$her4_v1nc[tmadf$her4_v1 == "strong staining, cytoplasmic and/or membranous"] <- 2

vmadf$her4_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$her4_v1nc[vmadf$her4_v1 == "no staining"] <- 0
vmadf$her4_v1nc[vmadf$her4_v1 == "weak staining, cytoplasmic and/or membranous"] <- 1
vmadf$her4_v1nc[vmadf$her4_v1 == "strong staining, cytoplasmic and/or membranous"] <- 2


tmadf$her4_b01v2_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$her4_b01v2_v1n[tmadf$her4_b01v2_v1 == "no staining or weak staining, cytoplasmic and/or membranous"] <- 0
tmadf$her4_b01v2_v1n[tmadf$her4_b01v2_v1 == "strong staining, cytoplasmic and/or membranous"] <- 1

vmadf$her4_b01v2_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$her4_b01v2_v1n[vmadf$her4_b01v2_v1 == "no staining or weak staining, cytoplasmic and/or membranous"] <- 0
vmadf$her4_b01v2_v1n[vmadf$her4_b01v2_v1 == "strong staining, cytoplasmic and/or membranous"] <- 1

##hsp27
tmadf$hsp27_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$hsp27_v1nc[tmadf$hsp27_v1 == "negative"] <- 0
tmadf$hsp27_v1nc[tmadf$hsp27_v1 == ">5%, weakly positive"] <- 1
tmadf$hsp27_v1nc[tmadf$hsp27_v1 == "moderately positive"] <- 2
tmadf$hsp27_v1nc[tmadf$hsp27_v1 == "strongly positive"] <- 3

vmadf$hsp27_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$hsp27_v1nc[vmadf$hsp27_v1 == "negative"] <- 0
vmadf$hsp27_v1nc[vmadf$hsp27_v1 == ">5%, weakly positive"] <- 1
vmadf$hsp27_v1nc[vmadf$hsp27_v1 == "moderately positive"] <- 2
vmadf$hsp27_v1nc[vmadf$hsp27_v1 == "strongly positive"] <- 3


tmadf$hsp27_b0v123_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$hsp27_b0v123_v1n[tmadf$hsp27_b0v123_v1 == "{0}"] <- 0
tmadf$hsp27_b0v123_v1n[tmadf$hsp27_b0v123_v1 == "{1,2,3}"] <- 1

vmadf$hsp27_b0v123_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$hsp27_b0v123_v1n[vmadf$hsp27_b0v123_v1 == "{0}"] <- 0
vmadf$hsp27_b0v123_v1n[vmadf$hsp27_b0v123_v1 == "{1,2,3}"] <- 1

##igf1r
tmadf$igf1r_v2.0.1.ppn <- rep(NA_integer_, nrow(tmadf))
tmadf$igf1r_v2.0.1.ppn[tmadf$igf1r_v2.0.1.percent == "none"] <- 0
tmadf$igf1r_v2.0.1.ppn[tmadf$igf1r_v2.0.1.percent == "1-10%"] <- 1
tmadf$igf1r_v2.0.1.ppn[tmadf$igf1r_v2.0.1.percent == "11-33%"] <- 2
tmadf$igf1r_v2.0.1.ppn[tmadf$igf1r_v2.0.1.percent == "34-66%"] <- 3
tmadf$igf1r_v2.0.1.ppn[tmadf$igf1r_v2.0.1.percent == "67-100%"] <- 4

vmadf$igf1r_v2.0.1.ppn <- rep(NA_integer_, nrow(vmadf))
vmadf$igf1r_v2.0.1.ppn[vmadf$igf1r_v2.0.1.percent == "none"] <- 0
vmadf$igf1r_v2.0.1.ppn[vmadf$igf1r_v2.0.1.percent == "1-10%"] <- 1
vmadf$igf1r_v2.0.1.ppn[vmadf$igf1r_v2.0.1.percent == "11-33%"] <- 2
vmadf$igf1r_v2.0.1.ppn[vmadf$igf1r_v2.0.1.percent == "34-66%"] <- 3
vmadf$igf1r_v2.0.1.ppn[vmadf$igf1r_v2.0.1.percent == "67-100%"] <- 4


tmadf$igf1r_b0v123_v2.1.1n <- rep(NA_integer_, nrow(tmadf))
tmadf$igf1r_b0v123_v2.1.1n[tmadf$igf1r_b0v123_v2.1.1 == "No staining is observed in invasive tumer cells"] <- 0
tmadf$igf1r_b0v123_v2.1.1n[tmadf$igf1r_b0v123_v2.1.1 == "1,2,3+"] <- 1

vmadf$igf1r_b0v123_v2.1.1n <- rep(NA_integer_, nrow(vmadf))
vmadf$igf1r_b0v123_v2.1.1n[vmadf$igf1r_b0v123_v2.1.1 == "No staining is observed in invasive tumer cells"] <- 0
vmadf$igf1r_b0v123_v2.1.1n[vmadf$igf1r_b0v123_v2.1.1 == "1,2,3+"] <- 1

##igfbp2
tmadf$igfbp2_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$igfbp2_v1nc[tmadf$igfbp2_v1 == "negative"] <- 0
tmadf$igfbp2_v1nc[tmadf$igfbp2_v1 == "weak positive"] <- 1
tmadf$igfbp2_v1nc[tmadf$igfbp2_v1 == "medium positive"] <- 2
tmadf$igfbp2_v1nc[tmadf$igfbp2_v1 == "strong positive"] <- 3

vmadf$igfbp2_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$igfbp2_v1nc[vmadf$igfbp2_v1 == "negative"] <- 0
vmadf$igfbp2_v1nc[vmadf$igfbp2_v1 == "weak positive"] <- 1
vmadf$igfbp2_v1nc[vmadf$igfbp2_v1 == "medium positive"] <- 2
vmadf$igfbp2_v1nc[vmadf$igfbp2_v1 == "strong positive"] <- 3


tmadf$igfbp2_b01v23_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$igfbp2_b01v23_v1n[tmadf$igfbp2_b01v23_v1 == "negative/weak positive"] <- 0
tmadf$igfbp2_b01v23_v1n[tmadf$igfbp2_b01v23_v1 == "medium/strong positive"] <- 1

vmadf$igfbp2_b01v23_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$igfbp2_b01v23_v1n[vmadf$igfbp2_b01v23_v1 == "negative/weak positive"] <- 0
vmadf$igfbp2_b01v23_v1n[vmadf$igfbp2_b01v23_v1 == "medium/strong positive"] <- 1

##inpp4b
tmadf$inpp4b_c_v1.ppn <- as.numeric(tmadf$inpp4b_c_v1.pp)

tmadf$inpp4b_v1.intn <- rep(NA_integer_, nrow(tmadf))
tmadf$inpp4b_v1.intn[tmadf$inpp4b_v1.int == "negative"] <- 0
tmadf$inpp4b_v1.intn[tmadf$inpp4b_v1.int == "weak staining"] <- 1
tmadf$inpp4b_v1.intn[tmadf$inpp4b_v1.int == "moderate staining"] <- 2
tmadf$inpp4b_v1.intn[tmadf$inpp4b_v1.int == "strong staining"] <- 3

tmadf$inpp4b_b0v123_v1 <- tmadf$inpp4b_v1.intn
tmadf$inpp4b_b0v123_v1[tmadf$inpp4b_v1.intn >= 1] <- 1

vmadf$inpp4b_c_v1.ppn <- as.numeric(vmadf$inpp4b_c_v1.pp)

vmadf$inpp4b_v1.intn <- rep(NA_integer_, nrow(vmadf))
vmadf$inpp4b_v1.intn[vmadf$inpp4b_v1.int == "negative"] <- 0
vmadf$inpp4b_v1.intn[vmadf$inpp4b_v1.int == "weak staining"] <- 1
vmadf$inpp4b_v1.intn[vmadf$inpp4b_v1.int == "moderate staining"] <- 2
vmadf$inpp4b_v1.intn[vmadf$inpp4b_v1.int == "strong staining"] <- 3

vmadf$inpp4b_b0v123_v1 <- vmadf$inpp4b_v1.intn
vmadf$inpp4b_b0v123_v1[vmadf$inpp4b_v1.intn >= 1] <- 1

##ki67
tmadf$ki67_c_v1.6nc <- as.numeric(tmadf$ki67_c_v1.6)
vmadf$ki67_c_v1.6nc <- as.numeric(vmadf$ki67_c_v1.6)

tmadf$ki67_c1v3_v1.6n <- rep(NA_integer_, nrow(tmadf))
tmadf$ki67_c1v3_v1.6n[tmadf$ki67_c1v3_v1.6 == "<14%"] <- 0
tmadf$ki67_c1v3_v1.6n[tmadf$ki67_c1v3_v1.6 == ">=14%"] <- 1

vmadf$ki67_c1v3_v1.6n <- rep(NA_integer_, nrow(vmadf))
vmadf$ki67_c1v3_v1.6n[vmadf$ki67_c1v3_v1.6 == "<14%"] <- 0
vmadf$ki67_c1v3_v1.6n[vmadf$ki67_c1v3_v1.6 == ">=14%"] <- 1

## kitm   Ckit mast cells
tmadf$kitm_c_v1n <- as.numeric(tmadf$kitm_c_v1)

vmadf$kitm_c_v1n <- as.numeric(vmadf$kitm_c_v1)

##kitt   Ckit tumour cells
tmadf$kitt_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$kitt_v1nc[tmadf$kitt_v1 == "negative"] <- 0
tmadf$kitt_v1nc[tmadf$kitt_v1 == "weak staining intensity"] <- 1
tmadf$kitt_v1nc[tmadf$kitt_v1 == "strong staining intensity"] <- 2
tmadf$kit_v1nc <- tmadf$kitt_v1nc

vmadf$kitt_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$kitt_v1nc[vmadf$kitt_v1 == "negative"] <- 0
vmadf$kitt_v1nc[vmadf$kitt_v1 == "weak staining intensity"] <- 1
vmadf$kitt_v1nc[vmadf$kitt_v1 == "strong staining intensity"] <- 2
vmadf$kit_v1nc <- vmadf$kitt_v1nc


##krt5
tmadf$krt5_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$krt5_v1nc[tmadf$krt5_v1 == "negative"] <- 0
tmadf$krt5_v1nc[tmadf$krt5_v1 == "any staining"] <- 1
tmadf$krt5_v1nc[tmadf$krt5_v1 == "strong stain in > 20%"] <- 2

vmadf$krt5_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$krt5_v1nc[vmadf$krt5_v1 == "negative"] <- 0
vmadf$krt5_v1nc[vmadf$krt5_v1 == "any staining"] <- 1
vmadf$krt5_v1nc[vmadf$krt5_v1 == "strong stain in > 20%"] <- 2


tmadf$krt5_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$krt5_b0v12_v1n[tmadf$krt5_b0v12_v1 == "{0}"] <- 0
tmadf$krt5_b0v12_v1n[tmadf$krt5_b0v12_v1 == "{1,2}"] <- 1

vmadf$krt5_b0v12_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$krt5_b0v12_v1n[vmadf$krt5_b0v12_v1 == "{0}"] <- 0
vmadf$krt5_b0v12_v1n[vmadf$krt5_b0v12_v1 == "{1,2}"] <- 1

##ku7080
tmadf$ku7080_c_v1nc <- as.numeric(tmadf$ku7080_c_v1)
vmadf$ku7080_c_v1nc <- as.numeric(vmadf$ku7080_c_v1)


tmadf$ku7080_c1v1_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$ku7080_c1v1_v1n[tmadf$ku7080_c1v1_v1 == "<=1%"] <- 0
tmadf$ku7080_c1v1_v1n[tmadf$ku7080_c1v1_v1 == ">1%"] <- 1

vmadf$ku7080_c1v1_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$ku7080_c1v1_v1n[vmadf$ku7080_c1v1_v1 == "<=1%"] <- 0
vmadf$ku7080_c1v1_v1n[vmadf$ku7080_c1v1_v1 == ">1%"] <- 1

##mdm2
tmadf$mdm2_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$mdm2_v1nc[tmadf$mdm2_v1 == "negative"] <- 0
tmadf$mdm2_v1nc[tmadf$mdm2_v1 == "weak positive"] <- 1
tmadf$mdm2_v1nc[tmadf$mdm2_v1 == "strong positive"] <- 2

vmadf$mdm2_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$mdm2_v1nc[vmadf$mdm2_v1 == "negative"] <- 0
vmadf$mdm2_v1nc[vmadf$mdm2_v1 == "weak positive"] <- 1
vmadf$mdm2_v1nc[vmadf$mdm2_v1 == "strong positive"] <- 2


tmadf$mdm2_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$mdm2_b0v12_v1n[tmadf$mdm2_b0v12_v1 == "{0}"] <- 0
tmadf$mdm2_b0v12_v1n[tmadf$mdm2_b0v12_v1 == "{1,2}"] <- 1

vmadf$mdm2_b0v12_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$mdm2_b0v12_v1n[vmadf$mdm2_b0v12_v1 == "{0}"] <- 0
vmadf$mdm2_b0v12_v1n[vmadf$mdm2_b0v12_v1 == "{1,2}"] <- 1

##nestin
tmadf$nestin_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$nestin_v1nc[tmadf$nestin_v1 == "<1% positive tumor cells"] <- 0
tmadf$nestin_v1nc[tmadf$nestin_v1 == "1 to <10% positive tumor cells"] <- 1
tmadf$nestin_v1nc[tmadf$nestin_v1 == ">= 10% positive tumor cells"] <- 2

vmadf$nestin_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$nestin_v1nc[vmadf$nestin_v1 == "<1% positive tumor cells"] <- 0
vmadf$nestin_v1nc[vmadf$nestin_v1 == "1 to <10% positive tumor cells"] <- 1
vmadf$nestin_v1nc[vmadf$nestin_v1 == ">= 10% positive tumor cells"] <- 2
  
tmadf$nestin_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$nestin_v1n[tmadf$nestin_v1 == "<1% positive tumor cells"] <- 0
tmadf$nestin_v1n[tmadf$nestin_v1 == "1 to <10% positive tumor cells"] <- 1
tmadf$nestin_v1n[tmadf$nestin_v1 == ">= 10% positive tumor cells"] <- 1

vmadf$nestin_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$nestin_v1n[vmadf$nestin_v1 == "<1% positive tumor cells"] <- 0
vmadf$nestin_v1n[vmadf$nestin_v1 == "1 to <10% positive tumor cells"] <- 1
vmadf$nestin_v1n[vmadf$nestin_v1 == ">= 10% positive tumor cells"] <- 1
  
##p16
tmadf$p16_c1v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$p16_c1v1nc[tmadf$p16_v1 == "negative"] <- 0
tmadf$p16_c1v1nc[tmadf$p16_v1 == "1-10%"] <- 1
tmadf$p16_c1v1nc[tmadf$p16_v1 == "11-50%"] <- 2
tmadf$p16_c1v1nc[tmadf$p16_v1 == "51-80%"] <- 3
tmadf$p16_c1v1nc[tmadf$p16_v1 == ">80%"] <- 4

vmadf$p16_c1v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$p16_c1v1nc[vmadf$p16_v1 == "negative"] <- 0
vmadf$p16_c1v1nc[vmadf$p16_v1 == "1-10%"] <- 1
vmadf$p16_c1v1nc[vmadf$p16_v1 == "11-50%"] <- 2
vmadf$p16_c1v1nc[vmadf$p16_v1 == "51-80%"] <- 3
vmadf$p16_c1v1nc[vmadf$p16_v1 == ">80%"] <- 4

tmadf$p16_c1v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$p16_c1v1n[tmadf$p16_v1 == "negative"] <- 0
tmadf$p16_c1v1n[tmadf$p16_v1 == "1-10%"] <- 1
tmadf$p16_c1v1n[tmadf$p16_v1 == "11-50%"] <- 1
tmadf$p16_c1v1n[tmadf$p16_v1 == "51-80%"] <- 1
tmadf$p16_c1v1n[tmadf$p16_v1 == ">80%"] <- 1

vmadf$p16_c1v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$p16_c1v1n[vmadf$p16_v1 == "negative"] <- 0
vmadf$p16_c1v1n[vmadf$p16_v1 == "1-10%"] <- 1
vmadf$p16_c1v1n[vmadf$p16_v1 == "11-50%"] <- 1
vmadf$p16_c1v1n[vmadf$p16_v1 == "51-80%"] <- 1
vmadf$p16_c1v1n[vmadf$p16_v1 == ">80%"] <- 1

##p27
tmadf$p27_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$p27_v1nc[tmadf$p27_v1 == "no positive nuclei"] <- 0
tmadf$p27_v1nc[tmadf$p27_v1 == "any positive nuclei, but <25%"] <- 1
tmadf$p27_v1nc[tmadf$p27_v1 == "25-50% nuclei"] <- 2
tmadf$p27_v1nc[tmadf$p27_v1 == ">50% nuclei"] <- 3

vmadf$p27_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$p27_v1nc[vmadf$p27_v1 == "no positive nuclei"] <- 0
vmadf$p27_v1nc[vmadf$p27_v1 == "any positive nuclei, but <25%"] <- 1
vmadf$p27_v1nc[vmadf$p27_v1 == "25-50% nuclei"] <- 2
vmadf$p27_v1nc[vmadf$p27_v1 == ">50% nuclei"] <- 3


tmadf$p27_b012v3_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$p27_b012v3_v1n[tmadf$p27_b012v3_v1 == "<=50% nuclei"] <- 0
tmadf$p27_b012v3_v1n[tmadf$p27_b012v3_v1 == ">50% nuclei"] <- 1

vmadf$p27_b012v3_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$p27_b012v3_v1n[vmadf$p27_b012v3_v1 == "<=50% nuclei"] <- 0
vmadf$p27_b012v3_v1n[vmadf$p27_b012v3_v1 == ">50% nuclei"] <- 1
  
##p53
tmadf$p53_c_v1.1.1.pnnc <- as.numeric(tmadf$p53_c_v1.1.1.pn)
vmadf$p53_c_v1.1.1.pnnc <- as.numeric(vmadf$p53_c_v1.1.1.pn)


tmadf$p53_c1v1_v1.1.1n <- rep(NA_integer_, nrow(tmadf))
tmadf$p53_c1v1_v1.1.1n[tmadf$p53_c1v1_v1.1.1 == "<=10% positive nuclei"] <- 0
tmadf$p53_c1v1_v1.1.1n[tmadf$p53_c1v1_v1.1.1 == ">10% positive nuclei"] <- 1

vmadf$p53_c1v1_v1.1.1n <- rep(NA_integer_, nrow(vmadf))
vmadf$p53_c1v1_v1.1.1n[vmadf$p53_c1v1_v1.1.1 == "<=10% positive nuclei"] <- 0
vmadf$p53_c1v1_v1.1.1n[vmadf$p53_c1v1_v1.1.1 == ">10% positive nuclei"] <- 1
  
##pcad
tmadf$pcad_v1.1.ppnc <- as.numeric(tmadf$pcad_v1.1.pp)
vmadf$pcad_v1.1.ppnc <- as.numeric(vmadf$pcad_v1.1.pp)

tmadf$pcad_ble50vgt50_v1.1.ppn <- rep(NA_integer_, nrow(tmadf))
tmadf$pcad_ble50vgt50_v1.1.ppn[tmadf$pcad_ble50vgt50_v1.1.pp == "0% to 50%"] <- 0
tmadf$pcad_ble50vgt50_v1.1.ppn[tmadf$pcad_ble50vgt50_v1.1.pp == "60% or higher"] <- 1

vmadf$pcad_ble50vgt50_v1.1.ppn <- rep(NA_integer_, nrow(vmadf))
vmadf$pcad_ble50vgt50_v1.1.ppn[vmadf$pcad_ble50vgt50_v1.1.pp == "0% to 50%"] <- 0
vmadf$pcad_ble50vgt50_v1.1.ppn[vmadf$pcad_ble50vgt50_v1.1.pp == "60% or higher"] <- 1

##pgp
tmadf$pgp_v1nc <- as.numeric(tmadf$PGP)
vmadf$pgp_v1nc <- as.numeric(vmadf$PGP)


tmadf$pgp_c1v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$pgp_c1v1n[tmadf$pgp_v1 == "negative"] <- 0
tmadf$pgp_c1v1n[tmadf$pgp_v1 == "any staining"] <- 1
tmadf$pgp_c1v1n[tmadf$pgp_v1 == "strong stain in >20%"] <- 1

vmadf$pgp_c1v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$pgp_c1v1n[vmadf$pgp_v1 == "negative"] <- 0
vmadf$pgp_c1v1n[vmadf$pgp_v1 == "any staining"] <- 1
vmadf$pgp_c1v1n[vmadf$pgp_v1 == "strong stain in >20%"] <- 1

##pipki
tmadf$pipki_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$pipki_v1nc[tmadf$pipki_v1 == "negative"] <- 0
tmadf$pipki_v1nc[tmadf$pipki_v1 == "incomplete membranous / weak membranous in <10% cells"] <- 1
tmadf$pipki_v1nc[tmadf$pipki_v1 == "strong / complete membranous in >10% cells"] <- 2

vmadf$pipki_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$pipki_v1nc[vmadf$pipki_v1 == "negative"] <- 0
vmadf$pipki_v1nc[vmadf$pipki_v1 == "incomplete membranous / weak membranous in <10% cells"] <- 1
vmadf$pipki_v1nc[vmadf$pipki_v1 == "strong / complete membranous in >10% cells"] <- 2


tmadf$pipki_b0v12_v1n <- rep(NA_integer_, nrow(tmadf))
tmadf$pipki_b0v12_v1n[tmadf$pipki_b0v12_v1 == "negative"] <- 0
tmadf$pipki_b0v12_v1n[tmadf$pipki_b0v12_v1 == "any positive"] <- 1

vmadf$pipki_b0v12_v1n <- rep(NA_integer_, nrow(vmadf))
vmadf$pipki_b0v12_v1n[vmadf$pipki_b0v12_v1 == "negative"] <- 0
vmadf$pipki_b0v12_v1n[vmadf$pipki_b0v12_v1 == "any positive"] <- 1

##plau
tmadf$plau_c1v1 <- rep(NA_integer_, nrow(tmadf))
tmadf$plau_c1v1[tmadf$plau_v1 == "negative"] <- 0 
tmadf$plau_c1v1[tmadf$plau_v1 == "positive"] <- 1

vmadf$plau_c1v1 <- rep(NA_integer_, nrow(vmadf))
vmadf$plau_c1v1[vmadf$plau_v1 == "negative"] <- 0 
vmadf$plau_c1v1[vmadf$plau_v1 == "positive"] <- 1

##podxl
tmadf$podxl_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$podxl_v1nc[tmadf$podxl_v1 == "negative"] <- 0
tmadf$podxl_v1nc[tmadf$podxl_v1 == "any weak membranous staining"] <- 1
tmadf$podxl_v1nc[tmadf$podxl_v1 == ">10 and <50% strong membranous staining"] <- 2
tmadf$podxl_v1nc[tmadf$podxl_v1 == ">50% strong membranous staining"] <- 3

vmadf$podxl_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$podxl_v1nc[vmadf$podxl_v1 == "negative"] <- 0
vmadf$podxl_v1nc[vmadf$podxl_v1 == "any weak membranous staining"] <- 1
vmadf$podxl_v1nc[vmadf$podxl_v1 == ">10 and <50% strong membranous staining"] <- 2
vmadf$podxl_v1nc[vmadf$podxl_v1 == ">50% strong membranous staining"] <- 3


tmadf$podxl_c1v1 <- rep(NA_integer_, nrow(tmadf))
tmadf$podxl_c1v1[tmadf$podxl_b01v23_v1 == "negative or any weak membranous staining"] <- 0
tmadf$podxl_c1v1[tmadf$podxl_b01v23_v1 == ">10% strong membranous staining"] <- 1

vmadf$podxl_c1v1 <- rep(NA_integer_, nrow(vmadf))
vmadf$podxl_c1v1[vmadf$podxl_b01v23_v1 == "negative or any weak membranous staining"] <- 0
vmadf$podxl_c1v1[vmadf$podxl_b01v23_v1 == ">10% strong membranous staining"] <- 1

##pr
tmadf$pr_v3nc <- rep(NA_integer_, nrow(tmadf))
tmadf$pr_v3nc[tmadf$pr_v3 == "<1% tumor nuclei stained"] <- 0 
tmadf$pr_v3nc[tmadf$pr_v3 == "1-25% tumor nuclei stained"] <- 1
tmadf$pr_v3nc[tmadf$pr_v3 == "25-75% tumor nuclei stained"] <- 2
tmadf$pr_v3nc[tmadf$pr_v3 == ">75% tumor nuclei stained"] <- 3

vmadf$pr_v3nc <- rep(NA_integer_, nrow(vmadf))
vmadf$pr_v3nc[vmadf$pr_v3 == "<1% tumor nuclei stained"] <- 0 
vmadf$pr_v3nc[vmadf$pr_v3 == "1-25% tumor nuclei stained"] <- 1
vmadf$pr_v3nc[vmadf$pr_v3 == "25-75% tumor nuclei stained"] <- 2
vmadf$pr_v3nc[vmadf$pr_v3 == ">75% tumor nuclei stained"] <- 3


tmadf$pr_b0v123_v3n <- rep(NA_integer_, nrow(tmadf))
tmadf$pr_b0v123_v3n[tmadf$pr_b0v123_v3 == "<1% tumor nuclei stained"] <- 0 
tmadf$pr_b0v123_v3n[tmadf$pr_b0v123_v3 == ">=1% tumor nuclei stained"] <- 1

vmadf$pr_b0v123_v3n <- rep(NA_integer_, nrow(vmadf))
vmadf$pr_b0v123_v3n[vmadf$pr_b0v123_v3 == "<1% tumor nuclei stained"] <- 0 
vmadf$pr_b0v123_v3n[vmadf$pr_b0v123_v3 == ">=1% tumor nuclei stained"] <- 1

##psf
tmadf$psf_c_v1nc <- as.numeric(tmadf$psf_c_v1)
vmadf$psf_c_v1nc <- as.numeric(vmadf$psf_c_v1)


tmadf$psf_c1v1_v2 <- rep(NA_integer_, nrow(tmadf))
tmadf$psf_c1v1_v2[tmadf$psf_c1v1_v1 == "<=35%"] <- 0
tmadf$psf_c1v1_v2[tmadf$psf_c1v1_v1 == ">35%"] <- 1

vmadf$psf_c1v1_v2 <- rep(NA_integer_, nrow(vmadf))
vmadf$psf_c1v1_v2[vmadf$psf_c1v1_v1 == "<=35%"] <- 0
vmadf$psf_c1v1_v2[vmadf$psf_c1v1_v1 == ">35%"] <- 1

##pten
tmadf$pten_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$pten_v1nc[tmadf$pten_v1 == "negative"] <- 0
tmadf$pten_v1nc[tmadf$pten_v1 == "weak positive"] <- 1
tmadf$pten_v1nc[tmadf$pten_v1 == "strong positive"] <- 2

vmadf$pten_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$pten_v1nc[vmadf$pten_v1 == "negative"] <- 0
vmadf$pten_v1nc[vmadf$pten_v1 == "weak positive"] <- 1
vmadf$pten_v1nc[vmadf$pten_v1 == "strong positive"] <- 2

tmadf$pten_c1v1 <- rep(NA_integer_, nrow(tmadf))
tmadf$pten_c1v1[tmadf$pten_v1 == "negative"] <- 0
tmadf$pten_c1v1[tmadf$pten_v1 == "weak positive"] <- 1
tmadf$pten_c1v1[tmadf$pten_v1 == "strong positive"] <- 1

vmadf$pten_c1v1 <- rep(NA_integer_, nrow(vmadf))
vmadf$pten_c1v1[vmadf$pten_v1 == "negative"] <- 0
vmadf$pten_c1v1[vmadf$pten_v1 == "weak positive"] <- 1
vmadf$pten_c1v1[vmadf$pten_v1 == "strong positive"] <- 1

##ret
tmadf$ret_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$ret_v1nc[tmadf$ret_v1 == "negative"] <- 0
tmadf$ret_v1nc[tmadf$ret_v1 == "weak staining"] <- 1
tmadf$ret_v1nc[tmadf$ret_v1 == "strong staining"] <- 2

vmadf$ret_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$ret_v1nc[vmadf$ret_v1 == "negative"] <- 0
vmadf$ret_v1nc[vmadf$ret_v1 == "weak staining"] <- 1
vmadf$ret_v1nc[vmadf$ret_v1 == "strong staining"] <- 2


tmadf$ret_c1v1 <- rep(NA_integer_, nrow(tmadf))
tmadf$ret_c1v1[tmadf$ret_b0v12_v1 == "{0}"] <- 0
tmadf$ret_c1v1[tmadf$ret_b0v12_v1 == "{1,2}"] <- 1

vmadf$ret_c1v1 <- rep(NA_integer_, nrow(vmadf))
vmadf$ret_c1v1[vmadf$ret_b0v12_v1 == "{0}"] <- 0
vmadf$ret_c1v1[vmadf$ret_b0v12_v1 == "{1,2}"] <- 1

##skp2
tmadf$skp2_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$skp2_v1nc[tmadf$skp2_v1 == "<10%"] <- 0
tmadf$skp2_v1nc[tmadf$skp2_v1 == "10-50%"] <- 1
tmadf$skp2_v1nc[tmadf$skp2_v1 == ">50%"] <- 2

vmadf$skp2_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$skp2_v1nc[vmadf$skp2_v1 == "<10%"] <- 0
vmadf$skp2_v1nc[vmadf$skp2_v1 == "10-50%"] <- 1
vmadf$skp2_v1nc[vmadf$skp2_v1 == ">50%"] <- 2


tmadf$skp2_c1v1 <- rep(NA_integer_, nrow(tmadf))
tmadf$skp2_c1v1[tmadf$skp2_b0v12_v1 == "<10%"] <- 0
tmadf$skp2_c1v1[tmadf$skp2_b0v12_v1 == ">=10%"] <- 1

vmadf$skp2_c1v1 <- rep(NA_integer_, nrow(vmadf))
vmadf$skp2_c1v1[vmadf$skp2_b0v12_v1 == "<10%"] <- 0
vmadf$skp2_c1v1[vmadf$skp2_b0v12_v1 == ">=10%"] <- 1

##trim29
tmadf$trim29_v1nc <- rep(NA_integer_, nrow(tmadf))
tmadf$trim29_v1nc[tmadf$trim29_v1 == "background or lower"] <- 0
tmadf$trim29_v1nc[tmadf$trim29_v1 == "above background in <100% of cancer cells"] <- 1
tmadf$trim29_v1nc[tmadf$trim29_v1 == "above background in 100% of cancer cells or strong stain in >20%"] <- 2

vmadf$trim29_v1nc <- rep(NA_integer_, nrow(vmadf))
vmadf$trim29_v1nc[vmadf$trim29_v1 == "background or lower"] <- 0
vmadf$trim29_v1nc[vmadf$trim29_v1 == "above background in <100% of cancer cells"] <- 1
vmadf$trim29_v1nc[vmadf$trim29_v1 == "above background in 100% of cancer cells or strong stain in >20%"] <- 2


tmadf$trim29_c1v1 <- rep(NA_integer_, nrow(tmadf))
tmadf$trim29_c1v1[tmadf$trim29_b0v12_v1 == "{0}"] <- 0
tmadf$trim29_c1v1[tmadf$trim29_b0v12_v1 == "{1,2}"] <- 1

vmadf$trim29_c1v1 <- rep(NA_integer_, nrow(vmadf))
vmadf$trim29_c1v1[vmadf$trim29_b0v12_v1 == "{0}"] <- 0
vmadf$trim29_c1v1[vmadf$trim29_b0v12_v1 == "{1,2}"] <- 1

##yb1
tmadf$yb1_v1.2nc <- rep(NA_integer_, nrow(tmadf))
tmadf$yb1_v1.2nc[tmadf$yb1_v1.2 == "negative"] <- 0
tmadf$yb1_v1.2nc[tmadf$yb1_v1.2 == "weakly positive in >=50% cells"] <- 1
tmadf$yb1_v1.2nc[tmadf$yb1_v1.2 == "moderately positive in >=50% cells"] <- 2
tmadf$yb1_v1.2nc[tmadf$yb1_v1.2 == "strongly positive"] <- 3

vmadf$yb1_v1.2nc <- rep(NA_integer_, nrow(vmadf))
vmadf$yb1_v1.2nc[vmadf$yb1_v1.2 == "negative"] <- 0
vmadf$yb1_v1.2nc[vmadf$yb1_v1.2 == "weakly positive in >=50% cells"] <- 1
vmadf$yb1_v1.2nc[vmadf$yb1_v1.2 == "moderately positive in >=50% cells"] <- 2
vmadf$yb1_v1.2nc[vmadf$yb1_v1.2 == "strongly positive"] <- 3


tmadf$yb1_c1v1.2 <- rep(NA_integer_, nrow(tmadf))
tmadf$yb1_c1v1.2[tmadf$yb1_b0v123_v1.2 == "YB-1 = {0}"] <- 0
tmadf$yb1_c1v1.2[tmadf$yb1_b0v123_v1.2 == "YB-1 = {1,2,3}"] <- 1

vmadf$yb1_c1v1.2 <- rep(NA_integer_, nrow(vmadf))
vmadf$yb1_c1v1.2[vmadf$yb1_b0v123_v1.2 == "YB-1 = {0}"] <- 0
vmadf$yb1_c1v1.2[vmadf$yb1_b0v123_v1.2 == "YB-1 = {1,2,3}"] <- 1

all.equal(names(tmadf), names(vmadf))

wmadf <- rbind(tmadf, vmadf)
wmadf$kit_v1nc <- wmadf$kitt_v1nc





inpp4b_ppint.ctreeout <-
  ctree(Surv(survyrs, brdeath_os) ~ inpp4b_c_v1.ppn + inpp4b_v1.intn, data = tmadf,
        controls = ctree_control(mincriterion = 0.01, maxdepth = 4),
        subset = !is.na(tmadf$inpp4b_c_v1.ppn))



### h3k27me3_v1_b012v3
pdf("./Plots/H3K27me3RatesByAgev01.pdf", width = 6, height = 4)
win.metafile("./Plots/H3K27me3RatesByAgev06.wmf", width = 6, height = 4)
par(mfrow = c(1, 2), oma = c(2, 2, 0, 0), mar = c(4, 2, 4, 1))
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

idxp <- !is.na(tmadf$h3k27me3_v1_b012v3)
set.seed(3123)
with(tmadf[idxp, ], plot(age_at_diagnosis, jitter(1.0 * (h3k27me3_v1_b012v3 == "Strong"), factor = 0.1),
                         ylab = "", xlab = "", xlim = xlims, ylim = ylims))
with(tmadf[idxp, ], lines(supsmu(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"), bass = 3), lwd = 3))
abline(h = 0.5, v = 50, lty = 2)
title(main = "Training set")

vidxp <- !is.na(vmadf$h3k27me3_v1_b012v3)
set.seed(3123)
with(vmadf[vidxp, ], plot(age_at_diagnosis, jitter(1.0 * (h3k27me3_v1_b012v3 == "Strong"), factor = 0.1),
                         ylab = "", xlab = "", xlim = xlims, ylim = ylims))
with(vmadf[vidxp, ], lines(supsmu(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"), bass = 3), lwd = 3))
abline(h = 0.5, v = 50, lty = 2)
title(main = "Validation set")
mtext(text = "Patient age", side = 1, outer = TRUE)
mtext(text = "Proportion H3K27me3+", side = 2, outer = TRUE)
dev.off()

### loess smooth

pdf("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesLoessv01.pdf", width = 8, height = 10)
### win.graph()
par(mfrow = c(1, 2), oma = c(2, 2, 0, 0), mar = c(4, 2, 4, 1))
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

idxp <- !is.na(tmadf$h3k27me3_v1_b012v3)
set.seed(3123)
with(tmadf[idxp, ], scatter.smooth(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"), lwd = 3))
abline(h = 0.5, v = 50, lty = 2)
title(main = "Training set")

vidxp <- !is.na(vmadf$h3k27me3_v1_b012v3)
set.seed(3123)
with(vmadf[vidxp, ], scatter.smooth(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"), lwd = 3))
abline(h = 0.5, v = 50, lty = 2)
title(main = "Validation set")
mtext(text = "Patient age", side = 1, outer = TRUE)
mtext(text = "Proportion H3K27me3+", side = 2, outer = TRUE)
dev.off()

### Do plots by BrCa subtype



### h3k27me3_v1_b012v3
pdf("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesSpan4Bass10v03.pdf", width = 6, height = 4)
##win.metafile("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesv01.wmf", width = 6, height = 4)
##par(mfrow = c(1, 2), oma = c(2, 2, 0, 0), mar = c(4, 2, 4, 1))
pdf("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesSpan4Bass10v04.pdf", width = 8, height = 10)
win.metafile("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesSpan4Bass10v01.wmf", width = 8, height = 10)
par(mfrow = c(5, 2), oma = c(4, 8, 2, 8), mar = c(2, 5, 3, 5))
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

for ( bi in seq(along = unique(tmadf$BrCa4)) ) {
  brcai <- unique(tmadf$BrCa4)[bi]
  idxp <- ( ( !is.na(tmadf$h3k27me3_v1_b012v3) ) & ( tmadf$BrCa4 == brcai ) )
  set.seed(3123)
  with(tmadf[idxp, ], plot(age_at_diagnosis, jitter(1.0 * (h3k27me3_v1_b012v3 == "Strong"), factor = 0.1),
                           ylab = "", xlab = "", xlim = xlims, ylim = ylims))
  with(tmadf[idxp, ], lines(supsmu(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"),
                                   span = 0.4, bass = 10), lwd = 3))
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Training:",  brcai) )
  
  vidxp <- ( ( !is.na(vmadf$h3k27me3_v1_b012v3) ) & ( vmadf$BrCa4 == brcai ) )
  set.seed(3123)
  with(vmadf[vidxp, ], plot(age_at_diagnosis, jitter(1.0 * (h3k27me3_v1_b012v3 == "Strong"), factor = 0.1),
                            ylab = "", xlab = "", xlim = xlims, ylim = ylims))
  with(vmadf[vidxp, ], lines(supsmu(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"),
                                    span = 0.4, bass = 10), lwd = 3))
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Validation:", brcai) )
  mtext(text = "Patient age", side = 1, outer = TRUE)
  mtext(text = "Proportion H3K27me3+", side = 2, outer = TRUE)

}

dev.off()

### Loess

pdf("./Plots/H3K27me3RatesWithinBrCaSubtypesLoessv01.pdf", width = 8, height = 10)
##win.metafile("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesv01.wmf", width = 6, height = 4)
##par(mfrow = c(1, 2), oma = c(2, 2, 0, 0), mar = c(4, 2, 4, 1))

par(mfrow = c(5, 2), oma = c(4, 8, 2, 8), mar = c(2, 5, 3, 5))
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

for ( bi in seq(along = unique(tmadf$BrCa4)) ) {
  brcai <- unique(tmadf$BrCa4)[bi]
  idxp <- ( ( !is.na(tmadf$h3k27me3_v1_b012v3) ) & ( tmadf$BrCa4 == brcai ) )
  set.seed(3123)
  with(tmadf[idxp, ], scatter.smooth( age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3 == "Strong"), span = 7/8,
                                     lwd = 3 ) )
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Training:",  brcai) )
  
  vidxp <- ( ( !is.na(vmadf$h3k27me3_v1_b012v3) ) & ( vmadf$BrCa4 == brcai ) )
  set.seed(3123)
  with(vmadf[vidxp, ], scatter.smooth(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3 == "Strong"), span = 7/8,
                                      lwd = 3) )
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Validation:", brcai) )
  mtext(text = "Patient age", side = 1, outer = TRUE)
  mtext(text = "Proportion H3K27me3+", side = 2, outer = TRUE)

}

dev.off()


### h3k27me3_v1_b012v3
pdf("./Plots/H3K27me3RatesByAgeLT80WithinBrCaSubtypesSpan4Bass10v03.pdf", width = 6, height = 4)
##win.metafile("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesv01.wmf", width = 6, height = 4)
##par(mfrow = c(1, 2), oma = c(2, 2, 0, 0), mar = c(4, 2, 4, 1))
pdf("./Plots/H3K27me3RatesByAgeLT75WithinBrCaSubtypesSpan4Bass10v01.pdf", width = 8, height = 10)
par(mfrow = c(5, 2), oma = c(4, 8, 2, 8), mar = c(2, 5, 3, 5))
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

for ( bi in seq(along = unique(tmadf$BrCa4)) ) {
  brcai <- unique(tmadf$BrCa4)[bi]
  idxp <- ( ( !is.na(tmadf$h3k27me3_v1_b012v3) ) & ( tmadf$BrCa4 == brcai ) & ( tmadf$age_at_diagnosis <= 75) )
  set.seed(3123)
  with(tmadf[idxp, ], plot(age_at_diagnosis, jitter(1.0 * (h3k27me3_v1_b012v3 == "Strong"), factor = 0.1),
                           ylab = "", xlab = "", xlim = xlims, ylim = ylims))
  with(tmadf[idxp, ], lines(supsmu(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"),
                                   span = 0.4, bass = 10), lwd = 3))
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Training:",  brcai) )
  
  vidxp <- ( ( !is.na(vmadf$h3k27me3_v1_b012v3) ) & ( vmadf$BrCa4 == brcai ) & (vmadf$age_at_diagnosis <= 75) )
  set.seed(3123)
  with(vmadf[vidxp, ], plot(age_at_diagnosis, jitter(1.0 * (h3k27me3_v1_b012v3 == "Strong"), factor = 0.1),
                            ylab = "", xlab = "", xlim = xlims, ylim = ylims))
  with(vmadf[vidxp, ], lines(supsmu(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"),
                                    span = 0.4, bass = 10), lwd = 3))
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Validation:", brcai) )
  mtext(text = "Patient age", side = 1, outer = TRUE)
  mtext(text = "Proportion H3K27me3+", side = 2, outer = TRUE)

}

dev.off()


pdf("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesWholeSeriesSpan3Bass10v01.pdf", width = 8, height = 10)

win.metafile("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesWholeSeriesSpan4Bass10v01.wmf", width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(4, 4, 2, 4), mar = c(2, 5, 3, 2))
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

for ( bi in seq(along = unique(wmadf$BrCa4)) ) {
  brcai <- unique(wmadf$BrCa4)[bi]
  idxp <- ( ( !is.na(wmadf$h3k27me3_v1_b012v3) ) & ( wmadf$BrCa4 == brcai ) )
  set.seed(3123)
  with(wmadf[idxp, ], plot(age_at_diagnosis, jitter(1.0 * (h3k27me3_v1_b012v3 == "Strong"), factor = 0.1),
                           ylab = "", xlab = "", xlim = xlims, ylim = ylims))
  with(wmadf[idxp, ], lines(supsmu(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"),
                                   span = 0.4, bass = 10), lwd = 3))
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Whole Set:",  brcai) )
  
 

}

dev.off()



pdf("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesWholeSeriesLoessv01.pdf", width = 8, height = 10)

###win.metafile("./Plots/H3K27me3RatesByAgeWithinBrCaSubtypesWholeSeriesLoessv01.wmf", width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(4, 4, 2, 4), mar = c(2, 5, 3, 2))
xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

for ( bi in seq(along = unique(wmadf$BrCa4)) ) {
  brcai <- unique(wmadf$BrCa4)[bi]
  idxp <- ( ( !is.na(wmadf$h3k27me3_v1_b012v3) ) & ( wmadf$BrCa4 == brcai ) )
  set.seed(3123)
  with(wmadf[idxp, ], scatter.smooth(age_at_diagnosis, 1.0 * (h3k27me3_v1_b012v3== "Strong"),
                                     family = "gaussian", lwd = 3))
  abline(h = 0.5, v = 50, lty = 2)
  title(main = paste("Whole Set:",  brcai) )
}

dev.off()


### Look at age relationship with other biomarkers.

### Age by EZH2


if(FALSE){

ar_c1v1_v1n                    ##ar_c1v1_v1                     ##ar
bcl2_bge60_c_v2.ppn            ## bcl2_bge60_c_v2.pp             ##bcl2_c1v1_v2
ca9_b0v123_v1n                 ## ca9_b0v123_v1                  ##ca9
##cd163  TAM count, not binarized.  Do scatterplot.
##cd8  CD8 positive lymphocyte count, not binarized.  Do scatterplot.
ck56_b0v12_v1n                 ## ck56_b0v12_v1   ## ck5/6  krt5
##cldn3  Make 2 binary var:  binarize at: neg vs pos;  strong stain in >20% vs else
cldn3_b0v12_v1n                ##cldn3
cldn3_b01v2_v1n                ##cldn3
##cox2  Skip - blocks A-J only
### cox2_b0v123_v1 <- NULL
cryab_b0v12_v1.1n              ## cryab_b0v12_v1.1 ##cryab alpha-basic crystallin (as opposed to alpha-acidic)
ecad_b0v12_v1n                 ##ecad_b0v12_v1  ##ecad
egfr_b0v12_v1n                 ##egfr_b0v12_v1 ##egfr
er_b0v123_v2n                  ##er_b0v123_v2  ##er
ezh2_blt10vge10_v1.ppn         ##ezh2_blt10vge10_v1.pp   ##ezh2
##foxa1:  Review intensity and percent - binarize at 0.
### tmadf$foxa1_v1.intf <- factor(tmadf$foxa1_v1.intensity, levels = as.character(0:3))
### tmadf$foxa1_v1.ppn <- as.numeric(tmadf$foxa1_v1.percent)
### foxa1_v1.intensity
### foxa1_v1.percent
foxa1_b0v123_v1n               ##foxa1_v1.intensity  ##foxa1
gata3_b0v12_v1n                ##gata3_b0v12_v1       ##gata3
h3k27me3_v1_b012v3n            ##h3k27me3_v1_b012v3   ##h3k27me3
her2_b012v23a_v2n              ##her2_b012v23a_v2     ##her2
her3_b0v12_v1.1n               ##her3_b0v12_v1.1      ##her3
her4_b01v2_v1n                 ##her4_b01v2_v1        ##her4
hsp27_b0v123_v1n               ##hsp27_b0v123_v1      ##hsp27
igf1r_b0v123_v2.1.1n           ## igf1r_v2.0.1.intensity, igf1r_v2.intensity ##igf1r
igfbp2_b01v23_v1n              ##igfbp2_b01v23_v1               ##igfbp2
inpp4b_b0v123_v1               ## inpp4b_c_v1.pp  inpp4b_v1.int  inpp4b
ki67_c1v3_v1.6n                ##ki67_c1v3_v1.6  ##ki67
kitm_c_v1n                     ## kitm_c_v1 kitm   Ckit mast cells count.  Do scatterplot.
kitt_c_v1n                     ## kitt   Ckit tumour cells count.  Do scatterplot.
##kitt_b0v12_v1                ##kitt   Ckit tumour cells count.  Few positives.  Skip.
krt5_b0v12_v1n                 ##krt5_b0v12_v1  ##krt5
ku7080_c1v1_v1n                ##ku7080_c1v1_v1                 ##ku7080
mdm2_b0v12_v1n                 ##mdm2_b0v12_v1                  ##mdm2  Blocks A to J only
nestin_v1n                     ##nestin_v1   ##nestin 
p16_c1v1n                      ##p16
p27_b012v3_v1n                 ##p27_b012v3_v1                  ##p27
p53_c1v1_v1.1.1n               ## p53_c1v1_v1.1.1 ##p53
pcad_ble50vgt50_v1.1.ppn       ##pcad_ble50vgt50_v1.1.pp        ##pcad_c1v1_v1.1
pgp_c1v1n                      ## pgp_v1 pgp
pipki_b0v12_v1                 ## pipki_b01v2_v1(too few pos cases)   pipki
plau_c1v1                      ##plau_v1 ##plau
podxl_c1v1                     ##podxl_b01v23_v1 ##podxl
pr_b0v123_v3n                  ##pr_b0v123_v3  ##pr
psf_c1v1_v2                    ##psf_c1v1_v1 ##psf
pten_c1v1                      ##pten_v1 ##pten
ret_c1v1                       ##ret_b0v12_v1 ##ret
skp2_c1v1                      ##skp2_b0v12_v1 ##skp2
trim29_c1v1                    ##trim29_b0v12_v1  ##trim29_b01v2_v1  ##trim29
yb1_c1v1.2                     ##yb1_b0v123_v1.2

} ## end if(FALSE)

## Binary scatterplot smooths

bscvars <- c(
  "ar_c1v1_v1n",       ##ar_c1v1_v1                     ##ar
  "bcl2_bge60_c_v2.ppn", ## bcl2_bge60_c_v2.pp             ##bcl2_c1v1_v2
  "ca9_b0v123_v1n",      ## ca9_b0v123_v1                  ##ca9
  ##cd163,  TAM count, not binarized.  Do scatterplot.
  ##cd8,  CD8 positive lymphocyte count, not binarized.  Do scatterplot.
  "ck56_b0v12_v1n", ## ck56_b0v12_v1   ## ck5/6  krt5
  ##cldn3,  Make 2 binary var:  binarize at: neg vs pos;  strong stain in >20% vs else
  "cldn3_b0v12_v1n", ##cldn3 ## Appears to agree with rates reported in other BrCa studies
  ## "cldn3_b01v2_v1n", ##cldn3  ## Appears lower than rates reported in other BrCa studies - drop
  ##cox2,  Skip - blocks A-J only
  ## cox2_b0v123_v1 <- NULL
  "cryab_b0v12_v1.1n", ## cryab_b0v12_v1.1 ##cryab alpha-basic crystallin (as opposed to alpha-acidic)
  "ecad_b0v12_v1n",    ##ecad_b0v12_v1  ##ecad
  "egfr_b0v12_v1n",    ##egfr_b0v12_v1 ##egfr
  "er_b0v123_v2n",     ##er_b0v123_v2  ##er
  "ezh2_blt10vge10_v1.ppn", ##ezh2_blt10vge10_v1.pp   ##ezh2
  ##foxa1:,  Review intensity and percent - binarize at 0.
  ##, tmadf$foxa1_v1.intf <- factor(tmadf$foxa1_v1.intensity, levels = as.character(0:3))
  ##, tmadf$foxa1_v1.ppn <- as.numeric(tmadf$foxa1_v1.percent)
  ##, foxa1_v1.intensity
  ##, foxa1_v1.percent
  "foxa1_b0v123_v1n",  ##foxa1_v1.intensity  ##foxa1
  "gata3_b0v12_v1n",   ##gata3_b0v12_v1       ##gata3
  "h3k27me3_v1_b012v3n", ##h3k27me3_v1_b012v3   ##h3k27me3
  "her2_b012v23a_v2n",   ##her2_b012v23a_v2     ##her2
  "her3_b0v12_v1.1n",    ##her3_b0v12_v1.1      ##her3
  "her4_b01v2_v1n",      ##her4_b01v2_v1        ##her4
  "hsp27_b0v123_v1n",    ##hsp27_b0v123_v1      ##hsp27
  "igf1r_b0v123_v2.1.1n", ## igf1r_v2.0.1.intensity, igf1r_v2.intensity ##igf1r
  "igfbp2_b01v23_v1n",    ##igfbp2_b01v23_v1               ##igfbp2
  "inpp4b_b0v123_v1",     ## inpp4b_c_v1.pp  inpp4b_v1.int  inpp4b
  "ki67_c1v3_v1.6n",      ##ki67_c1v3_v1.6  ##ki67
  ##,   kitm_c_v1n ## kitm_c_v1 kitm   Ckit mast cells count.  Do scatterplot.
  ##,   kitt_c_v1n ## kitt   Ckit tumour cells count.  Do scatterplot.
  ##kitt_b0v12_v1,                ##kitt   Ckit tumour cells count.  Few positives.  Skip.
  "krt5_b0v12_v1n", ##krt5_b0v12_v1  ##krt5
  "ku7080_c1v1_v1n", ##ku7080_c1v1_v1                 ##ku7080
  "mdm2_b0v12_v1n", ##mdm2_b0v12_v1                  ##mdm2  Blocks A to J only
  "nestin_v1n",     ##nestin_v1   ##nestin 
  "p16_c1v1n",      ##p16
  "p27_b012v3_v1n", ##p27_b012v3_v1                  ##p27
  "p53_c1v1_v1.1.1n",       ## p53_c1v1_v1.1.1 ##p53
  "pcad_ble50vgt50_v1.1.ppn", ##pcad_ble50vgt50_v1.1.pp        ##pcad_c1v1_v1.1
  "pgp_c1v1n",                ## pgp_v1 pgp
  "pipki_b0v12_v1n", ## pipki_b01v2_v1(too few pos cases)   pipki
  "plau_c1v1",      ##plau_v1 ##plau
  "podxl_c1v1",     ##podxl_b01v23_v1 ##podxl
  "pr_b0v123_v3n",  ##pr_b0v123_v3  ##pr
  "psf_c1v1_v2",    ##psf_c1v1_v1 ##psf
  "pten_c1v1",      ##pten_v1 ##pten
  "ret_c1v1",       ##ret_b0v12_v1 ##ret
  "skp2_c1v1",      ##skp2_b0v12_v1 ##skp2
  "trim29_c1v1",    ##trim29_b0v12_v1  ##trim29_b01v2_v1  ##trim29
  "yb1_c1v1.2"     ##yb1_b0v123_v1.2
  )

source("smCrossTable.R")
binchisq <- function(x,
                     xdf = tmadf,
                     agecat = c(20, 40, 55, 70, 100) ) {
  ##browser()
  xdf$agecat <- cut(xdf$age_at_diagnosis, breaks=agecat,
                    labels = c(paste("<=", agecat[2]),
                      paste(agecat[2], "-", agecat[3]),
                      paste(agecat[3], "-", agecat[4]),
                      paste(">", agecat[4])
                      ), right = FALSE,
                    include.lowest = TRUE)
  idxp <- !is.na(xdf[, x])
  
  cat(paste("\n\nTraining set:", x, "\n"))
  ## print(chiout <- chisq.test(xdf[idxp, x], xdf[idxp, "agecat"]))
  ## mtext(text = paste("Proportion", x), side = 2, outer = TRUE)
  ## print(chiout$observed)
  ## print(chiout$expected)
  ctout <- smCrossTable(xdf[idxp, "agecat"], xdf[idxp, x], expected = TRUE, prop.r = TRUE,
                        prop.c = TRUE, prop.t = FALSE, prop.chisq = FALSE, chisq = TRUE,
                        format = "SPSS", digits = 1,
                        dnn = c("Age",
                          paste(toupper(strsplit(x, split = "_")[[1]][1]), "(", x, ")")))
  feout <- if (min(ctout$chisq$expected) < 2) {
    fisher.test(xdf[idxp, "agecat"], xdf[idxp, x])
  } else { NULL }
  if ( !is.null(feout) ) { print(feout) }
  return(list(
    ## chiout = chiout,
              ctout = ctout,
              feout = feout))
  
}

bsctableout <- sapply(bscvars, binchisq, simplify = FALSE, USE.NAMES = TRUE)
vbsctableout <- sapply(bscvars, binchisq, xdf = vmadf, simplify = FALSE, USE.NAMES = TRUE)
wbscpvals <- c(sapply(bsctableout, function(x) x[[1]]$chisq$p.value),
               sapply(vbsctableout, function(x) x[[1]]$chisq$p.value))
wbscpvalsadj <- p.adjust(wbscpvals, method = "BH")
bsctableoutadj <-
  lapply(seq(length(bsctableout)),
         function(x) {
           bsctableout[[x]]$ctout$chisq$p.adjust <- wbscpvalsadj[x]
           bsctableout[[x]]
         })
names(bsctableoutadj) <- names(bsctableout)
vbsctableoutadj <-
  lapply(seq(length(vbsctableout)),
         function(x) {
           vbsctableout[[x]]$ctout$chisq$p.adjust <- wbscpvalsadj[length(bsctableout) + x]
           vbsctableout[[x]]
         })
names(vbsctableoutadj) <- names(vbsctableout)



## Scatterplots with smooths

binsupsmu <- function(x, xdf = tmadf, ydf = NULL,
                      binchisqlist = NULL, ybinchisqlist = NULL,
                      tlbl = "") {

  idxp <- !is.na(xdf[, x])
  if ( length( unique(xdf[idxp, x]) ) == 2 ) {
    plot(xdf[idxp, "age_at_diagnosis"], jitter(xdf[idxp, x], factor = 0.1),
         ylab = paste("Proportion", x), xlab = "Age (years)", xlim = xlims, ylim = ylims)
    lines(supsmu(xdf[idxp, "age_at_diagnosis"], xdf[idxp, x], bass = 3), lwd = 3)
    abline(h = 0.5, v = 50, lty = 2)
    xdfttl <- "Training set"
    if (tlbl != "") xdfttl <- paste(xdfttl, "(", tlbl, ")" )
    title(main = xdfttl)
    if ( !is.null(binchisqlist) && !is.null(binchisqlist[[x]]) ) {
      ypos <- ifelse(binchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
      xpos <- 75
      text(xpos, ypos, paste(toupper(strsplit(x, split = "_")[[1]][1]),
                             "by AgeCat: p =",
                             format(binchisqlist[[x]]$ctout$chisq$p.value, digits = 2)))
      if ( !is.null( binchisqlist[[x]]$ctout$chisq$p.adjust ) ) {
        ypos <- ypos - 0.06
        xpos <- xpos + 7
        text(xpos, ypos, paste("p.adjust =",
                               format(binchisqlist[[x]]$ctout$chisq$p.adjust, digits = 2)))
      }
    }
  } else {
    plot(0, 0, type = "n", axes = FALSE, frame.plot = TRUE, xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = "", ylab = "")
    text(-0.2, 0, labels = paste("Data only available for level ",
                    unique(xdf[idxp, x]), "\nfor variable", x ) )    
  }
  if ( !is.null(ydf) ) {
    idxp <- !is.na(ydf[, x])
    if ( length( unique(ydf[idxp, x]) ) == 2 ) {
      plot(ydf[idxp, "age_at_diagnosis"], jitter(ydf[idxp, x], factor = 0.1),
           ylab = paste("Proportion", x), xlab = "Age (years)", xlim = xlims, ylim = ylims)
      lines(supsmu(ydf[idxp, "age_at_diagnosis"], ydf[idxp, x], bass = 3), lwd = 3)
      abline(h = 0.5, v = 50, lty = 2)
      ydfttl <- "Validation set"
      if (tlbl != "") ydfttl <- paste(ydfttl, "(", tlbl, ")" )
      title(main = ydfttl)
      if ( !is.null(ybinchisqlist) && !is.null(ybinchisqlist[[x]]) ) {
        ypos <- ifelse(ybinchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
        xpos <- 75
        text(xpos, ypos, paste(toupper(strsplit(x, split = "_")[[1]][1]),
                               "by AgeCat: p =",
                               format(ybinchisqlist[[x]]$ctout$chisq$p.value, digits = 2)))
        if ( !is.null( ybinchisqlist[[x]]$ctout$chisq$p.adjust ) ) {
          ypos <- ypos - 0.06
          xpos <- xpos + 7
          text(xpos, ypos, paste("p.adjust =",
                                 format(ybinchisqlist[[x]]$ctout$chisq$p.adjust, digits = 2)))
        }
      }
    } else {
      plot(0, 0, type = "n", axes = FALSE, frame.plot = TRUE, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "", ylab = "")
      text(-0.2, 0, labels = paste("Data only available for level ",
                      unique(xdf[idxp, x]), "\nfor variable", x ) )    
    }
  }
  ## mtext(text = paste("Proportion", x), side = 2, outer = TRUE) 
}

binsupsmu <- function(x, xdf = tmadf, ydf = NULL,
                      binchisqlist = NULL, ybinchisqlist = NULL,
                      tlbl = "") {

  idxp <- !is.na(xdf[, x])
  if ( length( unique(xdf[idxp, x]) ) == 2 ) {
    plot(xdf[idxp, "age_at_diagnosis"], jitter(xdf[idxp, x], factor = 0.1),
         ylab = paste("Proportion positive"), xlab = "Age (years)", xlim = xlims, ylim = ylims)
    lines(supsmu(xdf[idxp, "age_at_diagnosis"], xdf[idxp, x], bass = 3), lwd = 3)
    abline(h = 0.5, v = 50, lty = 2)
    xdfttl <- "Training set"
    #if (tlbl != "") xdfttl <- paste(xdfttl, "(", tlbl, ")" )
    title(main = xdfttl)
    if ( !is.null(binchisqlist) && !is.null(binchisqlist[[x]]) ) {
      ypos <- ifelse(binchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
      xpos <- 75
      text(xpos, ypos, paste(#toupper(strsplit(x, split = "_")[[1]][1]),
                             "p =",
                             format(binchisqlist[[x]]$ctout$chisq$p.value, digits = 2) ), cex=0.9 )
      if ( !is.null( binchisqlist[[x]]$ctout$chisq$p.adjust ) ) {
        ypos <- ypos - 0.06
        xpos <- xpos + 7
        text(xpos, ypos, paste("p.adjust =",
                               format(binchisqlist[[x]]$ctout$chisq$p.adjust, digits = 2)))
      }
    }
  } else {
    plot(0, 0, type = "n", axes = FALSE, frame.plot = TRUE, xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = "", ylab = "")
    text(-0.2, 0, labels = paste("Data only available for level ",
                    unique(xdf[idxp, x]), "\nfor variable", x ) )    
  }
  if ( !is.null(ydf) ) {
    idxp <- !is.na(ydf[, x])
    if ( length( unique(ydf[idxp, x]) ) == 2 ) {
      plot(ydf[idxp, "age_at_diagnosis"], jitter(ydf[idxp, x], factor = 0.1),
           ylab = paste("Proportion positive"), xlab = "Age (years)", xlim = xlims, ylim = ylims)
      lines(supsmu(ydf[idxp, "age_at_diagnosis"], ydf[idxp, x], bass = 3), lwd = 3)
      abline(h = 0.5, v = 50, lty = 2)
      ydfttl <- "Validation set"
      #if (tlbl != "") ydfttl <- paste(ydfttl, "(", tlbl, ")" )
      title(main = ydfttl)
      if ( !is.null(ybinchisqlist) && !is.null(ybinchisqlist[[x]]) ) {
        ypos <- ifelse(ybinchisqlist[[x]]$ctout$prop.row[4, 2] < 0.5, 0.75, 0.25)
        xpos <- 75
        text(xpos, ypos, paste(#toupper(strsplit(x, split = "_")[[1]][1]),
                               "p =",
                               format(ybinchisqlist[[x]]$ctout$chisq$p.value, digits = 2) ), cex = 0.9 )
        if ( !is.null( ybinchisqlist[[x]]$ctout$chisq$p.adjust ) ) {
          ypos <- ypos - 0.06
          xpos <- xpos + 7
          text(xpos, ypos, paste("p.adjust =",
                                 format(ybinchisqlist[[x]]$ctout$chisq$p.adjust, digits = 2)))
        }
      }
    } else {
      plot(0, 0, type = "n", axes = FALSE, frame.plot = TRUE, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "", ylab = "")
      text(-0.2, 0, labels = paste("Data only available for level ",
                      unique(xdf[idxp, x]), "\nfor variable", x ) )    
    }
  }
  ## mtext(text = paste("Proportion", x), side = 2, outer = TRUE) 
}

binsupsmu(bscvars[1], binchisqlist = bsctableout)

pdf("./Plots/BiomarkerRatesByAgev04.pdf", width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 0, 0), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu, binchisqlist = bsctableout )
dev.off()


pdf("./Plots/tvBiomarkerRatesByAgev05.pdf", width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 0, 0), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu, ydf = vmadf,
       binchisqlist = bsctableoutadj, ybinchisqlist = vbsctableoutadj,
       tlbl = "All subtypes")
dev.off()

### Luminalp
BrCa4Subtype <- "Luminalp"
lump_t_tblout <- sapply(bscvars, binchisq, xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)

lump_v_tblout <- sapply(bscvars, binchisq, xdf = vmadf[vmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)

pdf(paste("./Plots/tv_", BrCa4Subtype, "_BiomarkerRatesByAgev05.pdf", sep = ""), width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 2, 2), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu,
       xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
       ydf = vmadf[vmadf$BrCa4 == BrCa4Subtype, ],
       binchisqlist = lump_t_tblout, ybinchisqlist = lump_v_tblout,
       tlbl = BrCa4Subtype )
dev.off()



### TNP
BrCa4Subtype <- "TNP"
lump_t_tblout <- sapply(bscvars, binchisq, xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)

lump_v_tblout <- sapply(bscvars, binchisq, xdf = vmadf[vmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)

pdf(paste("./Plots/tv_", BrCa4Subtype, "_BiomarkerRatesByAgev01.pdf", sep = ""), width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 2, 2), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu,
       xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
       ydf = vmadf[vmadf$BrCa4 == BrCa4Subtype, ],
       binchisqlist = lump_t_tblout, ybinchisqlist = lump_v_tblout,
       tlbl = BrCa4Subtype )
dev.off()


### HER2+       "Luminal/HER2+" "HER2+/ER-PR-"
BrCa4Subtype <- "HER2+"
lump_t_tblout <- sapply(bscvars, binchisq,
                        xdf = tmadf[((tmadf$BrCa4 == "Luminal/HER2+") |
                                     (tmadf$BrCa4 == "HER2+/ER-PR-") ), ],
                        simplify = FALSE, USE.NAMES = TRUE)

lump_v_tblout <- sapply(bscvars, binchisq,
                        xdf = vmadf[((vmadf$BrCa4 == "Luminal/HER2+") |
                                     (vmadf$BrCa4 == "HER2+/ER-PR-") ), ],
                        simplify = FALSE, USE.NAMES = TRUE)

pdf(paste("./Plots/tv_", BrCa4Subtype, "_BiomarkerRatesByAgev01.pdf", sep = ""), width = 8, height = 10)
par(mfrow = c(3, 2), oma = c(2, 2, 2, 2), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu,
       xdf = tmadf[((tmadf$BrCa4 == "Luminal/HER2+") |
                    (tmadf$BrCa4 == "HER2+/ER-PR-") ), ],
       ydf = vmadf[((vmadf$BrCa4 == "Luminal/HER2+") |
                    (vmadf$BrCa4 == "HER2+/ER-PR-") ), ],
       binchisqlist = lump_t_tblout, ybinchisqlist = lump_v_tblout,
       tlbl = BrCa4Subtype )
dev.off()

### Need to adjust the p-values for the tests of association of biomarker
### with age.  40 biomarkers, training and validation, so adjust the 80 p-values for multiple comparisons.


### Need beanplots boxplots of age by subtype

require("beanplot")
win.graph()
beanplot(age_at_diagnosis ~ BrCa4f, data = tmadf,
         col = list("purple", "cyan", "orange", "red", "green"))
binchisq("BrCa4f")
### $ctout$prop.col
###          y
### x           Luminalp Luminal/HER2+ HER2+/ER-PR-        TNP Unassigned
###   <= 40   0.05299539    0.08264463   0.13333333 0.15723270 0.04225352
###   40 - 55 0.27956989    0.34710744   0.35000000 0.36792453 0.23943662
###   55 - 70 0.37403994    0.38842975   0.37500000 0.30188679 0.40845070
###   > 70    0.29339478    0.18181818   0.14166667 0.17295597 0.30985915

binchisq("BrCa4f", xdf = vmadf)
### $ctout$prop.col
###          y
### x           Luminalp Luminal/HER2+ HER2+/ER-PR-        TNP Unassigned
###   <= 40   0.04887510    0.06796117   0.12307692 0.13782051 0.09090909
###   40 - 55 0.30178433    0.35922330   0.28461538 0.37820513 0.27272727
###   55 - 70 0.38789760    0.35922330   0.44615385 0.31410256 0.38961039
###   > 70    0.26144298    0.21359223   0.14615385 0.16987179 0.24675325


### Grant figs
### "ezh2_blt10vge10_v1.ppn"
### Luminalp
BrCa4Subtype <- "TNP" #"Luminalp"
bscvars <- "p16_c1v1n" #"ki67_c1v3_v1.6n" #"foxa1_b0v123_v1n" #"ezh2_blt10vge10_v1.ppn"
#bscvars <- "h3k27me3_v1_b012v3n" 
lump_t_tblout <- sapply(bscvars, binchisq, xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)

lump_v_tblout <- sapply(bscvars, binchisq, xdf = vmadf[vmadf$BrCa4 == BrCa4Subtype, ],
                        simplify = FALSE, USE.NAMES = TRUE)



xlims <- c(20, 100)
ylims <- c(-0.05, 1.05)

win.metafile(filename = paste("./Plots/tv_", bscvars, "_", BrCa4Subtype, "_BiomarkerRatesByAgev06.wmf",
               sep = ""), width = 7, height = 4)
par(mfrow = c(1, 2), oma = c(1,1,1,1), mar = c(4, 5, 4, 1))
set.seed(747)
sapply(bscvars, binsupsmu,
       xdf = tmadf[tmadf$BrCa4 == BrCa4Subtype, ],
       ydf = vmadf[vmadf$BrCa4 == BrCa4Subtype, ],
       binchisqlist = lump_t_tblout, ybinchisqlist = lump_v_tblout,
       tlbl = BrCa4Subtype )
dev.off()



### Continuous plots
simsmufit <- function(x, xdf = tmadf, ydf = NULL, xvar = "age_at_diagnosis",
                      tlbl = "", nsim = 2000, alpha = 0.05) {
  rn <- nsim
  alphalev <- alpha/2.0
  UCLn <- trunc(rn * (1 - alphalev))
  LCLn <- trunc(rn * alphalev)
  agesv <- sort(unique(xdf[, xvar]), na.last = NA)
  simmat <- matrix(NA_real_, nrow = rn, ncol = length(agesv))
  ylims <- range(xdf[, x], na.rm = TRUE)
  ylimd <- abs(diff(ylims))
###   ylimlo <- ylims[1] - (0.5 * ylimd)
  ylimlo <- ylims[1] - (0.2 * ylimd)
  ylimhi <- ylims[2] + (0.2 * ylimd)
  valdx <- !is.na(xdf[, x])
###   xlb <- strsplit(xvar, split = "_")[[1]][1]
###   xlbc <- paste(toupper(substring(xlb, 1, 1)), tolower(substring(xlb, 2)), sep = "")
  xlbc <- "Dx age (years)"
  ylbc <- toupper(strsplit(x, split = "_")[[1]][1])
  plot( xdf[valdx, xvar], xdf[valdx, x], xlab = xlbc, ylab = ylbc,
       xlim = c(15, 95), ylim = c(ylimlo, ylimhi), pch = ".", col = "grey60" )
###  axis(side = 2, at = c(0, 20, 40, 60, 80, 100))
  for (i in seq(rn) ) {
    ssout <- supsmu( sample(xdf[valdx, xvar]), xdf[valdx, x], bass = 10)
    spout <- predict(smooth.spline(ssout$x, ssout$y), x = agesv)
    simmat[i, ] <- spout$y
###     lines(ssout, lty = 3, lwd = 1, col = "lightgrey")
  }
  ##points( xdf[, xvar], xdf[, x])
  srtmat <- apply(simmat, 2, sort)
###   lines(supsmu(agesv, srtmat[LCLn, ], bass = 10), lty = 1, lwd = 3, type = "l", col = "grey80")
###   lines(supsmu(agesv, srtmat[UCLn, ], bass = 10), lty = 1, lwd = 3, type = "l", col = "grey80")
  polygon(c(agesv, rev(agesv)), c(srtmat[LCLn, ], rev(srtmat[UCLn, ])), col = "wheat3", border = "wheat3")
  ssout <- supsmu( xdf[valdx, xvar], xdf[valdx, x], bass = 10)
  lines(ssout, lwd = 2, type = "l")
  srtmatcols <- match(ssout$x, agesv)
  if ( length(ssout$x) <= length(agesv) ) {
    if ( any( ssout$y < srtmat[LCLn, srtmatcols] ) |
         any( ssout$y > srtmat[UCLn, srtmatcols] ) ) {
      idxs <- c( which( ssout$y < srtmat[LCLn, srtmatcols] ),
                 which( ssout$y > srtmat[UCLn, srtmatcols] ) )
      points(ssout$x[idxs], ssout$y[idxs], pch = 20)
      text(x = 20, y = ylimhi, labels = paste("Age-dependent p<", alpha, sep = ""), adj = 0)
    } else {
      text(x = 20, y = ylimhi, labels = paste("Not age-dependent p>", alpha, sep = ""), adj = 0)
    }
  } else {
    text(x = 20, y = ylimhi, labels = "FIX ME missing ages", adj = 0)
  }
###   legend("bottomright", legend = paste(trunc(100*(1-alpha)), "% Confidence Bounds", sep = ""),
###          lty = 1, lwd = 4, seg.len = 3, text.width = 60, col = "wheat3")
  if ( !is.null(ydf) ) {
    agesv <- sort(unique(ydf[, xvar]), na.last = NA)
    simmat <- matrix(NA_real_, nrow = rn, ncol = length(agesv))
    ylims <- range(ydf[, x], na.rm = TRUE)
    ylimd <- abs(diff(ylims))
    ylimlo <- ylims[1] - (0.2 * ylimd)
    ylimhi <- ylims[2] + (0.2 * ylimd)
    valdy <- !is.na(ydf[, x])
    plot( ydf[valdy, xvar], ydf[valdy, x], xlab = xvar, ylab = x,
         xlim = c(15, 95), ylim = c(ylimlo, ylimhi), pch = ".", col = "grey60" )
###  axis(side = 2, at = c(0, 20, 40, 60, 80, 100))
    for (i in seq(rn) ) {
      ssout <- supsmu( sample(ydf[valdy, xvar]), ydf[valdy, x], bass = 10)
      spout <- predict(smooth.spline(ssout$x, ssout$y), x = agesv)
      simmat[i, ] <- spout$y
###      lines(ssout, lty = 3, lwd = 1, col = "lightgrey")
    }
    ##points( ydf[, xvar], ydf[, x])
    srtmat <- apply(simmat, 2, sort)
    polygon(c(agesv, rev(agesv)), c(srtmat[LCLn, ], rev(srtmat[UCLn, ])), col = "wheat3", border = "wheat3")
###     lines(agesv, srtmat[LCLn, ], lty = 1, lwd = 3, type = "l", col = "grey80")
###     lines(agesv, srtmat[UCLn, ], lty = 1, lwd = 3, type = "l", col = "grey80")
    lines(supsmu( ydf[valdy, xvar], ydf[valdy, x], bass = 10), lwd = 3, type = "l")
###     legend("bottomright", legend = paste(trunc(100*(1-alpha)), "% Confidence Bounds", sep = ""),
###            lty = 1, lwd = 4, seg.len = 3, text.width = 60, col = "wheat3")
  }
}

cscvars <- c(
##  "ar_c_v1nc", ## Under study by another team, leave out.
  "bcl2_c_v2.ppnc",
  "ca9_b0v123_v1n",
##  "cd163_c_v1nc", ## Macrophage - not so relevant?
  "ck56_v1nc",
  "cldn3_v1nc",
  "cryab_b0v12_v1.1n",
### ADD IN Cyclin e CCNE1
  "cycline_v1.1.ppn",
  "ecad_v1nc",
  "egfr_v1nc",
  "ephb4_v1.1nc",
  "er_v2nc",
  "er_dcc_v1n",
  "ezh2_v1.ppnc",
  "foxa1_v1.ppnc",
  "gata3_v1nc",
  "h3k27me3_v1.1.ppnc",
  "her2_v2nc",
  "her3_v1.1nc",
  "her4_v1nc",
  "hsp27_v1nc",
  "igf1r_v2.0.1.ppn",
  "igfbp2_v1nc",
  "inpp4b_c_v1.ppn",
  "ki67_c_v1.6nc",
##  "kitm_c_v1n",
##  "kitt_v1nc",
  "kit_v1nc",
  "krt5_v1nc",
  "ku7080_c_v1nc",
  "mdm2_v1nc",
  "nestin_v1nc",
  "p16_c1v1nc",
  "p27_v1nc",
  "p53_c_v1.1.1.pnnc",
  "pcad_v1.1.ppnc",
  "pgp_v1nc",
  "pipki_v1nc",
  "plau_c1v1",
  "podxl_v1nc",
  "pr_v3nc",
  "psf_c_v1nc",
  "pten_v1nc",
  "ret_v1nc",
  "skp2_v1nc",
  "trim29_v1nc",
  "yb1_v1.2nc"
  )


cscvars4 <- c(
  "er_v2nc",
  "ezh2_v1.ppnc",
  "foxa1_v1.ppnc",
  "h3k27me3_v1.1.ppnc"
  )


### Remove MB09 cases
tmadfbackupall <- tmadf
vmadfbackupall <- vmadf

### tmadf <- tmadf[!(tmadf$bseries_id %in% mb09matchdf$bseries_id), ]
### vmadf <- vmadf[!(vmadf$bseries_id %in% mb09matchdf$bseries_id), ]
tmadf <- tmadf[tmadf$bseries_id %in% match3df$bseries_id, ]
vmadf <- vmadf[vmadf$bseries_id %in% match3df$bseries_id, ]

### sapply("ephb4_v1.1nc", simsmufit, ydf = vmadf)
pdf(file = "BigSeries_NoMB09_NoMBex_IHC_Trends_95pct_SimCIs_shaded_v05.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
sapply(cscvars, simsmufit, ydf = vmadf)
dev.off()

### Whole series - Need a p-value declaration on each plot.
date()
pdf(file = "BigSeries_WholeSeries_NoMB09_NoMBex_IHC_Trends_95pct_2kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = wmadf, nsim = 2000, alpha = 0.05)
dev.off()
date()
date()
pdf(file = "BigSeries_WholeSeries_NoMB09_NoMBex_IHC_Trends_99pct_10kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = wmadf, nsim = 10000, alpha = 0.01)
dev.off()
date()
date()
pdf(file = "BigSeries_WholeSeries_NoMB09_NoMBex_IHC_Trends_99p9pct_40kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = wmadf, nsim = 400000, alpha = 0.001)
dev.off()
date()
date()
pdf(file = "BigSeries_WholeSeries_NoMB09_NoMBex_IHC_Trends_99p99pct_400kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = wmadf, nsim = 400000, alpha = 0.0001)
dev.off()
date()


### Whole series - Triple negatives - Need a p-value declaration on each plot.
tnmadf <- wmadf[wmadf$BrCa4 == "TNP", ]
lpmadf <- wmadf[wmadf$BrCa4 == "Luminalp", ]
date()
pdf(file = "BigSeries_TNP_WholeSeries_NoMB09_NoMBex_IHC_Trends_99pct_10kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = tnmadf, nsim = 10000, alpha = 0.01)
dev.off()
date()
date()
pdf(file = "BigSeries_TNP_WholeSeries_NoMB09_NoMBex_IHC_Trends_99p99pct_400kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = tnmadf, nsim = 400000, alpha = 0.0001)
dev.off()
date()

date()
pdf(file = "BigSeries_Luminalp_WholeSeries_NoMB09_NoMBex_IHC_Trends_99pct_10kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = lpmadf, nsim = 10000, alpha = 0.01)
dev.off()
date()
date()
pdf(file = "BigSeries_Luminalp_WholeSeries_NoMB09_NoMBex_IHC_Trends_99p99pct_400kSimCIs_shaded_v01.pdf",
    width = 8, height = 8)
par(mfrow = c(2, 2))
set.seed(473); sapply(cscvars, simsmufit, xdf = lpmadf, nsim = 400000, alpha = 0.0001)
dev.off()
date()


pdf(file = "BigSeries_NoMB09_NoMBex_IHC_Trends_95pct_SimCIs_4Tgts_shaded_v01.pdf", width = 8, height = 8)
par(mfrow = c(2, 2), oma = c(1,1,3,1))
mtext("BigSeries", outer=TRUE)
sapply(cscvars4, simsmufit, ydf = vmadf)
mtext("BigSeries", outer=TRUE)
dev.off()

NULL
## Is there any age related trend to missing data?  Is missing data pattern contributing?


  ar_c_v1nc
  bcl2_c_v2.ppnc
  ca9_b0v123_v1n
  cd163_c_v1nc
  ck56_v1nc
  cldn3_v1nc
  cryab_b0v12_v1.1n
  ecad_v1nc
  egfr_v1nc
  er_v2nc
  er_dcc_v1n
  ezh2_v1.ppnc
  foxa1_v1.ppnc
  gata3_v1nc
  h3k27me3_v1.1.ppnc
  her2_v2nc
  her3_v1.1nc
  her4_v1nc
  hsp27_v1nc
  igf1r_v2.0.1.ppn
  igfbp2_v1nc
  inpp4b_c_v1.ppn
  ki67_c_v1.6nc
  kitm_c_v1n
  kitt_v1nc
  krt5_v1nc
  ku7080_c_v1nc
  mdm2_v1nc
  nestin_v1nc
  p16_c1v1nc
  p27_v1nc
  p53_c_v1.1.1.pnnc
  pcad_v1.1.ppnc
  pgp_v1nc
  pipki_v1nc
  plau_c1v1
  podxl_v1nc
  pr_v3nc
  psf_c_v1nc
  pten_v1nc
  ret_v1nc
  skp2_v1nc
  trim29_v1nc
  yb1_v1.2nc


### Patient characteristics
sddf <-  rbind(tmadf, vmadf)
quantile(sddf$age_at_diagnosis, na.rm = TRUE)
sddf$Agef <- cut(sddf$age_at_diagnosis, breaks = c(0, 30, 40, 50, 60, 70, 80, 100), right = FALSE)
table(sddf$Agef, useNA = "always")
table(sddf$Agef, useNA = "always")/sum(table(sddf$Agef, useNA = "always"))*100

### BrCa subtype

smCrossTable(x = sddf$pam50_subtype, y = sddf$BrCa4)
smCrossTable(y = sddf$pam50_subtype, x = sddf$BrCa4)

table(sddf$BrCa4, useNA = "always")
table(sddf$BrCa4, useNA = "always")/sum(table(sddf$BrCa4, useNA = "always"))*100

### HER2 status

table(sddf$her2_b012v23a_v2, useNA = "always")
table(sddf$her2_b012v23a_v2, useNA = "always")/sum(table(sddf$her2_b012v23a_v2, useNA = "always"))*100

### Clinical T N M

table(sddf$tnm_clin_t_4cat, useNA = "always")
table(sddf$tnm_clin_t_4cat, useNA = "always")/sum(table(sddf$tnm_clin_t_4cat, useNA = "always"))*100

table(sddf$tnm_clin_n, useNA = "always")
table(sddf$tnm_clin_n, useNA = "always")/sum(table(sddf$tnm_clin_n, useNA = "always"))*100

table(sddf$tnm_clin_m, useNA = "always")
table(sddf$tnm_clin_m, useNA = "always")/sum(table(sddf$tnm_clin_m, useNA = "always"))*100

### Pathological T N M

table(sddf$tnm_surg_t_4cat, useNA = "always")
table(sddf$tnm_surg_t_4cat, useNA = "always")/sum(table(sddf$tnm_surg_t_4cat, useNA = "always"))*100

table(sddf$tnm_surg_n, useNA = "always")
table(sddf$tnm_surg_n, useNA = "always")/sum(table(sddf$tnm_surg_n, useNA = "always"))*100

table(sddf$tnm_surg_m, useNA = "always")
table(sddf$tnm_surg_m, useNA = "always")/sum(table(sddf$tnm_surg_m, useNA = "always"))*100

## Tx

[143] "ad"                                    
[144] "partial"                               
[145] "complete"                              
[146] "finsurg"                               
[147] "have_complete"                         
[148] "localtx"                               
[149] "rtintent"                              
[150] "brchwrt"                               
[151] "nodalrt"                               
[152] "boost"                                 
[153] "finrt"                                 
[154] "rt"                                    
[155] "behavior"                              
[156] "chemflag"                              
[157] "any_chem"                              
[158] "chemtype"                              
[159] "chemtype_cat2"                         
[160] "init_chemo"                            
[161] "hormflag"                              
[162] "hormtype"                              
[163] "any_horm"                              
[164] "systemic"      
table(sddf$finsurg, useNA = "always")
table(sddf$localtx, useNA = "always")
table(sddf$localtx, useNA = "always")/sum(table(sddf$localtx, useNA = "always"))*100

table(sddf$systemic, useNA = "always")
table(sddf$systemic, useNA = "always")/sum(table(sddf$systemic, useNA = "always"))*100

### EZH2 and H3K27me3 association:

win.graph()
ccidxp <- !( is.na(sddf$h3k27me3_v1.1.ppnc) | is.na(sddf$ezh2_v1.ppnc) )
plot(sddf$h3k27me3_v1.1.ppnc[ccidxp], sddf$ezh2_v1.ppnc[ccidxp])
lines( supsmu(sddf$h3k27me3_v1.1.ppnc[ccidxp], sddf$ezh2_v1.ppnc[ccidxp],
              span = 0.4, bass = 10), lwd = 3)
cor(sddf$h3k27me3_v1.1.ppnc[ccidxp], sddf$ezh2_v1.ppnc[ccidxp]) # [1] -0.1161101
anova(lm(sddf$ezh2_v1.ppnc[ccidxp] ~ sddf$h3k27me3_v1.1.ppnc[ccidxp])) # p = 2.462e-09
dev.print(device = png, file = "EZH2vH3K27me3_PctPos_v01.png",
          width = 6, height = 6, units = "in", res = 600)

win.graph()
par(mfrow = c(3, 2))
plot(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminalp")[ccidxp]],
     sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Luminalp")[ccidxp]],
     xlab = "H3K27me3", ylab = "EZH2", main = "Luminalp")
lines( supsmu(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminalp")[ccidxp]],
              sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Luminalp")[ccidxp]],
              span = 0.4, bass = 10), lwd = 3)
cor(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminalp")[ccidxp]],
    sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Luminalp")[ccidxp]]) # -0.04526539
anova(lm(sddf$ezh2_v1.ppnc[ccidxp] [(sddf$BrCa5 == "Luminalp")[ccidxp]] ~
           sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminalp")[ccidxp]])) # p =0.05384
plot(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]],
     sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]],
     xlab = "H3K27me3", ylab = "EZH2", main = "Luminal/HER2+")
lines( supsmu(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]],
              sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]],
              span = 0.4, bass = 10), lwd = 3)
cor(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]],
    sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]]) # 
anova(lm(sddf$ezh2_v1.ppnc[ccidxp] [(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]] ~
           sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Luminal/HER2+")[ccidxp]])) # 
plot(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]],
     sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]],
     xlab = "H3K27me3", ylab = "EZH2", main = "HER2+/ER-PR-")
lines( supsmu(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]],
              sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]],
              span = 0.4, bass = 10), lwd = 3)
cor(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]],
    sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]]) # 
anova(lm(sddf$ezh2_v1.ppnc[ccidxp] [(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]] ~
           sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "HER2+/ER-PR-")[ccidxp]])) # 
plot(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Core Basal")[ccidxp]],
     sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Core Basal")[ccidxp]],
     xlab = "H3K27me3", ylab = "EZH2", main = "Core Basal")
lines( supsmu(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Core Basal")[ccidxp]],
              sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Core Basal")[ccidxp]],
              span = 0.4, bass = 10), lwd = 3)
cor(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Core Basal")[ccidxp]],
    sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "Core Basal")[ccidxp]]) # 
anova(lm(sddf$ezh2_v1.ppnc[ccidxp] [(sddf$BrCa5 == "Core Basal")[ccidxp]] ~
           sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "Core Basal")[ccidxp]])) # 
plot(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "5NP")[ccidxp]],
     sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "5NP")[ccidxp]],
     xlab = "H3K27me3", ylab = "EZH2", main = "5NP")
lines( supsmu(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "5NP")[ccidxp]],
              sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "5NP")[ccidxp]],
              span = 0.4, bass = 10), lwd = 3)
cor(sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "5NP")[ccidxp]],
    sddf$ezh2_v1.ppnc[ccidxp][(sddf$BrCa5 == "5NP")[ccidxp]]) # 
anova(lm(sddf$ezh2_v1.ppnc[ccidxp] [(sddf$BrCa5 == "5NP")[ccidxp]] ~
           sddf$h3k27me3_v1.1.ppnc[ccidxp][(sddf$BrCa5 == "5NP")[ccidxp]])) # 


dev.print(device = png, file = "EZH2vH3K27me3_PctPos_BrCa5_v01.png",
          width = 6, height = 6, units = "in", res = 600)


### AIC for AR cutpoint
