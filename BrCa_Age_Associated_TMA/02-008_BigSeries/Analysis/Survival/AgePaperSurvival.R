## Survival analysis for age association paper
## Assess EZH2 impact on survival adjusting for ER, age, grade, tumour size, node status
##


library(here)
library(fs)
library(gmodels)  ## for e.g. CrossTable()
library(tidyverse)
library(survival)
### Make directory to hold plots
plot_dirs <- path(here("02-008_BigSeries/Analysis/Survival/Plots"))
dir_create(plot_dirs)

## Utility function to extract p-value from a survdiff output object (sdifobj)
lrpv <- function(sdifobj){
  sdifdf <- length(sdifobj$obs) - 1
  pval <- pchisq(sdifobj$chisq, df = sdifdf, lower.tail = FALSE)
}


wmadf <- read.table("02-008_BigSeries/Analysis/BuildData/02-008_BigSeries_Whole.CSV",
                    stringsAsFactors = FALSE, sep = ',', header = TRUE, quote = "\"",
                    na.strings = " ", comment.char = "", colClasses = "character")
names(wmadf)[1] <- "bseries_id"
wmatb <- as.tibble(wmadf)

wmadf$age <- as.numeric(wmadf$age_at_diagnosis)
table(wmadf$age, useNA = "always")
quantile(wmadf$age)
table(wmadf$dxage_g)
wmadf$ageg <- cut(wmadf$age, breaks = c(0, 45, 60, 70, 99))
table(wmadf$ageg)
with(wmadf, CrossTable(ageg, er_b0v123_v2))

table(wmadf$ezh2_blt10vge10_v1.pp, useNA = "always")
with(wmadf, CrossTable(ageg, er_b0v123_v2, ezh2_blt10vge10_v1.pp))

with(wmadf, table(ageg, er_b0v123_v2, ezh2_blt10vge10_v1.pp, useNA = "always"))
with(wmadf, table(ageg, er_b0v123_v2, ezh2_blt10vge10_v1.pp))

#leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml)
#levels(wmadf$ageg)

is.na(wmadf$er_b0v123_v2) <- 
  wmadf$er_b0v123_v2 == "uninterpretable/missing"
head(wmadf$survyrs1999 )
head(wmadf$survstat1999 )
wmadf$survyrs1999 <- as.numeric(wmadf$survyrs1999)
wmadf$survstatc1999 <- wmadf$survstat1999
wmadf$survstat1999 <- rep(NA_integer_, nrow(wmadf))
wmadf$survstat1999[wmadf$survstatc1999 == "Alive"] <- 0
wmadf$survstat1999[wmadf$survstatc1999 == "Dead"] <- 1
table(wmadf$survstat1999 )
with(wmadf, table(er_b0v123_v2, erposneg))
a1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(0,45]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
a2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(45,60]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
a3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(60,70]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
a4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(70,99]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))

n1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(0,45]" & er_b0v123_v2 == "negative" ))
n2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(45,60]" & er_b0v123_v2 == "negative" ))
n3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(60,70]" & er_b0v123_v2 == "negative" ))
n4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                   subset = (ageg == "(70,99]" & er_b0v123_v2 == "negative" ))
par(mfcol=c(4,2))
plot(a1.surv, lty = 2:3)
plot(a2.surv, lty = 2:3)
plot(a3.surv, lty = 2:3)
plot(a4.surv, lty = 2:3)
plot(n1.surv, lty = 2:3)
plot(n2.surv, lty = 2:3)
plot(n3.surv, lty = 2:3)
plot(n4.surv, lty = 2:3)

par(mfrow=c(2,4))
plot(a1.surv, lty = 2:3)
plot(a2.surv, lty = 2:3)
plot(a3.surv, lty = 2:3)
plot(a4.surv, lty = 2:3)
plot(n1.surv, lty = 2:3)
plot(n2.surv, lty = 2:3)
plot(n3.surv, lty = 2:3)
plot(n4.surv, lty = 2:3)
legend("topright",legend = levels(as.factor(wmadf$ezh2_blt10vge10_v1.pp))  , lty = 2:3)
title(main = "EZH2 low vs high")

ccs <- with(wmadf, complete.cases(ageg, ezh2_blt10vge10_v1.pp, er_b0v123_v2))
cmff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg * er_b0v123_v2, data = wmadf[ccs, ])
cmff
cmfez <-  coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg , data = wmadf[ccs, ])
cmfez
anova(cmfez, cmff)
cmfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2 , data = wmadf[ccs, ])
cmfer
anova(cmfer, cmff)

is.na(wmadf$grade_b12v3) <- wmadf$grade_b12v3 == "unknown"
is.na(wmadf$nodestat) <- wmadf$nodestat == "nodal status unknown"
is.na(wmadf$size_tumor_grp) <- wmadf$size_tumor_grp == "unknown"

table(wmadf$grade_b12v3)
table(wmadf$nodestat)
table(wmadf$size_tumor_grp)


ccsp <- with(wmadf, complete.cases(ageg, ezh2_blt10vge10_v1.pp, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
cmpff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
               data = wmadf[ccsp, ])
cmpff
cmpfez <-  coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccsp, ])
cmpfez
anova(cmpfez, cmpff)
cmpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccsp, ])
cmpfer
anova(cmpfer, cmpff)

### HER2 in place of EZH2  her2_b012v23_v2 <> ezh2_blt10vge10_v1.pp

table(wmadf$her2_b012v23_v2)
is.na(wmadf$her2_b012v23_v2) <- wmadf$her2_b012v23_v2 == "uninterpretable"

ah1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(0,45]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ah2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(45,60]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ah3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(60,70]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ah4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(70,99]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))

nh1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(0,45]" & er_b0v123_v2 == "negative" ))
nh2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(45,60]" & er_b0v123_v2 == "negative" ))
nh3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(60,70]" & er_b0v123_v2 == "negative" ))
nh4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = wmadf, 
                    subset = (ageg == "(70,99]" & er_b0v123_v2 == "negative" ))
par(mfcol=c(4,2))
plot(ah1.surv, lty = 2:3)
plot(ah2.surv, lty = 2:3)
plot(ah3.surv, lty = 2:3)
plot(ah4.surv, lty = 2:3)
plot(nh1.surv, lty = 2:3)
plot(nh2.surv, lty = 2:3)
plot(nh3.surv, lty = 2:3)
plot(nh4.surv, lty = 2:3)

par(mfrow=c(2,4))
plot(ah1.surv, lty = 2:3)
plot(ah2.surv, lty = 2:3)
plot(ah3.surv, lty = 2:3)
plot(ah4.surv, lty = 2:3)
plot(nh1.surv, lty = 2:3)
plot(nh2.surv, lty = 2:3)
plot(nh3.surv, lty = 2:3)
plot(nh4.surv, lty = 2:3)
legend("topright",legend = levels(as.factor(wmadf$her2_b012v23_v2))  , lty = 2:3)
title(main = "HER2 low vs high")

cchs <- with(wmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2))
cmhff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_b0v123_v2, data = wmadf[cchs, ])
cmhff
cmhfez <-  coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg , data = wmadf[cchs, ])
cmhfez
anova(cmhfez, cmhff)
cmhfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2 , data = wmadf[cchs, ])
cmhfer
anova(cmhfer, cmhff)

is.na(wmadf$grade_b12v3) <- wmadf$grade_b12v3 == "unknown"
is.na(wmadf$nodestat) <- wmadf$nodestat == "nodal status unknown"
is.na(wmadf$size_tumor_grp) <- wmadf$size_tumor_grp == "unknown"

table(wmadf$grade_b12v3)
table(wmadf$nodestat)
table(wmadf$size_tumor_grp)


cchsp <- with(wmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
cmhpff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                data = wmadf[cchsp, ])
cmhpff
cmhpfez <-  coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                  data = wmadf[cchsp, ])
cmhpfez
anova(cmhpfez, cmhpff)
cmhpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                  data = wmadf[cchsp, ])
cmhpfer
anova(cmhpfer, cmhpff)


### HER2 in place of EZH2  her2_b012v23_v2 <> ezh2_blt10vge10_v1.pp
### Training set
table(wmadf$isValidationSet_4543)

tmadf <- wmadf[wmadf$isValidationSet_4543 == "training set" ,]
dim(tmadf)

table(tmadf$her2_b012v23_v2)
is.na(tmadf$her2_b012v23_v2) <- tmadf$her2_b012v23_v2 == "uninterpretable"

ath1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ath2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ath3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ath4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))

nth1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "negative" ))
nth2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "negative" ))
nth3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "negative" ))
nth4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "negative" ))

par(mfrow=c(2,4))
plot(ath1.surv, lty = 2:3)
plot(ath2.surv, lty = 2:3)
plot(ath3.surv, lty = 2:3)
plot(ath4.surv, lty = 2:3)
plot(nth1.surv, lty = 2:3)
plot(nth2.surv, lty = 2:3)
plot(nth3.surv, lty = 2:3)
plot(nth4.surv, lty = 2:3)
legend("topright",legend = levels(as.factor(tmadf$her2_b012v23_v2))  , lty = 2:3)
title(main = "HER2 low vs high")

ccths <- with(tmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2))
cmthff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_b0v123_v2, data = tmadf[ccths, ])
cmthff
cmthfez <-  coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg , data = tmadf[ccths, ])
cmthfez
anova(cmthfez, cmthff)
cmthfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2 , data = tmadf[ccths, ])
cmthfer
anova(cmthfer, cmthff)

is.na(tmadf$grade_b12v3) <- tmadf$grade_b12v3 == "unknown"
is.na(tmadf$nodestat) <- tmadf$nodestat == "nodal status unknown"
is.na(tmadf$size_tumor_grp) <- tmadf$size_tumor_grp == "unknown"

table(tmadf$grade_b12v3)
table(tmadf$nodestat)
table(tmadf$size_tumor_grp)


ccthsp <- with(tmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
cmthpff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = tmadf[ccthsp, ])
cmthpff
cmthpfez <-  coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[ccthsp, ])
cmthpfez
anova(cmthpfez, cmthpff)
cmthpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[ccthsp, ])
cmthpfer
anova(cmthpfer, cmthpff)


### EZH2 ezh2_blt10vge10_v1.pp   her2_b012v23_v2 <> 
### Training set
table(wmadf$isValidationSet_4543)

tmadf <- wmadf[wmadf$isValidationSet_4543 == "training set" ,]
dim(tmadf)


atz1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
atz2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
atz3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
atz4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))

ntz1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "negative" ))
ntz2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "negative" ))
ntz3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "negative" ))
ntz4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "negative" ))

par(mfrow=c(2,4))
plot(atz1.surv, lty = 2:3)
plot(atz2.surv, lty = 2:3)
plot(atz3.surv, lty = 2:3)
plot(atz4.surv, lty = 2:3)
plot(ntz1.surv, lty = 2:3)
plot(ntz2.surv, lty = 2:3)
plot(ntz3.surv, lty = 2:3)
plot(ntz4.surv, lty = 2:3)
legend("topright",legend = levels(as.factor(tmadf$ezh2_blt10vge10_v1.pp))  , lty = 2:3)
title(main = "EZH2 low vs high")

cctzs <- with(tmadf, complete.cases(ageg, ezh2_blt10vge10_v1.pp, er_b0v123_v2))
cmtzff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg * er_b0v123_v2, data = tmadf[cctzs, ])
cmtzff
cmtzfez <-  coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg , data = tmadf[cctzs, ])
cmtzfez
anova(cmtzfez, cmtzff)
cmtzfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2 , data = tmadf[cctzs, ])
cmtzfer
anova(cmtzfer, cmtzff)

is.na(tmadf$grade_b12v3) <- tmadf$grade_b12v3 == "unknown"
is.na(tmadf$nodestat) <- tmadf$nodestat == "nodal status unknown"
is.na(tmadf$size_tumor_grp) <- tmadf$size_tumor_grp == "unknown"

table(tmadf$grade_b12v3)
table(tmadf$nodestat)
table(tmadf$size_tumor_grp)


cctzsp <- with(tmadf, complete.cases(ageg, ezh2_blt10vge10_v1.pp, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
cmtzpff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = tmadf[cctzsp, ])
cmtzpff
cmtzpfez <-  coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[cctzsp, ])
cmtzpfez
anova(cmtzpfez, cmtzpff)
cmtzpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[cctzsp, ])
cmtzpfer
anova(cmtzpfer, cmtzpff)


### HER2 in place of EZH2  her2_b012v23_v2 <> ezh2_blt10vge10_v1.pp
### Validation set
table(wmadf$isValidationSet_4543)

tmadf <- wmadf[wmadf$isValidationSet_4543 == "training set" ,]
dim(tmadf)

table(tmadf$her2_b012v23_v2)
is.na(tmadf$her2_b012v23_v2) <- tmadf$her2_b012v23_v2 == "uninterpretable"

ath1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ath2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ath3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
ath4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))

nth1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "negative" ))
nth2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "negative" ))
nth3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "negative" ))
nth4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "negative" ))

par(mfrow=c(2,4))
plot(ath1.surv, lty = 2:3)
plot(ath2.surv, lty = 2:3)
plot(ath3.surv, lty = 2:3)
plot(ath4.surv, lty = 2:3)
plot(nth1.surv, lty = 2:3)
plot(nth2.surv, lty = 2:3)
plot(nth3.surv, lty = 2:3)
plot(nth4.surv, lty = 2:3)
legend("topright",legend = levels(as.factor(tmadf$her2_b012v23_v2))  , lty = 2:3)
title(main = "HER2 low vs high")

ccths <- with(tmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2))
cmthff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_b0v123_v2, data = tmadf[ccths, ])
cmthff
cmthfez <-  coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg , data = tmadf[ccths, ])
cmthfez
anova(cmthfez, cmthff)
cmthfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2 , data = tmadf[ccths, ])
cmthfer
anova(cmthfer, cmthff)

is.na(tmadf$grade_b12v3) <- tmadf$grade_b12v3 == "unknown"
is.na(tmadf$nodestat) <- tmadf$nodestat == "nodal status unknown"
is.na(tmadf$size_tumor_grp) <- tmadf$size_tumor_grp == "unknown"

table(tmadf$grade_b12v3)
table(tmadf$nodestat)
table(tmadf$size_tumor_grp)


ccthsp <- with(tmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
cmthpff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = tmadf[ccthsp, ])
cmthpff
cmthpfez <-  coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[ccthsp, ])
cmthpfez
anova(cmthpfez, cmthpff)
cmthpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[ccthsp, ])
cmthpfer
anova(cmthpfer, cmthpff)


### EZH2 ezh2_blt10vge10_v1.pp   her2_b012v23_v2 <> 
### Training set
table(wmadf$isValidationSet_4543)

tmadf <- wmadf[wmadf$isValidationSet_4543 == "training set" ,]
dim(tmadf)
vmadf <- wmadf[wmadf$isValidationSet_4543 == "validation set" ,]
dim(vmadf)

## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"
atz1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
atz2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
atz3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))
atz4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "any nuclei staining (score 1-3)" ))

ntz1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(0,45]" & er_b0v123_v2 == "negative" ))
ntz2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(45,60]" & er_b0v123_v2 == "negative" ))
ntz3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(60,70]" & er_b0v123_v2 == "negative" ))
ntz4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = tmadf, 
                     subset = (ageg == "(70,99]" & er_b0v123_v2 == "negative" ))

par(mfrow=c(2,4))
plot(atz1.surv, lty = 2:3)
plot(atz2.surv, lty = 2:3)
plot(atz3.surv, lty = 2:3)
plot(atz4.surv, lty = 2:3)
plot(ntz1.surv, lty = 2:3)
plot(ntz2.surv, lty = 2:3)
plot(ntz3.surv, lty = 2:3)
plot(ntz4.surv, lty = 2:3)
legend("topright",legend = levels(as.factor(tmadf$ezh2_blt10vge10_v1.pp))  , lty = 2:3)
title(main = "EZH2 low vs high")

cctzs <- with(tmadf, complete.cases(ageg, ezh2_blt10vge10_v1.pp, er_b0v123_v2))
cmtzff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg * er_b0v123_v2, data = tmadf[cctzs, ])
cmtzff
cmtzfez <-  coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg , data = tmadf[cctzs, ])
cmtzfez
anova(cmtzfez, cmtzff)
cmtzfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2 , data = tmadf[cctzs, ])
cmtzfer
anova(cmtzfer, cmtzff)

is.na(tmadf$grade_b12v3) <- tmadf$grade_b12v3 == "unknown"
is.na(tmadf$nodestat) <- tmadf$nodestat == "nodal status unknown"
is.na(tmadf$size_tumor_grp) <- tmadf$size_tumor_grp == "unknown"

table(tmadf$grade_b12v3)
table(tmadf$nodestat)
table(tmadf$size_tumor_grp)


cctzsp <- with(tmadf, complete.cases(ageg, ezh2_blt10vge10_v1.pp, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
cmtzpff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = tmadf[cctzsp, ])
cmtzpff
cmtzpfez <-  coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[cctzsp, ])
cmtzpfez
anova(cmtzpfez, cmtzpff)
cmtzpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[cctzsp, ])
cmtzpfer
anova(cmtzpfer, cmtzpff)

### Need more levels than neg and pos to show idea for ER  ezh2_blt10vge10_v1.pp  er_v2
table(wmadf$er_v2)

is.na(wmadf$er_v2) <- wmadf$er_v2 == "uninterpretable/missing"


erpw1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                      subset = (ageg == "(0,45]"))
erpw2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                      subset = (ageg == "(45,60]"))
erpw3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                      subset = (ageg == "(60,70]"))
erpw4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                      subset = (ageg == "(70,99]"))

par(mfrow=c(2,2))
plot(erpw1.surv, lty = 2:5, main = "(0,45]")
plot(erpw2.surv, lty = 2:5, main = "(45,60]")
plot(erpw3.surv, lty = 2:5, main = "(60,70]")
plot(erpw4.surv, lty = 2:5, main = "(70,99]")

legend("topright",legend = levels(as.factor(wmadf$er_v2))  , lty = 2:5)
title(main = "ER low to high", line=3)
### Separation of pos curves for 0-45??

cce4tsp <- with(tmadf, complete.cases(ageg, er_b0v123_v2, er_v2, grade_b12v3, nodestat, size_tumor_grp))

cme4tpff <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = tmadf[cce4tsp, ])
cme4tpff
cme4tpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = tmadf[cce4tsp, ])
cme4tpfer
anova(cme4tpfer, cme4tpff)

cce4vsp <- with(vmadf, complete.cases(ageg, er_b0v123_v2, er_v2, grade_b12v3, nodestat, size_tumor_grp))

cme4vpff <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = vmadf[cce4vsp, ])
cme4vpff
cme4vpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = vmadf[cce4vsp, ])
cme4vpfer
anova(cme4vpfer, cme4vpff)

cce4wsp <- with(wmadf, complete.cases(ageg, er_b0v123_v2, er_v2, grade_b12v3, nodestat, size_tumor_grp))

cme4wpff <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                   data = wmadf[cce4wsp, ])
cme4wpff
cme4wpfer <-  coxph(Surv(survyrs1999, survstat1999) ~ ageg * er_b0v123_v2  + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cce4wsp, ])
cme4wpfer
anova(cme4wpfer, cme4wpff)

### Split by EZH2

table(wmadf$ezh2_blt10vge10_v1.pp)


erzlpw1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(0,45]" & ezh2_blt10vge10_v1.pp == "0 to 5%"))
erzlpw2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(45,60]" & ezh2_blt10vge10_v1.pp == "0 to 5%"))
erzlpw3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(60,70]" & ezh2_blt10vge10_v1.pp == "0 to 5%"))
erzlpw4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(70,99]" & ezh2_blt10vge10_v1.pp == "0 to 5%"))
erzhpw1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(0,45]" & ezh2_blt10vge10_v1.pp == "10 to 100%"))
erzhpw2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(45,60]" & ezh2_blt10vge10_v1.pp == "10 to 100%"))
erzhpw3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(60,70]" & ezh2_blt10vge10_v1.pp == "10 to 100%"))
erzhpw4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v2 , data = wmadf, 
                        subset = (ageg == "(70,99]" & ezh2_blt10vge10_v1.pp == "10 to 100%"))

par(mfrow=c(2,4))
plot(erzlpw1.surv, lty = 1:4, main = "(0,45]")
plot(erzlpw2.surv, lty = 1:4, main = "(45,60]")
plot(erzlpw3.surv, lty = 1:4, main = "(60,70]")
plot(erzlpw4.surv, lty = 1:4, main = "(70,99]")
plot(erzhpw1.surv, lty = 1:4, main = "(0,45]")
plot(erzhpw2.surv, lty = 1:4, main = "(45,60]")
plot(erzhpw3.surv, lty = 1:4, main = "(60,70]")
plot(erzhpw4.surv, lty = 1:4, main = "(70,99]")

legend("topright",legend = levels(as.factor(wmadf$er_v2))  , lty = 1:4)
title(main = "ER low to high", line=3)

## er vs age

wmadf$ageg4 <- cut(wmadf$age, breaks = c(0, 35, 55, 75, 99))
table(wmadf$ageg4)
# 
# (0,35] (35,55] (55,75] (75,99] 
# 138    1449    1979     425 
wmadf$ageg4 <- cut(wmadf$age, breaks = c(0, 40, 55, 70, 99))
table(wmadf$ageg4)
# 
# (0,40] (40,55] (55,70] (70,99] 
# 376    1211    1527     877 

table(wmadf$er_v3n, useNA = "always")

wmadf$er_v3n_6c <- cut(wmadf$er_v3n, breaks = c(-0.5, 0.5, 10.5, 25.5, 50.5, 75.5, 100.5), 
                       labels = c("0%", "1-10%", "11-25%", "26-50%", "51-75%", "76-100%"))
table(wmadf$er_v3n_6c, useNA = "always")
with(wmadf, CrossTable(ageg4, er_v3n_6c, expected = TRUE, chisq = TRUE))

wmadf$ageg6 <- cut(wmadf$age, breaks = c(0, 35, 45, 55, 65, 75, 99))
table(wmadf$ageg6)
# 
# (0,35] (35,45] (45,55] (55,65] (65,75] (75,99] 
# 138     651     798     975    1004     425 
# 

with(wmadf, CrossTable(ageg, er_v2, expected = TRUE, chisq = TRUE))

with(wmadf[wmadf$isValidationSet_4543 == "training set" , ], CrossTable(ageg6, er_v2, expected = TRUE, chisq = TRUE))
with(wmadf[wmadf$isValidationSet_4543 == "validation set" , ], CrossTable(ageg6, er_v2, expected = TRUE, chisq = TRUE))
with(wmadf, CrossTable(ageg6, er_v2, expected = TRUE, chisq = TRUE))

wmadf$er_c4_v2 <- 
  factor(wmadf$er_v2, 
         levels = c("no nuclei stained or <1%",
                    "1-25% of the nuclei stained",
                    "25-75% of the nuclei stained",
                    ">75% of the nuclei stained"))
levels(wmadf$er_c4_v2) <-c("none or <1%",
                           "1-25%",
                           "25-75%",
                           ">75% nuclei stained")
table(wmadf$er_c4_v2,wmadf$er_v2)

with(wmadf[wmadf$isValidationSet_4543 == "training set" , ], CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))
with(wmadf[wmadf$isValidationSet_4543 == "validation set" , ], CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))
with(wmadf, CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))
## ER pos subset only
with(wmadf[wmadf$er_c4_v2 != "none or <1%", ], CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))


##wmadf$her2_b012v23_v2
with(wmadf[wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}", ], CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))

with(wmadf[wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}", ], CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))

table(wmadf$her2_b012v23_v2=="{2 w/ FISH pos, 3}")
table(wmadf$her2_b012v23_v2)

### Split by xxxx

table(wmadf$ezh2_blt10vge10_v1.pp)


erhloag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(0,35]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}"  ))
erhloag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(35,45]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(45,55]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(55,65]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag5.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(65,75]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag6.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(75,99]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhhiag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(0,35]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}"  ))
erhhiag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(35,45]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(45,55]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(55,65]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag5.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(65,75]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag6.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c4_v2 , data = wmadf, 
                         subset = (ageg6 == "(75,99]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))

par(mfrow=c(3,4))
plot(erhloag1.surv, lty = 1:4, main = "(0,35] HER2-")
plot(erhloag2.surv, lty = 1:4, main = "(35,45] HER2-")
plot(erhloag3.surv, lty = 1:4, main = "(45,55] HER2-")
plot(erhloag4.surv, lty = 1:4, main = "(55,65] HER2-")
plot(erhloag5.surv, lty = 1:4, main = "(65,75] HER2-")
plot(erhloag6.surv, lty = 1:4, main = "(75,99] HER2-")
plot(erhhiag1.surv, lty = c(1, 3), main = "(0,35] HER2+")
plot(erhhiag2.surv, lty = 1:4, main = "(35,45] HER2+")
plot(erhhiag3.surv, lty = 1:4, main = "(45,55] HER2+")
plot(erhhiag4.surv, lty = 1:4, main = "(55,65] HER2+")
plot(erhhiag5.surv, lty = 1:4, main = "(65,75] HER2+")
plot(erhhiag6.surv, lty = c(1, 3, 4), main = "(75,99] HER2+")

legend("topright",legend = levels(as.factor(wmadf$er_c4_v2))  , lty = 1:4)
title(main = "ER low to high", line=3)

### ER in 3 cats

wmadf$er_c3_v2c <- rep(NA_character_, nrow(wmadf))
wmadf$er_c3_v2c[wmadf$er_v2 == "no nuclei stained or <1%"] <- "no nuclei stained or <1%"
wmadf$er_c3_v2c[wmadf$er_v2 == "1-25% of the nuclei stained"] <- "1-75% of the nuclei stained"
wmadf$er_c3_v2c[wmadf$er_v2 == "25-75% of the nuclei stained"] <- "1-75% of the nuclei stained"
wmadf$er_c3_v2c[wmadf$er_v2 == ">75% of the nuclei stained"] <- ">75% of the nuclei stained"

wmadf$er_c3_v2 <- 
  factor(wmadf$er_c3_v2c, 
         levels = c("no nuclei stained or <1%", 
                    "1-75% of the nuclei stained", 
                    ">75% of the nuclei stained"))
table(wmadf$er_v2, wmadf$er_c3_v2)
levels(wmadf$er_c3_v2) <-c("none or <1%",
                           "1-75%",
                           ">75% nuclei stained")

table(wmadf$er_v2,wmadf$er_c3_v2)

with(wmadf[wmadf$isValidationSet_4543 == "training set" , ], CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))
with(wmadf[wmadf$isValidationSet_4543 == "validation set" , ], CrossTable(ageg6, er_c4_v2, expected = TRUE, chisq = TRUE))
with(wmadf, CrossTable(ageg6, er_c3_v2, expected = TRUE, chisq = TRUE))


##wmadf$her2_b012v23_v2
with(wmadf[wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}", ], CrossTable(ageg6, er_c3_v2, expected = TRUE, chisq = TRUE))

with(wmadf[wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}", ], CrossTable(ageg6, er_c3_v2, expected = TRUE, chisq = TRUE))

table(wmadf$her2_b012v23_v2=="{2 w/ FISH pos, 3}")
table(wmadf$her2_b012v23_v2)

### Split by xxxx

table(wmadf$ezh2_blt10vge10_v1.pp)


erhloag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(0,35]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}"  ))
erhloag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(35,45]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(45,55]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(55,65]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag5.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(65,75]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag6.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(75,99]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhhiag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(0,35]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}"  ))
erhhiag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(35,45]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(45,55]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(55,65]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag5.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(65,75]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag6.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg6 == "(75,99]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))

pdf(file="./02-008_BigSeries/Analysis/Survival/ERbyIHCvsAgevsHER2.pdf", useDingbats=FALSE, width=8, height = 10)
par(mfrow=c(4,3))

plot(erhloag1.surv, lty = 1:3, main = "(0,35] HER2-")

plot(erhloag2.surv, lty = 1:3, main = "(35,45] HER2-")
title(main = "ER low to high", line=3)
legend("bottomleft",legend = levels(as.factor(wmadf$er_c3_v2))  , lty = 1:3)

plot(erhloag3.surv, lty = 1:3, main = "(45,55] HER2-")
plot(erhloag4.surv, lty = 1:3, main = "(55,65] HER2-")
plot(erhloag5.surv, lty = 1:3, main = "(65,75] HER2-")
plot(erhloag6.surv, lty = 1:3, main = "(75,99] HER2-")
plot(erhhiag1.surv, lty = c(1, 2), main = "(0,35] HER2+")
plot(erhhiag2.surv, lty = 1:3, main = "(35,45] HER2+")
plot(erhhiag3.surv, lty = 1:3, main = "(45,55] HER2+")
plot(erhhiag4.surv, lty = 1:3, main = "(55,65] HER2+")
plot(erhhiag5.surv, lty = 1:3, main = "(65,75] HER2+")
plot(erhhiag6.surv, lty = 1:3, main = "(75,99] HER2+")

dev.off()



ccwhsp <- with(wmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
cmwhpff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg6 * er_c3_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccwhsp, ])
cmwhpff
cmwhpfr <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg6 * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccwhsp, ])
cmwhpfr
anova(cmwhpfr, cmwhpff)

cmwhpfrr <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg6 + grade_b12v3 + nodestat + size_tumor_grp, 
                  data = wmadf[ccwhsp, ])
cmwhpfrr
anova(cmwhpfrr, cmwhpff)
anova(cmwhpfrr, cmwhpfr)


###
table(wmadf$ageg)

erhloag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(0,45]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}"  ))
erhloag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(45,60]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(60,70]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhloag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(70,99]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erhhiag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(0,45]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}"  ))
erhhiag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(45,60]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(60,70]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erhhiag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_c3_v2 , data = wmadf, 
                         subset = (ageg == "(70,99]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))

pdf(file="./02-008_BigSeries/Analysis/Survival/ERbyIHCvsAge4vsHER2.pdf", useDingbats=FALSE, width=8, height = 10)
par(mfcol=c(4,2))

plot(erhloag1.surv, lty = 1:3, main = "(0,45] HER2-")
title(main = "ER low to high", line=3)
legend("bottomleft",legend = levels(as.factor(wmadf$er_c3_v2))  , lty = 1:3)

plot(erhloag2.surv, lty = 1:3, main = "(45,60] HER2-")
plot(erhloag3.surv, lty = 1:3, main = "(60,70] HER2-")
plot(erhloag4.surv, lty = 1:3, main = "(70,99] HER2-")
plot(erhhiag1.surv, lty = 1:3, main = "(0,45] HER2+")
plot(erhhiag2.surv, lty = 1:3, main = "(45,60] HER2+")
plot(erhhiag3.surv, lty = 1:3, main = "(60,70] HER2+")
plot(erhhiag4.surv, lty = 1:3, main = "(70,99] HER2+")

dev.off()



ccwhsp <- with(wmadf, complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp))
head(ccwhsp)
cmwhpff <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_c3_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccwhsp, ])
cmwhpff
cmwhpfr <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg * er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccwhsp, ])
cmwhpfr
anova(cmwhpfr, cmwhpff)

cmwhpfrr <- coxph(Surv(survyrs1999, survstat1999) ~ her2_b012v23_v2 + ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                  data = wmadf[ccwhsp, ])
cmwhpfrr
anova(cmwhpfrr, cmwhpff)
anova(cmwhpfrr, cmwhpfr)


##wmadf$her2_b012v23_v2
with(wmadf[wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}", ], CrossTable(ageg, er_c3_v2, expected = TRUE, chisq = TRUE))

with(wmadf[wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}", ], CrossTable(ageg, er_c3_v2, expected = TRUE, chisq = TRUE))

### 0-45 HER2-


ccwhsp <- with(wmadf, (complete.cases(ageg, her2_b012v23_v2, er_b0v123_v2, grade_b12v3, nodestat, size_tumor_grp) & 
                         (ageg == "(0,45]") &
                         (her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}") ) )
table(ccwhsp)
head(ccwhsp)
cmwhpff <- coxph(Surv(survyrs1999, survstat1999) ~ er_c3_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccwhsp, ])
cmwhpff
cmwhpfr <- coxph(Surv(survyrs1999, survstat1999) ~ er_b0v123_v2 + grade_b12v3 + nodestat + size_tumor_grp, 
                 data = wmadf[ccwhsp, ])
cmwhpfr
anova(cmwhpfr, cmwhpff)
cmwhpff <- coxph(Surv(survyrs1999, survstat1999) ~ er_c3_v2, 
                 data = wmadf[ccwhsp, ])
cmwhpff
cmwhpfr <- coxph(Surv(survyrs1999, survstat1999) ~ er_b0v123_v2, 
                 data = wmadf[ccwhsp, ])
cmwhpfr
anova(cmwhpfr, cmwhpff)

## Read in ER percent cells stained data from Samuel  ERpp
##  /Users/stevenmckinney/git/Dapple/TMA/02-008_BigSeries/Scoring/ESR1_ER/er_v3.1_labels.csv 

big_er31 <- 
  read.table(file = paste("/Users/stevenmckinney/git/Dapple/TMA/",
                          "02-008_BigSeries/Scoring/ESR1_ER/",
                          "er_v3.1_labels.csv", sep = ""), sep = ",", header = FALSE, 
             col.names = c("bseries_id", "er_v3.1c"), stringsAsFactors = FALSE  )
wmadf <- merge(wmadf, big_er31, all.x=TRUE, all.y=FALSE )
wmadf$er_v3n <- as.numeric(wmadf$er_v3.1c)

## ERpp density plot:
# 
# horder <- c(8,7,6,5,4,3,9,2,1)
# #horder <- c(1:9)
# unique(tmaagedf$datasource)[horder]
# #> unique(tmaagedf$datasource)[horder]
# #[1] "TCGA ProstateCa"    "TCGA LungCa"        "TCGA KidneyCa"      "TCGA BrCa"          "METABRIC"           "NormalBreast13 TMA" "TCGA ThyroidCa"    
# #[8] "MB09 TMA"
# 
# tmaagedf$datasourcef <- factor(tmaagedf$datasource, levels = unique(tmaagedf$datasource)[horder])
# 
# ggplot(data=tmaagedf,aes(x=age_at_diagnosis, group=datasourcef, fill=datasourcef)) + 
#   geom_density(adjust=1.5 , alpha=0.7) +
#   xlab("Age at diagnosis (years)") + ylab("Density") +
#   scale_fill_manual(name = "Data source", 
#                     values = c("#FF5500", "#EE82EE",  "#00EE76", "#CD3278", "#8B0000",
#                                "#0000CD", "#FFAA00", "#00C5CD", "#FFFF40")[horder])
# 
# tmaagedf$datasourcef <- factor(tmaagedf$datasource, levels = sort(unique(tmaagedf$datasource)))
# ggplot(data=tmaagedf, aes(x = age_at_diagnosis, y = datasourcef)) + 
#   geom_density_ridges() + 
#   theme_minimal(base_size = 14) + theme(axis.text.y = element_text(vjust = 0)) +
#   xlab("Age at diagnosis (years)") + ylab("Density") +
#   scale_x_continuous(expand = c(0.01, 0)) +
#   scale_y_discrete(expand = c(0.01, 0))
# 
library(ggplot2)
library(ggridges)

pdf(file=paste(plot_dirs, "Age_ER_DensityPlots.pdf", sep = "/"), useDingbats = FALSE, width = 5, height = 5, paper = "letter")

png(file=paste(plot_dirs, "Age_ER_DensityPlots.png", sep = "/"),  width = 700, height = 500, res = 150)
ggplot(data = wmadf[which(!is.na(wmadf$er_v3n)),], aes(x = er_v3n, y = ageg4)) + geom_density_ridges(rel_min_height=0.01)  +
  xlab("ER: Percent of tumour cells positively stained") + ylab("Density by age at diagnosis group") +
  scale_x_continuous(expand = c(0.01, 0))  + annotate("text", label = paste("N =", table(wmadf[which(!is.na(wmadf$er_v3n)), "ageg4"])), x = 55, y = 0.35+1:4) +
  scale_y_discrete(expand = c(0.01, 0),  
                   labels = c("(0,40]" = "40 years or younger","(40,55]" = "41-55","(55,70]" = "56-70","(70,99]" = "71 years or older"))
dev.off()
## Positive cases only
pdf(file=paste(plot_dirs, "Age_ERpositive_DensityPlots.pdf", sep = "/"), useDingbats = FALSE, width = 5, height = 5, paper = "letter")
png(file=paste(plot_dirs, "Age_ERpositive_DensityPlots.png", sep = "/"),  width = 700, height = 500, res = 150)
ggplot(data = wmadf[which(wmadf$er_v3n >= 1), ], aes(x = er_v3n, y = ageg4)) + geom_density_ridges(rel_min_height=0.01)  +
  xlab("ER positive cases:\n Percent of tumour cells positively stained") + ylab("Density by age at diagnosis group") +
  scale_x_continuous(expand = c(0.01, 0)) + annotate("text", label = paste("N =", table(wmadf[which(wmadf$er_v3n >= 1), "ageg4"])), x = 55, y = 0.65+1:4) +
  scale_y_discrete(expand = c(0.01, 0),  
                   labels = c("(0,40]" = "40 years or younger","(40,55]" = "41-55","(55,70]" = "56-70","(70,99]" = "71 years or older"))
dev.off()
table(wmadf[which(wmadf$er_v3n >= 1), "ageg4"])

ggplot(data = wmadf, aes(x = er_v3n, y = ageg4)) + stat_binline(breaks=c(0,1:20,c(3:8)*10,81:100)) 
ggplot(data = wmadf, aes(x = er_v3n, y = ageg4)) + stat_binline(breaks=c(0,1,5,c(1:9)*10,95,100)) 
ggplot(data = wmadf, aes(x = er_v3n, y = ageg4)) + stat_binline(breaks=c(0,5,c(1:9)*10,95,100)) 

wmadf$er_v3c4 <- cut(wmadf$er_v3n, breaks = c(0, 1, 11, 51, 100), right = FALSE, include=TRUE, labels = c("[0,1)","[1,10]","[11,50]","[51,100]")   )


erv3hloag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(0,45]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}"  ))
erv3hloag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(45,60]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erv3hloag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(60,70]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erv3hloag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(70,99]" & wmadf$her2_b012v23_v2=="Her2 = {0,1,2 w/ FISH neg}" ))
erv3hhiag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(0,45]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}"  ))
erv3hhiag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(45,60]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erv3hhiag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(60,70]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))
erv3hhiag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ er_v3c4 , data = wmadf, 
                           subset = (ageg == "(70,99]" & wmadf$her2_b012v23_v2=="Her2 = {2 w/ FISH pos, 3}" ))

pdf(file="./02-008_BigSeries/Analysis/Survival/ERv3byIHCvsAge4vsHER2.pdf", useDingbats=FALSE, width=8, height = 10)
par(mfcol=c(4,2))

plot(erv3hloag1.surv, lty = 1:4, main = "(0,45] HER2-")
title(main = "ER low to high", line=3)
legend("bottomleft",legend = levels(as.factor(wmadf$er_v3c4))  , lty = 1:4)

plot(erv3hloag2.surv, lty = 1:4, main = "(45,60] HER2-")
plot(erv3hloag3.surv, lty = 1:4, main = "(60,70] HER2-")
plot(erv3hloag4.surv, lty = 1:4, main = "(70,99] HER2-")
plot(erv3hhiag1.surv, lty = 1:4, main = "(0,45] HER2+")
plot(erv3hhiag2.surv, lty = 1:4, main = "(45,60] HER2+")
plot(erv3hhiag3.surv, lty = 1:4, main = "(60,70] HER2+")
plot(erv3hhiag4.surv, lty = 1:4, main = "(70,99] HER2+")

dev.off()

### EZH2 and H3K27me3
# 
# > names(wmadf)[which(grepl("h3k", names(wmadf), ignore.case = TRUE))]
# [1] "h3k27me3_v1"        "h3k27me3_v1_b01v23" "h3k27me3_v1_b012v3" "h3k27me3_v1_cat5"   "h3k27me3_v1_all"    "h3k27me3_v1.1"      "h3k27me3_v1.1.pp"   "h3k27me3_v1.1.int" 

# [241] "col"                                    "score"                                  "h3k27me3_v1"                            "h3k27me3_v1_b01v23"                    
# [245] "h3k27me3_v1_b012v3"                     "h3k27me3_v1_cat5"                       "h3k27me3_v1_all"                        "h3k27me3_v1.1"                         
# [249] "h3k27me3_v1.1.pp"                       "h3k27me3_v1.1.int"                      "her2_v1"                                "her2_bXvS_v1"                          
# [253] "her2_b01v3_v1"                          "her2_v2"                                "her2_bXvS_v2"                           "her2_b01v3_v2"                         
# 

wmadf$h3k27me3_v1.1n <- as.numeric(wmadf$h3k27me3_v1.1.pp)
quantile(wmadf$h3k27me3_v1.1n, na.rm=TRUE)


## Original H3K27me3 analysis done on the odd classification by image colour in h3k27me3_v1
## Get some info on the percent positive version.

splitpctvec <- c(60, 70, 80, 90, 95)
splitpct <- 90

wmadf$h3k27me3_v1.1.pp_b40c <- rep(NA_character_, nrow(wmadf))
wmadf$h3k27me3_v1.1.pp_b50c <- rep(NA_character_, nrow(wmadf))
wmadf$h3k27me3_v1.1.pp_b60c <- rep(NA_character_, nrow(wmadf))
wmadf$h3k27me3_v1.1.pp_b70c <- rep(NA_character_, nrow(wmadf))
wmadf$h3k27me3_v1.1.pp_b80c <- rep(NA_character_, nrow(wmadf))
wmadf$h3k27me3_v1.1.pp_b90c <- rep(NA_character_, nrow(wmadf))
wmadf$h3k27me3_v1.1.pp_b95c <- rep(NA_character_, nrow(wmadf))

wmadf$h3k27me3_v1.1.pp_b40c[wmadf$h3k27me3_v1.1n <  40] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_b40c[wmadf$h3k27me3_v1.1n >= 40] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_b40f <- factor(wmadf$h3k27me3_v1.1.pp_b40c, levels = c("H3K27me3 low", "H3K27me3 high"))
wmadf$h3k27me3_v1.1.pp_b50c[wmadf$h3k27me3_v1.1n <  50] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_b50c[wmadf$h3k27me3_v1.1n >= 50] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_b50f <- factor(wmadf$h3k27me3_v1.1.pp_b50c, levels = c("H3K27me3 low", "H3K27me3 high"))
wmadf$h3k27me3_v1.1.pp_b60c[wmadf$h3k27me3_v1.1n <  60] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_b60c[wmadf$h3k27me3_v1.1n >= 60] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_b60f <- factor(wmadf$h3k27me3_v1.1.pp_b60c, levels = c("H3K27me3 low", "H3K27me3 high"))
wmadf$h3k27me3_v1.1.pp_b70c[wmadf$h3k27me3_v1.1n <  70] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_b70c[wmadf$h3k27me3_v1.1n >= 70] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_b70f <- factor(wmadf$h3k27me3_v1.1.pp_b70c, levels = c("H3K27me3 low", "H3K27me3 high"))
wmadf$h3k27me3_v1.1.pp_b80c[wmadf$h3k27me3_v1.1n <  80] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_b80c[wmadf$h3k27me3_v1.1n >= 80] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_b80f <- factor(wmadf$h3k27me3_v1.1.pp_b80c, levels = c("H3K27me3 low", "H3K27me3 high"))
wmadf$h3k27me3_v1.1.pp_b90c[wmadf$h3k27me3_v1.1n <  90] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_b90c[wmadf$h3k27me3_v1.1n >= 90] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_b90f <- factor(wmadf$h3k27me3_v1.1.pp_b90c, levels = c("H3K27me3 low", "H3K27me3 high"))
wmadf$h3k27me3_v1.1.pp_b95c[wmadf$h3k27me3_v1.1n <  95] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_b95c[wmadf$h3k27me3_v1.1n >= 95] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_b95f <- factor(wmadf$h3k27me3_v1.1.pp_b95c, levels = c("H3K27me3 low", "H3K27me3 high"))

table(wmadf$h3k27me3_v1.1.pp_b40f , useNA = "always")
table(wmadf$h3k27me3_v1.1.pp_b50f , useNA = "always")
table(wmadf$h3k27me3_v1.1.pp_b60f , useNA = "always")
table(wmadf$h3k27me3_v1.1.pp_b70f , useNA = "always")
table(wmadf$h3k27me3_v1.1.pp_b80f , useNA = "always")
table(wmadf$h3k27me3_v1.1.pp_b90f , useNA = "always")
table(wmadf$h3k27me3_v1.1.pp_b95f , useNA = "always")
### Split by EZH2

table(wmadf$h3k27me3_v1_cat5)
table(wmadf$h3k27me3_v1.1.pp)
table(wmadf$h3k27me3_v1.1.pp_b60f)
table(wmadf$h3k27me3_v1.1.pp_b70f)
table(wmadf$h3k27me3_v1.1.pp_b80f)
table(wmadf$h3k27me3_v1.1.pp_b90f)
table(wmadf$h3k27me3_v1.1.pp_b95f)
table(wmadf$ezh2_blt10vge10_v1.pp)


h3v11b40.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b40f , data = wmadf)
h3v11b50.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b50f , data = wmadf)
h3v11b60.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf)
h3v11b70.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b70f , data = wmadf)
h3v11b80.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b80f , data = wmadf)
h3v11b90.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b90f , data = wmadf)
h3v11b95.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b95f , data = wmadf)

par(mfrow=c(3, 3))
plot(h3v11b40.surv, lty = 1:2, main = "40")
plot(h3v11b50.surv, lty = 1:2, main = "50")
plot(h3v11b60.surv, lty = 1:2, main = "60")
plot(h3v11b70.surv, lty = 1:2, main = "70")
plot(h3v11b80.surv, lty = 1:2, main = "80")
plot(h3v11b90.surv, lty = 1:2, main = "90")
plot(h3v11b95.surv, lty = 1:2, main = "95")

legend("topright",legend = levels(as.factor(wmadf$er_v2))  , lty = 1:4)
title(main = "ER low to high", line=3)

## Cox model fits
table(wmadf$ageg)

cch3c3ag4 <- with(wmadf, complete.cases(h3k27me3_v1.1.pp_b60f, ageg,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
cch3c3ag4 <- with(wmadf, complete.cases(h3k27me3_v1.1.pp_b90f, ageg,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(cch3c3ag4)
head(cch3c3ag4)

cmh3b40c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b40f * ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ag4, ])
summary(cmh3b40c3agff)
cmh3b40c3agff$loglik  ##  [1] -8467.863 -8251.090

cmh3b50c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b50f * ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ag4, ])
summary(cmh3b50c3agff)
cmh3b50c3agff$loglik  ##  [1] -8467.863 -8247.177

cmh3b60c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f * ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ag4, ])
summary(cmh3b60c3agff)
cmh3b60c3agff$loglik  ##  [1] -8467.863 -8247.011


cmh3b70c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b70f * ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ag4, ])
summary(cmh3b70c3agff)
cmh3b70c3agff$loglik  ##  [1] -8467.863 -8249.402


cmh3b80c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b80f * ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ag4, ])
summary(cmh3b80c3agff)
cmh3b80c3agff$loglik  ##  [1] -8467.863 -8249.540


cmh3b90c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b90f * ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ag4, ])
summary(cmh3b90c3agff)
cmh3b90c3agff$loglik ## [1] -8467.863 -8249.234


cmh3b95c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b95f * ageg + grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ag4, ])
summary(cmh3b95c3agff)
cmh3b95c3agff$loglik  ## [1] -8467.863 -8252.726

## 60 pct cells positively stained is binary optimum. <60 vs >=60  h3k27me3_v1.1.pp_b60f 

##

wmadf$h3k27me3_v1.1.pp_t6095c <- rep(NA_character_, nrow(wmadf))

wmadf$h3k27me3_v1.1.pp_t6095c[wmadf$h3k27me3_v1.1n <  60] <- "H3K27me3 low"
wmadf$h3k27me3_v1.1.pp_t6095c[wmadf$h3k27me3_v1.1n >= 60 & wmadf$h3k27me3_v1.1n < 90] <- "H3K27me3 mid"
wmadf$h3k27me3_v1.1.pp_t6095c[wmadf$h3k27me3_v1.1n >= 95] <- "H3K27me3 high"
wmadf$h3k27me3_v1.1.pp_t6095f <- factor(wmadf$h3k27me3_v1.1.pp_t6095c, levels = c("H3K27me3 low", "H3K27me3 mid", "H3K27me3 high"))
table(wmadf$h3k27me3_v1.1.pp_t6095f)

h3v11t6095.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf)
plot(h3v11t6095.surv, lty = 1:3, main = "[0,60),[60-95),[95,100]")

### Split by xxxx

table(wmadf$h3k27me3_v1.1.pp_t6095f)
names(wmadf)[which(grepl("ez", names(wmadf), ignore.case = TRUE))]

h3c3ag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                        subset = (ageg6 == "(0,35]"  ))
h3c3ag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                        subset = (ageg6 == "(35,45]" ))
h3c3ag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                        subset = (ageg6 == "(45,55]" ))
h3c3ag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                        subset = (ageg6 == "(55,65]" ))
h3c3ag5.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                        subset = (ageg6 == "(65,75]" ))
h3c3ag6.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                        subset = (ageg6 == "(75,99]" ))

par(mfrow=c(3,2))
plot(h3c3ag1.surv, lty = 1:3, main = "(0,35] ")
plot(h3c3ag2.surv, lty = 1:3, main = "(35,45]")
plot(h3c3ag3.surv, lty = 1:3, main = "(45,55]")
plot(h3c3ag4.surv, lty = 1:3, main = "(55,65]")
plot(h3c3ag5.surv, lty = 1:3, main = "(65,75]")
plot(h3c3ag6.surv, lty = 1:3, main = "(75,99]")


legend("topright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_t6095f))  , lty = 1:3)
title(main = "H3K27me3 low to high", line=3)

## Cox models


cch3c3ag <- with(wmadf, complete.cases(h3k27me3_v1.1.pp_t6095f, ageg6,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(cch3c3ag)
head(cch3c3ag)

cmh3c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f * ageg6 + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cch3c3ag, ])
cmh3c3agff
cmh3c3agfr <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f + ageg6 + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cch3c3ag, ])
cmh3c3agfr
anova(cmh3c3agfr, cmh3c3agff)
# 
# > anova(cmh3c3agfr, cmh3c3agff)
# Analysis of Deviance Table
# Cox model: response is  Surv(survyrs1999, survstat1999)
# Model 1: ~ h3k27me3_v1.1.pp_t6095f + ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# Model 2: ~ h3k27me3_v1.1.pp_t6095f * ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# loglik Chisq Df P(>|Chi|)
# 1 -7707.6                   
# 2 -7704.1 6.888 10     0.736
#
# No apparent difference in survival rate patterns among the 3 h3k27me3 groups across ages.

cmh3c3agrr <- coxph(Surv(survyrs1999, survstat1999) ~  ageg6 + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cch3c3ag, ])
cmh3c3agrr
anova(cmh3c3agrr, cmh3c3agff)
anova(cmh3c3agrr, cmh3c3agfr)
# 
# > anova(cmh3c3agrr, cmh3c3agff)
# Analysis of Deviance Table
# Cox model: response is  Surv(survyrs1999, survstat1999)
# Model 1: ~ ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# Model 2: ~ h3k27me3_v1.1.pp_t6095f * ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# loglik Chisq Df P(>|Chi|)   
# 1 -7717.5                      
# 2 -7704.1 26.75 12  0.008393 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > anova(cmh3c3agrr, cmh3c3agfr)
# Analysis of Deviance Table
# Cox model: response is  Surv(survyrs1999, survstat1999)
# Model 1: ~ ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# Model 2: ~ h3k27me3_v1.1.pp_t6095f + ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# loglik  Chisq Df P(>|Chi|)    
# 1 -7717.5                        
# 2 -7707.6 19.862  2 4.864e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > Definite differences in survival rates due to H3K27me3 percent of cells positively stained.

### Add ezh2 to subset

table(wmadf$ezh2_blt10vge10_v1.pp)
# > table(wmadf$ezh2_blt10vge10_v1.pp)
# 
# 0 to 5% 10 to 100% 
#   1891       1747 

### Split by xxxx

table(wmadf$h3k27me3_v1.1.pp_t6095f)
names(wmadf)[which(grepl("ez", names(wmadf), ignore.case = TRUE))]

h3c3enag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(0,35]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c3enag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(35,45]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c3enag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(45,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c3enag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(55,65]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c3enag5.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(65,75]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c3enag6.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(75,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c3epag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(0,35]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c3epag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(35,45]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c3epag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(45,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c3epag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(55,65]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c3epag5.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(65,75]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c3epag6.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f , data = wmadf, 
                          subset = ((ageg6 == "(75,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )

par(mfrow=c(3,4))
plot(h3c3enag1.surv, lty = 1:3, main = "(0,35] EZH2 neg")
plot(h3c3enag2.surv, lty = 1:3, main = "(35,45] EZH2 neg")
plot(h3c3enag3.surv, lty = 1:3, main = "(45,55] EZH2 neg")
plot(h3c3enag4.surv, lty = 1:3, main = "(55,65] EZH2 neg")
plot(h3c3enag5.surv, lty = 1:3, main = "(65,75] EZH2 neg")
plot(h3c3enag6.surv, lty = 1:3, main = "(75,99] EZH2 neg")
plot(h3c3epag1.surv, lty = 1:3, main = "(0,35] EZH2 pos")
plot(h3c3epag2.surv, lty = 1:3, main = "(35,45] EZH2 pos")
plot(h3c3epag3.surv, lty = 1:3, main = "(45,55] EZH2 pos")
plot(h3c3epag4.surv, lty = 1:3, main = "(55,65] EZH2 pos")
plot(h3c3epag5.surv, lty = 1:3, main = "(65,75] EZH2 pos")
plot(h3c3epag6.surv, lty = 1:3, main = "(75,99] EZH2 pos")


legend("topright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_t6095f))  , lty = 1:3)
title(main = "H3K27me3 low to high", line=3)

## Cox


cch3c3ezag <- with(wmadf, complete.cases(h3k27me3_v1.1.pp_t6095f, ezh2_blt10vge10_v1.pp, 
                                         ageg6,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(cch3c3ezag)
head(cch3c3ezag)

cmh3c3ezagff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp * h3k27me3_v1.1.pp_t6095f * ageg6 + 
                        grade_b12v3 + nodestat + size_tumor_grp, 
                      data = wmadf[cch3c3ezag, ])
cmh3c3ezagff
cmh3c3ezagfr <- coxph(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_t6095f * ageg6 + grade_b12v3 + nodestat + size_tumor_grp, 
                      data = wmadf[cch3c3ezag, ])
cmh3c3ezagfr
anova(cmh3c3ezagfr, cmh3c3ezagff)
# 
# > anova(cmh3c3ezagfr, cmh3c3ezagff)
# Analysis of Deviance Table
# Cox model: response is  Surv(survyrs1999, survstat1999)
# Model 1: ~ h3k27me3_v1.1.pp_t6095f * ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# Model 2: ~ ezh2_blt10vge10_v1.pp * h3k27me3_v1.1.pp_t6095f * ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# loglik  Chisq Df P(>|Chi|)    
# 1 -7551.1                        
# 2 -7525.6 50.949 18 5.419e-05 ***     ##### H3K27me3 differences
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 

cmh3c3ezagfre <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp * ageg6 + 
                         grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[cch3c3ezag, ])
cmh3c3ezagfre
anova(cmh3c3ezagfre, cmh3c3ezagff)
# 
# > anova(cmh3c3ezagfre, cmh3c3ezagff)
# Analysis of Deviance Table
# Cox model: response is  Surv(survyrs1999, survstat1999)
# Model 1: ~ ezh2_blt10vge10_v1.pp * ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# Model 2: ~ ezh2_blt10vge10_v1.pp * h3k27me3_v1.1.pp_t6095f * ageg6 + grade_b12v3 + nodestat + size_tumor_grp
# loglik  Chisq Df P(>|Chi|)   
# 1 -7547.7                       
# 2 -7525.6 44.208 24  0.007216 **    ##### ezh2 differences  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > 


# h3k27me3_v1.1.pp_b60f ezh2_blt10vge10_v1.pp  ageg4  
# Basic survival analysis and graphs.


### Split by xxxx

table(wmadf$h3k27me3_v1.1.pp_b60f)
table(wmadf$ezh2_blt10vge10_v1.pp)
table(wmadf$ageg4)

names(wmadf)[which(grepl("ez", names(wmadf), ignore.case = TRUE))]

h3c2enag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(0,35]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(35,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(55,75]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(75,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2epag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(0,35]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(35,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(55,75]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(75,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
#logrank and pvals
h3c2enag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,35]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag1.lrpv <- lrpv(h3c2enag1.sdif)
h3c2enag1.lrpvtxt <- paste("p =", format(h3c2enag1.lrpv, digits = 3))
h3c2enag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(35,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag2.lrpv <- lrpv(h3c2enag2.sdif)
h3c2enag2.lrpvtxt <- paste("p =", format(h3c2enag2.lrpv, digits = 3))
h3c2enag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,75]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag3.lrpv <- lrpv(h3c2enag3.sdif)
h3c2enag3.lrpvtxt <- paste("p =", format(h3c2enag3.lrpv, digits = 3))
h3c2enag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(75,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag4.lrpv <- lrpv(h3c2enag4.sdif)
h3c2enag4.lrpvtxt <- paste("p =", format(h3c2enag4.lrpv, digits = 3))
h3c2epag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,35]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag1.lrpv <- lrpv(h3c2epag1.sdif)
h3c2epag1.lrpvtxt <- paste("p =", format(h3c2epag1.lrpv, digits = 3))
h3c2epag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(35,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag2.lrpv <- lrpv(h3c2epag2.sdif)
h3c2epag2.lrpvtxt <- paste("p =", format(h3c2epag2.lrpv, digits = 3))
h3c2epag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,75]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag3.lrpv <- lrpv(h3c2epag3.sdif)
h3c2epag3.lrpvtxt <- paste("p =", format(h3c2epag3.lrpv, digits = 3))
h3c2epag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(75,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag4.lrpv <- lrpv(h3c2epag4.sdif)
h3c2epag4.lrpvtxt <- paste("p =", format(h3c2epag4.lrpv, digits = 3))


ezc2hnag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(0,35]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(35,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(55,75]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(75,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(0,35]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(35,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(55,75]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(75,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
#logrank and pvals
ezc2hnag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(0,35]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag1.lrpv <- lrpv(ezc2hnag1.sdif)
ezc2hnag1.lrpvtxt <- paste("p =", format(ezc2hnag1.lrpv, digits = 3))
ezc2hnag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(35,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag2.lrpv <- lrpv(ezc2hnag2.sdif)
ezc2hnag2.lrpvtxt <- paste("p =", format(ezc2hnag2.lrpv, digits = 3))
ezc2hnag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(55,75]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag3.lrpv <- lrpv(ezc2hnag3.sdif)
ezc2hnag3.lrpvtxt <- paste("p =", format(ezc2hnag3.lrpv, digits = 3))
ezc2hnag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(75,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag4.lrpv <- lrpv(ezc2hnag4.sdif)
ezc2hnag4.lrpvtxt <- paste("p =", format(ezc2hnag4.lrpv, digits = 3))
ezc2hpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(0,35]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag1.lrpv <- lrpv(ezc2hpag1.sdif)
ezc2hpag1.lrpvtxt <- paste("p =", format(ezc2hpag1.lrpv, digits = 3))
ezc2hpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(35,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag2.lrpv <- lrpv(ezc2hpag2.sdif)
ezc2hpag2.lrpvtxt <- paste("p =", format(ezc2hpag2.lrpv, digits = 3))
ezc2hpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(55,75]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag3.lrpv <- lrpv(ezc2hpag3.sdif)
ezc2hpag3.lrpvtxt <- paste("p =", format(ezc2hpag3.lrpv, digits = 3))
ezc2hpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(75,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag4.lrpv <- lrpv(ezc2hpag4.sdif)
ezc2hpag4.lrpvtxt <- paste("p =", format(ezc2hpag4.lrpv, digits = 3))

par(mfcol=c(4, 4))
plot(h3c2enag1.surv, lty = 1:3, main = "(0,35] EZH2 neg")
text(x=0.5, y=0.1, labels = h3c2enag1.lrpvtxt, adj=0)
plot(h3c2enag2.surv, lty = 1:3, main = "(35,55] EZH2 neg")
text(x=0.5, y=0.1, labels = h3c2enag2.lrpvtxt, adj=0)
plot(h3c2enag3.surv, lty = 1:3, main = "(55,75] EZH2 neg")
text(x=0.5, y=0.1, labels = h3c2enag3.lrpvtxt, adj=0)
plot(h3c2enag4.surv, lty = 1:3, main = "(75,99] EZH2 neg")
text(x=0.5, y=0.1, labels = h3c2enag4.lrpvtxt, adj=0)
plot(h3c2epag1.surv, lty = 1:3, main = "(0,35] EZH2 pos")
text(x=0.5, y=0.1, labels = h3c2epag1.lrpvtxt, adj=0)
plot(h3c2epag2.surv, lty = 1:3, main = "(35,55] EZH2 pos")
text(x=0.5, y=0.1, labels = h3c2epag2.lrpvtxt, adj=0)
plot(h3c2epag3.surv, lty = 1:3, main = "(55,75] EZH2 pos")
text(x=0.5, y=0.1, labels = h3c2epag3.lrpvtxt, adj=0)
plot(h3c2epag4.surv, lty = 1:3, main = "(75,99] EZH2 pos")
text(x=0.5, y=0.1, labels = h3c2epag4.lrpvtxt, adj=0)

legend("topright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2)
title(main = "H3K27me3 low to high", line=3)

plot(ezc2hnag1.surv, lty = 1:3, main = "(0,35] H3K27me3 low")
text(x=0.5, y=0.1, labels = ezc2hnag1.lrpvtxt, adj=0)
plot(ezc2hnag2.surv, lty = 1:3, main = "(35,55] H3K27me3 low")
text(x=0.5, y=0.1, labels = ezc2hnag2.lrpvtxt, adj=0)
plot(ezc2hnag3.surv, lty = 1:3, main = "(55,75] H3K27me3 low")
text(x=0.5, y=0.1, labels = ezc2hnag3.lrpvtxt, adj=0)
plot(ezc2hnag4.surv, lty = 1:3, main = "(75,99] H3K27me3 low")
text(x=0.5, y=0.1, labels = ezc2hnag4.lrpvtxt, adj=0)
plot(ezc2hpag1.surv, lty = 1:3, main = "(0,35] H3K27me3 high")
text(x=0.5, y=0.1, labels = ezc2hpag1.lrpvtxt, adj=0)
plot(ezc2hpag2.surv, lty = 1:3, main = "(35,55] H3K27me3 high")
text(x=0.5, y=0.1, labels = ezc2hpag2.lrpvtxt, adj=0)
plot(ezc2hpag3.surv, lty = 1:3, main = "(55,75] H3K27me3 high")
text(x=0.5, y=0.1, labels = ezc2hpag3.lrpvtxt, adj=0)
plot(ezc2hpag4.surv, lty = 1:3, main = "(75,99] H3K27me3 high")
text(x=0.5, y=0.1, labels = ezc2hpag4.lrpvtxt, adj=0)

legend("topright",legend = levels(as.factor(wmadf$ezh2_blt10vge10_v1.pp))  , lty = 1:2)
title(main = "EZH2 low to high", line=3)

## Cox


cch3c2ezag <- with(wmadf, complete.cases(h3k27me3_v1.1.pp_b60f, ezh2_blt10vge10_v1.pp, 
                                         ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(cch3c2ezag)
head(cch3c2ezag)

cmh3c2ezagff <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp * h3k27me3_v1.1.pp_b60f * ageg4 + 
                        grade_b12v3 + nodestat + size_tumor_grp, 
                      data = wmadf[cch3c2ezag, ])
cmh3c2ezagff


### Split by xxxx


wmadf$ageg4 <- cut(wmadf$age, breaks = c(0, 40, 55, 70, 99))
table(wmadf$ageg4)


h3c2ag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                        subset = ((ageg4 == "(0,40]"  ) ) )
h3c2ag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                        subset = ((ageg4 == "(40,55]"  ) ) )
h3c2ag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                        subset = ((ageg4 == "(55,70]"  ) ) )
h3c2ag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                        subset = ((ageg4 == "(70,99]"  ) ) )
#logrank and pvals
h3c2ag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                         subset = ((ageg4 == "(0,40]"  ) ) )
h3c2ag1.lrpv <- lrpv(h3c2ag1.sdif)
h3c2ag1.lrpvtxt <- paste("p =", format(h3c2ag1.lrpv, digits = 3))
h3c2ag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                         subset = ((ageg4 == "(40,55]"  ) ) )
h3c2ag2.lrpv <- lrpv(h3c2ag2.sdif)
h3c2ag2.lrpvtxt <- paste("p =", format(h3c2ag2.lrpv, digits = 3))
h3c2ag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                         subset = ((ageg4 == "(55,70]"  ) ) )
h3c2ag3.lrpv <- lrpv(h3c2ag3.sdif)
h3c2ag3.lrpvtxt <- paste("p =", format(h3c2ag3.lrpv, digits = 3))
h3c2ag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                         subset = ((ageg4 == "(70,99]"  ) ) )
h3c2ag4.lrpv <- lrpv(h3c2ag4.sdif)
h3c2ag4.lrpvtxt <- paste("p =", format(h3c2ag4.lrpv, digits = 3))


## Stratify by EZH2 low/high
h3c2enag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2epag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                          subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
#logrank and pvals
h3c2enag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag1.lrpv <- lrpv(h3c2enag1.sdif)
h3c2enag1.lrpvtxt <- paste("p =", format(h3c2enag1.lrpv, digits = 3))
h3c2enag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag2.lrpv <- lrpv(h3c2enag2.sdif)
h3c2enag2.lrpvtxt <- paste("p =", format(h3c2enag2.lrpv, digits = 3))
h3c2enag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag3.lrpv <- lrpv(h3c2enag3.sdif)
h3c2enag3.lrpvtxt <- paste("p =", format(h3c2enag3.lrpv, digits = 3))
h3c2enag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%")) )
h3c2enag4.lrpv <- lrpv(h3c2enag4.sdif)
h3c2enag4.lrpvtxt <- paste("p =", format(h3c2enag4.lrpv, digits = 3))
h3c2epag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag1.lrpv <- lrpv(h3c2epag1.sdif)
h3c2epag1.lrpvtxt <- paste("p =", format(h3c2epag1.lrpv, digits = 3))
h3c2epag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag2.lrpv <- lrpv(h3c2epag2.sdif)
h3c2epag2.lrpvtxt <- paste("p =", format(h3c2epag2.lrpv, digits = 3))
h3c2epag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag3.lrpv <- lrpv(h3c2epag3.sdif)
h3c2epag3.lrpvtxt <- paste("p =", format(h3c2epag3.lrpv, digits = 3))
h3c2epag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%")) )
h3c2epag4.lrpv <- lrpv(h3c2epag4.sdif)
h3c2epag4.lrpvtxt <- paste("p =", format(h3c2epag4.lrpv, digits = 3))


ezc2ag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                        subset = ((ageg4 == "(0,40]"  ) ) )
ezc2ag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                        subset = ((ageg4 == "(40,55]"  ) ) )
ezc2ag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                        subset = ((ageg4 == "(55,70]"  ) ) )
ezc2ag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                        subset = ((ageg4 == "(70,99]"  ) ) )
#logrank and pvals
ezc2ag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                         subset = ((ageg4 == "(0,40]"  ) ) )
ezc2ag1.lrpv <- lrpv(ezc2ag1.sdif)
ezc2ag1.lrpvtxt <- paste("p =", format(ezc2ag1.lrpv, digits = 3))
ezc2ag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                         subset = ((ageg4 == "(40,55]"  ) ) )
ezc2ag2.lrpv <- lrpv(ezc2ag2.sdif)
ezc2ag2.lrpvtxt <- paste("p =", format(ezc2ag2.lrpv, digits = 3))
ezc2ag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                         subset = ((ageg4 == "(55,70]"  ) ) )
ezc2ag3.lrpv <- lrpv(ezc2ag3.sdif)
ezc2ag3.lrpvtxt <- paste("p =", format(ezc2ag3.lrpv, digits = 3))
ezc2ag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                         subset = ((ageg4 == "(70,99]"  ) ) )
ezc2ag4.lrpv <- lrpv(ezc2ag4.sdif)
ezc2ag4.lrpvtxt <- paste("p =", format(ezc2ag4.lrpv, digits = 3))

## Stratify by H3K27me3
ezc2hnag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                          subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
#logrank and pvals
ezc2hnag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag1.lrpv <- lrpv(ezc2hnag1.sdif)
ezc2hnag1.lrpvtxt <- paste("p =", format(ezc2hnag1.lrpv, digits = 3))
ezc2hnag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag2.lrpv <- lrpv(ezc2hnag2.sdif)
ezc2hnag2.lrpvtxt <- paste("p =", format(ezc2hnag2.lrpv, digits = 3))
ezc2hnag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag3.lrpv <- lrpv(ezc2hnag3.sdif)
ezc2hnag3.lrpvtxt <- paste("p =", format(ezc2hnag3.lrpv, digits = 3))
ezc2hnag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low")) )
ezc2hnag4.lrpv <- lrpv(ezc2hnag4.sdif)
ezc2hnag4.lrpvtxt <- paste("p =", format(ezc2hnag4.lrpv, digits = 3))
ezc2hpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag1.lrpv <- lrpv(ezc2hpag1.sdif)
ezc2hpag1.lrpvtxt <- paste("p =", format(ezc2hpag1.lrpv, digits = 3))
ezc2hpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag2.lrpv <- lrpv(ezc2hpag2.sdif)
ezc2hpag2.lrpvtxt <- paste("p =", format(ezc2hpag2.lrpv, digits = 3))
ezc2hpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag3.lrpv <- lrpv(ezc2hpag3.sdif)
ezc2hpag3.lrpvtxt <- paste("p =", format(ezc2hpag3.lrpv, digits = 3))
ezc2hpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high")) )
ezc2hpag4.lrpv <- lrpv(ezc2hpag4.sdif)
ezc2hpag4.lrpvtxt <- paste("p =", format(ezc2hpag4.lrpv, digits = 3))


pdf(file = paste(plot_dirs, "EZH2_H3K27me3_ER_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")

pdf(file = paste(plot_dirs, "EZH2_H3K27me3_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "EZH2_H3K27me3_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  plot(h3c2enag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40 EZH2 neg", paste(": N =", paste(h3c2enag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] EZH2 neg", paste(": N =", paste(h3c2enag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enag2.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] EZH2 neg", paste(": N =", paste(h3c2enag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] EZH2 neg", paste(": N =", paste(h3c2enag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enag4.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2epag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] EZH2 pos", paste(": N =", paste(h3c2epag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2epag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2epag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] EZH2 pos", paste(": N =", paste(h3c2epag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2epag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2epag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] EZH2 pos", paste(": N =", paste(h3c2epag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2epag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2epag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] EZH2 pos", paste(": N =", paste(h3c2epag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2epag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(ezc2hnag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] H3K27me3 low", paste(": N =", paste(ezc2hnag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hnag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] H3K27me3 low", paste(": N =", paste(ezc2hnag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = ezc2hnag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(ezc2hnag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] H3K27me3 low", paste(": N =", paste(ezc2hnag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hnag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] H3K27me3 low", paste(": N =", paste(ezc2hnag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnag4.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] H3K27me3 high", paste(": N =", paste(ezc2hpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] H3K27me3 high", paste(": N =", paste(ezc2hpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpag2.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] H3K27me3 high", paste(": N =", paste(ezc2hpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] H3K27me3 high", paste(": N =", paste(ezc2hpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = paste("Strata Age, EZH2: N =", sum(c(h3c2enag1.surv$n, h3c2enag2.surv$n, h3c2enag3.surv$n, h3c2enag4.surv$n,
                                                     h3c2epag1.surv$n, h3c2epag2.surv$n, h3c2epag3.surv$n, h3c2epag4.surv$n), na.rm = TRUE)), 
                     line=1, outer = TRUE, adj = 0.17, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2enag1.surv$n, h3c2enag2.surv$n, h3c2enag3.surv$n, h3c2enag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.05, cex = 0.7)
  mtext(text = paste("N =", sum(c(h3c2epag1.surv$n, h3c2epag2.surv$n, h3c2epag3.surv$n, h3c2epag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.43, cex = 0.7)
  
  mtext(text = paste("Strata Age, H3K27me3: N =", sum(c(ezc2hnag1.surv$n, ezc2hnag2.surv$n, ezc2hnag3.surv$n, ezc2hnag4.surv$n,
                                                         ezc2hpag1.surv$n, ezc2hpag2.surv$n, ezc2hpag3.surv$n, ezc2hpag4.surv$n), na.rm = TRUE)),
        line=1, outer = TRUE, adj = 0.9, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2hnag1.surv$n, ezc2hnag2.surv$n, ezc2hnag3.surv$n, ezc2hnag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.58, cex = 0.7)
  mtext(text = paste("N =", sum(c(ezc2hpag1.surv$n, ezc2hpag2.surv$n, ezc2hpag3.surv$n, ezc2hpag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.97, cex = 0.7)
  
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()




## ER+  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

h3c2enerpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2enerpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2enerpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2enerpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
#logrank and pvals
h3c2enerpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2enerpag1.lrpv <- lrpv(h3c2enerpag1.sdif)
h3c2enerpag1.lrpvtxt <- paste("p =", format(h3c2enerpag1.lrpv, digits = 3))
h3c2enerpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2enerpag2.lrpv <- lrpv(h3c2enerpag2.sdif)
h3c2enerpag2.lrpvtxt <- paste("p =", format(h3c2enerpag2.lrpv, digits = 3))
h3c2enerpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2enerpag3.lrpv <- lrpv(h3c2enerpag3.sdif)
h3c2enerpag3.lrpvtxt <- paste("p =", format(h3c2enerpag3.lrpv, digits = 3))
h3c2enerpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2enerpag4.lrpv <- lrpv(h3c2enerpag4.sdif)
h3c2enerpag4.lrpvtxt <- paste("p =", format(h3c2enerpag4.lrpv, digits = 3))
h3c2eperpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag1.lrpv <- lrpv(h3c2eperpag1.sdif)
h3c2eperpag1.lrpvtxt <- paste("p =", format(h3c2eperpag1.lrpv, digits = 3))
h3c2eperpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag2.lrpv <- lrpv(h3c2eperpag2.sdif)
h3c2eperpag2.lrpvtxt <- paste("p =", format(h3c2eperpag2.lrpv, digits = 3))
h3c2eperpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag3.lrpv <- lrpv(h3c2eperpag3.sdif)
h3c2eperpag3.lrpvtxt <- paste("p =", format(h3c2eperpag3.lrpv, digits = 3))
h3c2eperpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2eperpag4.lrpv <- lrpv(h3c2eperpag4.sdif)
h3c2eperpag4.lrpvtxt <- paste("p =", format(h3c2eperpag4.lrpv, digits = 3))


ezc2hnerpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hnerpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hnerpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hnerpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
#logrank and pvals
ezc2hnerpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hnerpag1.lrpv <- lrpv(ezc2hnerpag1.sdif)
ezc2hnerpag1.lrpvtxt <- paste("p =", format(ezc2hnerpag1.lrpv, digits = 3))
ezc2hnerpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hnerpag2.lrpv <- lrpv(ezc2hnerpag2.sdif)
ezc2hnerpag2.lrpvtxt <- paste("p =", format(ezc2hnerpag2.lrpv, digits = 3))
ezc2hnerpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hnerpag3.lrpv <- lrpv(ezc2hnerpag3.sdif)
ezc2hnerpag3.lrpvtxt <- paste("p =", format(ezc2hnerpag3.lrpv, digits = 3))
ezc2hnerpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hnerpag4.lrpv <- lrpv(ezc2hnerpag4.sdif)
ezc2hnerpag4.lrpvtxt <- paste("p =", format(ezc2hnerpag4.lrpv, digits = 3))
ezc2hperpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag1.lrpv <- lrpv(ezc2hperpag1.sdif)
ezc2hperpag1.lrpvtxt <- paste("p =", format(ezc2hperpag1.lrpv, digits = 3))
ezc2hperpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag2.lrpv <- lrpv(ezc2hperpag2.sdif)
ezc2hperpag2.lrpvtxt <- paste("p =", format(ezc2hperpag2.lrpv, digits = 3))
ezc2hperpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag3.lrpv <- lrpv(ezc2hperpag3.sdif)
ezc2hperpag3.lrpvtxt <- paste("p =", format(ezc2hperpag3.lrpv, digits = 3))
ezc2hperpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2hperpag4.lrpv <- lrpv(ezc2hperpag4.sdif)
ezc2hperpag4.lrpvtxt <- paste("p =", format(ezc2hperpag4.lrpv, digits = 3))


pdf(file = paste(plot_dirs, "EZH2_H3K27me3_ERpos_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "EZH2_H3K27me3_ERpos_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  plot(h3c2enerpag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40 EZH2 neg", paste(": N =", paste(h3c2enerpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enerpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enerpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] EZH2 neg", paste(": N =", paste(h3c2enerpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enerpag2.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enerpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] EZH2 neg", paste(": N =", paste(h3c2enerpag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enerpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enerpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] EZH2 neg", paste(": N =", paste(h3c2enerpag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enerpag4.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2eperpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] EZH2 pos", paste(": N =", paste(h3c2eperpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2eperpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2eperpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] EZH2 pos", paste(": N =", paste(h3c2eperpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2eperpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2eperpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] EZH2 pos", paste(": N =", paste(h3c2eperpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2eperpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2eperpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] EZH2 pos", paste(": N =", paste(h3c2eperpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2eperpag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(ezc2hnerpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] H3K27me3 low", paste(": N =", paste(ezc2hnerpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnerpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hnerpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] H3K27me3 low", paste(": N =", paste(ezc2hnerpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = ezc2hnerpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(ezc2hnerpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] H3K27me3 low", paste(": N =", paste(ezc2hnerpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnerpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hnerpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] H3K27me3 low", paste(": N =", paste(ezc2hnerpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnerpag4.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hperpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] H3K27me3 high", paste(": N =", paste(ezc2hperpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hperpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hperpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] H3K27me3 high", paste(": N =", paste(ezc2hperpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hperpag2.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hperpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] H3K27me3 high", paste(": N =", paste(ezc2hperpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hperpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hperpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] H3K27me3 high", paste(": N =", paste(ezc2hperpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hperpag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = paste("Strata ER+, Age, EZH2: N =", sum(c(h3c2enerpag1.surv$n, h3c2enerpag2.surv$n, h3c2enerpag3.surv$n, h3c2enerpag4.surv$n,
                                                          h3c2eperpag1.surv$n, h3c2eperpag2.surv$n, h3c2eperpag3.surv$n, h3c2eperpag4.surv$n), na.rm = TRUE)), 
        line=1, outer = TRUE, adj = 0.12, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2enerpag1.surv$n, h3c2enerpag2.surv$n, h3c2enerpag3.surv$n, h3c2enerpag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.05, cex = 0.7)
  mtext(text = paste("N =", sum(c(h3c2eperpag1.surv$n, h3c2eperpag2.surv$n, h3c2eperpag3.surv$n, h3c2eperpag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.43, cex = 0.7)
  
  
  mtext(text = paste("Strata ER+, Age, H3K27me3: N =", sum(c(ezc2hnerpag1.surv$n, ezc2hnerpag2.surv$n, ezc2hnerpag3.surv$n, ezc2hnerpag4.surv$n,
                                                              ezc2hperpag1.surv$n, ezc2hperpag2.surv$n, ezc2hperpag3.surv$n, ezc2hperpag4.surv$n))), 
        line=1, outer = TRUE, adj = 0.96, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2hnerpag1.surv$n, ezc2hnerpag2.surv$n, ezc2hnerpag3.surv$n, ezc2hnerpag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.58, cex = 0.7)
  mtext(text = paste("N =", sum(c(ezc2hperpag1.surv$n, ezc2hperpag2.surv$n, ezc2hperpag3.surv$n, ezc2hperpag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.97, cex = 0.7)
  
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()


## ER-  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

h3c2enernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2enernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2enernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2enernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2epernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
h3c2epernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
h3c2epernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
h3c2epernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
#logrank and pvals
h3c2enernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2enernag1.lrpv <- lrpv(h3c2enernag1.sdif)
h3c2enernag1.lrpvtxt <- paste("p =", format(h3c2enernag1.lrpv, digits = 3))
h3c2enernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2enernag2.lrpv <- lrpv(h3c2enernag2.sdif)
h3c2enernag2.lrpvtxt <- paste("p =", format(h3c2enernag2.lrpv, digits = 3))
h3c2enernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2enernag3.lrpv <- lrpv(h3c2enernag3.sdif)
h3c2enernag3.lrpvtxt <- paste("p =", format(h3c2enernag3.lrpv, digits = 3))
h3c2enernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "0 to 5%") & (er_b0v123_v2 == "negative")) )
h3c2enernag4.lrpv <- lrpv(h3c2enernag4.sdif)
h3c2enernag4.lrpvtxt <- paste("p =", format(h3c2enernag4.lrpv, digits = 3))
h3c2epernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
h3c2epernag1.lrpv <- lrpv(h3c2epernag1.sdif)
h3c2epernag1.lrpvtxt <- paste("p =", format(h3c2epernag1.lrpv, digits = 3))
h3c2epernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
h3c2epernag2.lrpv <- lrpv(h3c2epernag2.sdif)
h3c2epernag2.lrpvtxt <- paste("p =", format(h3c2epernag2.lrpv, digits = 3))
h3c2epernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
h3c2epernag3.lrpv <- lrpv(h3c2epernag3.sdif)
h3c2epernag3.lrpvtxt <- paste("p =", format(h3c2epernag3.lrpv, digits = 3))
h3c2epernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (ezh2_blt10vge10_v1.pp == "10 to 100%") & (er_b0v123_v2 == "negative")) )
h3c2epernag4.lrpv <- lrpv(h3c2epernag4.sdif)
h3c2epernag4.lrpvtxt <- paste("p =", format(h3c2epernag4.lrpv, digits = 3))


ezc2hnernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hnernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hnernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hnernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hpernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
ezc2hpernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
ezc2hpernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
ezc2hpernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
#logrank and pvals
ezc2hnernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hnernag1.lrpv <- lrpv(ezc2hnernag1.sdif)
ezc2hnernag1.lrpvtxt <- paste("p =", format(ezc2hnernag1.lrpv, digits = 3))
ezc2hnernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hnernag2.lrpv <- lrpv(ezc2hnernag2.sdif)
ezc2hnernag2.lrpvtxt <- paste("p =", format(ezc2hnernag2.lrpv, digits = 3))
ezc2hnernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hnernag3.lrpv <- lrpv(ezc2hnernag3.sdif)
ezc2hnernag3.lrpvtxt <- paste("p =", format(ezc2hnernag3.lrpv, digits = 3))
ezc2hnernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 low") & (er_b0v123_v2 == "negative")) )
ezc2hnernag4.lrpv <- lrpv(ezc2hnernag4.sdif)
ezc2hnernag4.lrpvtxt <- paste("p =", format(ezc2hnernag4.lrpv, digits = 3))
ezc2hpernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
ezc2hpernag1.lrpv <- lrpv(ezc2hpernag1.sdif)
ezc2hpernag1.lrpvtxt <- paste("p =", format(ezc2hpernag1.lrpv, digits = 3))
ezc2hpernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
ezc2hpernag2.lrpv <- lrpv(ezc2hpernag2.sdif)
ezc2hpernag2.lrpvtxt <- paste("p =", format(ezc2hpernag2.lrpv, digits = 3))
ezc2hpernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
ezc2hpernag3.lrpv <- lrpv(ezc2hpernag3.sdif)
ezc2hpernag3.lrpvtxt <- paste("p =", format(ezc2hpernag3.lrpv, digits = 3))
ezc2hpernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (h3k27me3_v1.1.pp_b60f == "H3K27me3 high") & (er_b0v123_v2 == "negative")) )
ezc2hpernag4.lrpv <- lrpv(ezc2hpernag4.sdif)
ezc2hpernag4.lrpvtxt <- paste("p =", format(ezc2hpernag4.lrpv, digits = 3))


pdf(file = paste(plot_dirs, "EZH2_H3K27me3_ERneg_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "EZH2_H3K27me3_ERneg_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  plot(h3c2enernag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40 EZH2 neg", paste(": N =", paste(h3c2enernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] EZH2 neg", paste(": N =", paste(h3c2enernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enernag2.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] EZH2 neg", paste(": N =", paste(h3c2enernag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2enernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] EZH2 neg", paste(": N =", paste(h3c2enernag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2enernag4.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2epernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] EZH2 pos", paste(": N =", paste(h3c2epernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2epernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2epernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] EZH2 pos", paste(": N =", paste(h3c2epernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2epernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2epernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] EZH2 pos", paste(": N =", paste(h3c2epernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2epernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2epernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] EZH2 pos", paste(": N =", paste(h3c2epernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2epernag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(ezc2hnernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] H3K27me3 low", paste(": N =", paste(ezc2hnernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hnernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] H3K27me3 low", paste(": N =", paste(ezc2hnernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = ezc2hnernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(ezc2hnernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] H3K27me3 low", paste(": N =", paste(ezc2hnernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hnernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] H3K27me3 low", paste(": N =", paste(ezc2hnernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hnernag4.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] H3K27me3 high", paste(": N =", paste(ezc2hpernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] H3K27me3 high", paste(": N =", paste(ezc2hpernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpernag2.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] H3K27me3 high", paste(": N =", paste(ezc2hpernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2hpernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] H3K27me3 high", paste(": N =", paste(ezc2hpernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.75, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = ezc2hpernag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = paste("Strata ER- : Age, EZH2: N =", sum(c(h3c2enernag1.surv$n, h3c2enernag2.surv$n, h3c2enernag3.surv$n, h3c2enernag4.surv$n,
                                                          h3c2epernag1.surv$n, h3c2epernag2.surv$n, h3c2epernag3.surv$n, h3c2epernag4.surv$n), na.rm = TRUE)), 
        line=1, outer = TRUE, adj = 0.12, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2enernag1.surv$n, h3c2enernag2.surv$n, h3c2enernag3.surv$n, h3c2enernag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.05, cex = 0.7)
  mtext(text = paste("N =", sum(c(h3c2epernag1.surv$n, h3c2epernag2.surv$n, h3c2epernag3.surv$n, h3c2epernag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.43, cex = 0.7)
  
  
  mtext(text = paste("Strata ER- : Age, H3K27me3: N =", sum(c(ezc2hnernag1.surv$n, ezc2hnernag2.surv$n, ezc2hnernag3.surv$n, ezc2hnernag4.surv$n,
                                                              ezc2hpernag1.surv$n, ezc2hpernag2.surv$n, ezc2hpernag3.surv$n, ezc2hpernag4.surv$n))), 
        line=1, outer = TRUE, adj = 0.96, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2hnernag1.surv$n, ezc2hnernag2.surv$n, ezc2hnernag3.surv$n, ezc2hnernag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.58, cex = 0.7)
  mtext(text = paste("N =", sum(c(ezc2hpernag1.surv$n, ezc2hpernag2.surv$n, ezc2hpernag3.surv$n, ezc2hpernag4.surv$n), na.rm = TRUE)), line=-0.5, outer = TRUE, adj = 0.97, cex = 0.7)
  
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()

# erpag # Hexadecimal color specification 
# brewer.pal(n = 8, name = "Dark2")
# 
# ## [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
# ## [8] "#666666"
# 
# cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# p + scale_color_manual(values = cbp1)

## Colourblind panels of colours
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Dark2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
# p + scale_color_manual(values = cbp1)


library(ggplot2)
library(viridis)
library(hexbin)
barplot(1:8, col = Dark2)
barplot(1:5, col=topo.colors(5))
ggplot(data.frame(x = rnorm(10000), y = rnorm(10000)), aes(x = x, y = y)) +
  geom_hex() + coord_fixed() +
  scale_fill_manual(values = Dark2) + theme_bw()
#scale_color_manual(values = cbp1) + theme_bw()
#scale_fill_viridis() + theme_bw()


### Cox models

cccmehagntp <- with(wmadf, complete.cases(h3k27me3_v1.1.pp_t6095f, ezh2_blt10vge10_v1.pp, er_b0v123_v2, her2_b012v23_v2, # adding her2 yields fit problems in full interaction model
                                          ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(cccmehagntp)

cmehagntm1 <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp + h3k27me3_v1.1.pp_b60f + 
                      ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmehagntp, ])
print(cmehagntm1)
summary(cmehagntm1)

cmehagntm2 <- coxph(Surv(survyrs1999, survstat1999) ~ ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmehagntp, ])
anova(cmehagntm2, cmehagntm1)

cmehagntm3 <- coxph(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp * h3k27me3_v1.1.pp_b60f * ageg4 + 
                      grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmehagntp, ])
print(cmehagntm3)
summary(cmehagntm3)

anova(cmehagntm3, cmehagntm1)
anova(cmehagntm3, cmehagntm2)

#er_b0v123_v2

cmehagntm4 <- coxph(Surv(survyrs1999, survstat1999) ~ er_b0v123_v2 * ezh2_blt10vge10_v1.pp * h3k27me3_v1.1.pp_b60f * ageg4 + 
                      grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmehagntp, ])
print(cmehagntm4)
summary(cmehagntm4)

anova(cmehagntm4, cmehagntm3)
anova(cmehagntm4, cmehagntm2)

cmehagntm5 <- coxph(Surv(survyrs1999, survstat1999) ~ er_b0v123_v2 + ezh2_blt10vge10_v1.pp + h3k27me3_v1.1.pp_b60f + ageg4 + 
                      grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmehagntp, ])
summary(cmehagntm5)
anova(cmehagntm5, cmehagntm2)


# her2_b012v23_v2 can not additionally be accommodated with interactions

# cmehagntm5 <- coxph(Surv(survyrs1999, survstat1999) ~ 
#                       er_b0v123_v2 * her2_b012v23_v2 * ezh2_blt10vge10_v1.pp * h3k27me3_v1.1.pp_b60f * ageg4 + 
#                       grade_b12v3 + nodestat + size_tumor_grp, 
#                     data = wmadf[cccmehagntp, ])
# Warning message:
#   In fitter(X, Y, strats, offset, init, control, weights = weights,  :
#               Loglik converged before variable  14,30,33,40,41,42,52,53,54,55,58,62,63,64,65,66,67 ; beta may be infinite. 
#             
# print(cmehagntm5)
# summary(cmehagntm5)

# Do additive no interaction model with HER2

cmehagntm6 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2 + ezh2_blt10vge10_v1.pp + h3k27me3_v1.1.pp_b60f + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmehagntp, ])

summary(cmehagntm6)
anova(cmehagntm6, cmehagntm2)


cmehagntm7 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2 + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmehagntp, ])

summary(cmehagntm7)
anova(cmehagntm7, cmehagntm6)


cmehagntm8 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2 + ezh2_blt10vge10_v1.pp + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmehagntp, ])

summary(cmehagntm8)
anova(cmehagntm8, cmehagntm7)


cmehagntm9 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2 + h3k27me3_v1.1.pp_b60f + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmehagntp, ])

summary(cmehagntm9)
anova(cmehagntm9, cmehagntm7)






#####  FOXA1  ######

## wmadf$h3k27me3_v1.1n -->> wmadf$foxa1_v1.ppnc
## h3k27me3_v1.1.pp -->> foxa1_v1.ppnc
## H3K27me3 -->> FOXA1
## h3v11 -->> fx1v1


wmadf$foxa1_v1.ppnc <- as.numeric(wmadf$foxa1_v1.percent)
table(wmadf$foxa1_v1.ppnc, useNA = "always")
## Original FOXA1 analysis done on the odd classification by image colour in FOXA1_v1
## Get some info on the percent positive version.

splitpctvec <- (1:10)*10
splitpct <- 90

wmadf$foxa1_v1.ppnc_b10c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b20c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b30c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b40c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b50c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b60c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b70c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b80c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b90c <- rep(NA_character_, nrow(wmadf))
wmadf$foxa1_v1.ppnc_b100c <- rep(NA_character_, nrow(wmadf))

wmadf$foxa1_v1.ppnc_b10c[wmadf$foxa1_v1.ppnc <  10] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b10c[wmadf$foxa1_v1.ppnc >= 10] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b10f <- factor(wmadf$foxa1_v1.ppnc_b10c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b20c[wmadf$foxa1_v1.ppnc <  20] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b20c[wmadf$foxa1_v1.ppnc >= 20] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b20f <- factor(wmadf$foxa1_v1.ppnc_b20c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b30c[wmadf$foxa1_v1.ppnc <  30] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b30c[wmadf$foxa1_v1.ppnc >= 30] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b30f <- factor(wmadf$foxa1_v1.ppnc_b30c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b40c[wmadf$foxa1_v1.ppnc <  40] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b40c[wmadf$foxa1_v1.ppnc >= 40] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b40f <- factor(wmadf$foxa1_v1.ppnc_b40c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b50c[wmadf$foxa1_v1.ppnc <  50] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b50c[wmadf$foxa1_v1.ppnc >= 50] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b50f <- factor(wmadf$foxa1_v1.ppnc_b50c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b60c[wmadf$foxa1_v1.ppnc <  60] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b60c[wmadf$foxa1_v1.ppnc >= 60] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b60f <- factor(wmadf$foxa1_v1.ppnc_b60c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b70c[wmadf$foxa1_v1.ppnc <  70] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b70c[wmadf$foxa1_v1.ppnc >= 70] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b70f <- factor(wmadf$foxa1_v1.ppnc_b70c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b80c[wmadf$foxa1_v1.ppnc <  80] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b80c[wmadf$foxa1_v1.ppnc >= 80] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b80f <- factor(wmadf$foxa1_v1.ppnc_b80c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b90c[wmadf$foxa1_v1.ppnc <  90] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b90c[wmadf$foxa1_v1.ppnc >= 90] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b90f <- factor(wmadf$foxa1_v1.ppnc_b90c, levels = c("FOXA1 low", "FOXA1 high"))
wmadf$foxa1_v1.ppnc_b100c[wmadf$foxa1_v1.ppnc <  100] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_b100c[wmadf$foxa1_v1.ppnc >= 100] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_b100f <- factor(wmadf$foxa1_v1.ppnc_b100c, levels = c("FOXA1 low", "FOXA1 high"))

table(wmadf$foxa1_v1.ppnc_b10f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b20f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b30f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b40f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b50f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b60f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b70f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b80f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b90f , useNA = "always")
table(wmadf$foxa1_v1.ppnc_b100f , useNA = "always")
### Split by EZH2

table(wmadf$FOXA1_v1_cat5)
table(wmadf$foxa1_v1.ppnc)
table(wmadf$foxa1_v1.ppnc_b60f)
table(wmadf$foxa1_v1.ppnc_b70f)
table(wmadf$foxa1_v1.ppnc_b80f)
table(wmadf$foxa1_v1.ppnc_b90f)
table(wmadf$foxa1_v1.ppnc_b100f)
table(wmadf$ezh2_blt10vge10_v1.pp)


fx1v1b10.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b10f , data = wmadf)
fx1v1b20.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b20f , data = wmadf)
fx1v1b30.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b30f , data = wmadf)
fx1v1b40.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b40f , data = wmadf)
fx1v1b50.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b50f , data = wmadf)
fx1v1b60.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf)
fx1v1b70.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b70f , data = wmadf)
fx1v1b80.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b80f , data = wmadf)
fx1v1b90.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b90f , data = wmadf)
fx1v1b100.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b100f , data = wmadf)

{par(mfrow=c(3, 4))
  plot(fx1v1b10.surv, lty = 1:2, main = "FOXA1 < vs >= 10")
  title(main = "FOXA1 low to high", line=2.2)
  plot(fx1v1b20.surv, lty = 1:2, main = "20")
  plot(fx1v1b30.surv, lty = 1:2, main = "30")
  plot(fx1v1b40.surv, lty = 1:2, main = "40")
  plot(fx1v1b50.surv, lty = 1:2, main = "50")
  plot(fx1v1b60.surv, lty = 1:2, main = "60")
  plot(fx1v1b70.surv, lty = 1:2, main = "70")
  plot(fx1v1b80.surv, lty = 1:2, main = "80")
  plot(fx1v1b90.surv, lty = 1:2, main = "90")
  plot(fx1v1b100.surv, lty = 1:2, main = "100")
  
  legend("bottomright",legend = levels(as.factor(wmadf$foxa1_v1.ppnc_b100f))  , lty = 1:2)
}

## Cox model fits
table(wmadf$ageg4)

ccfx1c2ag4 <- with(wmadf, complete.cases(foxa1_v1.ppnc_b60f, ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
ccfx1c2ag4 <- with(wmadf, complete.cases(foxa1_v1.ppnc_b90f, ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(ccfx1c2ag4)
head(ccfx1c2ag4)

cmfx1b10c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b10f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b10c3agff)
cmfx1b10c3agff$loglik  ##   [1] -9680.396 -9449.596

cmfx1b20c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b20f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b20c3agff)
cmfx1b20c3agff$loglik  ##   [1] -9680.396 -9447.652

cmfx1b30c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b30f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b30c3agff)
cmfx1b30c3agff$loglik  ##   [1] -9680.396 -9446.971

cmfx1b40c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b40f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b40c3agff)
cmfx1b40c3agff$loglik  ##   [1] -9680.396 -9448.858

cmfx1b50c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b50f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b50c3agff)
cmfx1b50c3agff$loglik  ##   [1] -9680.396 -9447.717

cmfx1b60c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b60c3agff)
cmfx1b60c3agff$loglik  ##  [1] -9680.396 -9448.592

cmfx1b70c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b70f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b70c3agff)
cmfx1b70c3agff$loglik  ##  [1] -9680.396 -9453.184

cmfx1b80c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b80f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b80c3agff)
cmfx1b80c3agff$loglik  ##  [1] -9680.396 -9454.080

cmfx1b90c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b90f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b90c3agff)
cmfx1b90c3agff$loglik ## [1] -9680.396 -9455.887

cmfx1b100c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b100f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                         data = wmadf[ccfx1c2ag4, ])
summary(cmfx1b100c3agff)
cmfx1b100c3agff$loglik  ## [1] -9680.396 -9455.013

## 60 pct cells positively stained is binary optimum. <60 vs >=60  foxa1_v1.ppnc_b60f 

##

wmadf$foxa1_v1.ppnc_t60100c <- rep(NA_character_, nrow(wmadf))

wmadf$foxa1_v1.ppnc_t60100c[wmadf$foxa1_v1.ppnc <  60] <- "FOXA1 low"
wmadf$foxa1_v1.ppnc_t60100c[wmadf$foxa1_v1.ppnc >= 60 & wmadf$foxa1_v1.ppnc < 100] <- "FOXA1 mid"
wmadf$foxa1_v1.ppnc_t60100c[wmadf$foxa1_v1.ppnc >= 100] <- "FOXA1 high"
wmadf$foxa1_v1.ppnc_t60100f <- factor(wmadf$foxa1_v1.ppnc_t60100c, levels = c("FOXA1 low", "FOXA1 mid", "FOXA1 high"))
table(wmadf$foxa1_v1.ppnc_t60100f)

fx1v1t60100.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_t60100f , data = wmadf)
plot(fx1v1t60100.surv, lty = 1:3, main = "[0,60),[60-100),[100,100]")

cmfx1v1t60100c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_t60100f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                             data = wmadf[ccfx1c2ag4, ])
summary(cmfx1v1t60100c3agff)
cmfx1v1t60100c3agff$loglik  ##  [1] -9680.396 -9446.449

cmfx1v1agnsc3agff <- coxph(Surv(survyrs1999, survstat1999) ~ ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                           data = wmadf[ccfx1c2ag4, ])
summary(cmfx1v1agnsc3agff)
cmfx1v1agnsc3agff$loglik 

anova(cmfx1v1agnsc3agff, cmfx1b60c3agff)  ##    P(>|Chi|) = 3.859e-05   Binary highly predictive




#####  bcl2  ######

## wmadf$h3k27me3_v1.1n -->> wmadf$bcl2_v1.ppnc
## h3k27me3_v1.1.pp -->> bcl2_v1.ppnc
## H3K27me3 -->> bcl2
## h3v11 -->> bc2v1

wmadf$bcl2_c_v2.ppnc <- as.numeric(wmadf$bcl2_c_v2.pp)

table(wmadf$bcl2_c_v2.ppnc, useNA = "always")


## Original bcl2 analysis done on the odd classification by image colour in bcl2_v1
## Get some info on the percent positive version.

splitpctvec <- (1:10)*10
splitpct <- 90

wmadf$bcl2_c_v2.ppnc_b10c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b20c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b30c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b40c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b50c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b60c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b70c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b80c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b90c <- rep(NA_character_, nrow(wmadf))
wmadf$bcl2_c_v2.ppnc_b100c <- rep(NA_character_, nrow(wmadf))

wmadf$bcl2_c_v2.ppnc_b10c[wmadf$bcl2_c_v2.ppnc <  10] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b10c[wmadf$bcl2_c_v2.ppnc >= 10] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b10f <- factor(wmadf$bcl2_c_v2.ppnc_b10c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b20c[wmadf$bcl2_c_v2.ppnc <  20] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b20c[wmadf$bcl2_c_v2.ppnc >= 20] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b20f <- factor(wmadf$bcl2_c_v2.ppnc_b20c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b30c[wmadf$bcl2_c_v2.ppnc <  30] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b30c[wmadf$bcl2_c_v2.ppnc >= 30] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b30f <- factor(wmadf$bcl2_c_v2.ppnc_b30c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b40c[wmadf$bcl2_c_v2.ppnc <  40] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b40c[wmadf$bcl2_c_v2.ppnc >= 40] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b40f <- factor(wmadf$bcl2_c_v2.ppnc_b40c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b50c[wmadf$bcl2_c_v2.ppnc <  50] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b50c[wmadf$bcl2_c_v2.ppnc >= 50] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b50f <- factor(wmadf$bcl2_c_v2.ppnc_b50c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b60c[wmadf$bcl2_c_v2.ppnc <  60] <- "BCL2 low"
wmadf$bcl2_c_v2.ppnc_b60c[wmadf$bcl2_c_v2.ppnc >= 60] <- "BCL2 high"
wmadf$bcl2_c_v2.ppnc_b60f <- factor(wmadf$bcl2_c_v2.ppnc_b60c, levels = c("BCL2 low", "BCL2 high"))
wmadf$bcl2_c_v2.ppnc_b70c[wmadf$bcl2_c_v2.ppnc <  70] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b70c[wmadf$bcl2_c_v2.ppnc >= 70] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b70f <- factor(wmadf$bcl2_c_v2.ppnc_b70c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b80c[wmadf$bcl2_c_v2.ppnc <  80] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b80c[wmadf$bcl2_c_v2.ppnc >= 80] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b80f <- factor(wmadf$bcl2_c_v2.ppnc_b80c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b90c[wmadf$bcl2_c_v2.ppnc <  90] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b90c[wmadf$bcl2_c_v2.ppnc >= 90] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b90f <- factor(wmadf$bcl2_c_v2.ppnc_b90c, levels = c("bcl2 low", "bcl2 high"))
wmadf$bcl2_c_v2.ppnc_b100c[wmadf$bcl2_c_v2.ppnc <  100] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_b100c[wmadf$bcl2_c_v2.ppnc >= 100] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_b100f <- factor(wmadf$bcl2_c_v2.ppnc_b100c, levels = c("bcl2 low", "bcl2 high"))

table(wmadf$bcl2_c_v2.ppnc_b10f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b20f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b30f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b40f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b50f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b60f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b70f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b80f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b90f , useNA = "always")
table(wmadf$bcl2_c_v2.ppnc_b100f , useNA = "always")


bc2v2b10.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b10f , data = wmadf)
bc2v2b20.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b20f , data = wmadf)
bc2v2b30.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b30f , data = wmadf)
bc2v2b40.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b40f , data = wmadf)
bc2v2b50.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b50f , data = wmadf)
bc2v2b60.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf)
bc2v2b70.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b70f , data = wmadf)
bc2v2b80.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b80f , data = wmadf)
bc2v2b90.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b90f , data = wmadf)
bc2v2b100.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b100f , data = wmadf)

{par(mfrow=c(3, 4))
  plot(bc2v2b10.surv, lty = 1:2, main = "bcl2 < vs >= 10")
  title(main = "bcl2 low to high", line=2.2)
  plot(bc2v2b20.surv, lty = 1:2, main = "20")
  plot(bc2v2b30.surv, lty = 1:2, main = "30")
  plot(bc2v2b40.surv, lty = 1:2, main = "40")
  plot(bc2v2b50.surv, lty = 1:2, main = "50")
  plot(bc2v2b60.surv, lty = 1:2, main = "60")
  plot(bc2v2b70.surv, lty = 1:2, main = "70")
  plot(bc2v2b80.surv, lty = 1:2, main = "80")
  plot(bc2v2b90.surv, lty = 1:2, main = "90")
  plot(bc2v2b100.surv, lty = 1:2, main = "100")
  
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b100f))  , lty = 1:2)
}

## Cox model fits
table(wmadf$ageg4, wmadf$bcl2_c_v2.ppnc_b60f)

ccbc2c2ag4 <- with(wmadf, complete.cases(bcl2_c_v2.ppnc_b60f, ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
ccbc2c2ag4 <- with(wmadf, complete.cases(bcl2_c_v2.ppnc_b90f, ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(ccbc2c2ag4)
head(ccbc2c2ag4)

cmbc2b10c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b10f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b10c3agff)
cmbc2b10c3agff$loglik  ##  [1] -10096.683  -9846.529 

cmbc2b20c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b20f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b20c3agff)
cmbc2b20c3agff$loglik  ##   [1] -10096.683  -9840.979

cmbc2b30c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b30f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b30c3agff)
cmbc2b30c3agff$loglik  ##   [1] -10096.683  -9839.002

cmbc2b40c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b40f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b40c3agff)
cmbc2b40c3agff$loglik  ##   [1] -10096.683  -9839.194

cmbc2b50c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b50f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b50c3agff)
cmbc2b50c3agff$loglik  ##   [1] -10096.683  -9838.537

cmbc2b60c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b60c3agff)
cmbc2b60c3agff$loglik  ##  [1] -10096.683  -9839.724

cmbc2b70c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b70f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b70c3agff)
cmbc2b70c3agff$loglik  ##  [1] -10096.683  -9843.998

cmbc2b80c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b80f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b80c3agff)
cmbc2b80c3agff$loglik  ##  [1] -10096.683  -9843.683

cmbc2b90c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b90f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                        data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b90c3agff)
cmbc2b90c3agff$loglik ## [1] -10096.683  -9841.684

cmbc2b100c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b100f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                         data = wmadf[ccbc2c2ag4, ])
summary(cmbc2b100c3agff)
cmbc2b100c3agff$loglik  ## [1] -10096.683  -9848.654

## 60 pct cells positively stained is binary optimum. <60 vs >=60  bcl2_c_v2.ppnc_b60f 

##

wmadf$bcl2_c_v2.ppnc_t60100c <- rep(NA_character_, nrow(wmadf))

wmadf$bcl2_c_v2.ppnc_t60100c[wmadf$bcl2_c_v2.ppnc <  60] <- "bcl2 low"
wmadf$bcl2_c_v2.ppnc_t60100c[wmadf$bcl2_c_v2.ppnc >= 60 & wmadf$bcl2_c_v2.ppnc < 100] <- "bcl2 mid"
wmadf$bcl2_c_v2.ppnc_t60100c[wmadf$bcl2_c_v2.ppnc >= 100] <- "bcl2 high"
wmadf$bcl2_c_v2.ppnc_t60100f <- factor(wmadf$bcl2_c_v2.ppnc_t60100c, levels = c("bcl2 low", "bcl2 mid", "bcl2 high"))
table(wmadf$bcl2_c_v2.ppnc_t60100f)

bc2v2t60100.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_t60100f , data = wmadf)
plot(bc2v2t60100.surv, lty = 1:3, main = "[0,60),[60-100),[100,100]")

cmbc2v2t60100c3agff <- coxph(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_t60100f * ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                             data = wmadf[ccbc2c2ag4, ])
summary(cmbc2v2t60100c3agff)
cmbc2v2t60100c3agff$loglik  ##  [1] -9680.396 -9446.449

cmbc2v2agnsc3agff <- coxph(Surv(survyrs1999, survstat1999) ~ ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                           data = wmadf[ccbc2c2ag4, ])
summary(cmbc2v2agnsc3agff)
cmbc2v2agnsc3agff$loglik 

anova(cmbc2v2agnsc3agff, cmbc2b60c3agff)  ##    P(>|Chi|) = 3.859e-05   Binary highly predictive

######################################

table(wmadf$bcl2_c_v2.ppnc_b60f)
table(wmadf$foxa1_v1.ppnc_b60f)
table(wmadf$ageg4)


# bcl2_c_v2.ppnc_b60f foxa1_v1.ppnc_b60f  ageg4  
# Basic survival analysis and graphs.


### Split by xxxx

table(wmadf$bcl2_c_v2.ppnc_b60f)
table(wmadf$foxa1_v1.ppnc_b60f)
table(wmadf$ageg4)



bc2c2ag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(0,40]"  ) ) )
bc2c2ag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(40,55]"  ) ) )
bc2c2ag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(55,70]"  ) ) )
bc2c2ag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(70,99]"  ) ) )

#logrank and pvals
bc2c2ag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(0,40]"  ) ) )
bc2c2ag1.lrpv <- lrpv(bc2c2ag1.sdif)
bc2c2ag1.lrpvtxt <- paste("p =", format(bc2c2ag1.lrpv, digits = 3))
bc2c2ag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(40,55]"  ) ) )
bc2c2ag2.lrpv <- lrpv(bc2c2ag2.sdif)
bc2c2ag2.lrpvtxt <- paste("p =", format(bc2c2ag2.lrpv, digits = 3))
bc2c2ag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(55,70]"  ) ) )
bc2c2ag3.lrpv <- lrpv(bc2c2ag3.sdif)
bc2c2ag3.lrpvtxt <- paste("p =", format(bc2c2ag3.lrpv, digits = 3))
bc2c2ag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(70,99]"  ) ) )
bc2c2ag4.lrpv <- lrpv(bc2c2ag4.sdif)
bc2c2ag4.lrpvtxt <- paste("p =", format(bc2c2ag4.lrpv, digits = 3))


## Stratify by FOXA1
bc2c2fnag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
#logrank and pvals
bc2c2fnag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag1.lrpv <- lrpv(bc2c2fnag1.sdif)
bc2c2fnag1.lrpvtxt <- paste("p =", format(bc2c2fnag1.lrpv, digits = 3))
bc2c2fnag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag2.lrpv <- lrpv(bc2c2fnag2.sdif)
bc2c2fnag2.lrpvtxt <- paste("p =", format(bc2c2fnag2.lrpv, digits = 3))
bc2c2fnag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag3.lrpv <- lrpv(bc2c2fnag3.sdif)
bc2c2fnag3.lrpvtxt <- paste("p =", format(bc2c2fnag3.lrpv, digits = 3))
bc2c2fnag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag4.lrpv <- lrpv(bc2c2fnag4.sdif)
bc2c2fnag4.lrpvtxt <- paste("p =", format(bc2c2fnag4.lrpv, digits = 3))
bc2c2fpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag1.lrpv <- lrpv(bc2c2fpag1.sdif)
bc2c2fpag1.lrpvtxt <- paste("p =", format(bc2c2fpag1.lrpv, digits = 3))
bc2c2fpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag2.lrpv <- lrpv(bc2c2fpag2.sdif)
bc2c2fpag2.lrpvtxt <- paste("p =", format(bc2c2fpag2.lrpv, digits = 3))
bc2c2fpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag3.lrpv <- lrpv(bc2c2fpag3.sdif)
bc2c2fpag3.lrpvtxt <- paste("p =", format(bc2c2fpag3.lrpv, digits = 3))
bc2c2fpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag4.lrpv <- lrpv(bc2c2fpag4.sdif)
bc2c2fpag4.lrpvtxt <- paste("p =", format(bc2c2fpag4.lrpv, digits = 3))


fx1c2ag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(0,40]"  ) ) )
fx1c2ag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(40,55]"  ) ) )
fx1c2ag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(55,70]"  ) ) )
fx1c2ag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                         subset = ((ageg4 == "(70,99]"  ) ) )

#logrank and pvals
fx1c2ag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(0,40]"  ) ) )
fx1c2ag1.lrpv <- lrpv(fx1c2ag1.sdif)
fx1c2ag1.lrpvtxt <- paste("p =", format(fx1c2ag1.lrpv, digits = 3))
fx1c2ag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(40,55]"  ) ) )
fx1c2ag2.lrpv <- lrpv(fx1c2ag2.sdif)
fx1c2ag2.lrpvtxt <- paste("p =", format(fx1c2ag2.lrpv, digits = 3))
fx1c2ag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(55,70]"  ) ) )
fx1c2ag3.lrpv <- lrpv(fx1c2ag3.sdif)
fx1c2ag3.lrpvtxt <- paste("p =", format(fx1c2ag3.lrpv, digits = 3))
fx1c2ag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                          subset = ((ageg4 == "(70,99]"  ) ) )
fx1c2ag4.lrpv <- lrpv(fx1c2ag4.sdif)
fx1c2ag4.lrpvtxt <- paste("p =", format(fx1c2ag4.lrpv, digits = 3))


## Stratify by BCL2
fx1c2bnag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
#logrank and pvals
fx1c2bnag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag1.lrpv <- lrpv(fx1c2bnag1.sdif)
fx1c2bnag1.lrpvtxt <- paste("p =", format(fx1c2bnag1.lrpv, digits = 3))
fx1c2bnag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag2.lrpv <- lrpv(fx1c2bnag2.sdif)
fx1c2bnag2.lrpvtxt <- paste("p =", format(fx1c2bnag2.lrpv, digits = 3))
fx1c2bnag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag3.lrpv <- lrpv(fx1c2bnag3.sdif)
fx1c2bnag3.lrpvtxt <- paste("p =", format(fx1c2bnag3.lrpv, digits = 3))
fx1c2bnag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag4.lrpv <- lrpv(fx1c2bnag4.sdif)
fx1c2bnag4.lrpvtxt <- paste("p =", format(fx1c2bnag4.lrpv, digits = 3))
fx1c2bpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag1.lrpv <- lrpv(fx1c2bpag1.sdif)
fx1c2bpag1.lrpvtxt <- paste("p =", format(fx1c2bpag1.lrpv, digits = 3))
fx1c2bpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag2.lrpv <- lrpv(fx1c2bpag2.sdif)
fx1c2bpag2.lrpvtxt <- paste("p =", format(fx1c2bpag2.lrpv, digits = 3))
fx1c2bpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag3.lrpv <- lrpv(fx1c2bpag3.sdif)
fx1c2bpag3.lrpvtxt <- paste("p =", format(fx1c2bpag3.lrpv, digits = 3))
fx1c2bpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag4.lrpv <- lrpv(fx1c2bpag4.sdif)
fx1c2bpag4.lrpvtxt <- paste("p =", format(fx1c2bpag4.lrpv, digits = 3))

{par(mfcol=c(4, 4))
  plot(bc2c2fnag1.surv, lty = 1:3, main = "(0,40] FOXA1 neg")
  text(x=0.5, y=0.1, labels = bc2c2fnag1.lrpvtxt, adj=0)
  plot(bc2c2fnag2.surv, lty = 1:3, main = "(40,55] FOXA1 neg")
  text(x=0.5, y=0.1, labels = bc2c2fnag2.lrpvtxt, adj=0)
  plot(bc2c2fnag3.surv, lty = 1:3, main = "(55,70] FOXA1 neg")
  text(x=0.5, y=0.1, labels = bc2c2fnag3.lrpvtxt, adj=0)
  plot(bc2c2fnag4.surv, lty = 1:3, main = "(70,99] FOXA1 neg")
  text(x=0.5, y=0.1, labels = bc2c2fnag4.lrpvtxt, adj=0)
  plot(bc2c2fpag1.surv, lty = 1:3, main = "(0,40] FOXA1 pos")
  text(x=0.5, y=0.1, labels = bc2c2fpag1.lrpvtxt, adj=0)
  plot(bc2c2fpag2.surv, lty = 1:3, main = "(40,55] FOXA1 pos")
  text(x=0.5, y=0.1, labels = bc2c2fpag2.lrpvtxt, adj=0)
  plot(bc2c2fpag3.surv, lty = 1:3, main = "(55,70] FOXA1 pos")
  text(x=0.5, y=0.1, labels = bc2c2fpag3.lrpvtxt, adj=0)
  plot(bc2c2fpag4.surv, lty = 1:3, main = "(70,99] FOXA1 pos")
  text(x=0.5, y=0.1, labels = bc2c2fpag4.lrpvtxt, adj=0)
  
  legend("topright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2)
  title(main = "BCL2 low to high", line=3)
  
  plot(fx1c2bnag1.surv, lty = 1:3, main = "(0,40] BCL2 low")
  text(x=0.5, y=0.1, labels = fx1c2bnag1.lrpvtxt, adj=0)
  plot(fx1c2bnag2.surv, lty = 1:3, main = "(40,55] BCL2 low")
  text(x=0.5, y=0.1, labels = fx1c2bnag2.lrpvtxt, adj=0)
  plot(fx1c2bnag3.surv, lty = 1:3, main = "(55,70] BCL2 low")
  text(x=0.5, y=0.1, labels = fx1c2bnag3.lrpvtxt, adj=0)
  plot(fx1c2bnag4.surv, lty = 1:3, main = "(70,99] BCL2 low")
  text(x=0.5, y=0.1, labels = fx1c2bnag4.lrpvtxt, adj=0)
  plot(fx1c2bpag1.surv, lty = 1:3, main = "(0,40] BCL2 high")
  text(x=0.5, y=0.1, labels = fx1c2bpag1.lrpvtxt, adj=0)
  plot(fx1c2bpag2.surv, lty = 1:3, main = "(40,55] BCL2 high")
  text(x=0.5, y=0.1, labels = fx1c2bpag2.lrpvtxt, adj=0)
  plot(fx1c2bpag3.surv, lty = 1:3, main = "(55,70] BCL2 high")
  text(x=0.5, y=0.1, labels = fx1c2bpag3.lrpvtxt, adj=0)
  plot(fx1c2bpag4.surv, lty = 1:3, main = "(70,99] BCL2 high")
  text(x=0.5, y=0.1, labels = fx1c2bpag4.lrpvtxt, adj=0)
  
  legend("topright",legend = levels(as.factor(wmadf$foxa1_v1.ppnc_b60f))  , lty = 1:2)
  title(main = "FOXA1 low to high", line=3)
}
## Cox


ccbc2c2fzag <- with(wmadf, complete.cases(bcl2_c_v2.ppnc_b60f, foxa1_v1.ppnc_b60f, 
                                          ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(ccbc2c2fzag)
head(ccbc2c2fzag)

cmbc2c2fzagff <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f * bcl2_c_v2.ppnc_b60f * ageg4 + 
                         grade_b12v3 + nodestat + size_tumor_grp, 
                       data = wmadf[ccbc2c2fzag, ])
cmbc2c2fzagff


### Split by xxxx


wmadf$ageg4 <- cut(wmadf$age, breaks = c(0, 40, 55, 70, 99))
table(wmadf$ageg4)


bc2c2fnag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
#logrank and pvals
bc2c2fnag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag1.lrpv <- lrpv(bc2c2fnag1.sdif)
bc2c2fnag1.lrpvtxt <- paste("p =", format(bc2c2fnag1.lrpv, digits = 3))
bc2c2fnag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag2.lrpv <- lrpv(bc2c2fnag2.sdif)
bc2c2fnag2.lrpvtxt <- paste("p =", format(bc2c2fnag2.lrpv, digits = 3))
bc2c2fnag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag3.lrpv <- lrpv(bc2c2fnag3.sdif)
bc2c2fnag3.lrpvtxt <- paste("p =", format(bc2c2fnag3.lrpv, digits = 3))
bc2c2fnag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low")) )
bc2c2fnag4.lrpv <- lrpv(bc2c2fnag4.sdif)
bc2c2fnag4.lrpvtxt <- paste("p =", format(bc2c2fnag4.lrpv, digits = 3))
bc2c2fpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag1.lrpv <- lrpv(bc2c2fpag1.sdif)
bc2c2fpag1.lrpvtxt <- paste("p =", format(bc2c2fpag1.lrpv, digits = 3))
bc2c2fpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag2.lrpv <- lrpv(bc2c2fpag2.sdif)
bc2c2fpag2.lrpvtxt <- paste("p =", format(bc2c2fpag2.lrpv, digits = 3))
bc2c2fpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag3.lrpv <- lrpv(bc2c2fpag3.sdif)
bc2c2fpag3.lrpvtxt <- paste("p =", format(bc2c2fpag3.lrpv, digits = 3))
bc2c2fpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high")) )
bc2c2fpag4.lrpv <- lrpv(bc2c2fpag4.sdif)
bc2c2fpag4.lrpvtxt <- paste("p =", format(bc2c2fpag4.lrpv, digits = 3))


fx1c2bnag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
#logrank and pvals
fx1c2bnag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag1.lrpv <- lrpv(fx1c2bnag1.sdif)
fx1c2bnag1.lrpvtxt <- paste("p =", format(fx1c2bnag1.lrpv, digits = 3))
fx1c2bnag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag2.lrpv <- lrpv(fx1c2bnag2.sdif)
fx1c2bnag2.lrpvtxt <- paste("p =", format(fx1c2bnag2.lrpv, digits = 3))
fx1c2bnag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag3.lrpv <- lrpv(fx1c2bnag3.sdif)
fx1c2bnag3.lrpvtxt <- paste("p =", format(fx1c2bnag3.lrpv, digits = 3))
fx1c2bnag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low")) )
fx1c2bnag4.lrpv <- lrpv(fx1c2bnag4.sdif)
fx1c2bnag4.lrpvtxt <- paste("p =", format(fx1c2bnag4.lrpv, digits = 3))
fx1c2bpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag1.lrpv <- lrpv(fx1c2bpag1.sdif)
fx1c2bpag1.lrpvtxt <- paste("p =", format(fx1c2bpag1.lrpv, digits = 3))
fx1c2bpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag2.lrpv <- lrpv(fx1c2bpag2.sdif)
fx1c2bpag2.lrpvtxt <- paste("p =", format(fx1c2bpag2.lrpv, digits = 3))
fx1c2bpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag3.lrpv <- lrpv(fx1c2bpag3.sdif)
fx1c2bpag3.lrpvtxt <- paste("p =", format(fx1c2bpag3.lrpv, digits = 3))
fx1c2bpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high")) )
fx1c2bpag4.lrpv <- lrpv(fx1c2bpag4.sdif)
fx1c2bpag4.lrpvtxt <- paste("p =", format(fx1c2bpag4.lrpv, digits = 3))


pdf(file = paste(plot_dirs, "FOXA1_BCL2_ER_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")

pdf(file = paste(plot_dirs, "FOXA1_BCL2_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "FOXA1_BCL2_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  plot(bc2c2fnag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40 FOXA1 neg", paste(": N =", paste(bc2c2fnag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] FOXA1 neg", paste(": N =", paste(bc2c2fnag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnag2.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] FOXA1 neg", paste(": N =", paste(bc2c2fnag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] FOXA1 neg", paste(": N =", paste(bc2c2fnag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnag4.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] FOXA1 pos", paste(": N =", paste(bc2c2fpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] FOXA1 pos", paste(": N =", paste(bc2c2fpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = bc2c2fpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.6, col = Dark2[c(1, 2)])
  
  plot(bc2c2fpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] FOXA1 pos", paste(": N =", paste(bc2c2fpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] FOXA1 pos", paste(": N =", paste(bc2c2fpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fpag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(fx1c2bnag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] BCL2 low", paste(": N =", paste(fx1c2bnag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bnag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] BCL2 low", paste(": N =", paste(fx1c2bnag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2bnag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.6, col = Dark2[c(3, 4)])
  
  plot(fx1c2bnag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] BCL2 low", paste(": N =", paste(fx1c2bnag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bnag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] BCL2 low", paste(": N =", paste(fx1c2bnag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnag4.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] BCL2 high", paste(": N =", paste(fx1c2bpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] BCL2 high", paste(": N =", paste(fx1c2bpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpag2.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] BCL2 high", paste(": N =", paste(fx1c2bpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] BCL2 high", paste(": N =", paste(fx1c2bpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "Strata: Age, FOXA1", line=1, outer = TRUE, adj = 0.2, cex = 1)
  mtext(text = "Strata: Age, BCL2", line=1, outer = TRUE, adj = 0.8, cex = 1)
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()

## ER+  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

bc2c2fnerpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fnerpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fnerpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fnerpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
#logrank and pvals
bc2c2fnerpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fnerpag1.lrpv <- lrpv(bc2c2fnerpag1.sdif)
bc2c2fnerpag1.lrpvtxt <- paste("p =", format(bc2c2fnerpag1.lrpv, digits = 3))
bc2c2fnerpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fnerpag2.lrpv <- lrpv(bc2c2fnerpag2.sdif)
bc2c2fnerpag2.lrpvtxt <- paste("p =", format(bc2c2fnerpag2.lrpv, digits = 3))
bc2c2fnerpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fnerpag3.lrpv <- lrpv(bc2c2fnerpag3.sdif)
bc2c2fnerpag3.lrpvtxt <- paste("p =", format(bc2c2fnerpag3.lrpv, digits = 3))
bc2c2fnerpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fnerpag4.lrpv <- lrpv(bc2c2fnerpag4.sdif)
bc2c2fnerpag4.lrpvtxt <- paste("p =", format(bc2c2fnerpag4.lrpv, digits = 3))
bc2c2fperpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag1.lrpv <- lrpv(bc2c2fperpag1.sdif)
bc2c2fperpag1.lrpvtxt <- paste("p =", format(bc2c2fperpag1.lrpv, digits = 3))
bc2c2fperpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag2.lrpv <- lrpv(bc2c2fperpag2.sdif)
bc2c2fperpag2.lrpvtxt <- paste("p =", format(bc2c2fperpag2.lrpv, digits = 3))
bc2c2fperpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag3.lrpv <- lrpv(bc2c2fperpag3.sdif)
bc2c2fperpag3.lrpvtxt <- paste("p =", format(bc2c2fperpag3.lrpv, digits = 3))
bc2c2fperpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2fperpag4.lrpv <- lrpv(bc2c2fperpag4.sdif)
bc2c2fperpag4.lrpvtxt <- paste("p =", format(bc2c2fperpag4.lrpv, digits = 3))


fx1c2bnerpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bnerpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bnerpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bnerpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
#logrank and pvals
fx1c2bnerpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bnerpag1.lrpv <- lrpv(fx1c2bnerpag1.sdif)
fx1c2bnerpag1.lrpvtxt <- paste("p =", format(fx1c2bnerpag1.lrpv, digits = 3))
fx1c2bnerpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bnerpag2.lrpv <- lrpv(fx1c2bnerpag2.sdif)
fx1c2bnerpag2.lrpvtxt <- paste("p =", format(fx1c2bnerpag2.lrpv, digits = 3))
fx1c2bnerpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bnerpag3.lrpv <- lrpv(fx1c2bnerpag3.sdif)
fx1c2bnerpag3.lrpvtxt <- paste("p =", format(fx1c2bnerpag3.lrpv, digits = 3))
fx1c2bnerpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bnerpag4.lrpv <- lrpv(fx1c2bnerpag4.sdif)
fx1c2bnerpag4.lrpvtxt <- paste("p =", format(fx1c2bnerpag4.lrpv, digits = 3))
fx1c2bperpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag1.lrpv <- lrpv(fx1c2bperpag1.sdif)
fx1c2bperpag1.lrpvtxt <- paste("p =", format(fx1c2bperpag1.lrpv, digits = 3))
fx1c2bperpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag2.lrpv <- lrpv(fx1c2bperpag2.sdif)
fx1c2bperpag2.lrpvtxt <- paste("p =", format(fx1c2bperpag2.lrpv, digits = 3))
fx1c2bperpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag3.lrpv <- lrpv(fx1c2bperpag3.sdif)
fx1c2bperpag3.lrpvtxt <- paste("p =", format(fx1c2bperpag3.lrpv, digits = 3))
fx1c2bperpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2bperpag4.lrpv <- lrpv(fx1c2bperpag4.sdif)
fx1c2bperpag4.lrpvtxt <- paste("p =", format(fx1c2bperpag4.lrpv, digits = 3))


pdf(file = paste(plot_dirs, "FOXA1_BCL2_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "FOXA1_BCL2_ERpos_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  plot(bc2c2fnerpag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40 FOXA1 neg", paste(": N =", paste(bc2c2fnerpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnerpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnerpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] FOXA1 neg", paste(": N =", paste(bc2c2fnerpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnerpag2.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnerpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] FOXA1 neg", paste(": N =", paste(bc2c2fnerpag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnerpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnerpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] FOXA1 neg", paste(": N =", paste(bc2c2fnerpag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnerpag4.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fperpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] FOXA1 pos", paste(": N =", paste(bc2c2fperpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fperpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fperpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] FOXA1 pos", paste(": N =", paste(bc2c2fperpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = bc2c2fperpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, col = Dark2[c(1, 2)], cex = 0.6)
  
  plot(bc2c2fperpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] FOXA1 pos", paste(": N =", paste(bc2c2fperpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fperpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fperpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] FOXA1 pos", paste(": N =", paste(bc2c2fperpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fperpag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(fx1c2bnerpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] BCL2 low", paste(": N =", paste(fx1c2bnerpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnerpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bnerpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] BCL2 low", paste(": N =", paste(fx1c2bnerpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2bnerpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, col = Dark2[c(3, 4)], cex = 0.6)
  
  plot(fx1c2bnerpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] BCL2 low", paste(": N =", paste(fx1c2bnerpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnerpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bnerpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] BCL2 low", paste(": N =", paste(fx1c2bnerpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnerpag4.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bperpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] BCL2 high", paste(": N =", paste(fx1c2bperpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bperpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bperpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] BCL2 high", paste(": N =", paste(fx1c2bperpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bperpag2.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bperpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] BCL2 high", paste(": N =", paste(fx1c2bperpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bperpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bperpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] BCL2 high", paste(": N =", paste(fx1c2bperpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bperpag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "Strata ER+ : Age, FOXA1", line=1, outer = TRUE, adj = 0.2, cex = 1)
  mtext(text = "Strata ER+ : Age, BCL2", line=1, outer = TRUE, adj = 0.8, cex = 1)
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()


## ER-  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

bc2c2fnernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fnernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fnernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fnernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
#logrank and pvals
bc2c2fnernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fnernag1.lrpv <- lrpv(bc2c2fnernag1.sdif)
bc2c2fnernag1.lrpvtxt <- paste("p =", format(bc2c2fnernag1.lrpv, digits = 3))
bc2c2fnernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fnernag2.lrpv <- lrpv(bc2c2fnernag2.sdif)
bc2c2fnernag2.lrpvtxt <- paste("p =", format(bc2c2fnernag2.lrpv, digits = 3))
bc2c2fnernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fnernag3.lrpv <- lrpv(bc2c2fnernag3.sdif)
bc2c2fnernag3.lrpvtxt <- paste("p =", format(bc2c2fnernag3.lrpv, digits = 3))
bc2c2fnernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 low") & (er_b0v123_v2 == "negative")) )
bc2c2fnernag4.lrpv <- lrpv(bc2c2fnernag4.sdif)
bc2c2fnernag4.lrpvtxt <- paste("p =", format(bc2c2fnernag4.lrpv, digits = 3))
bc2c2fpernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag1.lrpv <- lrpv(bc2c2fpernag1.sdif)
bc2c2fpernag1.lrpvtxt <- paste("p =", format(bc2c2fpernag1.lrpv, digits = 3))
bc2c2fpernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag2.lrpv <- lrpv(bc2c2fpernag2.sdif)
bc2c2fpernag2.lrpvtxt <- paste("p =", format(bc2c2fpernag2.lrpv, digits = 3))
bc2c2fpernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag3.lrpv <- lrpv(bc2c2fpernag3.sdif)
bc2c2fpernag3.lrpvtxt <- paste("p =", format(bc2c2fpernag3.lrpv, digits = 3))
bc2c2fpernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (foxa1_v1.ppnc_b60f == "FOXA1 high") & (er_b0v123_v2 == "negative")) )
bc2c2fpernag4.lrpv <- lrpv(bc2c2fpernag4.sdif)
bc2c2fpernag4.lrpvtxt <- paste("p =", format(bc2c2fpernag4.lrpv, digits = 3))


fx1c2bnernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bnernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bnernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bnernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                              subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
#logrank and pvals
fx1c2bnernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bnernag1.lrpv <- lrpv(fx1c2bnernag1.sdif)
fx1c2bnernag1.lrpvtxt <- paste("p =", format(fx1c2bnernag1.lrpv, digits = 3))
fx1c2bnernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bnernag2.lrpv <- lrpv(fx1c2bnernag2.sdif)
fx1c2bnernag2.lrpvtxt <- paste("p =", format(fx1c2bnernag2.lrpv, digits = 3))
fx1c2bnernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bnernag3.lrpv <- lrpv(fx1c2bnernag3.sdif)
fx1c2bnernag3.lrpvtxt <- paste("p =", format(fx1c2bnernag3.lrpv, digits = 3))
fx1c2bnernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 low") & (er_b0v123_v2 == "negative")) )
fx1c2bnernag4.lrpv <- lrpv(fx1c2bnernag4.sdif)
fx1c2bnernag4.lrpvtxt <- paste("p =", format(fx1c2bnernag4.lrpv, digits = 3))
fx1c2bpernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(0,40]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag1.lrpv <- lrpv(fx1c2bpernag1.sdif)
fx1c2bpernag1.lrpvtxt <- paste("p =", format(fx1c2bpernag1.lrpv, digits = 3))
fx1c2bpernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(40,55]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag2.lrpv <- lrpv(fx1c2bpernag2.sdif)
fx1c2bpernag2.lrpvtxt <- paste("p =", format(fx1c2bpernag2.lrpv, digits = 3))
fx1c2bpernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(55,70]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag3.lrpv <- lrpv(fx1c2bpernag3.sdif)
fx1c2bpernag3.lrpvtxt <- paste("p =", format(fx1c2bpernag3.lrpv, digits = 3))
fx1c2bpernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                               subset = ((ageg4 == "(70,99]"  ) & (bcl2_c_v2.ppnc_b60f == "BCL2 high") & (er_b0v123_v2 == "negative")) )
fx1c2bpernag4.lrpv <- lrpv(fx1c2bpernag4.sdif)
fx1c2bpernag4.lrpvtxt <- paste("p =", format(fx1c2bpernag4.lrpv, digits = 3))


pdf(file = paste(plot_dirs, "FOXA1_BCL2_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "FOXA1_BCL2_ERneg_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)


{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  plot(bc2c2fnernag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40 FOXA1 neg", paste(": N =", paste(bc2c2fnernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] FOXA1 neg", paste(": N =", paste(bc2c2fnernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnernag2.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] FOXA1 neg", paste(": N =", paste(bc2c2fnernag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fnernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] FOXA1 neg", paste(": N =", paste(bc2c2fnernag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fnernag4.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fpernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] FOXA1 pos", paste(": N =", paste(bc2c2fpernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fpernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fpernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] FOXA1 pos", paste(": N =", paste(bc2c2fpernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = bc2c2fpernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, col = Dark2[c(1, 2)], cex = 0.6)
  
  plot(bc2c2fpernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] FOXA1 pos", paste(": N =", paste(bc2c2fpernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fpernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2fpernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] FOXA1 pos", paste(": N =", paste(bc2c2fpernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = bc2c2fpernag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(fx1c2bnernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] BCL2 low", paste(": N =", paste(fx1c2bnernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bnernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] BCL2 low", paste(": N =", paste(fx1c2bnernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2bnernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, col = Dark2[c(3, 4)], cex = 0.6)
  
  plot(fx1c2bnernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] BCL2 low", paste(": N =", paste(fx1c2bnernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bnernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] BCL2 low", paste(": N =", paste(fx1c2bnernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bnernag4.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40] BCL2 high", paste(": N =", paste(fx1c2bpernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55] BCL2 high", paste(": N =", paste(fx1c2bpernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpernag2.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70] BCL2 high", paste(": N =", paste(fx1c2bpernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2bpernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99] BCL2 high", paste(": N =", paste(fx1c2bpernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2bpernag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "Strata ER- : Age, FOXA1", line=1, outer = TRUE, adj = 0.2, cex = 1)
  mtext(text = "Strata ER- : Age, BCL2", line=1, outer = TRUE, adj = 0.8, cex = 1)
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()



### Cox models

# foxa1_v1.ppnc_b60f bcl2_c_v2.ppnc_b60f == "BCL2 high" foxa1_v1.ppnc_b60f == "FOXA1 high"
table(wmadf$foxa1_v1.ppnc_b60f)
table(wmadf$bcl2_c_v2.ppnc_b60f)

cccmfbagntp <- with(wmadf, complete.cases(foxa1_v1.ppnc_b60f, bcl2_c_v2.ppnc_b60f, er_b0v123_v2, her2_b012v23_v2, 
                                          ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(cccmfbagntp)

cmfbagntm1 <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f + bcl2_c_v2.ppnc_b60f + 
                      ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmfbagntp, ])
print(cmfbagntm1)
summary(cmfbagntm1)

cmfbagntm2 <- coxph(Surv(survyrs1999, survstat1999) ~ ageg4 + grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmfbagntp, ])
anova(cmfbagntm2, cmfbagntm1)

cmfbagntm3 <- coxph(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f * bcl2_c_v2.ppnc_b60f * ageg4 + 
                      grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmfbagntp, ])
print(cmfbagntm3)
summary(cmfbagntm3)

anova(cmfbagntm3, cmfbagntm1)
anova(cmfbagntm3, cmfbagntm2)

#er_b0v123_v2

cmfbagntm4 <- coxph(Surv(survyrs1999, survstat1999) ~ er_b0v123_v2 * foxa1_v1.ppnc_b60f * bcl2_c_v2.ppnc_b60f * ageg4 + 
                      grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmfbagntp, ])
print(cmfbagntm4)
summary(cmfbagntm4)

anova(cmfbagntm4, cmfbagntm3)
anova(cmfbagntm4, cmfbagntm2)

cmfbagntm5 <- coxph(Surv(survyrs1999, survstat1999) ~ er_b0v123_v2 + foxa1_v1.ppnc_b60f + bcl2_c_v2.ppnc_b60f + ageg4 + 
                      grade_b12v3 + nodestat + size_tumor_grp, 
                    data = wmadf[cccmfbagntp, ])
summary(cmfbagntm5)
anova(cmfbagntm5, cmfbagntm2)


# her2_b012v23_v2 can not additionally be accommodated with interactions

# cmfbagntm5 <- coxph(Surv(survyrs1999, survstat1999) ~ 
#                       er_b0v123_v2 * her2_b012v23_v2 * foxa1_v1.ppnc_b60f * bcl2_c_v2.ppnc_b60f * ageg4 + 
#                       grade_b12v3 + nodestat + size_tumor_grp, 
#                     data = wmadf[cccmfbagntp, ])
# Warning message:
#   In fitter(X, Y, strats, offset, init, control, weights = weights,  :
#               Loglik converged before variable  14,30,33,40,41,42,52,53,54,55,58,62,63,64,65,66,67 ; beta may be infinite. 
#             
# print(cmfbagntm5)
# summary(cmfbagntm5)

# Do additive no interaction model with HER2

cmfbagntm6 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2 + foxa1_v1.ppnc_b60f + bcl2_c_v2.ppnc_b60f + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmfbagntp, ])

summary(cmfbagntm6)
anova(cmfbagntm6, cmfbagntm2)

cmfbagntm7 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2  + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmfbagntp, ])

summary(cmfbagntm7)
anova(cmfbagntm7, cmfbagntm6)


cmfbagntm8 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2 + foxa1_v1.ppnc_b60f  + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmfbagntp, ])

summary(cmfbagntm8)
anova(cmfbagntm8, cmfbagntm7)


cmfbagntm9 <- coxph(Surv(survyrs1999, survstat1999) ~
                      er_b0v123_v2 + her2_b012v23_v2 +  bcl2_c_v2.ppnc_b60f + ageg4 +
                      grade_b12v3 + nodestat + size_tumor_grp,
                    data = wmadf[cccmfbagntp, ])

summary(cmfbagntm9)
anova(cmfbagntm9, cmfbagntm7)

## 4df omnibus test of EZH2, H3K27me3, FOXA1, BCL2


cccmehfbagntp <- 
  with(wmadf, 
       complete.cases(ezh2_blt10vge10_v1.pp, h3k27me3_v1.1.pp_b60f,
                      foxa1_v1.ppnc_b60f, bcl2_c_v2.ppnc_b60f, er_b0v123_v2, her2_b012v23_v2, 
                      ageg4,  grade_b12v3, nodestat, size_tumor_grp, survyrs1999, survstat1999) )
table(cccmehfbagntp)


cmehfbagntm1 <- coxph(Surv(survyrs1999, survstat1999) ~
                        er_b0v123_v2 + her2_b012v23_v2 + 
                        ezh2_blt10vge10_v1.pp + h3k27me3_v1.1.pp_b60f + 
                        foxa1_v1.ppnc_b60f + bcl2_c_v2.ppnc_b60f + 
                        ageg4 + grade_b12v3 + nodestat + size_tumor_grp,
                      data = wmadf[cccmehfbagntp, ])
summary(cmehfbagntm1)

cmehfbagntm2 <- coxph(Surv(survyrs1999, survstat1999) ~
                        er_b0v123_v2 + her2_b012v23_v2 + 
                        ageg4 + grade_b12v3 + nodestat + size_tumor_grp,
                      data = wmadf[cccmehfbagntp, ])

anova(cmehfbagntm2, cmehfbagntm1)

#### EZH2, H3K27me3, FOXA1, BCL2 stratified by age alone and within Age, ER+/-




## ER+  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

h3c2erpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2erpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2erpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2erpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )

#logrank and pvals
h3c2erpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2erpag1.lrpv <- lrpv(h3c2erpag1.sdif)
h3c2erpag1.lrpvtxt <- paste("p =", format(h3c2erpag1.lrpv, digits = 3))
h3c2erpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2erpag2.lrpv <- lrpv(h3c2erpag2.sdif)
h3c2erpag2.lrpvtxt <- paste("p =", format(h3c2erpag2.lrpv, digits = 3))
h3c2erpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2erpag3.lrpv <- lrpv(h3c2erpag3.sdif)
h3c2erpag3.lrpvtxt <- paste("p =", format(h3c2erpag3.lrpv, digits = 3))
h3c2erpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
h3c2erpag4.lrpv <- lrpv(h3c2erpag4.sdif)
h3c2erpag4.lrpvtxt <- paste("p =", format(h3c2erpag4.lrpv, digits = 3))



ezc2erpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2erpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2erpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2erpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )

#logrank and pvals
ezc2erpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2erpag1.lrpv <- lrpv(ezc2erpag1.sdif)
ezc2erpag1.lrpvtxt <- paste("p =", format(ezc2erpag1.lrpv, digits = 3))
ezc2erpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2erpag2.lrpv <- lrpv(ezc2erpag2.sdif)
ezc2erpag2.lrpvtxt <- paste("p =", format(ezc2erpag2.lrpv, digits = 3))
ezc2erpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2erpag3.lrpv <- lrpv(ezc2erpag3.sdif)
ezc2erpag3.lrpvtxt <- paste("p =", format(ezc2erpag3.lrpv, digits = 3))
ezc2erpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
ezc2erpag4.lrpv <- lrpv(ezc2erpag4.sdif)
ezc2erpag4.lrpvtxt <- paste("p =", format(ezc2erpag4.lrpv, digits = 3))




## ER-  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

h3c2ernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
h3c2ernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
h3c2ernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
h3c2ernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )

#logrank and pvals
h3c2ernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
h3c2ernag1.lrpv <- lrpv(h3c2ernag1.sdif)
h3c2ernag1.lrpvtxt <- paste("p =", format(h3c2ernag1.lrpv, digits = 3))
h3c2ernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
h3c2ernag2.lrpv <- lrpv(h3c2ernag2.sdif)
h3c2ernag2.lrpvtxt <- paste("p =", format(h3c2ernag2.lrpv, digits = 3))
h3c2ernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
h3c2ernag3.lrpv <- lrpv(h3c2ernag3.sdif)
h3c2ernag3.lrpvtxt <- paste("p =", format(h3c2ernag3.lrpv, digits = 3))
h3c2ernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ h3k27me3_v1.1.pp_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )
h3c2ernag4.lrpv <- lrpv(h3c2ernag4.sdif)
h3c2ernag4.lrpvtxt <- paste("p =", format(h3c2ernag4.lrpv, digits = 3))



ezc2ernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
ezc2ernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
ezc2ernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
ezc2ernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                           subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )

#logrank and pvals
ezc2ernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
ezc2ernag1.lrpv <- lrpv(ezc2ernag1.sdif)
ezc2ernag1.lrpvtxt <- paste("p =", format(ezc2ernag1.lrpv, digits = 3))
ezc2ernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
ezc2ernag2.lrpv <- lrpv(ezc2ernag2.sdif)
ezc2ernag2.lrpvtxt <- paste("p =", format(ezc2ernag2.lrpv, digits = 3))
ezc2ernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
ezc2ernag3.lrpv <- lrpv(ezc2ernag3.sdif)
ezc2ernag3.lrpvtxt <- paste("p =", format(ezc2ernag3.lrpv, digits = 3))
ezc2ernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ ezh2_blt10vge10_v1.pp , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )
ezc2ernag4.lrpv <- lrpv(ezc2ernag4.sdif)
ezc2ernag4.lrpvtxt <- paste("p =", format(ezc2ernag4.lrpv, digits = 3))




#################################################################


## ER+  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

bc2c2erpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2erpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2erpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2erpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )

#logrank and pvals
bc2c2erpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2erpag1.lrpv <- lrpv(bc2c2erpag1.sdif)
bc2c2erpag1.lrpvtxt <- paste("p =", format(bc2c2erpag1.lrpv, digits = 3))
bc2c2erpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2erpag2.lrpv <- lrpv(bc2c2erpag2.sdif)
bc2c2erpag2.lrpvtxt <- paste("p =", format(bc2c2erpag2.lrpv, digits = 3))
bc2c2erpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2erpag3.lrpv <- lrpv(bc2c2erpag3.sdif)
bc2c2erpag3.lrpvtxt <- paste("p =", format(bc2c2erpag3.lrpv, digits = 3))
bc2c2erpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
bc2c2erpag4.lrpv <- lrpv(bc2c2erpag4.sdif)
bc2c2erpag4.lrpvtxt <- paste("p =", format(bc2c2erpag4.lrpv, digits = 3))


fx1c2erpag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2erpag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2erpag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2erpag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )

#logrank and pvals
fx1c2erpag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2erpag1.lrpv <- lrpv(fx1c2erpag1.sdif)
fx1c2erpag1.lrpvtxt <- paste("p =", format(fx1c2erpag1.lrpv, digits = 3))
fx1c2erpag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2erpag2.lrpv <- lrpv(fx1c2erpag2.sdif)
fx1c2erpag2.lrpvtxt <- paste("p =", format(fx1c2erpag2.lrpv, digits = 3))
fx1c2erpag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2erpag3.lrpv <- lrpv(fx1c2erpag3.sdif)
fx1c2erpag3.lrpvtxt <- paste("p =", format(fx1c2erpag3.lrpv, digits = 3))
fx1c2erpag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "any nuclei staining (score 1-3)")) )
fx1c2erpag4.lrpv <- lrpv(fx1c2erpag4.sdif)
fx1c2erpag4.lrpvtxt <- paste("p =", format(fx1c2erpag4.lrpv, digits = 3))


## ER-  ## ER pos:  er_b0v123_v2 == "any nuclei staining (score 1-3)"   ER neg:   er_b0v123_v2 == "negative"

bc2c2ernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
bc2c2ernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
bc2c2ernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
bc2c2ernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )

#logrank and pvals
bc2c2ernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
bc2c2ernag1.lrpv <- lrpv(bc2c2ernag1.sdif)
bc2c2ernag1.lrpvtxt <- paste("p =", format(bc2c2ernag1.lrpv, digits = 3))
bc2c2ernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
bc2c2ernag2.lrpv <- lrpv(bc2c2ernag2.sdif)
bc2c2ernag2.lrpvtxt <- paste("p =", format(bc2c2ernag2.lrpv, digits = 3))
bc2c2ernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
bc2c2ernag3.lrpv <- lrpv(bc2c2ernag3.sdif)
bc2c2ernag3.lrpvtxt <- paste("p =", format(bc2c2ernag3.lrpv, digits = 3))
bc2c2ernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ bcl2_c_v2.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )
bc2c2ernag4.lrpv <- lrpv(bc2c2ernag4.sdif)
bc2c2ernag4.lrpvtxt <- paste("p =", format(bc2c2ernag4.lrpv, digits = 3))


fx1c2ernag1.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
fx1c2ernag2.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
fx1c2ernag3.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
fx1c2ernag4.surv <- survfit(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                            subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )

#logrank and pvals
fx1c2ernag1.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(0,40]"  ) & (er_b0v123_v2 == "negative")) )
fx1c2ernag1.lrpv <- lrpv(fx1c2ernag1.sdif)
fx1c2ernag1.lrpvtxt <- paste("p =", format(fx1c2ernag1.lrpv, digits = 3))
fx1c2ernag2.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(40,55]"  ) & (er_b0v123_v2 == "negative")) )
fx1c2ernag2.lrpv <- lrpv(fx1c2ernag2.sdif)
fx1c2ernag2.lrpvtxt <- paste("p =", format(fx1c2ernag2.lrpv, digits = 3))
fx1c2ernag3.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(55,70]"  ) & (er_b0v123_v2 == "negative")) )
fx1c2ernag3.lrpv <- lrpv(fx1c2ernag3.sdif)
fx1c2ernag3.lrpvtxt <- paste("p =", format(fx1c2ernag3.lrpv, digits = 3))
fx1c2ernag4.sdif <- survdiff(Surv(survyrs1999, survstat1999) ~ foxa1_v1.ppnc_b60f , data = wmadf, 
                             subset = ((ageg4 == "(70,99]"  ) & (er_b0v123_v2 == "negative")) )
fx1c2ernag4.lrpv <- lrpv(fx1c2ernag4.sdif)
fx1c2ernag4.lrpvtxt <- paste("p =", format(fx1c2ernag4.lrpv, digits = 3))


#### EZH2, H3K27me3, FOXA1, BCL2 stratified by age alone and Age, ER+/-

pdf(file = paste(plot_dirs, "EZH2_H3K27me3_FOXA1_BCL2_ER_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")

pdf(file = paste(plot_dirs, "EZH2_H3K27me3_FOXA1_BCL2_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "EZH2_H3K27me3_FOXA1_BCL2_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  ##   bc2c2fn  ezc2
  plot(ezc2ag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40", paste(": N =", paste(ezc2ag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2ag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(ezc2ag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = ezc2ag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(ezc2ag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(ezc2ag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2ag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(ezc2ag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ag4.lrpvtxt, adj=0, cex = 0.9)
  
  ##   bc2c2fp  h3c2
  plot(h3c2ag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(h3c2ag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2ag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(h3c2ag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2ag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2ag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(h3c2ag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2ag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(h3c2ag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ag4.lrpvtxt, adj=0, cex = 0.9)
  
  ## fx1c2bn  fx1c2
  plot(fx1c2ag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(fx1c2ag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2ag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(fx1c2ag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2ag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(fx1c2ag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(fx1c2ag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2ag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(fx1c2ag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ag4.lrpvtxt, adj=0, cex = 0.9)
  
  ## fx1c2bp  bc2c2
  plot(bc2c2ag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(bc2c2ag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2ag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(bc2c2ag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = bc2c2ag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(bc2c2ag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(bc2c2ag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2ag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(bc2c2ag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "EZH2", line=1, outer = TRUE, adj = 0.12, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2ag1.surv$n, ezc2ag2.surv$n, ezc2ag3.surv$n, ezc2ag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.117, cex = 0.7)
  mtext(text = "H3K27me3", line=1, outer = TRUE, adj = 0.37, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2ag1.surv$n, h3c2ag2.surv$n, h3c2ag3.surv$n, h3c2ag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.377, cex = 0.7)
  mtext(text = "FOXA1", line=1, outer = TRUE, adj = 0.638, cex = 1)
  mtext(text = paste("N =", sum(c(fx1c2ag1.surv$n, fx1c2ag2.surv$n, fx1c2ag3.surv$n, fx1c2ag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.637, cex = 0.7)
  mtext(text = "BCL2", line=1, outer = TRUE, adj = 0.9, cex = 1)
  mtext(text = paste("N =", sum(c(bc2c2ag1.surv$n, bc2c2ag2.surv$n, bc2c2ag3.surv$n, bc2c2ag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.903, cex = 0.7)
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()


#### EZH2, H3K27me3, FOXA1, BCL2 stratified by age and ER+


pdf(file = paste(plot_dirs, "EZH2_H3K27me3_FOXA1_BCL2_ERpos_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "EZH2_H3K27me3_FOXA1_BCL2_ERpos_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  ##   bc2c2fn  ezc2
  plot(ezc2erpag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40", paste(": N =", paste(ezc2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(ezc2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = ezc2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(ezc2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(ezc2erpag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(ezc2erpag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2erpag4.lrpvtxt, adj=0, cex = 0.9)
  
  ##   bc2c2fp  h3c2
  plot(h3c2erpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(h3c2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(h3c2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(h3c2erpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(h3c2erpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2erpag4.lrpvtxt, adj=0, cex = 0.9)
  
  ## fx1c2bn  fx1c2
  plot(fx1c2erpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(fx1c2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(fx1c2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(fx1c2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(fx1c2erpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(fx1c2erpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2erpag4.lrpvtxt, adj=0, cex = 0.9)
  
  ## fx1c2bp  bc2c2
  plot(bc2c2erpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(bc2c2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(bc2c2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = bc2c2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(bc2c2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(bc2c2erpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(bc2c2erpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2erpag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "ER+:  EZH2", line=1, outer = TRUE, adj = 0.095, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2erpag1.surv$n, ezc2erpag2.surv$n, ezc2erpag3.surv$n, ezc2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.117, cex = 0.7)
  mtext(text = "ER+:  H3K27me3", line=1, outer = TRUE, adj = 0.36, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2erpag1.surv$n, h3c2erpag2.surv$n, h3c2erpag3.surv$n, h3c2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.377, cex = 0.7)
  mtext(text = "ER+:  FOXA1", line=1, outer = TRUE, adj = 0.65, cex = 1)
  mtext(text = paste("N =", sum(c(fx1c2erpag1.surv$n, fx1c2erpag2.surv$n, fx1c2erpag3.surv$n, fx1c2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.637, cex = 0.7)
  mtext(text = "ER+:  BCL2", line=1, outer = TRUE, adj = 0.93, cex = 1)
  mtext(text = paste("N =", sum(c(bc2c2erpag1.surv$n, bc2c2erpag2.surv$n, bc2c2erpag3.surv$n, bc2c2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.903, cex = 0.7)
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()



#### EZH2, H3K27me3, FOXA1, BCL2 stratified by age and ER-



pdf(file = paste(plot_dirs, "EZH2_H3K27me3_FOXA1_BCL2_ERneg_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "EZH2_H3K27me3_FOXA1_BCL2_ERneg_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  ##   bc2c2fn  ezc2
  plot(ezc2ernag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40", paste(": N =", paste(ezc2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(ezc2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = ezc2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(ezc2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(ezc2ernag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(ezc2ernag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  ##   bc2c2fp  h3c2
  plot(h3c2ernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(h3c2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(h3c2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(h3c2ernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(h3c2ernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  ## fx1c2bn  fx1c2
  plot(fx1c2ernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(fx1c2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(fx1c2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(fx1c2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(fx1c2ernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(fx1c2ernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  ## fx1c2bp  bc2c2
  plot(bc2c2ernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(bc2c2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(bc2c2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = bc2c2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(bc2c2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(bc2c2ernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(bc2c2ernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "ER-:  EZH2", line=1, outer = TRUE, adj = 0.095, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2ernag1.surv$n, ezc2ernag2.surv$n, ezc2ernag3.surv$n, ezc2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.117, cex = 0.7)
  mtext(text = "ER-:  H3K27me3", line=1, outer = TRUE, adj = 0.36, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2ernag1.surv$n, h3c2ernag2.surv$n, h3c2ernag3.surv$n, h3c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.377, cex = 0.7)
  mtext(text = "ER-:  FOXA1", line=1, outer = TRUE, adj = 0.65, cex = 1)
  mtext(text = paste("N =", sum(c(fx1c2ernag1.surv$n, fx1c2ernag2.surv$n, fx1c2ernag3.surv$n, fx1c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.637, cex = 0.7)
  mtext(text = "ER-:  BCL2", line=1, outer = TRUE, adj = 0.93, cex = 1)
  mtext(text = paste("N =", sum(c(bc2c2ernag1.surv$n, bc2c2ernag2.surv$n, bc2c2ernag3.surv$n, bc2c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.903, cex = 0.7)
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()

############

#### EZH2, H3K27me3 stratified by age and ER+ (left 2 cols) and ER- (right 2 cols)


pdf(file = paste(plot_dirs, "EZH2_H3K27me3_ERposneg_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "EZH2_H3K27me3_ERposneg_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
  ## ER+
  ##   bc2c2fn  ezc2
  plot(ezc2erpag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40", paste(": N =", paste(ezc2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(ezc2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = ezc2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(ezc2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(ezc2erpag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(ezc2erpag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2erpag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(h3c2erpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(h3c2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(h3c2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(h3c2erpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(h3c2erpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2erpag4.lrpvtxt, adj=0, cex = 0.9)
  
  
  ## ER-
  
  plot(ezc2ernag1.surv, lty = 1:2, lwd = 2, main = paste("(20-40]", paste(": N =", paste(ezc2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(ezc2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = ezc2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("EZH2 neg", "EZH2 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(ezc2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(ezc2ernag3.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(ezc2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(ezc2ernag4.surv$n, collapse = ", "))),  
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = ezc2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(h3c2ernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(h3c2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(h3c2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.3, labels = h3c2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$h3k27me3_v1.1.pp_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(1, 2)])
  
  plot(h3c2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(h3c2ernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(h3c2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(h3c2ernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(1, 2)])
  text(x=0.5, y=0.1, labels = h3c2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "ER+:  EZH2", line=1, outer = TRUE, adj = 0.095, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2erpag1.surv$n, ezc2erpag2.surv$n, ezc2erpag3.surv$n, ezc2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.117, cex = 0.7)
  mtext(text = "ER+:  H3K27me3", line=1, outer = TRUE, adj = 0.36, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2erpag1.surv$n, h3c2erpag2.surv$n, h3c2erpag3.surv$n, h3c2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.377, cex = 0.7)
  mtext(text = "ER-:  EZH2", line=1, outer = TRUE, adj = 0.65, cex = 1)
  mtext(text = paste("N =", sum(c(ezc2ernag1.surv$n, ezc2ernag2.surv$n, ezc2ernag3.surv$n, ezc2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.637, cex = 0.7)
  mtext(text = "ER-:  H3K27me3", line=1, outer = TRUE, adj = 0.97, cex = 1)
  mtext(text = paste("N =", sum(c(h3c2ernag1.surv$n, h3c2ernag2.surv$n, h3c2ernag3.surv$n, h3c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.903, cex = 0.7)
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  # 
  # mtext(text = "ER-:  EZH2", line=1, outer = TRUE, adj = 0.095, cex = 1)
  # mtext(text = paste("N =", sum(c(ezc2ernag1.surv$n, ezc2ernag2.surv$n, ezc2ernag3.surv$n, ezc2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.117, cex = 0.7)
  # mtext(text = "ER-:  H3K27me3", line=1, outer = TRUE, adj = 0.36, cex = 1)
  # mtext(text = paste("N =", sum(c(h3c2ernag1.surv$n, h3c2ernag2.surv$n, h3c2ernag3.surv$n, h3c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.377, cex = 0.7)
  # mtext(text = "ER-:  FOXA1", line=1, outer = TRUE, adj = 0.65, cex = 1)
  # mtext(text = paste("N =", sum(c(fx1c2ernag1.surv$n, fx1c2ernag2.surv$n, fx1c2ernag3.surv$n, fx1c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.637, cex = 0.7)
  # mtext(text = "ER-:  BCL2", line=1, outer = TRUE, adj = 0.93, cex = 1)
  # mtext(text = paste("N =", sum(c(bc2c2ernag1.surv$n, bc2c2ernag2.surv$n, bc2c2ernag3.surv$n, bc2c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.903, cex = 0.7)
  # mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  # mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  # 
}
dev.off()

### FOXA1, BCL2 ER+ left two cols, ER- right two cols

#### EZH2, H3K27me3, FOXA1, BCL2 stratified by age and ER+


pdf(file = paste(plot_dirs, "FOXA1_BCL2_ERposneg_Age_SurvivalPlots.pdf", sep = "/"), useDingbats = FALSE, width = 7, height = 8, paper = "letter")
png(file = paste(plot_dirs, "FOXA1_BCL2_ERposneg_Age_SurvivalPlots.png", sep = "/"),  width = 2200, height = 2400, res = 300)

{
  par(mfcol=c(4, 4), mar = c(2, 2, 3, 1), oma = c(3, 3, 3, 1))
  
 
  ##  foxa1 bcl2 er+
  plot(fx1c2erpag1.surv, lty = 1:2, lwd = 2, main = paste("Age 20-40", paste(": N =", paste(fx1c2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(fx1c2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(fx1c2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(fx1c2erpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(fx1c2erpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2erpag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(bc2c2erpag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(bc2c2erpag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2erpag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2erpag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(bc2c2erpag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = bc2c2erpag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(bc2c2erpag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(bc2c2erpag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2erpag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2erpag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(bc2c2erpag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2erpag4.lrpvtxt, adj=0, cex = 0.9)

  ## foxa1 bcl2 er-
  plot(fx1c2ernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(fx1c2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(fx1c2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = fx1c2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright", legend = c("FOXA1 neg", "FOXA1 pos")  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(fx1c2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(fx1c2ernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(fx1c2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(fx1c2ernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = fx1c2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  plot(bc2c2ernag1.surv, lty = 1:2, lwd = 2, main = paste("(20,40]", paste(": N =", paste(bc2c2ernag1.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ernag1.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2ernag2.surv, lty = 1:2, lwd = 2, main = paste("(40,55]", paste(": N =", paste(bc2c2ernag2.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.3, labels = bc2c2ernag2.lrpvtxt, adj=0, cex = 0.9)
  legend("bottomright",legend = levels(as.factor(wmadf$bcl2_c_v2.ppnc_b60f))  , lty = 1:2, lwd = 2, seg.len = 4, cex = 0.7, col = Dark2[c(3, 4)])
  
  plot(bc2c2ernag3.surv, lty = 1:2, lwd = 2, main = paste("(55,70]", paste(": N =", paste(bc2c2ernag3.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ernag3.lrpvtxt, adj=0, cex = 0.9)
  plot(bc2c2ernag4.surv, lty = 1:2, lwd = 2, main = paste("(70,99]", paste(": N =", paste(bc2c2ernag4.surv$n, collapse = ", "))), 
       xlim = c(0, 14), mark.time = TRUE, cex=0.7, cex.axis = 0.8, cex.main = 0.9, col = Dark2[c(3, 4)])
  text(x=0.5, y=0.1, labels = bc2c2ernag4.lrpvtxt, adj=0, cex = 0.9)
  
  mtext(text = "ER+:  FOXA1", line=1, outer = TRUE, adj = 0.085, cex = 1)
  mtext(text = paste("N =", sum(c(fx1c2erpag1.surv$n, fx1c2erpag2.surv$n, fx1c2erpag3.surv$n, fx1c2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.112, cex = 0.7)
  mtext(text = "ER+:  BCL2", line=1, outer = TRUE, adj = 0.367, cex = 1)
  mtext(text = paste("N =", sum(c(bc2c2erpag1.surv$n, bc2c2erpag2.surv$n, bc2c2erpag3.surv$n, bc2c2erpag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.377, cex = 0.7)
  mtext(text = "ER-:  FOXA1", line=1, outer = TRUE, adj = 0.660, cex = 1)
  mtext(text = paste("N =", sum(c(fx1c2ernag1.surv$n, fx1c2ernag2.surv$n, fx1c2ernag3.surv$n, fx1c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.637, cex = 0.7)
  mtext(text = "ER-:  BCL2", line=1, outer = TRUE, adj = 0.935, cex = 1)
  mtext(text = paste("N =", sum(c(bc2c2ernag1.surv$n, bc2c2ernag2.surv$n, bc2c2ernag3.surv$n, bc2c2ernag4.surv$n), na.rm = TRUE)), line=0, outer = TRUE, adj = 0.908, cex = 0.7)
  
  
  mtext(text = "Years post diagnosis", side = 1, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
  mtext(text = "Overall survival rate", side = 2, line = 1.5, outer = TRUE, adj = 0.5, cex = 1)
}
dev.off()
