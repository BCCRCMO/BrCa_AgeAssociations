### Review patients, age distribution

ptadf <- read.table(file = "../../NormalBreast13Creation/Patient_Ages.csv",
                    sep = ",", stringsAsFactors = FALSE, header = TRUE)
ptadf$PatientID <- ptadf$PatID
ptcdf <- read.table(file = "../../NormalBreast13Creation/Normal_Breast_TMA_SMmods.csv",
                    sep = ",", stringsAsFactors = FALSE, header = TRUE)
ptcdf$BeatExcelIntoSubmission <- NULL

ptcdf$DiagnosisYear <- 2000 + ptcdf$YrDx
ptcdf$TMAnumber <- ptcdf$SampleinID

## Drop cases without FFPE block or H&E slide
ptcdf <- ptcdf[!(toupper(ptcdf$SlideID) == toupper("No slide")), ]
ptcdf <- ptcdf[!(toupper(ptcdf$BlockID) == toupper("No block")), ]
## Drop cases with no identifiable spots to take a core.
ptcdf <- ptcdf[!is.na(ptcdf$Ncores), ]
ptcdf <- ptcdf[!(ptcdf$Ncores == "0"), ]



### {set.seed(2431)
###  pids <- unique(ptcdf$PatientID)
###  for (pti in seq(along = pids) ) {
###    ptidxs <- which(ptcdf$PatientID %in% pids[pti])
###    npts <- length(ptidxs)
###    ptcdf$TMAnumber[ptidxs] <- sample.int(npts, npts)
###    ## If the patient has a sample that will admit four cores, that sample should be on the first TMA
###    if ( any( ptcdf$Ncores[ptidxs] == 4 ) && ( ptcdf$Ncores[ptidxs][ptcdf$TMAnumber[ptidxs] == 1] < 4 ) ) {
###      TMAnum <- ptcdf$TMAnumber[ptidxs]
###      tma1ptidx <- which( ( ( ptcdf$PatientID %in% pids[pti] ) & ( ptcdf$Ncores == 4 ) ) )[1]
###      ptcdf$TMAnumber[tma1ptidx] <- 1
###      ptcdf$TMAnumber[setdiff(ptidxs, tma1ptidx)] <- TMAnum[TMAnum != 1]
###    }
###  }
### }
### head(ptcdf, 15)
### 
### write.csv(ptcdf, file = "Randomized2431NormalBreast13_v02.csv")
### 
### ptcwdf <- reshape(data = ptcdf, timevar = "TMAnumber", idvar = "PatientID", direction = "wide")
### batch1 <- which(ptcwdf$Ncores.1 == 4)
### batch2 <- which(ptcwdf$Ncores.1 != 4)
### { set.seed(2431)
###   ptcwdf$NB13ID <- sample(nrow(ptcwdf))
###   ptcwdf[batch1, "NB13ID"] <- sample(length(batch1))
###   if ( length(batch2) ) { ptcwdf[batch2, "NB13ID"] <- ( length(batch1) + sample(length(batch2)) ) }
### }
### 
### varcols <- c("NB13ID", "PatientID", "DiagnosisYear.1", "Age.1", 
###              "BlockID.1", "Ncores.1", "Comments.1", "Slide.box.number.1", "Block.case.number.1", 
###              "BlockID.2", "Ncores.2", "Comments.2", "Slide.box.number.2", "Block.case.number.2", 
###              "BlockID.3", "Ncores.3", "Comments.3", "Slide.box.number.3", "Block.case.number.3", 
###              "BlockID.4", "Ncores.4", "Comments.4", "Slide.box.number.4", "Block.case.number.4"
###   )
### 
###    
### 
### write.csv(ptcwdf, file = "Randomized2431Wide_all_NormalBreast13_v01.csv")
### write.csv(ptcwdf[, varcols], file = "Randomized2431Wide_NormalBreast13_v02.csv")

### uptcdf <- ptcdf[( ( ptcdf$SampleinID == 0 ) & ( !is.na(ptcdf$Age) ) ), ]
### dim(uptcdf)
### 
### table(uptcdf$Age)
### 
### quartz()
### hist(uptcdf$Age)
### hist(uptcdf$Age, breaks = 14:87)

### ### Left - Right comparisons within case
### ### How many cases have a left and right sample?
### 
### ptcwdf$LeftSampleAvailable <-
###   sapply(seq(nrow(ptcwdf)), function(x) {
###     (length(grep(toupper("L"), toupper(ptcwdf[x, c("BlockID.1", "BlockID.2", "BlockID.3", "BlockID.4")]))) > 0) * 1
###   } )
### 
### 
### ptcwdf$RightSampleAvailable <- 
###   sapply(seq(nrow(ptcwdf)), function(x) {
###     (length(grep(toupper("R"), toupper(ptcwdf[x, c("BlockID.1", "BlockID.2", "BlockID.3", "BlockID.4")]))) > 0) * 1
###   } )
### 
### ptcwdf$LeftAndRightSampleAvailable <-   ptcwdf$LeftSampleAvailable * ptcwdf$RightSampleAvailable
### write.csv(ptcwdf, file = "Randomized2431Wide_allLR_NormalBreast13_v01.csv")

### Left vs right comparisons desired.  Need to randomize.
### Unit is patient left or patient right.
### Patients with a left and a right must be randomized to the same chip, but not side by side.
### 84 samples per slide, duplicate cores, so 42 patient samples per slide.
### Rows 2 through 8 (7 rows) Cols 1 to 12
### About 863 evaluable patient samples, will need about 21 TMA blocks.

### Tricky PatientIDs - Laterality Unknown or Unilateral, multiple samples. 
### "022-01","217-03", "261-00", "285-00", "50-01", "51-01", "19-00"

### Drop the SlideID "19(L)-00-2"
### Drop the SlideID "19(R)-00-2"
### Drop the SlideID "261-00-2"
### Drop the SlideID "285(R)-00"
### Drop the SlideID "022-01-1"
### Drop the SlideID "50-01-1"
### Drop the SlideID "51-01-2"
### Drop the SlideID "217(L)-03"

dim(ptcdf) ## [1] 863  36

SlideIDsToDrop <-
  c("19(L)-00-2", "19(R)-00-2", "261-00-2", "285(R)-00", "022-01-1", "50-01-1", "51-01-2", "217(L)-03")

ptcdf <- ptcdf[!(ptcdf$SlideID %in% SlideIDsToDrop), ]

dim(ptcdf) ## [1] 855  36

### Need all 4 core samples to go on first blocks

ptc4df <- ptcdf[ptcdf$Ncores != "2", ]

### Need all 2 core samples to go on last blocks

ptc2df <- ptcdf[ptcdf$Ncores == "2", ]


NULL
### 1.5mm core slide can hold 84 cores. 2 per patient => 42 patients per TMA slide/block.

### Need a list of all blocks that have L and and a list having R but not (LR)
### so must e.g. grep for (L) and (R) to avoid (LR).
### Need a list of slides that have L and a list having R
### Then merge all these lists by patient ID.
blockL.df <- ptcdf[grepl(toupper("\\(L\\)"), toupper(ptcdf$BlockID)), c("PatientID", "BlockID")]
names(blockL.df)[names(blockL.df) == "BlockID"] <- "BlockID.L"
## Drop the second block for PatientID "19-00"
blockL.df <- blockL.df[!(blockL.df$BlockID.L == "19(L)-00-2"), ]

blockR.df <- ptcdf[grepl(toupper("\\(R\\)"), toupper(ptcdf$BlockID)), c("PatientID", "BlockID")]
names(blockR.df)[names(blockR.df) == "BlockID"] <- "BlockID.R"
## Drop the second block for PatientID "19-00"
blockR.df <- blockR.df[!(blockR.df$BlockID.R == "19(R)-00-2"), ]

blockLR.df <- merge(blockL.df, blockR.df)

slideL.df <- ptcdf[grepl(toupper("\\(L\\)"), toupper(ptcdf$SlideID)), c("PatientID", "SlideID")]
names(slideL.df)[names(slideL.df) == "SlideID"] <- "SlideID.L"
## Drop the second slide for PatientID "19-00"
slideL.df <- slideL.df[!(slideL.df$SlideID.L == "19(L)-00-2"), ]

slideR.df <- ptcdf[grepl(toupper("\\(R\\)"), toupper(ptcdf$SlideID)), c("PatientID", "SlideID")]
names(slideR.df)[names(slideR.df) == "SlideID"] <- "SlideID.R"
## Drop the second slide for PatientID "19-00"
slideR.df <- slideR.df[!(slideR.df$SlideID.R == "19(R)-00-2"), ]

slideLR.df <- merge(slideL.df, slideR.df)

blockslideLR.df <- merge(blockLR.df, slideLR.df)

### 318 cases have left and right blocks.

### Need to pre-specify which sample is evaluated per patient when only one sample per patient.
### is used in analysis.  20130627:  This will need to be redone for the second set of blocks
set.seed(719); blockslideLR.df$SingleSample <- c("L", "R")[sample(c(1, 2), nrow(blockslideLR.df), replace = TRUE)]
head(blockslideLR.df)
## Assign 16 cases to first 18 blocks and 15 cases to 19 and 20.
TMABlockNos <- c(rep(seq(18), each = 16), rep(19:20, each = 15))
table(TMABlockNos)
### TMABlockNos
###  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
### 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 15 15 
length(TMABlockNos)
### [1] 318
set.seed(985); blockslideLR.df$TMABlockn <- sample(TMABlockNos)

### Now need list of cases without laterality
blockslideU.df <- ptcdf[!(ptcdf$PatientID %in% blockslideLR.df$PatientID), c("PatientID", "BlockID", "SlideID")]
dim(blockslideU.df)  ## [1] 219   3   ## Was [1] 225   3 before culling SlideIDsToDrop above.
blockslideU.df$SingleSample <- rep("U", nrow(blockslideU.df))  ## Unilateral or Unknown.
## Assign 10 cases to first 18 blocks, 12 cases to 19 and 20 and 15 (was 21) cases to slide 21.
### TMABlockUNos <- c(rep(seq(18), each = 10), rep(19:20, each = 12), rep(21, 21))
TMABlockUNos <- c(rep(seq(18), each = 10), rep(19:20, each = 12), rep(21, 15))
table(TMABlockUNos)
set.seed(763); blockslideU.df$TMABlockn <- sample(TMABlockUNos)

ptcdf[!(ptcdf$BlockID %in% c(blockslideU.df$BlockID, blockslideLR.df$BlockID.L, blockslideLR.df$BlockID.R)), ]
### <0 rows> (or 0-length row.names) after deleting SlideIDs in SlideIDsToDrop above.
### Was as below before deleting.
###      BlockID    SlideID Ncores Comments Slide.box.number Block.tray.number
### 5 19(L)-00-2 19(L)-00-2      4                         1                 1
### 7 19(R)-00-2 19(R)-00-2      4                         1                 1

### Now rbind samples to go on TMA
nbtmaIDL.df <- blockslideLR.df[, c("PatientID", "BlockID.L", "SlideID.L", "SingleSample", "TMABlockn")]
names(nbtmaIDL.df)[names(nbtmaIDL.df) == "BlockID.L"] <- "BlockID"
names(nbtmaIDL.df)[names(nbtmaIDL.df) == "SlideID.L"] <- "SlideID"
nbtmaIDL.df$Laterality <- rep("L", nrow(nbtmaIDL.df))

nbtmaIDR.df <- blockslideLR.df[, c("PatientID", "BlockID.R", "SlideID.R", "SingleSample", "TMABlockn")]
names(nbtmaIDR.df)[names(nbtmaIDR.df) == "BlockID.R"] <- "BlockID"
names(nbtmaIDR.df)[names(nbtmaIDR.df) == "SlideID.R"] <- "SlideID"
nbtmaIDR.df$Laterality <- rep("R", nrow(nbtmaIDR.df))

nbtmaIDU.df <- blockslideU.df[, c("PatientID", "BlockID", "SlideID", "SingleSample", "TMABlockn")]
nbtmaIDU.df$Laterality <- rep("U", nrow(nbtmaIDU.df))

head(nbtmaIDL.df)
head(nbtmaIDR.df)
head(nbtmaIDU.df)

nbtmaID.df <- rbind(nbtmaIDL.df, nbtmaIDR.df, nbtmaIDU.df)
dim(nbtmaID.df)
### Need also Ncores data available
nbtmaIDbkup.df  <- nbtmaID.df
nbtmaID.df <- merge(nbtmaIDbkup.df, ptcdf[, c("BlockID", "Ncores")], by = "BlockID", all.x = TRUE, all.y = FALSE, suffixes = c("", ".Y"))

NULL
### ## 16 pairs on first 18 slides.  15 pairs on slides 19,20.  318 paired.                   ## 16*18 + 15*2       [1] 318
### ## 10 unpaired on first 18 slides. 12 unpaired on slides 19, 20. 21 unpaired on slide 21. ## 10*18 + 12*2 + 21  [1] 225
## 16 pairs on first 18 slides.  15 pairs on slides 19,20.  318 paired.                   ## 16*18 + 15*2       [1] 318
## 10 unpaired on first 18 slides. 12 unpaired on slides 19, 20. 15 unpaired on slide 21. ## 10*18 + 12*2 + 15  [1] 219
### > table(nbtmaID.df$Laterality)
### ### 
### ###   L   R   U 
### ### 318 318 225  ## Before culling
###   L   R   U 
### 318 318 219   ## After culling

### Now need to place samples on the TMA.  Then assign nb13_id
withinTMABlockID <- rep(seq(42), times = 21)[seq(nrow(nbtmaID.df))]
{ nbtmaID.df$withinTMABlockID <- withinTMABlockID
  nbtmaID.df[order(nbtmaID.df$TMABlockn), "withinTMABlockID"] <- withinTMABlockID }
with(nbtmaID.df, table(withinTMABlockID, TMABlockn))

nbtmaID.df[order(nbtmaID.df$PatientID), ][nbtmaID.df[order(nbtmaID.df$PatientID), ]$TMABlockn == 1, ]
### Too structured within block.  Need to randomize within TMAblock.
### 20130627:  Need to put all the 2-core cases at the end of each TMA block so that
###            there are not a bunch of holes in the second set of blocks.
nbtmaID.df[nbtmaID.df$TMABlockn == 1, "withinTMABlockID"]
{ set.seed(5301)
  tids <- unique(nbtmaID.df$TMABlockn)
  for ( ti in seq(along = tids) ) {
    tmablki <- tids[ti]
###     tmablkiposns <- nbtmaID.df[nbtmaID.df$TMABlockn == tmablki, "withinTMABlockID"]
###     nbtmaID.df[nbtmaID.df$TMABlockn == tmablki, "withinTMABlockID"] <- sample(tmablkiposns)
    tmablki4cposns <- nbtmaID.df[((nbtmaID.df$TMABlockn == tmablki) & (nbtmaID.df$Ncores != "2")), "withinTMABlockID"]
    numtmablki4cposns <- length(tmablki4cposns)
    if ( numtmablki4cposns > 0 ) {
      nbtmaID.df[((nbtmaID.df$TMABlockn == tmablki) & (nbtmaID.df$Ncores != "2")),
                 "withinTMABlockID"] <- sample(seq(numtmablki4cposns))
    }
    tmablki2cposns <- nbtmaID.df[((nbtmaID.df$TMABlockn == tmablki) & (nbtmaID.df$Ncores == "2")), "withinTMABlockID"]
    numtmablki2cposns <- length(tmablki2cposns)
    if ( numtmablki2cposns > 0 ) {
      nbtmaID.df[((nbtmaID.df$TMABlockn == tmablki) & (nbtmaID.df$Ncores == "2")),
                 "withinTMABlockID"] <- sample(seq(numtmablki2cposns) + numtmablki4cposns)
    }    
  }
  ## Sort and assign ID
  nbtmaID.df <- nbtmaID.df[order(nbtmaID.df$TMABlockn, nbtmaID.df$withinTMABlockID), ]
  nbtmaID.df$nb13_id <- seq(nrow(nbtmaID.df))
}
nbtmaID.df$TMABlock <- LETTERS[nbtmaID.df$TMABlockn]
nbtmaID.df$isUnpairedp <- (nbtmaID.df$SingleSample == nbtmaID.df$Laterality)
nbtmaID.df$isUnpaired <- 1.0 * (nbtmaID.df$SingleSample == nbtmaID.df$Laterality)
### Merge data in ptcdf back with nbtmaID.df

nb13tmadf <- merge(x = nbtmaID.df, y = ptcdf, by = "BlockID", all.x = TRUE, all.y = FALSE, suffixes = c("", ".Y"))
### Find all pairs - if one has 4 cores and other has 2 cores, make sure the 4 core sample is TRUE for isUnpairedp

{
  nb13PtIDs <- unique(nb13tmadf$PatientID)
  for ( pti in seq(along = nb13PtIDs) ) {
    ptIDi <- nb13PtIDs[pti]
    ptIDi.idxs <- which(nb13tmadf$PatientID %in% ptIDi)
    if ( length(ptIDi.idxs) == 2 && any(nb13tmadf$Ncores[ptIDi.idxs] == "2") && any(nb13tmadf$Ncores[ptIDi.idxs] != "2")) {
      nb13tmadf[((nb13tmadf$PatientID %in% ptIDi) & (nb13tmadf$Ncores == "2")), "isUnpairedp"] <- FALSE
      nb13tmadf[((nb13tmadf$PatientID %in% ptIDi) & (nb13tmadf$Ncores != "2")), "isUnpairedp"] <- TRUE
    }
  }
  nb13tmadf$isUnpaired <- 1.0 * nb13tmadf$isUnpairedp
}

dim(nb13tmadf)  ## [1] 855  46 after culling
head(nb13tmadf)

### Tricky PatientIDs - Laterality Unknown or Unilateral, multiple samples. 
### "022-01","217-03", "261-00", "285-00", "50-01", "51-01", "19-00"
table(table(nb13tmadf[nb13tmadf$isUnpairedp, "PatientID"]))
###   1 
### 537   ## After culling
table(table(nb13tmadf[!nb13tmadf$isUnpairedp, "PatientID"]))
###   1 
### 318

write.csv(nb13tmadf[order(nb13tmadf$TMABlockn,nb13tmadf$withinTMABlockID),
                    c("BlockID", "Age", "Ncores", "Comments", "Slide.box.number",
                      "Block.tray.number", "Block.tray.column.number", "PatientID", "SlideID", "TMABlock",
                      "withinTMABlockID", "nb13_id", "TMABlockn", "DiagnosisYear", "Laterality",
                      "isUnpairedp", "isUnpaired"),],
          file = "../../NormalBreast13Creation/RandomizedLat_NormalBreastTMA_NB13_v03.csv")
                    


NULL

###--------------------------------------------------------------------------------------------------------
### 4 core samples
### Code below revealed that attempting to place all 2 punch cases at end of whole TMA would yield a loss
### of 31% of pairings.  Code above modified to place 2 punch cases at end of each slide.

### Code below not used to set up TMA.

NULL
### 1.5mm core slide can hold 84 cores. 2 per patient => 42 patients per TMA slide/block.

### Need a list of all blocks that have L and and a list having R but not (LR)
### so must e.g. grep for (L) and (R) to avoid (LR).
### Need a list of slides that have L and a list having R
### Then merge all these lists by patient ID.
blockL.4df <- ptc4df[grepl(toupper("\\(L\\)"), toupper(ptc4df$BlockID)), c("PatientID", "BlockID")]
names(blockL.4df)[names(blockL.4df) == "BlockID"] <- "BlockID.L"
## Drop the second block for PatientID "19-00"
blockL.4df <- blockL.4df[!(blockL.4df$BlockID.L == "19(L)-00-2"), ]

blockR.4df <- ptc4df[grepl(toupper("\\(R\\)"), toupper(ptc4df$BlockID)), c("PatientID", "BlockID")]
names(blockR.4df)[names(blockR.4df) == "BlockID"] <- "BlockID.R"
## Drop the second block for PatientID "19-00"
blockR.4df <- blockR.4df[!(blockR.4df$BlockID.R == "19(R)-00-2"), ]

blockLR.4df <- merge(blockL.4df, blockR.4df)

slideL.4df <- ptc4df[grepl(toupper("\\(L\\)"), toupper(ptc4df$SlideID)), c("PatientID", "SlideID")]
names(slideL.4df)[names(slideL.4df) == "SlideID"] <- "SlideID.L"
## Drop the second slide for PatientID "19-00"
slideL.4df <- slideL.4df[!(slideL.4df$SlideID.L == "19(L)-00-2"), ]

slideR.4df <- ptc4df[grepl(toupper("\\(R\\)"), toupper(ptc4df$SlideID)), c("PatientID", "SlideID")]
names(slideR.4df)[names(slideR.4df) == "SlideID"] <- "SlideID.R"
## Drop the second slide for PatientID "19-00"
slideR.4df <- slideR.4df[!(slideR.4df$SlideID.R == "19(R)-00-2"), ]

slideLR.4df <- merge(slideL.4df, slideR.4df)

blockslideLR.4df <- merge(blockLR.4df, slideLR.4df)

### 201 cases have left and right blocks.

### Need to pre-specify which sample is evaluated per patient when only one sample per patient
### is used in analysis.
set.seed(719); blockslideLR.4df$SingleSample <- c("L", "R")[sample(c(1, 2), nrow(blockslideLR.4df), replace = TRUE)]
head(blockslideLR.4df)
## Assign 13 cases to first 15 blocks and 6 cases to block 16.
TMABlock4Nos <- c(rep(seq(15), each = 13), rep(16, each = 6))
table(TMABlock4Nos)
### TMABlock4Nos
###  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
### 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 15 15 
length(TMABlock4Nos)
### [1] 201
set.seed(985); blockslideLR.4df$TMABlockn <- sample(TMABlock4Nos)

### Now need list of cases without laterality
blockslideU.4df <- ptc4df[!(ptc4df$PatientID %in% blockslideLR.4df$PatientID), c("PatientID", "BlockID", "SlideID")]
dim(blockslideU.4df)  ## [1] 219   3   ## Was [1] 225   3 before culling SlideIDsToDrop above.
blockslideU.4df$SingleSample <- rep("U", nrow(blockslideU.4df))  ## Unilateral or Unknown.
## Assign 10 cases to first 18 blocks, 12 cases to 19 and 20 and 15 (was 21) cases to slide 21.
### TMABlockU4Nos <- c(rep(seq(18), each = 10), rep(19:20, each = 12), rep(21, 21))
TMABlockU4Nos <- c(rep(seq(18), each = 10), rep(19:20, each = 12), rep(21, 15))
table(TMABlockU4Nos)
set.seed(763); blockslideU.4df$TMABlockn <- sample(TMABlockU4Nos)

ptc4df[!(ptc4df$BlockID %in% c(blockslideU.4df$BlockID, blockslideLR.4df$BlockID.L, blockslideLR.4df$BlockID.R)), ]
### <0 rows> (or 0-length row.names) after deleting SlideIDs in SlideIDsToDrop above.
### Was as below before deleting.
###      BlockID    SlideID Ncores Comments Slide.box.number Block.tray.number
### 5 19(L)-00-2 19(L)-00-2      4                         1                 1
### 7 19(R)-00-2 19(R)-00-2      4                         1                 1

### Now rbind samples to go on TMA
nbtmaIDL.4df <- blockslideLR.4df[, c("PatientID", "BlockID.L", "SlideID.L", "SingleSample", "TMABlockn")]
names(nbtmaIDL.4df)[names(nbtmaIDL.4df) == "BlockID.L"] <- "BlockID"
names(nbtmaIDL.4df)[names(nbtmaIDL.4df) == "SlideID.L"] <- "SlideID"
nbtmaIDL.4df$Laterality <- rep("L", nrow(nbtmaIDL.4df))

nbtmaIDR.4df <- blockslideLR.4df[, c("PatientID", "BlockID.R", "SlideID.R", "SingleSample", "TMABlockn")]
names(nbtmaIDR.4df)[names(nbtmaIDR.4df) == "BlockID.R"] <- "BlockID"
names(nbtmaIDR.4df)[names(nbtmaIDR.4df) == "SlideID.R"] <- "SlideID"
nbtmaIDR.4df$Laterality <- rep("R", nrow(nbtmaIDR.4df))

nbtmaIDU.4df <- blockslideU.4df[, c("PatientID", "BlockID", "SlideID", "SingleSample", "TMABlockn")]
nbtmaIDU.4df$Laterality <- rep("U", nrow(nbtmaIDU.4df))

head(nbtmaIDL.4df)
head(nbtmaIDR.4df)
head(nbtmaIDU.4df)

nbtmaID.4df <- rbind(nbtmaIDL.4df, nbtmaIDR.4df, nbtmaIDU.4df)
dim(nbtmaID.4df)

### ## 16 pairs on first 18 slides.  15 pairs on slides 19,20.  318 paired.                   ## 16*18 + 15*2       [1] 318
### ## 10 unpaired on first 18 slides. 12 unpaired on slides 19, 20. 21 unpaired on slide 21. ## 10*18 + 12*2 + 21  [1] 225
## 16 pairs on first 18 slides.  15 pairs on slides 19,20.  318 paired.                   ## 16*18 + 15*2       [1] 318
## 10 unpaired on first 18 slides. 12 unpaired on slides 19, 20. 15 unpaired on slide 21. ## 10*18 + 12*2 + 15  [1] 219
### > table(nbtmaID.4df$Laterality)
### ### 
### ###   L   R   U 
### ### 318 318 225  ## Before culling
###   L   R   U 
### 318 318 219   ## After culling

### Now need to place samples on the TMA.  Then assign nb13_id
withinTMABlockID <- rep(seq(42), times = 21)[seq(nrow(nbtmaID.4df))]
{ nbtmaID.4df$withinTMABlockID <- withinTMABlockID
  nbtmaID.4df[order(nbtmaID.4df$TMABlockn), "withinTMABlockID"] <- withinTMABlockID }
with(nbtmaID.4df, table(withinTMABlockID, TMABlockn))

nbtmaID.4df[order(nbtmaID.4df$PatientID), ][nbtmaID.4df[order(nbtmaID.4df$PatientID), ]$TMABlockn == 1, ]
### Too structured within block.  Need to randomize within TMAblock.
nbtmaID.4df[nbtmaID.4df$TMABlockn == 1, "withinTMABlockID"]
{ set.seed(5301)
  tids <- unique(nbtmaID.4df$TMABlockn)
  for ( ti in seq(along = tids) ) {
    tmablki <- tids[ti]
    tmablkiposns <- nbtmaID.4df[nbtmaID.4df$TMABlockn == tmablki, "withinTMABlockID"]
    nbtmaID.4df[nbtmaID.4df$TMABlockn == tmablki, "withinTMABlockID"] <- sample(tmablkiposns)
  }
  ## Sort and assign ID
  nbtmaID.4df <- nbtmaID.4df[order(nbtmaID.4df$TMABlockn, nbtmaID.4df$withinTMABlockID), ]
  nbtmaID.4df$nb13_id <- seq(nrow(nbtmaID.4df))
}
nbtmaID.4df$TMABlock <- LETTERS[nbtmaID.4df$TMABlockn]
nbtmaID.4df$isUnpairedp <- (nbtmaID.4df$SingleSample == nbtmaID.4df$Laterality)
nbtmaID.4df$isUnpaired <- 1.0 * (nbtmaID.4df$SingleSample == nbtmaID.4df$Laterality)
### Merge data in ptc4df back with nbtmaID.4df

nb13tma4df <- merge(x = nbtmaID.4df, y = ptc4df, by = "BlockID", all.x = TRUE, all.y = FALSE, suffixes = c("", ".Y"))
dim(nb13tma4df)  ## [1] 855  46 after culling
head(nb13tma4df)

### Tricky PatientIDs - Laterality Unknown or Unilateral, multiple samples. 
### "022-01","217-03", "261-00", "285-00", "50-01", "51-01", "19-00"
table(table(nb13tma4df[nb13tma4df$isUnpairedp, "PatientID"]))
###   1 
### 537   ## After culling
table(table(nb13tma4df[!nb13tma4df$isUnpairedp, "PatientID"]))
###   1 
### 318

write.csv(nb13tma4df[order(nb13tma4df$TMABlockn,nb13tma4df$withinTMABlockID),
                    c("BlockID", "Age", "Ncores", "Comments", "Slide.box.number",
                      "Block.tray.number", "Block.tray.column.number", "PatientID", "SlideID", "TMABlock",
                      "withinTMABlockID", "nb13_id", "TMABlockn", "DiagnosisYear", "Laterality",
                      "isUnpairedp", "isUnpaired"),],
          file = "../../NormalBreast13Creation/RandomizedLat4core_NormalBreastTMA_NB13_v01.csv")
                    

### 2 core samples
NULL
### 1.5mm core slide can hold 84 cores. 2 per patient => 42 patients per TMA slide/block.

### Need a list of all blocks that have L and and a list having R but not (LR)
### so must e.g. grep for (L) and (R) to avoid (LR).
### Need a list of slides that have L and a list having R
### Then merge all these lists by patient ID.
blockL.2df <- ptc2df[grepl(toupper("\\(L\\)"), toupper(ptc2df$BlockID)), c("PatientID", "BlockID")]
names(blockL.2df)[names(blockL.2df) == "BlockID"] <- "BlockID.L"
## Drop the second block for PatientID "19-00"
blockL.2df <- blockL.2df[!(blockL.2df$BlockID.L == "19(L)-00-2"), ]

blockR.2df <- ptc2df[grepl(toupper("\\(R\\)"), toupper(ptc2df$BlockID)), c("PatientID", "BlockID")]
names(blockR.2df)[names(blockR.2df) == "BlockID"] <- "BlockID.R"
## Drop the second block for PatientID "19-00"
blockR.2df <- blockR.2df[!(blockR.2df$BlockID.R == "19(R)-00-2"), ]

blockLR.2df <- merge(blockL.2df, blockR.2df)

slideL.2df <- ptc2df[grepl(toupper("\\(L\\)"), toupper(ptc2df$SlideID)), c("PatientID", "SlideID")]
names(slideL.2df)[names(slideL.2df) == "SlideID"] <- "SlideID.L"
## Drop the second slide for PatientID "19-00"
slideL.2df <- slideL.2df[!(slideL.2df$SlideID.L == "19(L)-00-2"), ]

slideR.2df <- ptc2df[grepl(toupper("\\(R\\)"), toupper(ptc2df$SlideID)), c("PatientID", "SlideID")]
names(slideR.2df)[names(slideR.2df) == "SlideID"] <- "SlideID.R"
## Drop the second slide for PatientID "19-00"
slideR.2df <- slideR.2df[!(slideR.2df$SlideID.R == "19(R)-00-2"), ]

slideLR.2df <- merge(slideL.2df, slideR.2df)

blockslideLR.2df <- merge(blockLR.2df, slideLR.2df)

### 19 cases have left and right blocks.

### Need to pre-specify which sample is evaluated per patient when only one sample per patient
### is used in analysis.
set.seed(719); blockslideLR.2df$SingleSample <- c("L", "R")[sample(c(1, 2), nrow(blockslideLR.2df), replace = TRUE)]
head(blockslideLR.2df)
## Assign 13 cases to first 15 blocks and 6 cases to block 16.
TMABlock2Nos <- c(rep(seq(15), each = 13), rep(16, each = 6))
table(TMABlock2Nos)
### TMABlock2Nos
###  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
### 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 15 15 
length(TMABlock2Nos)
### [1] 201
set.seed(985); blockslideLR.2df$TMABlockn <- sample(TMABlock2Nos)

### Now need list of cases without laterality
blockslideU.2df <- ptc2df[!(ptc2df$PatientID %in% blockslideLR.2df$PatientID), c("PatientID", "BlockID", "SlideID")]
dim(blockslideU.2df)  ## [1] 219   3   ## Was [1] 225   3 before culling SlideIDsToDrop above.
blockslideU.2df$SingleSample <- rep("U", nrow(blockslideU.2df))  ## Unilateral or Unknown.
## Assign 10 cases to first 18 blocks, 12 cases to 19 and 20 and 15 (was 21) cases to slide 21.
### TMABlockU2Nos <- c(rep(seq(18), each = 10), rep(19:20, each = 12), rep(21, 21))
TMABlockU2Nos <- c(rep(seq(18), each = 10), rep(19:20, each = 12), rep(21, 15))
table(TMABlockU2Nos)
set.seed(763); blockslideU.2df$TMABlockn <- sample(TMABlockU2Nos)

ptc2df[!(ptc2df$BlockID %in% c(blockslideU.2df$BlockID, blockslideLR.2df$BlockID.L, blockslideLR.2df$BlockID.R)), ]
### <0 rows> (or 0-length row.names) after deleting SlideIDs in SlideIDsToDrop above.
### Was as below before deleting.
###      BlockID    SlideID Ncores Comments Slide.box.number Block.tray.number
### 5 19(L)-00-2 19(L)-00-2      4                         1                 1
### 7 19(R)-00-2 19(R)-00-2      4                         1                 1

### Now rbind samples to go on TMA
nbtmaIDL.2df <- blockslideLR.2df[, c("PatientID", "BlockID.L", "SlideID.L", "SingleSample", "TMABlockn")]
names(nbtmaIDL.2df)[names(nbtmaIDL.2df) == "BlockID.L"] <- "BlockID"
names(nbtmaIDL.2df)[names(nbtmaIDL.2df) == "SlideID.L"] <- "SlideID"
nbtmaIDL.2df$Laterality <- rep("L", nrow(nbtmaIDL.2df))

nbtmaIDR.2df <- blockslideLR.2df[, c("PatientID", "BlockID.R", "SlideID.R", "SingleSample", "TMABlockn")]
names(nbtmaIDR.2df)[names(nbtmaIDR.2df) == "BlockID.R"] <- "BlockID"
names(nbtmaIDR.2df)[names(nbtmaIDR.2df) == "SlideID.R"] <- "SlideID"
nbtmaIDR.2df$Laterality <- rep("R", nrow(nbtmaIDR.2df))

nbtmaIDU.2df <- blockslideU.2df[, c("PatientID", "BlockID", "SlideID", "SingleSample", "TMABlockn")]
nbtmaIDU.2df$Laterality <- rep("U", nrow(nbtmaIDU.2df))

head(nbtmaIDL.2df)
head(nbtmaIDR.2df)
head(nbtmaIDU.2df)

nbtmaID.2df <- rbind(nbtmaIDL.2df, nbtmaIDR.2df, nbtmaIDU.2df)
dim(nbtmaID.2df)

### ## 16 pairs on first 18 slides.  15 pairs on slides 19,20.  318 paired.                   ## 16*18 + 15*2       [1] 318
### ## 10 unpaired on first 18 slides. 12 unpaired on slides 19, 20. 21 unpaired on slide 21. ## 10*18 + 12*2 + 21  [1] 225
## 16 pairs on first 18 slides.  15 pairs on slides 19,20.  318 paired.                   ## 16*18 + 15*2       [1] 318
## 10 unpaired on first 18 slides. 12 unpaired on slides 19, 20. 15 unpaired on slide 21. ## 10*18 + 12*2 + 15  [1] 219
### > table(nbtmaID.2df$Laterality)
### ### 
### ###   L   R   U 
### ### 318 318 225  ## Before culling
###   L   R   U 
### 318 318 219   ## After culling

### Now need to place samples on the TMA.  Then assign nb13_id
withinTMABlockID <- rep(seq(42), times = 21)[seq(nrow(nbtmaID.2df))]
{ nbtmaID.2df$withinTMABlockID <- withinTMABlockID
  nbtmaID.2df[order(nbtmaID.2df$TMABlockn), "withinTMABlockID"] <- withinTMABlockID }
with(nbtmaID.2df, table(withinTMABlockID, TMABlockn))

nbtmaID.2df[order(nbtmaID.2df$PatientID), ][nbtmaID.2df[order(nbtmaID.2df$PatientID), ]$TMABlockn == 1, ]
### Too structured within block.  Need to randomize within TMAblock.
nbtmaID.2df[nbtmaID.2df$TMABlockn == 1, "withinTMABlockID"]
{ set.seed(5301)
  tids <- unique(nbtmaID.2df$TMABlockn)
  for ( ti in seq(along = tids) ) {
    tmablki <- tids[ti]
    tmablkiposns <- nbtmaID.2df[nbtmaID.2df$TMABlockn == tmablki, "withinTMABlockID"]
    nbtmaID.2df[nbtmaID.2df$TMABlockn == tmablki, "withinTMABlockID"] <- sample(tmablkiposns)
  }
  ## Sort and assign ID
  nbtmaID.2df <- nbtmaID.2df[order(nbtmaID.2df$TMABlockn, nbtmaID.2df$withinTMABlockID), ]
  nbtmaID.2df$nb13_id <- seq(nrow(nbtmaID.2df))
}
nbtmaID.2df$TMABlock <- LETTERS[nbtmaID.2df$TMABlockn]
nbtmaID.2df$isUnpairedp <- (nbtmaID.2df$SingleSample == nbtmaID.2df$Laterality)
nbtmaID.2df$isUnpaired <- 1.0 * (nbtmaID.2df$SingleSample == nbtmaID.2df$Laterality)
### Merge data in ptc2df back with nbtmaID.2df

nb13tma2df <- merge(x = nbtmaID.2df, y = ptc2df, by = "BlockID", all.x = TRUE, all.y = FALSE, suffixes = c("", ".Y"))
dim(nb13tma2df)  ## [1] 855  46 after culling
head(nb13tma2df)

### Tricky PatientIDs - Laterality Unknown or Unilateral, multiple samples. 
### "022-01","217-03", "261-00", "285-00", "50-01", "51-01", "19-00"
table(table(nb13tma2df[nb13tma2df$isUnpairedp, "PatientID"]))
###   1 
### 537   ## After culling
table(table(nb13tma2df[!nb13tma2df$isUnpairedp, "PatientID"]))
###   1 
### 318

write.csv(nb13tma2df[order(nb13tma2df$TMABlockn,nb13tma2df$withinTMABlockID),
                    c("BlockID", "Age", "Ncores", "Comments", "Slide.box.number",
                      "Block.tray.number", "Block.tray.column.number", "PatientID", "SlideID", "TMABlock",
                      "withinTMABlockID", "nb13_id", "TMABlockn", "DiagnosisYear", "Laterality",
                      "isUnpairedp", "isUnpaired"),],
          file = "../../NormalBreast13Creation/RandomizedLat4core_NormalBreastTMA_NB13_v01.csv")
                    
