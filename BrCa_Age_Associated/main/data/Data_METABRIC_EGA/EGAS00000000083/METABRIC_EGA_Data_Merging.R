
# Setup -------------------------------------------------------------------

# Load packages and set dir path for loading data
library(tidyverse)
library(magrittr)
ega.dir <- "main/data/Data_METABRIC_EGA/EGAS00000000083/"

# Merging -----------------------------------------------------------------

# Copy number data METABRIC tumour cases
mbcnaddf <- read_rds(file.path(ega.dir, "EGAD00010000213/discovery_CNA_CBS.rds"))  # d_cna
mbcnavdf <- read_rds(file.path(ega.dir, "EGAD00010000215/validation_CNA_CBS.rds"))  # v_cna
mbcnvddf <- read_rds(file.path(ega.dir, "EGAD00010000214/discovery_CNV_CBS.rds"))  # d_cnv
mbcnvvdf <- read_rds(file.path(ega.dir, "EGAD00010000216/validation_CNV_CBS.rds"))  # v_cnv

# Add MBid vars and discovery set indicator
mbcnaddf <- mutate(mbcnaddf, MBid = setMBid(METABRIC_ID), isDiscovery = TRUE)
mbcnavdf <- mutate(mbcnavdf, MBid = setMBid(METABRIC_ID), isDiscovery = FALSE)
mbcnvddf <- mutate(mbcnvddf, MBid = setMBid(METABRIC_ID), isDiscovery = TRUE)
mbcnvvdf <- mutate(mbcnvvdf, MBid = setMBid(METABRIC_ID), isDiscovery = FALSE)

# Row-bind discovery and validation data for CNA and CNV
mbcnadf <- select(rbind(mbcnaddf, mbcnavdf), MBid, everything())
mbcnvdf <- select(rbind(mbcnvddf, mbcnvvdf), MBid, everything())


# Expression Data ---------------------------------------------------------

# Discovery, Validation, Normals data from European Genome-Phenome Archive from the METABRIC study
mbexddf <- read_rds(file.path(ega.dir, "EGAD00010000210/discovery_ExpressionMatrix.rds"))
mbexvdf <- read_rds(file.path(ega.dir, "EGAD00010000211/validation_ExpressionMatrix.rds"))
mbexndf <- read_rds(file.path(ega.dir, "EGAD00010000212/normals_ExpressionMatrix.rds"))

# Note: dplyr::arrange() removes row.names
mbexddf <- mbexddf %>% 
  extract(order(row.names(.)), )
mbexvdf <- mbexvdf %>% 
  extract(order(row.names(.)), )
mbexndf <- mbexndf %>% 
  extract(order(row.names(.)), ) %>% 
  set_names(setMBid(names(.)))

# Extract expression data into matrix for faster R processing.

# Discovery set
mbexdmat <- t(as.matrix(mbexddf))
dtrtcentre <- MBid2Ctr(dimnames(mbexdmat)[[1]])

# Validation set
mbexvmat <- t(as.matrix(mbexvdf))
vtrtcentre <- MBid2Ctr(dimnames(mbexvmat)[[1]])

# Create standardized values using robust location/scale
DatasetI.rzs <- rzscore.strat.transform(mbexdmat, dtrtcentre) %>% 
  extract(, order(colnames(.)))
DatasetII.rzs <- rzscore.strat.transform(mbexvmat, vtrtcentre) %>% 
  extract(, order(colnames(.)))

Dataset.r <- mbexdvmat <- rbind(mbexdmat, mbexvmat)
Dataset.rzs <- rbind(DatasetI.rzs, DatasetII.rzs)
