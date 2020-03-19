# Setup packages and directory for saving data
library(readr)
library(fs)
load_dir <- "/Volumes/KilroyHD/kilroy/Projects/MOlab/METABRIC/ExpressionDataExampleShare/EGA_Data/EGAS00000000083"
save_dir <- "main/data/Data_METABRIC_EGA/EGAS00000000083"

# Load, Compress & Save
# Discovery expression matrix
load_path <- path(load_dir, "EGAD00010000210/discovery_ExpressionMatrix.RData")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
load(load_path)
write_rds(mbexddf, path = save_path, compress = "xz")

# Validation expression matrix
load_path <- path(load_dir, "EGAD00010000211/validation_ExpressionMatrix.RData")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
load(load_path)
write_rds(mbexvdf, path = save_path, compress = "xz")

# Normalized expression matrix
load_path <- path(load_dir, "EGAD00010000212/normals_ExpressionMatrix.RData")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
load(load_path)
write_rds(mbexndf, path = save_path, compress = "xz")

# Discovery Copy Number Abberation CBS
load_path <- path(load_dir, "EGAD00010000213/discovery_CNA_CBS.txt")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
d_cna_cbs <- read_tsv(load_path, col_types = cols(loc.end = "d"))
write_rds(d_cna_cbs, path = save_path, compress = "xz")

# Discovery Copy Number Variant CBS
load_path <- path(load_dir, "EGAD00010000214/discovery_CNV_CBS.txt")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
d_cnv_cbs <- read_tsv(load_path)
write_rds(d_cnv_cbs, path = save_path, compress = "xz")

# Validation Copy Number Abberation CBS
load_path <- path(load_dir, "EGAD00010000215/validation_CNA_CBS.txt")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
v_cna_cbs <- read_tsv(load_path)
write_rds(v_cna_cbs, path = save_path, compress = "xz")

# Validation Copy Number Variant CBS
load_path <- path(load_dir, "EGAD00010000216/validation_CNV_CBS.txt")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
v_cnv_cbs <- read_tsv(load_path)
write_rds(v_cnv_cbs, path = save_path, compress = "xz")

# Discovery Copy Number Abberation HMM
load_path <- path(load_dir, "EGAD00010000217/discovery_CNA_HMM.txt")
save_path <- path_ext_set(gsub(load_dir, save_dir, load_path), "rds")
d_cna_hmm <- read_tsv(load_path, col_names = c(
  "METABRIC_ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"))
write_rds(d_cna_hmm, path = save_path, compress = "xz")
