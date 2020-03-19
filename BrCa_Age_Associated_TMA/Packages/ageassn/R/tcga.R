#' Load TCGA data
#'
#' Load TCGA data pertaining to a specific organ into the current environment.
#'
#' The data directories need to have stable, constant paths because the function
#' hard-codes the directory and file name for each data object. The data objects
#' are:
#' * `cldf`: clinical data
#' * `ealldf`: expression data
#' * `annodf`: annotation data
#' * `ardf`: age-related genes
#' * `ezdf`: EZH2-related genes
#' * `ebdf`: ER binding genes frmo ChIP-Seq
#' * `ebdf_t1`: NKI Tier 1 ER binding genes
#' * `ebdf_t2`: NKI Tier 2 ER binding genes
#' * `icdf`: IntClust groups (only for breast cancer)
#'
#' Only `cldf` and `ealldf` are dependent on `organ`. An additional `icdf` data
#' object is loaded when `organ == "Br"`.
#'
#' @param organ organ to obtain TCGA data
#' @param spec logical; if `TRUE` and objects do not already exist, then column
#'   specification parsing output from `readr` is shown.
#' @param overwrite logical; only load if set to `TRUE` or if the object doesn't
#'   already exist in the current environment. `organ`-specific objects `cldf`
#'   and `ealldf` are always loaded.
#' @return The function is returned invisibly. All data objects are loaded into
#'   the current environment.
#' @author Derek Chiu
#' @export
tcga_load <- function(organ = c("Br", "Kidney", "Lung", "Prostate", "Thyroid"),
                      spec = FALSE, overwrite = FALSE) {
  # Choose organ
  org <- match.arg(organ)

  # Whether to show readr parsing output of column specs
  col_types <- if (spec) NULL else readr::cols()

  # Set directories
  dat_dir <- "main/data"
  exp_dir <- fs::path(dat_dir, glue::glue("Data_TCGA_{org}Ca_RNASeq_Expression_Trends"))
  arez_dir <- fs::path(dat_dir, "Data_TomoOsako_BrCa_AgeRelated_EZH2Related")
  an_dir <- fs::path(dat_dir, "Data_METABRIC_Expression_Trends")
  eb_dir <- fs::path(dat_dir, "Data_JasonCarroll")
  nki_dir <- fs::path(dat_dir, "Data_StaceyJoostenNKI_ERbinding")

  # Load organ-specific objects
  cldf <- readr::read_tsv(fs::path(exp_dir, "data_bcr_clinical_data_patient.txt"),
                         col_types = readr::cols(INITIAL_PATHOLOGIC_DX_YEAR = "c", # Missing data coded as characters
                                                 AGE = "c",
                                                 DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS = "c"),
                         skip = 4)  # Skip extra prepending metadata
  ealldf <- readr::read_rds(fs::path(exp_dir, "data_RNA_Seq_v2_expression_median.rds"))
  assign("cldf", cldf, pos = pos)
  assign("ealldf", ealldf, pos = pos)

  # Store organ-nonspecific objects in vector, conditionally add IntClust data if BrCa requested
  obj_names <- c("annodf", "ardf", "ezdf", "ebdf", "ebdf_t1", "ebdf_t2")
  if (org == "Br") obj_names <- c(obj_names, "icdf")

  purrr::walk(obj_names, ~ {
    # Only load if object doesn't exist or specify to overwrite existing object
    if (!exists(.) | overwrite) {
      obj <- switch(
        .,
        annodf = readr::read_rds(fs::path(an_dir, "Annotation_Illumina_Human-WG-V3_hg18_V1.0.0_Aug09.rds")),
        ardf = readr::read_csv(fs::path(arez_dir, "AgeRelated_Accession_Numbersv02.csv"),
                               col_types = col_types),
        ezdf = readr::read_csv(fs::path(arez_dir, "EZH2relatedGenesv08.csv"),
                               col_types = col_types),
        ebdf = utils::read.delim(fs::path(eb_dir, "ER_binding_gene_hg19_threshold_1e-5.txt"),
                                 stringsAsFactors = FALSE), # row names can't be parsed by read_tsv
        ebdf_t1 = readr::read_tsv(fs::path(nki_dir, "Tier_1_EntrezID_and_Symbol_ER_at_promotors.txt"),
                                  col_types = col_types),
        ebdf_t2 = readr::read_tsv(fs::path(nki_dir, "Tier_2_EntrezID_and_Symbol_total_ER_binding_via_20kb.txt"),
                                  col_types = col_types),
        icdf = readr::read_tsv(fs::path(exp_dir, "TCGAIntClust1100.txt"), col_types = col_types)
      )
      # Assign objects to current environment
      assign(., obj, pos = pos)
    }
  })
}

#' Process TCGA data
#'
#' Process TCGA data sources to produce outputs and results ready for analysis.
#'
#' The data sources loaded from `tcga_load` are processed so that we may use the
#' outputs and results for age-associated biomarker analysis.
#'
#' @param cldf clinical data
#' @param ealldf expression data
#' @param annodf annotation data
#' @param ardf age related genes
#' @param ezdf EZH2 related genes
#' @param ebdf ER binding genes from ChIP-Seq
#' @param ebdf_t1 NKI Tier 1 ER binding genes
#' @param ebdf_t2 NKI Tier 2 ER binding genes
#' @param icdf IntClust groups (leave `NULL` unless for breast cancer)
#' @param gender Either "all", "male", or "female" subsets
#'
#' @return A list with elements
#' * `exandf`: expression data restricted to annotation data
#' * `sedf`: data for regression analysis
#' * `erbnms`: ER binding gene names
#' * `arnms`: Age-related gene names
#' * `tcgaarnms`: TCGA gene names
#' @export
tcga_process <- function(cldf, ealldf, annodf, ardf, ezdf, ebdf,
                         ebdf_t1, ebdf_t2 ,icdf = NULL,
                         gender = c("all", "male", "female")) {
  # Choose gender subset
  gender <- toupper(match.arg(gender))
  if (gender == "ALL") gender <- c("MALE", "FEMALE")

  # Add variables to annotation data
  annodf <- annodf %>%
    dplyr::mutate(
      !!"Entrez_Gene_Id" := .data$Entrez_Gene_ID_0,
      !!"GL_split" := strsplit(.data$Genomic_location, split = ":"),
      !!"Chrcr" := vapply(.data$GL_split, "[", character(1), 1),
      !!"Chrc" := gsub("_qbl_hap2|_cox_hap1|_h2_hap1|_random", "", .data$Chrcr),
      !!"Chrn" := as.numeric(stringr::str_replace_all(.data$Chrc, c(
        "chr" = "", "X" = "23", "Y" = "24"
      ))),
      !!"CHR" := .data$Chrn,
      !!"Startcr" := vapply(.data$GL_split, "[", character(1), 2),
      !!"Stopcr" := vapply(.data$GL_split, "[", character(1), 3),
      !!"Startn" := as.numeric(.data$Startcr),
      !!"Stopn" := as.numeric(.data$Stopcr),
      !!"BP" := trunc((.data$Startn + .data$Stopn) / 2),
      !!"SNP" := .data$Entrez_Gene_Id
    ) %>%
    dplyr::select(-.data$GL_split)

  # Process IntClust variables
  if (!is.null(icdf)) {
    icdf <- icdf %>%
      dplyr::filter(!substring(.data$ID, 14, 15) == "06") %>%
      dplyr::mutate(!!"tcgaid" := paste0(substring(.data$ID, 1, 12), "-01"))
  }

  # Many entries 0, jitter medians so models are not singular, add small noise
  annoColNames <- c("Hugo_Symbol", "Entrez_Gene_Id")
  dataColidxs <- setdiff(names(ealldf), annoColNames)
  neDataRows <- nrow(ealldf)
  neDataCols <- length(dataColidxs)
  set.seed(719)
  jitterMat <- matrix(runif(neDataRows * neDataCols, min = 1e-06, max = 1e-02),
                      nrow = neDataRows, ncol = neDataCols)
  ealldf[, dataColidxs] <- ealldf[, dataColidxs] + jitterMat

  # Drop probesets with no variation
  ealldfSDs <- apply(log2(16 + ealldf[, dataColidxs]), 1, sd, na.rm = TRUE)
  edf <- ealldf[ealldfSDs > 1e-6, ]

  # Merge edf and annodf, drop unknown chromosome and chromosome 24 entries
  eandf <- edf[, annoColNames] %>%
    cbind(annodf[match(.$Entrez_Gene_Id, annodf$Entrez_Gene_Id),
                 c("Entrez_Gene_ID_0", "RefSeq_ID_0", "Accession_0",
                   "Gene_symbol", "SNP", "CHR", "BP")]) %>%
    dplyr::filter(.data$Entrez_Gene_Id %in% annodf$Entrez_Gene_Id,
                  !is.na(.data$CHR),
                  .data$CHR != 24) %>%
    dplyr::mutate(
      !!"TCGAsymbolEQMETABRICsymbol" := .data$Hugo_Symbol == .data$Gene_symbol,
      !!"refseqID" := sub("^(.*?)(\\..*|$)", "\\1", .data$RefSeq_ID_0)
    )

  # Get unique gene names from age-related and EZH2 gene lists
  arnms <- unique(unlist(strsplit(c(ardf$GeneName, ardf$IHCName), "/")))
  eznms <- unique(unlist(strsplit(c(ezdf$GeneName, ezdf$IHCName), "/")))
  areznms <- unique(unlist(lapply(list(arnms, eznms), toupper)))

  # TCGA age-related gene names
  tcgaarnms <- edf$Hugo_Symbol %>%
    `[`(. %in% areznms) %>%
    sort()

  # Get ER binding gene names
  erbnms <- tidy_genes(ebdf, "hgnc_symbol")
  erbnms_t1 <- tidy_genes(ebdf_t1, "Entrez ID")
  erbnms_t2 <- tidy_genes(ebdf_t2, "Entrez ID")

  # Match clinical ID to RNASeq IDs
  clKeepVars <- c("tcgaid", "AGE", "GENDER", "HISTORY_OTHER_MALIGNANCY",
                  "HISTORY_NEOADJUVANT_TRTYN", "RADIATION_TREATMENT_ADJUVANT",
                  "LYMPH_NODES_EXAMINED", "LYMPH_NODE_EXAMINED_COUNT",
                  "AJCC_TUMOR_PATHOLOGIC_PT", "AJCC_NODES_PATHOLOGIC_PN",
                  "AJCC_METASTASIS_PATHOLOGIC_PM", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                  "HISTOLOGICAL_DIAGNOSIS", "OS_STATUS", "OS_MONTHS")
  if (!is.null(icdf)) {
    clKeepVars <- c("tcgaid", "AGE", "GENDER", "MENOPAUSE_STATUS",
                    "ER_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS",
                    "HISTORY_OTHER_MALIGNANCY", "HISTORY_NEOADJUVANT_TRTYN",
                    "RADIATION_TREATMENT_ADJUVANT", "PHARMACEUTICAL_TX_ADJUVANT",
                    "SURGICAL_PROCEDURE_FIRST", "LYMPH_NODES_EXAMINED",
                    "LYMPH_NODE_EXAMINED_COUNT", "LYMPH_NODES_EXAMINED_HE_COUNT",
                    "AJCC_TUMOR_PATHOLOGIC_PT", "AJCC_NODES_PATHOLOGIC_PN",
                    "AJCC_METASTASIS_PATHOLOGIC_PM", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                    "METASTATIC_SITE", "HISTOLOGICAL_DIAGNOSIS", "OS_STATUS",
                    "OS_MONTHS")
  }

  # Filtered analysis data
  sdf <- cldf %>%
    dplyr::mutate(!!"tcgaid" := paste(.data$PATIENT_ID, "01", sep = "-")) %>%
    dplyr::select(dplyr::one_of(clKeepVars)) %>%
    dplyr::mutate(!!"Age" := purrr::map_dbl(
      .data$AGE, ~ ifelse(grepl("[[:alpha:]]", .), NA_real_, as.numeric(.))
    )) %>%
    dplyr::filter(!is.na(.data$Age) &
                    !is.na(match(.data$tcgaid, names(edf))) &
                    .data$GENDER %in% gender)
  if (!is.null(icdf)) {
    sdf <- dplyr::left_join(sdf, icdf, by = "tcgaid")
  }

  # Restrict expression data to annotated data in eandf
  exandf <- edf %>%
    dplyr::filter(.data$Entrez_Gene_Id %in% eandf$Entrez_Gene_Id) %>%
    dplyr::select(dplyr::one_of(setdiff(c(annoColNames, sdf$tcgaid),
                                        names(eandf)))) %>%
    cbind(eandf, .)

  # Data used for plotting and analysis
  sedf <- sdf %>%
    cbind(log2(16 + t(exandf[, match(.$tcgaid, names(exandf))]))) %>%
    purrr::set_names(c(names(sdf), exandf$Hugo_Symbol))
  if (!is.null(icdf)) {
    sedf <- sedf %>%
      dplyr::mutate(
        !!"BrCaEH" := dplyr::case_when(
          .data$ER_STATUS_BY_IHC == "Positive" &
            (.data$IHC_HER2 == "Negative" | (.data$IHC_HER2 == "Equivocal" & .data$HER2_FISH_STATUS == "Negative")) ~ "ER+/HER2-",
          .data$ER_STATUS_BY_IHC == "Positive" &
            (.data$IHC_HER2 == "Positive" | (.data$IHC_HER2 == "Equivocal" & .data$HER2_FISH_STATUS == "Positive")) ~ "ER+/HER2+",
          .data$ER_STATUS_BY_IHC == "Negative" &
            (.data$IHC_HER2 == "Positive" | (.data$IHC_HER2 == "Equivocal" & .data$HER2_FISH_STATUS == "Positive")) ~ "ER-/HER2+",
          .data$ER_STATUS_BY_IHC == "Negative" &
            (.data$IHC_HER2 == "Negative" | (.data$IHC_HER2 == "Equivocal" & .data$HER2_FISH_STATUS == "Negative")) ~ "ER-/HER2-",
          TRUE ~ NA_character_
        ),
        !!"BrCaEHf" := factor(.data$BrCaEH, levels = c("ER+/HER2-", "ER+/HER2+", "ER-/HER2+", "ER-/HER2-"))
      ) %>%
      magrittr::set_rownames(sdf$tcgaid)
  }

  # Assign objects to current environment
  dplyr::lst(exandf, sedf, erbnms, erbnms_t1, erbnms_t2, arnms, tcgaarnms) %>%
    purrr::iwalk(~ assign(.y, .x, pos = pos))
}

#' Create aroutdf object indicating p-values, estimates, age-dependent trend
#' indicators from the regression result
#'
#' @param exandf expression data
#' @param sedf analytical data
#' @param erbnms ER binding names
#' @param root.dir root directory of aroutdf object
#' @param age age threshold for calculating biological significance
#' @param fdr numeric vector of false discovery rates to obtain results for
#' @param fc numeric vector of fold change values to obtain results for
#' @param subset if `TRUE`, also save subsetted aroutdf results
#' @param overwrite if `TRUE`, save result even if it already exists
#' @param file.name core file name
#' @param suffix suffix for file name
#'
#' @return The object named aroutdf, with columns "Hugo_Symbol",
#'   "Entrez_Gene_Id", and the statistical results
#' @export
tcga_aroutdf <- function(exandf, sedf, erbnms, root.dir, age = 60,
                         fdr = 0.01, fc = 1.25, subset = FALSE,
                         overwrite = FALSE, file.name = NULL, suffix = NULL) {
  file.name <- file.name %||% "aroutdf"
  path_args <- list(root.dir, fdr, fc)
  df.path <- purrr::invoke(
    tcga_path, path_args, file.name = file.name, suffix = suffix)
  if (file.exists(df.path) & !overwrite) {
    aroutdf <- readr::read_csv(df.path)
  } else {
    # Set biological significance and multiple comparison p-value constants
    age.max <- max(sedf$Age)
    age.min <- min(sedf$Age)
    biosigslopeLE <- log2(fc) / (age - age.min + 1)
    biosigslopeGT <- log2(fc) / (age.max - age)
    biosigslopeAllAges <- log2(fc) / (age.max - age.min + 1)
    multcompPval <- fdr / (2 * nrow(exandf))
    co_names <- paste0(c("LE", "GT"), age)

    # Get below, above, overall lm estimates and p-values, plus age
    # dependence indicator (takes ~ 3 mins to run)
    arlmdf <- data.frame(Age = sedf$Age)
    aroutdf <- vapply(seq_len(nrow(exandf)), function(x) {
      psi <- exandf$Hugo_Symbol[x]
      arlmdf$probesetni <- sedf[, match(psi, names(sedf))]
      lm.all <- with(arlmdf, lapply(list(Age <= age, Age > age, NULL), function(s)
        lm(probesetni ~ Age, data = arlmdf, na.action = na.exclude, subset = s))) %>%
        lapply(broom::tidy) %>%
        dplyr::bind_rows() %>%
        magrittr::extract(.$term %in% "Age", c("estimate", "p.value")) %>%
        purrr::flatten_dbl() %>%
        c((abs(.[1]) > biosigslopeLE ||
             abs(.[2]) > biosigslopeGT ||
             abs(.[3]) > biosigslopeAllAges) &&
            (any(.[4:6] < multcompPval))) %>%
        purrr::set_names(
          purrr::cross2(c(co_names, "AllAges"), c("slope", "pval")) %>%
            purrr::map_chr(paste, collapse = "_") %>%
            c("AgeDependentp")
        )
      return(lm.all)
    }, numeric(7)) %>%
      t() %>%
      cbind(exandf[, seq_len(2)], .)

    # THIS ASSIGNS OBJECTS TO GLOBAL ENVIRONMENT
    # age_vars <- c("pval", "BHadj_pval", "slope", "log2FC", "Bestp")
    # sym_list <- purrr::cross2(paste0(c("LE", "GT"), age), age_vars) %>%
    #   purrr::map_chr(paste, collapse = "_") %>%
    #   purrr::set_names(gsub(age, "", .))
    # purrr::iwalk(sym_list, ~ assign(.y, rlang::sym(.x), pos = pos))

    # SPECIFY OBJECTS INDIVIDUALLY
    LE_pval <- sym_var(age, "LE", "pval")
    GT_pval <- sym_var(age, "GT", "pval")
    LE_BHadj_pval <- sym_var(age, "LE", "BHadj_pval")
    GT_BHadj_pval <- sym_var(age, "GT", "BHadj_pval")
    LE_slope <- sym_var(age, "LE", "slope")
    GT_slope <- sym_var(age, "GT", "slope")
    LE_log2FC <- sym_var(age, "LE", "log2FC")
    GT_log2FC <- sym_var(age, "GT", "log2FC")
    LE_Bestp <- sym_var(age, "LE", "Bestp")
    GT_Bestp <- sym_var(age, "GT", "Bestp")

    # Adjust all p-values and redo age-dependent selection on adjusted
    # p-values < fdr and find the largest significant fold change
    aroutdf <- aroutdf %>%
      dplyr::mutate(
        !!LE_BHadj_pval := p.adjust(!!LE_pval, method = "BH"),
        !!GT_BHadj_pval := p.adjust(!!GT_pval, method = "BH"),
        !!"AllAges_BHadj_pval" := p.adjust(.data$AllAges_pval, method = "BH"),
        !!"BHadj_and_AgeDependentp" := (abs(!!LE_slope) > biosigslopeLE & !!LE_BHadj_pval < fdr) |
          (abs(!!GT_slope) > biosigslopeGT & !!GT_BHadj_pval < fdr) |
          (abs(.data$AllAges_slope) > biosigslopeAllAges & .data$AllAges_BHadj_pval < fdr),
        !!"BHadj_signifp" := !!LE_BHadj_pval < fdr | !!GT_BHadj_pval < fdr | .data$AllAges_BHadj_pval < fdr,
        !!LE_log2FC := !!LE_slope * (age - age.min + 1),
        !!GT_log2FC := !!GT_slope * (age.max - age),
        !!"AllAges_log2FC" := .data$AllAges_slope * (age.max - age.min + 1),
        !!LE_Bestp := abs(!!LE_slope) > biosigslopeLE & !!LE_BHadj_pval < fdr,
        !!GT_Bestp := abs(!!GT_slope) > biosigslopeGT & !!GT_BHadj_pval < fdr &
          abs(!!GT_slope) > abs(!!LE_slope),
        !!"AllAges_Bestp" := abs(.data$AllAges_slope) > biosigslopeAllAges & .data$AllAges_BHadj_pval < fdr &
          (abs(.data$AllAges_log2FC) > abs(!!LE_log2FC) | abs(.data$AllAges_log2FC) > abs(!!GT_log2FC)),
        !!"Best_log2FC" := dplyr::case_when(
          .data$AllAges_Bestp ~ .data$AllAges_log2FC,
          !!GT_Bestp ~ !!GT_log2FC,
          !!LE_Bestp ~ !!LE_log2FC,
          TRUE ~ .data$AllAges_log2FC
        ),
        !!"Best_BHadj_pval" := dplyr::case_when(
          .data$AllAges_Bestp ~ .data$AllAges_BHadj_pval,
          !!GT_Bestp ~ !!GT_BHadj_pval,
          !!LE_Bestp ~ !!LE_BHadj_pval,
          TRUE ~ .data$AllAges_BHadj_pval
        ),
        !!"ERbinding" := .data$Hugo_Symbol %in% erbnms) %>%
      dplyr::select(-!!LE_Bestp, -!!GT_Bestp, -.data$AllAges_Bestp)

    # Write results to file
    readr::write_csv(aroutdf, path = df.path)
    # Write subsetted results
    if (subset)
      tcga_subset(aroutdf, path_args, file.name)
  }
  aroutdf
}

#' Final processing step to augment aroutdf with additional columns
#' @param erbnms ER binding gene names
#' @param erbnms_t1 NKI Tier 1 ER binding genes
#' @param erbnms_t2 NKI Tier 2 ER binding genes
#' @param arnms Age-related gene names
#' @inheritParams tcga_process
#' @inheritParams tcga_aroutdf
#' @export
tcga_augment <- function(root.dir, annodf, erbnms, erbnms_t1, erbnms_t2, arnms,
                         fdr = 0.01, fc = 1.25) {
  df.path <- tcga_path(root.dir, "aroutdf", fdr, fc, "tables/")
  aroutdf <- readr::read_csv(df.path)
  aroutdf <- aroutdf %>%
    dplyr::mutate(
      !!"PBHadj" := .data$Best_BHadj_pval,
      !!"PBHadj_log" := -log10(.data$PBHadj),
      !!"symbol_match" := match(.data$Entrez_Gene_Id, annodf$Original_Entrez),
      !!"gene_symbol" := annodf$Gene_symbol[.data$symbol_match],
      !!"entrez" := annodf$Entrez[.data$symbol_match],
      !!"ERbinding" := .data$gene_symbol %in% erbnms,
      !!"ERbinding1" := .data$entrez %in% erbnms_t1,
      !!"ERbinding2" := .data$entrez %in% erbnms_t2,
      !!"vcp" := abs(.data$Best_log2FC) > log2(fc) & .data$PBHadj < fdr,
      !!"binding_group" := dplyr::case_when(
        .data$vcp & .data$ERbinding ~ "ER binding",
        .data$vcp & !.data$ERbinding ~ "Non binding",
        TRUE ~ " "
      ),
      !!"binding_group_NKI" := dplyr::case_when(
        .data$vcp & .data$ERbinding1  & .data$ERbinding2 ~ "Tier 1",
        .data$vcp & !.data$ERbinding1 & .data$ERbinding2 ~ "Tier 2 Only",
        .data$vcp & !.data$ERbinding1 & !.data$ERbinding2 ~ "Non binding",
        TRUE ~ " "
      ),
      !!"evcp" := (abs(.data$Best_log2FC) > log2(fc) & -log10(.data$PBHadj) > 2) |
        -log10(.data$PBHadj) > 4,
      !!"location_sortable" := stringr::str_replace_all(
        annodf$Cytoband[.data$symbol_match], c("X" = "23", "Y" = "24")),
      !!"CHR" := stringr::str_pad(gsub("[_*A-z].+", "", .data$location_sortable),
                                  width = 2, side = "left", pad = "0") %>%
        as.numeric(),
      !!"Chrpq" := stringr::str_pad(gsub("(.*_*[A-z])|(.)|(\\.[0-9]*)", "\\2",
                                         .data$location_sortable),
                                    width = 2, side = "left", pad = "0") %>%
        ifelse(. == "00", NA, .),
      !!"BP" := as.numeric(.data$Chrpq) * 1e6,
      !!"SNP" := .data$gene_symbol,
      !!"arGeneSet_BHadj_and_AgeDependentp" := .data$BHadj_and_AgeDependentp &
        .data$gene_symbol %in% arnms)
  assign(paste("aroutdf", str_decimal(fdr), str_decimal(fc), sep = "_"),
         aroutdf, pos = pos)
}

#' Modify the biological significance conditions for the aroutdf object
#'
#' @inheritParams tcga_aroutdf
#' @param aroutdf age-related output data frame
#' @export
tcga_modify <- function(aroutdf, sedf, root.dir, age = 60, fdr = 0.01,
                        fc = 1.25, subset = FALSE, file.name = NULL,
                        suffix = NULL) {
  # Set biological significance and multiple comparison p-value constants
  file.name <- file.name %||% "aroutdf"
  path_args <- list(root.dir, fdr, fc)
  df.path <- purrr::invoke(tcga_path, path_args, file.name = file.name,
                           suffix = suffix)
  age.max <- max(sedf$Age)
  age.min <- min(sedf$Age)
  biosigslopeLE <- log2(fc) / (age - age.min + 1)
  biosigslopeGT <- log2(fc) / (age.max - age)
  biosigslopeAllAges <- log2(fc) / (age.max - age.min + 1)
  multcompPval <- fdr / (2 * nrow(aroutdf))

  LE_pval <- sym_var(age, "LE", "pval")
  GT_pval <- sym_var(age, "GT", "pval")
  LE_BHadj_pval <- sym_var(age, "LE", "BHadj_pval")
  GT_BHadj_pval <- sym_var(age, "GT", "BHadj_pval")
  LE_slope <- sym_var(age, "LE", "slope")
  GT_slope <- sym_var(age, "GT", "slope")

  # Modify conditions based on fdr and fc
  aroutdf <- aroutdf %>%
    dplyr::mutate(
      !!"AgeDependentp" := ifelse(
        (abs(!!LE_slope) > biosigslopeLE |
           abs(!!GT_slope) > biosigslopeGT |
           abs(.data$AllAges_slope) > biosigslopeAllAges) &
          pmin(!!LE_pval, !!GT_pval, .data$AllAges_pval) < multcompPval, 1L, 0L
      ),
      !!"BHadj_and_AgeDependentp" := (abs(!!LE_slope) > biosigslopeLE &
                                        !!LE_BHadj_pval < fdr) |
        (abs(!!GT_slope) > biosigslopeGT & !!GT_BHadj_pval < fdr) |
        (abs(.data$AllAges_slope) > biosigslopeAllAges &
           .data$AllAges_BHadj_pval < fdr),
      !!"BHadj_signifp" := !!LE_BHadj_pval < fdr |
        !!GT_BHadj_pval < fdr | .data$AllAges_BHadj_pval < fdr
    )
  # Write results to file
  readr::write_csv(aroutdf, path = df.path)
  # Write subsetted results
  if (subset)
    tcga_subset(aroutdf, path_args, file.name)
  aroutdf
}

#' @noRd
tcga_subset <- function(aroutdf, path_args, file.name) {
  # AgeDependent
  aroutdf_ad <- dplyr::filter(aroutdf, .data$AgeDependentp == 1)
  df.path.ad <- purrr::invoke(
    tcga_path,
    path_args,
    file.name = paste0("AgeDependent_", file.name))
  readr::write_csv(aroutdf_ad, path = df.path.ad)
  # BHadj_and_AgeDependentp
  aroutdf_ad_bhadj <- dplyr::filter(aroutdf, .data$BHadj_and_AgeDependentp)
  df.path.ad.bhadj <- purrr::invoke(
    tcga_path,
    path_args,
    file.name = paste0("BHadj_and_AgeDependent_", file.name))
  readr::write_csv(aroutdf_ad_bhadj, path = df.path.ad.bhadj)
}

#' Run Fisher's Exact Test on TCGA studies for specified FDR level
#' @inheritParams tcga_aroutdf
#' @param nki logical; if `TRUE` use NKI ER Binding gene list specifications
#' @export
tcga_fet <- function(root.dir, fdr = 0.01, fc = 1.25, file.name = NULL,
                     suffix = NULL, nki = FALSE) {
  file.name <- file.name %||% "aroutdf"
  if (nki) {
    elab <- c("Tier 1", "Tier 2 Only", "Not ER binding")
  } else {
    elab <- c("ER binding", "Not ER binding")
  }
  alab <- c("Age associated", "Not age associated")

  dplyr::lst(suffix, fc, fdr) %>%
    purrr::compact() %>%
    purrr::cross() %>%
    purrr::map(rev) %>%
    purrr::set_names(purrr::map(., paste, collapse = "_")) %>%
    purrr::map(~ {
      df.path <- tcga_path(root.dir, file.name, .$fdr, .$fc, suffix = .$suffix)
      fetdf <- readr::read_csv(file = df.path)
      fetmat <- fetdf %>%
        dplyr::transmute(
          `ER binding` = if (nki) {
            factor(ERbinding_NKI, elab)
          } else {
            factor(ifelse(ERbinding, elab[1], elab[2]), elab)
          },
          `Age association` = factor(ifelse(BHadj_and_AgeDependentp,
                                            alab[1], alab[2]), alab)
        ) %>%
        table() %>%
        tidy_ct()
    })
}

#' TCGA linear regression model scatterplots
#'
#' @inheritParams tcga_aroutdf
#' @param tcgaarnms TCGA age-related probe names
#' @param erbnms ER-binding probe names
#' @export
tcga_lm <- function(root.dir, file.name, sedf, tcgaarnms, erbnms, age = 60,
                    fdr = 0.01, fc = 1.25) {
  df.path <- tcga_path(root.dir, "aroutdf", fdr, fc, "tables/")
  aroutdf <- readr::read_csv(df.path)

  # Filter for genes in TCGA, add ER Binding Status
  LE_pval <- sym_var(age, "LE", "pval")
  GT_pval <- sym_var(age, "GT", "pval")

  aroutdf.lm <- aroutdf %>%
    dplyr::filter(.data$Hugo_Symbol %in% tcgaarnms) %>%
    dplyr::select(Gene = .data$Hugo_Symbol, !!LE_pval, !!GT_pval, .data$AllAges_pval) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::matches("pval")),
                     dplyr::funs(vapply(., format, digits = 5, character(1)))) %>%
    dplyr::mutate(
      !!LE_pval := paste0("[<=", age, "] p = ", !!LE_pval),
      !!GT_pval := paste0("[>", age, "] p = ",  !!GT_pval),
      !!"AllAges_pval" := paste0("[All] p = ", .data$AllAges_pval),
      !!"Trend" := ifelse(
        .data$Gene %in% with(aroutdf, Hugo_Symbol[BHadj_and_AgeDependentp]),
        paste0("Age-dependent trend:\n\n|FC|>", fc, " and adjPval<", fdr),
        paste0("No detectable trend:\n\n|FC|<", fc, " or adjPval>", fdr)
      ),
      !!"ERBS" := ifelse(.data$Gene %in% erbnms, "(ER binding)", "(non-ER binding)")
    )
  # Create plotting data: list for each gene with p-values and trend indicator
  lmdf <- sedf %>%
    dplyr::select(.data$Age, dplyr::one_of(tcgaarnms)) %>%
    tidyr::gather("Gene", "Expression", -1) %>%
    merge(aroutdf.lm) %>%
    split(.$Gene)

  # Apply function to list of all gene expression data
  all.plots <- purrr::map(lmdf, scatterfit_plot, conf.level = 1 - fdr)

  # All plots in one file
  grobs <- suppressWarnings(
    gridExtra::marrangeGrob(all.plots, nrow = 2, ncol = 2)
  )
  ggsave(tcga_path(
    root.dir, file.name, fdr, fc, "figures/", "loess", "pdf"),
    plot = grobs, width = 8, height = 10)

  # One plot per file
  suppressWarnings(
    purrr::iwalk(
      all.plots,
      ~ ggsave(tcga_path(root.dir, paste0(file.name, "_", .y), fdr, fc,
                         "figures/single_figures/", "loess", "pdf"),
               .x, width = 5, height = 6)
    )
  )
}

#' Read in a TCGA object from path specified by four parameters
#' @inheritParams tcga_aroutdf
#' @param sub.dir subdirectory
#' @param suffix suffix for full name, right before extension
#' @param extension file extension
#' @export
tcga_path <- function(root.dir, file.name, fdr, fc, sub.dir = NULL,
                      suffix = NULL, extension = "csv") {
  if (!is.null(suffix)) suffix <- paste0("_", str_sign(suffix))
  paste0(root.dir, sub.dir, file.name, "_FDR_", str_decimal(fdr),
         "_FC_", str_decimal(fc), suffix, ".", extension)
}
