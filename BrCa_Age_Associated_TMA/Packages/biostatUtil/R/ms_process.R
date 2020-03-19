#' Process mass spectrometry data
#'
#' Process mass spectrometry data for filtering out contaminant samples,
#' manipulating variables, removing duplicates, and more.
#'
#' @param psm PSM file
#' @param protein Protein file
#' @param treatment character vector of treatment groups
#' @param samples character vector of sample names to keep
#' @param sample.id character vector for sample IDs. Order of samples must match
#'   that in the psm raw data.
#' @param path file path to save return element `pep`
#' @param ... additional arguments to `ms_condition`
#' @return A list with the following elements
#' \item{pep}{processed data frame to be used by `ms_summarize`}
#' \item{raw}{raw data values}
#' \item{l2}{log2 raw data values}
#' \item{vsn}{vsn raw data values}
#' @family Mass Spectrometry functions
#' @author Derek Chiu
#' @export
ms_process <- function(psm, protein, treatment, samples = NULL,
                       sample.id = NULL, path = NULL, ...) {
  # Make raw data file column names into R column names
  psm <- psm %>% magrittr::set_colnames(make.names(colnames(.)))
  protein <- protein %>% magrittr::set_colnames(make.names(colnames(.)))

  # Variables to keep
  samples <- samples %||% grep("^X[0-9]", names(psm), value = TRUE)
  if (is.null(sample.id)) {
    ns <- seq_along(samples)
    sample.id <- paste0("X", ceiling(ns * 2 / max(ns)),
                        "_", seq_len(max(ns) / 2))
    assertthat::assert_that(all(purrr::map_lgl(treatment,
                                               ~ any(grepl(.x, sample.id)))))
  }
  psmKeepVars <-
    c("Annotated.Sequence", "Modifications", "Number.of.Protein.Groups",
      "Number.of.Proteins", "Master.Protein.Accessions", "Protein.Accessions",
      "Confidence", "Reporter.Quan.Result.ID", "Quan.Info", "PSM.Ambiguity",
      samples)

  # Select relevant columns for analysis
  # Rename original sample labels to reflect sample grouping
  # Filter out:
  #   - non-unique peptides
  #   - no quan values
  #   - insufficient data: use ms_condition()
  #   - Master Protein Accession missing or "sp"
  pep <- psm %>%
    select(one_of(psmKeepVars)) %>%
    rename_(.dots = stats::setNames(samples, sample.id)) %>%
    filter_(.dots = list(lazyeval::interp(
      ~NOP == 1 & !grepl("NoQuanValues", QI) & (is.na(MPA) | MPA != "sp"),
      NOP = quote(Number.of.Proteins), QI = quote(Quan.Info),
      MPA = quote(Master.Protein.Accessions)
    ))) %>%
    magrittr::extract(ms_condition(., treatment = treatment, ...), )

  pro <- protein %>%
    select_("Accession", "Description", "MW.in.kDa")

  # For each Reporter.Quan.Result.ID in pep, remove duplicates for 4 vars
  pep <- pep %>%
    group_by_("Reporter.Quan.Result.ID") %>%
    do(remove_dup(., c("Annotated.Sequence", "Modifications",
                       "Master.Protein.Accessions", "Protein.Accessions")))

  # Parse peptide accession and merge the protein descriptions with the peptide file
  pep <- pep %>%
    mutate_(.dots = stats::setNames(list(lazyeval::interp(
      ~purrr::map_chr(strsplit(as.character(MPA), ";"), `[`, 1),  # Strip ";"
      MPA = quote(Master.Protein.Accessions))), "Accession")) %>%
    mutate_(.dots = stats::setNames(list(lazyeval::interp(
      ~purrr::map_chr(strsplit(as.character(A), " | "), `[`, 1),  # Strip " | "
      A = quote(Accession))), "Accession")) %>%
    merge(pro, by = "Accession") %>%  # Merge with protein set on Accession
    mutate_(.dots = stats::setNames(list(
      lazyeval::interp(~sub(".*?GN=(.*?)( .*|$)", "\\1", D),
                       D = quote(Description)),  # Get the gene name out
      lazyeval::interp(~toupper(sub(".*?\\.(.*?)(\\..*|$)", "\\1", AS)),
                       AS = quote(Annotated.Sequence)),   # Parse the peptide column for amino acids
      lazyeval::interp(~sub("(.*?)( OS=.*|$)", "\\1", D),
                       D = quote(Description))),  # Filter information from Description
      c("Gene", "Sequence", "Descriptions"))) %>%
    filter_(.dots = list(lazyeval::interp(
      ~!grepl("Keratin", quote(D)) &
        !grepl("sp", quote(A), ignore.case = FALSE) &
        !grepl("ribosomal", quote(D)),  # Remove specific proteins (typically contaminants from other sources)
      D = quote(Descriptions), A = quote(Accession))))
  if (!is.null(path))
    readr::write_csv(pep, path = path)

  # Raw, log2, and vsn transformed expression data
  raw <- pep[, sample.id]
  l2 <- log2(raw) %>%
    magrittr::set_colnames(paste("l2", names(.), sep = "_"))
  vsn <- limma::normalizeVSN(raw, verbose = FALSE) %>%
    magrittr::set_colnames(paste("vsn", colnames(.), sep = "_"))
  return(list(pep = pep, raw = raw, l2 = l2, vsn = vsn))
}
