#' Compare age association across cancer datasets
#'
#' Compare age associated status (presence/absence) between TCGA and METABRIC
#' datasets.
#'
#' @param metabric METABRIC dataset object
#' @param tcga list of TCGA dataset objects. Must have specific naming syntax.
#' @param fdr false discovery rate
#' @param fc fold change
#' @param organ TCGA organ determines how `tcga` is subsetted for the
#'   corresponding `fdr` and `fc` levels.
#'
#' @return An object of class "CrossTable" comparing the age association across
#'   different cancer datasets with a contingency table and different hypothesis
#'   tests of independence.
#' @export
comp_metabric_tcga <- function(metabric, tcga, fdr, fc, organ = "Br") {
  # METABRIC Entrez Gene IDs and age association absence/presence
  METABRIC_Entrez <- metabric %>%
    dplyr::transmute(
      Entrez = as.numeric(.data$Entrez_Gene_ID_0),
      METABRIC = ifelse(.data$Best_BHadj_pval < fdr &
                          .data$Abs_FoldChange > fc, 1, 0)
    )
  # TCGA Entrez Gene IDs and age association absence/presence
  TCGA_Entrez <- ageassn_tcga(tcga, organ, fdr, fc, "TCGA")
  # Relabel factor labels
  comp <- dplyr::full_join(METABRIC_Entrez, TCGA_Entrez, by = "Entrez") %>%
    tidyr::replace_na(list(TCGA = 0)) %>%
    dplyr::mutate_at(
      2:3,
      factor,
      levels = c("1", "0"),
      labels = c("Age Associated", "Not Age Associated")
    )
  # Contingency table and tidy cross tabulation results
  tab <- table(comp[[2]], comp[[3]],
               dnn = list("METABRIC", paste0("TCGA_", organ, "Ca")))
  suppressWarnings(tidy_ct(tab))
}

#' Compare pairs of TCGA data on age association
#'
#' @param x first list of TCGA data
#' @param y second list of TCGA other
#' @param x_organ first TCGA organ
#' @param y_organ second TCGA organ
#' @param N total number of cases
#' @inheritParams comp_metabric_tcga
#' @export
comp_pair_tcga <- function(x, y, x_organ, y_organ, fdr, fc, N = 18950) {
  # TCGA Entrez Gene IDs and age association absence/presence for both data
  TCGA_x <- ageassn_tcga(x, x_organ, fdr, fc, "TCGA_x")
  TCGA_y <- ageassn_tcga(y, y_organ, fdr, fc, "TCGA_y")
  dnn <- list(paste0("TCGA_", x_organ, "Ca"), paste0("TCGA_", y_organ, "Ca"))
  if (nrow(TCGA_x) == 0 & nrow(TCGA_y) == 0) {
    # In case both x and y have no cases
    naa <- factor(rep("Not Age Associated", N),
                  levels = c("Age Associated", "Not Age Associated"))
    tab <- table(naa, naa, dnn = dnn)
  } else {
    # Relabel factor labels
    comp <- dplyr::full_join(TCGA_x, TCGA_y, by = "Entrez") %>%
      tidyr::replace_na(list(TCGA_x = 0, TCGA_y = 0)) %>%
      dplyr::mutate_at(
        2:3,
        factor,
        levels = c("1", "0"),
        labels = c("Age Associated", "Not Age Associated")
      )
    # Contingency table and tidy cross tabulation results
    tab <- table(comp[[2]], comp[[3]], dnn = dnn)
    if (sum(tab) != N) tab[2, 2] <- N - sum(tab)
  }
  suppressWarnings(tidy_ct(tab))
}

#' Depending on `organ`, extract a tibble with Entrez Gene ID and age
#' association identifier columns
#' @inheritParams comp_metabric_tcga
#' @param .id column name of age association identifier
#' @noRd
ageassn_tcga <- function(tcga, organ, fdr, fc, .id = NULL) {
  # Set identifier name
  id <- .id %||% "id"
  if (organ == "Br") {
    # TCGA BrCa results needs to be summarized
    tcga[[as.character(fdr)]] %>%
      dplyr::transmute(
        Entrez = .data$Entrez_Gene_Id,
        !!rlang::sym(id) := ifelse(.data$PBHadj < fdr &
                                     .data$Abs_FoldChange > fc, 1, 0)
      )
  } else {
    # TCGA Other results already filtered, so all are age associated
    tcga[[paste(fdr, fc)]] %>%
      dplyr::transmute(
        Entrez = as.numeric(.data$Entrez_Gene_Id),
        !!rlang::sym(id) := 1
      )
  }
}
