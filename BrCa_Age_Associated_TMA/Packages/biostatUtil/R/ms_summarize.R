#' Summarize mass spectrometry data
#'
#' Summarize mass spectrometry data for expression differences between sample
#' groups. Apply Benjamini-Hochberg p-value adjustment.
#'
#' Gene-level analysis is performed on the "Gene" variable. Peptide-level
#' analysis is performed on distinct combinations of the Accession, Gene,
#' Descriptions, Sequence, and Modifications variables. This combined variable
#' is renamed to "AGDSM".
#'
#' @param x data frame with expression values
#' @param g vector of factor levels for samples
#' @param level analysis is on the "Gene" level or "Peptide" level
#' @param col.names vector of column names for output data frame
#' @param info.vars vector of column names containing metadata information.
#'   These variables are collapsed if not unique.
#' @param path file path to save result object
#' @return A data frame of statistics from analyzing mass spec data. Includes
#'   t-values, Wald p-values, effect sizes, fold change, absolute fold change.
#' @family Mass Spectrometry functions
#' @author Derek Chiu
#' @export
ms_summarize <- function(x, g, level = c("Gene", "Peptide"), col.names = NULL,
                         info.vars = NULL, path = NULL) {

  # Determine level of variable split
  level <- match.arg(level)
  var.split <- switch(level, Gene = "Gene", Peptide = "AGDSM")

  # NULL defaults
  if (is.null(info.vars))
    info.vars <- c(var.split, "Accession", "Sequence", "Annotated.Sequence",
                   "Descriptions", "Modifications", "Reporter.Quan.Result.ID")

  # Combine values, add/select vars, per-level analyses, numeric coercion
  expr <- x %>%
    magrittr::extract(c("pep", "l2", "vsn")) %>%
    unname() %>%
    purrr::invoke(cbind, .) %>%
    mutate_(.dots = stats::setNames(list(
      lazyeval::interp(~paste(A, G, D, sep = " || "),
                       A = quote(Accession), G = quote(Gene),
                       D = quote(Description)),
      lazyeval::interp(~paste(A, G, D, S, M, sep = " || "),
                       A = quote(Accession), G = quote(Gene),
                       D = quote(Descriptions), S = quote(Sequence),
                       M = quote(Modifications)),
      lazyeval::interp(~as.character(R), R = quote(Reporter.Quan.Result.ID))),
      c("AGD", "AGDSM", "Block"))) %>%
    select(one_of(c(info.vars, "Block", colnames(x$vsn))))
  res <- purrr::map(sort(unique(expr[[var.split]])), ~
                      ms_analyze(expr[expr[[var.split]] == .x, ], g = g,
                                 level = level, col.names = col.names,
                                 info.vars = info.vars[-1])) %>%
    purrr::invoke(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate_at(.cols = vars(-one_of(info.vars), matches("dir|adj")),
              .funs = as.numeric)

  # Data only with adjusted p-values (BH), synced column names
  adj <- res %>%
    mutate_at(.cols = vars(matches("p-*val"), -matches("adj")),
              .funs = funs(Temp = stats::p.adjust(., method = "BH"))) %>%
    select(contains("Temp")) %>%
    magrittr::set_names(grep("adj", col.names, value = TRUE))

  # Replace placeholder columns using new columns with adjusted p-values
  adj.ind <- match(names(adj), names(res))
  all.ind <- order(match(c(names(res), names(adj)), col.names))
  final <- cbind(res, adj)[, all.ind[!(all.ind %in% adj.ind)]]

  if (!is.null(path))
    readr::write_csv(final, path)
  return(final)
}

#' Analyze mass spectrometry data for expression differences between sample
#' groups.
#' @noRd
ms_analyze <- function(x, g, level, col.names, info.vars) {
  # Create factor variable for different treatments in gene-specific data frame
  # Modify `mutate_()` call depending on variable coding
  adf <- tidyr::gather(data = x, key = "Sample", value = "tInt",
                       matches("vsn")) %>%
    mutate_(.dots = stats::setNames(list(lazyeval::interp(
      ~factor(gsub("vsn_(.*)_.+", "\\1", var), levels = g),
      var = quote(Sample))), "Trtf"))

  # Test interaction model if more than one block, otherwise test null model
  if (length(unique(adf$Block)) > 1) {
    mod1 <- stats::lm(tInt ~ Block, data = adf)
    mod2 <- stats::lm(tInt ~ Block * Trtf, data = adf)
  } else {
    mod1 <- stats::lm(tInt ~ 1, data = adf)
    mod2 <- stats::lm(tInt ~ Trtf, data = adf)
  }

  # Variable level
  Level <- switch(level,
                  Gene = unique(adf$Gene),
                  Peptide = unique(adf$AGDSM))

  # Extract F value, df (num), df (den), p-value from Omnibus test
  OmnibusTrt <- stats::anova(mod1, mod2)[2, c("F", "Df", "Res.Df", "Pr(>F)")]
  Omnibus_obj <- c(as.character(OmnibusTrt), "BHadj_OmnibusTrtPvalHere")

  # Fitted means and standard error
  # Set up contrasts: trt effect = 1, otherwise = 0
  Trtf_idx <- grepl(paste0("Trtf", g[2]), names(stats::coef(mod2)))
  ate_contr <- ifelse(Trtf_idx, 1, 0)

  # Transform any interaction coefficient contrasts to 1 / # of Trt coefs
  if (any(Trtf_idx)) {
    Intn_idx <- grepl(paste0(":Trtf", g[2]), names(stats::coef(mod2)))
    ate_contr[Intn_idx] <- 1 / sum(Trtf_idx)
  }

  coefs <- stats::coef(mod2)
  Intcp <- coefs[grep("Intercept", names(coefs))]
  Trtmt <- Intcp + coefs[grep("Intercept|Block", names(coefs), invert = TRUE)]
  Stdev <- sqrt(ate_contr %*% stats::vcov(mod2) %*% ate_contr)
  MeanVar_obj <- unname(c(Intcp, Trtmt, Stdev))

  # Compute statistics for treatment levels (omit first group: reference)
  Trt_stats <- unlist(lapply(g[-1], treatment_stats, data = adf, mod = mod2))

  # Descriptive information
  Desc_obj <- adf %>%
    select(one_of(info.vars)) %>%
    mutate_all(funs(collapse_var)) %>%
    unique() %>%
    as.character()

  # Combine data objects and set column names
  return(stats::setNames(c(Level, Omnibus_obj, MeanVar_obj, Trt_stats, Desc_obj),
                         col.names))
}

#' Compute statistics for treatment-specific comparisons
#' @param trt treatment level of grouping variable
#' @param data data frame for per-level analysis
#' @param mod model object; usually full or larger model
#' @param alpha significance level; defaults to 0.05
#' @noRd
treatment_stats <- function(trt, data, mod, alpha = 0.05) {

  # Set up contrasts: trt effect = 1, otherwise = 0
  Trtf_idx <- grepl(paste0("Trtf", trt), names(stats::coef(mod)))
  ate_contr <- ifelse(Trtf_idx, 1, 0)

  # Transform any interaction coefficient contrasts to 1 / # of Trt coefs
  if (any(Trtf_idx)) {
    Intn_idx <- grepl(paste0(":Trtf", trt), names(stats::coef(mod)))
    ate_contr[Intn_idx] <- 1 / sum(Trtf_idx)
  }

  # Effect, means, and standard error (from covariance matrix)
  Effect <- ate_contr %*% stats::coef(mod)
  Stdev <- sqrt(ate_contr %*% stats::vcov(mod) %*% ate_contr)

  # Compute t value, (adj) p-values, (log2) effect, (absolute) fold change
  tcrit <- stats::qt(1 - alpha / 2, df = mod$df.residual, lower.tail = TRUE)
  tval <- Effect / Stdev
  Waldpval <- 2 * stats::pt(abs(tval), df = mod$df.residual, lower.tail = FALSE)
  Effect_obj <- c(tval, Waldpval, "BHadj_WaldPvalHere")
  Log2Effect_obj <- Effect + c(-1:1) * tcrit * Stdev
  FC_obj <- 2 ^ Log2Effect_obj
  Effect_dir <- rep(Effect > 0, 4)
  AbsFC_obj <- ifelse(Effect_dir, c("Up", FC_obj), c("Down", 1 / rev(FC_obj)))
  return(c(Effect_obj, Log2Effect_obj, FC_obj, AbsFC_obj))
}
