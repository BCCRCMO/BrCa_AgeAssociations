#' Manhattan Plot
#'
#' @param x data object
#' @param CHR chromosome numbers
#' @param BP base pairs
#' @param PBHadj BH-adjusted p-values
#' @param snp SNP gene names
#' @param logp logical; if `TRUE` use -log10 p-values instead of raw
#' @param chrlabs labels for chromosomes
#' @param col colours to use for chromosomes; can be recycled
#' @param suggestiveline line for suggestive association
#' @param genomewideline line for genomewide association
#' @param highlight logical; if `TRUE`, points annotated by shape and label
#' @param highlight.var variable name to use highlighting
#' @param plotpointhiliteidxp whether to highlight points by shape
#' @param hilitelbls labels to use for highlighting
#' @param ... additional parameters to `plot`
#'
#' @return A manhattan plot of chromosome number vs adjusted p-value
#' @author Derek Chiu
#' @references https://github.com/stephenturner/qqman/blob/master/R/manhattan.R
#' @export
manhattan_plot <- function(x, CHR = "CHR", BP = "BP", PBHadj = "PBHadj",
                           snp = "SNP", logp = TRUE, chrlabs = NULL,
                           col = c("black", "red3", "green4", "blue",
                                   "cyan3", "magenta3", "orange", "gray40"),
                           suggestiveline = 0.05, genomewideline = FALSE,
                           highlight = FALSE, highlight.var = NULL,
                           plotpointhiliteidxp = FALSE, hilitelbls = NULL,
                           ...) {
  # Ensure columns exist in data
  assertthat::assert_that(CHR %in% names(x), BP %in% names(x),
                          PBHadj %in% names(x), snp %in% names(x))
  assertthat::assert_that(is.numeric(x[[CHR]]), is.numeric(x[[BP]]),
                          is.numeric(x[[PBHadj]]))

  # Plot ER binding status as different shapes
  if (plotpointhiliteidxp) {
    pch.idx <- x$ERbinding
  } else {
    pch.idx <- rep(FALSE, nrow(x))
  }

  # Variable to use for labelling
  if (highlight) {
    if (is.null(highlight.var)) {
      highlight.var <- "BHadj_and_AgeDependentp"
    }
    gene.hilite <- x[[snp]][x[[highlight.var]]]
  }

  # Choose/add variables, filter for complete cases, and set p-value type
  d <- x %>%
    dplyr::select(CHR, BP, PBHadj) %>%
    dplyr::mutate(
      ppch = ifelse(pch.idx, 8, 20),
      SNP = ifelse(rep(!is.null(x[[snp]]), nrow(.)), x[[snp]], NA)
    ) %>%
    dplyr::select_if(~ !inherits(., "logical")) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::mutate(
      logp = ifelse(rep(logp, nrow(.)), -log10(PBHadj), PBHadj),
      pos = NA,
      index = CHR
    )

  # For single chromosome and general case
  nchr <- length(unique(d$CHR))
  if (nchr == 1) {
    options(scipen = 999)
    d$pos <- d$BP / 1e6
    ticks <- floor(length(d$pos)) / 2 + 1
    xlabel <- paste("Chromosome", unique(d$CHR), "position(Mb)")
    labs <- ticks
  } else {
    lastbase <-
      with(d, purrr::map_dbl(unique(index), ~ dplyr::last(BP[index == . - 1]))) %>%
      ifelse(is.na(.), 0, .) %>%
      cumsum()
    d$pos <-
      with(d, purrr::map(unique(index), ~ BP[index == .] + lastbase[.])) %>%
      unlist()
    ticks <-
      with(d, purrr::map_dbl(unique(index), ~ mean(range(pos[CHR == .])) + 1))
    xlabel <- "Chromosome"
    labs <- unique(d$CHR)
  }

  # Plotting parameters
  xmax <- ceiling(max(d$pos) * 1.03)
  xmin <- floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax),
                   ylim = c(0, 1.1 * ceiling(max(d$logp))),
                   xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                            names(dotargs)]))

  # Chromosome labels with appropriate warnings
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }

  # Axes for single chr and general case (iteratively to avoid skipping)
  if (nchr == 1) {
    axis(1, ...)
  } else {
    for (i in seq_along(ticks))
      axis(1, at = ticks[i], labels = labs[i], las = 2, ...)
  }

  # Plotting points for single chromosome and general case
  ccol = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = ccol[1], ...))

  } else {
    for (i in unique(d$index)) {
      with(d[d$index == i, ], points(pos, logp, col = ccol[i], pch = ppch, ...))
      with(d[d$index == i, ], abline(v = max(pos), lty = "dashed", col = "grey80"))
    }
  }

  # Add lines of interest
  if (suggestiveline) {
    if (logp) {
      abline(h = -log10(suggestiveline), col = "blue")
    } else {
      abline(h = suggestiveline, col = "blue")
    }
  }
  if (genomewideline)  {
    if (logp) {
      abline(h = -log10(genomewideline), col = "red")
    } else {
      abline(h = genomewideline, col = "blue")
    }
  }

  # Add gene labelling
  if (highlight) {
    if (length(gene.hilite) == 0) {
      warning("There are no significant SNPs to highlight in your results.")
    } else if (all(!gene.hilite %in% d$SNP)) {
      warning("You're trying to highlight SNPs that don't exist in your results.")
    } else {
      d.highlight = d[which(d$SNP %in% gene.hilite), ]
      with(d.highlight, maptools::pointLabel(pos, logp, SNP, cex = 0.8))
    }
    if (plotpointhiliteidxp && !is.null(hilitelbls) && any(pch.idx)) {
      legend("topleft", legend = hilitelbls, pch = c(8, 20), lwd = 2, lty = 0)
    }
  }
}

#' Make a manhattan plot specifically for data structures in TCGA studies
#'
#' @param data list of data objects for plotting
#' @param root.dir root directory to save plots
#' @param plot.suffix suffix at end of plot file name
#' @param title plot title
#' @export
tcga_manhattan <- function(data, root.dir, plot.suffix, title) {
  # store file names
  filenames <- data %>%
    purrr::imap(~ {
      fn_split <- stringr::str_split_fixed(.y, "_", 3)
      tcga_path(root.dir = root.dir, file.name = "AgeRelated_Manhattan",
                fdr = fn_split[, 2], fc = fn_split[, 3], sub.dir = "figures/",
                suffix = plot.suffix, extension = "pdf")
    })
  # store suggestive lines
  sl <- data %>%
    purrr::imap_chr(~ stringr::str_split_fixed(.y, "_", 3)[, 2]) %>%
    purrr::map_chr(~ gsub("p", "\\.", .)) %>%
    as.numeric()
  # common arguments to invoke plotting function over
  args <- list(chrlabs = c(1:22, "X"),
               ylab = "-log10(Benjamini-Hochberg adjusted P-values)",
               main = title, highlight = TRUE, plotpointhiliteidxp = TRUE,
               hilitelbls = c("ER binding", "No binding"))
  # invoke plotting function on each of the objects, saving to different files
  purrr::pwalk(list(filenames, data, sl), ~ {
    pdf(..1)
    purrr::invoke(manhattan_plot, args, x = ..2, suggestiveline = ..3)
    dev.off()
  })
}
