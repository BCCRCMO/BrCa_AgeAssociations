#' Volcano Plot
#'
#' @param data data object with regression fit result variables
#' @param x log2 fold change
#' @param y adjusted p-value (preferrably in -log10 scale)
#' @param group grouping variable of interest (e.g. ER binding group)
#' @param label gene variable for labelling points
#' @param subset subset of `data` used for labelling points
#' @param legend character vector for legend key items
#' @param color character vector for colour of each legend key item
#' @param title title of plot
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param seed reproducible text repel direction
#' @param filename file path to save figure to (PDF)
#' @param style use [graphics::smoothScatter()] or [ggplot2::ggplot()] plotting
#'   device?
#' @param thin logical; if `TRUE`, a clustering is applied on the non
#'   age-dependent points and only certain low/high density points are plotted.
#'   Smoothing with [graphics::smoothScatter()] is not used.
#'
#' @return a ggplot object storing the volcano plot
#' @author Derek Chiu
#' @export
volcano_plot <- function(data, x, y, group, label, subset,
                         legend = NULL, color = rainbow(3), title = "",
                         xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                         seed = 1, filename = NULL, style = c("ss", "gg"),
                         thin = FALSE, legloc="bottomleft") {
  # Choose plotting style
  style <- match.arg(style)
  # Set axis label and limit defaults
  xlab <- xlab %||% "(Down)        FC        (Up)   "
  ylab <- ylab %||% "-log10(Benjamini-Hochberg adjusted P-values)"
  xlim <- xlim %||% max(abs(range(data[[x]], na.rm = TRUE))) * c(-1, 1) * 1.1
  ylim <- ylim %||% c(0, 1.1 * max(data[[y]]))
  gglegloc <- switch(legloc,
                     "bottomleft"=c(0, 0),
                     "topleft"=c(0, 1),
                     "topright"=c(1, 1),
                     "bottomright"=c(1, 0),
                     "bottom"=c(0.5, 0),
                     "left"=c(0, 0.5),
                     "top"=c(0.5, 1),
                     "right"=c(1, 0.5),
                     "center"=c(0.5, 0.5),
                     c(0, 0)
                     )

  if (style == "ss") {
    if (!is.null(filename)) {
      pdf(filename)
    }
    xy_args <- dplyr::lst(xaxt = "n", xlab, ylab, xlim, ylim)
    if (thin) {
      # Thin points to avoid smudges when zooming in on plot
      args <- c(xy_args, pch = 19, col = col2rgba("black"))
      data %>%
        thin_plot(x, y, group, subset) %>%
        with(purrr::invoke(plot, args, x = .[[x]], y = .[[y]]))
    } else {
      # Scatterplot for all points without ER binding status
      args <- c(xy_args,
                colramp = colorRampPalette(c("white", col2rgba("black"))),
                nrpoints = 0)
      data %>%
        dplyr::filter(!.[[subset]]) %>%
        with(purrr::invoke(smoothScatter, args, x = .[[x]], y = .[[y]]))
    }
    # Add coloured points with ER binding status
    purrr::walk2(legend, color, ~ {
      data %>%
        dplyr::filter(.data[[group]] == .x) %>%
        with(points(x = .[[x]], y = .[[y]], pch = 19, col = col2rgba(.y)))
    })
    # Add back title, axis tick marks, legend box
    title(title, line = 1)
    axis(side = 1, at = -5:5, labels = 2 ^ abs(-5:5))
    legend(x=legloc, legend = legend, col = color, lwd = 3, pch = 1)
    # Add annotated symbol names for genes filtered in subset
    if (any(data[[group]][data[[subset]]] != " ")) {
      data %>%
        dplyr::filter(.[[subset]] & .[[group]] != " ") %>%
        with(maptools::pointLabel(
          x = .[[x]], y = .[[y]], labels = .[[label]], cex = 0.6,
          col = dplyr::recode(.[[group]], !!!as.list(purrr::set_names(color, legend)))
        ))
    }
    if (!is.null(filename)) {
      dev.off()
    }
  } else if (style == "gg") {
    # Use seed to make sure the geom_text_repel directions don't change
    set.seed(seed)
    # By default, obtain legend key names from unique elements of grouping var
    legend <- legend %||% unique(data[[group]])
    # Volcano plot uses a positive-only x-axis label, 30% transparency
    vp <- ggplot(data, aes_string(x, y)) +
      geom_point(aes_string(color = group), size = 1, alpha = 0.3) +
      ggrepel::geom_text_repel(data = data[data[[subset]], ],
                               aes_string(x, y, label = label, color = group),
                               segment.colour = NA, show.legend = FALSE) +
      scale_color_manual(breaks = legend, values = c("grey70", color)) +
      scale_x_continuous(breaks = -3:3, labels = 2 ^ abs(-3:3), limits = xlim) +
      ylim(ylim) +
      labs(title = title, x = xlab, y = ylab) +
      theme_linedraw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = gglegloc,
            legend.justification = c(0, 0),
            legend.background = element_rect(colour = "black"),
            legend.title = element_blank())
    # If filename given, save vp
    if (!is.null(filename)) {
      ggsave(filename = filename, plot = vp)
    }
    return(vp)
  }
}

#' Thin a volcano plot
#'
#' Thin the plot by keeping some low density points and representative high
#' density points.
#'
#' @inheritParams volcano_plot
#' @param prop_lo proportion of low density points
#' @param prop_hi proportion of high density points
#' @param max_pts maximum number of thinned points to keep
#' @param max_size maximum size of data subsets to apply clustering on
#' @return filtered `data` keeping indices of thinned points
#' @noRd
thin_plot <- function(data, x, y, group, subset, prop_lo = 0.66,
                      prop_hi = 1 - prop_lo, max_pts = 1500, max_size = 15000) {
  # Split data into equal parts before clustering
  np <- sum(data[[group]] == " ")
  n_split <- ifelse(np > max_size, np %/% max_size, 1)
  split_data <- data %>%
    dplyr::filter(.[[group]] == " ") %>%
    dplyr::select(x, y) %>%
    dplyr::mutate_at(y, jitter) %>%
    split_n(n = n_split)
  # Low and high density indices
  density_idxs <- suppressWarnings(
    split_data %>%
      purrr::map(
        amap::hcluster, # HC: euclidean distance, single linkage on x & y
        link = "single",
        nbproc = parallel::detectCores()
      ) %>%
      purrr::map("merge") %>%
      purrr::map(~ {
        # Number of low and high density points based on clustering results
        nNotAgeDep <- sum(. < 0, na.rm = TRUE)
        nLowDens <- trunc(prop_lo * min(max_pts / n_split, nNotAgeDep))
        nHiDens <- trunc(prop_hi * min(max_pts / n_split, nNotAgeDep))
        # Pick out points to plot from both density regions
        singletons <- -1 * purrr::keep(., ~ . < 0)
        lowDensityidxs <- rev(singletons)[seq_len(nLowDens)]
        hiDensityidxs <-
          sample(setdiff(singletons, lowDensityidxs), size = nHiDens)
        c(lowDensityidxs, hiDensityidxs)
      }) %>%
      list_index(split_data, .) # Combine indices to keep in each data subset
  )
  # Filter data, including points from `subset`
  NotAgeDepidxs <- which(data[[group]] == " ")
  MBlblidxs <- purrr::keep(which(data[[subset]]), ~ . %in% NotAgeDepidxs)
  plotNotAgeDepidxs <- c(density_idxs, MBlblidxs)
  dplyr::slice(data, NotAgeDepidxs[plotNotAgeDepidxs])
}

#' Randomly split a data frame into n equal parts
#' @param data data frame
#' @param n number of splits
#' @noRd
split_n <- function(data, n = 3) {
  nr <- nrow(data)
  s <- nr %/% n
  split(data, sample(rep(1:trunc(nr / s), each = s, length.out = nr)))
}

#' Return unique vector of indices for a combined list of data frames
#' given indices for each individual data frame
#' @param .x list of data frames
#' @param idx list of indices
#' @noRd
list_index <- function(.x, idx) {
  rows <- purrr::map_int(.x, nrow)
  lasts <- cumsum(rows)
  firsts <- lasts - rows + 1
  all_idx <- purrr::map2(firsts, lasts, seq)
  purrr::flatten_int(purrr::map2(all_idx, idx, `[`))
}

#' Maximum axes limits
#'
#' Calculate maximum x-axis/y-axis limits across list of datasets for
#' plotting equal axes across a subgroup.
#'
#' @param .x list of datasets
#' @param var variable name
#' @name max_limits
#' @export
max_xlim <- function(.x, var = "Best_log2FC") {
  .x %>%
    purrr::map(var) %>%
    purrr::map_dbl(~ max(abs(range(., na.rm = TRUE)))) %>%
    max() %>%
    magrittr::multiply_by(1.1) %>%
    rep(2)  # +ve to +ve for xlim since abs value taken on log2FC
}

#' @name max_limits
#' @export
max_ylim <- function(.x, var = "PBHadj_log") {
  .x %>%
    purrr::map(var) %>%
    purrr::map_dbl(max) %>%
    max() %>%
    magrittr::multiply_by(1.1) %>%
    c(0, .)
}

#' Make a volcano plot specifically for data structures in TCGA studies
#'
#' @param data list of data objects for plotting
#' @param root.dir root directory to save plots
#' @param plot.suffix suffix at end of plot file name
#' @param title plot title
#' @param nki if `TRUE`, use ER binding gene lists for NKI Tier 1/2 instead of
#'   ChIP-Seq
#' @param ... named arguments to pass to `volcano_plot`
#' @export
tcga_volcano <- function(data, root.dir, plot.suffix, title, nki = FALSE, ...) {
  # Store file names
  filenames <- data %>%
    purrr::imap(~ {
      fn_split <- stringr::str_split_fixed(.y, "_", 3)
      tcga_path(root.dir = root.dir, file.name = "AgeRelated_Volcano",
                fdr = fn_split[, 2], fc = fn_split[, 3], sub.dir = "figures/",
                suffix = plot.suffix, extension = "pdf")
    })
  # Choose which ER binding gene lists to highlight
  if (nki) {
    group <- "binding_group_NKI"
    legend <- c("Tier 1", "Tier 2 Only", "Non binding")
    color <- c("red", "orange", "blue")
  } else {
    group <- "binding_group"
    legend <- c("ER binding", "Non binding")
    color <- c("red", "blue")
  }
  # Common arguments to invoke plotting function over
  args <- c(dplyr::lst(x = "Best_log2FC", y = "PBHadj_log", group,
                       label = "Hugo_Symbol", subset = "evcp",
                       legend, color, title),
            list(...))

  # Invoke plotting function on each of the objects, saving to different files
  data %>%
    purrr::set_names(filenames) %>%
    purrr::iwalk(~ purrr::invoke(volcano_plot, args, data = .x, filename = .y))
}

#' Make a volcano plot specifically for data structures in METABRIC studies
#'
#' @inheritParams tcga_volcano
#' @param ... named arguments to pass to `volcano_plot`
#' @export
metabric_volcano <- function(data, root.dir, plot.suffix, title, ...) {
  # Store file names
  filenames <- file.path(root.dir, plot.suffix)
  # Common arguments to invoke plotting function over
  args <- c(list(x = "Best_log2FC", y = "PBHadj_log", group = "binding_group",
                 label = "Gene_symbol", subset = "evcp",
                 legend = c("Tier 1", "Tier 2 Only", "Non binding"),
                 color = c("red", "orange", "blue")),
            list(...))

  # Invoke plotting function on each of the objects, saving to different files
  list(data, filenames, title) %>%
    purrr::pwalk(~ purrr::invoke(volcano_plot, args,
                                 data = ..1, filename = ..2, title = ..3))
}
