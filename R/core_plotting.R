# core_plotting.R
# Pure R functions for visualising protein and peptide quantitative data.
# No Shiny dependencies.
#
# Functions in this file:
#   plot_intensity_distributions()  -- violin/box plots of sample intensities
#   plot_pca()                      -- interactive PCA plot
#   plot_correlation_heatmap()      -- sample correlation heatmap
#   plot_missing_values()           -- missing value pattern heatmap
#   plot_volcano()                  -- volcano plot for a single Limma result
#   plot_protein_heatmap()          -- protein-level interactive heatmap
#   plot_peptide_heatmap()          -- peptide-level interactive heatmap


# Exported functions -----------------------------------------------------------

# plot_intensity_distributions -------------------------------------------------

#' Plot sample intensity distributions
#'
#' @title Plot sample intensity distributions as violin plots
#'
#' @description Produces a ggplot2 violin-plus-boxplot figure showing the
#'   distribution of \code{Intensity} values for each sample. Samples are
#'   arranged on the x-axis in the order they appear in \code{metadata} and
#'   are optionally coloured by a grouping variable. A horizontal dashed
#'   reference line is drawn at the median of all intensities to aid visual
#'   comparison between samples.
#'
#'   The input \code{data} is the long-form tibble produced by
#'   \code{clean_protein_data()} — one row per protein-sample observation with
#'   columns \code{id}, \code{Sample_name}, and \code{Intensity}.
#'
#' @param data A long-form tibble with at minimum the columns \code{Sample_name}
#'   and \code{Intensity}. Typically the output of \code{clean_protein_data()}.
#' @param metadata A tibble describing sample metadata. Must contain a
#'   \code{Sample_name} column whose values match those in \code{data}. The
#'   row order of \code{metadata} determines the left-to-right order of samples
#'   on the plot x-axis.
#' @param color_by A single character string giving the name of a column in
#'   \code{metadata} to use for point/violin colour. Defaults to
#'   \code{"Group"}. Pass \code{NULL} to suppress colour mapping.
#'
#' @return A \code{ggplot} object. Print it directly or save with
#'   \code{ggplot2::ggsave()}.
#'
#' @examples
#' \dontrun{
#' p <- plot_intensity_distributions(
#'   data     = protein_long,
#'   metadata = sample_table,
#'   color_by = "Group"
#' )
#' print(p)
#' ggplot2::ggsave("intensity_distributions.tiff", p,
#'                 width = 14, height = 6, dpi = 300)
#' }
#'
#' @export
plot_intensity_distributions <- function(data,
                                         metadata,
                                         color_by = "Group") {

  sample_levels <- metadata[["Sample_name"]]

  plot_data <- data |>
    dplyr::filter(Sample_name %in% sample_levels) |>
    dplyr::mutate(Sample_name = factor(Sample_name, levels = sample_levels))

  if (!is.null(color_by) &&
      color_by %in% colnames(metadata) &&
      !color_by %in% colnames(plot_data)) {
    plot_data <- plot_data |>
      dplyr::left_join(
        dplyr::select(metadata, Sample_name, dplyr::all_of(color_by)),
        by = "Sample_name"
      )
  }

  global_median <- stats::median(plot_data[["Intensity"]], na.rm = TRUE)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Sample_name, y = Intensity)
  )

  if (!is.null(color_by) && color_by %in% colnames(plot_data)) {
    p <- p +
      ggplot2::geom_violin(
        ggplot2::aes(fill = .data[[color_by]]),
        alpha = 0.6, trim = TRUE
      ) +
      ggplot2::geom_boxplot(
        ggplot2::aes(colour = .data[[color_by]]),
        width = 0.15, outlier.size = 0.5, fill = "white", alpha = 0.8
      )
  } else {
    p <- p +
      ggplot2::geom_violin(fill = "steelblue", alpha = 0.6, trim = TRUE) +
      ggplot2::geom_boxplot(
        width = 0.15, outlier.size = 0.5, fill = "white", alpha = 0.8
      )
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = global_median,
      linetype   = "dashed",
      colour     = "grey40",
      linewidth  = 0.4
    ) +
    ggplot2::labs(
      x      = NULL,
      y      = "log\u2082 Intensity",
      fill   = color_by,
      colour = color_by
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "right"
    )

  p
}


# plot_pca ---------------------------------------------------------------------

#' Plot an interactive PCA of sample intensities
#'
#' @title Interactive PCA plot of samples
#'
#' @description Runs principal component analysis on the protein intensity
#'   matrix (samples as rows, proteins as columns) and produces an interactive
#'   \code{plotly} scatter plot of PC1 vs PC2. Each point represents a sample,
#'   labelled with \code{Sample_name} and coloured by \code{Group}. The
#'   percentage of variance explained by each PC is shown on the axis labels.
#'   Proteins with any missing value across samples are excluded before PCA.
#'
#' @param data A long-form tibble with columns \code{id}, \code{Sample_name},
#'   and \code{Intensity}. Typically the output of \code{normalize_intensities()}.
#' @param metadata A tibble with at minimum the columns \code{Sample_name},
#'   \code{Group}, and \code{Batch} used for hover text and colour.
#'
#' @return A \code{plotly} object.
#'
#' @examples
#' \dontrun{
#' pca_plot <- plot_pca(protein_long, sample_table)
#' pca_plot
#' }
#'
#' @export
plot_pca <- function(data, metadata) {

  # Pivot to wide, drop proteins with any NA, transpose so samples are rows
  mat <- data |>
    dplyr::select(id, Sample_name, Intensity) |>
    tidyr::pivot_wider(names_from = id, values_from = Intensity) |>
    tibble::column_to_rownames("Sample_name") |>
    as.matrix()

  mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]

  pca_res    <- stats::prcomp(mat, scale. = TRUE, center = TRUE)
  pct_var    <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  scores     <- as.data.frame(pca_res$x[, 1:2])
  scores[["Sample_name"]] <- rownames(scores)

  plot_data <- scores |>
    dplyr::left_join(metadata, by = "Sample_name")

  plotly::plot_ly(
    data   = plot_data,
    x      = ~PC1,
    y      = ~PC2,
    color  = ~Group,
    text   = ~paste0(Sample_name, "<br>Batch: ", Batch),
    type   = "scatter",
    mode   = "markers",
    marker = list(size = 10)
  ) |>
    plotly::layout(
      xaxis = list(title = paste0("PC1 (", pct_var[1], "% variance)")),
      yaxis = list(title = paste0("PC2 (", pct_var[2], "% variance)"))
    )
}


# plot_correlation_heatmap -----------------------------------------------------

#' Plot a sample-level Pearson correlation heatmap
#'
#' @title Sample correlation heatmap
#'
#' @description Computes pairwise Pearson correlations between samples and
#'   displays the result as a \code{ComplexHeatmap} heatmap with a red colour
#'   scale. Proteins with any missing value across samples are excluded before
#'   computing correlations.
#'
#' @param data A long-form tibble with columns \code{id}, \code{Sample_name},
#'   and \code{Intensity}.
#' @param metadata A tibble with a \code{Sample_name} column used to order
#'   samples on both axes of the heatmap.
#' @param scale_rows Logical. Currently unused (retained for API consistency
#'   with the original \code{make_cor_plot}). Defaults to \code{FALSE}.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object. Print it to display.
#'
#' @examples
#' \dontrun{
#' cor_plot <- plot_correlation_heatmap(protein_long, sample_table)
#' print(cor_plot)
#' }
#'
#' @export
plot_correlation_heatmap <- function(data, metadata, scale_rows = FALSE) {

  sample_levels <- metadata[["Sample_name"]]

  mat <- data |>
    dplyr::filter(Sample_name %in% sample_levels) |>
    dplyr::select(id, Sample_name, Intensity) |>
    tidyr::pivot_wider(names_from = Sample_name, values_from = Intensity) |>
    tibble::column_to_rownames("id") |>
    as.matrix()

  # Drop proteins with any NA so cor() works without use = "pairwise"
  mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
  mat <- mat[, sample_levels[sample_levels %in% colnames(mat)], drop = FALSE]

  cor_mat <- stats::cor(mat, method = "pearson")

  ComplexHeatmap::Heatmap(
    cor_mat,
    name              = "Pearson r",
    col               = circlize::colorRamp2(c(0.8, 1), c("white", "#E41A1C")),
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = grid::gpar(fontsize = 8),
    column_names_gp   = grid::gpar(fontsize = 8),
    column_names_rot  = 45
  )
}


# plot_missing_values ----------------------------------------------------------

#' Plot missing value patterns across samples
#'
#' @title Missing value pattern heatmap
#'
#' @description Produces a \code{ComplexHeatmap} binary heatmap where each cell
#'   indicates whether a protein was observed (black) or missing (white) in a
#'   given sample. Proteins are sorted by the number of missing values
#'   (descending) so that the most-incomplete proteins appear at the top.
#'
#' @param data A long-form tibble with columns \code{id}, \code{Sample_name},
#'   and \code{Intensity}.
#' @param metadata A tibble with a \code{Sample_name} column used to order
#'   samples on the x-axis.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object.
#'
#' @examples
#' \dontrun{
#' na_plot <- plot_missing_values(protein_long, sample_table)
#' print(na_plot)
#' }
#'
#' @export
plot_missing_values <- function(data, metadata) {

  sample_levels <- metadata[["Sample_name"]]

  mat <- data |>
    dplyr::filter(Sample_name %in% sample_levels) |>
    dplyr::select(id, Sample_name, Intensity) |>
    tidyr::pivot_wider(names_from = Sample_name, values_from = Intensity) |>
    tibble::column_to_rownames("id") |>
    as.matrix()

  mat <- mat[, sample_levels[sample_levels %in% colnames(mat)], drop = FALSE]

  # Binary presence/absence matrix; sort rows by missingness
  binary_mat  <- (!is.na(mat)) * 1L
  row_order   <- order(rowSums(binary_mat))
  binary_mat  <- binary_mat[row_order, , drop = FALSE]

  ComplexHeatmap::Heatmap(
    binary_mat,
    name              = "Observed",
    col               = c("0" = "white", "1" = "black"),
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = FALSE,
    show_column_names = TRUE,
    column_names_gp   = grid::gpar(fontsize = 8),
    column_names_rot  = 45,
    heatmap_legend_param = list(
      at     = c(0, 1),
      labels = c("Missing", "Observed")
    )
  )
}


# plot_volcano -----------------------------------------------------------------

#' Plot a volcano plot for a single Limma result table
#'
#' @title Volcano plot of differential expression results
#'
#' @description Produces an interactive \code{plotly} scatter plot with
#'   \eqn{-\log_{10}(\text{adjusted p-value})} on the y-axis and
#'   \eqn{\log_2} fold change on the x-axis. Points are coloured by
#'   significance category: \emph{Up}, \emph{Down}, or \emph{NS}. Dashed
#'   threshold lines are drawn as reference guides.
#'
#' @param results A tibble (one row per protein) with at minimum the columns
#'   \code{logFC}, \code{adj.P.Val}, and the column named by \code{label_col}.
#'   Typically one element of the list returned by \code{run_limma_protein()}.
#' @param fc_threshold Numeric. Absolute log2 fold-change threshold. Defaults
#'   to \code{1}.
#' @param p_threshold Numeric. Adjusted p-value threshold. Defaults to
#'   \code{0.05}.
#' @param label_col Character. Column in \code{results} to use as hover-text
#'   label. Defaults to \code{"protein"}.
#'
#' @return A \code{plotly} object.
#'
#' @examples
#' \dontrun{
#' pv <- plot_volcano(
#'   results      = limma_results[["A - B"]],
#'   fc_threshold = 1,
#'   p_threshold  = 0.05,
#'   label_col    = "Gene_name"
#' )
#' pv
#' }
#'
#' @export
plot_volcano <- function(results,
                         fc_threshold = 1,
                         p_threshold  = 0.05,
                         label_col    = "protein") {

  plot_data <- results |>
    dplyr::mutate(
      neg_log10_p  = -log10(adj.P.Val),
      significance = dplyr::case_when(
        logFC >  fc_threshold & adj.P.Val < p_threshold ~ "Up",
        logFC < -fc_threshold & adj.P.Val < p_threshold ~ "Down",
        TRUE                                             ~ "NS"
      ),
      significance = factor(significance, levels = c("Up", "Down", "NS"))
    )

  sig_colours <- c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey70")

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x      = logFC,
      y      = neg_log10_p,
      colour = significance,
      text   = .data[[label_col]]
    )
  ) +
    ggplot2::geom_point(size = 1.5, alpha = 0.7) +
    ggplot2::geom_vline(
      xintercept = c(-fc_threshold, fc_threshold),
      linetype   = "dashed", colour = "grey50", linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_threshold),
      linetype   = "dashed", colour = "grey50", linewidth = 0.4
    ) +
    ggplot2::scale_colour_manual(values = sig_colours, name = NULL) +
    ggplot2::labs(
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "Adjusted p-value")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "right")

  plotly::ggplotly(p, tooltip = c("text", "x", "y"))
}


# plot_protein_heatmap ---------------------------------------------------------

#' Plot an interactive protein-level heatmap
#'
#' @title Interactive heatmap of protein intensities across samples
#'
#' @description Produces an interactive \code{heatmaply} heatmap showing
#'   normalised intensities for a selected set of proteins. Rows are optionally
#'   z-score scaled. Groups listed in \code{exclude_groups} are removed before
#'   plotting.
#'
#' @param data A long-form tibble with columns \code{id}, \code{Sample_name},
#'   and \code{Intensity}.
#' @param metadata A tibble with columns \code{Sample_name} and \code{Group}.
#' @param proteins A character vector of \code{id} values to include.
#' @param scale_rows Logical. If \code{TRUE} (default), rows are z-score scaled.
#' @param exclude_groups A character vector of group labels to exclude from the
#'   heatmap. Defaults to \code{NULL}.
#'
#' @return A \code{plotly}/\code{heatmaply} object.
#'
#' @examples
#' \dontrun{
#' sig_proteins <- limma_results[["A - B"]] |>
#'   dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
#'   dplyr::pull(protein)
#'
#' ph <- plot_protein_heatmap(
#'   data       = normalised_long,
#'   metadata   = sample_table,
#'   proteins   = sig_proteins
#' )
#' ph
#' }
#'
#' @export
plot_protein_heatmap <- function(data,
                                 metadata,
                                 proteins,
                                 scale_rows     = TRUE,
                                 exclude_groups = NULL) {

  if (!is.null(exclude_groups) && length(exclude_groups) > 0) {
    keep_samples <- metadata |>
      dplyr::filter(!Group %in% exclude_groups) |>
      dplyr::pull(Sample_name)
    data     <- dplyr::filter(data, Sample_name %in% keep_samples)
    metadata <- dplyr::filter(metadata, Sample_name %in% keep_samples)
  }

  mat <- data |>
    dplyr::filter(id %in% proteins) |>
    dplyr::select(id, Sample_name, Intensity) |>
    tidyr::pivot_wider(names_from = Sample_name, values_from = Intensity) |>
    as.data.frame() |>
    tibble::column_to_rownames("id") |>
    as.matrix()

  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

  col_annotation <- metadata |>
    dplyr::filter(Sample_name %in% colnames(mat)) |>
    dplyr::select(Sample_name, Group) |>
    as.data.frame() |>
    tibble::column_to_rownames("Sample_name")
  col_annotation <- col_annotation[colnames(mat), , drop = FALSE]

  heatmaply::heatmaply(
    mat,
    scale           = if (scale_rows) "row" else "none",
    col_side_colors = col_annotation,
    showticklabels  = c(TRUE, FALSE),
    xlab            = "Sample",
    ylab            = "Protein",
    main            = "",
    plot_method     = "plotly",
    colors          = grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    )(256)
  )
}


# plot_peptide_heatmap ---------------------------------------------------------

#' Plot an interactive peptide-level heatmap
#'
#' @title Interactive heatmap of peptide intensities for selected proteins
#'
#' @description Produces an interactive \code{heatmaply} heatmap showing
#'   individual peptide intensities for one or more selected proteins. Peptide
#'   rows are labelled with the precursor sequence when available.
#'
#' @param peptide_data A long-form tibble with columns \code{id},
#'   \code{Sample_name}, and \code{Intensity}.
#' @param metadata A tibble with columns \code{Sample_name} and \code{Group}.
#' @param proteins A character vector of protein accession identifiers.
#' @param peptide_metadata A tibble of peptide annotations with columns
#'   \code{id} and \code{PG.ProteinAccessions}. Optionally a \code{Sequence}
#'   column provides row labels.
#'
#' @return A \code{plotly}/\code{heatmaply} object.
#'
#' @examples
#' \dontrun{
#' pph <- plot_peptide_heatmap(
#'   peptide_data     = peptide_long,
#'   metadata         = sample_table,
#'   proteins         = c("P04406"),
#'   peptide_metadata = peptide_meta
#' )
#' pph
#' }
#'
#' @export
plot_peptide_heatmap <- function(peptide_data,
                                 metadata,
                                 proteins,
                                 peptide_metadata) {

  protein_pattern <- paste(proteins, collapse = "|")

  matching_ids <- peptide_metadata |>
    dplyr::filter(grepl(protein_pattern, PG.ProteinAccessions, fixed = FALSE)) |>
    dplyr::pull(id) |>
    as.character()

  if (length(matching_ids) == 0) {
    stop(
      "No peptides found for the requested protein(s): ",
      paste(proteins, collapse = ", ")
    )
  }

  pep_subset <- peptide_data |>
    dplyr::filter(as.character(id) %in% matching_ids) |>
    dplyr::select(id, Sample_name, Intensity) |>
    tidyr::pivot_wider(names_from = Sample_name, values_from = Intensity) |>
    as.data.frame() |>
    tibble::column_to_rownames("id") |>
    as.matrix()

  pep_subset <- pep_subset[rowSums(!is.na(pep_subset)) > 0, , drop = FALSE]

  if ("Sequence" %in% colnames(peptide_metadata)) {
    label_map <- peptide_metadata |>
      dplyr::filter(as.character(id) %in% rownames(pep_subset)) |>
      dplyr::mutate(id = as.character(id)) |>
      dplyr::select(id, Sequence) |>
      dplyr::distinct()

    new_rownames <- label_map[["Sequence"]][
      match(rownames(pep_subset), label_map[["id"]])
    ]
    new_rownames[is.na(new_rownames)] <- rownames(pep_subset)[is.na(new_rownames)]
    rownames(pep_subset) <- new_rownames
  }

  col_annotation <- metadata |>
    dplyr::filter(Sample_name %in% colnames(pep_subset)) |>
    dplyr::select(Sample_name, Group) |>
    as.data.frame() |>
    tibble::column_to_rownames("Sample_name")
  col_annotation <- col_annotation[colnames(pep_subset), , drop = FALSE]

  heatmaply::heatmaply(
    pep_subset,
    scale           = "row",
    col_side_colors = col_annotation,
    showticklabels  = c(TRUE, TRUE),
    xlab            = "Sample",
    ylab            = "Peptide",
    main            = "",
    plot_method     = "plotly",
    colors          = grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    )(256)
  )
}
