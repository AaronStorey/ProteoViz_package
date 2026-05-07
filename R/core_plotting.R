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
#'   \eqn{-\log_{10}(\text{p-value})} on the y-axis and \eqn{\log_2} fold
#'   change on the x-axis. All points are drawn in black. Dashed threshold
#'   lines are drawn as reference guides. Rich hover labels (black background,
#'   white text) display protein description, gene name, UniProt ID, logFC,
#'   p.value, and adj.P.Val when those columns are present in \code{results}.
#'
#'   Setting \code{plotly_source} enables \code{plotly::event_data()} so that
#'   downstream Shiny outputs can react to clicks and drag-selections on the
#'   plot. The \code{key_col} value is returned in \code{event_data()$key}.
#'
#' @param results A tibble (one row per protein) with at minimum the columns
#'   \code{logFC}, \code{P.Value}, and \code{adj.P.Val}. Optional metadata
#'   columns \code{PG.ProteinDescriptions}, \code{PG.Genes}, and
#'   \code{PG.UniProtIds} are included in the hover label when present.
#'   Typically the output of \code{combine_limma_results()}, optionally
#'   left-joined with a protein metadata table before calling this function.
#' @param fc_threshold Numeric. Absolute log2 fold-change threshold for the
#'   vertical dashed guide lines. Defaults to \code{1}.
#' @param p_threshold Numeric. P-value threshold for the horizontal dashed
#'   guide line. Defaults to \code{0.05}.
#' @param use_adj_p Logical. If \code{TRUE} (default), the y-axis shows
#'   \eqn{-\log_{10}(\text{adj.P.Val})}. If \code{FALSE}, raw \code{P.Value}
#'   is used instead.
#' @param plotly_source Character. Passed as the \code{source} argument to
#'   \code{plotly::plot_ly()}, enabling \code{plotly::event_data(source = ...)}
#'   in Shiny. Defaults to \code{"proteinVolcano"}.
#' @param key_col Character. Column in \code{results} whose values are used as
#'   the plotly \code{key} — the identifier returned by
#'   \code{event_data()$key} on click or selection. Defaults to
#'   \code{"protein"}.
#' @param color_var Character or \code{NULL}. Name of a numeric column in
#'   \code{results} to use for continuous point colouring (Viridis scale with
#'   a colour bar). When \code{NULL} (default) points are drawn in solid black.
#' @param color_label Character. Label shown on the colour bar. Defaults to
#'   \code{color_var}.
#'
#' @return A \code{plotly} object.
#'
#' @examples
#' \dontrun{
#' pv <- plot_volcano(
#'   results      = limma_results[["A - B"]],
#'   fc_threshold = 1,
#'   p_threshold  = 0.05,
#'   use_adj_p    = TRUE
#' )
#' pv
#' }
#'
#' @export
plot_volcano <- function(results,
                         fc_threshold  = 1,
                         p_threshold   = 0.05,
                         use_adj_p     = TRUE,
                         plotly_source = "proteinVolcano",
                         key_col       = "protein",
                         color_var     = NULL,
                         color_label   = color_var) {

  p_col   <- if (use_adj_p) "adj.P.Val" else "P.Value"
  y_label <- if (use_adj_p) "-log10 adj.p.value" else "-log10 p.value"

  plot_data <- results |>
    dplyr::mutate(neg_log10_p = -log10(.data[[p_col]]))

  # Build per-row hover text from whatever metadata columns are present
  hover_lines <- list(
    paste0("Log2 fold change: ", round(plot_data$logFC, 3)),
    paste0("p.value: ",          signif(plot_data$P.Value,   3)),
    paste0("adj.p.value: ",      signif(plot_data$adj.P.Val, 3))
  )
  if (!is.null(color_var) && color_var %in% names(plot_data))
    hover_lines <- c(hover_lines,
                     list(paste0(color_label, ": ", plot_data[[color_var]])))
  if ("PG.UniProtIds" %in% names(plot_data))
    hover_lines <- c(list(paste0("Uniprot ID: ", plot_data$PG.UniProtIds)), hover_lines)
  if ("PG.Genes" %in% names(plot_data))
    hover_lines <- c(list(paste0("Gene name: ", plot_data$PG.Genes)), hover_lines)
  if ("PG.ProteinDescriptions" %in% names(plot_data))
    hover_lines <- c(list(paste0("Protein: ", plot_data$PG.ProteinDescriptions)), hover_lines)

  plot_data$hover_text <- do.call(paste, c(hover_lines, sep = "<br>"))

  use_color_var <- !is.null(color_var) && color_var %in% names(plot_data)

  marker_spec <- if (use_color_var) {
    color_vals <- pmin(plot_data[[color_var]], 5L)
    list(
      color        = color_vals,
      colorscale   = "Viridis",
      showscale    = TRUE,
      reversescale = TRUE,
      cmin         = 1,
      cmax         = 5,
      colorbar     = list(title = color_label),
      size         = 5,
      opacity      = 0.7
    )
  } else {
    list(color = "black", size = 5, opacity = 0.7)
  }

  plotly::plot_ly(
    data         = plot_data,
    x            = ~logFC,
    y            = ~neg_log10_p,
    key          = plot_data[[key_col]],
    text         = ~hover_text,
    type         = "scatter",
    mode         = "markers",
    marker       = marker_spec,
    source       = plotly_source,
    hovertemplate = "%{text}<extra></extra>"
  ) |>
    plotly::layout(
      dragmode   = "select",
      xaxis      = list(title = "log2 FC"),
      yaxis      = list(title = y_label),
      hoverlabel = list(
        bgcolor  = "black",
        font     = list(color = "white", size = 12)
      ),
      shapes = list(
        list(type = "line",
             x0 = fc_threshold,  x1 = fc_threshold,
             y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "grey50", width = 1, dash = "dash")),
        list(type = "line",
             x0 = -fc_threshold, x1 = -fc_threshold,
             y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "grey50", width = 1, dash = "dash")),
        list(type = "line",
             x0 = 0, x1 = 1, xref = "paper",
             y0 = -log10(p_threshold), y1 = -log10(p_threshold),
             line = list(color = "grey50", width = 1, dash = "dash"))
      )
    )
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
#' @param protein_anno Optional tibble with columns \code{id}, \code{Description},
#'   and \code{Gene_name} used to label heatmap rows. Up to 200 characters of
#'   \code{Description} are used, followed by the \code{Gene_name}.
#' @param scale_rows Logical. If \code{TRUE} (default), rows are z-score scaled
#'   and an RdBu diverging palette is used; if \code{FALSE} a viridis palette
#'   is applied to raw intensities.
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
#'   data         = normalised_long,
#'   metadata     = sample_table,
#'   proteins     = sig_proteins,
#'   protein_anno = protein_metadata
#' )
#' ph
#' }
#'
#' @export
plot_protein_heatmap <- function(data,
                                 metadata,
                                 proteins,
                                 protein_anno   = NULL,
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

  if (!is.null(protein_anno) && nrow(mat) > 0) {
    id_vec   <- rownames(mat)
    anno_sub <- protein_anno |>
      dplyr::mutate(id = as.character(id)) |>
      dplyr::filter(id %in% id_vec) |>
      dplyr::select(id,
                    dplyr::any_of(c("PG.ProteinDescriptions", "Description")),
                    dplyr::any_of(c("Gene_name", "PG.Genes"))) |>
      dplyr::distinct(id, .keep_all = TRUE)

    desc_col <- intersect(c("PG.ProteinDescriptions", "Description"), names(anno_sub))[1]
    gene_col <- intersect(c("Gene_name", "PG.Genes"),                 names(anno_sub))[1]

    label_map <- if (!is.na(desc_col) && !is.na(gene_col)) {
      stats::setNames(
        paste0(substr(anno_sub[[desc_col]], 1, 75), " | ",
               substr(anno_sub[[gene_col]], 1, 15)),
        anno_sub$id
      )
    } else if (!is.na(gene_col)) {
      stats::setNames(substr(anno_sub[[gene_col]], 1, 15), anno_sub$id)
    } else {
      stats::setNames(anno_sub$id, anno_sub$id)
    }

    rownames(mat) <- ifelse(
      id_vec %in% names(label_map), label_map[id_vec], id_vec
    )
  }

  # Order columns by design table row order, then build annotation in that order
  design_col_order <- metadata$Sample_name[metadata$Sample_name %in% colnames(mat)]
  mat <- mat[, design_col_order, drop = FALSE]

  col_annotation <- metadata |>
    dplyr::filter(Sample_name %in% colnames(mat)) |>
    dplyr::select(Sample_name, Group) |>
    as.data.frame() |>
    tibble::column_to_rownames("Sample_name")
  col_annotation <- col_annotation[colnames(mat), , drop = FALSE]

  make_heatmap <- function(Rowv = TRUE, Colv = FALSE, main = "") {
    if (scale_rows) {
      heatmaply::heatmaply(
        mat,
        scale                   = "row",
        col_side_colors         = col_annotation,
        showticklabels          = c(TRUE, TRUE),
        margins                 = c(60, 300, 40, 20),
        main                    = main,
        plot_method             = "ggplot",
        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
          low  = "blue",
          mid  = "white",
          high = "red"
        ),
        Rowv                    = Rowv,
        Colv                    = Colv
      )
    } else {
      heatmaply::heatmaply(
        mat,
        scale           = "none",
        col_side_colors = col_annotation,
        showticklabels  = c(TRUE, TRUE),
        margins         = c(60, 300, 40, 20),
        main            = main,
        plot_method     = "plotly",
        colors          = viridisLite::viridis(256),
        Rowv            = Rowv,
        Colv            = Colv
      )
    }
  }

  tryCatch(
    make_heatmap(),
    error = function(e) {
      if (!grepl("NA/NaN/Inf", conditionMessage(e), fixed = TRUE)) stop(e)
      # Drop all-NA rows and columns then retry with clustering
      mat            <<- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
      mat            <<- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
      col_annotation <<- col_annotation[colnames(mat), , drop = FALSE]
      tryCatch(
        make_heatmap(Rowv = TRUE),
        error = function(e2) {
          if (!grepl("NA/NaN/Inf", conditionMessage(e2), fixed = TRUE)) stop(e2)
          make_heatmap(Rowv = FALSE, Colv = FALSE,
                       main = "(clustering disabled — too many missing values)")
        }
      )
    }
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
#' @param scale_rows Logical. If \code{TRUE} rows are z-score scaled and the
#'   blue-white-red gradient is used; if \code{FALSE} raw intensities are shown
#'   with the viridis palette. Default \code{TRUE}.
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
                                 peptide_metadata,
                                 scale_rows = TRUE) {

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

  make_heatmap <- function(Rowv = TRUE, Colv = FALSE, main = "") {
    if (scale_rows) {
      heatmaply::heatmaply(
        pep_subset,
        scale                   = "row",
        col_side_colors         = col_annotation,
        showticklabels          = c(TRUE, TRUE),
        xlab                    = "Sample",
        ylab                    = "Peptide",
        main                    = main,
        plot_method             = "ggplot",
        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
          low  = "blue",
          mid  = "white",
          high = "red"
        ),
        Rowv                    = Rowv,
        Colv                    = Colv
      )
    } else {
      heatmaply::heatmaply(
        pep_subset,
        scale           = "none",
        col_side_colors = col_annotation,
        showticklabels  = c(TRUE, TRUE),
        xlab            = "Sample",
        ylab            = "Peptide",
        main            = main,
        plot_method     = "plotly",
        colors          = viridisLite::viridis(256),
        Rowv            = Rowv,
        Colv            = Colv
      )
    }
  }

  tryCatch(
    make_heatmap(),
    error = function(e) {
      if (!grepl("NA/NaN/Inf", conditionMessage(e), fixed = TRUE)) stop(e)
      # Drop all-NA rows and columns then retry with clustering
      pep_subset     <<- pep_subset[rowSums(!is.na(pep_subset)) > 0, , drop = FALSE]
      pep_subset     <<- pep_subset[, colSums(!is.na(pep_subset)) > 0, drop = FALSE]
      col_annotation <<- col_annotation[colnames(pep_subset), , drop = FALSE]
      tryCatch(
        make_heatmap(Rowv = TRUE),
        error = function(e2) {
          if (!grepl("NA/NaN/Inf", conditionMessage(e2), fixed = TRUE)) stop(e2)
          make_heatmap(Rowv = FALSE, Colv = FALSE,
                       main = "(clustering disabled — too many missing values)")
        }
      )
    }
  )
}
