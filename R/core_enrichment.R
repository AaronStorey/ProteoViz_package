# ── pick_populated_column ────────────────────────────────────────────────────

#' Pick the first candidate column that actually has non-NA data
#'
#' \code{intersect(candidates, names(df))[1]} only tests whether a column is
#' present, which silently prefers an all-NA column (e.g. an export where
#' \code{PG.Genes} exists but was never filled in) over one that is actually
#' populated, such as \code{Gene_name}. This checks for data, not just
#' presence.
#'
#' @param df A data frame.
#' @param candidates Character vector of column names, in preference order.
#'
#' @return The first candidate in \code{df} with at least one non-NA value,
#'   or \code{NA_character_} if none qualify.
#'
#' @keywords internal
pick_populated_column <- function(df, candidates) {
  for (col in intersect(candidates, names(df))) {
    if (any(!is.na(df[[col]]))) return(col)
  }
  NA_character_
}


# ── run_go_enrichment ────────────────────────────────────────────────────────

#' Run GO over-representation analysis on proteins selected from the volcano plot
#'
#' Tests the box-selected proteins (the foreground, "group A") for GO term
#' enrichment against a background universe ("group B") of every protein
#' identified in the experiment — not the full genome annotation database.
#' Because the hypergeometric test draws the foreground as a sample of the
#' universe, restricting the universe to the identified proteome already
#' compares group A against "every identified protein not selected"; the
#' statistics are equivalent to a direct group-A-vs-group-B contingency table.
#'
#' @param quant_data Wide tibble: column \code{id} + one column per sample.
#'   Every \code{id} here is treated as "identified" and forms the background
#'   universe.
#' @param metadata Protein metadata with \code{id} and a gene symbol column
#'   (\code{PG.Genes} or \code{Gene_name}).
#' @param selected_data Subset of \code{quant_data} for the box-selected
#'   proteins (the foreground gene set).
#' @param ont GO ontology: \code{"BP"}, \code{"MF"}, \code{"CC"}, or
#'   \code{"ALL"}.
#' @param pvalue_cutoff Adjusted p-value cutoff.
#'
#' @return A list with elements \code{ego} (an \code{enrichResult}) and
#'   \code{gene_list} (a named numeric vector of median-centred intensities
#'   for the selected genes, used only to colour downstream plots), or
#'   \code{NULL} when no terms pass the cutoff.
#'
#' @importFrom dplyr filter select mutate inner_join distinct pull
#' @importFrom tibble tibble
#' @importFrom stringr str_extract str_trim
#' @importFrom clusterProfiler enrichGO
#' @export
run_go_enrichment <- function(quant_data, metadata, selected_data,
                              ont           = "BP",
                              pvalue_cutoff = 0.05) {

  selected_ids <- selected_data$id

  gene_col <- pick_populated_column(metadata, c("PG.Genes", "Gene_name"))
  if (is.na(gene_col)) {
    stop("No gene symbol column found in metadata (expected PG.Genes or Gene_name).")
  }

  id_to_gene <- metadata |>
    dplyr::filter(id %in% quant_data$id) |>
    dplyr::mutate(gene = stringr::str_trim(stringr::str_extract(.data[[gene_col]], "^[^;]+"))) |>
    dplyr::filter(!is.na(gene), gene != "") |>
    dplyr::distinct(id, gene)

  universe <- unique(id_to_gene$gene)

  selected_genes <- id_to_gene |>
    dplyr::filter(id %in% selected_ids) |>
    dplyr::pull(gene) |>
    unique()

  if (length(selected_genes) < 5) {
    stop("Fewer than 5 gene symbols matched in the selected proteins — select more proteins before running enrichment.")
  }

  quant_mat     <- dplyr::select(quant_data, -id)
  cohort_means  <- tibble::tibble(
    id             = quant_data$id,
    mean_intensity = rowMeans(quant_mat, na.rm = TRUE)
  )
  cohort_median <- stats::median(cohort_means$mean_intensity, na.rm = TRUE)

  gene_list_df <- cohort_means |>
    dplyr::filter(id %in% selected_ids) |>
    dplyr::mutate(mean_intensity = mean_intensity - cohort_median) |>
    dplyr::inner_join(id_to_gene, by = "id") |>
    dplyr::distinct(gene, .keep_all = TRUE)

  gene_list        <- gene_list_df$mean_intensity
  names(gene_list) <- gene_list_df$gene

  ego <- clusterProfiler::enrichGO(
    gene          = selected_genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    minGSSize     = 5,
    maxGSSize     = 500,
    pvalueCutoff  = pvalue_cutoff,
    qvalueCutoff  = 1,
    pAdjustMethod = "BH"
  )

  if (is.null(ego) || nrow(ego) == 0) return(NULL)

  list(ego = ego, gene_list = gene_list)
}


# ── plot_go_dotplot ───────────────────────────────────────────────────────────

#' Interactive plotly dotplot of GO enrichment results
#'
#' @param ego An \code{enrichResult} object, or \code{NULL}.
#' @param show_category Maximum GO terms to display, ordered by adjusted p-value.
#' @param plotly_source Source string for \code{plotly::event_data}.
#'
#' @return A \code{plotly} figure.
#'
#' @importFrom dplyr arrange slice_head mutate
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient
#'   scale_size_continuous labs theme element_text theme_void
#' @importFrom plotly ggplotly
#' @export
plot_go_dotplot <- function(ego, show_category = 20, plotly_source = "goDotplot") {

  # ggplotly() measures text via the current graphics device; inside Shiny,
  # a prior renderPlot() call elsewhere in the session can leave a broken
  # device current, which makes that measurement fail with "invalid 'width'
  # or 'height'". Opening (and auto-closing) a throwaway device here keeps
  # the conversion independent of whatever device state preceded it.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (is.null(ego) || nrow(ego) == 0) {
    return(plotly::ggplotly(
      ggplot2::ggplot() +
        ggplot2::labs(title    = "No significant GO terms found",
                      subtitle = "Try selecting more proteins or relaxing the p-value cutoff.") +
        ggplot2::theme_void()
    ))
  }

  plot_df <- ego@result |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = show_category) |>
    dplyr::mutate(
      GeneRatio   = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      Description = factor(substr(Description, 1, 100),
                           levels = rev(substr(Description, 1, 100))),
      label = paste0(
        "<b>", Description, "</b>",
        "<br>GO ID: ", ID,
        "<br>Gene ratio: ", round(GeneRatio, 3),
        "<br>p.adjust: ", signif(p.adjust, 3),
        "<br>Count: ", Count
      )
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = GeneRatio, y = Description, size = Count,
                 color = p.adjust, key = ID, text = label)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "red", high = "blue") +
    ggplot2::scale_size_continuous(range = c(3, 10)) +
    ggplot2::labs(x = "Gene ratio", y = NULL, color = "p.adjust", size = "Count") +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

  plotly::ggplotly(p, tooltip = "text", source = plotly_source)
}


# ── plot_go_cnetplot ──────────────────────────────────────────────────────────

#' Category-gene network plot of GO enrichment results
#'
#' @param ego An \code{enrichResult} object, or \code{NULL}.
#' @param gene_list Named numeric vector of median-centred intensities used as
#'   fold-change colours on the gene nodes.
#' @param show_category Number of GO terms to display.
#'
#' @return A \code{ggplot2} figure.
#'
#' @importFrom enrichplot cnetplot
#' @importFrom ggplot2 ggplot labs theme_void
#' @export
plot_go_cnetplot <- function(ego, gene_list = NULL, show_category = 20) {

  if (is.null(ego) || nrow(ego) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "No significant GO terms found") +
        ggplot2::theme_void()
    )
  }

  enrichplot::cnetplot(
    ego,
    foldChange    = gene_list,
    showCategory  = show_category,
    node_label    = "all"
  )
}


# ── plot_go_heatplot ──────────────────────────────────────────────────────────

#' Heatmap of genes across enriched GO terms
#'
#' @param ego An \code{enrichResult} object, or \code{NULL}.
#' @param gene_list Named numeric vector used to colour gene tiles.
#' @param show_category Number of GO terms to display.
#'
#' @return A \code{ggplot2} figure.
#'
#' @importFrom enrichplot heatplot
#' @importFrom ggplot2 ggplot labs theme_void theme element_text
#' @export
plot_go_heatplot <- function(ego, gene_list = NULL, show_category = 20) {

  if (is.null(ego) || nrow(ego) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "No significant GO terms found") +
        ggplot2::theme_void()
    )
  }

  enrichplot::heatplot(
    ego,
    foldChange   = gene_list,
    showCategory = show_category
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7, angle = 90, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8))
}
