# ── run_go_enrichment ────────────────────────────────────────────────────────

#' Run GO enrichment analysis on proteins selected from the volcano plot
#'
#' Returns a named list with the raw \code{gseaResult} object and the ranked
#' gene list so that downstream plotting functions can colour nodes by
#' expression level.
#'
#' @param quant_data Wide tibble: column \code{id} + one column per sample.
#' @param metadata Protein metadata with \code{id} and a gene symbol column
#'   (\code{PG.Genes} or \code{Gene_name}).
#' @param selected_data Subset of \code{quant_data} for the selected proteins.
#' @param ont GO ontology: \code{"BP"}, \code{"MF"}, \code{"CC"}, or
#'   \code{"ALL"}.
#' @param score_type GSEA score type: \code{"std"}, \code{"pos"}, or
#'   \code{"neg"}.
#' @param pvalue_cutoff Adjusted p-value cutoff.
#'
#' @return A list with elements \code{ego} (a \code{gseaResult}) and
#'   \code{gene_list} (the named numeric ranking vector), or \code{NULL} when
#'   no terms pass the cutoff.
#'
#' @importFrom dplyr filter select mutate inner_join arrange distinct any_of desc
#' @importFrom tibble tibble
#' @importFrom stringr str_extract str_trim
#' @importFrom clusterProfiler gseGO
#' @export
run_go_enrichment <- function(quant_data, metadata, selected_data,
                              ont           = "BP",
                              score_type    = "std",
                              pvalue_cutoff = 0.05) {

  selected_ids <- selected_data$id

  gene_col <- intersect(c("PG.Genes", "Gene_name"), names(metadata))[1]
  if (is.na(gene_col)) {
    stop("No gene symbol column found in metadata (expected PG.Genes or Gene_name).")
  }

  meta_slim <- metadata |>
    dplyr::filter(id %in% selected_ids) |>
    dplyr::select(id, gene = !!rlang::sym(gene_col))

  quant_mat     <- dplyr::select(quant_data, -id)
  cohort_means  <- tibble::tibble(
    id             = quant_data$id,
    mean_intensity = rowMeans(quant_mat, na.rm = TRUE)
  )
  cohort_median <- stats::median(cohort_means$mean_intensity, na.rm = TRUE)

  mean_intensity <- cohort_means |>
    dplyr::filter(id %in% selected_ids) |>
    dplyr::mutate(mean_intensity = mean_intensity - cohort_median)

  gene_list_df <- mean_intensity |>
    dplyr::inner_join(meta_slim, by = "id") |>
    dplyr::mutate(gene = stringr::str_trim(stringr::str_extract(gene, "^[^;]+"))) |>
    dplyr::filter(!is.na(gene), gene != "") |>
    dplyr::arrange(dplyr::desc(mean_intensity)) |>
    dplyr::distinct(gene, .keep_all = TRUE)

  if (nrow(gene_list_df) < 5) {
    stop("Fewer than 5 gene symbols matched — select more proteins before running enrichment.")
  }

  gene_list        <- gene_list_df$mean_intensity
  names(gene_list) <- gene_list_df$gene
  gene_list        <- sort(gene_list, decreasing = TRUE)

  ego <- clusterProfiler::gseGO(
    geneList     = gene_list,
    OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = ont,
    minGSSize    = 5,
    maxGSSize    = 500,
    pvalueCutoff = pvalue_cutoff,
    scoreType    = score_type,
    eps          = 0,
    verbose      = FALSE,
    BPPARAM      = BiocParallel::SerialParam()
  )

  if (nrow(ego) == 0) return(NULL)

  list(ego = ego, gene_list = gene_list)
}


# ── plot_go_dotplot ───────────────────────────────────────────────────────────

#' Interactive plotly dotplot of GO enrichment results
#'
#' @param ego A \code{gseaResult} object, or \code{NULL}.
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
      Description = factor(substr(Description, 1, 100),
                           levels = rev(substr(Description, 1, 100))),
      label = paste0(
        "<b>", Description, "</b>",
        "<br>GO ID: ", ID,
        "<br>NES: ", round(NES, 3),
        "<br>p.adjust: ", signif(p.adjust, 3),
        "<br>Set size: ", setSize
      )
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = NES, y = Description, size = setSize,
                 color = p.adjust, key = ID, text = label)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "red", high = "blue") +
    ggplot2::scale_size_continuous(range = c(3, 10)) +
    ggplot2::labs(x = "NES", y = NULL, color = "p.adjust", size = "Set size") +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

  plotly::ggplotly(p, tooltip = "text", source = plotly_source)
}


# ── plot_go_cnetplot ──────────────────────────────────────────────────────────

#' Category-gene network plot of GO enrichment results
#'
#' @param ego A \code{gseaResult} object, or \code{NULL}.
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
#' @param ego A \code{gseaResult} object, or \code{NULL}.
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


# ── plot_go_gseaplot ──────────────────────────────────────────────────────────

#' GSEA running-score plot for a selected GO term
#'
#' Shows the enrichment plot for \code{go_id} if supplied and present in
#' \code{ego}, otherwise defaults to the top-ranked term.
#'
#' @param ego A \code{gseaResult} object, or \code{NULL}.
#' @param go_id GO term ID string, or \code{NULL} to use the top term.
#'
#' @return A \code{ggplot2} figure.
#'
#' @importFrom enrichplot gseaplot2
#' @importFrom ggplot2 ggplot labs theme_void
#' @export
plot_go_gseaplot <- function(ego, go_id = NULL) {

  if (is.null(ego) || nrow(ego) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "No significant GO terms found") +
        ggplot2::theme_void()
    )
  }

  gene_set_id <- if (!is.null(go_id) && go_id %in% ego@result$ID) go_id else 1

  enrichplot::gseaplot2(ego, geneSetID = gene_set_id)
}
