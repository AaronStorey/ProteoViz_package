# Profiling script for run_go_enrichment()
# Uses the real LCMD example data (4021 proteins x 32 samples)
# Run interactively — timed output prints to console, profvis opens in Viewer.
#
# Install if needed:
#   install.packages("profvis")
#   BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"))

library(profvis)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# ── Load real example data ────────────────────────────────────────────────────

data_dir <- "C:/Users/aaron.storey/Desktop/ProteoViz_refactor/example_inputs/StoreyAJ_20260323_02_DIA"

quant_data <- read_tsv(file.path(data_dir, "StoreyAJ_20260323_02_DIA_protein_quantitative_data.tsv"),
                       show_col_types = FALSE)
metadata   <- read_tsv(file.path(data_dir, "StoreyAJ_20260323_02_DIA_protein_metadata.tsv"),
                       show_col_types = FALSE)

# Simulate a volcano selection: take 200 random proteins
set.seed(42)
N_SELECTED    <- 200
selected_ids  <- sample(quant_data$id, N_SELECTED)
selected_data <- filter(quant_data, id %in% selected_ids)

cat(sprintf("Proteins: %d | Samples: %d | Selected: %d\n",
            nrow(quant_data), ncol(quant_data) - 1L, N_SELECTED))
cat(sprintf("Gene name column: %s\n",
            intersect(names(metadata), c("PG.Genes", "Gene_name"))))

# ── Inline copy with per-step timing ─────────────────────────────────────────

run_go_enrichment_timed <- function(quant_data, metadata, selected_data) {

  t0 <- proc.time()

  # Step 1: locate gene column and slim metadata to selected proteins
  selected_ids <- selected_data$id
  gene_col <- intersect(c("PG.Genes", "Gene_name"), names(metadata))[1]
  if (is.na(gene_col)) stop("No gene symbol column found in metadata.")

  meta_slim <- metadata |>
    filter(id %in% selected_ids) |>
    dplyr::select(id, gene = !!rlang::sym(gene_col))

  t1 <- proc.time()
  cat(sprintf("Step 1 (metadata join):      %.2f s\n", (t1 - t0)[["elapsed"]]))

  # Step 2: cohort-median-centred row means via rowMeans (fast, no rowwise)
  quant_mat     <- dplyr::select(quant_data, -id)
  cohort_means  <- tibble(id = quant_data$id,
                          mean_intensity = rowMeans(quant_mat, na.rm = TRUE))
  cohort_median <- median(cohort_means$mean_intensity, na.rm = TRUE)

  mean_intensity <- cohort_means |>
    filter(id %in% selected_ids) |>
    mutate(mean_intensity = mean_intensity - cohort_median)

  t2 <- proc.time()
  cat(sprintf("Step 2 (row means):          %.2f s\n", (t2 - t1)[["elapsed"]]))

  # Step 3: join, expand semicolons, sort, deduplicate
  gene_list_df <- mean_intensity |>
    inner_join(meta_slim, by = "id") |>
    mutate(gene = str_trim(str_extract(gene, "^[^;]+"))) |>
    filter(!is.na(gene), gene != "") |>
    arrange(desc(mean_intensity)) |>
    distinct(gene, .keep_all = TRUE)

  t3 <- proc.time()
  cat(sprintf("Step 3 (join + sort):        %.2f s\n", (t3 - t2)[["elapsed"]]))
  cat(sprintf("Gene list length:            %d\n", nrow(gene_list_df)))

  gene_list <- setNames(gene_list_df$mean_intensity, gene_list_df$gene)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Step 4: gseGO
  ego <- gseGO(
    geneList     = gene_list,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = "BP",
    minGSSize    = 5,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    scoreType    = "std",
    eps          = 0,
    verbose      = TRUE,
    BPPARAM      = BiocParallel::SerialParam()
  )

  t4 <- proc.time()
  cat(sprintf("Step 4 (gseGO):              %.2f s\n", (t4 - t3)[["elapsed"]]))
  cat(sprintf("Enriched terms found:        %d\n", nrow(ego)))

  # Step 5: dotplot
  p <- if (nrow(ego) > 0) {
    dotplot(ego, showCategory = 20)
  } else {
    message("No significant GO terms — try relaxing pvalueCutoff or selecting more proteins.")
    ggplot2::ggplot() + ggplot2::labs(title = "No significant GO terms found")
  }

  t5 <- proc.time()
  cat(sprintf("Step 5 (dotplot):            %.2f s\n", (t5 - t4)[["elapsed"]]))
  cat(sprintf("Total:                       %.2f s\n", (t5 - t0)[["elapsed"]]))

  p
}

# ── Run timed version ─────────────────────────────────────────────────────────
cat("\n--- Timed run ---\n")
p <- run_go_enrichment_timed(quant_data, metadata, selected_data)
print(p)

# ── profvis flame graph ───────────────────────────────────────────────────────
cat("\n--- profvis run (opens flame graph in Viewer) ---\n")
profvis({

  selected_ids <- selected_data$id

  gene_col  <- intersect(c("PG.Genes", "Gene_name"), names(metadata))[1]
  meta_slim <- metadata |>
    filter(id %in% selected_ids) |>
    dplyr::select(id, gene = !!rlang::sym(gene_col))

  quant_mat     <- dplyr::select(quant_data, -id)
  cohort_means  <- tibble(id = quant_data$id,
                          mean_intensity = rowMeans(quant_mat, na.rm = TRUE))
  cohort_median <- median(cohort_means$mean_intensity, na.rm = TRUE)

  mean_intensity <- cohort_means |>
    filter(id %in% selected_ids) |>
    mutate(mean_intensity = mean_intensity - cohort_median)

  gene_list_df <- mean_intensity |>
    inner_join(meta_slim, by = "id") |>
    mutate(gene = str_trim(str_extract(gene, "^[^;]+"))) |>
    filter(!is.na(gene), gene != "") |>
    arrange(desc(mean_intensity)) |>
    distinct(gene, .keep_all = TRUE)

  gene_list <- setNames(gene_list_df$mean_intensity, gene_list_df$gene)
  gene_list <- sort(gene_list, decreasing = TRUE)

  gseGO(
    geneList     = gene_list,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = "BP",
    minGSSize    = 5,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    scoreType    = "std",
    eps          = 0,
    verbose      = FALSE,
    BPPARAM      = BiocParallel::SerialParam()
  )
})
