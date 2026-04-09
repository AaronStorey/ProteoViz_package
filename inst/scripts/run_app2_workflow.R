# =============================================================================
# run_app2_workflow.R
# ProteoViz App2 — Statistical Analysis CLI Workflow (TASK-021)
#
# Performs a complete limma differential expression analysis on DIA proteomics
# data entirely from the command line, with no Shiny dependency.
#
# Workflow sections:
#   1. Setup          — libraries, file paths, output directory
#   2. Load data      — read TSV inputs
#   3. Prepare matrix — pivot protein quant to matrix and long form
#   4. Build matrices — design and contrast matrices via core_limma.R helpers
#   5. Run Limma      — run_limma_protein(); print significant hits per contrast
#   6. Save summary   — write_limma_summary() wide-format TSV
#   7. Volcano plots  — plot_volcano() per contrast, saved as HTML widgets
#   8. Protein heatmap — top 50 significant proteins across all contrasts
#   9. Peptide heatmap — peptide-level plot for the top protein by |logFC|
# =============================================================================


# =============================================================================
# 1. SETUP
# =============================================================================

library(ProteoViz)       # exposes build_design_matrix, build_contrast_matrix,
                         # run_limma_protein, write_limma_summary,
                         # plot_volcano, plot_protein_heatmap, plot_peptide_heatmap
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(purrr)
library(htmlwidgets)

# ---- Project identity -------------------------------------------------------

project_name <- "StoreyAJ_20260323_02_DIA"

# ---- Input file paths -------------------------------------------------------

input_dir <- "C:/Users/aaron.storey/Desktop/ProteoViz_refactor/example_inputs/StoreyAJ_20260323_02_DIA"

path_protein_quant  <- file.path(
  input_dir,
  "StoreyAJ_20260323_02_DIA_protein_quantitative_data.tsv"
)
path_protein_meta   <- file.path(
  input_dir,
  "StoreyAJ_20260323_02_DIA_protein_metadata.tsv"
)
path_peptide_quant  <- file.path(
  input_dir,
  "StoreyAJ_20260323_02_DIA_peptide_quantitative_data.tsv"
)
path_peptide_meta   <- file.path(
  input_dir,
  "StoreyAJ_20260323_02_DIA_peptide_metadata.tsv"
)
path_design_table   <- file.path(
  input_dir,
  "StoreyAJ_20260323_02_DIA_design_table.tsv"
)
path_contrast_list  <- file.path(
  input_dir,
  "StoreyAJ_20260323_02_DIA_Sample_contrast_list.tsv"
)

# ---- Output directory -------------------------------------------------------

output_dir <- file.path("output", project_name)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Output will be written to: ", normalizePath(output_dir, mustWork = FALSE))


# =============================================================================
# 2. LOAD DATA
# =============================================================================

message("\n--- Loading data ---")

# Protein quantitative data — wide format (id + one column per sample)
protein_quant <- read_tsv(path_protein_quant, show_col_types = FALSE)
message("protein_quant: ", nrow(protein_quant), " proteins x ", ncol(protein_quant) - 1, " samples")

# Protein metadata — one row per protein id
protein_meta  <- read_tsv(path_protein_meta, show_col_types = FALSE)
message("protein_meta:  ", nrow(protein_meta), " rows, columns: ",
        paste(colnames(protein_meta), collapse = ", "))

# Peptide quantitative data — wide format (id + one column per sample)
peptide_quant <- read_tsv(path_peptide_quant, show_col_types = FALSE)
message("peptide_quant: ", nrow(peptide_quant), " peptides x ", ncol(peptide_quant) - 1, " samples")

# Peptide metadata — one row per peptide id, includes PG.ProteinAccessions
peptide_meta  <- read_tsv(path_peptide_meta, show_col_types = FALSE)
message("peptide_meta:  ", nrow(peptide_meta), " rows, columns: ",
        paste(colnames(peptide_meta), collapse = ", "))

# Experimental design table — Sample_name, Group, Batch, Block
design_table  <- read_tsv(path_design_table, show_col_types = FALSE)
message("design_table:  ", nrow(design_table), " samples")

# Contrast list — single column Contrast_name containing contrast strings
contrast_list <- read_tsv(path_contrast_list, show_col_types = FALSE)
message("contrast_list: ", nrow(contrast_list), " contrast(s): ",
        paste(contrast_list$Contrast_name, collapse = " | "))


# =============================================================================
# 3. PREPARE PROTEIN MATRIX
# =============================================================================

message("\n--- Preparing protein matrix ---")

# --- Matrix form: proteins (rows) x samples (columns) -----------------------
# column_to_rownames requires a data.frame; as.matrix converts to numeric
protein_matrix <- protein_quant |>
  as.data.frame() |>
  tibble::column_to_rownames("id") |>
  as.matrix()

message("protein_matrix dimensions: ", nrow(protein_matrix), " x ", ncol(protein_matrix))

# --- Long form: used by plot_protein_heatmap ---------------------------------
# pivot_longer over all sample columns, retaining the id column as an
# identifier; produces columns id, Sample_name, Intensity
protein_long <- protein_quant |>
  tidyr::pivot_longer(
    cols      = -id,
    names_to  = "Sample_name",
    values_to = "Intensity"
  )

# --- Long form for peptide data: used by plot_peptide_heatmap ----------------
peptide_long <- peptide_quant |>
  tidyr::pivot_longer(
    cols      = -id,
    names_to  = "Sample_name",
    values_to = "Intensity"
  )


# =============================================================================
# 4. BUILD DESIGN AND CONTRAST MATRICES
# =============================================================================

message("\n--- Building design and contrast matrices ---")

# design_matrix: one row per sample, one column per group (~0 + Group)
design_matrix   <- build_design_matrix(design_table)

message("Design matrix groups: ", paste(colnames(design_matrix), collapse = ", "))

# contrast_matrix: encodes each comparison as a numeric vector over group means
contrast_matrix <- build_contrast_matrix(design_matrix, contrast_list)

message("Contrasts defined: ", paste(colnames(contrast_matrix), collapse = " | "))


# =============================================================================
# 5. RUN LIMMA
# =============================================================================

message("\n--- Running limma ---")

# run_limma_protein() returns a named list of tibbles, one per contrast.
# Each tibble has columns: protein, logFC, AveExpr, t, P.Value, adj.P.Val
limma_results <- run_limma_protein(
  data           = protein_matrix,
  design_matrix  = design_matrix,
  contrast_matrix = contrast_matrix,
  sample_table   = design_table
)

# --- Print number of significant hits per contrast --------------------------
# Significance criteria: adj.P.Val < 0.05 AND |logFC| > 1
message("\nSignificant hits per contrast (adj.P.Val < 0.05, |logFC| > 1):")

sig_summary <- purrr::imap_dfr(limma_results, function(tbl, contrast_name) {
  n_up   <- tbl |> dplyr::filter(adj.P.Val < 0.05,  logFC >  1) |> nrow()
  n_down <- tbl |> dplyr::filter(adj.P.Val < 0.05,  logFC < -1) |> nrow()
  tibble::tibble(
    Contrast = contrast_name,
    Up       = n_up,
    Down     = n_down,
    Total    = n_up + n_down
  )
})

print(sig_summary)


# =============================================================================
# 6. SAVE LIMMA SUMMARY
# =============================================================================

message("\n--- Saving limma summary ---")

# Combine all per-contrast tibbles into a single long tibble for write_limma_summary
combined_results <- purrr::imap_dfr(
  limma_results,
  function(tbl, contrast_name) dplyr::mutate(tbl, Comparison = contrast_name)
)

# write_limma_summary() writes <prefix>_Summarized_output.tsv
# It expects quant_data with an "id" column and metadata with an "id" column
summary_prefix <- file.path(output_dir, project_name)

write_limma_summary(
  quant_data       = protein_quant,
  combined_results = combined_results,
  metadata         = protein_meta,
  prefix           = summary_prefix
)

message("Summary TSV written to: ", summary_prefix, "_Summarized_output.tsv")


# =============================================================================
# 7. VOLCANO PLOTS
# =============================================================================

message("\n--- Generating volcano plots ---")

# One interactive plotly volcano plot per contrast, saved as a self-contained
# HTML file via htmlwidgets::saveWidget()
purrr::iwalk(limma_results, function(tbl, contrast_name) {
  # Sanitise contrast name for use in a filename (replace spaces and special
  # characters with underscores)
  safe_name  <- gsub("[^A-Za-z0-9_.-]", "_", contrast_name)
  out_file   <- file.path(output_dir, paste0(project_name, "_volcano_", safe_name, ".html"))

  volcano_plot <- plot_volcano(
    results      = tbl,
    fc_threshold = 1,
    p_threshold  = 0.05,
    label_col    = "protein"
  )

  htmlwidgets::saveWidget(
    widget   = volcano_plot,
    file     = normalizePath(out_file, mustWork = FALSE),
    selfcontained = TRUE
  )

  message("  Volcano saved: ", out_file)
})


# =============================================================================
# 8. PROTEIN HEATMAP — top 50 significant proteins
# =============================================================================

message("\n--- Generating protein heatmap ---")

# Collect significant proteins from all contrasts, ranked by |logFC|, keeping
# up to the top 50 unique identifiers
top_proteins <- combined_results |>
  dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
  dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC))) |>
  dplyr::distinct(protein) |>
  dplyr::slice_head(n = 50) |>
  dplyr::pull(protein) |>
  as.character()

message("  Number of proteins for heatmap: ", length(top_proteins))

if (length(top_proteins) > 0) {
  protein_heatmap <- plot_protein_heatmap(
    data           = protein_long,
    metadata       = design_table,
    proteins       = top_proteins,
    scale_rows     = TRUE,
    exclude_groups = NULL
  )

  heatmap_file <- file.path(output_dir, paste0(project_name, "_protein_heatmap.html"))

  htmlwidgets::saveWidget(
    widget        = protein_heatmap,
    file          = normalizePath(heatmap_file, mustWork = FALSE),
    selfcontained = TRUE
  )

  message("  Protein heatmap saved: ", heatmap_file)
} else {
  message("  No significant proteins found — skipping protein heatmap.")
}


# =============================================================================
# 9. PEPTIDE HEATMAP — top protein by |logFC|
# =============================================================================

message("\n--- Generating peptide heatmap ---")

# Identify the single protein with the largest absolute log2 fold change
# across all contrasts, among those meeting the significance thresholds
top_protein_row <- combined_results |>
  dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
  dplyr::slice_max(abs(logFC), n = 1, with_ties = FALSE)

if (nrow(top_protein_row) > 0) {
  top_protein_id <- as.character(top_protein_row$protein)

  message(
    "  Top protein by |logFC|: ", top_protein_id,
    " (logFC = ", round(top_protein_row$logFC, 3),
    ", contrast: ", top_protein_row$Comparison, ")"
  )

  # plot_peptide_heatmap() matches on PG.ProteinAccessions in peptide_meta.
  # The protein_meta Uniprot_ID or PG.ProteinGroups column maps numeric ids to
  # accession strings; resolve the accession for the top protein id here.
  top_protein_accession <- protein_meta |>
    dplyr::filter(as.character(id) == top_protein_id) |>
    dplyr::pull(PG.ProteinGroups) |>
    dplyr::first()

  if (!is.na(top_protein_accession) && length(top_protein_accession) > 0) {
    message("  Resolved accession: ", top_protein_accession)

    peptide_heatmap <- plot_peptide_heatmap(
      peptide_data     = peptide_long,
      metadata         = design_table,
      proteins         = top_protein_accession,
      peptide_metadata = peptide_meta
    )

    pep_heatmap_file <- file.path(
      output_dir,
      paste0(project_name, "_peptide_heatmap_", gsub("[^A-Za-z0-9_.-]", "_", top_protein_accession), ".html")
    )

    htmlwidgets::saveWidget(
      widget        = peptide_heatmap,
      file          = normalizePath(pep_heatmap_file, mustWork = FALSE),
      selfcontained = TRUE
    )

    message("  Peptide heatmap saved: ", pep_heatmap_file)
  } else {
    message("  Could not resolve protein accession for id '", top_protein_id,
            "' — skipping peptide heatmap.")
  }
} else {
  message("  No significant proteins found — skipping peptide heatmap.")
}

message("\n=== App2 workflow complete ===")
message("All outputs written to: ", normalizePath(output_dir, mustWork = FALSE))
