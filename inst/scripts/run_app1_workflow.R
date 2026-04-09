# =============================================================================
# run_app1_workflow.R
# App1 scripted workflow: DIA data cleaning and normalisation
#
# Mirrors the interactive steps performed in runApp1DIA() but runs entirely
# from the R terminal without launching Shiny.  Adjust the parameters in
# each section to match your experiment before sourcing the script.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 1: Setup
# -----------------------------------------------------------------------------

library(ProteoViz)   # loads all core_* functions via the package namespace
library(dplyr)
library(readr)
library(ggplot2)

# Project name — used as the output sub-directory and file prefix
project_name <- "StoreyAJ_20260323_02_DIA"

# Data type understood by the readers ("DIA" or "Spectronaut")
data_type <- "DIA"

# Normalisation method — one of:
#   "None", "Mean", "Median", "Quantile",
#   "Cyclic Loess", "Cyclic Loess then Median", "Rlr", "Gi"
norm_method <- "Median"

# log2-transform positive intensities before filtering / normalisation?
log2_transform <- TRUE

# Remove zero / NA intensities?
remove_zeros <- TRUE

# Create the output directory if it does not yet exist
if (!dir.exists(project_name)) dir.create(project_name, recursive = TRUE)


# -----------------------------------------------------------------------------
# Section 2: Load data
# -----------------------------------------------------------------------------

protein_path <- "C:/Users/aaron.storey/Desktop/ProteoViz_refactor/example_inputs/StoreyAJ_20260323_02_DIA/20260325_134643_StoreyAJ_20260323_02_DIA_Protein_Report.tsv"
peptide_path <- "C:/Users/aaron.storey/Desktop/ProteoViz_refactor/example_inputs/StoreyAJ_20260323_02_DIA/20260327_025131_StoreyAJ_20260323_02_DIA_Peptide_Report.tsv"
sample_path  <- "C:/Users/aaron.storey/Desktop/ProteoViz_refactor/example_inputs/StoreyAJ_20260323_02_DIA/StoreyAJ_20260323_02_DIA_Sample_name_table.tsv"

# Load protein data (raw wide-form) and the sample name table
protein_loaded <- load_protein_data(
  protein_path = protein_path,
  sample_path  = sample_path,
  type         = data_type
)

protein_raw  <- protein_loaded$protein_data   # long-form, already cleaned
sample_table <- protein_loaded$sample_table

# Load peptide data
peptide_loaded <- load_peptide_data(
  peptide_path = peptide_path,
  sample_path  = sample_path,
  type         = data_type
)

peptide_raw <- peptide_loaded$peptide_data    # long-form, already cleaned

message("Proteins loaded: ", length(unique(protein_raw$id)))
message("Peptides loaded: ", length(unique(peptide_raw$id)))
message("Samples in table: ", nrow(sample_table))


# -----------------------------------------------------------------------------
# Section 3: Remove samples
# -----------------------------------------------------------------------------
# List any Sample_name values to drop.  Leave as character(0) to keep all.

samples_to_remove <- character(0)
# Example: samples_to_remove <- c("Sample_03", "Sample_07")

groups_to_remove <- character(0)
# Example: groups_to_remove <- c("QC")

# Filter sample table
sample_table_filtered <- sample_table |>
  dplyr::filter(
    !Sample_name %in% samples_to_remove,
    !Group       %in% groups_to_remove
  )

# Drop those samples from protein and peptide long-form data
protein_removed <- remove_samples(
  protein_data      = protein_raw,
  sample_table      = sample_table_filtered,
  samples_to_remove = samples_to_remove
)

peptide_removed <- remove_samples(
  protein_data      = peptide_raw,
  sample_table      = sample_table_filtered,
  samples_to_remove = samples_to_remove
)

protein_clean    <- protein_removed$protein_data
peptide_clean    <- peptide_removed$protein_data
sample_table_use <- protein_removed$sample_table

message("Samples retained after exclusion: ", nrow(sample_table_use))


# -----------------------------------------------------------------------------
# Section 4: Build observation requirements
# -----------------------------------------------------------------------------
# group_counts lists each group and how many samples it contains.
# For each group, set the minimum number of valid observations a protein must
# have in that group to be retained.  The simplest approach is a single
# numeric value applied uniformly across all groups (e.g. 2 out of n).

group_counts <- sample_table_use |>
  dplyr::count(Group, name = "Number_of_samples") |>
  dplyr::ungroup()

# Uniform per-group minimum — adjust as needed
per_group_min <- 2L

group_obs_table <- make_observation_required_table(group_counts, per_group_min)

# Global requirement: a protein must have >= global_min_obs observations in
# >= global_min_groups groups (set both to 0 to apply only the per-group rule)
global_min_obs    <- 0L
global_min_groups <- 0L

message("Groups detected: ", paste(group_counts$Group, collapse = ", "))
print(group_obs_table)


# -----------------------------------------------------------------------------
# Section 5: Filter missing values
# -----------------------------------------------------------------------------

protein_filtered <- filter_missing_values(
  protein_data  = protein_clean,
  group_table   = group_obs_table,
  min_per_group = global_min_obs,
  min_groups    = global_min_groups
)

peptide_filtered <- filter_missing_values(
  protein_data  = peptide_clean,
  group_table   = group_obs_table,
  min_per_group = global_min_obs,
  min_groups    = global_min_groups
)

message(
  "Proteins after missing-value filter: ",
  length(unique(protein_filtered$id)),
  " (removed ",
  length(unique(protein_clean$id)) - length(unique(protein_filtered$id)),
  ")"
)
message(
  "Peptides after missing-value filter: ",
  length(unique(peptide_filtered$id))
)


# -----------------------------------------------------------------------------
# Section 6: Normalise
# -----------------------------------------------------------------------------

protein_normalised <- normalize_intensities(
  protein_data = protein_filtered,
  method       = norm_method
)

peptide_normalised <- normalize_intensities(
  protein_data = peptide_filtered,
  method       = norm_method
)

message("Normalisation method applied: ", norm_method)


# -----------------------------------------------------------------------------
# Section 7: QC plots
# -----------------------------------------------------------------------------

# -- 7a. Pre-normalisation intensity distributions ---------------------------
p_pre <- plot_intensity_distributions(
  data     = protein_clean,
  metadata = sample_table_use,
  color_by = "Group"
)

ggplot2::ggsave(
  filename = file.path(project_name, paste0(project_name, "_pre_norm_boxplot.tiff")),
  plot     = p_pre,
  width    = 14,
  height   = 6,
  dpi      = 300
)

# -- 7b. Post-normalisation intensity distributions --------------------------
p_post <- plot_intensity_distributions(
  data     = protein_normalised,
  metadata = sample_table_use,
  color_by = "Group"
)

ggplot2::ggsave(
  filename = file.path(project_name, paste0(project_name, "_post_norm_boxplot.tiff")),
  plot     = p_post,
  width    = 14,
  height   = 6,
  dpi      = 300
)

# -- 7c. PCA (interactive; saved as HTML) ------------------------------------
pca_plt <- plot_pca(
  data     = protein_normalised,
  metadata = sample_table_use
)

htmlwidgets::saveWidget(
  widget = pca_plt,
  file   = file.path(project_name, paste0(project_name, "_pca.html")),
  selfcontained = TRUE
)

# -- 7d. Correlation heatmap --------------------------------------------------
cor_plt <- plot_correlation_heatmap(
  data      = protein_normalised,
  metadata  = sample_table_use,
  scale_rows = FALSE
)

tiff(
  filename = file.path(project_name, paste0(project_name, "_correlation_heatmap.tiff")),
  width    = 10,
  height   = 9,
  units    = "in",
  res      = 300
)
print(cor_plt)
dev.off()

# -- 7e. Missing value pattern (pre-filter, on cleaned data) -----------------
na_plt <- plot_missing_values(
  data     = protein_clean,
  metadata = sample_table_use
)

tiff(
  filename = file.path(project_name, paste0(project_name, "_missing_values.tiff")),
  width    = 10,
  height   = 8,
  units    = "in",
  res      = 300
)
print(na_plt)
dev.off()

message("QC plots saved to: ", project_name)


# -----------------------------------------------------------------------------
# Section 8: Save outputs
# -----------------------------------------------------------------------------

file_prefix <- project_name

# -- Protein quantitative table (wide: proteins x samples) -------------------
protein_normalised |>
  dplyr::select(id, Sample_name, Intensity) |>
  dplyr::mutate(
    Sample_name = factor(Sample_name, levels = unique(sample_table_use$Sample_name))
  ) |>
  tidyr::pivot_wider(names_from = Sample_name, values_from = Intensity) |>
  dplyr::mutate(id = as.integer(id)) |>
  dplyr::arrange(id) |>
  readr::write_tsv(
    file.path(project_name, paste0(file_prefix, "_protein_quantitative_data.tsv"))
  )

# -- Peptide quantitative table ----------------------------------------------
peptide_normalised |>
  dplyr::select(id, Sample_name, Intensity) |>
  dplyr::mutate(
    Sample_name = factor(Sample_name, levels = unique(sample_table_use$Sample_name))
  ) |>
  tidyr::pivot_wider(names_from = Sample_name, values_from = Intensity) |>
  dplyr::mutate(id = as.integer(id)) |>
  dplyr::arrange(id) |>
  readr::write_tsv(
    file.path(project_name, paste0(file_prefix, "_peptide_quantitative_data.tsv"))
  )

# -- Sample design table (for App2 / Limma) ----------------------------------
sample_table_use |>
  dplyr::distinct(Sample_name, Group, Batch) |>
  dplyr::mutate(Block = "1") |>
  readr::write_tsv(
    file.path(project_name, paste0(file_prefix, "_design_table.tsv"))
  )

# -- Contrast list -----------------------------------------------------------
groups_for_contrasts <- unique(sample_table_use$Group)
groups_for_contrasts <- groups_for_contrasts[groups_for_contrasts != "Pool"]

expand.grid(x1 = groups_for_contrasts, y1 = groups_for_contrasts) |>
  dplyr::mutate(Contrast_name = paste0(y1, " - ", x1)) |>
  dplyr::filter(x1 != y1) |>
  dplyr::select(Contrast_name) |>
  readr::write_tsv(
    file.path(project_name, paste0(file_prefix, "_Sample_contrast_list.tsv"))
  )

message("All outputs saved to: ", project_name)
message("Workflow complete.")
