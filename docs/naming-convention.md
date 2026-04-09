# ProteoViz Naming Convention Reference

Generated during refactor audit — 2026-04-09. Applies to all new and refactored code.

---

## 1. Function Naming

**Pattern:** `verb_noun` in `snake_case`

- Start with an active verb that describes what the function does
- Follow with a noun describing what it acts on
- All lowercase, words separated by underscores
- No camelCase

| Verb | Use for |
|------|---------|
| `load_` | Reading data from files into R |
| `build_` | Constructing matrices or structured objects (design matrix, contrast matrix) |
| `run_` | Executing a multi-step statistical pipeline |
| `plot_` | Returning a plot object (ggplot2 or plotly) |
| `normalize_` | Applying a normalization transformation |
| `filter_` | Subsetting rows/columns based on criteria |
| `write_` | Saving data to a file |
| `extract_` | Pulling a subset of results from a larger object |
| `make_` | Avoid — use a more specific verb from the list above |

---

## 2. Variable Naming

**Pattern:** `snake_case`, descriptive noun phrases

- No single letters or numbers alone (`x1`, `df2`, `P7`)
- No generic sequences (`Box1`, `Box2`, `Box3`)
- Temporary data pipeline steps: suffix with `_clean`, `_filtered`, `_normalized`, `_matrix` to indicate stage
- Shiny UI layout objects: suffix with `_box`, `_tab`, `_panel`, `_row`

### Data pipeline variable convention (App server code)

| Stage | Suffix | Example |
|-------|--------|---------|
| Raw import | *(none)* | `protein_data` |
| After cleaning | `_clean` | `protein_data_clean` |
| After filtering | `_filtered` | `protein_data_filtered` |
| As matrix | `_matrix` | `protein_data_matrix` |
| After normalization | `_normalized` | `protein_data_normalized` |
| Back to long form | `_long` | `protein_data_long` |

---

## 3. File Naming

**Pattern:** `role_scope.R` in `snake_case`

| Role prefix | Purpose | Examples |
|-------------|---------|---------|
| `core_` | Pure R functions, no Shiny dependencies | `core_data_loading.R`, `core_normalization.R`, `core_filtering.R`, `core_limma.R`, `core_plotting.R` |
| `app1_` | App1-specific Shiny UI or server | `app1_ui.R`, `app1_server.R` |
| `app2_` | App2-specific Shiny UI or server | `app2_ui.R`, `app2_server.R` |
| `run_` | Top-level launcher functions | `run_app1.R`, `run_app2.R` |

Files in `inst/scripts/` use descriptive names: `run_app1_workflow.R`, `run_app2_workflow.R`.

---

## 4. Before / After Examples

### Functions

| Before | After | Reason |
|--------|-------|--------|
| `readSamplesReport()` | `load_samples_report()` | `load_` for file I/O; snake_case |
| `makeDesignTable()` | `build_design_matrix()` | `build_` for constructing structured objects; more precise noun |
| `makeContrastTable()` | `build_contrast_matrix()` | Same pattern; aligns with limma terminology |
| `runLimma()` | `run_limma()` | Already correct verb; just snake_case |
| `fitContrasts()` | `fit_contrasts()` | snake_case; verb already correct |
| `combineLimmaData()` | `combine_limma_results()` | snake_case; `results` clearer than `data` |
| `makeProteinHeatHeatmap()` | `plot_protein_heatmap()` | `plot_` for functions returning plot objects; remove redundant `Heat` |
| `makePeptideHeatHeatmap()` | `plot_peptide_heatmap()` | Same |
| `PCA_plotly()` | `plot_pca()` | `plot_` prefix; implementation detail (plotly) not in name |
| `make_boxplot()` | `plot_intensity_boxplot()` | More descriptive noun; consistent `plot_` prefix |

### UI Variables (Shiny server code)

| Before | File | After |
|--------|------|-------|
| `Box1` | App1 | `import_files_box` |
| `Box3` | App1 | `protein_qc_controls_box` |
| `Box4` | App1 | `group_observation_criteria_box` |
| `Box4_5` | App1 | `global_observation_criteria_box` |
| `Box5` | App1 | `normalization_controls_box` |
| `Box6` | App1 | `protein_export_box` |
| `Box7` | App1 | `peptide_export_box` |
| `Box0_1` | App2 | `file_upload_box` |
| `Box2` | App2 | `experimental_design_box` |
| `Box5` | App2 | `limma_plot_box` |
| `BoxP8` | App2 | `upset_diagram_box` |
| `BoxP9` | App2 | `set_selection_box` |
| `Box6` | App2 | `data_export_box` |

### Shiny reactive variables (server pipeline)

| Before | After |
|--------|-------|
| `protein_df` | `protein_data` |
| `protein_df2` | `protein_data_clean` |
| `protein_df3` | `protein_data_filtered` |
| `protein_df4` | `protein_data_matrix` |
| `protein_df5` | `protein_data_normalized` |
| `protein_df6` | `protein_data_long` |
| `sampleNameTable` | `sample_table` |
| `sampleNameTable2` | `sample_table_df` |
| `sampleNameTable3` | `sample_table_filtered` |
| `groupTable` | `group_counts` |
| `groupTable2` | `group_observation_requirements` |
| `design_table2` | `design_matrix_df` |
| `contrast_list` | `contrast_table_df` |

---

## 5. Rules Summary

1. All names: `snake_case`, no camelCase
2. Functions: start with an active verb from the approved list
3. Files: `role_scope.R` pattern
4. UI box variables: end with `_box`, `_tab`, or `_panel`
5. Data pipeline reactives: noun + stage suffix (`_clean`, `_filtered`, `_matrix`, `_normalized`, `_long`)
6. No numbered suffixes (`df2`, `Box3`) — use descriptive stage or purpose names
