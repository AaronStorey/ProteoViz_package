# ProteoViz

ProteoViz is an R package for analyzing quantitative proteomics data (DIA/EncyclopeDIA, Spectronaut, and TMT workflows). It provides two interactive Shiny applications plus a set of scriptable core functions, all built around the same three-stage workflow:

1. **Clean** — load raw quantitative and metadata tables, filter low-confidence or sparsely-detected peptides/proteins, and normalize intensities across samples (e.g. median, quantile, or VSN normalization via `NormalyzerDE`/`preprocessCore`).
2. **Test for differential abundance** — build a design matrix from the experimental groups and contrasts of interest, then fit a moderated t-test using `limma` (with optional blocking/duplicate-correlation for repeated-measures designs) to estimate log fold-changes and adjusted p-values for each contrast.
3. **Visualize** — summarize and explore the results interactively: volcano plots, heatmaps, PCA, UpSet plots, and gene ontology enrichment, plus export of a combined summary table joining metadata, quantitative values, and statistical results.

## Applications

- **App 1** (`Run_app1*_shiny.R`) — data cleaning, filtering, and normalization for DIA, Spectronaut, and TMT inputs.
- **App 2** (`Run_app2_shiny.R`) — statistical testing (limma) and downstream visualization of an already-cleaned dataset.

The underlying functions (`core_data_loading.R`, `core_data_cleaning.R`, `core_filtering.R`, `core_normalization.R`, `core_limma.R`, `core_plotting.R`, `core_enrichment.R`) can also be called directly for reproducible, scripted analyses outside of the Shiny interface.
