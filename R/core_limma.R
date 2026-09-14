#' @title Build a limma design matrix
#'
#' @description Constructs a limma-compatible design matrix from a sample table.
#'   Groups are encoded as a factor (preserving the order in which they first
#'   appear in the table) and the intercept is dropped so that each column
#'   represents one group mean — the standard \code{~0 + Groups}
#'   parameterisation used with contrast matrices in limma.
#'
#' @param sample_table A data frame with at minimum two columns:
#'   \describe{
#'     \item{Sample_name}{Character. Unique sample identifiers used as row
#'       names of the returned matrix.}
#'     \item{Group}{Character or factor. Group membership of each sample.
#'       The order of levels is taken from the first occurrence of each value
#'       in the column.}
#'   }
#'   Additional columns (e.g. \code{Batch}, \code{Block}) are accepted but
#'   are not used by this function.
#'
#' @return A numeric matrix with one row per sample and one column per group.
#'   Row names are set to \code{sample_table$Sample_name} and column names are
#'   set to the group levels.  The matrix is directly compatible with
#'   \code{limma::lmFit()}.  Returns \code{NULL} if fewer than two distinct
#'   groups are present.
#'
#' @export
build_design_matrix <- function(sample_table) {
  sample_table_clean <- sample_table |>
    dplyr::mutate(Group = factor(Group, levels = unique(Group)))

  groups <- sample_table_clean$Group

  if (length(unique(groups)) < 2) return(NULL)

  design_matrix <- model.matrix(~0 + groups)

  rownames(design_matrix) <- sample_table_clean$Sample_name
  colnames(design_matrix) <- levels(groups)

  design_matrix
}


#' @title Build a limma contrast matrix
#'
#' @description Constructs a limma-compatible contrast matrix from a set of
#'   contrast strings and an existing design matrix.  Each contrast string
#'   should be an arithmetic expression over the column names of
#'   \code{design_matrix}, e.g. \code{"A - B"} or
#'   \code{"Treatment - Control"}.
#'
#' @param design_matrix A numeric matrix as returned by
#'   \code{\link{build_design_matrix}}.  Its column names define the valid group
#'   names that may appear in contrast expressions.
#' @param contrasts A data frame with a column named \code{Contrast_name}
#'   containing contrast strings, \strong{or} a character vector of contrast
#'   strings directly.  Blank/empty strings are silently ignored.
#'
#' @return A numeric contrast matrix compatible with
#'   \code{limma::contrasts.fit()}.  Returns \code{NULL} if no non-empty
#'   contrast strings are found.
#'
#' @export
build_contrast_matrix <- function(design_matrix, contrasts) {
  if (is.data.frame(contrasts)) {
    contrast_strings <- contrasts$Contrast_name
  } else {
    contrast_strings <- as.character(contrasts)
  }

  contrast_strings_clean <- contrast_strings[nchar(trimws(contrast_strings)) > 0]

  if (length(contrast_strings_clean) < 1) return(NULL)

  limma::makeContrasts(
    contrasts = as.list(contrast_strings_clean),
    levels    = colnames(design_matrix)
  )
}


#' @title Fit a limma linear model
#'
#' @description Fits a linear model to a numeric protein-by-sample matrix using
#'   \code{limma::lmFit()}.  When \code{block = TRUE} the function first
#'   estimates the intra-block correlation via
#'   \code{limma::duplicateCorrelation()} using the \code{Block} column of
#'   \code{sample_table}, then passes that consensus correlation to
#'   \code{lmFit()} alongside the block vector.
#'
#' @param data_matrix A numeric matrix with proteins in rows and samples in
#'   columns.  Column names must match \code{sample_table$Sample_name}; columns
#'   are reordered automatically.
#' @param design_matrix A numeric design matrix as returned by
#'   \code{\link{build_design_matrix}}.
#' @param sample_table A data frame with a \code{Sample_name} column.  When
#'   \code{block = TRUE} a \code{Block} column must also be present.
#' @param block Logical (default \code{FALSE}).  When \code{TRUE}, applies
#'   duplicate-correlation blocking using \code{sample_table$Block}.
#'
#' @return A fitted \code{MArrayLM} object ready to be passed to
#'   \code{\link{fit_contrasts}}.
#'
#' @export
run_limma <- function(data_matrix, design_matrix, sample_table, block = FALSE) {
  quant_data <- data_matrix[, sample_table$Sample_name]

  if (block) {
    sample_table <- sample_table |>
      dplyr::mutate(Block = factor(Block, levels = unique(Block)))

    corfit <- limma::duplicateCorrelation(
      quant_data,
      design_matrix,
      block = sample_table$Block
    )

    limma::lmFit(
      quant_data,
      design_matrix,
      block       = sample_table$Block,
      correlation = corfit$consensus.correlation
    )
  } else {
    limma::lmFit(quant_data, design_matrix)
  }
}


#' @title Fit contrasts and apply empirical Bayes moderation
#'
#' @description Applies a contrast matrix to a fitted limma model using
#'   \code{limma::contrasts.fit()}, then applies empirical Bayes variance
#'   moderation via \code{limma::eBayes()}.
#'
#' @param limma_fit A fitted \code{MArrayLM} object as returned by
#'   \code{\link{run_limma}}.
#' @param contrast_matrix A contrast matrix as returned by
#'   \code{\link{build_contrast_matrix}}.
#'
#' @return A moderated \code{MArrayLM} object containing t-statistics,
#'   p-values, and log-odds for each contrast.
#'
#' @export
fit_contrasts <- function(limma_fit, contrast_matrix) {
  limma_fit |>
    limma::contrasts.fit(contrast_matrix) |>
    limma::eBayes()
}


# ---------------------------------------------------------------------------
# Helper: extract results for a single contrast
# ---------------------------------------------------------------------------
extract_limma_results <- function(contrast_fit, comparison) {
  limma::topTable(contrast_fit, coef = comparison, number = Inf) |>
    tibble::rownames_to_column("protein") |>
    dplyr::mutate(Comparison = comparison) |>
    tibble::as_tibble()
}


# ---------------------------------------------------------------------------
# Helper: combine results across all contrasts into one long tibble
# ---------------------------------------------------------------------------
combine_limma_results <- function(contrast_fit, design_matrix, contrast_matrix) {
  comparisons <- colnames(contrast_matrix)

  purrr::map(comparisons, \(cmp) extract_limma_results(contrast_fit, cmp)) |>
    purrr::list_rbind() |>
    dplyr::filter(!is.na(adj.P.Val))
}


# ---------------------------------------------------------------------------
# Helper: pivot decideTests output to a tidy long tibble with Up/Down labels
# ---------------------------------------------------------------------------
make_tidy_test_results <- function(test_results) {
  test_results@.Data |>
    tibble::as_tibble(rownames = "protein") |>
    tidyr::pivot_longer(
      cols      = -protein,
      names_to  = "Test",
      values_to = "Result"
    ) |>
    dplyr::filter(Result %in% c(1L, -1L)) |>
    dplyr::mutate(
      Result = dplyr::if_else(Result == 1L, "Up", "Down")
    ) |>
    dplyr::group_by(Test, Result) |>
    tidyr::nest() |>
    tidyr::unite("Test", Test, Result, sep = " ")
}


# ---------------------------------------------------------------------------
# Helper: convert the nested tidy tibble to a named list of protein id vectors
# ---------------------------------------------------------------------------
make_test_results_list <- function(tidy_test_results) {
  protein_lists <- purrr::map(tidy_test_results$data, "protein")
  stats::setNames(protein_lists, tidy_test_results$Test)
}


# ---------------------------------------------------------------------------
# Exported: find proteins exclusively detected in one side of a contrast
# ---------------------------------------------------------------------------

#' @title Find proteins exclusively detected on one side of a contrast
#'
#' @description Identifies proteins for which every sample on one side of a
#'   contrast is missing (\code{NA}) while at least one sample on the other
#'   side has a value. These are the proteins responsible for limma's
#'   \code{"Partial NA coefficients"} warning — a moderated t-test is
#'   undefined for them because one group has no quantitative values at all.
#'   Rather than forcing a p-value, this function reports detection-based
#'   statistics (mean intensity and number of replicates detected) so these
#'   proteins can be ranked and visualised separately from the standard
#'   volcano plot.
#'
#' @param data_matrix A numeric matrix with proteins in rows and samples in
#'   columns, as passed to \code{\link{run_limma}}. Row names are used as
#'   protein identifiers.
#' @param sample_table A data frame with \code{Sample_name} and \code{Group}
#'   columns (as used by \code{\link{build_design_matrix}}).
#' @param comparison A single contrast string of the form
#'   \code{"GroupA - GroupB"}, matching the group names in
#'   \code{sample_table$Group}.
#' @param peptide_counts An optional data frame with columns \code{protein}
#'   and \code{n_peptides} to join onto the result. Defaults to \code{NULL}.
#'
#' @return A tibble with one row per exclusively-detected protein and the
#'   columns \code{protein}, \code{Comparison}, \code{contrast_side}
#'   (\code{"A"} if detected in the first group named in \code{comparison},
#'   \code{"B"} if the second), \code{exclusive_group} (the group the protein
#'   was detected in), \code{other_group}, \code{mean_intensity} (mean of the
#'   non-missing values), \code{n_detected}, and \code{n_total_in_group}.
#'   Proteins detected (or missing) on both sides of the contrast are
#'   excluded. Returns an empty tibble if \code{comparison} cannot be split
#'   into two group names.
#'
#' @examples
#' \dontrun{
#' exclusive <- find_exclusive_proteins(
#'   data_matrix    = limma_input_matrix,
#'   sample_table   = sample_table,
#'   comparison     = "S20_15105 - Pregnancy",
#'   peptide_counts = peptide_counts
#' )
#' }
#'
#' @export
find_exclusive_proteins <- function(data_matrix, sample_table, comparison,
                                    peptide_counts = NULL) {
  parts <- trimws(strsplit(comparison, " - ", fixed = TRUE)[[1]])
  if (length(parts) != 2) return(tibble::tibble())

  group_a <- parts[1]
  group_b <- parts[2]

  samples_a <- sample_table$Sample_name[sample_table$Group == group_a]
  samples_b <- sample_table$Sample_name[sample_table$Group == group_b]

  mat_a <- data_matrix[, colnames(data_matrix) %in% samples_a, drop = FALSE]
  mat_b <- data_matrix[, colnames(data_matrix) %in% samples_b, drop = FALSE]

  n_detected_a <- rowSums(!is.na(mat_a))
  n_detected_b <- rowSums(!is.na(mat_b))

  exclusive_a <- n_detected_a > 0 & n_detected_b == 0
  exclusive_b <- n_detected_b > 0 & n_detected_a == 0

  mean_a <- rowMeans(mat_a, na.rm = TRUE)
  mean_b <- rowMeans(mat_b, na.rm = TRUE)

  result <- tibble::tibble(
    protein          = rownames(data_matrix),
    Comparison       = comparison,
    contrast_side    = dplyr::case_when(
      exclusive_a ~ "A",
      exclusive_b ~ "B",
      TRUE        ~ NA_character_
    ),
    exclusive_group  = dplyr::case_when(
      exclusive_a ~ group_a,
      exclusive_b ~ group_b,
      TRUE        ~ NA_character_
    ),
    other_group      = dplyr::case_when(
      exclusive_a ~ group_b,
      exclusive_b ~ group_a,
      TRUE        ~ NA_character_
    ),
    mean_intensity   = dplyr::case_when(
      exclusive_a ~ mean_a,
      exclusive_b ~ mean_b,
      TRUE        ~ NA_real_
    ),
    n_detected       = dplyr::case_when(
      exclusive_a ~ n_detected_a,
      exclusive_b ~ n_detected_b,
      TRUE        ~ NA_integer_
    ),
    n_total_in_group = dplyr::case_when(
      exclusive_a ~ ncol(mat_a),
      exclusive_b ~ ncol(mat_b),
      TRUE        ~ NA_integer_
    )
  ) |>
    dplyr::filter(!is.na(exclusive_group))

  if (!is.null(peptide_counts)) {
    result <- dplyr::left_join(result, peptide_counts, by = "protein")
  }

  result
}


# ---------------------------------------------------------------------------
# Exported: write a wide-format summary TSV joining quant + limma results
# ---------------------------------------------------------------------------

#' @title Write a limma summary TSV
#'
#' @description Joins quantitative data and metadata with limma differential
#'   expression results (spread to wide format: one column per comparison per
#'   statistic) and writes the result to a tab-separated file named
#'   \code{<prefix>_Summarized_output.tsv}.
#'
#' @param quant_data A data frame of quantitative intensities that contains an
#'   \code{id} column matching the row identifiers used in the limma analysis.
#' @param combined_results A long-format tibble as returned by
#'   \code{combine_limma_results}, with columns \code{protein},
#'   \code{logFC}, \code{P.Value}, \code{adj.P.Val}, and \code{Comparison}.
#' @param metadata A data frame with an \code{id} column and at minimum protein
#'   annotation columns to include in the output.
#' @param prefix A character string prepended to the output file name, e.g.
#'   \code{"Experiment1"} produces \code{Experiment1_Summarized_output.tsv}.
#'
#' @return The joined wide-format data frame, invisibly.  The primary side
#'   effect is writing the TSV file.
#'
#' @export
write_limma_summary <- function(quant_data, combined_results, metadata, prefix) {
  spread_limma <- combined_results |>
    dplyr::select(protein, logFC, P.Value, adj.P.Val, Comparison) |>
    tidyr::pivot_longer(
      cols      = c(logFC, P.Value, adj.P.Val),
      names_to  = "Type",
      values_to = "Value"
    ) |>
    tidyr::unite("Type", Comparison, Type, sep = " ") |>
    tidyr::pivot_wider(names_from = Type, values_from = Value)

  out <- metadata |>
    dplyr::right_join(quant_data, by = "id") |>
    dplyr::left_join(spread_limma, by = c("id" = "protein"))

  readr::write_tsv(out, paste0(prefix, "_Summarized_output.tsv"))
  invisible(out)
}


# ---------------------------------------------------------------------------
# Exported: high-level convenience wrapper — the main user-facing function
# ---------------------------------------------------------------------------

#' @title Run a complete limma differential expression analysis on protein data
#'
#' @description Convenience wrapper that chains \code{\link{run_limma}},
#'   \code{\link{fit_contrasts}}, and \code{combine_limma_results} into
#'   a single call, then splits the combined results into a named list of
#'   per-contrast tidy data frames ready for downstream use or export.
#'
#' @param data A numeric matrix with proteins in rows and samples in columns.
#'   Column names must match \code{sample_table$Sample_name}.
#' @param design_matrix A numeric design matrix as returned by
#'   \code{\link{build_design_matrix}}.
#' @param contrast_matrix A contrast matrix as returned by
#'   \code{\link{build_contrast_matrix}}.
#' @param sample_table A data frame with a \code{Sample_name} column.  When
#'   \code{block = TRUE} a \code{Block} column must also be present.
#' @param block Logical (default \code{FALSE}).  When \code{TRUE}, applies
#'   duplicate-correlation blocking via \code{limma::duplicateCorrelation()}.
#'
#' @return A named list of tibbles, one per contrast (names taken from
#'   \code{colnames(contrast_matrix)}).  Each tibble contains the columns:
#'   \describe{
#'     \item{protein}{Row identifier from the input matrix (character).}
#'     \item{logFC}{Log2 fold change.}
#'     \item{AveExpr}{Average log2 expression across all samples.}
#'     \item{t}{Moderated t-statistic.}
#'     \item{P.Value}{Raw p-value.}
#'     \item{adj.P.Val}{Benjamini-Hochberg adjusted p-value.}
#'   }
#'   Rows with \code{NA} adjusted p-values (proteins with insufficient data for
#'   a given contrast) are dropped.
#'
#' @examples
#' \dontrun{
#' design   <- build_design_matrix(sample_table)
#' contrasts <- build_contrast_matrix(design, c("GroupA - GroupB"))
#' results  <- run_limma_protein(
#'   data            = protein_matrix,
#'   design_matrix   = design,
#'   contrast_matrix = contrasts,
#'   sample_table    = sample_table
#' )
#' # results[["GroupA - GroupB"]] is a tibble of per-protein statistics
#' }
#'
#' @export
run_limma_protein <- function(data, design_matrix, contrast_matrix,
                              sample_table, block = FALSE) {
  contrast_fit <- run_limma(data, design_matrix, sample_table, block = block) |>
    fit_contrasts(contrast_matrix)

  combined <- combine_limma_results(contrast_fit, design_matrix, contrast_matrix)

  keep_cols <- c("protein", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")

  comparisons <- colnames(contrast_matrix)
  results <- purrr::map(
    comparisons,
    \(cmp) combined |>
      dplyr::filter(Comparison == cmp) |>
      dplyr::select(dplyr::any_of(keep_cols))
  )
  stats::setNames(results, comparisons)
}
