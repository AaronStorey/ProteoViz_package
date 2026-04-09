# core_data_cleaning.R
# Pure R functions for cleaning and reshaping protein quantitative data.
# No Shiny dependencies.


# Exported functions -----------------------------------------------------------

#' Clean and reshape protein quantitative data
#'
#' @title Clean and reshape wide-form protein data into long form
#'
#' @description Takes a wide-form protein (or peptide) tibble as produced by
#'   the internal readers in \code{core_data_loading.R}, selects the sample
#'   columns identified in \code{sample_table}, pivots to long form, optionally
#'   applies a log2 transformation to positive intensity values, optionally
#'   removes zero and missing observations, and left-joins the full sample
#'   metadata so that every row carries group, batch, and other annotation
#'   columns.
#'
#'   The \code{normToReferenceChannel} parameter found in the original
#'   \code{clean_protein_df()} is intentionally omitted — that branch is
#'   TMT-specific and is out of scope for this module.
#'
#' @param protein_raw A wide-form tibble with an \code{id} column and one
#'   column per sample whose names match the \code{Data_name} column of
#'   \code{sample_table}.  Typically the output of \code{read_quant_report()}
#'   or \code{read_samples_report()}.
#' @param sample_table A tibble describing sample metadata. Must contain at
#'   minimum two columns: \code{Data_name} (values matching the sample column
#'   names of \code{protein_raw}) and \code{Sample_name} (the human-readable
#'   label used for downstream analysis). Additional columns (e.g.
#'   \code{Group}, \code{Batch}, \code{Reference}) are preserved and attached
#'   to every row of the returned tibble.
#' @param log2_transform Logical. If \code{TRUE} (default), intensity values
#'   greater than zero are replaced by their log2. Values that are zero or
#'   negative are left unchanged; they will be removed by the
#'   \code{remove_zeros} step if that is also \code{TRUE}.
#' @param remove_zeros Logical. If \code{TRUE} (default), rows whose
#'   \code{Intensity} is not strictly greater than zero (i.e. zeros, negative
#'   values, and \code{NA}s) are dropped. When both \code{log2_transform} and
#'   \code{remove_zeros} are \code{TRUE} the effect is: transform positive
#'   values to log2 scale, then discard everything that did not transform.
#'
#' @return A long-form tibble with one row per protein-sample observation. The
#'   columns are:
#'   \describe{
#'     \item{\code{id}}{Protein/peptide identifier, carried over from
#'       \code{protein_raw}.}
#'     \item{\code{Sample_name}}{Human-readable sample label from
#'       \code{sample_table}.}
#'     \item{\code{Intensity}}{Quantitative intensity value, log2-transformed
#'       if \code{log2_transform = TRUE}.}
#'     \item{additional columns}{All other columns present in
#'       \code{sample_table} (e.g. \code{Data_name}, \code{Group},
#'       \code{Batch}, \code{Reference}).}
#'   }
#'
#' @export
clean_protein_data <- function(protein_raw,
                               sample_table,
                               log2_transform = TRUE,
                               remove_zeros   = TRUE) {

  # Resolve the Data_name column regardless of what the first column is called
  # in the supplied sample_table (mirrors the original defensive rename).
  colnames(sample_table)[1] <- "Data_name"
  data_names <- sample_table[["Data_name"]]

  # Pivot from wide to long, selecting only the sample columns present in
  # sample_table. This silently ignores any samples listed in sample_table that
  # are absent from protein_raw (dplyr::any_of behaviour), and it drops metadata
  # columns that are not sample columns.
  df <- protein_raw |>
    dplyr::select(id, dplyr::any_of(data_names)) |>
    tidyr::pivot_longer(
      cols      = dplyr::any_of(data_names),
      names_to  = "Data_name",
      values_to = "Intensity"
    ) |>
    dplyr::left_join(sample_table, by = "Data_name")

  if (log2_transform) {
    df <- df |>
      dplyr::mutate(
        Intensity = dplyr::if_else(Intensity > 0, log2(Intensity), Intensity)
      )
  }

  if (remove_zeros) {
    df <- df |>
      dplyr::filter(Intensity > 0)
  }

  df
}


# Internal helpers -------------------------------------------------------------

# protein_df_to_matrix ---------------------------------------------------------
# Converts a long-form protein tibble (columns: id, Sample_name, Intensity,
# ...) to a numeric matrix with proteins as rows and samples as columns.
# Used by normalization functions in core_normalization.R.
protein_df_to_matrix <- function(protein_data) {
  protein_data |>
    dplyr::select(id, Sample_name, Intensity) |>
    tidyr::pivot_wider(
      names_from  = Sample_name,
      values_from = Intensity
    ) |>
    as.data.frame() |>
    tibble::column_to_rownames("id") |>
    as.matrix()
}


# matrix_to_long_form ----------------------------------------------------------
# Converts a wide numeric matrix (rows=proteins, cols=samples) back to a
# long-form tibble with columns id, Sample_name, Intensity. Rows where
# Intensity is NA are removed.
# Used by normalization functions in core_normalization.R to return results in
# the standard long form expected by the rest of the pipeline.
matrix_to_long_form <- function(protein_matrix) {
  protein_matrix |>
    tibble::as_tibble(rownames = "id") |>
    tidyr::pivot_longer(
      cols      = -id,
      names_to  = "Sample_name",
      values_to = "Intensity"
    ) |>
    dplyr::filter(!is.na(Intensity))
}
