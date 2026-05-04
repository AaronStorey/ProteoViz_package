# core_data_loading.R
# Pure R functions for loading protein and peptide quantitative data.
# No Shiny dependencies.


# Exported functions -----------------------------------------------------------

#' Load protein quantitative data
#'
#' @title Load protein quantitative data from a report file
#'
#' @description Reads a protein report file and a sample name table TSV, then
#'   returns a named list containing the protein data in long form and the
#'   sample table. Supports Spectronaut and DIA (EncyclopeDIA) input formats.
#'
#' @param protein_path Character. Path to the protein report file.
#' @param sample_path Character. Path to the sample name table TSV. The file
#'   must contain the columns: \code{Data_name}, \code{Sample_name},
#'   \code{Type}, \code{Group}, \code{Batch}, \code{Reference}.
#' @param type Character. Input format: \code{"Spectronaut"} or \code{"DIA"}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{protein_data}}{A tibble in long form with columns
#'       \code{id}, \code{Sample_name}, \code{Intensity}, and any additional
#'       sample-table columns joined from \code{sample_path}.}
#'     \item{\code{sample_table}}{A tibble with the parsed sample name table.}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- load_protein_data(
#'   protein_path = "data/protein_report.tsv",
#'   sample_path  = "data/sample_names.tsv",
#'   type         = "Spectronaut"
#' )
#' protein_long <- result$protein_data
#' sample_table <- result$sample_table
#' }
#'
#' @export
load_protein_data <- function(protein_path, sample_path, type) {

  if (!type %in% c("Spectronaut", "DIA")) {
    stop("'type' must be \"Spectronaut\" or \"DIA\". Got: ", type)
  }

  sample_table <- read_sample_name_table(sample_path)

  protein_raw <- switch(
    type,
    Spectronaut = read_quant_report(protein_path, type),
    DIA         = read_samples_report(protein_path, type)
  )

  protein_data <- clean_protein_data(protein_raw, sample_table)

  list(
    protein_data = protein_data,
    sample_table = sample_table
  )
}


#' Load peptide quantitative data
#'
#' @title Load peptide quantitative data from a report file
#'
#' @description Reads a peptide report file and a sample name table TSV, then
#'   returns a named list containing the peptide data in long form and the
#'   sample table. Supports Spectronaut and DIA (EncyclopeDIA) input formats.
#'
#' @param peptide_path Character. Path to the peptide report file.
#' @param sample_path Character. Path to the sample name table TSV. The file
#'   must contain the columns: \code{Data_name}, \code{Sample_name},
#'   \code{Type}, \code{Group}, \code{Batch}, \code{Reference}.
#' @param type Character. Input format: \code{"Spectronaut"} or \code{"DIA"}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{peptide_data}}{A tibble in long form with columns
#'       \code{id}, \code{Sample_name}, \code{Intensity}, and any additional
#'       sample-table columns joined from \code{sample_path}.}
#'     \item{\code{sample_table}}{A tibble with the parsed sample name table.}
#'   }
#'
#' @export
load_peptide_data <- function(peptide_path, sample_path, type) {

  if (!type %in% c("Spectronaut", "DIA")) {
    stop("'type' must be \"Spectronaut\" or \"DIA\". Got: ", type)
  }

  sample_table <- read_sample_name_table(sample_path)

  peptide_raw <- switch(
    type,
    Spectronaut = read_peptide_spectronaut(peptide_path, type),
    DIA         = read_peptide_report(peptide_path, type)
  )

  peptide_data <- clean_protein_data(peptide_raw, sample_table)

  list(
    peptide_data = peptide_data,
    sample_table = sample_table
  )
}


# Internal helpers -------------------------------------------------------------

# read_samples_report ----------------------------------------------------------
# Reads a ScaffoldDIA / EncyclopeDIA samples report.
# The file has a multi-row preamble; the data begins at the row containing
# "Visible". Rows without a value in the "Visible" column and DECOY rows are
# discarded.
read_samples_report <- function(file_path, type) {
  if (type == "DIA") {
    first_read  <- readLines(file_path)
    file_length <- length(first_read)
    skip_number <- grep("Visible", first_read) - 1

    readr::read_csv(
      file_path,
      skip  = skip_number,
      na    = c("", "Missing Value"),
      n_max = file_length - skip_number - 2
    ) %>%
      purrr::set_names(gsub(" ", "_", names(.))) %>%
      dplyr::filter(!is.na(Visible)) %>%
      dplyr::rename(id = `#`) %>%
      dplyr::filter(!grepl("\\.", id)) %>%
      dplyr::filter(!grepl("DECOY", Protein_Name))
  } else {
    stop("read_samples_report() supports type = \"DIA\" only.")
  }
}


# read_quant_report ------------------------------------------------------------
# Reads a Spectronaut protein groups report (tab-separated).
# Selects metadata columns plus all sample-level raw.PG.Quantity columns,
# assigns a numeric row id, then strips the ".raw.PG.Quantity" suffix from
# sample column names so they equal the bare R.FileName value — the same as
# the column names produced by read_peptide_spectronaut().  This allows a
# single sample name table (Data_name = R.FileName) to be shared across both
# protein and peptide reports from the same Spectronaut run.
read_quant_report <- function(file_path, type) {
  df <- readr::read_tsv(file_path, guess_max = 10000) %>%
    purrr::set_names(gsub(" ", "_", names(.))) %>%
    dplyr::select(
      PG.ProteinGroups, PG.Genes, PG.ProteinDescriptions, PG.UniProtIds,
      dplyr::contains("raw.PG.Quantity")
    ) %>%
    dplyr::distinct() %>%
    tibble::rownames_to_column("id")

  # "StoreyAJ_20260324_01_DIA_01.raw.PG.Quantity" → "StoreyAJ_20260324_01_DIA_01"
  quant_cols <- grep("raw\\.PG\\.Quantity", names(df), value = TRUE)
  names(df)[match(quant_cols, names(df))] <- sub("\\.raw\\.PG\\.Quantity$", "", quant_cols)

  df
}


# read_peptide_report ----------------------------------------------------------
# Reads a DIA (EncyclopeDIA) peptide report.
# The file uses a two-row header: row 1 contains friendly column names, row 2
# contains data-type codes. Intensity columns are identified by the value
# "Quant. Intensity" in row 1.
read_peptide_report <- function(file_path, type) {
  header_row  <- read.csv(file_path, nrows = 1)
  quant_cols  <- names(header_row)[header_row == "Quant. Intensity"]
  col_indices <- c(1:6, grep("Intensity", header_row))
  col_names   <- c(as.character(header_row[1:6]), quant_cols)

  readr::read_csv(file_path, skip = 2, col_names = FALSE) %>%
    dplyr::select(dplyr::all_of(col_indices)) %>%
    purrr::set_names(col_names) %>%
    purrr::set_names(gsub(" ", "_", names(.))) %>%
    dplyr::mutate(id = as.double(rownames(.)))
}


# read_peptide_spectronaut -----------------------------------------------------
# Reads a Spectronaut peptide report (tab-separated).
# Renames EG.TotalQuantity to Intensity, removes NaN rows, and pivots to wide
# form so that each sample becomes a column named by its bare R.FileName value.
# Column names therefore match those produced by read_quant_report(), allowing
# a single sample name table to be shared across protein and peptide reports.
read_peptide_spectronaut <- function(file_path, type) {
  readr::read_tsv(file_path, guess_max = 10000) |>
    dplyr::rename(Intensity = `EG.TotalQuantity (Settings)`) |>
    dplyr::filter(!is.nan(Intensity)) |>
    dplyr::select(
      R.FileName, EG.PrecursorId,
      PG.ProteinAccessions, PG.ProteinDescriptions, PG.ProteinNames,
      Intensity
    ) |>
    tidyr::pivot_wider(names_from = R.FileName, values_from = Intensity) |>
    tibble::rowid_to_column("id")
}


# read_sample_name_table -------------------------------------------------------
# Reads a sample name table TSV. Forces the first column to be named
# "Data_name" regardless of what the file header says.
read_sample_name_table <- function(sample_file) {
  sample_table <- readr::read_tsv(sample_file)
  colnames(sample_table)[1] <- "Data_name"
  sample_table
}


# make_sample_name_table -------------------------------------------------------
# Auto-generates a minimal sample name table from the column names of a wide
# protein/peptide data frame. Useful when no hand-curated sample file exists.
# Supports Spectronaut and DIA formats.
#
# For Spectronaut, tries protein-style columns (raw.PG.Quantity) first, then
# falls back to peptide-style columns (raw_<FileName>) when the data came from
# a peptide report pivot.  Group information is NOT available from column names
# alone; use extract_spectronaut_sample_table() when the original long-format
# file is available to get real Group values from R.Condition.
make_sample_name_table <- function(protein_df, type) {
  if (type == "Spectronaut") {
    # Sample columns are everything except the known metadata columns.
    # Works for both protein reports (read_quant_report strips the
    # .raw.PG.Quantity suffix) and peptide reports (read_peptide_spectronaut
    # uses bare R.FileName values), so the resulting Data_name values are
    # identical across both report types.
    meta_cols <- c("id", "PG.ProteinGroups", "PG.Genes",
                   "PG.ProteinDescriptions", "PG.UniProtIds",
                   "EG.PrecursorId", "PG.ProteinAccessions", "PG.ProteinNames")
    data_cols <- setdiff(names(protein_df), meta_cols)

    return(
      tibble::tibble(
        Data_name   = data_cols,
        Sample_name = data_cols,
        Type        = "Lysate",
        Group       = "A",
        Batch       = 1L,
        Reference   = "N"
      )
    )
  } else if (type == "DIA") {
    col_indices <- grep("\\.", names(protein_df))
  } else {
    stop("make_sample_name_table() supports type = \"Spectronaut\" or \"DIA\" only.")
  }

  tibble::tibble(
    Data_name   = names(protein_df)[col_indices],
    Sample_name = paste0("Sample_", seq_along(col_indices)),
    Type        = "Lysate",
    Group       = "A",
    Batch       = 1L,
    Reference   = "N"
  )
}


# extract_spectronaut_sample_table ---------------------------------------------
# Builds a sample name table directly from the embedded metadata in a
# long-format Spectronaut report (protein group report OR peptide report).
# Reads only the R.Condition and R.FileName columns, so it is fast even on
# large files.
#
# The resulting Data_name values ("raw_<FileName>") match the column names
# produced by read_peptide_spectronaut() after it pivots the data wide,
# allowing the sample table to be joined correctly without any manual editing.
# R.Condition is mapped to Group, preserving real group labels rather than the
# placeholder "A" that make_sample_name_table() assigns.
#
# Rows are ordered by Group then FileName so the handsontable displays in a
# logical order.
#
#' @title Build a sample table from a long-format Spectronaut report
#'
#' @description Reads only the \code{R.Condition} and \code{R.FileName}
#'   metadata columns from a long-format Spectronaut protein or peptide report
#'   and returns a sample name table ready for use in the Shiny apps or CLI
#'   workflows. The \code{Data_name} values are prefixed with \code{"raw_"} to
#'   match the column names produced when the same report is loaded with
#'   \code{read_peptide_spectronaut()}.
#'
#' @param file_path Character. Path to a long-format Spectronaut TSV report
#'   that contains \code{R.Condition} and \code{R.FileName} columns.
#'
#' @return A tibble with columns \code{Data_name}, \code{Sample_name},
#'   \code{Type}, \code{Group}, \code{Batch}, \code{Reference}. One row per
#'   unique sample. \code{Group} is populated from \code{R.Condition}.
#'   \code{Sample_name} is set to the raw file base name (without "raw_"
#'   prefix) for readability. \code{Type}, \code{Batch}, and \code{Reference}
#'   default to \code{"Lysate"}, \code{1}, and \code{"N"} respectively and can
#'   be edited in the Shiny handsontable before downstream analysis.
#'
#' @export
extract_spectronaut_sample_table <- function(file_path) {
  readr::read_tsv(
    file_path,
    col_select  = c("R.Condition", "R.FileName"),
    guess_max   = 10000,
    show_col_types = FALSE
  ) |>
    dplyr::distinct() |>
    dplyr::arrange(R.Condition, R.FileName) |>
    dplyr::mutate(
      Data_name   = R.FileName,
      Sample_name = R.FileName,
      Group       = R.Condition,
      Type        = "Lysate",
      Batch       = 1L,
      Reference   = "N"
    ) |>
    dplyr::select(Data_name, Sample_name, Type, Group, Batch, Reference)
}
