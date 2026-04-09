# core_normalization.R
# Exported normalization entry point and unexported per-method helpers.
# All functions are stateless; no Shiny dependencies.


# Exported ----------------------------------------------------------------

#' Normalize protein intensities
#'
#' @title Normalize protein intensities
#'
#' @description Applies a sample-level normalization method to a long-form
#'   protein intensity tibble and returns a tibble in the same format.
#'   Internally the data are pivoted to a proteins-by-samples matrix,
#'   normalized, and pivoted back to long form; rows that were \code{NA}
#'   (i.e. missing values) before normalization remain absent from the output.
#'
#' @param protein_data A long-form \code{\link[tibble]{tibble}} with at least
#'   three columns:
#'   \describe{
#'     \item{id}{Character or numeric protein/row identifier.}
#'     \item{Sample_name}{Character sample label.}
#'     \item{Intensity}{Numeric intensity value (log2-scale expected for all
#'       methods except \code{"None"}).}
#'   }
#'   Additional columns (e.g. \code{Group}, \code{Batch}) are preserved
#'   through a join after normalization.
#' @param method Character scalar. One of:
#'   \itemize{
#'     \item \code{"None"} — no normalization (default)
#'     \item \code{"Mean"} — mean normalization
#'     \item \code{"Median"} — median normalization
#'     \item \code{"Quantile"} — quantile normalization via
#'       \pkg{preprocessCore}
#'     \item \code{"Cyclic Loess"} — cyclic LOESS normalization via
#'       \pkg{limma}
#'     \item \code{"Rlr"} — robust linear-regression normalization via
#'       \pkg{NormalyzerDE}
#'     \item \code{"Gi"} — global-intensity normalization via
#'       \pkg{NormalyzerDE}
#'     \item \code{"Cyclic Loess then Median"} — cyclic LOESS followed by
#'       median normalization
#'   }
#'
#' @return A long-form \code{\link[tibble]{tibble}} with the same columns as
#'   \code{protein_data}.  The \code{Intensity} column reflects the
#'   normalized values.  Row order may differ from the input.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # protein_data is a long-form tibble with columns id, Sample_name, Intensity
#' norm_data <- normalize_intensities(protein_data, method = "Median")
#' }
normalize_intensities <- function(protein_data, method = "None") {

  valid_methods <- c(
    "None", "Mean", "Median", "Quantile",
    "Cyclic Loess", "Rlr", "Gi", "Cyclic Loess then Median"
  )
  if (!method %in% valid_methods) {
    stop(
      "`method` must be one of: ",
      paste(paste0('"', valid_methods, '"'), collapse = ", "),
      call. = FALSE
    )
  }

  required_cols <- c("id", "Sample_name", "Intensity")
  missing_cols <- setdiff(required_cols, colnames(protein_data))
  if (length(missing_cols) > 0) {
    stop(
      "`protein_data` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------
  # Convert long form -> matrix (proteins x samples)
  # ------------------------------------------------------------------
  mat <- protein_data |>
    dplyr::select(id, Sample_name, Intensity) |>
    tidyr::pivot_wider(names_from = Sample_name, values_from = Intensity) |>
    as.data.frame() |>
    tibble::column_to_rownames("id") |>
    as.matrix()

  # ------------------------------------------------------------------
  # Apply normalization
  # ------------------------------------------------------------------
  norm_mat <- switch(
    method,
    "None"                    = mat,
    "Mean"                    = mean_norm(mat),
    "Median"                  = median_norm(mat),
    "Quantile"                = quant_norm(mat),
    "Cyclic Loess"            = cyc_loess_norm(mat),
    "Rlr"                     = rlr_norm(mat),
    "Gi"                      = gi_norm(mat),
    "Cyclic Loess then Median" = median_norm(cyc_loess_norm(mat))
  )

  # ------------------------------------------------------------------
  # Convert matrix -> long form, drop NAs (missing values)
  # ------------------------------------------------------------------
  norm_long <- norm_mat |>
    tibble::as_tibble(rownames = "id") |>
    tidyr::pivot_longer(
      cols      = -id,
      names_to  = "Sample_name",
      values_to = "Intensity"
    ) |>
    dplyr::filter(!is.na(Intensity))

  # ------------------------------------------------------------------
  # Re-attach any extra columns from protein_data (e.g. Group, Batch)
  # ------------------------------------------------------------------
  extra_cols <- setdiff(colnames(protein_data), c("id", "Sample_name", "Intensity"))

  if (length(extra_cols) > 0) {
    meta <- protein_data |>
      dplyr::select(dplyr::all_of(c("id", "Sample_name", extra_cols))) |>
      dplyr::distinct()

    # coerce id to character for safe join
    meta     <- meta     |> dplyr::mutate(id = as.character(id))
    norm_long <- norm_long |> dplyr::mutate(id = as.character(id))

    norm_long <- norm_long |>
      dplyr::left_join(meta, by = c("id", "Sample_name"))
  }

  norm_long
}


# Internal helpers --------------------------------------------------------
# These are unexported; called only by normalize_intensities().
# Each function accepts a numeric matrix (proteins x samples) and returns
# a numeric matrix of the same dimensions.

log_norm <- function(dat) {
  log_int <- log2(dat)
  log_int[is.infinite(as.matrix(log_int))] <- NA
  as.matrix(log_int)
}

median_norm <- function(log_dat) {
  sample_med <- apply(log_dat, 2, median, na.rm = TRUE)
  mean_med   <- mean(sample_med, na.rm = TRUE)
  out        <- t(t(log_dat) / sample_med) * mean_med
  as.matrix(out)
}

mean_norm <- function(log_dat) {
  sample_mean <- apply(log_dat, 2, mean, na.rm = TRUE)
  mean_mean   <- mean(sample_mean, na.rm = TRUE)
  out         <- t(t(log_dat) / sample_mean) * mean_mean
  as.matrix(out)
}

vsn_norm <- function(dat) {
  vsn_normed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
  colnames(vsn_normed) <- colnames(dat)
  rownames(vsn_normed) <- rownames(dat)
  as.matrix(vsn_normed)
}

quant_norm <- function(log_dat) {
  quant_normed <- preprocessCore::normalize.quantiles(as.matrix(log_dat), copy = TRUE)
  colnames(quant_normed) <- colnames(log_dat)
  rownames(quant_normed) <- rownames(log_dat)
  as.matrix(quant_normed)
}

cyc_loess_norm <- function(log_dat) {
  cyc_normed <- limma::normalizeCyclicLoess(as.matrix(log_dat), method = "fast")
  colnames(cyc_normed) <- colnames(log_dat)
  rownames(cyc_normed) <- rownames(log_dat)
  as.matrix(cyc_normed)
}

rlr_norm <- function(log_dat) {
  rlr_normed <- NormalyzerDE::performGlobalRLRNormalization(as.matrix(log_dat), noLogTransform = TRUE)
  colnames(rlr_normed) <- colnames(log_dat)
  rownames(rlr_normed) <- rownames(log_dat)
  as.matrix(rlr_normed)
}

gi_norm <- function(log_dat) {
  gi_normed <- NormalyzerDE::globalIntensityNormalization(as.matrix(log_dat), noLogTransform = TRUE)
  colnames(gi_normed) <- colnames(log_dat)
  rownames(gi_normed) <- rownames(log_dat)
  as.matrix(gi_normed)
}
