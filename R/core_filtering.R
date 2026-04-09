# core_filtering.R
# Pure R functions for filtering protein data by missing-value thresholds and
# by sample exclusion.
# No Shiny dependencies.


# Exported functions -----------------------------------------------------------

#' Filter proteins by missing-value thresholds
#'
#' @title Filter proteins by per-group and global observation thresholds
#'
#' @description Removes proteins (rows identified by \code{id}) that do not
#'   meet minimum-observation requirements. Two complementary criteria are
#'   applied in sequence:
#'
#'   \enumerate{
#'     \item \strong{Per-group threshold}: each protein must have at least
#'       \code{req_observations} valid intensity values within every group
#'       listed in \code{group_table}.
#'     \item \strong{Global threshold}: after the per-group filter, a protein
#'       must have at least \code{min_per_group} observations in at least
#'       \code{min_groups} groups.
#'   }
#'
#' @param protein_data A long-form tibble with at minimum the columns
#'   \code{id}, \code{Sample_name}, \code{Intensity}, and \code{Group}.
#' @param group_table A tibble with columns \code{Group} and
#'   \code{req_observations}. Use \code{make_observation_required_table()} to
#'   build this from a scalar or per-group lookup tibble.
#' @param min_per_group Numeric. Global minimum observations per group for a
#'   group to count toward \code{min_groups}. Defaults to \code{0}.
#' @param min_groups Numeric. Minimum number of groups that must each have at
#'   least \code{min_per_group} observations. Defaults to \code{0}.
#'
#' @return A long-form tibble with the same columns as \code{protein_data},
#'   containing only rows whose \code{id} values satisfy both threshold
#'   criteria.
#'
#' @export
filter_missing_values <- function(protein_data,
                                  group_table,
                                  min_per_group = 0,
                                  min_groups = 0) {

  group_tbl <- group_table |>
    dplyr::mutate(Group = factor(Group, levels = group_table$Group))

  obs_df <- protein_data |>
    dplyr::mutate(Group = factor(Group, levels = levels(group_tbl$Group))) |>
    dplyr::group_by(id) |>
    dplyr::count(Group, name = "Observations", .drop = FALSE) |>
    dplyr::ungroup()

  per_group_pass <- group_tbl |>
    dplyr::left_join(obs_df, by = "Group") |>
    dplyr::group_by(id) |>
    dplyr::filter(all(Observations >= req_observations)) |>
    dplyr::ungroup()

  global_pass <- per_group_pass |>
    dplyr::filter(Observations >= min_per_group) |>
    dplyr::count(id, name = "Number_of_groups") |>
    dplyr::filter(Number_of_groups >= min_groups)

  protein_data |>
    dplyr::filter(id %in% global_pass$id)
}


#' Remove samples from protein data and sample table
#'
#' @title Exclude specified samples from protein data and sample metadata
#'
#' @description Drops rows corresponding to one or more named samples from
#'   both the long-form protein data tibble and the sample metadata table.
#'
#' @param protein_data A long-form tibble with at minimum the columns
#'   \code{id}, \code{Sample_name}, and \code{Intensity}.
#' @param sample_table A tibble of sample metadata with a \code{Sample_name}
#'   column.
#' @param samples_to_remove Character vector of \code{Sample_name} values to
#'   exclude. Values not present in the data are silently ignored. Defaults to
#'   \code{character(0)} (no samples removed).
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{protein_data}}{Input \code{protein_data} with excluded
#'       sample rows removed.}
#'     \item{\code{sample_table}}{Input \code{sample_table} with excluded
#'       sample rows removed.}
#'   }
#'
#' @export
remove_samples <- function(protein_data,
                           sample_table,
                           samples_to_remove = character(0)) {

  list(
    protein_data = protein_data |>
      dplyr::filter(!Sample_name %in% samples_to_remove),
    sample_table = sample_table |>
      dplyr::filter(!Sample_name %in% samples_to_remove)
  )
}


# Internal helpers -------------------------------------------------------------

# make_observation_required_table ----------------------------------------------
# Builds a per-group minimum-observation table for use with
# filter_missing_values().
#
# group_table: tibble with at minimum a Group column.
# x1:          either a length-1 numeric (applied uniformly to all groups) or
#              a tibble with columns Group and req_observations.
make_observation_required_table <- function(group_table, x1) {

  if (tibble::is_tibble(x1)) {
    group_table |>
      dplyr::left_join(x1, by = "Group") |>
      dplyr::mutate(
        req_observations = dplyr::if_else(
          is.na(req_observations), 0, req_observations
        )
      )
  } else if (is.numeric(x1) && length(x1) == 1L) {
    group_table |>
      dplyr::mutate(req_observations = x1)
  } else {
    stop(
      "'x1' must be a length-1 numeric or a tibble with columns ",
      "'Group' and 'req_observations'."
    )
  }
}
