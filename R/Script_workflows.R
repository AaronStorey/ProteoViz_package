#' Processes a MaxQuant PhosphoTMT data set and outputs Limma-ready tables
#'
#'
#' @export
#'
#'
ScriptApp1PhosphoTMT <- function(project_name, protein_file, samples_file, phospho_file, P1exclude, P1groupexclude,
                                 GlobalGroupX, GlobalGroupY, required_observations, normalization_option,
                                 log2_check, remove_0_values, normalize_to_reference, normalize_to_protein,
                                 type) {

  #Sample table
  sampleNameTable2 <- readSampleNameTable(samples_file)
  sampleNameTable3 <- sampleNameTable2 %>%
    dplyr::filter(!Sample_name %in% P1exclude,!Group %in% P1groupexclude) %>%
    dplyr::filter(Type == "Lysate")

  #Group table
  groupTable <- sampleNameTable3 %>%
    dplyr::count(Group, name = "Number_of_samples") %>%
    dplyr::ungroup()


  #Table that lists required number of observations for each group,
  #as specified from the input sliders
  groupTable2 <- make_observation_required_table(groupTable, required_observations)


  # Protein dfs -------------------------------------------------------------

  protein_df <- readSamplesReport(protein_file, type)
  protein_df2 <- clean_protein_df(protein_df, sampleNameTable3, log2_check, remove_0_values, normalize_to_reference)
  protein_df3 <- filter_protein_df(protein_df2, groupTable2, GlobalGroupX, GlobalGroupY)
  protein_df4 <- protein_df_to_matrix(protein_df3)
  protein_df5 <- normalize_protein_df(protein_df4, normalization_option)
  protein_df6 <- matrix_to_long_form(protein_df5)


  # Phospho section ---------------------------------------------------------

  phosphoSampleNameTable3 <- sampleNameTable2 %>%
    dplyr::filter(!Sample_name %in% P1exclude,!Group %in% P1groupexclude) %>%
    dplyr::filter(Type == "Phospho")
  phosphogroupTable <- phosphoSampleNameTable3 %>%
    dplyr::count(Group, name = "Number_of_samples") %>%
    dplyr::ungroup()
  phosphogroupTable2 <- make_observation_required_table(phosphogroupTable, required_observations)



  phospho_df <- readPhosphoSites(phospho_file)
  phospho_to_protein <- phospho_df %>%
    dplyr::select(id, Protein_group_IDs)

  phospho_df1 <- preclean_phospho_df(phospho_df)
  phospho_df2 <- clean_phospho_df(phospho_df1, phosphoSampleNameTable3, log2_check, remove_0_values, normalize_to_reference,
                                  normalize_to_protein, protein_df3)
  protein_match_df <- makeProteinMatchdf(phospho_df2, normalize_to_protein)

  #These use the protein df functions
  phospho_df3 <- filter_protein_df(phospho_df2, phosphogroupTable2, GlobalGroupX, GlobalGroupY)
  phospho_df4 <- protein_df_to_matrix(phospho_df3)
  phospho_df5 <- normalize_protein_df(phospho_df4, normalization_option)
  phospho_df6 <- matrix_to_long_form(phospho_df5)


  # Save_dfs ----------------------------------------------------------------


  save_quantitative_data(protein_df6, project_name, sampleNameTable3, "protein")
  save_protein_meta(protein_df, project_name, type)
  save_quantitative_data(phospho_df6, project_name, phosphoSampleNameTable3, "phospho")
  save_phospho_meta(phospho_df, project_name, normalize_to_protein, protein_match_df)
  save_design_table(sampleNameTable3, project_name)
  make_contrast_list(sampleNameTable3, project_name)

}


#' Processes a ScaffoldDIA Samples Report and outputs Limma-ready tables
#'
#'
#' @export
#'
#'
ScriptApp1DIA <- function(project_name, protein_file, samples_file, P1exclude, P1groupexclude,
                                 GlobalGroupX, GlobalGroupY, required_observations, normalization_option,
                                 log2_check, remove_0_values, normalize_to_reference, type) {


  #Sample table
  sampleNameTable2 <- readSampleNameTable(samples_file)
  sampleNameTable3 <- sampleNameTable2 %>%
    dplyr::filter(!Sample_name %in% P1exclude,!Group %in% P1groupexclude)

  #Group table
  groupTable <- sampleNameTable3 %>%
    dplyr::count(Group, name = "Number_of_samples") %>%
    dplyr::ungroup()


  #Table that lists required number of observations for each group,
  #as specified from the input sliders
  groupTable2 <- make_observation_required_table(groupTable, required_observations)


  #Protein dfs
  protein_df <- readSamplesReport(protein_file, type)
  protein_df2 <- clean_protein_df(protein_df, sampleNameTable3, log2_check, remove_0_values, normalize_to_reference)
  protein_df3 <- filter_protein_df(protein_df2, groupTable2, GlobalGroupX, GlobalGroupY)
  protein_df4 <- protein_df_to_matrix(protein_df3)
  protein_df5 <- normalize_protein_df(protein_df4, normalization_option)
  protein_df6 <- matrix_to_long_form(protein_df5)



  # Save_dfs ----------------------------------------------------------------



  save_quantitative_data(protein_df6, project_name, sampleNameTable3, "protein")
  save_protein_meta(protein_df, project_name, type)
  save_design_table(sampleNameTable3, project_name)
  make_contrast_list(sampleNameTable3, project_name)

}
