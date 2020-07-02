# Parameters --------------------------------------------------------------

project_name <- "test_output/Chambers_110518_Splitphospho"
protein_file <- "data/TMTphos_Chambers_110518/proteinGroups.txt"
samples_file <- "data/TMTphos_Chambers_110518/Chambers_110518_Sample_name_table.txt"
phospho_file <- "data/TMTphos_Chambers_110518/Phospho (STY)Sites.txt"


#Samples to exclude
P1exclude <- ""
P1groupexclude <- "Blank"

#Require at least X observations in Y groups
GlobalGroupX <-  0
GlobalGroupY <-  0

#Required observations for specific groups. Either a single number for all groups or a tibble specifying Group and number
required_observations <- 0
normalization_option <- "Cyclic Loess"

log2_check <- TRUE
remove_0_values <- TRUE
normalize_to_reference <- TRUE
normalize_to_protein <- TRUE
type <- "TMT"


# Script ------------------------------------------------------------------

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
phospho_df3 <- filter_protein_df(phospho_df2, phosphogroupTable2, GlobalGroupX, GlobalGroupY) %>%
  split_protein_match()
phospho_df4 <- protein_df_to_matrix(phospho_df3)
phospho_df5 <- normalize_protein_df(phospho_df4, normalization_option)
phospho_df6 <- matrix_to_long_form(phospho_df5) %>%
  reunite_protein_match()


# Save_dfs ----------------------------------------------------------------


save_quantitative_data(protein_df6, project_name, sampleNameTable3, "protein")
save_protein_meta(protein_df, project_name, type)
save_quantitative_data(phospho_df6, project_name, phosphoSampleNameTable3, "phospho")
save_phospho_meta(phospho_df, project_name, normalize_to_protein, protein_match_df)
save_design_table(sampleNameTable3, project_name)
make_contrast_list(sampleNameTable3, project_name)

