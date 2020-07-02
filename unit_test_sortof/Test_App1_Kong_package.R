# Parameters --------------------------------------------------------------


project_name <- "test_output/Scripting_test_Kong"
protein_file <- "data/DIA_Kong_042920/Samples Report of Kong_042920.csv"
samples_file <- "data/DIA_Kong_042920/Kong_042920_Sample_table_fixed.txt"


#Samples to exclude
P1exclude <- ""
P1groupexclude <- ""

#Require at least X observations in Y groups
GlobalGroupX <-  0
GlobalGroupY <-  0

#Required observations. Either a single number or a tibble specifying Group and number
required_observations <- tibble::tibble(
  Group = "Pool",
  req_observations = 13
)

normalization_option <- "Cyclic Loess"
log2_check <- TRUE
remove_0_values <- TRUE
normalize_to_reference <- FALSE
type <- "DIA"

# Script ------------------------------------------------------------------

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

