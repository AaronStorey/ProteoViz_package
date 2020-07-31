


# Import ------------------------------------------------------------------
#General
readSamplesReport <- function(filePath, type){

  if (type == "DIA"){
    #The file will contain a variable number of uncommented metadata.  This reads the first 100 lines,
    #Finds the column row and sets skip value to that - 1.
    first_read <- readLines(filePath)
    file_length <- length(first_read)
    skip_number <- grep("Visible", first_read) - 1

    read_csv(filePath, skip = skip_number, na = c("", "Missing Value"), n_max = file_length - skip_number - 2) %>%
      {purrr::set_names(., gsub(" ", "_", names(.)))} %>%
      filter(!is.na(Visible)) %>% #Removes END OF FILE line
      rename(id = `#`) %>%
      filter(!grepl("\\.", id)) %>% #Removes protein clusters, keeps protein group
      filter(!grepl("DECOY", Protein_Name))
  } else if (type == "TMT"){
    read_tsv(filePath, guess_max = 10000) %>%
      {purrr::set_names(., gsub(" ", "_", names(.)))}
  }

}
make_observation_required_table <- function(group_table, x1){
  if(tibble::is_tibble(x1)){
    df1 <-  group_table %>%
      left_join(x1) %>%
      mutate(req_observations = ifelse(is.na(req_observations), 0, req_observations))
  }

  if(is.numeric(x1) & length(x1) == 1) {
    df1 <- group_table %>%
      mutate(req_observations = x1)
  }

  df1

}



# Data --------------------------------------------------------------------

#General
clean_protein_df <- function(dfx, sample_table, log2check = TRUE, removecheck = TRUE, normToReferenceChannel = FALSE){

  column_name <- colnames(sample_table)[[1]]
  samples <- sample_table[[column_name]]
  #Fix this later
  colnames(sample_table)[[1]] <- "Data_name"

  df <- dfx %>%
    select(id, one_of(samples)) %>%
    gather(Data_name, Intensity, one_of(samples)) %>%
    left_join(sample_table, by = "Data_name")

  if(log2check){
    df <- df %>%
      mutate(Intensity = if_else(Intensity > 0, log2(Intensity), Intensity))
  }

  if(removecheck){
    df <- df %>%
      filter(Intensity > 0)
  }

  if(normToReferenceChannel){
    df <- normalizeReferenceChannels(df)
  }

  df
}
filter_protein_df <- function(dfx, group_table, GlobalGroupX = 0, GlobalGroupY = 0){

  group_tbl <- group_table %>%
    mutate(Group = factor(Group, levels = group_table$Group))

  df <- dfx %>%
    #Group needs to be a factor so that missing groups are not dropped
    mutate(Group = factor(Group, levels = group_tbl$Group)) %>%
    group_by(id) %>%
    count(Group, name = "Observations", .drop = FALSE) %>%
    ungroup()


  # Filters for proteins which meet all required number of observations
  df2 <- group_tbl %>%
    left_join(df, by = "Group") %>% #Check
    group_by(id) %>%
    filter(all(Observations >= req_observations)) %>%
    ungroup()

  #Then filters for proteins that meet the global requirements
  obs_min <- GlobalGroupX
  group_min <- GlobalGroupY

  df3 <- df2 %>%
    filter(Observations >= obs_min) %>%
    count(id, name = "Number_of_groups") %>%
    filter(Number_of_groups >= group_min)

  dfx %>%
    filter(id %in% df3$id)

}
protein_df_to_matrix <- function(dfx){
  dfx %>%
    select(id, Sample_name, Intensity) %>%
    spread(Sample_name, Intensity) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("id") %>%
    as.matrix()
}
normalize_protein_df <- function(dfx, normOption = "None"){

  df <- dfx

  if(normOption == "None"){
    df1 <- df
  } else if(normOption == "Mean"){
    df1 <- meanNorm(df)
  } else if(normOption == "Median"){
    df1 <- medianNorm(df)
  } else if(normOption == "Quantile"){
    df1 <- quantNorm(df)
  } else if(normOption == "Cyclic Loess"){
    df1 <- cycLoessNorm(df)
  } else if(normOption == "Rlr"){
    df1 <- rlrNorm(df)
  } else if(normOption == "Gi"){
    df1 <- giNorm(df)
  } else if(normOption == "Cyclic Loess then Median"){
    df1 <- medianNorm(cycLoessNorm(df))
  }

  df1


}
matrix_to_long_form <- function(dfx) {

  dfx %>%
    tibble::as_tibble(rownames = "id") %>%
    gather(Sample_name, Intensity, -id) %>%
    filter(!is.na(Intensity))

}
normalizeReferenceChannels <- function(protein_df){
  #Must be log2 transformed
  protein_df %>%
    group_by(id) %>%
    mutate(Ref_Intensity = mean(Intensity[which(Reference == "Y")])) %>%
    group_by(id, Batch) %>%
    mutate(Batch_Ref_Intensity = mean(Intensity[which(Reference == "Y")])) %>%
    ungroup() %>%
    mutate(Intensity = Intensity - (Batch_Ref_Intensity - Ref_Intensity)) %>%
    select(-Ref_Intensity, -Batch_Ref_Intensity)
}



# Phospho -----------------------------------------------------------------

readPhosphoSites <- function(filePath){
  read_tsv(filePath, guess_max = 20000) %>%
    {purrr::set_names(., gsub(" ", "_", names(.)))}
}
preclean_phospho_df <- function(dfx){
  df <- dfx %>%
    filter(is.na(Reverse),
                  is.na(Potential_contaminant)) %>%
    filter(Localization_prob > 0.75) %>%
    select(id,Protein_group_IDs, matches("\\_{3}[[:digit:]]$")) %>%
    gather(Data_name, Intensity, matches("\\_{3}[[:digit:]]$"))

  #Summarizes all phos_number columns into single Intensity for each phosphosite
  df2 <- df %>%
    filter(Intensity > 0) %>%
    separate(Data_name, into = c("Data_name", "Phos_number"), sep = "___") %>%
    group_by(id, Protein_group_IDs, Data_name) %>%
    summarize(Intensity = sum(Intensity), .groups = "drop") %>%
    ungroup()
}
clean_phospho_df <-function(dfx,
                            sampleTable,
                            log2check = TRUE,
                            removecheck = TRUE,
                            normToReferenceChannel = FALSE,
                            normToProtein = TRUE,
                            proteinData = NULL){
    df2 <- dfx

    df3 <- df2 %>%
      right_join(sampleTable, by = "Data_name")


    if (log2check) {
      df3 <- df3 %>%
        mutate(Intensity = if_else(Intensity > 0, log2(Intensity), Intensity))
    }

    if (removecheck) {
      df3 <- df3 %>%
        filter(Intensity > 0)
    }

    if (normToReferenceChannel) {
      df3 <- normalizeReferenceChannels(df3)
    }

    if (normToProtein) {
      df3 <- normalizePhosphoToProtein(df3, proteinData)
    }

    df3
  }
normalizePhosphoToProtein <- function(phospho_df, protein_df = NULL){
  # Must be log2 transformed
  # protein_df and phospho_df must be in tidy(long) format
  # Metadata table must contain a "Sample_name" column which matches protein samples to phospho samples
  if (is.null(protein_df)){
    stop("No protein data supplied")
  }
  protein_estimate_df <- protein_df %>%
    select(id, Sample_name, Batch, Intensity) %>%
    group_by(id, Batch) %>%
    mutate(Protein_intensity = Intensity - mean(Intensity, na.rm = TRUE)) %>%
    ungroup() %>%
    select(Protein_group_IDs = id, Sample_name, Batch, Protein_intensity) %>%
    mutate(Protein_group_IDs = as.character(Protein_group_IDs))

  phospho_df %>%
    left_join(protein_estimate_df,
                     by = c("Protein_group_IDs", "Batch", "Sample_name")) %>%
    mutate(
      Protein_group_match = ifelse(is.na(Protein_intensity), "0", "1"),
      Protein_group_match = ifelse(grepl("\\;", Protein_group_IDs), "Many", Protein_group_match)
    ) %>%
    mutate(Intensity = ifelse(
      Protein_group_match == 1,
      Intensity - Protein_intensity,
      Intensity
    )) %>%
    select(-Protein_intensity)

}
makeProteinMatchdf <- function(dfx, normToProtein = TRUE){
  if(normToProtein){
    dfx %>%
      distinct(id, Protein_group_IDs, Protein_group_match) %>%
      group_by(id, Protein_group_IDs) %>%
      summarize(Protein_group_match = paste0(Protein_group_match, collapse = "_")) %>%
      ungroup() %>%
      mutate(Protein_group_match = ifelse(Protein_group_match %in% c("0", "1", "Many"), Protein_group_match, "Some"))
  } else return(NULL)
}


# "Separates phosphosites" that don't match protein groups.
# Attaches Protein_group_match value to Sample name
split_protein_match <- function(dfx){
  dfx %>%
    mutate(Protein_group_match = sub("0", "Nomatch", Protein_group_match) %>%
             gsub("Many", "Multimatch", .) %>%
             sub("1", "Match", .)) %>%
    unite(Sample_name, Sample_name, Protein_group_match)
}

reunite_protein_match <- function(dfx){
  dfx %>%
    mutate(Sample_name = gsub("\\_[^\\_]+$", "", Sample_name))
}

reorder_split_phospho <- function(sample_names){
  paste0(rep(sample_names, each = 3), c("_Nomatch", "_Match", "_Multimatch"))
}


# Save --------------------------------------------------------------------

makeProteinMeta <- function(x, type){

  if (type == "DIA"){
    x %>%
      select(id, Protein_Name, Accession_Number, Identified_Peptide_Count, Protein_Group_Score) %>%
      mutate(
        Description = stringr::str_extract(Protein_Name, "(?<= )[^\\|]+(?= OS\\=)"),
        Gene_name = stringr::str_extract(Protein_Name, "(?<=GN\\=)[^\\|]+(?= PE\\=)"),
        Uniprot_ID = stringr::str_extract(Protein_Name, "(?<=\\|)[^\\|]+(?=\\|)")
      )
  } else if (type == "TMT"){
    #Selects Phosphosite IDs only if present in the MaxQuant data
    column_test <- c("Majority_protein_IDs", "Fasta_headers", "Score", "id", "Phospho_(STY)_site_IDs")
    column_test <- column_test[column_test %in% names(x)]

    x %>%
      select(one_of(column_test)) %>%
      mutate(Description = stringr::str_extract(Fasta_headers, "(?<= )[^\\|]+(?= OS\\=)"),
                    Gene_name = stringr::str_extract(Fasta_headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)"),
                    Uniprot_ID = stringr::str_extract(Fasta_headers, "(?<=\\|)[^\\|]+(?=\\|)"))
  }

}
makePhosphoMeta <- function(x){

  x %>%
    select(Proteins:Score, Amino_acid, Sequence_window, `Phospho_(STY)_Probabilities`, Charge, id:Evidence_IDs) %>%
    mutate(Description = stringr::str_extract(Fasta_headers, "(?<= )[^\\|]+(?= OS\\=)"),
                  Gene_name = stringr::str_extract(Fasta_headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)"),
                  Uniprot_ID = stringr::str_extract(Fasta_headers, "(?<=\\|)[^\\|]+(?=\\|)"),
                  Flanking = gsub("\\;.*$", "", Sequence_window) %>%
                    stringr::str_sub(9,23) %>%
                    paste0("-p"))
}
save_protein_meta <- function(dfx, projectName, type) {
  if(!dir.exists(projectName)) dir.create(projectName, recursive = TRUE)
  file_prefix <- gsub("^.*\\/", "", projectName)

  makeProteinMeta(dfx, type)  %>%
    write_tsv(paste0(projectName, "/", file_prefix, "_", "protein_metadata.tsv"))
}
save_phospho_meta <- function(dfx, projectName, normToProtein = TRUE, proteinmatchdf = NULL){

  df <- makePhosphoMeta(dfx)
  if(!dir.exists(projectName)) dir.create(projectName, recursive = TRUE)
  file_prefix <- gsub("^.*\\/", "", projectName)

  if(normToProtein){
    df <- left_join(df, proteinmatchdf, by = c("id", "Protein_group_IDs"))
  }

  df %>%
    write_tsv(paste0(projectName, "/", file_prefix, "_", "phospho_metadata.tsv"))
}
save_quantitative_data <- function(dfx, projectName, sample_table, dataType = "protein"){
  if(!dir.exists(projectName)) dir.create(projectName, recursive = TRUE)
  file_prefix <- gsub("^.*\\/", "", projectName)
  dfx %>%
    select(id, Sample_name, Intensity) %>%
    mutate(Sample_name = factor(Sample_name, levels = unique(sample_table$Sample_name))) %>%
    spread(Sample_name, Intensity) %>%
    mutate(id = as.integer(id)) %>%
    arrange(id) %>%
    write_tsv(paste0(projectName, "/", file_prefix, "_", dataType, "_", "quantitative_data.tsv"))
}
save_design_table <- function(sampleTable, projectName){
  if(!dir.exists(projectName)) dir.create(projectName, recursive = TRUE)
  file_prefix <- gsub("^.*\\/", "", projectName)

  sampleTable %>%
    distinct(Sample_name, Group, Batch) %>%
    mutate(Block = "1") %>%
    write_tsv(paste0(projectName, "/", file_prefix, "_", "design_table.tsv"))

}
make_contrast_list <- function(sampleTable, projectName){
  if(!dir.exists(projectName)) dir.create(projectName, recursive = TRUE)
  file_prefix <- gsub("^.*\\/", "", projectName)


  #Removes Pool from list of groups.  For convenience
  Groups <- unique(sampleTable$Group)
  Groups <- Groups[Groups != "Pool"]

  expand.grid(x1 = Groups, y1 = Groups) %>%
    mutate(Contrast_name = paste0(y1, " - ", x1)) %>%
    filter(x1 != y1) %>%
    select(Contrast_name) %>%
    write_tsv(paste0(projectName, "/", file_prefix, "_", "Sample_contrast_list.tsv"))

}




# Shiny -------------------------------------------------------------------


#Only used in Shiny
makeSampleNameTable <- function(protein_df, type){


  if(type == "DIA"){
    x1 <- grep("\\.", names(protein_df))
  } else if (type == "TMT"){
    x1 <- grep("Reporter\\_intensity\\_corrected\\_.{3,}", names(protein_df))
  }

  tibble::tibble(Data_name = names(protein_df)[x1],
                 Sample_name = paste0("Sample_", seq_along(x1)),
                 Type = "Lysate",
                 Group = "A",
                 Batch = 1,
                 Reference = "N"
  )
}

readSampleNameTable <- function(sampleFile){
  x1 <- read_tsv(sampleFile)
  colnames(x1)[1] <- "Data_name"
  x1
}

make_boxplot <- function(dfx, sampleTable){

  #Fixes order for multisite phospho.
  if (all(grepl("Nomatch|Multimatch|Match", dfx$Sample_name))) {
    sample_levels <- reorder_split_phospho(unique(sampleTable$Sample_name))
  } else {sample_levels <- unique(sampleTable$Sample_name) }

  dfx %>%
    mutate(
      Sample_name = factor(Sample_name,
                           levels = sample_levels)
    ) %>%
    ggplot(aes(x = Sample_name, y = Intensity, fill = Group)) + geom_boxplot() +
    scale_fill_viridis_d(end = 0.9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9)) +
    labs(title = "Before normalization")

}

make_boxplot2 <- function(dfx, sampleTable){

  #Fixes order for multisite phospho. What a fucking mess.
  if (all(grepl("Nomatch|Multimatch|Match", dfx$Sample_name))) {
    sample_levels <- reorder_split_phospho(unique(sampleTable$Sample_name))
    sample_groups <- rep(sampleTable$Sample_name, each = 3)
    sampleTable <- tibble(Sample_name = sample_levels,
                          Group = rep(sampleTable$Group, each = 3))
  } else {sample_levels <- unique(sampleTable$Sample_name) }

  df <- dfx %>%
    left_join(sampleTable, by = "Sample_name")

  df %>%
    mutate(
      Sample_name = factor(Sample_name,
                           levels = sample_levels)
    ) %>%
    ggplot(aes(x = Sample_name, y = Intensity, fill = Group)) + geom_boxplot() +
    scale_fill_viridis_d(end = 0.9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))
}

make_CV_histogram <- function(x1, samples = NULL, unlog = FALSE) {
  if(unlog) {
    x1 <- x1 %>%
      mutate(Intensity = 2^Intensity)
  }
  if(!is.null(samples)){
    x1 <- x1[x1$Sample_name %in% samples,]
  }
  
  x1 %>%
    group_by(id) %>%
    summarize(CV = sd(Intensity, na.rm = TRUE) / mean(Intensity, na.rm = TRUE)) %>%
    ggplot(aes(x = CV)) + geom_histogram(bins = 100)
}


PCA_plotly <- function(dfx, sampleTable){
  a1 <- na.omit(dfx)
  sample_df <- sampleTable

  PCA_df <- prcomp(t(a1))
  #Percent variance explained
  PC1_var <- summary(PCA_df)$importance[2,1] %>% {signif(. * 100, 3)} # Pretty format
  PC2_var <- summary(PCA_df)$importance[2,2] %>% {signif(. * 100, 3)} #

  df4 <- PCA_df$x %>%
    as_tibble(rownames = "Sample_name") %>%
    right_join(sample_df, ., by = "Sample_name")

  p1 <- df4 %>%
    ggplot(aes(x = PC1, y = PC2, color = Group,
               text = paste(
                 "Sample name:", Sample_name, "\n",
                 "Data name:", Data_name, "\n",
                 "Batch:", Batch, "\n"
               ))) +
    geom_point(size = 5) +
    scale_color_viridis_d(end = 0.9) +
    labs(x = paste0("PC1 (", PC1_var, "%)"),
         y = paste0("PC2 (", PC2_var, "%)")) +
    theme(legend.position = "bottom")

  ggplotly(p1, source = "prot_PCA", tooltip = "text") %>%
    layout(dragmode = "select",
           hoverlabel=list(bgcolor="black"))
}

make_cor_plot <- function(dfx, sampleTable, scaleRows = FALSE){


  a1 <- dfx %>%
    mutate(
      Sample_name = factor(Sample_name,
                           levels = unique(sampleTable$Sample_name))
    ) %>%
    select(id, Sample_name, Intensity) %>%
    spread(Sample_name, Intensity) %>%
    as.data.frame() %>%
    column_to_rownames("id") %>%
    as.matrix()

  sample_labels <- colnames(a1)

  if(scaleRows){
    a1 <- a1 - rowMeans(a1, na.rm = TRUE)
  }


  cor_mat <- cor(a1, use = "pairwise.complete.obs")

  hm_corr = Heatmap(cor_mat,
                    col = circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7)),RColorBrewer::brewer.pal(8, "Reds")),
                    heatmap_legend_param = list(color_bar = "continuous",
                                                legend_direction = "horizontal",
                                                legend_width = unit(5, "cm"),
                                                title_position = "topcenter"),
                    name = "Pearson correlation",
                    column_names_gp = gpar(fontsize = 12),
                    row_names_gp = gpar(fontsize = 12),
                    # top_annotation = ColAnn,
                    # left_annotation = RowAnn,
                    column_labels = sample_labels, row_labels = sample_labels
  )

  draw(hm_corr, heatmap_legend_side = "top")
}

make_NA_plot <- function(dfx, sampleTable){


  a1 <- dfx %>%
    select(id, Sample_name, Intensity) %>%
    mutate(
      Sample_name = factor(Sample_name,
                           levels = unique(sampleTable$Sample_name))
    ) %>%
    spread(Sample_name, Intensity) %>%
    as.data.frame() %>%
    column_to_rownames("id") %>%
    as.matrix()

  sample_labels <- colnames(a1)


  missing = !is.na(a1)
  complete = apply(missing, 1, all)
  completeNA = apply(!missing, 1, all)
  missing = missing[!complete & !completeNA,]

  hm_clust = Heatmap(
    missing + 0,
    col = c("white", "black"),
    column_names_side = "top",
    show_row_names = FALSE,
    show_column_names = TRUE,
    name = "Missing values pattern",
    column_names_gp = gpar(fontsize = 16),
    heatmap_legend_param = list(
      at = c(0, 1),
      labels = c("Missing value", "Valid value")
    ),
    column_labels = sample_labels
  )

  draw(hm_clust)
}
