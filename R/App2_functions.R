


###########Import functions

#General
makeDesignTable <- function(sampleTable){
  df1 <- sampleTable %>%
    mutate(Group = factor(Group, levels = unique(Group)))

  Groups <- df1$Group

  if(length(unique(Groups)) < 2) return(NULL)

  df2 <- model.matrix(~0 + Groups)

  rownames(df2) <- df1$Sample_name
  colnames(df2) <- levels(Groups)

  df2
}
makeContrastTable <- function(design_table, contrast_data){
  df <- contrast_data %>%
    filter(Contrast_name != "")

  if(nrow(df) < 1) return(NULL)

  makeContrasts(contrasts = as.list(df$Contrast_name), levels = colnames(design_table))
}
runLimma <- function(limmaInput, designTable, sampleTable, block = FALSE){
  #Rearrange quant data to match order of sampleNameTable
  quant_data <- limmaInput[,sampleTable$Sample_name]

  if(block){

    sampleTable <- sampleTable %>%
      mutate(Block = factor(Block, levels = unique(Block)))

    corfit <- duplicateCorrelation(quant_data, designTable, block=sampleTable$Block)

    fit1 <- lmFit(quant_data, designTable, block = sampleTable$Block, correlation = corfit$consensus.correlation)
  } else {
    fit1 <- lmFit(quant_data, designTable)
  }

  fit1
}
fitContrasts <- function(limmaOutput, contrastTable){
  fit2 <- contrasts.fit(limmaOutput, contrastTable)
  eBayes(fit2)
}
extractLimmaData <- function(limmaData, comparison){

  a1 <- topTable(limmaData, comparison, number = Inf) %>%
    rownames_to_column("id") %>%
    mutate(Comparison = comparison) %>%
    mutate(id = as.integer(id)) %>%
    as_tibble()

  a1
}
combineLimmaData <- function(contrastFit, designTable, contrastTable){

  req(contrastFit)
  req(designTable)
  req(contrastTable)

  Comparisons <- contrastTable$Contrast_name

  map_df(Comparisons, ~extractLimmaData(contrastFit, .x)) %>%
    filter(!is.na(adj.P.Val))
}
makeTidyTestResults <- function(testResults){
  #Probably an easier way to get to a named list of up/downregulated hits
  #Filters for up/downregulated hits, makes nested tibble
    df1 <- testResults@.Data %>%
    as_tibble(rownames = "id") %>%
    gather(Test, Result, -id) %>%
    filter(Result %in% c(1, -1))

    df1 %>%
      mutate(Result = as.integer(Result),
             Result = ifelse(Result == 1, "Up", "Down")) %>%
      group_by(Test, Result) %>%
      nest() %>%
      unite(Test, Test, Result, sep = " ")
}
makeTidyTestResultsList <- function(tidyTestResults){
  a2 <- map(tidyTestResults$data, "id")
  names(a2) <- tidyTestResults$Test
  a2
}
makeSummary <- function(quantData, combinedData, metaData){
  spread_protein_limma <- combinedData %>%
    dplyr::select(id, logFC, P.Value, adj.P.Val, Comparison) %>%
    gather(Type, Value, logFC, P.Value, adj.P.Val) %>%
    unite(Type, Comparison, Type, sep = " ") %>%
    spread(Type, Value) %>%
    mutate(id = as.integer(id))

  meta_df <- metaData %>%
    mutate(id = as.integer(id))

  meta_df %>%
    right_join(quantData) %>%
    left_join(spread_protein_limma) %>%
    write_tsv("Summarized_output.tsv")
}

makeSummary2 <- function(quantData, combinedData, metaData, prefix){
  spread_protein_limma <- combinedData %>%
    dplyr::select(id, logFC, P.Value, adj.P.Val, Comparison) %>%
    gather(Type, Value, logFC, P.Value, adj.P.Val) %>%
    unite(Type, Comparison, Type, sep = " ") %>%
    spread(Type, Value) %>%
    mutate(id = as.integer(id))

  meta_df <- metaData %>%
    mutate(id = as.integer(id))

  meta_df %>%
    right_join(quantData) %>%
    left_join(spread_protein_limma) %>%
    write_tsv(paste0(prefix, "_Summarized_output.tsv"))
}

makeSummary3 <- function(quantData, combinedData, metaData, prefix){
  spread_protein_limma <- combinedData %>%
    dplyr::select(id, logFC, P.Value, adj.P.Val, Comparison) %>%
    gather(Type, Value, logFC, P.Value, adj.P.Val) %>%
    unite(Type, Comparison, Type, sep = " ") %>%
    spread(Type, Value) %>%
    mutate(id = as.integer(id))

  meta_df <- metaData %>%
    mutate(id = as.integer(id))

  meta_df %>%
    right_join(quantData) %>%
    right_join(spread_protein_limma) %>%
    write_tsv(paste0(prefix, "_Summarized_output.tsv"))
}

makeProteinHeatHeatmap <- function(selectData, metaData, sampleOrder, groupOrder, groupexclude = NULL) {

  protein_df <- metaData

  p <- selectData %>%
    left_join(select(protein_df, id, Description, Gene_name), by = "id") %>%
    mutate(Description = str_sub(Description, 1, 75)) %>%
    unite(Protein_name, Description, Gene_name, id, sep = " ") %>%
    as.data.frame() %>%
    column_to_rownames("Protein_name")

  #Exclude samples belonging to indicated groups
  sample_order1 <- sampleOrder
  group_order1 <- groupOrder

  if(isTruthy(groupexclude)){
    sample_order1 <- sample_order1[!group_order1 %in% groupexclude]
    group_order1 <- group_order1[!group_order1 %in% groupexclude]
  }

  #Rearrange based on sample table order
  p1 <- p[,sample_order1]

  #Column color annotations

  p1 %>%
    heatmaply(
      scale = "row",
      scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                                     midpoint = 0, name = "log2 FC"),
      Colv = FALSE,
      col_side_colors = group_order1

    )
}
