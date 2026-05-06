#' Runs the Shiny App2 Statistical Analysis
#' @param options A list of options passed to \code{shiny::shinyApp()}.
#' @import stringr
#' @import shinydashboard
#' @import rhandsontable
#' @import limma
#' @import shiny
#' @import cowplot
#' @importFrom plotly ggplotly layout plotlyOutput renderPlotly
#' @import heatmaply
#' @import UpSetR
#' @import ggsci
#' @export
runApp2 <- function(options = list()){
  ggplot2::theme_set(theme_cowplot())
  options(shiny.maxRequestSize = 5000*1024^2)
  
  
  
  

# Modules -----------------------------------------------------------------

  downloadablePlotUI <- function(id, label1) {
    ns <- NS(id)
    tabPanel(label1,
             fluidRow(column(3, textInput(ns("DLWidth"), "Width", value = 8)),
                      column(3, textInput(ns("DLHeight"), "Height", value = 6)),
                      column(3, textInput(ns("DLDpi"), "dpi", value = 300)),
                      column(3, textInput(ns("DLFileName"), "Filename", value = "Plot"))),
             fluidRow(column(3, sliderInput(ns("DLXLabels"), "X-axis label size", min = 4, max = 16, value = 10, step = 1)),
                      column(3, sliderInput(ns("DLXAngle"), "X-axis label angle", min = 0, max = 90, value = 45, step =5)),
                      column(3, sliderInput(ns("DLYLabels"), "Y-axis label size", min = 4, max = 16, value = 10, step = 1))),
             fluidRow(column(3, sliderInput(ns("DLXTitle"), "X-axis title size", min = 0, max = 25, value = 12, step = 1)),
                      column(3, sliderInput(ns("DLYTitle"), "Y-axis title size", min = 0, max = 25, value = 12, step = 1)),
                      column(3, sliderInput(ns("DLLegendTitle"), "Legend title size", min = 0, max = 25, value = 12, step = 1)),
                      column(3, sliderInput(ns("DLLegendText"), "Legend text size", min = 0, max = 25, value = 12, step = 1))),
             fluidRow(column(3, downloadButton(ns("DLButton")))),
             plotOutput(ns("DLPlot"), height = "900px"))
    
  }
  
  
  downloadablePlotServer <- function(id, plotfun) {
    moduleServer(
      id,
      function(input, output, session) {
        
        
        adjusted_plot <- reactive({
          x_size = as.numeric(input$DLXLabels)
          x_angle = as.numeric(input$DLXAngle)
          y_size = as.numeric(input$DLYLabels)
          x_title = as.numeric(input$DLXTitle)
          y_title = as.numeric(input$DLYTitle)
          legend_title_size = as.numeric(input$DLLegendTitle)
          legend_text_size = as.numeric(input$DLLegendText)
          
          plotfun() +
            theme(axis.text.x = element_text(size = x_size, angle = x_angle, hjust = 1),
                  axis.text.y = element_text(size = y_size),
                  axis.title.x = element_text(size = x_title),
                  axis.title.y = element_text(size = y_title),
                  legend.title = element_text(size = legend_title_size),
                  legend.text = element_text(size = legend_text_size))
        })
        
        output$DLPlot <- renderPlot({
          plotfun()
        })
        
        output$DLTest <- renderText({
          paste(input$DLFileName, '.tiff', sep='')
        })
        
        output$DLButton <- downloadHandler(
          filename = function() {paste(input$DLFileName, '.tiff', sep='') },
          content = function(file) {
            width1 = as.numeric(input$DLWidth)
            height1 = as.numeric(input$DLHeight)
            dpi1 = as.numeric(input$DLDpi)
            
            ggsave(file,
                   plot = adjusted_plot(),
                   width = width1,
                   height = height1,
                   dpi = dpi1,
                   device = "tiff",
                   compression = "lzw")
            
          },
          contentType = "image/tiff"
        )
        
      }
    )
  }  



  # Limma tab --------------------------------------------------------------

  #Deprecated
  # Box0 <- {
  #   fluidRow(
  #     box(
  #       width = 12,
  #       fileInput("allUpload", "Upload multiple files", multiple = TRUE)
  #     )
  #   )
  # }
  # 
  # Box1 <- {
  #   fluidRow(
  #     box(
  #       width = 12,
  #       fileInput("quantData", "Import quantitative data"),
  #       fileInput("metadata", "Import metadata"),
  #       fileInput("peptideQuantData", "Import peptide quantitative data"),
  #       fileInput("peptideMetaData", "Import peptide metadata")
  #     )
  #   )
  # }
  

# Module testing ----------------------------------------------------------


  file_upload_box <- {
    fluidRow(
      tabBox(id = "fileInputTab", width = 12,
             tabPanel("Upload multiple files",
                      fileInput("allUpload", "Upload multiple files", multiple = TRUE)),
             tabPanel("Upload individual files",
                      fileInput("quantData", "Import quantitative data"),
                      fileInput("metadata", "Import metadata"),
                      fileInput("peptideQuantData", "Import peptide quantitative data"),
                      fileInput("peptideMetaData", "Import peptide metadata"))
             )
    )
  }




  experimental_design_box <- {
    fluidRow(
      tabBox(
        id = "Samplegrouping",
        width = 5,
        title = "Sample grouping",
        tabPanel(
          "Table",
          fileInput("sampleFile", "Import sample file"),
          rHandsontableOutput("sample_table")
        ),
        tabPanel(
          "Buttons",
          selectizeInput("sampleGroupA", "Group A", multiple = TRUE, choices = character(0)),
          selectizeInput("sampleGroupB", "Group B", multiple = TRUE, choices = character(0)),
          actionButton("sampleActionButton", "All other samples")
        )

      ),
      box(
        width = 5,
        title = "Design matrix",
        checkboxInput("blockCheck", "Blocked design"),
        rHandsontableOutput("design_matrix_table")
      ),
      box(
        width = 2,
        title = "Contrast matrix",
        fileInput("contrastFile", "Import contrast list"),
        rHandsontableOutput("contrast_table")
      )

    )
  }

  test_results_table_box <- {
    fluidRow(
      box(
        width = 6,
        title = "test_results_table_1",
        tableOutput("test_results_table_1")
      ),
      box(
        width = 6,
        title = "test_results_table_2",
        tableOutput("test_results_table_2")
      )
    )
  }

  test_results_table_box2 <- {
    fluidRow(
      box(
        width = 6,
        title = "test_results_table_3",
        tableOutput("test_results_table_3")
      ),
      box(
        width = 6,
        title = "test_results_table_4",
        tableOutput("test_results_table_4")
      )
    )
  }

  limma_plot_box <- {
    fluidRow(
      box(
        width = 12,
        title = "test_summary_plot",
        plotOutput("test_summary_plot")
      )
    )
  }

  upset_diagram_box <- {
    fluidRow(
      box(
        width = 12,
        title = "P8 Venn",
        plotOutput("upset_plot",
                   height = "900px")
      )
    )
  }

  set_selection_box <- {
    fluidRow(
      box(
        width = 12,
        title = "P9 Set Selection",
        selectInput("set_selector", "Select sets", choices = NULL, selected = NULL,
                    multiple = TRUE, width = "50%"),
        tableOutput("selected_set_table"),
        actionButton("saveSets", "Extract and save all sets")
      )
    )
  }

  data_export_box <- {
    fluidRow(
      box(
        width = 12,
        title = "Save_output",
        column(4, actionButton("Save_output", "Save output")),
        column(4, actionButton("Save_summary", "Save limma and metadata"))

      )
    )
  }

  Limma_tab <- {
    tabItem("Limma",
            #Box0,
            # Box1,
            file_upload_box,
            experimental_design_box,
            # test_results_table_box,
            # test_results_table_box2,
            limma_plot_box,
            upset_diagram_box,
            set_selection_box,
            data_export_box
    )
  }



  # Visualize tab -----------------------------------------------------------

  Viz_volcano <- {
    fluidRow(
      box(
        width = 12,
        fluidRow(
          column(3, selectInput("ProteinVolcanoComparison", label = "Comparison",
                                choices = NULL,
                                multiple = TRUE,
                                selected = NULL)),
          column(1,
                 textInput("pv1FC", "FC cutoff:",
                           value = 0.2)),
          column(1,
                 textInput("pv1pvalue", "p-value cutoff:",
                           value = 1)),
          column(2,
                 checkboxInput("pv1useAdjP", "Use adjusted p.values", value = TRUE))

        ),
        fluidRow(
          column(6,
                 plotlyOutput("ProteinVolcanoplot"),
                 selectInput("ProteinSearch", label = "Highlight a protein:",
                             choices = character(0), multiple = TRUE)
          ),
          column(6,
                 box(
                   width = 12,
                   title = "GO Enrichment (selected proteins)",
                   checkboxInput("goEnrichCheck",
                                 "Run GO enrichment on selected proteins",
                                 value = FALSE),
                   fluidRow(
                     column(4,
                            selectInput("goOnt", "Ontology",
                                        choices  = c("Biological Process" = "BP",
                                                     "Molecular Function" = "MF",
                                                     "Cellular Component" = "CC",
                                                     "All"               = "ALL"),
                                        selected = "BP")),
                     column(4,
                            selectInput("goScoreType", "Score type",
                                        choices  = c("Standard" = "std",
                                                     "Positive" = "pos",
                                                     "Negative" = "neg"),
                                        selected = "std")),
                     column(4,
                            textInput("goPvalueCutoff", "p-value cutoff", value = "0.05"))
                   ),
                   fluidRow(
                     column(4,
                            textInput("goShowCategory", "Categories to show", value = "20")),
                     column(4,
                            selectInput("goGeneSet", "Heatmap gene set",
                                        choices  = c("Core enrichment" = "core",
                                                     "Full gene set"   = "full"),
                                        selected = "core"))
                   ),
                   tabBox(
                     id    = "goPlotTabs",
                     width = 12,
                     tabPanel("Dot plot",
                              plotlyOutput("goEnrichPlot", height = "500px")),
                     tabPanel("Network (cnetplot)",
                              plotOutput("goCnetplot", height = "600px")),
                     tabPanel("Heat plot",
                              plotOutput("goHeatplot", height = "500px")),
                     tabPanel("GSEA plot",
                              plotOutput("goGseaplot", height = "500px"))
                   )
                 )
          ))
      )
    )

  }

  Viz_heatmap <- {
    fluidRow(
      box(
        title = "Heat map of selected proteins",
        collapsible = TRUE,
        collapsed = FALSE,
        width = 12,
        fluidRow(
          column(2, checkboxInput("phcheck", label = "Plot heatmap")),
          column(2, checkboxInput("phscalecheck", label = "Scale Rows")),
          column(4, selectInput("proteinHeatGroupExclude", label = "Exclude groups:",
                                choices = character(0), multiple = TRUE))
        ),
        fluidRow(
          style = "height:900px",
          tabBox(title = "Selected proteins", id = "tab1", width = 12,
                 tabPanel("Select proteins",
                          plotlyOutput("proteinHeat", height = "900px")),
                 downloadablePlotUI("proteinDLPlot", "Protein downloadable heatmap"),
                 downloadablePlotUI("proteinDLPoint", "Protein downloadable PointPlot"),
                 tabPanel("Summarized",
                          plotlyOutput("ProteinHeatAverage", height = "900px")),
                 tabPanel("Peptide plotly heatmap",
                          plotlyOutput("peptideHeat", height = "900px"),
                          tableOutput("peptide1"),
                          tableOutput("peptide2")),
                 downloadablePlotUI("peptideDLPlot", "Peptide ggplot heatmap")
                 
          )
        )
        # fluidRow(
        #   tableOutput("Test1")
        # )
      )
    )
  }

  Viz_tab <- {
    tabItem("Viz",
            Viz_volcano,
            Viz_heatmap
    )
  }


  # UI ----------------------------------------------------------------------


  ui <- dashboardPage(
    dashboardHeader(title = "Proteoviz Model"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Limma", tabName = "Limma", icon = icon("th")),
        menuItem("Vizualize", tabName = "Viz", icon = icon("th"))
      )
    ),
    dashboardBody(
      tabItems(
        Limma_tab,
        Viz_tab
      )
    )
  )


  # Server ------------------------------------------------------------------

  server <- function(input, output, session) {


    ##################################### Limma ----------------------------------
    # Normalized data ---------------------------------------------------------
    

    quant_data <- reactive({
      if(isTruthy(input$allUpload)){
        a <- as_tibble(input$allUpload) %>%
          filter(grepl("protein\\_quantitative\\_data", name))
        if(nrow(a) == 1){
          return(read_tsv(a$datapath))
        }
        
      } 
      
      req(input$quantData)
      df <- input$quantData$datapath
      read_tsv(df)
      
    })
    

    
    peptide_quant_data <- reactive({
      if(isTruthy(input$allUpload)){
        a <- as_tibble(input$allUpload) %>%
          filter(grepl("peptide\\_quantitative\\_data", name))
        if(nrow(a) == 1){
          return(read_tsv(a$datapath))
        }
      } 
      
      req(input$peptideQuantData)
      df <- input$peptideQuantData$datapath
      read_tsv(df)
      
    })

    limma_input <- reactive({
      req(quant_data())
      quant_data() %>%
        as.data.frame() %>%
        column_to_rownames("id")
    })



    # Metadata ----------------------------------------------------------------

    # metadata <- reactive({
    #   req(input$metadata)
    #   if(isTruthy(input$metadata)){
    #     df <- input$metadata
    #     read_tsv(df$datapath)
    #   }
    # })
    # 
    # peptide_metadata <- reactive({
    #   req(input$peptideMetaData)
    #   if(isTruthy(input$peptideMetaData)){
    #     df <- input$peptideMetaData
    #     read_tsv(df$datapath)
    #   }
    # })

    metadata <- reactive({
      if(isTruthy(input$allUpload)){
        a <- as_tibble(input$allUpload) %>%
          filter(grepl("protein\\_metadata", name))
        if(nrow(a) == 1){
         return(read_tsv(a$datapath)) 
        }
        
      }
      
      req(input$metadata)
      df <- input$metadata
      read_tsv(df$datapath)

    })
    
    peptide_metadata <- reactive({
      if(isTruthy(input$allUpload)){
        a <- as_tibble(input$allUpload) %>%
          filter(grepl("peptide\\_metadata", name))
        if(nrow(a) == 1){
          return(read_tsv(a$datapath))
          }

      } 
      
      req(input$peptideMetaData)
      df <- input$peptideMetaData
      read_tsv(df$datapath)
      
    })
    
    observe({
      if(isTruthy(metadata())){
        x1 <- metadata() %>%
          select(Protein_Name)
        updateSelectInput(session, "ProteinSearch",
                          choices = x1)
      }
    })


    # Sample table ------------------------------------------------------------


    # sampleNameTable <- reactive({
    #   
    #   if (isTruthy(input$sampleFile)) {
    #     col_check <- c("Sample_name", "Group", "Batch", "Block")
    #     
    #     read_tsv(input$sampleFile$datapath) %>%
    #       select(one_of(col_check)) %>%
    #       mutate(Batch = as.integer(Batch)) %>%
    #       rhandsontable() %>%
    #       hot_col("Sample_name", readOnly = T)
    #   } else if(isTruthy(quant_data())){
    #     quant_column_names <- quant_data() %>%
    #       dplyr::select(-id) %>%
    #       colnames()
    #     tibble(Sample_name = quant_column_names,
    #            Group = "",
    #            Batch = "",
    #            Block = ""
    #     ) %>%
    #       as.data.frame() %>%
    #       rhandsontable() %>%
    #       hot_col("Sample_name", readOnly = T)
    #   }
    # 
    # 
    # })
    
    sample_table <- reactive({
      col_check <- c("Sample_name", "Group", "Batch", "Block")
      
      if(isTruthy(input$allUpload)){
        a <- as_tibble(input$allUpload) %>%
          filter(grepl("design\\_table", name))
        if(nrow(a) > 0){
          a1 <- read_tsv(a$datapath) %>%
            select(one_of(col_check)) %>%
            mutate(Batch = as.integer(Batch)) %>%
            rhandsontable() %>%
            hot_col("Sample_name", readOnly = T)
          return(a1)
        }
      }
      
      if (isTruthy(input$sampleFile)) {
        
        
        read_tsv(input$sampleFile$datapath) %>%
          select(one_of(col_check)) %>%
          mutate(Batch = as.integer(Batch)) %>%
          rhandsontable() %>%
          hot_col("Sample_name", readOnly = T)
      } else if(isTruthy(quant_data())){
        quant_column_names <- quant_data() %>%
          dplyr::select(-id) %>%
          colnames()
        tibble(Sample_name = quant_column_names,
               Group = "",
               Batch = "",
               Block = ""
        ) %>%
          as.data.frame() %>%
          rhandsontable() %>%
          hot_col("Sample_name", readOnly = T)
      }
      
      
    })

    # Outputs the Rhandsontable
    output$sample_table <- renderRHandsontable({
      req(sample_table())
      sample_table()
    })

    # Changes the handsontable back into a dataframe
    sample_table_df <- reactive({
      if(is.null(input$sample_table)) return (NULL)
      df <- hot_to_r(input$sample_table)
    })
    
    
    observe({
      if(isTruthy(sample_table_df())){
        x1 <- sample_table_df() %>%
          select(Sample_name)
        updateSelectInput(session, "sampleGroupA",
                          choices = x1$Sample_name)
        updateSelectInput(session, "sampleGroupB",
                          choices = x1$Sample_name)
      }
    })
    
    observe({
      req(input$sampleGroupA)
      req(sample_table_df())
      
      x1 <- sample_table_df() %>%
        select(Sample_name)
      x2 <- input$sampleGroupA
      x3 <- x1 %>%
        filter(!Sample_name %in% x2)
      
      updateSelectInput(session, "sampleGroupB",
                        choices = x3$Sample_name)
    })
    
    observeEvent(input$sampleActionButton, {
      req(sample_table_df())
      x1 <- sample_table_df() %>%
        select(Sample_name)
      x2 <- input$sampleGroupA
      x3 <- x1 %>%
        filter(!Sample_name %in% x2)
      
      updateSelectizeInput(session, "sampleGroupB",
                        choices = x3,
                        selected = x3$Sample_name)
      
    })




    # Design table ------------------------------------------------------------
    # design_table <- reactive({
    #   req(sample_table_df())
    #   makeDesignTable(sample_table_df())
    # })
    
    # design_table <- reactive({
    #   req(sample_table_df())
    #   if(isTruthy(input$sampleGroupA)){
    #     a <- sample_table_df() %>%
    #       mutate(Group = ifelse(Sample_name %in% input$sampleGroupA, "A", "B"))
    #     return(makeDesignTable(a))
    #   } else {
    #     makeDesignTable(sample_table_df())
    #   }
    #   
    # }) %>% throttle(500)
    
    design_table <- reactive({
      req(sample_table_df())
      a <- sample_table_df()
      
      if(isTruthy(input$sampleGroupA)){
        a <- a %>%
          mutate(Group = ifelse(Sample_name %in% input$sampleGroupA, "A", Group))
      }
      
      if(isTruthy(input$sampleGroupB)){
        a <- a %>%
          mutate(Group = ifelse(Sample_name %in% input$sampleGroupB, "B", Group))
      }
      
      if(all(isTruthy(input$sampleGroupB) & isTruthy(input$sampleGroupB))){
        a1 <- c(input$sampleGroupA, input$sampleGroupB)
        a <- a %>%
          mutate(Group = ifelse(!Sample_name %in% a1, "C", Group))
      }
      
      return(build_design_matrix(a))

    }) %>% debounce(200)

    sample_order <- reactive({
      req(sample_table_df())
      sample_table_df()$Sample_name
    })

    group_order <- reactive({
      req(sample_table_df())
      sample_table_df()$Group

    })


    #Not sure whether to make the design matrix editable or not.
    design_matrix_df <- reactive({
      if(is.null(input$design_matrix_table)) return (NULL)
      df <- hot_to_r(input$design_matrix_table)
      df
    })



    # Contrast table ----------------------------------------------------------

    contrast_table <- reactive({
      if(is.null(input$contrast_table)) return (NULL)
      df <- hot_to_r(input$contrast_table) %>%
        #Removes blank entries
        filter(Contrast_name != "")
      df
    })


    contrast_table_df <- reactive({
      req(design_matrix_df())
      req(contrast_table())

      try(build_contrast_matrix(design_matrix_df(), contrast_table()), silent = TRUE)
    })





    # Limma -------------------------------------------------------------------
    limma_output <- reactive({
      req(limma_input())
      req(design_matrix_df())
      req(sample_table_df())

      run_limma(limma_input(), design_matrix_df(), sample_table_df(), block = input$blockCheck)

    })

    contrast_fit <- reactive({
      req(limma_output())
      req(contrast_table_df())

      fit_contrasts(limma_output(), contrast_table_df())
    })

    combined_data <- reactive({

      req(contrast_fit())
      req(design_matrix_df())
      req(contrast_table())

      combine_limma_results(contrast_fit(), design_matrix_df(), contrast_table_df())
    })

    test_results <- reactive({
      req(contrast_fit())

      decideTests(contrast_fit(), lfc = log2(1.5))

    })

    tidy_test_results <- reactive({
      req(test_results())

      make_tidy_test_results(test_results())
    })

    test_results_list <- reactive({
      req(tidy_test_results())
      make_test_results_list(tidy_test_results())

    })


    # Table plots -------------------------------------------------------------
    output$sample_table_preview <- renderTable({
      req(quant_data())
      quant_data() %>%
        dplyr::slice(1:10)
    })

    output$design_matrix_table <- renderRHandsontable({
      req(design_table())

      sample_names <- rownames(design_table())

      df3 <- apply(design_table(), 2, as.integer) %>%
        as.data.frame()

      rownames(df3) <- sample_names

      rhandsontable(df3)

    })


    output$contrast_table <- renderRHandsontable({
      
      if(isTruthy(input$allUpload)){
        a <- as_tibble(input$allUpload) %>%
          filter(grepl("contrast\\_list", name))
        if(nrow(a) == 1){
          a1 <- read_tsv(a$datapath) %>%
            as.data.frame() %>%
            rhandsontable()
          return(a1)
        }
      }
      if(isTruthy(input$contrastFile)){
        #File import
        df <- input$contrastFile$datapath
        read_tsv(df) %>%
          as.data.frame() %>%
          rhandsontable()
      } else {
        #Manual entry
        req(sample_table_df())
        tibble(Contrast_name = "",
        ) %>%
          as.data.frame() %>%
          rhandsontable()
      }

    })



    # Save output -------------------------------------------------------------
    observeEvent(input$Save_output, {
      req(quant_data())
      req(combined_data())
      req(metadata())
      write_limma_summary(quant_data(), combined_data(), metadata(), prefix = "")
    })

    observeEvent(input$Save_summary, {
      req(quant_data())
      req(combined_data())
      req(metadata())
      write_limma_summary(quant_data(), combined_data(), metadata(), prefix = "Summarized")
    })

    observeEvent(input$saveSets, {
      req(tidy_test_results())
      req(metadata())

      meta_df <- metadata() %>%
        mutate(id = as.character(id))

      unnest(tidy_test_results(), data) %>%
        group_by(protein) %>%
        summarize(Sets = str_c(Test, collapse = ", "),
                  Number_of_groups = length(Test)) %>%
        dplyr::select(Number_of_groups, Sets, protein) %>%
        ungroup() %>%
        mutate(protein = as.character(protein)) %>%
        left_join(meta_df, by = c("protein" = "id")) %>%
        arrange(desc(Number_of_groups)) %>%
        write_tsv("Extracted_sets.tsv")
    })





    # Debugging ---------------------------------------------------------------

    output$test_results_table_1 <- renderTable({
      req(contrast_table_df())
      contrast_table_df()
    }, rownames = TRUE)

    output$test_results_table_2 <- renderTable({
      req(limma_output())

      head(limma_output()$coefficients)
    })

    output$test_results_table_3 <- renderTable({
      req(combined_data())

      head(combined_data())
    })

    output$test_results_table_4 <- renderTable({
      req(contrast_fit())
      req(contrast_table())
      x1 <- contrast_table()$Contrast_name %>%
        .[[1]]
      head(topTable(contrast_fit(), x1, number = Inf))

    })

    output$test_summary_plot <- renderPlot({
      req(combined_data())

      df <- combined_data()

      df %>%
        ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point() +
        geom_hline(yintercept = -log10(.05), linetype = "dotted") +
        geom_vline(xintercept = -log2(1.5), linetype = "dotted") +
        geom_vline(xintercept = log2(1.5), linetype = "dotted") +
        facet_wrap(~Comparison, scales = "free_y")
    })

    output$upset_plot <- renderPlot({

      req(tidy_test_results())

      a2 <- map(tidy_test_results()$data, "protein")
      names(a2) <- tidy_test_results()$Test


      upset(fromList(a2), nsets = length(a2),
            text.scale = 1.5,
            point.size = 4)


    })

    observe({
      req(tidy_test_results())

      df <- tidy_test_results()

      x2 <- unique(df$Test)

      updateSelectInput(session, "set_selector",
                        choices = x2)
    })
    

    output$selected_set_table <- renderTable({
      req(input$set_selector)
      req(tidy_test_results())
      req(metadata())

      queries <- input$set_selector

      a3 <- tidy_test_results() %>%
        unnest(data)

      # a4 <- a3 %>%
      #   group_by(id) %>%
      #   mutate(Number_of_groups = length(Test)) %>%
      #   filter(Test %in% queries,
      #          Number_of_groups == length(queries)) %>%
      #   ungroup() %>%
      #   distinct(id)

      a4 <- a3 %>%
        group_by(protein) %>%
        mutate(Number_of_groups = length(Test),
               Matches_query = Test %in% queries) %>%
        filter(all(Matches_query),
               Number_of_groups == length(queries)) %>%
        ungroup() %>%
        distinct(protein)

      meta_df <- metadata() %>%
        mutate(id = as.character(id))

      a5 <- a4 %>%
        mutate(protein = as.character(protein)) %>%
        left_join(meta_df, by = c("protein" = "id"))

      a5


    })





    # Viz ---------------------------------------------------------------------

    observe({
      req(contrast_table())

      x1 <- contrast_table()$Contrast_name

      updateSelectInput(session, "ProteinVolcanoComparison",
                        choices = x1)
    })

    observe({
      req(group_order())
      x1 <- unique(group_order())
      updateSelectInput(session, "proteinHeatGroupExclude",
                        choices = x1)
    })

    excluded_groups <- reactive({
      inp <- input$proteinHeatGroupExclude
      if (is.null(inp)) character(0) else inp
    })

    output$ProteinVolcanoplot <- renderPlotly({
      req(combined_data())
      req(input$ProteinVolcanoComparison)

      input$pv1button
      volcanoComparison <- isolate(input$ProteinVolcanoComparison)
      FCcutoff  <- isolate(as.numeric(input$pv1FC))
      Pcutoff   <- isolate(as.numeric(input$pv1pvalue))
      useAdjP   <- isTRUE(input$pv1useAdjP)

      # Join protein metadata so hover text can show description/gene/uniprot.
      # Match on protein (character rowname from matrix) == id (from metadata).
      plot_data <- combined_data() |>
        dplyr::filter(Comparison %in% volcanoComparison)

      if (isTruthy(metadata())) {
        meta_slim <- metadata() |>
          dplyr::mutate(protein = as.character(id)) |>
          dplyr::select(protein,
                        dplyr::any_of(c("PG.ProteinDescriptions",
                                        "PG.Genes",
                                        "PG.UniProtIds")))
        plot_data <- dplyr::left_join(plot_data, meta_slim, by = "protein")
      }

      plot_volcano(
        plot_data,
        fc_threshold  = FCcutoff,
        p_threshold   = Pcutoff,
        use_adj_p     = useAdjP,
        plotly_source = "proteinVolcano",
        key_col       = "protein"
      )
    })

    selected_data <- reactive({
      req(quant_data())
      d <- plotly::event_data("plotly_selected", source = "proteinVolcano")
      req(!is.null(d))

      quant_data() %>%
        filter(as.character(id) %in% d$key)
    })

    protein_click_data <- reactive({
      req(quant_data())
      d <- plotly::event_data("plotly_click", source = "proteinVolcano")
      req(!is.null(d))

      quant_data() %>%
        filter(as.character(id) %in% d$key)
    })

    peptide_click_data <- reactive({
      req(quant_data())
      req(metadata())
      req(peptide_metadata())
      d <- plotly::event_data("plotly_click", source = "proteinVolcano")
      req(!is.null(d))

      protein_id <- metadata() %>%
        filter(as.character(id) %in% d$key)
      
      #If DIA peptide metadata
      if("Protein_Accessions" %in% names(peptide_metadata())){
        peptide_metadata() %>%
          filter(Protein_Accessions %in% protein_id$Accession_Number)
      } else if("PG.ProteinAccessions" %in% names(peptide_metadata())){
        #For Spectronaut data
        # peptide_metadata() %>%
        #   filter(PG.ProteinAccessions %in% protein_id$PG.ProteinAccessions)
        
        peptide_metadata() %>%
          filter(PG.ProteinAccessions %in% protein_id$PG.ProteinGroups)
        
      }else {
        #Else, assumes MaxQuant peptide metadata
        peptide_metadata() %>%
          filter(Protein_group_IDs %in% protein_id$id)
      }
    })
    
    
    
    
    peptide_click_table <- reactive({
      req(peptide_click_data())
      req(peptide_quant_data())
      
      peptide_quant_data() %>%
        filter(id %in% peptide_click_data()$id)
      
    })
    
    output$peptide1 <- renderTable({
      req(peptide_click_data())
      req(metadata())
      peptide_click_data()
      
    })
    
    output$peptide2 <- renderTable({
      req(peptide_click_table())
      peptide_click_table()
      
    })

    output$proteinHeat <- renderPlotly({
      req(input$phcheck)
      req(quant_data())
      req(sample_table_df())

      # GO dotplot click takes priority; fall back to volcano selection
      heat_ids <- if (isTruthy(go_click_proteins()) && nrow(go_click_proteins()) > 0) {
        go_click_proteins()$id
      } else {
        req(selected_data())
        selected_data()$id
      }

      quant_long <- quant_data() |>
        tidyr::pivot_longer(-id, names_to = "Sample_name", values_to = "Intensity")

      plot_protein_heatmap(
        quant_long,
        sample_table_df(),
        proteins       = heat_ids,
        protein_anno   = metadata(),
        scale_rows     = input$phscalecheck,
        exclude_groups = excluded_groups()
      )
    })

    output$peptideHeat <- renderPlotly({
      req(input$phcheck)
      req(peptide_click_table())
      req(peptide_metadata())
      req(protein_click_data())
      req(metadata())
      req(sample_table_df())

      peptide_long <- peptide_click_table() |>
        tidyr::pivot_longer(-id, names_to = "Sample_name", values_to = "Intensity")

      protein_names <- metadata() |>
        dplyr::filter(id %in% protein_click_data()$id) |>
        dplyr::pull(Protein_Name)

      plot_peptide_heatmap(
        peptide_long,
        sample_table_df(),
        proteins         = protein_names,
        peptide_metadata = peptide_metadata()
      )
    })

    downloadablePlotServer("peptideDLPlot", reactive({
      req(input$phcheck)
      req(peptide_click_table())
      req(peptide_metadata())
      req(protein_click_data())
      req(metadata())
      req(sample_table_df())

      peptide_long <- peptide_click_table() |>
        tidyr::pivot_longer(-id, names_to = "Sample_name", values_to = "Intensity")

      protein_names <- metadata() |>
        dplyr::filter(id %in% protein_click_data()$id) |>
        dplyr::pull(Protein_Name)

      plot_peptide_heatmap(
        peptide_long,
        sample_table_df(),
        proteins         = protein_names,
        peptide_metadata = peptide_metadata()
      )
    }))


    downloadablePlotServer("proteinDLPlot", reactive({
      req(input$phcheck)
      req(quant_data())
      req(sample_table_df())
      req(selected_data())

      quant_long <- quant_data() |>
        tidyr::pivot_longer(-id, names_to = "Sample_name", values_to = "Intensity")

      plot_protein_heatmap(
        quant_long,
        sample_table_df(),
        proteins       = selected_data()$id,
        protein_anno   = metadata(),
        scale_rows     = input$phscalecheck,
        exclude_groups = excluded_groups()
      )
    }))
    
    
    downloadablePlotServer("proteinDLPoint", reactive({
      req(input$phcheck)
      req(protein_click_data())
      req(metadata())
      req(sample_order())
      req(group_order())
      
      makeProteinPointPlot(protein_click_data(), metadata(), sample_order(), group_order(), excluded_groups(), input$phscalecheck)
    }))
    
    # output$Test1 <- renderTable({
    #   req(input$phcheck)
    #   req(protein_click_data())
    #   req(metadata())
    #   req(sample_order())
    #   req(group_order())
    #   
    #   f_Test1(protein_click_data(), metadata(), sample_order(), group_order(), excluded_groups(), input$phscalecheck)
    # })
    
    
    
    
    
    
    


    # GO Enrichment -----------------------------------------------------------

    go_result <- reactive({
      req(input$goEnrichCheck)
      req(quant_data())
      req(metadata())
      req(selected_data())

      run_go_enrichment(
        quant_data(), metadata(), selected_data(),
        ont           = input$goOnt,
        score_type    = input$goScoreType,
        pvalue_cutoff = as.numeric(input$goPvalueCutoff)
      )
    })

    # Unpack ego and gene_list from go_result(); return NULL if not yet computed
    go_ego <- reactive({
      r <- tryCatch(go_result(), error = function(e) NULL)
      if (is.null(r)) NULL else r$ego
    })

    go_gene_list <- reactive({
      r <- tryCatch(go_result(), error = function(e) NULL)
      if (is.null(r)) NULL else r$gene_list
    })

    go_show_category <- reactive({
      n <- suppressWarnings(as.integer(input$goShowCategory))
      if (is.na(n) || n < 1) 20L else n
    })

    output$goEnrichPlot <- renderPlotly({
      plot_go_dotplot(go_ego(),
                      show_category = go_show_category(),
                      plotly_source = "goDotplot")
    })

    output$goCnetplot <- renderPlot({
      plot_go_cnetplot(go_ego(), go_gene_list(),
                       show_category = go_show_category())
    })

    output$goHeatplot <- renderPlot({
      plot_go_heatplot(go_ego(), go_gene_list(),
                       show_category = go_show_category())
    })

    output$goGseaplot <- renderPlot({
      d     <- plotly::event_data("plotly_click", source = "goDotplot")
      go_id <- if (!is.null(d)) d$key else NULL
      plot_go_gseaplot(go_ego(), go_id)
    })

    go_click_proteins <- reactive({
      result <- tryCatch(go_result(), error = function(e) NULL)
      if (is.null(result)) return(NULL)
      ego <- result$ego

      d <- plotly::event_data("plotly_click", source = "goDotplot")
      if (is.null(d)) return(NULL)

      req(quant_data())
      req(metadata())

      go_id    <- d$key
      gene_col <- intersect(c("PG.Genes", "Gene_name"), names(metadata()))[1]

      gene_symbols <- if (input$goGeneSet == "core") {
        ego@result |>
          dplyr::filter(ID == go_id) |>
          dplyr::pull(core_enrichment) |>
          stringr::str_split("/") |>
          unlist()
      } else {
        ego@geneSets[[go_id]]
      }

      protein_ids <- metadata() |>
        dplyr::mutate(gene = stringr::str_trim(
          stringr::str_extract(!!rlang::sym(gene_col), "^[^;]+")
        )) |>
        dplyr::filter(gene %in% gene_symbols) |>
        dplyr::pull(id)

      quant_data() |>
        dplyr::filter(id %in% protein_ids)
    })

    # -------------------------------------------------------------------------

# Debugging/scraps/tests --------------------------------------------------




  }

  shinyApp(ui, server, options = options)
}

