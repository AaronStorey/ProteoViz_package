#' Runs the Shiny DIA App
#' @param options A list of options passed to \code{shiny::shinyApp()}.
#' @import dplyr
#' @import purrr
#' @import readr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import shinydashboard
#' @import rhandsontable
#' @import limma
#' @import shiny
#' @import cowplot
#' @importFrom plotly ggplotly layout plotlyOutput renderPlotly
#' @export
runApp1DIA <- function(options = list()){
  theme_set(theme_cowplot())
  options(shiny.maxRequestSize = 5000*1024^2)
  type <- "DIA"


  # UI ----------------------------------------------------------------------


  # Import tab --------------------------------------------------------------




  import_files_box <- {
    fluidRow(
      box(
        width = 12,
        textInput("projectName", "Project name", value = "Demo", width = "50%"),
        fileInput("proteinFile", "Protein Samples Report"),
        fileInput("peptideFile", "Peptide Quant Report"),
        fileInput("sampleFile", "Import sample file"),
        rHandsontableOutput("sample_table"),
        br(),
        actionButton("saveSampleNameTable", "Save Table"),
        tableOutput("Output1"),
        textOutput("Output2"),
        tableOutput("Output3")
      )
    )
  }

  # Box2 <- {
  #   fluidRow(
  #     box(
  #       width = 12,
  #       fileInput("peptidesFile", "Peptides (TODO)")
  #     )
  #   )
  # }

  Import_tab <- {
    tabItem("Import",
            import_files_box
            #Box2
    )
  }


  # Protein tab -------------------------------------------------------------



  protein_qc_controls_box <- {
    fluidRow(
      box(
        width = 12,
        fluidRow(
          column(width= 3,
                 selectInput("samples_to_exclude", "Samples to exclude:", choices = NULL, multiple = TRUE, selectize = TRUE)
          ),
          column(width = 3,
                 selectInput("groups_to_exclude", "Groups to exclude:", choices = NULL, multiple = TRUE, selectize = TRUE)
          ),
          column(width = 2,
                 checkboxInput("P1log2", "log2 transform", value = TRUE),
                 checkboxInput("P1remove", "Remove 0 values", value = TRUE)
          ),
          column(width = 2,
                 checkboxInput("P1Normcheckbox", "Normalize to Reference channels"),
                 checkboxInput("P1checkbox", "Plot output", value = TRUE))
        ),
        plotOutput("pre_norm_boxplot")
      )
    )
  }


  group_observation_criteria_box <- {
    fluidRow(
      box(
        title = "Required number of observations for each group:",
        width = 12,
        fluidRow(column(12, checkboxInput("NaP3checkbox", "Plot missing values", value = FALSE))),
        fluidRow(column(12, uiOutput("group_sliders"))),
        plotOutput("missing_value_plot")
      )
    )
  }

  global_observation_criteria_box <- {
    fluidRow(
      box(
        width = 6,
        title = "Global group requirements",
        fluidRow(
          column(3, sliderInput("GlobalGroupX", "At least X observations", min = 0, max = 10, value = 0, step = 1)),
          column(3, sliderInput("GlobalGroupY", "In at least Y groups", min = 0, max = 10, value = 0, step = 1))
        )
      )
    )
  }




  normalization_controls_box <- {
    fluidRow(
      box(
        width = 12,
        fluidRow(
          column(
            width = 4,
            selectInput("normOption", "Normalization method",
                        choices = c("None", "Mean", "Median", "Quantile", "Cyclic Loess", "Cyclic Loess then Median", "Rlr", "Gi")
            )
          ),
          column(
            width = 2,
            checkboxInput("normCheckbox", "Run normalization",
            )
          )

        ),
        tabBox(
          width = 12,
          tabPanel(
            "Boxplot",
            plotOutput("post_norm_boxplot")
          ),
          tabPanel(
            "PCA",
            plotlyOutput("pca_plot")
          ),
          tabPanel(
            "Correlation",
            checkboxInput("CorP3check", "Scale by row"),
            plotOutput("correlation_plot")
          )
        )
      )
    )
  }



  protein_export_box <- {
    fluidRow(
      box(
        column(3, actionButton("P1dfbutton", "Save quantitative table")),
        column(3, actionButton("saveProteinMeta", "Save protein metadata")),
        column(3, actionButton("saveDesignTable", "Save sample design table")),
        column(3, actionButton("makeContrastList", "Make contrast list"))
      )
    )
  }
  
  peptide_export_box <- {
    fluidRow(
      box(
        width = 6,
        column(6, actionButton("P1peptidedfbutton", "Save peptide quantitative table")),
        column(6, actionButton("savePeptideMeta", "Save peptide metadata"))
      )
    )
  }

  ProteinNorm_tab <- {
    tabItem("ProteinNorm",
            protein_qc_controls_box,
            group_observation_criteria_box,
            global_observation_criteria_box,
            normalization_controls_box,
            protein_export_box,
            peptide_export_box)
  }


  ui <- dashboardPage(
    dashboardHeader(title = "Import ScaffoldDIA"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Import", tabName = "Import", icon = icon("th")),
        menuItem("Protein Normalization", tabName = "ProteinNorm", icon = icon("th"))
      )
    ),
    dashboardBody(
      tabItems(
        Import_tab,
        ProteinNorm_tab
      )
    )
  )


  # Server ------------------------------------------------------------------

  server <- function(input, output, session) {


    # Protein Server Section --------------------------------------------------
    # Sample table ------------------------------------------------------------
    # Credit: Steffan's ProteoNorm app
    # Creates handsontable where metaData1() will be edited

    sample_table <- reactive({
      if (isTruthy(input$sampleFile)) {
        read_sample_name_table(input$sampleFile$datapath) %>%
          rhandsontable() %>%
          hot_col("Data_name", readOnly = T)
      } else if(isTruthy(protein_data())){
        x1 <- make_sample_name_table(protein_data(), type = type) %>%
          as.data.frame() %>%
          rhandsontable() %>%
          hot_col("Data_name", readOnly = T)
        return(x1)
      }


    })

    # Outputs the Rhandsontable
    output$sample_table <- renderRHandsontable({
      req(sample_table())
      if(isTruthy(sample_table())){
        sample_table()
      }
    })

    # Changes the handsontable back into a dataframe
    sample_table_df <- reactive({
      if(is.null(input$sample_table)) return (NULL)
      req(input$sample_table)
      df <- hot_to_r(input$sample_table)
      df
    })

    # Removes excluded samples and groups from table
    sample_table_filtered <- reactive({
      req(sample_table_df())
      req(!is.null(sample_table_df()))
      req(length(sample_table_df()$Group) > 0)
      req(any(!is.na(sample_table_df()$Group)))
      x1 <- input$samples_to_exclude
      x2 <- input$groups_to_exclude

      sample_table_df() %>%
        filter(!Sample_name %in% x1,
               !Group %in% x2)
    })


    # Observers to update selection box
    observe({
      req(sample_table_df())

      df <- sample_table_df()

      x1 <- df$Sample_name
      x2 <- unique(df$Group)

      updateSelectInput(session, "samples_to_exclude",
                        choices = x1)
      updateSelectInput(session, "groups_to_exclude",
                        choices = x2)
    })

    # Save table
    observeEvent(input$saveSampleNameTable, {
      req(input$projectName)
      project_name <- input$projectName
      if(!dir.exists(project_name)) dir.create(project_name)

      sample_table_df() %>%
        write_tsv(paste0(project_name, "/", project_name, "_", "Sample_name_table.tsv"))
    })


    # Group tables ------------------------------------------------------------------
    #Counts the number of Samples in each group, after excluding specified groups/samples
    group_counts <- reactive({
      req(sample_table_filtered())
      req(!is.null(sample_table_filtered()))
      req(length(unique(sample_table_filtered()$Group)) > 0)
      sample_table_filtered() %>%
        count(Group, name = "Number_of_samples") %>%
        ungroup()
    })

    #Table that lists required number of observations for each group,
    #as specified from the input sliders

    group_observation_requirements <- reactive({
      req(group_counts())
      req(!is.null(group_counts()))

      group_tbl <- group_counts()

      slider_data <- map(group_tbl$Group, ~input[[paste0("slider", .x)]]) %>%
        unlist() %>%
        as.integer()

      if(length(slider_data) != nrow(group_tbl)) return (NULL)

      x1 <- tibble::tibble(Group = group_tbl$Group, req_observations = slider_data)

      make_observation_required_table(group_tbl, x1)

    })


    #Update global requirement slider
    observe({
      req(group_counts())

      Biggest_group <- max(group_counts()$Number_of_samples)
      Number_of_groups <- length(unique(group_counts()$Group))

      updateSliderInput(session, "GlobalGroupX",
                        max = Biggest_group)
      updateSliderInput(session, "GlobalGroupY",
                        max = Number_of_groups)
    })


    # Protein dfs -------------------------------------------------------------

    #Import
    protein_data <- reactive({
      if(isTruthy(input$proteinFile)){
        x1 <- input$proteinFile
        req(input$proteinFile)
        if(grepl("quant\\_report", x1$name)){
          read_quant_report(x1$datapath, type)
        } else read_samples_report(x1$datapath, type)

      } else if(isTruthy(input$peptideFile)){

        x1 <- input$peptideFile
        req(input$peptideFile)
        read_peptide_report(x1$datapath, type)
      }

    })
    
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

    #Exclude samples, tidy, transform (log2, remove 0s, NormalizeRefs)
    protein_data_clean <- reactive({
      req(protein_data())
      req(sample_table_filtered())

      clean_protein_data(protein_data(), sample_table_filtered(),
                         log2_transform = input$P1log2,
                         remove_zeros   = input$P1remove)

    })

    #Filter for proteins that meet observation criteria
    protein_data_filtered <- reactive({
      req(protein_data_clean())
      req(group_observation_requirements())
      req(input$GlobalGroupX)
      req(input$GlobalGroupY)

      filter_missing_values(protein_data_clean(), group_observation_requirements(),
                            min_per_group = input$GlobalGroupX,
                            min_groups    = input$GlobalGroupY)

    })

    #To matrix (retained for any downstream use; not passed to normalization)
    protein_data_matrix <- reactive({
      req(protein_data_filtered())
      protein_df_to_matrix(protein_data_filtered())
    })

    #Normalize — normalize_intensities() takes and returns long-form data
    protein_data_normalized <- reactive({
      req(protein_data_filtered())
      req(input$normOption)
      req(input$normCheckbox == TRUE)

      normalize_intensities(protein_data_filtered(), method = input$normOption)

    })

    #Long-form normalized data — same reactive, aliased for clarity
    protein_data_long <- reactive({
      req(protein_data_normalized())
      protein_data_normalized()

    })



    # Plots -------------------------------------------------------------------
    output$pre_norm_boxplot <- renderPlot({
      req(protein_data_clean())
      req(input$P1checkbox == TRUE)
      req(sample_table_filtered())

      plot_intensity_distributions(protein_data_clean(), sample_table_filtered())
    })

    output$post_norm_boxplot <- renderPlot({
      req(protein_data_long())
      req(sample_table_filtered())
      plot_intensity_distributions(protein_data_long(), sample_table_filtered())
    })

    output$pca_plot <- renderPlotly({
      req(protein_data_long())
      req(sample_table_filtered())

      plot_pca(protein_data_long(), sample_table_filtered())

    })

    output$correlation_plot <- renderPlot({
      req(protein_data_long())
      req(sample_table_filtered())
      plot_correlation_heatmap(protein_data_long(), sample_table_filtered(), scale_rows = input$CorP3check)
    })

    output$missing_value_plot <- renderPlot({
      req(protein_data_clean())
      req(sample_table_filtered())
      if(input$NaP3checkbox){
        plot_missing_values(protein_data_clean(), sample_table_filtered())
      }
    })


    # Extra -------------------------------------------------------------------

    #Reactive rendering of slider inputs
    output$group_sliders <- renderUI({
      req(group_counts())
      df <- group_counts()
      a <- list()
      for(i in seq_along(df$Group)){
        a[[i]] <- column(3,
                         sliderInput(paste0("slider", df$Group[[i]]), df$Group[[i]],
                                     min = 0, max = df$Number_of_samples[[i]],
                                     value = 0,
                                     step = 1,
                                     width = "100%"))
      }
      tagList(a)

    })




    #Save df
    observeEvent(input$P1dfbutton, {
      req(input$projectName)
      req(protein_data_long())
      req(sample_table_filtered())
      save_quantitative_data(protein_data_long(), input$projectName, sample_table_filtered(), "protein")
    })
    
    observeEvent(input$P1peptidedfbutton, {
      req(input$projectName)
      req(protein_data_long())
      req(sample_table_filtered())
      save_quantitative_data(protein_data_long(), input$projectName, sample_table_filtered(), "peptide")
    })

    observeEvent(input$saveProteinMeta, {
      req(protein_data())
      req(input$projectName)

      save_protein_meta(protein_data(), input$projectName, type)
    })
    
    observeEvent(input$savePeptideMeta, {
      req(protein_data())
      req(input$projectName)
      
      save_peptide_meta(protein_data(), input$projectName, type)
    })

    observeEvent(input$saveDesignTable, {
      req(sample_table_filtered())
      req(input$projectName)
      save_design_table(sample_table_filtered(), input$projectName)
    })

    observeEvent(input$makeContrastList, {
      req(sample_table_filtered())
      req(input$projectName)
      make_contrast_list(sample_table_filtered(), input$projectName)
    })



    output$Output1 <- renderTable(head(protein_data()))
    output$Output2 <- renderText({
      req(input$proteinFile)
      input$proteinFile$datapath
    })
    output$Output3 <- renderTable({
      req(input$proteinFile)
      input$proteinFile
    })
  }

  shinyApp(ui, server, options = options)

}
