#' Runs the Shiny DIA App
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
#' @import plotly
#' @export
runApp1DIA <- function(options = list()){
  theme_set(theme_cowplot())
  options(shiny.maxRequestSize = 5000*1024^2)
  type <- "DIA"


  # UI ----------------------------------------------------------------------


  # Import tab --------------------------------------------------------------




  Box1 <- {
    fluidRow(
      box(
        width = 12,
        textInput("projectName", "Project name", value = "Demo", width = "50%"),
        fileInput("proteinFile", "Protein Samples Report"),
        fileInput("peptideFile", "Peptide Quant Report"),
        fileInput("sampleFile", "Import sample file"),
        rHandsontableOutput("sampleNameTable"),
        br(),
        actionButton("saveSampleNameTable", "Save Table")
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
            Box1
            #Box2
    )
  }


  # Protein tab -------------------------------------------------------------



  Box3 <- {
    fluidRow(
      box(
        width = 12,
        fluidRow(
          column(width= 3,
                 selectInput("P1exclude", "Samples to exclude:", choices = NULL, multiple = TRUE, selectize = TRUE)
          ),
          column(width = 3,
                 selectInput("P1groupexclude", "Groups to exclude:", choices = NULL, multiple = TRUE, selectize = TRUE)
          ),
          column(width = 2,
                 checkboxInput("P1log2", "log2 transform", value = TRUE),
                 checkboxInput("P1remove", "Remove 0 values", value = TRUE)
          ),
          column(width = 2,
                 checkboxInput("P1Normcheckbox", "Normalize to Reference channels"),
                 checkboxInput("P1checkbox", "Plot output", value = TRUE))
        ),
        plotOutput("P1")
      )
    )
  }


  Box4 <- {
    fluidRow(
      box(
        title = "Required number of observations for each group:",
        width = 12,
        fluidRow(column(12, checkboxInput("NaP3checkbox", "Plot missing values", value = FALSE))),
        fluidRow(column(12, uiOutput("groupSlider"))),
        plotOutput("NaP3")
      )
    )
  }

  Box4_5 <- {
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




  Box5 <- {
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
            plotOutput("P2")
          ),
          tabPanel(
            "PCA",
            plotlyOutput("P3")
          ),
          tabPanel(
            "Correlation",
            checkboxInput("CorP3check", "Scale by row"),
            plotOutput("CorP3")
          )
        )
      )
    )
  }



  Box6 <- {
    fluidRow(
      box(
        column(3, actionButton("P1dfbutton", "Save quantitative table")),
        column(3, actionButton("saveProteinMeta", "Save protein metadata")),
        column(3, actionButton("saveDesignTable", "Save sample design table")),
        column(3, actionButton("makeContrastList", "Make contrast list"))
      )
    )
  }
  
  Box7 <- {
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
            Box3,
            Box4,
            Box4_5,
            Box5,
            Box6,
            Box7)
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

    sampleNameTable <- reactive({

      if (isTruthy(input$sampleFile)) {
        readSampleNameTable(input$sampleFile$datapath) %>%
          rhandsontable() %>%
          hot_col("Data_name", readOnly = T)
      } else if (isTruthy(protein_df())){
        makeSampleNameTable(protein_df(), type = type) %>%
          as.data.frame() %>%
          rhandsontable() %>%
          hot_col("Data_name", readOnly = T)
      }
    })

    # Outputs the Rhandsontable
    output$sampleNameTable <- renderRHandsontable({
      req(sampleNameTable())
      sampleNameTable()
    })

    # Changes the handsontable back into a dataframe
    sampleNameTable2 <- reactive({
      if(is.null(input$sampleNameTable)) return (NULL)
      df <- hot_to_r(input$sampleNameTable)
      df
    })

    # Removes excluded samples and groups from table
    sampleNameTable3 <- reactive({
      req(sampleNameTable2())
      req(!is.null(sampleNameTable2()))
      req(length(sampleNameTable2()$Group) > 0)
      req(any(!is.na(sampleNameTable2()$Group)))
      x1 <- input$P1exclude
      x2 <- input$P1groupexclude

      sampleNameTable2() %>%
        filter(!Sample_name %in% x1,
               !Group %in% x2)
    })


    # Observers to update selection box
    observe({
      req(sampleNameTable2())

      df <- sampleNameTable2()

      x1 <- df$Sample_name
      x2 <- unique(df$Group)

      updateSelectInput(session, "P1exclude",
                        choices = x1)
      updateSelectInput(session, "P1groupexclude",
                        choices = x2)
    })

    # Save table
    observeEvent(input$saveSampleNameTable, {
      req(input$projectName)
      project_name <- input$projectName
      if(!dir.exists(project_name)) dir.create(project_name)

      sampleNameTable2() %>%
        write_tsv(paste0(project_name, "/", project_name, "_", "Sample_name_table.tsv"))
    })


    # Group tables ------------------------------------------------------------------
    #Counts the number of Samples in each group, after excluding specified groups/samples
    groupTable <- reactive({
      req(sampleNameTable3())
      req(!is.null(sampleNameTable3()))
      req(length(unique(sampleNameTable3()$Group)) > 0)
      sampleNameTable3() %>%
        count(Group, name = "Number_of_samples") %>%
        ungroup()
    })

    #Table that lists required number of observations for each group,
    #as specified from the input sliders

    groupTable2 <- reactive({
      req(groupTable())
      req(!is.null(groupTable()))

      group_tbl <- groupTable()

      slider_data <- map(group_tbl$Group, ~input[[paste0("slider", .x)]]) %>%
        unlist() %>%
        as.integer()

      if(length(slider_data) != nrow(group_tbl)) return (NULL)

      group_tbl$req_observations <- slider_data

      #group_tbl
      group_tbl

    })


    #Update global requirement slider
    observe({
      req(groupTable())

      Biggest_group <- max(groupTable()$Number_of_samples)
      Number_of_groups <- length(unique(groupTable()$Group))

      updateSliderInput(session, "GlobalGroupX",
                        max = Biggest_group)
      updateSliderInput(session, "GlobalGroupY",
                        max = Number_of_groups)
    })


    # Protein dfs -------------------------------------------------------------

    #Import
    protein_df <- reactive({
      if(isTruthy(input$proteinFile)){
        x1 <- input$proteinFile
        req(input$proteinFile)
        readSamplesReport(x1$datapath, type)
      } else if(isTruthy(input$peptideFile)){

        x1 <- input$peptideFile
        req(input$peptideFile)
        readPeptideReport(x1$datapath, type)
      }

      
      
      # if(isTruthy(x1)){
      #   #Import, replace spaces with underscores
      #   x2 <- readSamplesReport(x1$datapath, type)
      #   x2
      # }
      
    })

    #Exclude samples, tidy, transform (log2, remove 0s, NormalizeRefs)
    protein_df2 <- reactive({
      req(protein_df())
      req(sampleNameTable3())

      clean_protein_df(protein_df(), sampleNameTable3())

    })

    #Filter for proteins that meet observation criteria
    protein_df3 <- reactive({
      req(protein_df2())
      req(groupTable2())
      req(input$GlobalGroupX)
      req(input$GlobalGroupY)

      filter_protein_df(protein_df2(), groupTable2(), input$GlobalGroupX, input$GlobalGroupY)

    })

    #To matrix
    protein_df4 <- reactive({
      req(protein_df3())
      protein_df_to_matrix(protein_df3())
    })

    #Normalize
    protein_df5 <- reactive({
      req(protein_df4())
      req(input$normOption)
      req(input$normCheckbox == TRUE)

      normalize_protein_df(protein_df4(), input$normOption)

    })

    #Back to tidy
    protein_df6 <- reactive({
      req(protein_df5())
      matrix_to_long_form(protein_df5())

    })



    # Plots -------------------------------------------------------------------
    output$P1 <- renderPlot({
      req(protein_df2())
      req(input$P1checkbox == TRUE)
      req(sampleNameTable2())

      make_boxplot(protein_df2(), sampleNameTable2())
    })

    output$P2 <- renderPlot({
      req(protein_df6())
      req(sampleNameTable3())
      make_boxplot2(protein_df6(), sampleNameTable3())
    })

    output$P3 <- renderPlotly({
      req(protein_df5())
      req(sampleNameTable3())

      PCA_plotly(protein_df5(), sampleNameTable3())

    })

    output$CorP3 <- renderPlot({
      req(protein_df6())
      req(sampleNameTable3())
      make_cor_plot(protein_df6(), sampleNameTable3(), scaleRows = input$CorP3check)
    })

    output$NaP3 <- renderPlot({
      req(protein_df3())
      req(sampleNameTable3())
      if(input$NaP3checkbox){
        make_NA_plot(protein_df3(), sampleNameTable3())
      }
    })


    # Extra -------------------------------------------------------------------

    #Reactive rendering of slider inputs
    output$groupSlider <- renderUI({
      req(groupTable())
      df <- groupTable()
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
      req(protein_df6())
      req(sampleNameTable3())
      save_quantitative_data(protein_df6(), input$projectName, sampleNameTable3(), "protein")
    })
    
    observeEvent(input$P1peptidedfbutton, {
      req(input$projectName)
      req(protein_df6())
      req(sampleNameTable3())
      save_quantitative_data(protein_df6(), input$projectName, sampleNameTable3(), "peptide")
    })

    observeEvent(input$saveProteinMeta, {
      req(protein_df())
      req(input$projectName)

      save_protein_meta(protein_df(), input$projectName, type)
    })
    
    observeEvent(input$savePeptideMeta, {
      req(protein_df())
      req(input$projectName)
      
      save_peptide_meta(protein_df(), input$projectName, type)
    })

    observeEvent(input$saveDesignTable, {
      req(sampleNameTable3())
      req(input$projectName)
      save_design_table(sampleNameTable3(), input$projectName)
    })

    observeEvent(input$makeContrastList, {
      req(sampleNameTable3())
      req(input$projectName)
      make_contrast_list(sampleNameTable3(), input$projectName)
    })



  }

  shinyApp(ui, server, options = options)

}
