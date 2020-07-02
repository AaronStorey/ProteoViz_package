#' Runs the Shiny TMT App
#' @import tidyverse
#' @import shinydashboard
#' @import rhandsontable
#' @import limma
#' @import shiny
#' @import cowplot
#' @import plotly
#' @import heatmaply
#' @import UpSetR
#' @export
#'

runApp2 <- function(){
  theme_set(theme_cowplot())
  options(shiny.maxRequestSize = 5000*1024^2)



  # Limma tab --------------------------------------------------------------




  Box1 <- {
    fluidRow(
      box(
        width = 12,
        fileInput("quantData", "Import quantitative data"),
        fileInput("metadata", "Import metadata")
      )
    )
  }

  Box2 <- {
    fluidRow(
      box(
        width = 5,
        title = "Sample grouping",
        fileInput("sampleFile", "Import sample file"),
        rHandsontableOutput("sampleNameTable")
      ),
      box(
        width = 5,
        title = "Design matrix",
        checkboxInput("blockCheck", "Blocked design"),
        rHandsontableOutput("T2")
      ),
      box(
        width = 2,
        title = "Contrast matrix",
        fileInput("contrastFile", "Import contrast list"),
        rHandsontableOutput("T3")
      )

    )
  }

  Box3 <- {
    fluidRow(
      box(
        width = 6,
        title = "T4",
        tableOutput("T4")
      ),
      box(
        width = 6,
        title = "T5",
        tableOutput("T5")
      )
    )
  }

  Box4 <- {
    fluidRow(
      box(
        width = 6,
        title = "T6",
        tableOutput("T6")
      ),
      box(
        width = 6,
        title = "T7",
        tableOutput("T7")
      )
    )
  }

  Box5 <- {
    fluidRow(
      box(
        width = 12,
        title = "P7",
        plotOutput("P7")
      )
    )
  }

  BoxP8 <- {
    fluidRow(
      box(
        width = 12,
        title = "P8 Venn",
        plotOutput("P8",
                   height = "900px")
      )
    )
  }

  BoxP9 <- {
    fluidRow(
      box(
        width = 12,
        title = "P9 Set Selection",
        selectInput("P9select", "Select sets", choices = NULL, selected = NULL,
                    multiple = TRUE, width = "50%"),
        tableOutput("P9"),
        actionButton("saveSets", "Extract and save all sets")
      )
    )
  }

  Box6 <- {
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
            Box1,
            Box2,
            # Box3,
            # Box4,
            Box5,
            BoxP8,
            BoxP9,
            Box6
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
                 actionButton("pv1button", "Update volcano plot"))

        ),
        fluidRow(
          column(6,
                 plotlyOutput("ProteinVolcanoplot"),
                 selectInput("ProteinSearch", label = "Highlight a protein:",
                             choices = character(0), multiple = TRUE)
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
          column(4, selectInput("proteinHeatGroupExclude", label = "Exclude groups:",
                                choices = character(0), multiple = TRUE))
        ),
        fluidRow(
          style = "height:900px",
          tabBox(title = "Selected proteins", id = "tab1", width = 12,
                 tabPanel("Select proteins",
                          plotlyOutput("proteinHeat", height = "900px")),
                 tabPanel("Summarized",
                          plotlyOutput("ProteinHeatAverage", height = "900px"))
          )
        )
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
      req(input$quantData)
      df <- input$quantData$datapath
      read_tsv(df)
    })

    limma_input <- reactive({
      req(quant_data())
      quant_data() %>%
        as.data.frame() %>%
        column_to_rownames("id")
    })



    # Metadata ----------------------------------------------------------------

    metadata <- reactive({
      req(input$metadata)
      if(isTruthy(input$metadata)){
        df <- input$metadata
        read_tsv(df$datapath)
      }
    })




    # Sample table ------------------------------------------------------------


    sampleNameTable <- reactive({

      if (isTruthy(input$sampleFile)) {

        col_check <- c("Sample_name", "Group", "Batch", "Block")

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
    output$sampleNameTable <- renderRHandsontable({
      req(sampleNameTable())
      sampleNameTable()
    })

    # Changes the handsontable back into a dataframe
    sampleNameTable2 <- reactive({
      if(is.null(input$sampleNameTable)) return (NULL)
      df <- hot_to_r(input$sampleNameTable)
    })




    # Design table ------------------------------------------------------------
    design_table <- reactive({
      req(sampleNameTable2())
      makeDesignTable(sampleNameTable2())
    })

    sample_order <- reactive({
      req(sampleNameTable2())
      sampleNameTable2()$Sample_name
    })

    group_order <- reactive({
      req(sampleNameTable2())
      sampleNameTable2()$Group

    })


    #Not sure whether to make the design matrix editable or not.
    design_table2 <- reactive({
      if(is.null(input$T2)) return (NULL)
      df <- hot_to_r(input$T2)
      df
    })



    # Contrast table ----------------------------------------------------------

    contrast_table <- reactive({
      if(is.null(input$T3)) return (NULL)
      df <- hot_to_r(input$T3) %>%
        #Removes blank entries
        filter(Contrast_name != "")
      df
    })


    contrast_list <- reactive({
      req(design_table2())
      req(contrast_table())

      makeContrastTable(design_table2(), contrast_table())
    })





    # Limma -------------------------------------------------------------------
    limma_output <- reactive({
      req(limma_input())
      req(design_table2())
      req(sampleNameTable2())

      runLimma(limma_input(), design_table2(), sampleNameTable2(), block = input$blockCheck)

    })

    contrast_fit <- reactive({
      req(limma_output())
      req(contrast_list())

      fitContrasts(limma_output(), contrast_list())
    })

    combined_data <- reactive({

      req(contrast_fit())
      req(design_table2())
      req(contrast_table())

      combineLimmaData(contrast_fit(), design_table(), contrast_table())
    })

    test_results <- reactive({
      req(contrast_fit())

      decideTests(contrast_fit(), lfc = log2(1.5))

    })

    tidy_test_results <- reactive({
      req(test_results())

      makeTidyTestResults(test_results())
    })

    test_results_list <- reactive({
      req(tidy_test_results())
      makeTidyTestResultsList(tidy_test_results())

    })


    # Table plots -------------------------------------------------------------
    output$T1 <- renderTable({
      req(quant_data())
      quant_data() %>%
        dplyr::slice(1:10)
    })

    output$T2 <- renderRHandsontable({
      req(design_table())

      sample_names <- rownames(design_table())

      df3 <- apply(design_table(), 2, as.integer) %>%
        as.data.frame()

      rownames(df3) <- sample_names

      rhandsontable(df3)

    })


    output$T3 <- renderRHandsontable({

      if(isTruthy(input$contrastFile)){
        #File import
        df <- input$contrastFile$datapath
        read_tsv(df) %>%
          as.data.frame() %>%
          rhandsontable()
      } else {
        #Manual entry
        req(sampleNameTable2())
        tibble(Contrast_name = "",
        ) %>%
          as.data.frame() %>%
          rhandsontable()
      }

    })



    # Save output -------------------------------------------------------------
    observeEvent(input$Save_output, {
      req(combined_data())
      combined_data() %>%
        write_tsv("Limma_output.tsv")
    })

    observeEvent(input$Save_summary, {
      req(quant_data())
      req(combined_data())
      req(metadata())

      spread_protein_limma <- combined_data() %>%
        dplyr::select(id, logFC, P.Value, adj.P.Val, Comparison) %>%
        gather(Type, Value, logFC, P.Value, adj.P.Val) %>%
        unite(Type, Comparison, Type, sep = " ") %>%
        spread(Type, Value) %>%
        mutate(id = as.integer(id))

      meta_df <- metadata() %>%
        mutate(id = as.integer(id))

      meta_df %>%
        right_join(quant_data()) %>%
        left_join(spread_protein_limma) %>%
        write_tsv("Summarized_output.tsv")
    })

    observeEvent(input$saveSets, {
      req(tidy_test_results())
      req(metadata())

      meta_df <- metadata() %>%
        mutate(id = as.character(id))

      unnest(tidy_test_results(), data) %>%
        group_by(id) %>%
        summarize(Sets = str_c(Test, collapse = ", "),
                  Number_of_groups = length(Test)) %>%
        dplyr::select(Number_of_groups, Sets, id) %>%
        ungroup() %>%
        left_join(meta_df, by = "id") %>%
        arrange(desc(Number_of_groups)) %>%
        write_tsv("Extracted_sets.tsv")
    })





    # Debugging ---------------------------------------------------------------

    output$T4 <- renderTable({
      req(contrast_list())
      contrast_list()
    }, rownames = TRUE)

    output$T5 <- renderTable({
      req(limma_output())

      head(limma_output()$coefficients)
    })

    output$T6 <- renderTable({
      req(combined_data())

      head(combined_data())
    })

    output$T7 <- renderTable({
      req(contrast_fit())
      req(contrast_table())
      x1 <- contrast_table()$Contrast_name %>%
        .[[1]]
      head(topTable(contrast_fit(), x1, number = Inf))

    })

    output$P7 <- renderPlot({
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

    output$P8 <- renderPlot({

      req(tidy_test_results())

      a2 <- map(tidy_test_results()$data, "id")
      names(a2) <- tidy_test_results()$Test


      upset(fromList(a2), nsets = length(a2),
            text.scale = 1.5,
            point.size = 4)


    })

    observe({
      req(tidy_test_results())

      df <- tidy_test_results()

      x2 <- unique(df$Test)

      updateSelectInput(session, "P9select",
                        choices = x2)
    })

    output$P9 <- renderTable({
      req(input$P9select)
      req(tidy_test_results())
      req(metadata())

      queries <- input$P9select

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
        group_by(id) %>%
        mutate(Number_of_groups = length(Test),
               Matches_query = Test %in% queries) %>%
        filter(all(Matches_query),
               Number_of_groups == length(queries)) %>%
        ungroup() %>%
        distinct(id)

      meta_df <- metadata() %>%
        mutate(id = as.character(id))

      a5 <- a4 %>%
        mutate(id = as.character(id)) %>%
        left_join(meta_df)

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
      input$proteinHeatGroupExclude
    }) %>%
      debounce(750)

    output$ProteinVolcanoplot <- renderPlotly({
      req(combined_data)
      req(input$ProteinVolcanoComparison)
      req(metadata())

      input$pv1button
      volcanoComparison <- isolate(input$ProteinVolcanoComparison)
      FCcutoff <- isolate(as.numeric(input$pv1FC))
      Pcutoff <- isolate(as.numeric(input$pv1pvalue))

      protein_df <- metadata()


      p1 <- combined_data() %>%
        filter(Comparison %in% volcanoComparison,
               abs(logFC) > FCcutoff,
               adj.P.Val < Pcutoff) %>%
        left_join(protein_df, by = "id")

      p2 <- p1 %>%
        ggplot(aes(x = logFC, y = -log10(adj.P.Val), key = id,
                   text = paste0(
                     "Protein: ",
                     Description,
                     "\n",
                     "Gene name: ",
                     Gene_name,
                     "\n",
                     "Uniprot ID: ",
                     Uniprot_ID,
                     "\n",
                     "Log2 fold change: ",
                     signif(logFC, 4),
                     "\n",
                     "p.value: ",
                     signif(adj.P.Val, 4),
                     "\n",
                     "id: ",
                     id,
                     "\n"
                   )
        )) +
        geom_point() +
        geom_hline(yintercept = -log10(.05), linetype = "dotted") +
        geom_vline(xintercept = log2(1.5), linetype = "dotted") +
        geom_vline(xintercept = -log2(1.5), linetype = "dotted") +
        labs(x = "log2 FC", y= "-log10 adj.p.value")

      if(isTruthy(input$ProteinSearch)){
        p_search <- protein_df %>%
          filter(Protein_name2 %in% input$ProteinSearch,
          ) %>%
          left_join(protein_limma, by = "id") %>%
          filter(Comparison %in% volcanoComparison)

        if(nrow(p_search) > 0){
          p2 <- p2 +
            geom_point(data = p_search, color = "red", size = 5)
        }

      }

      ggplotly(p2, source = "proteinVolcano", tooltip = "text") %>%
        layout(dragmode = "select")
    })

    selected_data <- reactive({
      req(quant_data())
      d <- event_data("plotly_selected", source = "proteinVolcano")
      req(!is.null(d))

      quant_data() %>%
        filter(id %in% d$key)
    })

    output$proteinHeat <- renderPlotly({
      req(input$phcheck)
      req(quant_data())
      req(metadata())
      req(sample_order())
      req(group_order())
      req(selected_data())

      makeProteinHeatHeatmap(selected_data(), metadata(), sample_order(), group_order(), excluded_groups())
    })



    # -------------------------------------------------------------------------


  }

  shinyApp(ui, server)
}
