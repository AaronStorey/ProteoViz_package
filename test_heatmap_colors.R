# Minimal test for heatmaply color scale with scale_rows = TRUE
# Run with: shiny::runApp("test_heatmap_colors.R")
# or source() it — opens in the Viewer/browser

library(shiny)
library(heatmaply)
library(ggplot2)
library(viridisLite)

# ---- synthetic data --------------------------------------------------------
set.seed(42)
mat <- matrix(rnorm(80), nrow = 10, ncol = 8,
              dimnames = list(
                paste0("Protein", 1:10),
                paste0("Sample", 1:8)
              ))

col_anno <- data.frame(
  Group = rep(c("A", "B"), each = 4),
  row.names = colnames(mat)
)

bwr_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(256)

# ---- four approaches to test -----------------------------------------------
make_plot <- function(approach) {
  switch(approach,

    "1_scale_fun_colors_NULL" = heatmaply(
      mat,
      scale                   = "row",
      col_side_colors         = col_anno,
      plot_method             = "plotly",
      colors                  = NULL,
      scale_fill_gradient_fun = scale_fill_gradient2(
        low = "blue", mid = "white", high = "red"
      )
    ),

    "2_scale_fun_no_colors_arg" = heatmaply(
      mat,
      scale                   = "row",
      col_side_colors         = col_anno,
      plot_method             = "plotly",
      scale_fill_gradient_fun = scale_fill_gradient2(
        low = "blue", mid = "white", high = "red"
      )
    ),

    "3_colors_bwr_palette" = heatmaply(
      mat,
      scale           = "row",
      col_side_colors = col_anno,
      plot_method     = "plotly",
      colors          = bwr_palette
    ),

    "4_ggplot_method_scale_fun" = heatmaply(
      mat,
      scale                   = "row",
      col_side_colors         = col_anno,
      plot_method             = "ggplot",
      scale_fill_gradient_fun = scale_fill_gradient2(
        low = "blue", mid = "white", high = "red"
      )
    )
  )
}

# ---- app -------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("heatmaply color scale test — scale_rows = TRUE"),
  selectInput("approach", "Approach:",
              choices = c(
                "1: scale_fill_gradient_fun + colors=NULL"   = "1_scale_fun_colors_NULL",
                "2: scale_fill_gradient_fun (no colors arg)" = "2_scale_fun_no_colors_arg",
                "3: colors = bwr colorRampPalette"           = "3_colors_bwr_palette",
                "4: plot_method='ggplot' + scale_fill_gradient_fun" = "4_ggplot_method_scale_fun"
              )),
  plotly::plotlyOutput("heatmap", height = "600px")
)

server <- function(input, output, session) {
  output$heatmap <- plotly::renderPlotly({
    make_plot(input$approach)
  })
}

shinyApp(ui, server)
