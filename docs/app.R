library(shiny)

ui <- fluidPage(
  titlePanel("Scatterplot Example - Base R & ShinyLive Compatible"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "Number of points:", 10, 1000, 200),
      downloadButton("downloadData", "Download Data (CSV)")
    ),
    
    mainPanel(
      plotOutput("scatterPlot", height = "500px")
    )
  )
)

server <- function(input, output, session) {
  
  # Generate reactive dataset
  data_reactive <- reactive({
    set.seed(123)
    x <- rnorm(input$n)
    y <- rnorm(input$n)
    data.frame(x = x, y = y)
  })
  
  # Base R plot
  output$scatterPlot <- renderPlot({
    d <- data_reactive()
    plot(
      d$x, d$y,
      pch = 19,
      main = "Base R Scatterplot (ShinyLive Compatible)",
      xlab = "X",
      ylab = "Y"
    )
  })
  
  # Download CSV (works in ShinyLive)
  output$downloadData <- downloadHandler(
    filename = function() { "scatter_data.csv" },
    content = function(file) {
      write.csv(data_reactive(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
