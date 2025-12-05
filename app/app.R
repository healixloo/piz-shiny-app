# app.R (ShinyLive compatible - Base R only)

# --- 1. Load Packages (Base R + Shiny only) ---
library(shiny)

# --- 2. Load Preprocessed RDS Data ---
d_biopsies_noexcl <- readRDS("data/d_biopsies_noexcl.rds")
meta_biopsies     <- readRDS("data/meta_biopsies.rds")
meta_pg           <- readRDS("data/meta_pg.rds")

# --- 3. String Cleaning (Base R) ---
extract_well <- function(x) {
  x <- sub(".*__", "", x)
  x <- sub("_.*$", "", x)
  x
}

# --- 4. Preprocessing (Base R only) ---

# Convert to long format (base R)
evo12_cols <- grep("Evo12", names(d_biopsies_noexcl), value = TRUE)

d_long_list <- lapply(evo12_cols, function(col) {
  data.frame(
    Protein.Group = d_biopsies_noexcl$Protein.Group,
    ms_id = col,
    int = as.numeric(d_biopsies_noexcl[[col]]),
    stringsAsFactors = FALSE
  )
})

d_long_noexcl <- do.call(rbind, d_long_list)
d_long_noexcl <- d_long_noexcl[d_long_noexcl$int != 0 & !is.na(d_long_noexcl$int), ]

# Compute low coverage threshold
n_per_sample <- tapply(d_long_noexcl$int, d_long_noexcl$ms_id, length)
stats_mean <- mean(n_per_sample)
stats_sd <- sd(n_per_sample)
stats_lower_n <- stats_mean - 1.5 * stats_sd

# Determine included samples
samples_pass <- names(n_per_sample)[n_per_sample > stats_lower_n]
well_ids <- extract_well(samples_pass)

# Filter meta_biopsies
meta_filtered <- meta_biopsies[
  meta_biopsies$well_id %in% well_ids &
    (meta_biopsies$tech_rep %in% c(1, NA) | is.na(meta_biopsies$tech_rep)) &
    meta_biopsies$include == TRUE,
]

SA_incl <- samples_pass[well_ids %in% meta_filtered$well_id]

# Final filtered dataset
d_long <- d_long_noexcl[d_long_noexcl$ms_id %in% SA_incl, ]
d_long$well_id <- extract_well(d_long$ms_id)

# Merge with meta_pg to get Genes
processed_data <- merge(d_long, meta_pg[, c("Protein.Group", "Genes")], 
                        by = "Protein.Group", all.x = TRUE)
processed_data <- processed_data[!is.na(processed_data$Genes), ]

# Get available genes
available_genes <- sort(unique(processed_data$Genes))

# --- 5. UI ---
ui <- fluidPage(
  titlePanel("Gene Correlation Explorer (Base R Version)"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("x_gene", "X-Axis Gene:", 
                  choices = available_genes, 
                  selected = if("SERPINA1" %in% available_genes) "SERPINA1" else available_genes[1]),
      selectInput("y_gene", "Y-Axis Gene:", 
                  choices = available_genes, 
                  selected = if("CAPN2" %in% available_genes) "CAPN2" else available_genes[2]),
      textInput("plot_title", "Plot Title (Optional):", value = ""),
      
      hr(),
      
      h4("Download Options"),
      radioButtons("download_format", "Format:",
                   choices = c("PNG (High Res)" = "png",
                               "SVG (Vector)" = "svg"),
                   selected = "png"),
      downloadButton("downloadPlot", "Download Plot"),
      
      hr(),
      
      verbatimTextOutput("correlation_stats")
    ),
    
    mainPanel(
      plotOutput("gene_plot", height = "600px")
    )
  )
)

# --- 6. Server ---
server <- function(input, output) {
  
  # Prepare data for plot (Base R)
  gene_plot_data <- reactive({
    req(input$x_gene, input$y_gene)
    
    # Filter for selected genes
    df <- processed_data[processed_data$Genes %in% c(input$x_gene, input$y_gene), ]
    
    # Aggregate by ms_id and Genes (mean intensity)
    df_agg <- aggregate(int ~ ms_id + Genes, data = df, FUN = mean)
    
    # Convert to wide format (base R)
    x_data <- df_agg[df_agg$Genes == input$x_gene, ]
    y_data <- df_agg[df_agg$Genes == input$y_gene, ]
    
    df_wide <- merge(x_data[, c("ms_id", "int")], 
                     y_data[, c("ms_id", "int")], 
                     by = "ms_id", suffixes = c("_x", "_y"))
    
    # Remove NAs
    df_wide <- df_wide[complete.cases(df_wide), ]
    
    if (nrow(df_wide) < 2) return(NULL)
    
    # Compute correlation
    x <- log2(df_wide$int_x)
    y <- log2(df_wide$int_y)
    
    cor_result <- cor.test(x, y, method = "pearson")
    
    list(
      x = x,
      y = y,
      cor = cor_result$estimate,
      pval = cor_result$p.value,
      n = nrow(df_wide)
    )
  })
  
  # Generate Plot (Base R graphics)
  output$gene_plot <- renderPlot({
    dat <- gene_plot_data()
    
    if (is.null(dat)) {
      plot.new()
      text(0.5, 0.5, "Not enough data to plot", cex = 1.5)
      return()
    }
    
    x <- dat$x
    y <- dat$y
    
    # Set up plot
    plot(x, y,
         xlab = paste("log2 Intensity:", input$x_gene),
         ylab = paste("log2 Intensity:", input$y_gene),
         main = if(input$plot_title == "") {
           paste(input$y_gene, "vs", input$x_gene)
         } else {
           input$plot_title
         },
         pch = 16,
         col = rgb(0, 0, 0, 0.5),
         cex = 1.2,
         cex.lab = 1.2,
         cex.main = 1.3,
         font.main = 2
    )
    
    # Add loess smooth line
    loess_fit <- loess(y ~ x)
    x_pred <- seq(min(x), max(x), length.out = 100)
    y_pred <- predict(loess_fit, newdata = data.frame(x = x_pred))
    lines(x_pred, y_pred, col = "blue", lwd = 2)
    
    # Add confidence band (approximate)
    pred_obj <- predict(loess_fit, newdata = data.frame(x = x_pred), se = TRUE)
    lines(x_pred, pred_obj$fit + 1.96 * pred_obj$se.fit, 
          col = "lightblue", lty = 2, lwd = 1.5)
    lines(x_pred, pred_obj$fit - 1.96 * pred_obj$se.fit, 
          col = "lightblue", lty = 2, lwd = 1.5)
    
    # Add correlation text
    cor_text <- sprintf("R = %.3f\nP = %.3g\nN = %d", 
                        dat$cor, dat$pval, dat$n)
    
    # Position text in upper left
    x_range <- par("usr")[1:2]
    y_range <- par("usr")[3:4]
    text_x <- x_range[1] + 0.05 * diff(x_range)
    text_y <- y_range[2] - 0.05 * diff(y_range)
    
    text(text_x, text_y, cor_text, adj = c(0, 1), cex = 1.2, font = 2)
    
    # Add grid
    grid(col = "gray80", lty = 3)
  })
  
  # Display correlation statistics
  output$correlation_stats <- renderText({
    dat <- gene_plot_data()
    
    if (is.null(dat)) {
      return("Select genes to see statistics")
    }
    
    paste(
      "Correlation Statistics:",
      "----------------------",
      sprintf("R: %.4f", dat$cor),
      sprintf("R²: %.4f", dat$cor^2),
      sprintf("P-value: %.4e", dat$pval),
      sprintf("N samples: %d", dat$n),
      sep = "\n"
    )
  })
  
  # Create plot function for download
  create_download_plot <- function() {
    dat <- gene_plot_data()
    if (is.null(dat)) return(NULL)
    
    function() {
      x <- dat$x
      y <- dat$y
      
      plot(x, y,
           xlab = paste("log2 Intensity:", input$x_gene),
           ylab = paste("log2 Intensity:", input$y_gene),
           main = if(input$plot_title == "") {
             paste(input$y_gene, "vs", input$x_gene)
           } else {
             input$plot_title
           },
           pch = 16,
           col = rgb(0, 0, 0, 0.5),
           cex = 1.5,
           cex.lab = 1.3,
           cex.main = 1.4,
           font.main = 2
      )
      
      loess_fit <- loess(y ~ x)
      x_pred <- seq(min(x), max(x), length.out = 100)
      y_pred <- predict(loess_fit, newdata = data.frame(x = x_pred))
      lines(x_pred, y_pred, col = "blue", lwd = 3)
      
      pred_obj <- predict(loess_fit, newdata = data.frame(x = x_pred), se = TRUE)
      lines(x_pred, pred_obj$fit + 1.96 * pred_obj$se.fit, 
            col = "lightblue", lty = 2, lwd = 2)
      lines(x_pred, pred_obj$fit - 1.96 * pred_obj$se.fit, 
            col = "lightblue", lty = 2, lwd = 2)
      
      cor_text <- sprintf("R = %.3f\nP = %.3g\nN = %d", 
                          dat$cor, dat$pval, dat$n)
      
      x_range <- par("usr")[1:2]
      y_range <- par("usr")[3:4]
      text_x <- x_range[1] + 0.05 * diff(x_range)
      text_y <- y_range[2] - 0.05 * diff(y_range)
      
      text(text_x, text_y, cor_text, adj = c(0, 1), cex = 1.3, font = 2)
      grid(col = "gray80", lty = 3)
    }
  }
  
  # Download handler
  output$downloadPlot <- downloadHandler(
    filename = function() {
      ext <- if(input$download_format == "png") ".png" else ".svg"
      paste0(input$y_gene, "_vs_", input$x_gene, "_", Sys.Date(), ext)
    },
    content = function(file) {
      plot_func <- create_download_plot()
      if (is.null(plot_func)) return()
      
      if (input$download_format == "png") {
        png(file, width = 800, height = 800, res = 150)
      } else {
        svg(file, width = 10, height = 10)
      }
      
      plot_func()
      dev.off()
    }
  )
}

# Run app
shinyApp(ui, server)