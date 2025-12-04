# app.R (ShinyLive compatible)

# --- 1. Load Packages (ShinyLive Safe) ---
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 2. Load Preprocessed RDS Data ---
# (RDS loads reliably in WebR – TSV does not)
d_biopsies_noexcl <- readRDS("data/d_biopsies_noexcl.rds")
meta_biopsies     <- readRDS("data/meta_biopsies.rds")
meta_pg           <- readRDS("data/meta_pg.rds")

# --- 3. Minimal String Cleaning Without stringr ---
extract_well <- function(x) {
  x <- sub(".*__", "", x)     # remove prefix
  x <- sub("_.*$", "", x)     # remove suffix
  x
}

# --- 4. Preprocessing (WebR-safe) ---

# Convert to long format
d_long_noexcl <- d_biopsies_noexcl %>%
  pivot_longer(cols = contains("Evo12"),
               names_to = "ms_id",
               values_to = "int") %>%
  mutate(int = as.numeric(int)) %>%
  filter(int != 0)

# Compute low coverage threshold
stats <- d_long_noexcl %>%
  count(ms_id) %>%
  summarise(mean = mean(n), sd = sd(n))

stats_lower_n <- stats$mean - 1.5 * stats$sd

# Determine included samples
SA_incl <- d_long_noexcl %>%
  count(ms_id) %>%
  filter(n > stats_lower_n) %>%
  mutate(well_id = extract_well(ms_id)) %>%
  inner_join(meta_biopsies, by = "well_id") %>%
  filter(tech_rep %in% c(1, NA)) %>%
  filter(include == TRUE) %>%
  pull(ms_id)

# Final filtered long dataset
d_long <- d_long_noexcl %>%
  filter(ms_id %in% SA_incl) %>%
  mutate(well_id = extract_well(ms_id))

# All available genes
available_genes <- d_long %>%
  left_join(meta_pg, by = "Protein.Group") %>%
  filter(!is.na(Genes)) %>%
  distinct(Genes) %>%
  arrange(Genes) %>%
  pull(Genes)

# Processed dataset for plotting
processed_data <- d_long %>%
  left_join(meta_pg, by = "Protein.Group") %>%
  select(ms_id, Genes, int) %>%
  filter(!is.na(Genes))


# --- 5. UI ---
ui <- fluidPage(
  titlePanel("Gene correlation explorer (ShinyLive Version)"),

  sidebarLayout(
    sidebarPanel(
      selectInput("x_gene", "X-Axis Gene:", choices = available_genes, selected = "SERPINA1"),
      selectInput("y_gene", "Y-Axis Gene:", choices = available_genes, selected = "CAPN2"),
      textInput("plot_title", "Plot Title (Optional):", value = ""),
      downloadButton("downloadPlot", "Download Plot (.png)")
    ),
    mainPanel(
      plotOutput("gene_plot")
    )
  )
)

# --- 6. Server ---
server <- function(input, output) {

  # Prepare data for plot (WebR safe: NO pivot_wider(values_fn))
  gene_plot_data <- reactive({
    req(input$x_gene, input$y_gene)

    df <- processed_data %>%
      filter(Genes %in% c(input$x_gene, input$y_gene)) %>%
      group_by(ms_id, Genes) %>%
      summarise(int = mean(int), .groups = "drop") %>%
      pivot_wider(names_from = Genes, values_from = int) %>%
      na.omit()

    if (nrow(df) < 2) return(NULL)

    # Compute correlation
    x <- log2(df[[input$x_gene]])
    y <- log2(df[[input$y_gene]])

    cor_result <- cor.test(x, y, method = "pearson")

    label <- sprintf("R = %.3f, P = %.3g", 
                     cor_result$estimate, cor_result$p.value)

    # Label position
    x_range <- range(x, na.rm = TRUE)
    y_range <- range(y, na.rm = TRUE)

    label_df <- data.frame(
      x = x_range[1] + 0.05 * diff(x_range),
      y = y_range[2] - 0.05 * diff(y_range),
      label = label
    )

    list(df = df, label_df = label_df)
  })

  # Generate Plot
  output$gene_plot <- renderPlot({
    dat <- gene_plot_data()
    if (is.null(dat)) {
      return(ggplot() + 
               labs(title = "Not enough data to plot") +
               theme_void())
    }

    xg <- input$x_gene
    yg <- input$y_gene

    ggplot(dat$df, aes(x = log2(.data[[xg]]), y = log2(.data[[yg]]))) +
      geom_point() +
      geom_smooth(method = "loess", se = TRUE) +
      geom_text(data = dat$label_df, aes(x = x, y = y, label = label),
                hjust = 0, vjust = 1, size = 4) +
      theme_classic(base_size = 13) +
      labs(
        x = paste("log2 Intensity:", xg),
        y = paste("log2 Intensity:", yg),
        title = ifelse(input$plot_title == "", 
                       paste(yg, "vs", xg), 
                       input$plot_title)
      )
  })



gene_plot_reactive <- reactive({
  dat <- gene_plot_data()
  if (is.null(dat)) {
    return(ggplot() + labs(title = "Not enough data to plot") + theme_void())
  }

  xg <- input$x_gene
  yg <- input$y_gene

  ggplot(dat$df, aes(x = log2(.data[[xg]]), y = log2(.data[[yg]]))) +
    geom_point() +
    geom_smooth(method = "loess", se = TRUE) +
    geom_text(data = dat$label_df, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 4) +
    theme_classic(base_size = 13) +
    labs(
      x = paste("log2 Intensity:", xg),
      y = paste("log2 Intensity:", yg),
      title = ifelse(input$plot_title == "", paste(yg, "vs", xg), input$plot_title)
    )
})

  # Download PNG
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("correlation_plot_", input$y_gene, "_vs_", input$x_gene, ".png")
    },
    content = function(file) {
      png(file, width = 7, height = 7, units = "in", res = 600)
      print(gene_plot_reactive())
      dev.off()
    }
  )

}

# Run app
shinyApp(ui, server)
