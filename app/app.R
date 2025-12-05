library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Load preprocessed data (RDS only) ---
d_biopsies_noexcl <- readRDS("data/d_biopsies_noexcl.rds")
meta_biopsies     <- readRDS("data/meta_biopsies.rds")
meta_pg           <- readRDS("data/meta_pg.rds")

# --- Helper: extract well_id ---
extract_well <- function(x) {
  x <- sub(".*__", "", x)
  sub("_.*$", "", x)
}

# --- Preprocessing ---
d_long_noexcl <- d_biopsies_noexcl %>%
  pivot_longer(cols = contains("Evo12"),
               names_to = "ms_id",
               values_to = "int") %>%
  mutate(int = as.numeric(int)) %>%
  filter(int != 0)

stats_lower_n <- d_long_noexcl %>%
  count(ms_id) %>%
  summarise(mean = mean(n), sd = sd(n)) %>%
  { .$mean - 1.5 * .$sd }

SA_incl <- d_long_noexcl %>%
  count(ms_id) %>%
  filter(n > stats_lower_n) %>%
  mutate(well_id = extract_well(ms_id)) %>%
  inner_join(meta_biopsies, by = "well_id") %>%
  filter(tech_rep %in% c(1, NA), include == TRUE) %>%
  pull(ms_id)

d_long <- d_long_noexcl %>%
  filter(ms_id %in% SA_incl) %>%
  mutate(well_id = extract_well(ms_id))

available_genes <- d_long %>%
  left_join(meta_pg, by = "Protein.Group") %>%
  filter(!is.na(Genes)) %>%
  distinct(Genes) %>%
  arrange(Genes) %>%
  pull(Genes)

processed_data <- d_long %>%
  left_join(meta_pg, by = "Protein.Group") %>%
  select(ms_id, Genes, int) %>%
  filter(!is.na(Genes))

# --- UI ---
ui <- fluidPage(
  titlePanel("Gene Correlation Explorer (Shinylive)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("x_gene", "X-Axis Gene:", choices = available_genes, selected = available_genes[1]),
      selectInput("y_gene", "Y-Axis Gene:", choices = available_genes, selected = available_genes[2]),
      textInput("plot_title", "Plot Title (Optional):", value = "")
    ),
    mainPanel(
      plotOutput("gene_plot")
    )
  )
)

# --- Server ---
server <- function(input, output) {
  gene_plot_data <- reactive({
    df <- processed_data %>%
      filter(Genes %in% c(input$x_gene, input$y_gene)) %>%
      group_by(ms_id, Genes) %>%
      summarise(int = mean(int), .groups = "drop") %>%
      pivot_wider(names_from = Genes, values_from = int) %>%
      na.omit()
    if (nrow(df) < 2) return(NULL)

    cor_result <- cor.test(log2(df[[input$x_gene]]), log2(df[[input$y_gene]]))
    label_df <- data.frame(
      x = min(log2(df[[input$x_gene]])) + 0.05 * diff(range(log2(df[[input$x_gene]]))),
      y = max(log2(df[[input$y_gene]])) - 0.05 * diff(range(log2(df[[input$y_gene]]))),
      label = sprintf("R = %.3f, P = %.3g", cor_result$estimate, cor_result$p.value)
    )
    list(df = df, label_df = label_df)
  })

  output$gene_plot <- renderPlot({
    dat <- gene_plot_data()
    if (is.null(dat)) {
      return(ggplot() + labs(title = "Not enough data") + theme_void())
    }
    ggplot(dat$df, aes(x = log2(.data[[input$x_gene]]), y = log2(.data[[input$y_gene]]))) +
      geom_point() +
      geom_smooth(method = "loess", se = TRUE) +
      geom_text(data = dat$label_df, aes(x = x, y = y, label = label),
                hjust = 0, vjust = 1, size = 4) +
      theme_classic(base_size = 13) +
      labs(
        x = paste("log2 Intensity:", input$x_gene),
        y = paste("log2 Intensity:", input$y_gene),
        title = ifelse(input$plot_title == "", 
                       paste(input$y_gene, "vs", input$x_gene), 
                       input$plot_title)
      )
  })
}

shinyApp(ui, server)
