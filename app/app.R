# app.R

# --- 1. Load Libraries ---
library(shiny)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr) # For str_replace

# --- 2. Global Data Loading and Preprocessing (ROBUST CONVERSION FOR SHINYLIVE) ---

# Set up data paths assuming the 'data' folder is next to app.R
data_dir <- "data/"

# Load data - Part 1: Intensities
d_biopsies_noexcl <- read_tsv(
    paste0(data_dir, "Biopsies_PiZ_report.tsv"), 
    show_col_types = FALSE # Let readr guess, but we will force conversion below
) %>%
  # Select the protein ID and all intensity columns
  select(Protein.Group, contains("Evo12"))

# Load data - Part 2: Sample Metadata
meta_biopsies <- read_tsv(paste0(data_dir, "meta_biopsies.txt"), show_col_types = FALSE)

# Load data - Part 3: Protein Metadata
meta_pg <- read_tsv(paste0(data_dir, "Biopsies_PiZ_report.tsv"), show_col_types = FALSE) %>%
  select(Protein.Group, Genes) 


# --- Data cleaning and filtering (CRITICAL FIX FOR TYPE CONSISTENCY)
d_biopsies_noexcl %>%
  
  # CRITICAL FIX: Explicitly force all Evo12 columns to be numeric before gathering. 
  # This resolves inconsistent type guessing in the WebAssembly environment.
  mutate(across(starts_with("Evo12"), as.numeric)) %>%

  # Gather the data
  gather(ms_id, int, contains("Evo12")) %>%
  
  # Filter out zero intensity 
  filter(int != 0) -> d_long_noexcl

# Calculate exclusion threshold (stats_lower_n)
stats <- d_long_noexcl %>%
  group_by(ms_id) %>%
  summarise(n = n(), .groups = 'drop') %>%
  ungroup() %>%
  summarise(mean = mean(n), sd = sd(n))
stats_lower_n <- stats$mean - 1.5 * stats$sd

# Determine included samples (SA_incl)
SA_incl <- d_long_noexcl %>%
  group_by(ms_id) %>%
  summarise(n = n(), .groups = 'drop') %>%
  filter(n > stats_lower_n) %>%
  mutate(well_id = str_replace(str_replace(ms_id, ".*__", ""), "_.*", "")) %>%
  filter(well_id %in% meta_biopsies$well_id) %>%
  left_join(meta_biopsies, by = "well_id") %>%
  filter(tech_rep %in% c(1, NA)) %>%
  filter(include == TRUE) %>%
  pull(ms_id)

# Filter main long data frame
d_long <- d_long_noexcl %>%
  filter(ms_id %in% SA_incl) %>%
  mutate(well_id = str_replace(str_replace(ms_id, ".*__", ""), "_.*", ""))

# Get unique list of available genes for input selectors
available_genes <- d_long %>% 
  left_join(meta_pg, by = "Protein.Group") %>% 
  distinct(Genes) %>% 
  pull(Genes) %>% 
  na.omit() %>% 
  unique() %>%
  sort()

# Pre-join and process data for efficient reactive filtering
processed_data <- d_long %>%
  left_join(meta_pg, by = "Protein.Group") %>%
  select(ms_id, Genes, int) %>%
  na.omit()


# --- 3. UI (User Interface) ---
ui <- fluidPage(
    titlePanel("Gene Correlation Explorer (Biopsies Data)"),

    sidebarLayout(
        sidebarPanel(
            selectInput("x_gene", "X-Axis Gene:", 
                        choices = available_genes, 
                        selected = "SERPINA1"),
            
            selectInput("y_gene", "Y-Axis Gene:", 
                        choices = available_genes, 
                        selected = "CAPN2"),
            
            textInput("plot_title", "Plot Title (Optional):", 
                      value = ""),
            
            downloadButton("downloadPlot", "Download Plot (.pdf)")
        ),

        mainPanel(
           plotOutput("gene_plot")
        )
    )
)

# --- 4. Server Logic ---
server <- function(input, output) {

    # Reactive function to prepare data and generate the plot
    gene_plot_data <- reactive({
        req(input$x_gene, input$y_gene) # Ensure genes are selected

        x_gene <- input$x_gene
        y_gene <- input$y_gene
        
        # 1. --- Prepare Data ---
        df_plot <- processed_data %>%
            filter(Genes %in% c(x_gene, y_gene)) %>%
            tidyr::pivot_wider(
                id_cols = ms_id, 
                names_from = Genes, 
                values_from = int, 
                values_fn = mean
            ) %>%
            na.omit() 
        
        if(nrow(df_plot) < 2) {
            return(NULL)
        }

        # 2. --- Calculate Correlation Statistics ---
        cor_result <- cor.test(log2(df_plot[[x_gene]]), log2(df_plot[[y_gene]]), method = "pearson")
        
        r_label <- paste0("R = ", format(cor_result$estimate, digits = 3))
        p_label <- paste0("P = ", format.pval(cor_result$p.value, digits = 3, eps = 0.001))
        label_text <- paste(r_label, p_label, sep = ", ")

        # 3. --- Determine Label Position (Top Left) ---
        x_log <- log2(df_plot[[x_gene]])
        y_log <- log2(df_plot[[y_gene]])
        
        x_range <- range(x_log, na.rm = TRUE)
        y_range <- range(y_log, na.rm = TRUE)
        
        x_pos <- x_range[1] + 0.05 * (x_range[2] - x_range[1])
        y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])

        label_df <- data.frame(
            x = x_pos,
            y = y_pos,
            label = label_text
        )
        
        list(
            df_plot = df_plot, 
            label_df = label_df
        )
    })

    # Reactive function to generate the actual plot
    gene_plot_reactive <- reactive({
        plot_data_list <- gene_plot_data()
        
        if (is.null(plot_data_list)) {
            return(
                ggplot() + 
                labs(title = paste("Insufficient data for:", input$y_gene, "vs", input$x_gene)) +
                theme_void()
            )
        }
        
        df_plot <- plot_data_list$df_plot
        label_df <- plot_data_list$label_df
        x_gene <- input$x_gene
        y_gene <- input$y_gene
        
        # 4. --- Plot Generation ---
        plot <- ggplot(df_plot, aes(x = log2(.data[[x_gene]]), y = log2(.data[[y_gene]]))) +
            
            geom_point() +
            
            geom_smooth(method = "loess", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.5) +
            
            geom_text(
                data = label_df, 
                aes(x = x, y = y, label = label), 
                inherit.aes = FALSE,
                hjust = 0, 
                vjust = 1,
                size = 4
            ) +
            
            theme_classic(base_size = 13) +
            labs(
                x = paste("log2 Intensity:", x_gene),
                y = paste("log2 Intensity:", y_gene),
                title = ifelse(
                    input$plot_title == "", 
                    paste(y_gene, "vs", x_gene), 
                    input$plot_title
                )
            )
        
        return(plot)
    })


    # Render the reactive plot to the UI
    output$gene_plot <- renderPlot({
        gene_plot_reactive()
    })

    # Download Handler for PDF (SYNTAX IS CORRECT FOR SHINYLIVE)
    output$downloadPlot <- downloadHandler(
        filename = function() {
            # Ensure file extension is explicitly .pdf
            paste0("correlation_plot_", input$y_gene, "_vs_", input$x_gene, ".pdf")
        },
        content = function(file) {
            # Shiny's Wasm bindings intercept this and handle the stream download
            ggsave(file, plot = gene_plot_reactive(), device = "pdf", width = 7, height = 7)
        }
    )
}

# --- 5. Run the application ---
shinyApp(ui = ui, server = server)
