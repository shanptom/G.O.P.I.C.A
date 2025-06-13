library(shiny)
library(phyloseq)
library(ggplot2)
library(vegan)
library(ranacapa)
library(phylosmith)

ui <- navbarPage("Microbiome Data Explorer",
                 tabPanel("Upload Data",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("asv", "Upload ASV Table (CSV)", accept = ".csv"),
                              fileInput("tax", "Upload Taxonomy Table (CSV)", accept = ".csv"),
                              fileInput("meta", "Upload Metadata Table (CSV)", accept = ".csv"),
                              fileInput("phylo", "Or Upload Phyloseq Object (.rds)", accept = ".rds"),
                              checkboxInput("doRarefy", "Apply rarefaction", value = FALSE)
                            ),
                            mainPanel(
                              verbatimTextOutput("upload_status")
                            )
                          )
                 ),
                 tabPanel("Rarefaction Plot",
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("rarefaction_color_selector"),
                              uiOutput("rarefaction_facet_selector")
                            ),
                            mainPanel(
                              plotOutput("rarefactionPlot")
                            )
                          )
                 ),
                 
                 tabPanel("Alpha Diversity",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("alpha_index", "Select Diversity Index:",
                                          choices = c("shannon", "simpson", "invsimpson"), selected = "shannon"),
                              uiOutput("alpha_group_selector"),
                              textInput("alpha_order", "Custom order (comma-separated values)", value = ""),
                            ),
                            mainPanel(
                              plotOutput("alphaPlot")
                            )
                          )
                 ),
                 
                 tabPanel("Beta Diversity", plotOutput("betaPlot")),
                 tabPanel("Abundance",
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("tax_selector")
                            ),
                            mainPanel(
                              plotOutput("abundancePlot")
                            )
                          )
                 ),
                 tabPanel("Ordination", plotOutput("ordinationPlot"))
)

server <- function(input, output) {
  data <- reactive({
    if (!is.null(input$phylo)) {
      ps <- readRDS(input$phylo$datapath)
    } else {
      req(input$asv, input$tax, input$meta)
      otu <- read.csv(input$asv$datapath, row.names = 1)
      tax <- as.matrix(read.csv(input$tax$datapath, row.names = 1))
      meta <- read.csv(input$meta$datapath, row.names = 1)
      
      ps <- phyloseq(
        otu_table(as.matrix(otu), taxa_are_rows = TRUE), 
        tax_table(tax), 
        sample_data(meta)
      )
    }
    
    # Optional rarefaction
    if (input$doRarefy) {
      ps <- rarefy_even_depth(
        ps, 
        sample.size = min(sample_sums(ps)), 
        rngseed = 123, 
        replace = TRUE, 
        trimOTUs = TRUE, 
        verbose = TRUE
      )
    }
    
    return(ps)
  })
  
  output$rarefaction_color_selector <- renderUI({
    ps <- data()
    meta_cols <- colnames(sample_data(ps))
    selectInput("rare_color", "Color by (metadata column):", choices = meta_cols, selected = meta_cols[1])
  })
  
  output$rarefaction_facet_selector <- renderUI({
    ps <- data()
    meta_cols <- colnames(sample_data(ps))
    selectInput("rare_facet", "Facet by (metadata column):", 
                choices = c("None", meta_cols), selected = "None")
  })
  
  output$alpha_group_selector <- renderUI({
    ps <- data()
    meta_cols <- colnames(sample_data(ps))
    selectInput("alpha_group", "Group by (metadata column):",
                choices = meta_cols, selected = meta_cols[1])
  })
  
  
  output$tax_selector <- renderUI({
    req(data())  # wait until full phyloseq object is built
    
    mdf <- tryCatch(psmelt(data()), error = function(e) return(NULL))
    if (is.null(mdf)) return(NULL)
    
    # Extract taxonomy columns from psmelt (exclude known non-taxonomic fields)
    tax_cols <- setdiff(colnames(mdf), c("OTU", "Sample", "Abundance", "Count", "value", "variable"))
    
    if (length(tax_cols) == 0) return(NULL)
    
    selectInput("tax_level", "Select Taxonomy Level", choices = tax_cols, selected = tax_cols[1])
  })
  
  
  
  output$upload_status <- renderPrint({
    if (!is.null(input$phylo)) {
      cat("RDS file uploaded successfully")
    } else if (!is.null(input$asv) && !is.null(input$tax) && !is.null(input$meta)) {
      cat("CSV files uploaded successfully")
    } else {
      cat("Waiting for file uploads...")
    }
  })
  
  
  reorder_factor_column <- function(p, group_var, order_string) {
    if (order_string != "" && group_var %in% colnames(p$data)) {
      custom_order <- trimws(unlist(strsplit(order_string, ",")))
      p$data[[group_var]] <- factor(p$data[[group_var]], levels = custom_order)
    }
    return(p)
  }
  
  
  output$rarefactionPlot <- renderPlot({
    ps <- data()
    req(input$rare_color)
    
    p <- ggrare(ps, step = 1000, color = input$rare_color, label = "Sample", se = FALSE) +
      theme_minimal()
    
    if (!is.null(input$rare_facet) && input$rare_facet != "None") {
      p <- p + facet_wrap(as.formula(paste("~", input$rare_facet)))
    }
    
    p
  })
  
  
  
  output$alphaPlot <- renderPlot({
    req(input$alpha_index, input$alpha_group)
    ps <- data()
    
    # Create plot
    p <- alpha_diversity_graph(ps, index = input$alpha_index,
                               treatment = input$alpha_group, subset = NULL, colors = "default")
    p <- reorder_factor_column(p, input$alpha_group, input$alpha_order)
    
    # Reorder factor levels if specified
    if (input$alpha_order != "") {
      custom_order <- trimws(unlist(strsplit(input$alpha_order, ",")))  # comma-separated input
      group_var <- input$alpha_group
      
      if (group_var %in% colnames(p$data)) {
        p$data[[group_var]] <- factor(p$data[[group_var]], levels = custom_order)
      }
    }
    
    p
  })
  
  
  
  output$betaPlot <- renderPlot({
    ps <- data()
    dist <- phyloseq::distance(ps, method = "bray")
    ord <- ordinate(ps, method = "PCoA", distance = dist)
    plot_ordination(ps, ord, color = "SampleID") + theme_minimal()
  })
  
  output$abundancePlot <- renderPlot({
    req(input$tax_level)  # Ensure input is ready
    
    ps <- data()
    ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
    mdf <- psmelt(ps_rel)
    
    # Check if input$tax_level exists in the melted dataframe
    if (!(input$tax_level %in% colnames(mdf))) {
      showNotification("Selected taxonomy level not found in data.", type = "error")
      return(NULL)
    }
    
    # Build the plot using aes_string
    ggplot(mdf, aes_string(x = "Sample", y = "Abundance", fill = input$tax_level)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(title = paste("Relative Abundance at", input$tax_level))
  })
  
  
  output$ordinationPlot <- renderPlot({
    ps <- data()
    ord <- ordinate(ps, method = "NMDS", distance = "bray")
    plot_ordination(ps, ord, color = "SampleID") + theme_minimal()
  })
}

shinyApp(ui = ui, server = server)