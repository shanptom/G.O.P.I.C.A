library(shiny)
library(bslib)
library(phyloseq)
library(ggplot2)
library(vegan)
library(ranacapa)
library(phylosmith)
library(microeco)
library(file2meco)
library(RColorBrewer)

ui <- navbarPage("MetaX",
                 id = "main_nav",
                 theme = bs_theme(
                   bootswatch = "cerulean",
                   primary = "#006699",
                   base_font = font_google("Inter")
                 ),
                 tabPanel("Upload Data",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("asv", "Upload ASV Table (CSV)", accept = ".csv"),
                              fileInput("tax", "Upload Taxonomy Table (CSV)", accept = ".csv"),
                              fileInput("meta", "Upload Metadata Table (CSV)", accept = ".csv"),
                              fileInput("phylo", "Or Upload Phyloseq Object (.rds)", accept = ".rds")
                            ),
                            mainPanel(
                              verbatimTextOutput("upload_status")
                            )
                          )
                 ),
                 tabPanel("Filtering",
                          sidebarLayout(
                            sidebarPanel(
                              actionButton("skip_filter", "Skip Filtering"),
                              actionButton("go_analysis", "Go to Analysis"),
                              checkboxInput("doRarefy", "Apply rarefaction", value = FALSE),
                              uiOutput("taxa_filters"),
                              actionButton("apply_filter", "Apply Filtering"),
                              actionButton("go_analysis", "Go to Analysis")
                            ),
                            mainPanel(
                              verbatimTextOutput("filter_status")
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
                 tabPanel("Abundance",
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("tax_rank_selector"),
                              sliderInput("ntaxa", "Number of Top Taxa:", min = 5, max = 15, value = 8, step = 1)),
                            mainPanel(
                              plotOutput("abundancePlot")
                            )
                          )
                 ),
                 tabPanel("Alpha Diversity",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("alpha_index", "Select Diversity Index:",
                                          choices = c("shannon", "simpson", "invsimpson"), selected = "shannon"),
                              uiOutput("alpha_group_selector"),
                              textInput("alpha_order", "Custom order (comma-separated values)", value = "")
                            ),
                            mainPanel(
                              plotOutput("alphaPlot")
                            )
                          )
                 ),
                 tabPanel("Beta Diversity", plotOutput("betaPlot")),
                 tabPanel("Ordination", plotOutput("ordinationPlot"))
)

server <- function(input, output, session) {
  final_physeq <- reactiveVal()
  ordering_rules <- reactiveValues()
  
  raw_physeq <- reactive({
    if (!is.null(input$phylo)) {
      readRDS(input$phylo$datapath)
    } else {
      req(input$asv, input$tax, input$meta)
      otu <- read.csv(input$asv$datapath, row.names = 1)
      tax <- as.matrix(read.csv(input$tax$datapath, row.names = 1))
      meta <- read.csv(input$meta$datapath, row.names = 1)
      phyloseq(
        otu_table(as.matrix(otu), taxa_are_rows = TRUE),
        tax_table(tax),
        sample_data(meta)
      )
    }
  })
  
  output$upload_status <- renderPrint({
    if (!is.null(input$phylo)) {
      cat("RDS file uploaded successfully")
      updateNavbarPage(session, "main_nav", selected = "Filtering")
    } else if (!is.null(input$asv) && !is.null(input$tax) && !is.null(input$meta)) {
      cat("CSV files uploaded successfully")
      updateNavbarPage(session, "main_nav", selected = "Filtering")
    } else {
      cat("Waiting for file uploads...")
    }
  })
  
  output$taxa_filters <- renderUI({
    req(raw_physeq())
    ranks <- colnames(as.data.frame(tax_table(raw_physeq())))
    lapply(ranks, function(rank) {
      textInput(inputId = paste0("filter_", rank),
                label = paste("Exclude", rank, "(comma-separated):"),
                value = "")
    })
  })
  
  observeEvent(input$apply_filter, {
    ps <- raw_physeq()
    taxdf <- as.data.frame(tax_table(ps))
    ranks <- colnames(taxdf)
    to_remove <- rep(FALSE, ntaxa(ps))
    
    for (rank in ranks) {
      input_vals <- input[[paste0("filter_", rank)]]
      if (nzchar(input_vals)) {
        exclude_list <- trimws(unlist(strsplit(input_vals, ",")))
        to_remove <- to_remove | taxdf[[rank]] %in% exclude_list
      }
    }
    
    ps <- prune_taxa(!to_remove, ps)
    
    if (input$doRarefy) {
      ps <- rarefy_even_depth(ps,
                              sample.size = min(sample_sums(ps)),
                              rngseed = 123,
                              replace = TRUE,
                              trimOTUs = TRUE,
                              verbose = FALSE)
    }
    
    df <- as.data.frame(sample_data(ps))
    for (var in colnames(df)) {
      if (is.character(df[[var]]) || is.factor(df[[var]])) {
        ordering_rules[[var]] <- unique(df[[var]])
      }
    }
    ordering_rules$Sample <- sample_names(ps)
    final_physeq(ps)
    
    if (any(sample_sums(ps) == 0)) {
      showNotification("Warning: some samples have zero reads after filtering.", type = "warning")
    }
    
  })
  
  
  
  observeEvent(input$skip_filter, {
    ps <- raw_physeq()
    if (input$doRarefy) {
      ps <- rarefy_even_depth(ps,
                              sample.size = min(sample_sums(ps)),
                              rngseed = 123,
                              replace = TRUE,
                              trimOTUs = TRUE,
                              verbose = FALSE)
    }
    df <- as.data.frame(sample_data(ps))
    for (var in colnames(df)) {
      if (is.character(df[[var]]) || is.factor(df[[var]])) {
        ordering_rules[[var]] <- unique(df[[var]])
      }
    }
    ordering_rules$Sample <- sample_names(ps)
    final_physeq(ps)
  })
  
  observeEvent(input$go_analysis, {
    req(final_physeq())
    updateNavbarPage(session, "main_nav", selected = "Rarefaction Plot")
  })
  
  output$filter_status <- renderPrint({
    req(final_physeq())
    cat("Current number of ASVs:", ntaxa(final_physeq()))
  })
  
  output$rarefactionPlot <- renderPlot({
    req(final_physeq(), input$rare_color)
    
    ps <- final_physeq()
    
    # Ensure counts are integers
    otu_table(ps) <- otu_table(round(otu_table(ps)), taxa_are_rows = TRUE)
    
    # Double-check sample_sums
    if (any(sample_sums(ps) == 0)) {
      showNotification("Some samples have 0 counts. Rarefaction plot may not work.", type = "error")
      return(NULL)
    }
    
    p <- ggrare(ps, step = 1000, color = input$rare_color, label = "Sample", se = FALSE) +
      theme_minimal()
    
    if (!is.null(input$rare_facet) && input$rare_facet != "None") {
      p <- p + facet_wrap(as.formula(paste("~", input$rare_facet)))
    }
    
    p
  })
  
  
  output$rarefaction_color_selector <- renderUI({
    req(final_physeq())
    meta_cols <- colnames(sample_data(final_physeq()))
    selectInput("rare_color", "Color by (metadata column):", choices = meta_cols, selected = meta_cols[1])
  })
  
  output$rarefaction_facet_selector <- renderUI({
    req(final_physeq())
    meta_cols <- colnames(sample_data(final_physeq()))
    selectInput("rare_facet", "Facet by (metadata column):", choices = c("None", meta_cols), selected = "None")
  })
  
  output$alpha_group_selector <- renderUI({
    req(final_physeq())
    meta_cols <- colnames(sample_data(final_physeq()))
    selectInput("alpha_group", "Group by (metadata column):", choices = meta_cols, selected = meta_cols[1])
  })
  
  
  output$tax_rank_selector <- renderUI({
    req(final_physeq())
    tax_ranks <- colnames(as.data.frame(tax_table(final_physeq())))
    selectInput("tax_rank", "Select Taxonomic Rank:",
                choices = tax_ranks,
                selected = tail(tax_ranks, 1))  # default to most specific
  })
  
  output$abundancePlot <- renderPlot({
    req(final_physeq())
    
    dataset <- phyloseq2meco(final_physeq())
    dataset$tidy_dataset()
    dataset$cal_abund()
    
    t1 <- trans_abund$new(
      dataset = dataset,
      taxrank = input$tax_rank,
      ntaxa = input$ntaxa
    )
    p4 <- t1$plot_bar(
      bar_type = "notfull",
      use_alluvium = TRUE,
      clustering = TRUE,
      xtext_angle = 90,
      xtext_size = 6,
      color_values = RColorBrewer::brewer.pal(8, "Set2")
    )
    
    for (var in names(ordering_rules)) {
      if (var %in% colnames(p4$data)) {
        p4$data[[var]] <- factor(p4$data[[var]], levels = ordering_rules[[var]])
      }
    }
    
    p4
  })
  
  output$betaPlot <- renderPlot({
    req(final_physeq())
    dist <- phyloseq::distance(final_physeq(), method = "bray")
    ord <- ordinate(final_physeq(), method = "PCoA", distance = dist)
    plot_ordination(final_physeq(), ord, color = "SampleID") + theme_minimal()
  })
  
  output$alphaPlot <- renderPlot({
    req(final_physeq(), input$alpha_index, input$alpha_group)
    p <- alpha_diversity_graph(final_physeq(), index = input$alpha_index,
                               treatment = input$alpha_group, subset = NULL, colors = "default")
    p <- reorder_factor_column(p, input$alpha_group, input$alpha_order)
    p
  })
  
  output$ordinationPlot <- renderPlot({
    req(final_physeq())
    ord <- ordinate(final_physeq(), method = "NMDS", distance = "bray")
    plot_ordination(final_physeq(), ord, color = "SampleID") + theme_minimal()
  })
  
  output$tax_selector <- renderUI({
    req(final_physeq())
    mdf <- tryCatch(psmelt(final_physeq()), error = function(e) return(NULL))
    if (is.null(mdf)) return(NULL)
    tax_cols <- setdiff(colnames(mdf), c("OTU", "Sample", "Abundance", "Count", "value", "variable"))
    if (length(tax_cols) == 0) return(NULL)
    selectInput("tax_level", "Select Taxonomy Level", choices = tax_cols, selected = tax_cols[1])
  })
  
  reorder_factor_column <- function(p, group_var, order_string) {
    if (order_string != "" && group_var %in% colnames(p$data)) {
      custom_order <- trimws(unlist(strsplit(order_string, ",")))
      p$data[[group_var]] <- factor(p$data[[group_var]], levels = custom_order)
    }
    return(p)
  }
}

shinyApp(ui = ui, server = server)
