

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

ui <- fluidPage(
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "favicon.png")
  ),
  navbarPage("MetaX",
                 id = "main_nav",
                 theme = bs_theme(
                   bootswatch = "morph", # https://bootswatch.com/
                   primary = "#006699",
                   base_font = font_google("Inter")
                 ),
                 tabPanel("Home",
                          fluidRow(
                            column(10, offset = 1,
                                   h2("Welcome to MetaX"),
                                   p("This application allows you to explore microbial community data using various visualizations and analyses."),
                                   tags$ul(
                                     tags$li("💾 Start by uploading your ASV, taxonomy, and metadata tables(.csv) or phyloseq object(.rds) under 'Upload Data'."),
                                     tags$li("🧪 Remove unwanted taxa and rarefy at filtering tab"),
                                     tags$li("📈 Explore abundance, diversity, and dendrograms in the respective tabs."),
                                     tags$li("🎯 Customize plots with the sidebar controls."),
                                     tags$li("🧬 All plots support dynamic interaction based on your metadata."),
                                     tags$li("🧑‍💻❌ No coding experience required — just upload your files and explore!")
                                   ),
                                   br(),
                                   p("For questions, contact the ",
                                     tags$a(href = "https://shanptom.github.io", target = "_blank", "developer on GitHub"),
                                     "."),
                                   
                                   p("This application was built using the following R packages: ",
                                     strong("shiny"), ", ",
                                     strong("phyloseq"), ", ",
                                     strong("microeco"), ", ",
                                     strong("phylosmith"), ", ",
                                     strong("vegan"), ", ",
                                     strong("ggplot2"), ", ",
                                     strong("RColorBrewer"), ", and ",
                                     strong("bslib"), ".")
                                   
                            )
                          )
                 ),
                 
                 tabPanel("Upload Data",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("asv", "Upload ASV Table (CSV)", accept = ".csv"),
                              fileInput("tax", "Upload Taxonomy Table (CSV)", accept = ".csv"),
                              fileInput("meta", "Upload Metadata Table (CSV)", accept = ".csv"),
                              fileInput("phylo", "Or Upload Phyloseq Object (.rds)", accept = ".rds")
                            ),
                            mainPanel(verbatimTextOutput("upload_status"))
                          )
                 ),
                 tabPanel("Filtering",
                          sidebarLayout(
                            sidebarPanel(
                              checkboxInput("skip_filter", "Skip Filtering"),
                              checkboxInput("doRarefy", "Apply rarefaction", value = FALSE),
                              uiOutput("taxa_filters"),
                              actionButton("apply_filter", "Apply Filtering"),
                              actionButton("go_analysis", "Go to Analysis")
                            ),
                            mainPanel(verbatimTextOutput("filter_status"))
                          )
                 ),
                 tabPanel("Rarefaction Plot",
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("rarefaction_color_selector"),
                              uiOutput("rarefaction_facet_selector"),
                              sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12)
                            ),
                            mainPanel(plotOutput("rarefactionPlot",height = "1000px", width = "100%"))
                          )
                 ),
                 tabPanel("Abundance",
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("abund_plot_type", "Plot Type:",
                                           choices = c("Bar" = "bar", "Line" = "line","Heatmap" = "heat"),
                                           selected = "bar"),
                              uiOutput("tax_rank_selector"),
                              sliderInput("ntaxa", "Number of Top Taxa:", min = 5, max = 15, value = 8, step = 1),
                              uiOutput("abundance_facet_selector"),
                              textInput("abund_order", "Custom Order (comma-separated):", value = ""),
                              sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12),
                              checkboxInput("flip_abundance", "Flip axes (horizontal plot)", value = FALSE),
                            ),
                            mainPanel(plotOutput("abundancePlot", height = "1000px", width = "100%"))
                          )
                 ),
                 tabPanel("Dendrogram",
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("dend_treatment_selector"),
                              selectInput("dend_method", "Select distance method:",
                                          choices = c("euclidian", "manhattan", "canberra", "clark", "bray",
                                                      "kulczynski", "jaccard", "gower", "altGower", "morisita",
                                                      "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"),
                                          selected = "bray"),
                              sliderInput("dend_label_size", "Label size:", min = 3, max = 10, value = 5, step = 1)
                            ),

                            mainPanel(plotOutput("dendrogramPlot", height = "950px", width = "100%")))
                          ),
                 tabPanel("Alpha Diversity",
                          sidebarLayout(
                            sidebarPanel(
                              checkboxGroupInput("alpha_index", "Select Diversity Index:",
                                                 choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"),
                                                 selected = c("Shannon")),
                              uiOutput("alpha_group_selector"),
                              uiOutput("alpha_colour_selector"),
                              checkboxInput("flip_alpha", "Flip axes (horizontal plot)", value = FALSE),
                              textInput("alpha_order", "Custom order (comma-separated values)", value = ""),
                              sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12)
                            ),
                            mainPanel(plotOutput("alphaPlot",height = "950px", width = "100%"))
                          )
                 ),
                 tabPanel("Beta Diversity",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("beta_dist", "Distance Method:",
                                          choices = c("bray", "unifrac", "wunifrac", "jaccard", "dpcoa", "jsd", "manhattan",
                                                      "euclidean", "canberra", "binomial"),
                                          selected = "bray"),
                              selectInput("beta_ord", "Ordination Method:",
                                          choices = c("NMDS", "MDS", "PCoA", "DCA", "CCA", "RDA", "CAP", "DPCoA"),
                                          selected = "NMDS"),
                              uiOutput("beta_color_selector"),
                              uiOutput("beta_shape_selector"),
                              uiOutput("beta_label_selector"),
                              uiOutput("beta_facet_selector"),
                              sliderInput("beta_label_size", "Axis Text Size:", min = 6, max = 20, value = 12),
                              sliderInput("beta_label_text_size", "Label Text Size:", min = 2, max = 15, value = 3),
                              sliderInput("beta_shape_size", "Shape Size:", min = 1, max = 10, value = 4),
                              tags$hr(),
                              h4("PERMANOVA"),
                              uiOutput("permanova_group_selector"),
                              actionButton("run_permanova", "Run PERMANOVA"),
                              verbatimTextOutput("permanova_result")
                            ),
                            mainPanel(
                              plotOutput("betaPlot", height = "950px", width = "100%"))
                          )
                 )

))

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
      textInput(paste0("filter_", rank), paste("Exclude", rank, "(comma-separated):"), "")
    })
  })
  
  observeEvent(input$apply_filter, {
    ps <- raw_physeq()
    taxdf <- as.data.frame(tax_table(ps))
    ranks <- colnames(taxdf)
    to_remove <- rep(FALSE, ntaxa(ps))
    
    for (rank in ranks) {
      vals <- input[[paste0("filter_", rank)]]
      if (nzchar(vals)) {
        exclude_list <- trimws(unlist(strsplit(vals, ",")))
        to_remove <- to_remove | taxdf[[rank]] %in% exclude_list
      }
    }
    
    ps <- prune_taxa(!to_remove, ps)
    df <- as.data.frame(sample_data(ps))
    for (var in colnames(df)) {
      if (is.character(df[[var]]) || is.factor(df[[var]])) ordering_rules[[var]] <- unique(df[[var]])
    }
    ordering_rules$Sample <- sample_names(ps)
    final_physeq(ps)
    
    if (any(sample_sums(ps) == 0)) {
      showNotification("Warning: some samples have zero reads after filtering.", type = "warning")
    }
  })
  
  observeEvent(input$skip_filter, {
    ps <- raw_physeq()
    df <- as.data.frame(sample_data(ps))
    for (var in colnames(df)) {
      if (is.character(df[[var]]) || is.factor(df[[var]])) ordering_rules[[var]] <- unique(df[[var]])
    }
    ordering_rules$Sample <- sample_names(ps)
    final_physeq(ps)
  })
  
  observeEvent(input$go_analysis, {
    req(final_physeq())
    if (input$doRarefy) {
      ps <- rarefy_even_depth(final_physeq(),
                              sample.size = min(sample_sums(final_physeq())),
                              rngseed = 123, replace = TRUE,
                              trimOTUs = TRUE, verbose = FALSE
      )
      final_physeq(ps)
    }
    updateNavbarPage(session, "main_nav", selected = "Rarefaction Plot")
  })
  
  output$filter_status <- renderPrint({
    req(final_physeq())
    cat("Current number of ASVs:", ntaxa(final_physeq()))
  })
  
  output$rarefaction_color_selector <- renderUI({
    req(final_physeq())
    cols <- colnames(sample_data(final_physeq()))
    selectInput("rare_color", "Color by:", choices = cols, selected = cols[1])
  })
  
  output$rarefaction_facet_selector <- renderUI({
    req(final_physeq())
    cols <- colnames(sample_data(final_physeq()))
    selectInput("rare_facet", "Facet by:", choices = c("None", cols), selected = "None")
  })
  
  reorder_factor_column <- function(p, group_var, order_string) {
    if (order_string != "" && group_var %in% colnames(p$data)) {
      custom_order <- trimws(unlist(strsplit(order_string, ",")))
      p$data[[group_var]] <- factor(p$data[[group_var]], levels = custom_order)
    }
    return(p)
  }
  
  output$rarefactionPlot <- renderPlot({
    req(final_physeq(), input$rare_color)
    ps <- final_physeq()
    otu_table(ps) <- otu_table(round(otu_table(ps)), taxa_are_rows = TRUE)
    if (any(sample_sums(ps) == 0)) {
      showNotification("Some samples have 0 counts. Rarefaction plot may not work.", type = "error")
      return(NULL)
    }
    p <- ggrare(ps, step = 1000, color = input$rare_color, label = "Sample", se = FALSE) + theme_minimal()
    if (!is.null(input$rare_facet) && input$rare_facet != "None") {
      p <- p + facet_wrap(as.formula(paste("~", input$rare_facet)))
    }
    p +
      theme(
        axis.text = element_text(size = input$beta_label_size),
        axis.title = element_text(size = input$beta_label_size),
        legend.text = element_text(size = input$beta_label_size),
        legend.title = element_text(size = input$beta_label_size)
      )
  })
  
  output$tax_rank_selector <- renderUI({
    req(final_physeq())
    tax_ranks <- colnames(as.data.frame(tax_table(final_physeq())))
    selectInput("tax_rank", "Select Taxonomic Rank:", choices = tax_ranks, selected = tail(tax_ranks, 1))
  })
  
  output$abundance_facet_selector <- renderUI({
    req(final_physeq())
    meta_cols <- colnames(sample_data(final_physeq()))
    selectInput("abund_facet", "Facet by:", choices = c("None", meta_cols), selected = "None")
  })
  
  output$abundancePlot <- renderPlot({
    req(final_physeq(), input$tax_rank, input$ntaxa)
    dataset <- phyloseq2meco(final_physeq())
    dataset$tidy_dataset()
    dataset$cal_abund()
    t1 <- trans_abund$new(dataset = dataset, taxrank = input$tax_rank, ntaxa = input$ntaxa)
    scale_type <- if (input$flip_abundance) "free_y" else "free_x"
    
    if (input$abund_plot_type == "line") {
      p4 <- t1$plot_bar(
        bar_type = "notfull",
        use_alluvium = TRUE,
        clustering = TRUE,
        xtext_angle = 90,
        color_values = RColorBrewer::brewer.pal(8, "Set2")
      )
    } else if (input$abund_plot_type == "heat") {
      p4 <- t1$plot_heatmap(
        xtext_keep = FALSE,
        withmargin = FALSE,
        plot_breaks = c(0.01, 0.1, 1, 10)
      )
    } else {
      p4 <- t1$plot_bar(
        others_color = "grey70",
        xtext_angle = 90,
        legend_text_italic = FALSE
      )
    }
    
    for (var in names(ordering_rules)) {
      if (var %in% colnames(p4$data)) {
        p4$data[[var]] <- factor(p4$data[[var]], levels = ordering_rules[[var]])
      }
    }
    
    if (!is.null(input$abund_facet) && input$abund_facet != "None") {
      facet_formula <- as.formula(paste("~", input$abund_facet))
      p4 <- p4 + facet_wrap(facet_formula, scales = scale_type)
    }
    # Apply custom ordering if specified
    if (nzchar(input$abund_order)) {
      custom_order <- trimws(unlist(strsplit(input$abund_order, ",")))
      
      if (input$abund_facet %in% colnames(p4$data)) {
        p4$data[[input$abund_facet]] <- factor(p4$data[[input$abund_facet]], levels = custom_order)
      } else if ("Sample" %in% colnames(p4$data)) {
        p4$data$Sample <- factor(p4$data$Sample, levels = custom_order)
      }
    }
    
    if (input$flip_abundance) {
      p4 <- p4 + coord_flip()
    }
    
    p4 +
      theme(
        axis.text = element_text(size = input$beta_label_size),
        axis.title = element_text(size = input$beta_label_size),
        legend.text = element_text(size = input$beta_label_size),
        legend.title = element_text(size = input$beta_label_size)
      )
  })
  
  output$alpha_group_selector <- renderUI({
    req(final_physeq())
    meta_cols <- colnames(sample_data(final_physeq()))
    selectInput("alpha_group", "Group by:",
                choices = c("None", meta_cols),
                selected = "None")
  })
  
  output$alpha_colour_selector <- renderUI({
    req(final_physeq())
    meta_cols <- colnames(sample_data(final_physeq()))
    selectInput("alpha_colour", "Colour",  choices = c("None", meta_cols),
                selected = "None")
  })
  
  
  output$alphaPlot <- renderPlot({
    req(final_physeq())
    validate(need(length(input$alpha_index) > 0, "Please select at least one diversity index."))
    
    x_var <- if (input$alpha_group == "None") "samples" else input$alpha_group
    scale_type <- if (input$flip_alpha) "free_x" else "free_y"
    
    # Base plot
    p <- plot_richness(
      final_physeq(),
      x = x_var,
      measures = input$alpha_index,
      scales = scale_type
    )
    
    # Apply color only if valid
    if (input$alpha_colour != "None") {
      p <- p + aes_string(color = input$alpha_colour)
    }
    
    # Reorder if custom order given
    p <- reorder_factor_column(p, x_var, input$alpha_order)
    
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    if (input$flip_alpha) {
      p <- p + coord_flip()
    }
    
    p +
      theme(
        axis.text = element_text(size = input$beta_label_size),
        axis.title = element_text(size = input$beta_label_size),
        legend.text = element_text(size = input$beta_label_size),
        legend.title = element_text(size = input$beta_label_size),
        strip.text = element_text(size = input$beta_label_size)
      )
  })
  
  output$beta_color_selector <- renderUI({
    req(final_physeq())
    cols <- colnames(sample_data(final_physeq()))
    selectInput("beta_color", "Color by:", choices = c("None", cols), selected = "None")
  })
  
  output$beta_shape_selector <- renderUI({
    req(final_physeq())
    cols <- colnames(sample_data(final_physeq()))
    selectInput("beta_shape", "Shape by:", choices = c("None", cols), selected = "None")
  })
  
  output$beta_label_selector <- renderUI({
    req(final_physeq())
    cols <- colnames(sample_data(final_physeq()))
    selectInput("beta_label", "Label points by:", choices = c("None", cols), selected = "None")
  })
  
  output$beta_facet_selector <- renderUI({
    req(final_physeq())
    cols <- colnames(sample_data(final_physeq()))
    selectInput("beta_facet", "Facet by:", choices = c("None", cols), selected = "None")
  })
  
  
  
  output$betaPlot <- renderPlot({
    req(final_physeq())
    
    dist <- distance(final_physeq(), method = input$beta_dist)
    ord <- ordinate(final_physeq(), method = input$beta_ord, distance = dist)
    
    p <- plot_ordination(final_physeq(), ord, 
                         color = input$beta_color, 
                         shape = if (input$beta_shape != "None") input$beta_shape else NULL) +
      theme_minimal() +
      geom_point(size = input$beta_shape_size)
    
    #  labels only if selected
    if (!is.null(input$beta_label) && input$beta_label != "None") {
      p <- p + geom_text(aes_string(label = input$beta_label), 
                         size = input$beta_label_text_size, 
                         vjust = -1)
    }
    
    # Apply theme and text size
    p <- p + theme(
      axis.text = element_text(size = input$beta_label_size),
      axis.title = element_text(size = input$beta_label_size),
      legend.text = element_text(size = input$beta_label_size),
      legend.title = element_text(size = input$beta_label_size),
      strip.text = element_text(size = input$beta_label_size)
    )
    

    if (!is.null(input$beta_facet) && input$beta_facet != "None") {
      p <- p + facet_wrap(as.formula(paste("~", input$beta_facet)))
    }
    
    p
  })
  
  output$permanova_group_selector <- renderUI({
    req(final_physeq())
    selectInput("permanova_group", "Select grouping variable:",
                choices = names(sample_data(final_physeq())),
                selected = names(sample_data(final_physeq()))[1])
  })
  
  permanova_result <- reactiveVal()
  
  observeEvent(input$run_permanova, {
    req(final_physeq(), input$permanova_group)
    
    dataset <- phyloseq2meco(final_physeq())
    dataset$tidy_dataset()
    dataset$cal_betadiv()
    
    t1 <- trans_beta$new(
      dataset = dataset,
      group = input$permanova_group,
      measure = "bray"
    )
    
    t1$cal_manova(manova_all = FALSE)
    
    permanova_result(t1$res_manova)
  })

  output$permanova_result <- renderPrint({
    req(permanova_result())
    permanova_result()
  })
  
  dendrogram_phyloseq_custom <- function(phyloseq_obj, treatment = NULL, method = "bray",
                                         colors = "default", label_size = 2.5) {
    dend <- phylosmith::dendrogram_phyloseq(
      phyloseq_obj = phyloseq_obj,
      treatment = treatment,
      method = method,
      colors = colors
    )
    
    # Identify and replace the geom_label layer with new size
    label_data <- NULL
    label_mapping <- NULL
    
    for (layer in dend$layers) {
      if (inherits(layer$geom, "GeomLabel") || inherits(layer$geom, "GeomText")) {
        label_data <- layer$data
        label_mapping <- layer$mapping
      }
    }
    
    if (!is.null(label_data)) {
      dend$layers <- Filter(function(l) {
        !inherits(l$geom, "GeomLabel") && !inherits(l$geom, "GeomText")
      }, dend$layers)
      
      dend <- dend + 
        ggplot2::geom_label(
          data = label_data,
          mapping = label_mapping,
          size = label_size,
          color = "white",
          label.padding = unit(0.2, "lines"),
          fontface = "bold",
          hjust = 1.05
        )
    }
    
    return(dend)
  }
  
  output$dend_treatment_selector <- renderUI({
    req(final_physeq())
    selectInput("dend_treatment", "Select metadata column for grouping:",
                choices = names(sample_data(final_physeq())),
                selected = names(sample_data(final_physeq()))[1])
  })
  
  output$dendrogramPlot <- renderPlot({
    req(final_physeq(), input$dend_method, input$dend_treatment, input$dend_label_size)
    
    dend <- dendrogram_phyloseq_custom(
      phyloseq_obj = final_physeq(),
      treatment = input$dend_treatment,
      method = input$dend_method,
      label_size = input$dend_label_size
    )
    
    dend
  })
  
  
  
}

shinyApp(ui = ui, server = server)
