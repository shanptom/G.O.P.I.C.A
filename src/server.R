library(phyloseq)
library(ggplot2)
library(vegan)
library(ranacapa)
library(phylosmith)
library(microeco)
library(file2meco)
library(GUniFrac)
library(RColorBrewer)
library(ggalluvial)
library(dplyr)
library(ggcor)
library(ggpubr)
library(plotly)

server <- function(input, output, session) {
  final_physeq <- reactiveVal()
  ordering_rules <- reactiveValues()
  reactiveValues_envfit <- reactiveValues(transenv = NULL)
  selected_analysis <- reactiveVal(NULL)
  show_tsne <- reactiveVal(FALSE)
  
  # Initially disable all analysis tabs
  analysis_tabs <- c("filter_tab", "rarefaction_tab", "abundance_tab", "alpha_tab", "dendrogram_tab", "ordination_tab", "metadata_tab", "regression_tab", "indicator_tab")
  
  observe({
    for(tab in analysis_tabs) {
      shinyjs::disable(selector = paste0("a[data-value='", tab, "']"))
    }
  })
  
  observeEvent(input$run_rda, {
    selected_analysis("rda")
  })
  
  observeEvent(input$run_corr, {
    selected_analysis("corr")
  })
  
  observeEvent(input$run_mantel, {
    selected_analysis("mantel")
  })
  
  observeEvent(input$run_tsne, {
    show_tsne(TRUE)
  })
  
  output$show_tsne_flag <- reactive({
    show_tsne()
  })
  outputOptions(output, "show_tsne_flag", suspendWhenHidden = FALSE)
  
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
  
  output$upload_status_ui <- renderUI({
    if (!is.null(input$phylo) || (!is.null(input$asv) && !is.null(input$tax) && !is.null(input$meta))) {
      # Enable analysis tabs
      for(tab in analysis_tabs) {
        shinyjs::enable(selector = paste0("a[data-value='", tab, "']"))
      }
      
      # Show success notification
      showNotification("Files uploaded successfully!", type = "message", duration = 5)
      
      # Update the navbar to the filter tab
      updateNavbarPage(session, "main_nav", selected = "filter_tab")
      
      # Process the data
      if (!is.null(input$phylo)) {
        phy <- readRDS(input$phylo$datapath)
        final_physeq(phy)
        df <- as.data.frame(sample_data(phy))
      } else {
        otu <- read.csv(input$asv$datapath, row.names = 1)
        tax <- as.matrix(read.csv(input$tax$datapath, row.names = 1))
        meta <- read.csv(input$meta$datapath, row.names = 1)
        ps <- phyloseq(
          otu_table(as.matrix(otu), taxa_are_rows = TRUE),
          tax_table(tax),
          sample_data(meta)
        )
        final_physeq(ps)
        df <- as.data.frame(sample_data(ps))
      }
      
      # Reset ordering rules
      for (var in colnames(df)) {
        if (is.character(df[[var]]) || is.factor(df[[var]])) ordering_rules[[var]] <- unique(df[[var]])
      }
      ordering_rules$Sample <- sample_names(final_physeq())
      
      return(h4("Upload successful. Please proceed to the next tabs."))
    } else {
      return(h4("Waiting for file uploads or demo data loading..."))
    }
  })
  
  # Handle loading of demo data
  observeEvent(input$load_demo, {
    demo_path <- "data/"
    if (input$demo_file == "rds") {
      phy <- readRDS(file.path(demo_path, "demo_ps.rds"))
      final_physeq(phy)
      df <- as.data.frame(sample_data(phy))
      showNotification("Demo Phyloseq RDS loaded successfully!", type = "message", duration = 5)
    } else if (input$demo_file == "csv") {
      otu <- read.csv(file.path(demo_path, "demo_asv.csv"), row.names = 1)
      tax <- as.matrix(read.csv(file.path(demo_path, "demo_tax.csv"), row.names = 1))
      meta <- read.csv(file.path(demo_path, "demo_meta.csv"), row.names = 1)
      ps <- phyloseq(
        otu_table(as.matrix(otu), taxa_are_rows = TRUE),
        tax_table(tax),
        sample_data(meta)
      )
      final_physeq(ps)
      df <- as.data.frame(sample_data(ps))
      showNotification("Demo CSV files loaded successfully!", type = "message", duration = 5)
    }
    
    # Enable analysis tabs
    for(tab in analysis_tabs) {
      shinyjs::enable(selector = paste0("a[data-value='", tab, "']"))
    }
    
    # Update the navbar to the filter tab
    updateNavbarPage(session, "main_nav", selected = "filter_tab")
    
    # Reset ordering rules
    for (var in colnames(df)) {
      if (is.character(df[[var]]) || is.factor(df[[var]])) ordering_rules[[var]] <- unique(df[[var]])
    }
    ordering_rules$Sample <- sample_names(final_physeq())
  })
  
  
  output$taxa_filters <- renderUI({
    req(final_physeq())
    ranks <- colnames(as.data.frame(tax_table(final_physeq())))
    lapply(ranks, function(rank) {
      textInput(paste0("filter_", rank), paste("Exclude", rank, "(comma-separated):"), "")
    })
  })
  
  observeEvent(input$apply_filter, {
    ps <- final_physeq()
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
    } else {
      showNotification("Filtering applied successfully.", type = "message")
    }
  })
  
  observeEvent(input$go_analysis, {
    req(final_physeq())
    ps <- final_physeq()
    
    if (input$doRarefy) {
      ps <- rarefy_even_depth(ps,
                              sample.size = min(sample_sums(ps)),
                              rngseed = 123, replace = TRUE,
                              trimOTUs = TRUE, verbose = FALSE
      )
      showNotification("Rarefaction applied.", type = "message")
    } else if (input$doTSS) {
      ps <- transform_sample_counts(ps, function(x) x / sum(x))
      showNotification("TSS normalization applied.", type = "message")
    }
    
    final_physeq(ps)
    updateNavbarPage(session, "main_nav", selected = "rarefaction_tab")
  })
  
  observeEvent(input$doRarefy, {
    if (input$doRarefy && input$doTSS) {
      updateCheckboxInput(session, "doTSS", value = FALSE)
    }
  })
  
  observeEvent(input$doTSS, {
    if (input$doTSS && input$doRarefy) {
      updateCheckboxInput(session, "doRarefy", value = FALSE)
    }
  })
  
  
  output$filter_status <- renderPrint({
    req(final_physeq())
    cat("Current number of ASVs:", ntaxa(final_physeq()))
  })
  
  output$rarefaction_color_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("rare_color", "Color by:", choices = categorical_cols, selected = categorical_cols[1] %||% "None")
  })
  
  output$rarefaction_facet_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("rare_facet", "Facet by:", choices = c("None", categorical_cols), selected = "None")
  })

  output$rarefaction_facet_order_selector <- renderUI({
    req(final_physeq(), input$rare_facet)
    if (!is.null(input$rare_facet) && input$rare_facet != "None") {
      df <- as.data.frame(sample_data(final_physeq()))
      choices <- unique(df[[input$rare_facet]])
      selectizeInput("rarefaction_facet_order", "Custom Facet Order (drag to reorder):", 
                     choices = choices, selected = choices, multiple = TRUE, 
                     options = list(plugins = list('drag_drop')))
    } else {
      return(NULL)
    }
  })
  
  reorder_factor_column <- function(p, group_var, order_string) {
    if (order_string != "" && group_var %in% colnames(p$data)) {
      custom_order <- trimws(unlist(strsplit(order_string, ",")))
      p$data[[group_var]] <- factor(p$data[[group_var]], levels = custom_order)
    }
    return(p)
  }
  
  # Store the base rarefaction plot object to avoid recalculation
  base_rarefaction_plot <- reactiveVal(NULL)
  
  observe({
    req(final_physeq())
    ps <- final_physeq()
    otu_table(ps) <- otu_table(round(otu_table(ps)), taxa_are_rows = TRUE)
    if (any(sample_sums(ps) == 0)) {
      showNotification("Some samples have 0 counts. Rarefaction plot may not work.", type = "error")
      base_rarefaction_plot(NULL)
    } else {
      # Create the base plot without color or other aesthetics that might change
      base_plot <- ggrare(ps, step = 100, label = "Sample", se = FALSE) + theme_minimal()
      base_rarefaction_plot(base_plot)
    }
  })
  
  output$rarefactionPlot <- renderPlot({
    req(base_rarefaction_plot(), input$rare_color)
    p <- base_rarefaction_plot()
    if (is.null(p)) {
      return(NULL)
    }
    # Apply color aesthetic dynamically
    p <- p + aes(color = !!sym(input$rare_color))
    
    if (!is.null(input$rare_facet) && input$rare_facet != "None") {
      if (!is.null(input$rarefaction_facet_order) && length(input$rarefaction_facet_order) > 0) {
        # Apply custom facet order using factor with specified levels
        levels_str <- paste0("c(\"", paste(input$rarefaction_facet_order, collapse = "\", \""), "\")")
        p <- p + facet_wrap(as.formula(paste0("~factor(", input$rare_facet, ", levels = ", levels_str, ")")), scales = "free")
      } else {
        # Default facet without custom order
        p <- p + facet_wrap(as.formula(paste("~", input$rare_facet)), scales = "free")
      }
    }
    # Remove default geom_text layer if labels are not to be shown
    if (!input$show_rarefaction_labels) {
      p$layers <- lapply(p$layers, function(layer) {
        if (inherits(layer$geom, "GeomText")) {
          return(NULL)
        }
        return(layer)
      })
      p$layers <- Filter(Negate(is.null), p$layers)
    } else {
      # Override the default label size for sample labels
      p$layers <- lapply(p$layers, function(layer) {
        if (inherits(layer$geom, "GeomText")) {
          layer$aes_params$size <- input$rarefaction_label_size
        }
        return(layer)
      })
    }
    p +
      theme(
        axis.text = element_text(size = input$beta_label_size),
        axis.title = element_text(size = input$beta_label_size),
        strip.text = element_text(size = input$beta_label_size),
        legend.position = "none"
      )
  })
  
  output$tax_rank_selector <- renderUI({
    req(final_physeq())
    tax_ranks <- colnames(as.data.frame(tax_table(final_physeq())))
    selectInput("tax_rank", "Select Taxonomic Rank:", choices = tax_ranks, selected = tail(tax_ranks, 1))
  })
  
  output$abundance_facet_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("abund_facet", "Facet by:", choices = c("None", categorical_cols), selected = "None")
  })

  output$abundance_order_selector <- renderUI({
    req(final_physeq(), input$abund_facet)
    if (!is.null(input$abund_facet) && input$abund_facet != "None") {
      df <- as.data.frame(sample_data(final_physeq()))
      choices <- unique(df[[input$abund_facet]])
      selectizeInput("abund_order", "Custom Order (drag to reorder):", 
                     choices = choices, selected = choices, multiple = TRUE, 
                     options = list(plugins = list('drag_drop')))
    } else {
      choices <- sample_names(final_physeq())
      selectizeInput("abund_order", "Custom Sample Order (drag to reorder):", 
                     choices = choices, selected = choices, multiple = TRUE, 
                     options = list(plugins = list('drag_drop')))
    }
  })
  
  output$abundance_plot_output <- renderUI({
    if (input$abund_plot_type == "line" || input$abund_plot_type == "heat") {
      plotOutput("abundancePlotStatic", height = "770px", width = "100%")
    } else {
      plotlyOutput("abundancePlotly", height = "770px", width = "100%")
    }
  })

  abundance_plot_obj <- reactive({
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
      xtext_size = input$beta_label_size,
      color_values = RColorBrewer::brewer.pal(8, "Set2")
    )
  } else if (input$abund_plot_type == "bar") {
    p4 <- t1$plot_bar(
      others_color = "grey70",
      xtext_angle = 90,
      xtext_size = input$beta_label_size,
      legend_text_italic = FALSE
    )
  } else if (input$abund_plot_type == "heat") {
    p4 <- t1$plot_heatmap(
      xtext_keep = TRUE,
      xtext_angle = 90,
      xtext_size = input$beta_label_size,
      ytext_size = input$beta_label_size,
      withmargin = FALSE,
      plot_breaks = c(0.01, 0.1, 1, 10)
    )
  }

  for (var in names(ordering_rules)) {
    if (var %in% colnames(p4$data)) {
      p4$data[[var]] <- factor(p4$data[[var]], levels = ordering_rules[[var]])
    }
  }

  if (!is.null(input$abund_facet) && input$abund_facet != "None") {
    facet_formula <- as.formula(paste("~", input$abund_facet))
    if (input$abund_plot_type == "heat") {
      if (input$flip_abundance) {
        p4 <- p4 + facet_wrap(facet_formula, scales = "free_y")
      } else {
        p4 <- p4 + facet_wrap(facet_formula, scales = "free_x")
      }
    } else {
      p4 <- p4 + facet_wrap(facet_formula, scales = scale_type)
    }
  }

  if (!is.null(input$abund_order) && length(input$abund_order) > 0) {
    if (input$abund_facet != "None" && input$abund_facet %in% colnames(p4$data)) {
      p4$data[[input$abund_facet]] <- factor(p4$data[[input$abund_facet]], levels = input$abund_order)
    } else if ("Sample" %in% colnames(p4$data)) {
      p4$data$Sample <- factor(p4$data$Sample, levels = input$abund_order)
    }
  }

  if (input$flip_abundance) {
    p4 <- p4 + coord_flip()
  }

  if (input$abund_plot_type == "heat" && input$flip_abundance) {
  # Identify the variable used as x-axis (was y-axis before flipping)
  # In most microeco heatmaps, y = taxa, x = sample (or vice versa depending on transposition)
  axis_var <- names(p4$data)[1]  # or explicitly set e.g., "Sample" or "Taxon"
  
  # Drop unused levels in that axis_var
p4$data$Sample <- factor(p4$data$Sample)
p4$data$Sample <- droplevels(p4$data$Sample)
}

  p4 +
    theme(
      axis.text = element_text(size = input$beta_label_size),
      axis.title = element_text(size = input$beta_label_size),
      legend.text = element_text(size = input$beta_label_size),
      legend.title = element_text(size = input$beta_label_size),
      strip.text = element_text(size = input$beta_label_size)
    )
})

  output$abundancePlotly <- renderPlotly({
    p <- abundance_plot_obj()
    ggplotly(p)
  })
  
  output$abundancePlotStatic <- renderPlot({
    abundance_plot_obj()
  })
  
  output$alpha_group_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("alpha_group", "Group by:",
                choices = c("None", categorical_cols),
                selected = "None")
  })

  output$alpha_order_selector <- renderUI({
    req(final_physeq(), input$alpha_group)
    if (!is.null(input$alpha_group) && input$alpha_group != "None") {
      df <- as.data.frame(sample_data(final_physeq()))
      choices <- unique(df[[input$alpha_group]])
      selectizeInput("alpha_order", "Custom Order (drag to reorder):", 
                     choices = choices, selected = choices, multiple = TRUE, 
                     options = list(plugins = list('drag_drop')))
    } else {
      choices <- sample_names(final_physeq())
      selectizeInput("alpha_order", "Custom Sample Order (drag to reorder):", 
                     choices = choices, selected = choices, multiple = TRUE, 
                     options = list(plugins = list('drag_drop')))
    }
  })
  
  output$alpha_colour_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("alpha_colour", "Colour", choices = c("None", categorical_cols),
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
    ) + geom_point(size = 7)

    p + geom_point(size=5, alpha=0.7)
    
    # Apply color only if valid
    if (input$alpha_colour != "None") {
      p <- p + aes_string(color = input$alpha_colour)
      p + geom_point(size=7, alpha=0.7)
    }
    
    # Reorder if custom order given
    if (!is.null(input$alpha_order) && length(input$alpha_order) > 0) {
      if (x_var != "samples" && x_var %in% colnames(p$data)) {
        p$data[[x_var]] <- factor(p$data[[x_var]], levels = input$alpha_order)
      } else if ("samples" %in% colnames(p$data)) {
        p$data$samples <- factor(p$data$samples, levels = input$alpha_order)
      }
    }
    
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
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("beta_color", "Color by:", choices = c("None", categorical_cols), selected = "None")
  })
  
  output$beta_shape_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("beta_shape", "Shape by:", choices = c("None", categorical_cols), selected = "None")
  })
  
  output$beta_label_selector <- renderUI({
    req(final_physeq())
    cols <- colnames(sample_data(final_physeq()))
    selectInput("beta_label", "Label points by:", choices = c("None", cols), selected = "None")
  })
  
  output$beta_facet_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("beta_facet", "Facet by:", choices = c("None", categorical_cols), selected = "None")
  })

  output$beta_facet_order_selector <- renderUI({
    req(final_physeq(), input$beta_facet)
    if (!is.null(input$beta_facet) && input$beta_facet != "None") {
      df <- as.data.frame(sample_data(final_physeq()))
      choices <- unique(df[[input$beta_facet]])
      selectizeInput("beta_facet_order", "Custom Facet Order (drag to reorder):", 
                     choices = choices, selected = choices, multiple = TRUE, 
                     options = list(plugins = list('drag_drop')))
    } else {
      return(NULL)
    }
  })
  
  
  
  output$betaPlot <- renderPlot({
    req(final_physeq())
    
    dist <- distance(final_physeq(), method = input$beta_dist)
    ord <- ordinate(final_physeq(), method = input$beta_ord, distance = dist)
    
    p <- plot_ordination(final_physeq(), ord, 
                         color = if (input$beta_color != "None") input$beta_color else NULL, 
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
      if (!is.null(input$beta_facet_order) && length(input$beta_facet_order) > 0) {
        p$data[[input$beta_facet]] <- factor(p$data[[input$beta_facet]], levels = input$beta_facet_order)
      }
      p <- p + facet_wrap(as.formula(paste("~", input$beta_facet)))
    }
    
    p
  })
  
  output$permanova_group_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("permanova_group", "Select grouping variable:",
                choices = categorical_cols,
                selected = categorical_cols[1] %||% "None")
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
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("dend_treatment", "Select metadata column for grouping:",
                choices = categorical_cols,
                selected = categorical_cols[1] %||% "None")
  })
  
  output$dendrogramPlot <- renderPlot({
    req(final_physeq(), input$dend_method, input$dend_treatment, input$dend_label_size)
    
    dend <- dendrogram_phyloseq_custom(
      phyloseq_obj = final_physeq(),
      treatment = input$dend_treatment,
      method = input$dend_method,
      label_size = input$dend_label_size
    )
    
    dend <- dend + theme(
      axis.text = element_text(size = input$dend_text_size),
      axis.title = element_text(size = input$dend_text_size)
    )
    
    dend
  })
  
  output$tsne_group_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selected_val <- if (length(categorical_cols) >= 2) categorical_cols[2] else categorical_cols[1] %||% "None"
    selectInput("tsne_group", "Group samples by:",
                choices = categorical_cols,
                selected = selected_val)
  })
  
  
  output$tsne_perplexity_selector <- renderUI({
    req(final_physeq())
    
    n_samples <- nsamples(final_physeq())
    max_perplexity <- floor((n_samples - 1) / 3)
    max_perplexity <- max(5, min(max_perplexity, 50))  # reasonable bounds
    
    sliderInput("tsne_perplexity", "Perplexity:",
                min = 0, max = max_perplexity, value = 2)
  })
  
  output$tsne_label_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("tsne_label", "Label samples by:",
                choices = c("None", categorical_cols),
                selected = "None")
  })
  
  output$tsne_plot <- renderPlotly({
    req(input$run_tsne)
    req(final_physeq(), input$tsne_group)
    
    p <- tsne_phyloseq(
      phyloseq_obj = final_physeq(),
      treatment = input$tsne_group,
      perplexity = input$tsne_perplexity,
      circle = input$tsne_circle,
      labels = if (input$tsne_label != "None") input$tsne_label else NULL,
      colors = "default"
    )
    ggplotly(p)
  })
  
  observeEvent(input$reset_tsne, {
    show_tsne(FALSE)
  })
  
  output$numeric_column_selector_ui <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    selectInput("numeric_cols", "Select numeric metadata columns:",
                choices = numeric_cols, multiple = TRUE)
  })
  
  observeEvent(input$create_transenv, {
    req(final_physeq(), input$numeric_cols)
    
    dataset <- phyloseq2meco(final_physeq())
    dataset$tidy_dataset()
    dataset$cal_abund()
    dataset$cal_alphadiv()
    dataset$cal_betadiv()
    
    # Get the full metadata data frame
    meta_df <- as.data.frame(sample_data(final_physeq()))
    
    # Find the numeric indices of the selected column names
    selected_indices <- which(names(meta_df) %in% input$numeric_cols)
    
    # Create the trans_env object using the indices
    env_obj <- trans_env$new(dataset = dataset, env_cols = selected_indices)
    reactiveValues_envfit$transenv <- env_obj
    
    output$transenv_display <- renderPrint({
      env_obj
    })
  })
  
  output$visualization_sidebar <- renderUI({
    req(reactiveValues_envfit$transenv)
    
    sample_meta <- as.data.frame(sample_data(final_physeq()))
    factor_vars <- names(sample_meta)[sapply(sample_meta, function(x) is.factor(x) || is.character(x))]
    
    tagList(
      h4("Select Analysis"),
      div(
        style = "margin-top: 10px; margin-bottom: 10px;",
        fluidRow(
          column(4, div(style = "text-align: center;", actionButton("run_rda", "RDA", class = "btn-primary"))),
          column(4, div(style = "text-align: center;", actionButton("run_corr", "Correlation", class = "btn-secondary"))),
          column(4, div(style = "text-align: center;", actionButton("run_mantel", "Mantel", class = "btn-warning")))
        )
      ),
      
      
      # Only show RDA controls when RDA is selected
      conditionalPanel(
        condition = "input.run_rda % 2 == 1",
        wellPanel(
          checkboxInput("adjust_arrow_length", "Adjust Arrow Length", TRUE),
          sliderInput("max_perc_env", "Max Percentage of Explained Env Fit (arrows)", min = 0.05, max = 1, value = 0.3, step = 0.05),
          selectInput("rda_color", "Color by:", choices = factor_vars, selected = factor_vars[1] %||% "None"),
          selectInput("rda_shape", "Point Shape", choices = factor_vars, selected = factor_vars[1] %||% "None"),
          selectInput("rda_label", "Sample Labels:", choices = c("None", names(sample_meta)), selected = "None"),
          sliderInput("rda_textsize", "Text Size", value = 6, min = 6, max = 15, step = 1)
        )
      ),
      conditionalPanel(
        condition = "input.run_corr % 2 == 1",
        wellPanel(
          selectInput("input_method", "Correlation Method", choices = c("pearson", "spearman", "kendall"), selected = "spearman"),
          numericInput("input_threshold", "Abundance Threshold", value = 0.001, min = 0.00001, max = 0.9, step = 0.001),
          selectInput("p_type", "P-value Adjustment Type", choices = c("All", "Taxa", "Env"), selected = "All"),
          selectInput("group_by", "Group By", choices = c("None", names(sample_meta)), selected = "None"),
          sliderInput("Cortext_size", "Text Size", value = 10, min = 7, max = 20, step = 2),
          sliderInput("xtextangle", "Text angle",  value = 0,min=0,max = 90, step = 10)
        )
      ),
      conditionalPanel(
        condition = "output.analysis_mode == 'mantel'",
        wellPanel(
          selectInput("mantel_group", "Group By:", choices = factor_vars, selected = factor_vars[1] %||% "None"),
          actionButton("run_mantel_analysis", "Run Mantel Test")
        )
      )
    )
  })
  
  
  output$analysis_mode <- reactive({
    selected_analysis()
  })
  outputOptions(output, "analysis_mode", suspendWhenHidden = FALSE)
  
  output$visualization_main <- renderUI({
    req(selected_analysis())
    req(reactiveValues_envfit$transenv)
    req(selected_analysis())
    
    if (selected_analysis() == "rda") {
      plotOutput("rda_plot", height = "90vh")
    } else if (selected_analysis() == "corr") {
      plotOutput("corr_plot", height = "90vh")
    } else if (selected_analysis() == "mantel") {
      plotOutput("mantel_plot", height = "90vh")
    }
    
  })
  
  
  output$rda_plot <- renderPlot({
    req(selected_analysis())
    req(reactiveValues_envfit$transenv)
    
    e1 <- reactiveValues_envfit$transenv
    
    e1$cal_ordination(
      method = "RDA",
      feature_sel = FALSE
    )
    
    e1$trans_ordination(
      adjust_arrow_length = input$adjust_arrow_length,
      max_perc_env = input$max_perc_env
    )
    
    label_input <- if (input$rda_label == "None") NULL else input$rda_label
    shape_input <- if (input$rda_shape == "None") NULL else input$rda_shape
    

    
    
     e1$plot_ordination(
      plot_color = input$rda_color,
      plot_shape = shape_input,
      env_text_size = input$rda_textsize,
      taxa_text_size = input$rda_textsize,
      add_sample_label = label_input,
      point_size = input$rda_textsize
    )

 
  })
  
  output$corr_plot <- renderPlot({
    req(reactiveValues_envfit$transenv)
    
    e1 <- reactiveValues_envfit$transenv
    group_by <- if (input$group_by == "None") NULL else input$group_by
    
    e1$cal_cor(
      method = input$input_method,
      add_abund_table = NULL,
      filter_thres = input$input_threshold,
      p_adjust_method = "fdr",
      p_adjust_type = input$p_type,
      by_group = group_by,
      group_use = NULL
    )
    
    corr_plot <- e1$plot_cor(
      xtext_angle = input$xtextangle,
      xtext_size = input$Cortext_size,
      ytext_size = input$Cortext_size
    )
    print(corr_plot)
  })
  
  observeEvent(input$run_mantel_analysis, {
    req(final_physeq(), input$numeric_cols)
    
    dataset <- phyloseq2meco(final_physeq())
    dataset$tidy_dataset()
    
    # Get group variable
    group_col <- input$mantel_group
    group_vals <- unique(sample_data(final_physeq())[[group_col]])
    
    # Get the full metadata data frame
    meta_df <- as.data.frame(sample_data(final_physeq()))
    
    # Find the numeric indices of the selected column names
    selected_indices <- which(names(meta_df) %in% input$numeric_cols)
    
    if (length(group_vals) < 2) {
      showNotification("You need at least two groups to compare.", type = "error")
      return(NULL)
    }
    
    # Pre-allocate results
    mantel_tables <- list()
    mantel_objs <- list()
    
    for (g in group_vals) {
      d <- clone(dataset)
      d$sample_table <- d$sample_table[d$sample_table[[group_col]] == g, ]
      d$tidy_dataset()
      d$cal_betadiv()
      
      t <- trans_env$new(dataset = d, env_cols = selected_indices)
      t$cal_mantel(use_measure = "bray", partial_mantel = TRUE)
      
      x <- data.frame(spec = g, t$res_mantel)[, c(1, 3, 6, 8)]
      colnames(x) <- c("spec", "env", "r", "p.value")
      x <- x %>% mutate(
        rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
        pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
      )
      
      mantel_tables[[g]] <- x
      mantel_objs[[g]] <- t
    }
    
    combined_table <- do.call(rbind, mantel_tables)
    

    output$mantel_plot <- renderPlot({
      req(mantel_objs[[1]])
      
      MantCorr.Sn <- quickcor(mantel_objs[[1]]$data_env, type = "upper", cor.test = TRUE, show.diag = TRUE) +
        geom_square() +scale_fill_distiller(palette = "RdBu", direction = 1)+
        #geom_mark(sig.thres = 0.05, color = "black", size = 0) +
        anno_link(aes(colour = pd, size = rd), data = combined_table) +
        scale_size_manual(values = c(0.5, 1.5, 3)) +
        scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
        guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
               colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
               fill = guide_colorbar(title = "Pearson's r", order = 3))+
        theme(
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
      
      MantCorr.Sn
    })
  })
  
  observe({
    req(final_physeq())
    updateSelectInput(session, "env_var",
                      choices = names(sample_data(final_physeq())))
  })
  

  # Step 1: Dynamic taxonomic rank selector
  output$tax_rank_selector_regression <- renderUI({
    req(final_physeq())
    ranks <- colnames(as.data.frame(tax_table(final_physeq())))
    selectInput("tax_rank_regression", "Select Taxonomic Rank:", 
                choices = ranks, selected = tail(ranks, 1))
  })
  
  # Step 2: Taxa list updates based on selected rank
  output$taxa_selector_regression <- renderUI({
    req(final_physeq(), input$tax_rank_regression)
    
    dataset <- phyloseq2meco(final_physeq())
    dataset$tidy_dataset()
    dataset$cal_abund()
    
    rank_table <- dataset$taxa_abund[[input$tax_rank_regression]]
    validate(need(!is.null(rank_table), "Abundance table for selected rank not found."))
    
    taxa_names <- rownames(rank_table)
    selectInput("selected_taxon", "Select Taxon (Lineage):",
                choices = taxa_names, selected = taxa_names[1])
  })
  
  
  output$regression_group_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("group", "Group by:", choices = c("None", categorical_cols), selected = "None")
  })
  
  
  
  output$regression_plot <- renderPlot({
    req(input$run_scatter)
    req(input$tax_rank_regression, input$selected_taxon, input$env_var)
    req(final_physeq())
    
    # Prepare dataset and env object
    dataset <- phyloseq2meco(final_physeq())
    dataset$tidy_dataset()
    dataset$cal_abund()
    
    # Create trans_env object if not already stored
    env_obj <- trans_env$new(dataset = dataset, env_cols = input$env_var)
    
    # Get abundance vector for selected taxon
    rank_table <- dataset$taxa_abund[[input$tax_rank_regression]]
    validate(need(input$selected_taxon %in% rownames(rank_table), "Selected taxon not found."))
    lineage_vector <- as.numeric(rank_table[input$selected_taxon, ])
    
    group_val <- if (input$group == "None") NULL else input$group
    
    # Generate plot using plot_scatterfit()
    # Extract the last taxon name from input$selected_taxon (format: K_Taxa1|p_Taxa2|...)
    taxon_parts <- unlist(strsplit(input$selected_taxon, "|", fixed = TRUE))
    last_taxon <- tail(taxon_parts, 1)
    # Remove prefix like x__ from the last taxon (e.g., c__ from c__Bacilli)
    last_taxon <- sub("^[a-zA-Z]__", "", last_taxon)
    
    p <- env_obj$plot_scatterfit(
      x = lineage_vector,
      y = input$env_var,
      group = group_val,
      type = "lm",
      point_size = input$point_size,
      point_alpha = 1,
      line_color = "#2A0E3C",
      line_se_color = "#A87CA0",
      label.x.npc = "left", label.y.npc = "top",
      x_axis_title = last_taxon,
      y_axis_title = input$env_var
    ) + theme_classic() + 
    theme(
      axis.text = element_text(size = input$text_size),
      axis.title = element_text(size = input$text_size),
      legend.text = element_text(size = input$text_size),
      legend.title = element_text(size = input$text_size),
      text = element_text(size = input$text_size)  # This should cover annotations like equations
    )
    
    print(p)
  })
  
  # Indicator Species Analysis Server Logic
  source("src/shap-phyloseq.R", local = TRUE)

  output$indicator_variable_selector <- renderUI({
    req(final_physeq())
    df <- as.data.frame(sample_data(final_physeq()))
    categorical_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    selectInput("indicator_var", "Select Metadata Variable:", choices = categorical_cols)
  })

  output$indicator_group1_selector <- renderUI({
    req(final_physeq(), input$indicator_var)
    choices <- unique(as.data.frame(sample_data(final_physeq()))[[input$indicator_var]])
    selectInput("indicator_group1", "Select Group 1 (will be coded as 1):", choices = choices, selected = choices[1])
  })

  output$indicator_group2_selector <- renderUI({
    req(final_physeq(), input$indicator_var)
    choices <- unique(as.data.frame(sample_data(final_physeq()))[[input$indicator_var]])
    selectInput("indicator_group2", "Select Group 2 (will be coded as 0):", choices = choices, selected = choices[2])
  })

  indicator_results <- eventReactive(input$run_indicator_analysis, {
    req(final_physeq(), input$indicator_var, input$indicator_group1, input$indicator_group2)
    
    ps <- final_physeq()
    
    # Run the SHAP analysis function from the sourced script
    results <- run_shap_analysis(
      phyloseq_obj = ps,
      variable = input$indicator_var,
      group1 = input$indicator_group1,
      top_n = input$top_n_taxa,
      font_size = input$indicator_font_size
    )
    
    return(results)
  })

  output$indicator_plot <- renderPlot({
    results <- indicator_results()
    req(results$plot)
    results$plot
  })

  output$indicator_table <- renderDataTable({
    results <- indicator_results()
    req(results$table)
    results$table
  })
  
}
