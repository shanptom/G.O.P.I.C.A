library(shiny)
library(shinyBS)
library(shinyjs)
library(shinycssloaders)
library(bslib)
# library(plotly) # Commented out as we're using static ggplot2 plots

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "favicon.png"),
    tags$script(HTML("
      $(document).ready(function(){
        $('[data-toggle=\"tooltip\"]').tooltip();
      });
    "))
  ),
  navbarPage("goPICa",
                 id = "main_nav",
                 theme = bs_theme(
                   bootswatch = "morph", # https://bootswatch.com/
                   primary = "#006699",
                   base_font = font_google("Inter")
                 ),
                 tabPanel("Home",
                      fluidRow(
                        # Column for Logo (left side)
                        column(
                          width = 2,
                          div(style = "padding-top: 30px; text-align: center;",
                              img(src = "logo.png", height = "120px", style = "max-width: 100%;")
                          )
                        ),
                        # Column for Text (right side)
                        column(
                          width = 10,
                          h2("Graphical Operations Platform for Interactive Community Analysis"),
                          p("GOPICA allows you to explore microbial community data using various visualizations and analyses."),
                          tags$ul(
                            tags$li(" Start by uploading your ASV, taxonomy, and metadata tables (.csv) or phyloseq object (.rds) under 'Upload Data'."),
                            tags$li(" Remove unwanted taxa and rarefy at the filtering tab."),
                            tags$li(" Explore abundance, diversity, and dendrograms in the respective tabs."),
                            tags$li(" Customize plots with the sidebar controls."),
                            tags$li(" All plots support dynamic interaction based on your metadata."),
                            tags$li(" No coding experience required — just upload your files and explore!")
                          ),
                          br(),
                          
                          h3("Contact"),
                          p("For questions, contact the ",
                            tags$a(href = "https://shanptom.github.io", target = "_blank", "developer.")
                          ),
                          
                          h3("References"),
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
                tabPanel("User Manual",
                      fluidRow(
                        column(10, offset = 1,
                               includeMarkdown("docs/user_guide.md")
                        )
                      )
             ),
                 tabPanel("Upload Data",
                          sidebarLayout(
                            sidebarPanel( width = 3,
                              fileInput("asv", "Upload Count Table (CSV)", accept = ".csv"),
                              fileInput("tax", "Upload Taxonomy Table (CSV)", accept = ".csv"),
                              fileInput("meta", "Upload Metadata Table (CSV)", accept = ".csv"),
                              fileInput("phylo", "Or Upload Phyloseq Object (.rds)", accept = ".rds"),
                              hr(),
                              h4("Load Demo Data"),
                              selectInput("demo_file", "Select Demo Dataset:", 
                                          choices = c("Phyloseq RDS (demo_ps.rds)" = "rds", 
                                                      "CSV Set (demo_asv.csv, demo_meta.csv, demo_tax.csv)" = "csv")),
                              actionButton("load_demo", "Load Demo Data")
                            ),
                            mainPanel(width = 9,
                                      uiOutput("upload_status_ui")
                                      )
                          )
                 ),
                 tabPanel("Filter", value = "filter_tab",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                              checkboxInput("doRarefy", "Apply rarefaction", value = FALSE),
                              checkboxInput("doTSS", "Normalize by TSS", value = FALSE),
                              uiOutput("taxa_filters"),
                              actionButton("apply_filter", "Apply Filtering"),
                              actionButton("go_analysis", "Go to Analysis")
                            ),
                            mainPanel(width = 9,verbatimTextOutput("filter_status"))
                          )
                 ),
                 tabPanel("Rarefaction Plot", value = "rarefaction_tab",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                              uiOutput("rarefaction_color_selector"),
                              uiOutput("rarefaction_facet_selector"),
                              sliderInput("beta_label_size", "Axis Text Size:", min = 6, max = 20, value = 12),
                              sliderInput("rarefaction_label_size", "Sample Label Size:", min = 2, max = 10, value = 4),
                              checkboxInput("show_rarefaction_labels", "Show Sample Labels", value = TRUE),
                              div(
                                style = "display: flex; align-items: center;",
                                uiOutput("rarefaction_facet_order_selector"),
                                tags$span(
                                  tags$i(class = "fas fa-info-circle", style = "margin-left: 5px; cursor: pointer;"),
                                  `data-toggle` = "tooltip",
                                  `data-placement` = "right",
                                  `title` = "Drag to reorder facets as they appear on the plot."
                                )
                              )
                            ),
                            mainPanel(width = 9,withSpinner(plotOutput("rarefactionPlot",height = "770px", width = "100%")))
                          )
                 ),
                 tabPanel("Abundance", value = "abundance_tab",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          radioButtons("abund_plot_type", "Plot Type:",
                                       choices = c("Bar" = "bar", "Line" = "line", "Heatmap" = "heat"),
                                       selected = "bar"),
                          div(
                            style = "display: flex; align-items: center;",
                            uiOutput("tax_rank_selector"),
                            tags$span(
                              tags$i(class = "fas fa-info-circle", style = "margin-left: 5px; cursor: pointer;"),
                              `data-toggle` = "tooltip",
                              `data-placement` = "right",
                              `title` = "Select the taxonomic level for abundance analysis (e.g., Phylum, Genus)."
                            )
                          ),
                          sliderInput("ntaxa", "Number of Top Taxa:", min = 5, max = 15, value = 8, step = 1),
                          uiOutput("abundance_facet_selector"),
                          
                          div(
                            style = "display: flex; align-items: center;",
                            uiOutput("abundance_order_selector"),
                            tags$span(
                              tags$i(class = "fas fa-info-circle", style = "margin-left: 5px; cursor: pointer;"),
                              `data-toggle` = "tooltip",
                              `data-placement` = "right",
                              `title` = "Drag to reorder facets or samples as they appear on the plot."
                            )
                          ),
                          
                          sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12),
                          checkboxInput("flip_abundance", "Flip axes (horizontal plot)", value = FALSE)
                        ),
                        mainPanel(width = 9, withSpinner(plotOutput("abundancePlotly", height = "770px", width = "100%")))
                      )
                 ),

                 tabPanel("Alpha Diversity", value = "alpha_tab",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                              checkboxGroupInput("alpha_index", "Select Diversity Index:",
                                                 choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"),
                                                 selected = c("Shannon")),
                              uiOutput("alpha_group_selector"),
                              uiOutput("alpha_colour_selector"),
                              checkboxInput("flip_alpha", "Flip axes (horizontal plot)", value = FALSE),
                              div(
                                style = "display: flex; align-items: center;",
                                uiOutput("alpha_order_selector"),
                                tags$span(
                                  tags$i(class = "fas fa-info-circle", style = "margin-left: 5px; cursor: pointer;"),
                                  `data-toggle` = "tooltip",
                                  `data-placement` = "right",
                                  `title` = "Drag to reorder groups or samples as they appear on the plot."
                                )
                              ),
                              sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12)
                            ),
                            mainPanel(width = 9,withSpinner(plotOutput("alphaPlot",height = "770px", width = "100%")))
                          )
                 ),
                   tabPanel("Dendrogram", value = "dendrogram_tab",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     uiOutput("dend_treatment_selector"),
                                     selectInput("dend_method", "Select distance method:",
                                                 choices = c("euclidian", "manhattan", "canberra", "clark", "bray",
                                                             "kulczynski", "jaccard", "gower", "altGower", "morisita",
                                                             "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"),
                                                 selected = "bray"),
                                     sliderInput("dend_label_size", "Label Size:", min = 3, max = 10, value = 5, step = 1),
                                     sliderInput("dend_text_size", "Text Size:", min = 6, max = 20, value = 12)
                        ),
                        
                        mainPanel(width = 9,withSpinner(plotOutput("dendrogramPlot", height = "770px", width = "100%"))))
             ),
                tabPanel("Ordination", value = "ordination_tab",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          bsCollapse(id = "ordination_collapse_panel",
                            bsCollapsePanel("Ordination Settings",
                              selectInput("beta_dist", "Distance Method:",
                                          choices = c("bray", "unifrac", "wunifrac", "jaccard", "dpcoa", "jsd", "manhattan",
                                                      "euclidean", "canberra", "binomial"),
                                          selected = "bray"),
                              selectInput("beta_ord", "Ordination Method:",
                                          choices = c("NMDS", "MDS", "PCoA", "DCA", "CCA", "RDA", "DPCoA"),
                                          selected = "NMDS"),
                              style = "info"
                            ),
                            bsCollapsePanel("Aesthetics",
                              uiOutput("beta_color_selector"),
                              uiOutput("beta_shape_selector"),
                              uiOutput("beta_label_selector"),
                              uiOutput("beta_facet_selector"),
                              uiOutput("beta_facet_order_selector"),
                              sliderInput("beta_label_size", "Axis Text Size:", min = 6, max = 20, value = 12),
                              sliderInput("beta_label_text_size", "Label Text Size:", min = 2, max = 15, value = 3),
                              sliderInput("beta_shape_size", "Shape Size:", min = 1, max = 10, value = 4),
                              style = "info"
                            ),
                            bsCollapsePanel("PERMANOVA",
                              uiOutput("permanova_group_selector"),
                              actionButton("run_permanova", "Run PERMANOVA"),
                              verbatimTextOutput("permanova_result"),
                              style = "info"
                            ),
                            bsCollapsePanel("t-SNE Analysis",
                              uiOutput("tsne_group_selector"),
                              uiOutput("tsne_perplexity_selector"),
                              checkboxInput("tsne_circle", "Draw circles", value = FALSE),
                              uiOutput("tsne_label_selector"),
                              actionButton("run_tsne", "Run tSNE"),
                              conditionalPanel(
                                condition = "output.show_tsne_flag",
                                actionButton("reset_tsne", "Back to Ordination")
                              ),
                              style = "info"
                            )
                          )
                        ),
                        mainPanel(width = 9,
                                  conditionalPanel(
                                    condition = "!output.show_tsne_flag",
                                    withSpinner(plotOutput("betaPlot", height = "770px", width = "100%"))
                                  ),
                                  conditionalPanel(
                                    condition = "output.show_tsne_flag",
                                    withSpinner(plotOutput("tsne_plot", height = "770px", width = "100%"))
                                  )
                        
                        )
                      )
             ),
             tabPanel("Metadata Analysis", value = "metadata_tab",
                      fluidRow(
                        column(
                          width = 3,
                          # Initial controls (trans_env creation)
                          uiOutput("numeric_column_selector_ui"),
                          actionButton("create_transenv", "Create trans_env Object"),
                          br(), br(),
                          verbatimTextOutput("transenv_display"),
                          uiOutput("continue_button_ui"),
                          
                          # Collapsible RDA controls (only visible after "Continue")
                          uiOutput("visualization_sidebar")
                        ),
                        column(
                          width = 9,
                          withSpinner(uiOutput("visualization_main"))  # RDA plot output
                        )
                      )
                      
             ),
             tabPanel("Regression", value = "regression_tab",
                      fluidRow(
                        column(
                          width = 3,
                          uiOutput("tax_rank_selector_regression"),
                          uiOutput("taxa_selector_regression"),
                          selectInput("env_var", "Select Environmental Variable:",
                                      choices = NULL),
                          uiOutput("regression_group_selector"),
                          actionButton("run_scatter", "Run Scatter Plot"),
                          sliderInput("point_size", "Point Size:", min = 1, max = 6, value = 3),
                          sliderInput("text_size", "Text Size:", min = 6, max = 20, value = 12),
                          
                        ),
                        column(
                          width = 9,
                          withSpinner(plotOutput("regression_plot", height = "770px"))
                        )
                      )
             ),
             tabPanel("Indicator Species", value = "indicator_tab",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          h4("SHAP Analysis for Indicator Species"),
                          uiOutput("indicator_variable_selector"),
                          uiOutput("indicator_group1_selector"),
                          uiOutput("indicator_group2_selector"),
                          sliderInput("top_n_taxa", "Number of Top Taxa to Display:", min = 5, max = 30, value = 10, step = 1),
                          sliderInput("indicator_font_size", "Font Size:", min = 6, max = 20, value = 10, step = 1),
                          actionButton("run_indicator_analysis", "Run Analysis")
                        ),
                        mainPanel(width = 9,
                          withSpinner(plotOutput("indicator_plot", height = "770px")),
                          dataTableOutput("indicator_table")
                        )
                      )
             )
             
             
  ))
