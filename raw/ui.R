library(shiny)
library(shinyBS)
library(shinyjs)
library(shinycssloaders)
library(bslib)
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

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "favicon.png")
  ),
  navbarPage("MetaPix",
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
                          h2("Welcome to MetaPix"),
                          p("This application allows you to explore microbial community data using various visualizations and analyses."),
                          tags$ul(
                            tags$li("💾 Start by uploading your ASV, taxonomy, and metadata tables (.csv) or phyloseq object (.rds) under 'Upload Data'."),
                            tags$li("🧪 Remove unwanted taxa and rarefy at the filtering tab."),
                            tags$li("📈 Explore abundance, diversity, and dendrograms in the respective tabs."),
                            tags$li("🎯 Customize plots with the sidebar controls."),
                            tags$li("🧬 All plots support dynamic interaction based on your metadata."),
                            tags$li("👨‍💻 No coding experience required — just upload your files and explore!")
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
                tabPanel("User Guide",
                      fluidRow(
                        column(10, offset = 1,
                               includeMarkdown("user_guide.md")
                        )
                      )
             ),
                 tabPanel("Upload Data",
                          sidebarLayout(
                            sidebarPanel( width = 3,
                              fileInput("asv", "Upload Count Table (CSV)", accept = ".csv"),
                              fileInput("tax", "Upload Taxonomy Table (CSV)", accept = ".csv"),
                              fileInput("meta", "Upload Metadata Table (CSV)", accept = ".csv"),
                              fileInput("phylo", "Or Upload Phyloseq Object (.rds)", accept = ".rds")
                            ),
                            mainPanel(width = 9,verbatimTextOutput("upload_status"),
                                      )
                          )
                 ),
                 tabPanel("Filter",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                              checkboxInput("skip_filter", "Skip Filtering"),
                              checkboxInput("doRarefy", "Apply rarefaction", value = FALSE),
                              checkboxInput("doTSS", "Normalize by TSS", value = FALSE),
                              uiOutput("taxa_filters"),
                              actionButton("apply_filter", "Apply Filtering"),
                              actionButton("go_analysis", "Go to Analysis")
                            ),
                            mainPanel(width = 9,verbatimTextOutput("filter_status"))
                          )
                 ),
                 tabPanel("Rarefaction Plot",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                              uiOutput("rarefaction_color_selector"),
                              uiOutput("rarefaction_facet_selector"),
                              sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12)
                            ),
                            mainPanel(width = 9,plotOutput("rarefactionPlot",height = "770px", width = "100%"))
                          )
                 ),
                 tabPanel("Abundance",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          radioButtons("abund_plot_type", "Plot Type:",
                                       choices = c("Bar" = "bar", "Line" = "line", "Heatmap" = "heat"),
                                       selected = "bar"),
                          uiOutput("tax_rank_selector"),
                          sliderInput("ntaxa", "Number of Top Taxa:", min = 5, max = 15, value = 8, step = 1),
                          uiOutput("abundance_facet_selector"),
                          
                          div(
                            style = "display: flex; align-items: center;",
                            tags$label(
                              "Custom Order (comma-separated) ℹ ",
                              `title` = "Specify sample names in the order you want them to appear on the plot. Example: Sample3, Sample2, Sample1"
                            ),
                
                          ),
                          textInput("abund_order", label = NULL, value = ""),
                          
                          sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12),
                          checkboxInput("flip_abundance", "Flip axes (horizontal plot)", value = FALSE)
                        ),
                        mainPanel(width = 9,plotOutput("abundancePlot", height = "770px", width = "100%"))
                      )
                 ),

                 tabPanel("Alpha Diversity",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                              checkboxGroupInput("alpha_index", "Select Diversity Index:",
                                                 choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"),
                                                 selected = c("Shannon")),
                              uiOutput("alpha_group_selector"),
                              uiOutput("alpha_colour_selector"),
                              checkboxInput("flip_alpha", "Flip axes (horizontal plot)", value = FALSE),
                              textInput("alpha_order", "Custom order (comma-separated values) ", value = ""),
                              sliderInput("beta_label_size", "Text Label Size:", min = 6, max = 20, value = 12)
                            ),
                            mainPanel(width = 9,plotOutput("alphaPlot",height = "770px", width = "100%"))
                          )
                 ),
                   tabPanel("Dendrogram",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     uiOutput("dend_treatment_selector"),
                                     selectInput("dend_method", "Select distance method:",
                                                 choices = c("euclidian", "manhattan", "canberra", "clark", "bray",
                                                             "kulczynski", "jaccard", "gower", "altGower", "morisita",
                                                             "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"),
                                                 selected = "bray"),
                                     sliderInput("dend_label_size", "Label size:", min = 3, max = 10, value = 5, step = 1)
                        ),
                        
                        mainPanel(width = 9,plotOutput("dendrogramPlot", height = "770px", width = "100%")))
             ),
                tabPanel("Ordination",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          selectInput("beta_dist", "Distance Method:",
                                      choices = c("bray", "unifrac", "wunifrac", "jaccard", "dpcoa", "jsd", "manhattan",
                                                  "euclidean", "canberra", "binomial"),
                                      selected = "bray"),
                          selectInput("beta_ord", "Ordination Method:",
                                      choices = c("NMDS", "MDS", "PCoA", "DCA", "CCA", "RDA", "DPCoA"),
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
                          verbatimTextOutput("permanova_result"),
                          tags$hr(),
                          h4("tSNE Analysis"),
                          uiOutput("tsne_group_selector"),
                          uiOutput("tsne_perplexity_selector"),
                          checkboxInput("tsne_circle", "Draw circles", value = FALSE),
                          uiOutput("tsne_label_selector"),
                          actionButton("run_tsne", "Run tSNE"),
                          conditionalPanel(
                            condition = "output.show_tsne_flag",
                            actionButton("reset_tsne", "Back to Ordination")
                          )
                          
                        ),
                        mainPanel(width = 9,
                                  conditionalPanel(
                                    condition = "!output.show_tsne_flag",
                                    plotOutput("betaPlot", height = "770px", width = "100%")
                                  ),
                                  conditionalPanel(
                                    condition = "output.show_tsne_flag",
                                    plotOutput("tsne_plot", height = "770px", width = "100%")
                                  )
                        
                        )
                      )
             ),
             tabPanel("Metadata Analysis",
                      fluidRow(
                        column(
                          width = 3,
                          # Initial controls (trans_env creation)
                          numericInput("colx", "Enter the column number where numerical data starts:", value = 1, min = 1),
                          numericInput("coly", "Enter the column number where numerical data ends:", value = 1, min = 1),
                          actionButton("create_transenv", "Create trans_env Object"),
                          br(), br(),
                          verbatimTextOutput("transenv_display"),
                          uiOutput("continue_button_ui"),
                          
                          # Collapsible RDA controls (only visible after "Continue")
                          uiOutput("visualization_sidebar")
                        ),
                        column(
                          width = 9,
                          uiOutput("visualization_main")  # RDA plot output
                        )
                      )
                      
             ),
             tabPanel("Regression",
                      fluidRow(
                        column(
                          width = 3,
                          uiOutput("tax_rank_selector_regression"),
                          uiOutput("taxa_selector_regression"),
                          selectInput("env_var", "Select Environmental Variable:",
                                      choices = NULL),
                          uiOutput("regression_group_selector"),
                          actionButton("run_scatter", "Run Scatter Plot"),
                          sliderInput("point_alpha", "Point Alpha:", min = 0.1, max = 1, value = 0.5, step = 0.1),
                          sliderInput("point_size", "Point Size:", min = 1, max = 6, value = 3),
                          sliderInput("regression_text_size", "Text Size:", min = 6, max = 20, value = 12),
                          
                        ),
                        column(
                          width = 9,
                          plotOutput("regression_plot", height = "770px")
                        )
                      )
             )
             
             
  ))
