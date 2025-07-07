
cran_pkgs <- c("shiny", "shinyBS", "shinyjs", "shinycssloaders", "bslib",
               "ggplot2", "vegan", "microeco", "file2meco", "GUniFrac",
               "ragg",
               "RColorBrewer", "ggalluvial", "dplyr", "ggpubr", "xgboost", "SHAPforxgboost")


installed <- rownames(installed.packages())
to_install <- setdiff(cran_pkgs, installed)
if (length(to_install)) install.packages(to_install)


if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
github_pkgs <- c("gauravsk/ranacapa", "schuyler-smith/phylosmith", "joey711/phyloseq")
lapply(github_pkgs, devtools::install_github)