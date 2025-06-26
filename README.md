# `MetaPiX`: Interactive Microbial Community Analysis in Shiny

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![R Shiny](https://img.shields.io/badge/built%20with-R%20Shiny-blue)](https://shiny.rstudio.com/)
[![Status](https://img.shields.io/badge/status-active-brightgreen)]()

`MetaPiX` is an interactive R Shiny application for microbial community data exploration. Built for researchers analyzing high-throughput sequencing data (e.g., 16S/18S rRNA metabarcoding), `MetaPiX` offers a no-code interface to perform filtering, diversity analysis, regression, ordination, and reporting — all in your browser.

---

## 🚀 Features

- **Upload Data**: Load ASV, taxonomy, and metadata tables or a pre-built phyloseq `.rds` object.
- **Filtering**: Remove taxa, apply rarefaction, or normalize using TSS.
- **Alpha & Beta Diversity**: Plot Shannon, Simpson, NMDS, PCoA, PERMANOVA, etc.
- **Abundance Analysis**: Bar plots, heatmaps, and alluvial plots by taxonomic rank.
- **Dendrograms**: Visualize hierarchical clustering of samples.
- **Rarefaction**: Assess sequencing depth across samples.
- **Metadata Analysis**: Perform constrained ordination (RDA/CCA/dbRDA), correlations, and Mantel tests.
- **Regression**: Relate taxa abundance to environmental variables using `microeco` models.

---

## 📂 Input Requirements

You can upload:

### Option 1: CSVs
- **ASV Table**: rows = taxa, columns = samples.
- **Taxonomy Table**: rows = taxa, columns = ranks (Kingdom → Species).
- **Metadata Table**: rows = samples, columns = metadata variables.

### Option 2: `.rds` file
- A valid `phyloseq` object created using the `phyloseq` package.

---

## 🔧 Normalization Options

- **None**: Use raw counts.
- **Rarefaction**: Subsample each sample to a common read depth.
- **TSS**: Total Sum Scaling normalization.

---



## 📦 Dependencies

`MetaPiX` relies on the following R packages:

- `shiny`, `shinyjs`, `shinyBS`, `shinycssloaders`
- `phyloseq`, `microeco`, `phylosmith`
- `ggplot2`, `RColorBrewer`, `ggcor`
- `vegan`, `ranacapa`, `microbiome`
- `rmarkdown`, `knitr`

---

## 🧠 How It Works

`MetaPiX` wraps around powerful microbiome R tools using a modular backend:
- Uses `microeco`'s object-oriented classes for ordination, correlation, regression, and visualization.
- Converts all CSVs into a `phyloseq` object, which is passed through your workflow.
- Stores user inputs and plots reactively.

---

## 👨‍💻 Contact



For feedback or feature suggestions, feel free to open an issue or contact the developer.

---

## 📜 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## 📚 References

If you use `MetaPiX` in a publication, please cite the underlying R packages that power the app:

- McMurdie & Holmes (2013). *phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data*. PLoS ONE.
- Liu et al. (2021). *microeco: an R package for data mining in microbial community ecology*. bioRxiv.
- Dixon (2003). *VEGAN, a package of R functions for community ecology*. Journal of Vegetation Science.
- Wickham (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer.

We also acknowledge the use of:
- `shiny`, `shinyjs`, `ggcor`, `RColorBrewer`, `phylosmith`, `ranacapa`,  and other CRAN/Bioconductor packages associated with the above ones.


