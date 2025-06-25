# User Guide

**MetaPiX** is an interactive Shiny app for exploring microbiome data through visual analytics and statistical tools. This guide walks you through each major feature and explains the expected input formats.

---

##  1. Upload Data

You can upload your data in one of two ways:

a.  **CSV Files**:

  * **ASV Table**: Rows = ASVs, Columns = Samples
  * **Taxonomy Table**: Rows = ASVs, Columns = Taxonomic ranks
  * **Metadata Table**: Rows = Samples  
    

b.  **`phyloseq` Object** (`.rds`): A ready-to-use object containing all necessary components.  


> ⚠️ *To use phylogeny-based tools, your `phyloseq` object must include a phylogenetic tree.* Follow [this guide](https://github.com/shanptom/metaX/blob/main/Phyloseq.md) to create `phyloseq` object from csv and to construct phylogenetic trees.

---

##  2. Filter

* Remove unwanted taxa by name at any taxonomic level (e.g., *Chloroplast*, *Mitochondria*)
* Apply **rarefaction** (optional) to normalize sequencing depth
* Click **Apply Filtering** to save and continue
>  After completing the **Filter** step, you can freely explore other analysis tabs. The filtered object is stored and reused across all downstream modules.

---

##  3. Rarefaction Plot

* Visualize sequencing depth across all samples
* Customize using metadata-based **color** and **facet** controls

---

##  4. Abundance

Visualize taxonomic profiles across your samples:

*  Bar plots
*  Line plots
*  Heatmaps

Customize with:

* Taxonomic rank (e.g., Genus, Family)
* Top N most abundant taxa
* Sample order and facet settings
* Axis flip for stratigraphic or temporal data

---

##  5. Dendrogram

* Hierarchical clustering using various distance metrics
* Group samples using metadata
* Customize label size and alignment

---

##  6. Alpha Diversity

* Select diversity indices (e.g., Shannon, Simpson)
* Group samples by metadata
* Adjust ordering, axis flip, color, and label size

---

##  7. Beta Diversity & Ordination

* Choose ordination methods: **NMDS**, **PCoA**, **tSNE**, etc.
* Select distance measures: **Bray-Curtis**, **Jaccard**, etc.
* Customize aesthetics: color, shape, labels, facets
* Perform **PERMANOVA** to test group-level differences

---

##  8. Metadata Analysis

Explore relationships between microbial communities and **numerical metadata** (e.g., pH, temperature, age, BMI):

Steps:

1. Select numeric columns by their start and end positions
2. Click **"Create trans\_env"**
3. Choose one of the following analyses:

   * **RDA**: Redundancy Analysis
   * **Correlation**: Genus vs. metadata
   * **Mantel Test**: Distance-based correlation

> ⚠️ *Ensure all selected columns are numeric. Factor/text data will cause errors.*

> Powered by the **`microeco`** package.

---

##  9. Regression

* Select a **taxon** (e.g., Genus)
* Select a **numeric environmental variable**
* Optionally group by metadata (e.g., site, treatment)
* Generate linear model visualizations

---

##  Tips

* Hover over inputs for helpful tooltips
* Use comma-separated values for manual ordering of samples
* All plots dynamically update based on your selections
* For Large Files: The Shiny app may become slow or unresponsive. To improve performance, run the app locally by following these steps:
    1. Download or clone the repository from [GitHub](https://github.com/shanptom/metaX).
    2. Open the MetaPiX.R file in RStudio and run the app
* A set of demo files are also available on Github for reference 

---
## Running MetaPiX Locally

To run MetaPiX locally, follow these steps:


### Get MetaPix repo from GiHub

You can download the MetaPiX repository from GitHub using the following command. Make sure to replace `~/Path/to/your/folder` with the actual path where you want to clone the repository.
```
git clone https://github.com/shanptom/MetaPiX.git ~/Path/to/your/folder
```
When you run the app for first time open the MetaPiX.R file and install all the necessary packages

If you are using a HPC cluster, you can use the following commands to load the necessary modules and run the app:(The below commands are for the HPC cluster at the Ohio Supercomputer Center, please modify them according to your HPC cluster)

```bash
ml gcc/12.3.0
ml R/4.4.0
ml gdal/3.7.3
ml proj/9.2.1
ml geos/3.12.0
R
```
If you are using a local machine, you can run the app directly in RStudio or R console and can ignore the module loading commands above.

```r
library(shiny)
runApp('/Path/to/MetaPiX.R')
```