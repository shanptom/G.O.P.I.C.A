# MetaPiX User Guide

**MetaPiX** is an interactive Shiny app designed for exploring microbiome data through powerful visual analytics and statistical tools. This user guide provides a comprehensive walkthrough of each major feature, along with detailed information on expected input formats to ensure a smooth experience.

---

## 1. Upload Data

MetaPiX allows you to upload your microbiome data in one of two formats for analysis:

- **CSV Files**: Upload separate tables for different data components.
  - **Count Table**: Rows represent ASVs (Taxa), and columns represent Samples.
  - **Taxonomy Table**: Rows represent ASVs (Taxa), and columns represent taxonomic ranks.
  - **Metadata Table**: Rows represent Samples, with columns for associated metadata.

- **`phyloseq` Object** (`.rds`): A pre-assembled R object that integrates all necessary components (count table, taxonomy, metadata, and optionally a phylogenetic tree) into a single file for direct use in MetaPiX.

> **⚠️ Note**: To utilize phylogeny-based tools, your `phyloseq` object must include a phylogenetic tree. For instructions on creating a `phyloseq` object from CSV files and constructing phylogenetic trees, refer to [this guide](https://github.com/shanptom/metaX/blob/main/Phyloseq.md).

---

## 2. Filter Data

Refine your dataset by removing noise and normalizing data for consistent analysis:

- **Remove Unwanted Taxa**: Filter out specific taxa by name at any taxonomic level (e.g., *Chloroplast*, *Mitochondria*).
- **Apply Rarefaction** (Optional): Normalize sequencing depth across samples to reduce bias.
- **Save Changes**: Click **Apply Filtering** to update your dataset and proceed.

> **Note**: Once the **Filter** step is complete, the filtered dataset is stored and automatically applied to all downstream analysis tabs, allowing seamless exploration across modules.

---

## 3. Rarefaction Plot

Assess the sequencing depth of your samples to ensure adequate coverage:

- **Visualize Depth**: Generate plots to see sequencing depth distribution across all samples.
- **Customize Display**: Use metadata to control **color** coding and **facet** organization for clearer insights.

---

## 4. Abundance Analysis

Explore the taxonomic composition of your samples with various visualization options:

- **Plot Types**:
  - **Bar Plots**: Show relative abundance stacked by taxon.
  - **Line Plots**: Track abundance trends across samples or conditions.
  - **Heatmaps**: Display abundance patterns as a color-coded matrix.

- **Customization Options**:
  - **Taxonomic Rank**: Select the level of classification (e.g., Genus, Family).
  - **Top N Taxa**: Focus on the most abundant taxa by specifying a number.
  - **Sample Order & Facets**: Arrange samples manually or by metadata categories.
  - **Axis Flip**: Adjust orientation for stratigraphic or temporal data representation.

---

## 5. Dendrogram

Analyze sample relationships through hierarchical clustering:

- **Clustering Methods**: Choose from various distance metrics for grouping.
- **Group by Metadata**: Organize samples based on metadata categories.
- **Customize Appearance**: Adjust label size and alignment for readability.

---

## 6. Alpha Diversity

Measure within-sample diversity to understand community richness and evenness:

- **Diversity Indices**: Choose metrics like Shannon, Simpson, or others.
- **Group by Metadata**: Compare diversity across different sample categories.
- **Customize Visualization**: Adjust sample ordering, axis orientation, color schemes, and label sizes.

---

## 7. Beta Diversity & Ordination

Examine between-sample diversity to uncover community structure and patterns:

- **Ordination Methods**: Select from options like **NMDS**, **PCoA**, **tSNE**, and more.
- **Distance Measures**: Choose metrics such as **Bray-Curtis**, **Jaccard**, or others.
- **Visual Customization**: Adjust colors, shapes, labels, and facet layouts.
- **Statistical Testing**: Use **PERMANOVA** to assess significant differences between groups.

---

## 8. Metadata Analysis

Investigate how microbial communities correlate with **numerical metadata** (e.g., pH, temperature, age, BMI):

**Steps to Analyze**:
1. **Select Numeric Columns**: Choose columns by specifying their start and end positions in the metadata table.
2. **Create Environment Dataset**: Click **"Create trans_env"** to prepare data for analysis.
3. **Choose Analysis Type**:
   - **RDA (Redundancy Analysis)**: Assess how much variation in community structure is explained by metadata.
   - **Correlation**: Examine direct relationships between genera and metadata variables.
   - **Mantel Test**: Evaluate correlations based on distance matrices.

> **⚠️ Warning**: Ensure all selected columns contain numeric data. Including factor or text data will result in errors.

> **Note**: This module is powered by the **`microeco`** package for robust statistical analysis.

---

## 9. Regression Analysis

Model relationships between microbial taxa and environmental variables:

- **Select Taxon**: Choose a specific taxon (e.g., a Genus) for analysis.
- **Select Variable**: Pick a numeric environmental variable to test against.
- **Group by Metadata** (Optional): Stratify analysis by categories like site or treatment.
- **Visualize Results**: Generate plots based on linear models to interpret relationships.

---

## Tips for Using MetaPiX

Maximize your experience with these helpful pointers:

- **Tooltips**: Hover over input fields and buttons for contextual help and explanations.
- **Manual Ordering**: Use comma-separated values to customize the order of samples in visualizations.
- **Dynamic Updates**: All plots automatically refresh based on your input selections.
- **Handling Large Files**: The online Shiny app may slow down with large datasets. For better performance, consider running MetaPiX locally:
  1. Download or clone the repository from [GitHub](https://github.com/shanptom/metaX).
  2. Open the `MetaPiX.R` file in RStudio and launch the app.
- **Demo Files**: Access sample datasets on GitHub for testing and learning purposes.

---
## Running MetaPiX Locally

For optimal performance or to handle large datasets, you can run MetaPiX on your local machine. Follow these steps to set up and launch the application:

### Step 1: Clone the MetaPiX Repository

Download the MetaPiX source code from GitHub. Replace `~/Path/to/your/folder` with the desired directory path on your system:

```bash
git clone https://github.com/shanptom/MetaPiX.git ~/Path/to/your/folder
```

### Step 2: Install Dependencies

When launching MetaPiX for the first time, install the required R packages by executing this command in your R environment. Adjust the path to point to the `install_dep.R` script in your cloned repository:

```r
source("Path/to/install_dep.R")
```

### Step 3: Environment Setup for HPC Clusters (Optional)

If you're using a High-Performance Computing (HPC) cluster, load the necessary modules before running R. The commands below are tailored for the Ohio Supercomputer Center (OSC). Modify them based on your specific HPC environment:

```bash
ml gcc/12.3.0
ml R/4.4.0
ml gdal/3.7.3
ml proj/9.2.1
ml geos/3.12.0
R
```

> **Note**: If you're on a local machine, skip this step and proceed to launching the app.

### Step 4: Launch the Application

Run MetaPiX directly from RStudio or an R console. Replace `/Path/to/MetaPiX.R` with the actual path to the `MetaPiX.R` file in your cloned repository:

```r
library(shiny)
runApp('/Path/to/MetaPiX.R')
```
