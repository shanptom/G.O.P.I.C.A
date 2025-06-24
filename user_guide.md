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
* Use comma-separated values for manual ordering of samples.
* All plots dynamically update based on your selections

---

 You can return to this tab anytime for guidance.
