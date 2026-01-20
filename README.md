# Divergent mechanisms of active antibiotic resistome enrichment in soil driven by pesticide diversity

This repository provides the workflow (Shell, Python, and R) for analyzing the active soil resistome.

## 🛠️ Data Processing
* **`main_pipeline.sh`**: A comprehensive SLURM-ready shell script for raw data processing (fastp, MEGAHIT, Prodigal, Diamond).
* **`summarize_metagenome.py`**: A Python utility to calculate **TPM** and generate the final Master Matrix.

## 📊 Statistical Analysis (R)

### 1. Community & Diversity
* **`01_Beta_Diversity_PCoA.R`**: PCoA and PERMANOVA for community shifts.
* **`08_Procrustes_PCoA.R`**: Congruence analysis between ARG and MGE profiles.

### 2. Composition & Profiling
* **`02_Feature_Abundance.R`**: Stacked bar charts with T-tests.
* **`03_Target_Overall_ARG.R`**: Violin plots for specific gene abundance comparison.
* **`04_Abundance_Heatmap.R`**: Log10-transformed heatmaps with significance stars.

### 3. Drivers & Mechanisms
* **`05_Pesticide_Trend_Clustering.R`**: **Mfuzz** clustering for gradient response patterns.
* **`06_WGCNA_Genus_vs_ARG.R`**: **WGCNA** for microbial module-trait correlations.
* **`07_Mantel_Correlation.R`**: **Mantel tests** for ARG, MGE, and environmental links.

### 4. Interaction & Response
* **`09_Nested_ANOVA_BCa.R`**: **Nested ANOVA** with **BCa Bootstrap CI** for interaction effects.
* **`09_Transfer_ROS_Analysis.R`**: Analysis of ROS and horizontal transfer potential.

---
*For detailed parameters and usage, please refer to the headers of each script.*
