# Single-Cell RNA-Seq QC and Clustering

**scRNA-seq quality control, normalization, and cell type identification**

**Stack:** R, Seurat

## Overview

This project implements comprehensive single-cell RNA-sequencing analysis workflows including quality control, normalization, dimensionality reduction, clustering, and cell type identification. It demonstrates best practices for scRNA-seq data processing and interpretation using bone marrow mononuclear cell (BMMC) and CD34+ samples.

## Problem & Approach

**Problem:** Process raw scRNA-seq data to identify distinct cell populations and their characteristics in bone marrow samples.

**Approach:**
- Quality control filtering (genes, cells, mitochondrial content)
- Data normalization and scaling
- Highly variable gene identification
- PCA and dimensionality reduction
- Clustering (Louvain/Leiden algorithms)
- Cell type annotation with reference datasets
- Doublet detection and removal
- Batch effect correction for multi-sample integration
- Differential expression analysis
- Gene ontology enrichment
- Trajectory inference for lymphoid lineage

## Tech Stack

- **R** (≥4.0)
- **Seurat** - scRNA-seq analysis framework
- **DoubletFinder** - Doublet detection
- **SingleR** - Cell type annotation
- **monocle3** - Trajectory analysis
- **enrichR** - Gene ontology enrichment
- **CellChat** - Cell-cell communication
- Additional QC and visualization packages (ggplot2, patchwork, dplyr)

## Data Requirements

This analysis requires single-cell RNA-seq count matrices in RDS format. The script expects data in a `data/` subdirectory containing:

- `GSM4138872_scRNA_BMMC_D1T1.rds` - BMMC sample, Donor 1, Technical replicate 1
- `GSM4138873_scRNA_BMMC_D1T2.rds` - BMMC sample, Donor 1, Technical replicate 2
- `GSM4138874_scRNA_CD34_D2T1.rds` - CD34+ sample, Donor 2, Technical replicate 1
- `GSM4138875_scRNA_CD34_D3T1.rds` - CD34+ sample, Donor 3, Technical replicate 1

**To obtain the data:**
1. Download the samples from GEO (Gene Expression Omnibus) if publicly available
2. Place the RDS files in a `data/` folder in the repository root
3. Ensure the RDS files contain raw count matrices compatible with Seurat's `CreateSeuratObject()`

**Alternative:** If you have your own scRNA-seq data, modify lines 348-357 in `scrna_analysis.r` to point to your data files and update sample metadata accordingly.

## Setup / Installation

### Install Required R Packages

```r
# Install CRAN packages
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "tidyverse"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment"))

# Install GitHub packages
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
remotes::install_github("cole-trapnell-lab/monocle3")
remotes::install_github("satijalab/seurat-wrappers")
```

### Alternative: Using renv (Recommended)

For reproducible package management:

```r
# Install renv if not already installed
install.packages("renv")

# Restore the project environment (if renv.lock exists)
renv::restore()
```

## How to Run

```bash
# Run from repository root
Rscript scrna_analysis.r
```

**Note:** The script performs intensive computations including:
- Quality control across multiple filtering strategies
- PCA, UMAP, and t-SNE dimensionality reduction
- Doublet detection with parameter sweeps
- Batch correction and integration
- Cell type annotation
- Differential expression analysis
- Trajectory inference

Execution time will vary based on system resources (expect 30-60+ minutes on a typical workstation).

## Outputs

The script generates numerous intermediate and final outputs:

### RDS Objects (saved at various stages)
- `list_sample_read.rds` - Initial loaded samples
- `sample_info_table.rds` - Sample metadata summary
- `outlier_then_threshold.rds` - QC-filtered samples
- `samples_seurats_list_qc_finished.rds` - Post-QC and doublet removal
- `merged_samples_uncorrected.rds` - Merged without batch correction
- `merged_samples_corrected.rds` - Integrated data with batch correction
- `final_annotated.rds` - Final annotated dataset
- Various enrichment and annotation results

### Visualizations (PNG files)
- QC violin plots (pre/post filtering)
- Correlation plots (nCount vs nFeature)
- PCA, UMAP, and t-SNE plots
- Feature plots for marker genes
- Elbow plots
- Cell type proportion plots
- Volcano plots for differential expression
- Enrichment analysis plots

### Analysis Results
- Cell type annotations (manual and automated)
- Differentially expressed genes
- GO enrichment results
- Trajectory analysis outputs

## Reproducibility Notes

- **Paths:** Uses `getwd()` for portability - script works from the repository root
- **Random seed:** Set at line 3: `set.seed(42)` for reproducible clustering
- **Parallel processing:** Configurable via `detectCores() - 1` (line 24)
- **Data provenance:** Original data from GEO (bone marrow samples)

## Key Analysis Steps

1. **Data Loading & Initial QC** (lines 337-425)
   - Load 4 samples (2 BMMC, 2 CD34+)
   - Calculate QC metrics (mt%, hb%, ribo%, ERCC%)
   
2. **Quality Control Strategies** (lines 427-573)
   - Compare multiple filtering approaches
   - Apply outlier detection (MAD-based)
   - Threshold-based filtering
   - Visualize before/after QC

3. **Sample Processing** (lines 575-586)
   - Normalization, variable feature selection, scaling
   - PCA with automatic dimensionality selection
   - Clustering and doublet removal
   - UMAP and t-SNE visualization

4. **Data Integration** (lines 587-636)
   - Merge samples with/without batch correction
   - Compare uncorrected vs corrected integration
   - Visualize batch effects

5. **Cell Type Annotation** (lines 639-794)
   - Automated annotation with SingleR
   - Manual annotation with marker genes
   - Cell type proportion analysis

6. **Differential Expression** (lines 859-969)
   - B-cells vs T-cells
   - T-cells vs Monocytes
   - BMMC vs CD34+ samples
   - Volcano plots and top DEGs

7. **Functional Enrichment** (lines 1065-1125)
   - GO Biological Process
   - GO Cellular Component
   - GO Molecular Function

8. **Trajectory Analysis** (lines 1127-1200)
   - Lymphoid lineage trajectory
   - Pseudotime ordering
   - Root cell identification (HSCs)

## Notes

- Comprehensive scRNA-seq QC pipeline with multiple filtering strategies
- Automated doublet detection with parameter optimization
- Reference-based cell type annotation
- Multiple clustering resolutions explored
- Batch effect correction demonstrated
- Rich visualization outputs for publication-quality figures
- Modular function design for reusability
