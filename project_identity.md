# Project Identity: Single-Cell RNA-Seq QC and Clustering

## Professional Project Identity

### Display Title
Single-Cell RNA-Seq QC and Clustering with Seurat

### Repo Slug (Recommended)
`scrna-qc-clustering-annotation-seurat`

### Tagline
scRNA-seq quality control, normalization, and cell type identification

### GitHub Description
Comprehensive single-cell RNA-sequencing analysis pipeline implementing quality control, normalization, dimensionality reduction, clustering, and cell type identification for bone marrow samples using Seurat and related tools.

### Primary Stack
- R (≥4.0)
- Seurat
- DoubletFinder
- SingleR
- monocle3
- enrichR
- CellChat

### Topics/Keywords
- scrna-seq
- single-cell-rna-seq
- seurat
- bioinformatics
- cell-type-annotation
- quality-control
- clustering
- doublet-detection
- batch-correction
- trajectory-analysis
- differential-expression
- r

### What Problem It Solves
Single-cell RNA sequencing generates complex datasets requiring rigorous quality control, normalization, and analysis to identify distinct cell populations and their characteristics. This project provides a comprehensive, reproducible pipeline for processing raw scRNA-seq data from bone marrow samples (BMMC and CD34+) to identify cell types, detect doublets, correct batch effects, and perform downstream analyses including differential expression, gene ontology enrichment, and trajectory inference.

### Approach
The pipeline implements industry-standard best practices:
1. Multi-strategy quality control with outlier detection
2. Seurat-based normalization and variable feature selection
3. PCA with automatic dimensionality selection
4. DoubletFinder for doublet removal
5. Batch effect correction using Seurat integration
6. Automated and manual cell type annotation with SingleR
7. Comprehensive differential expression analysis
8. Gene ontology enrichment analysis
9. Trajectory inference for lineage reconstruction

### Inputs
- Single-cell RNA-seq count matrices in RDS format
- Expected 4 samples: 2 BMMC replicates, 2 CD34+ samples
- Sample metadata (donor, replicate, sex)

### Outputs
- Quality-controlled and annotated Seurat objects
- Dimensionality reduction visualizations (PCA, UMAP, t-SNE)
- Cell type annotations and proportions
- Differentially expressed genes between cell types and conditions
- Gene ontology enrichment results
- Trajectory analysis for lymphoid lineage
- Publication-ready PNG figures throughout
