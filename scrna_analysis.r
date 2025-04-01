rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])
# initializing
set.seed(42)

 # detectCores()


suppressPackageStartupMessages({
  library(dplyr)
  library(spatstat.core)
  library(Seurat)
  library(patchwork)
  library(DoubletFinder)
  library(SingleR)
  library(enrichR)
  library(CellChat)
  library(SingleCellExperiment)
  library(tidyverse)
  library(celldex)
  library(parallel)
  library(monocle3)
  library(SeuratWrappers)
})
options(mc.cores = detectCores() - 1)




work_sample_seurat <- function(test_sample,name = "NA",doublet= T,Normal=T,Feature=T,Scale=T) {
  if (name == "NA"){
  name <- test_sample$orig.ident[1]
  }
  # normalize data finding feature variables and scaling the data
  
  if (Normal)
  test_sample <- test_sample %>% NormalizeData()
  if (Feature)
  test_sample <- test_sample %>% FindVariableFeatures() 
  if (Scale)
  test_sample <- test_sample %>% ScaleData()
  
  test_sample <- RunPCA(test_sample)
  
  stdv <- test_sample[["pca"]]@stdev
  sum.stdv <- sum(test_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc

  test_sample <- RunPCA(test_sample, npcs = min.pc)
  pc<-DimPlot(test_sample)
  ggsave(paste0(name, "_pca.png"), pc, width = 20, height = 20)
  
  
  dim_loadings <- test_sample[["pca"]]@feature.loadings
  
  # Extract top genes for dimensions 1 and 2
  top_genes_dim1 <- head(sort(dim_loadings[, 1], decreasing = TRUE), n = 10)  # Top 10 for Dim 1
  top_genes_dim2 <- head(sort(dim_loadings[, 2], decreasing = TRUE), n = 10)  # Top 10 for Dim 2
  
  
  gn <- VizDimLoadings(test_sample,dims = 1:2)
  ggsave(paste0(name, "_VizDimLoadings.png"), gn, width = 20, height = 20)
  
  fp <- FeaturePlot(test_sample, features = names(top_genes_dim1))
  ggsave(paste0(name, "_FeaturePlot_dim1.png"), fp, width = 20, height = 20)
  
  fb2 <- FeaturePlot(test_sample, features = names(top_genes_dim2))
  ggsave(paste0(name, "_FeaturePlot_dim2.png"), fb2, width = 20, height = 20)
  el <- ElbowPlot(test_sample,ndims = 50)
  ggsave(paste0(name, "_ElbowPlot.png"), el, width = 20, height = 20)
  
  test_sample <- FindNeighbors(object = test_sample, dims = 1:min.pc)              
  test_sample <- FindClusters(object = test_sample, resolution = 0.1)
  test_sample <- RunUMAP(test_sample, dims = 1:min.pc)
  
  umap_plot <-DimPlot(test_sample, reduction = "umap")
  ggsave(paste0(name, "_umap.png"), umap_plot, width = 20, height = 20)
  
  test_sample <- RunTSNE(test_sample, dims = 1:min.pc)
  tsne_plot <-DimPlot(test_sample, reduction = "tsne")
  ggsave(paste0(name, "_tsne.png"), tsne_plot, width = 20, height = 20)
  
  if (doublet){

  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(test_sample, PCs = 1:min.pc)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- test_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(test_sample@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  test_sample <- doubletFinder(seu = test_sample,
                               pN = 0.25,
                               PCs = 1:min.pc, 
                               pK = optimal.pk,
                               nExp = nExp.poi.adj,
                               reuse.pANN = F,
                               sct = F)
  
  
  metadata <- test_sample@meta.data
  colnames(metadata)[length(colnames(metadata))] <- "doublet_finder"
  test_sample@meta.data <- metadata 
  
  
  
  # subset and save
  test_sample <- subset(test_sample, doublet_finder == "Singlet")
  }
  return(test_sample)
}

add_metrics <- function(pbmc) {
  #   Mitochondrial Genes (percent.mt)
  #      High Values: Suggest stressed or dying cells; may result from cell lysis.
  #      Low Values: Indicate healthy cells with intact membranes.
  #   Hemoglobin Genes (percent.hb)
  #      High Values:
  #          Expected in erythroid cells.
  #          Unexpected in other cell types; may indicate contamination.
  #      Low Values: Suggest minimal contamination from blood cells.
  #   Ribosomal Genes (percent.ribo)
  #      High Values:
  #          May reflect technical artifacts or high protein synthesis activity.
  #          Could be associated with certain cell cycle phases.
  #      Low to Moderate Values: Typically acceptable and expected in most cells.
  #   ERCC Spike-Ins (percent.ERCC)
  #      High Values:
  #          May indicate issues with RNA input amounts or library prep.
  #          Could suggest low endogenous RNA content.
  #      Low Values: Reflect appropriate levels of spike-in controls.
  
  
  # percentage of hemoglobin genes
  # reason : hemoglobin genes are expressed in blood cells (filter no blood contamination)
  # but in our case, sicne in bone marrow, erythroid cells are an essential component. High percent.hb is expected and represents genuine biological signal. so we wont consider this metrics in future process
  pbmc[["percent.hb"]] <- PercentageFeatureSet(pbmc, pattern = "^HB[^(P)]")
  # percentage of mitochondrial genes
  # reason : mitochondrial genes are expressed in dying cells (filter no dying cells)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  # percentage of ERCC spike-in genes
  # reason : ERCC spike-in genes are used to estimate technical noise (filter no technical noise)
  # sinc in our case we dont have ERCC spike-in genes, we can remove this metric in future
  pbmc[["percent.errc"]] <- PercentageFeatureSet(pbmc, pattern = "^ERCC-")
  # percentage of ribosomal genes
  # reason : ribosomal genes are expressed in all cells
  # in thbis context (bone marrow), Actively dividing cells or progenitor cells might have higher ribosomal content. but since we assume that all cells are in the same stage of cell cycle, we can remove this metric in future
  pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RPS|^RPL|^MRPS|^MRPL")
  return(pbmc)
}

add_metrics_all <- function(list_sample_read) {
  for (i in 1:length(list_sample_read)) {
    list_sample_read[[i]] <- add_metrics(list_sample_read[[i]])
  }
  return(list_sample_read)
}

# Remove outlier MAD approach:
is_outlier <- function(scdata, metric, nmads) {
  M <- scdata@meta.data[, metric]
  outlier <- (M < median(M) - nmads * mad(M)) | (median(M) + nmads * mad(M) < M)
  return(outlier)
}

add_is_outlier_all <- function(list_sample_read) {
  for (i in 1:length(list_sample_read)) {
    list_sample_read[[i]]@meta.data$outlier  <- (is_outlier(list_sample_read[[i]], "nCount_RNA", 5) | is_outlier(list_sample_read[[i]], "nFeature_RNA", 5))
  }
  return(list_sample_read)
}

plot_violin <- function(pbmc, features) {
  all_plots <- list()
  for (feat in features) {
    print(feat)
    vln_plot <- VlnPlot(pbmc, features = feat, pt.size = 0.1, ncol = 1) +
      NoLegend() + labs(x = pbmc$orig.ident) +
      theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(3, 1, 1, 1), "cm"))
    
    vln_plot <- vln_plot + scale_y_continuous(breaks = scales::pretty_breaks(n = 20))
    all_plots[[feat]] <- vln_plot
  }
  
  wrapped_plots <- wrap_plots(plotlist = all_plots, ncol = length(features)) +
    plot_annotation(title = paste("QC metrics", pbmc$orig.ident)) &
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  return(wrapped_plots)
}


plot_violin_all <- function(list_sample_read, qc_features) {
  violin_plots_pre_filter <- c()
  for (i in 1:length(list_sample_read)) {
    violin_plots_pre_filter[[i]] <- plot_violin(list_sample_read[[i]], qc_features)
  }
  violin_plots_pre_filter_wraped <- wrap_plots(violin_plots_pre_filter, ncol = 2)
  # return both wraped and singles
  return(list(violin_plots_pre_filter_wraped, violin_plots_pre_filter))
}

# corelation plot function
plot_corelate <- function(seurat_obj) {
  # Extract metadata
  data <- seurat_obj@meta.data
  
  # Calculate the correlation coefficient
  cor_value <- cor(data$nCount_RNA, data$nFeature_RNA)
  
  # Create the scatter plot
  p <- ggplot(data, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(alpha = 0.3, color = "blue") +
    labs(
      title = paste("Correlation for sample : ", seurat_obj$orig.ident),
      subtitle = paste("Pearson correlation:", round(cor_value, 2)),
      x = "UMI Counts (nCount_RNA)",
      y = "Number of Features (nFeature_RNA)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    ) +
    geom_smooth(method = "lm", color = "red", se = FALSE)
  
  # Return the plot
  return(p)
}

plot_corelate_all <- function(list_sample_read) {
  cor_plots_pre_filter <- c()
  
  for (i in 1:length(list_sample_read)) {
    cor_plots_pre_filter[[i]] <- plot_corelate(list_sample_read[[i]])
  }
  cor_plots_pre_filter_wraped <- wrap_plots(cor_plots_pre_filter, ncol = 2)
  # return both wraped and singles
  return(list(cor_plots_pre_filter_wraped, cor_plots_pre_filter))
}

filter_seurat_object_thresholds <- function(seurat_obj, thresholds) {
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > thresholds$nFeature_RNA_bottom &
      nFeature_RNA < thresholds$nFeature_RNA_top &
      nCount_RNA > thresholds$nCount_RNA_bottom &
      nCount_RNA < thresholds$nCount_RNA_top &
      percent.mt < thresholds$percent_mt
  )
  return(seurat_obj)
}

# Function to filter a list of Seurat objects based on individual threshold values
filter_seurat_list_thresholds <- function(seurat_list, thresholds_list) {
  # Apply the filtering function to each object in the list
  for (i in seq_along(seurat_list)) {
    seurat_list[[i]] <- filter_seurat_object_thresholds(seurat_list[[i]], thresholds_list[[i]])
  }
  return(seurat_list)
}
filter_seurat_outliers <- function(seurat_obj) {
  non_outlier_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$outlier == FALSE]
  
  # Subset the Seurat object by specifying the cells to keep
  seurat_obj <- subset(seurat_obj, cells = non_outlier_cells)
  return(seurat_obj)
}

filter_seurat_list_outliers <- function(seurat_list) {
  for (i in 1:length(seurat_list)) {
    seurat_list[[i]] <- filter_seurat_outliers(seurat_list[[i]])
  }
  return(seurat_list)
}

# Extract sample information from Seurat object
get_sample_info <- function(seurat_obj) {
  name_file <- seurat_obj$orig.ident[1]
  # Get the number of cells
  num_cells <- ncol(seurat_obj)
  
  # Get the number of genes
  num_genes <- nrow(seurat_obj)
  
  # Get the names of metadata columns
  metadata_cols <- colnames(seurat_obj@meta.data)
  
  # Create a list to store the information
  sample_info <- list(
    name_file = name_file,
    num_cells = num_cells,
    num_genes = num_genes,
    metadata_columns = metadata_cols
  )
  
  return(sample_info)
}

# load data sets
loadDataSet <- function(sample_item, i) {
  filename <- sample_item$file[i]
  raw_counts <- readRDS(file = filename)
  pbmc <- CreateSeuratObject(counts = raw_counts, project = sample_item$names[i], assay = "RNA")
  return(pbmc)
}


# add metadata
addMetaData <- function(pbmc, samples, i) {
  pbmc$orig.ident <- samples$names[i]
  pbmc$donor <- samples$donor[i]
  pbmc$replicate <- samples$replicate[i]
  pbmc$sex <- samples$sex[i]
  return(pbmc)
}



getwd()
# Use project root directory instead of hardcoded path
dir <- getwd()
relative_path_data <- "data"
data_folder <- file.path(dir, relative_path_data)


# setwd(dir) # Not needed if already in project directory


# Data loading configuration
data_frame_samples <- data.frame(
    file = c("GSM4138872_scRNA_BMMC_D1T1.rds", "GSM4138873_scRNA_BMMC_D1T2.rds", "GSM4138874_scRNA_CD34_D2T1.rds", "GSM4138875_scRNA_CD34_D3T1.rds"),
    names = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1"),
    donor = c("D1", "D1", "D2", "D3"),
    replicate = c("T1", "T2", "T1", "T1"),
    sex = c("F", "F", "M", "F")
)

# add data_folder to each file
data_frame_samples$file <- file.path(data_folder, data_frame_samples$file)


list_sample_read <- c()

for (i in 1:length(data_frame_samples$file)) {
    list_sample_read <- c(list_sample_read, loadDataSet(data_frame_samples, i))
}

list_sample_read <- lapply(c(1:4), function(i) addMetaData(list_sample_read[[i]], data_frame_samples, i))



sample_info_df <- lapply(list_sample_read, get_sample_info)
# Convert the list of results into a data.frame
sample_info_table <- do.call(rbind, lapply(sample_info_df, function(info) {
    data.frame(
        name_file = info$name_file,
        num_cells = info$num_cells,
        num_genes = info$num_genes,
        metadata_columns = paste(info$metadata_columns, collapse = ", "),
        row.names = NULL
    )
}))

# View the result (use View() in interactive R or print)
print(sample_info_table)
saveRDS(sample_info_table, file = "sample_info_table.rds")
saveRDS(list_sample_read, file = "list_sample_read.rds")

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])

list_sample_read <- readRDS(file = "list_sample_read.rds")
sample_info_table <- readRDS(file = "sample_info_table.rds")

# QC metrics

# add metrics to each sample


qc_features <- c("percent.mt", "percent.hb", "percent.ribo", "percent.errc", "nCount_RNA", "nFeature_RNA")
# since percent.errc are all zero and High expression of ribo can be a normal feature of certain cell types in our context , we can remove them
# Note: This dataset has omitted the mitochondrial genes
# also since hb genes are expressed in blood cells, and in our context, we have bone marrow, we can remove it too
qc_features <- qc_features[!qc_features %in% c("percent.ribo", "percent.errc", "percent.mt", "percent.hb")]



# original samples
orig_samples <- list_sample_read
threshold_filter_only <- list_sample_read
outlier_only <- list_sample_read
threshold_then_outlier <- list_sample_read
outlier_then_threshold <- list_sample_read

# add metrics
orig_samples <- add_metrics_all(orig_samples)
threshold_filter_only <- add_metrics_all(threshold_filter_only)
outlier_only <- add_metrics_all(outlier_only)
threshold_then_outlier <- add_metrics_all(threshold_then_outlier)
outlier_then_threshold <- add_metrics_all(outlier_then_threshold)


# plot original samples
violin_plots_pre_filter_wraped <- plot_violin_all(orig_samples, qc_features)[[1]]
violin_plots_pre_filter_wraped

cor_plots_pre_filter_wraped <- plot_corelate_all(orig_samples)[[1]]
cor_plots_pre_filter_wraped

# #these threshold values are based on the results of the corelation plot and the violin

threshold_list <- list(
  list(
    nFeature_RNA_top = 3400,
    nFeature_RNA_bottom = 700,
    nCount_RNA_top = 8500,
    nCount_RNA_bottom = 1100,
    percent_mt = 5
  ),
  list(
    nFeature_RNA_top = 3400,
    nFeature_RNA_bottom = 700,
    nCount_RNA_top = 9100,
    nCount_RNA_bottom = 1100,
    percent_mt = 5
  ),
  list(
    nFeature_RNA_top = 3500,
    nFeature_RNA_bottom = 700,
    nCount_RNA_top = 9500,
    nCount_RNA_bottom = 1100,
    percent_mt = 5
  ),
  list(
    nFeature_RNA_top = 2900,
    nFeature_RNA_bottom = 700,
    nCount_RNA_top = 7500,
    nCount_RNA_bottom = 1100,
    percent_mt = 5
  )
)

# threshold filter only
threshold_filter_only <- filter_seurat_list_thresholds(threshold_filter_only, threshold_list)

# violin plot
violin_plots_threshold_filter_only_wraped <- plot_violin_all(threshold_filter_only, qc_features)[[1]]
violin_plots_threshold_filter_only_wraped

# corelation plot
cor_plots_threshold_filter_only_wraped <- plot_corelate_all(threshold_filter_only)[[1]]
cor_plots_threshold_filter_only_wraped

# outlier filter only
outlier_only <- add_is_outlier_all(outlier_only)
outlier_only <- filter_seurat_list_outliers(outlier_only)
outlier_only

# violin plot
violin_plots_outlier_only_wraped <- plot_violin_all(outlier_only, qc_features)[[1]]
violin_plots_outlier_only_wraped

# corelation plot
cor_plots_outlier_only_wraped <- plot_corelate_all(outlier_only)[[1]]
cor_plots_outlier_only_wraped

# threshold then outlier
threshold_then_outlier <- filter_seurat_list_thresholds(threshold_then_outlier, threshold_list)
threshold_then_outlier <- add_is_outlier_all(threshold_then_outlier)
threshold_then_outlier <- filter_seurat_list_outliers(threshold_then_outlier)

# violin plot
violin_plots_threshold_then_outlier_wraped <- plot_violin_all(threshold_then_outlier, qc_features)[[1]]
violin_plots_threshold_then_outlier_wraped

# corelation plot
cor_plots_threshold_then_outlier_wraped <- plot_corelate_all(threshold_then_outlier)[[1]]
cor_plots_threshold_then_outlier_wraped

# outlier then threshold
outlier_then_threshold <- add_is_outlier_all(outlier_then_threshold)
outlier_then_threshold <- filter_seurat_list_outliers(outlier_then_threshold)

# plot violin
violin_plots_outlier_then_threshold_wraped <- plot_violin_all(outlier_then_threshold, qc_features)[[1]]
violin_plots_outlier_then_threshold_wraped

# plot corelation
cor_plots_outlier_then_threshold_wraped <- plot_corelate_all(outlier_then_threshold)[[1]]
cor_plots_outlier_then_threshold_wraped

# new threshold

threshold_list_outlier <- list(
  list(
    nFeature_RNA_top = 2800,
    nFeature_RNA_bottom = 620,
    nCount_RNA_top = 8500,
    nCount_RNA_bottom = 1100,
    percent_mt = 5
  ),
  list(
    nFeature_RNA_top = 3000,
    nFeature_RNA_bottom = 650,
    nCount_RNA_top = 9100,
    nCount_RNA_bottom = 1100,
    percent_mt = 5
  ),
  list(
    nFeature_RNA_top = 3600,
    nFeature_RNA_bottom = 600,
    nCount_RNA_top = 9900,
    nCount_RNA_bottom = 1100,
    percent_mt = 5
  ),
  list(
    nFeature_RNA_top = 3100,
    nFeature_RNA_bottom = 600,
    nCount_RNA_top = 7650,
    nCount_RNA_bottom = 1000,
    percent_mt = 5
  )
)

outlier_then_threshold <- filter_seurat_list_thresholds(outlier_then_threshold, threshold_list_outlier)

# plot violin
violin_plots_outlier_then_threshold_wraped <- plot_violin_all(outlier_then_threshold, qc_features)[[1]]
violin_plots_outlier_then_threshold_wraped

# plot corelation
cor_plots_outlier_then_threshold_wraped <- plot_corelate_all(outlier_then_threshold)[[1]]
cor_plots_outlier_then_threshold_wraped

# save all plots
ggsave("violin_plots_pre_filter_wraped.png", violin_plots_pre_filter_wraped, width = 20, height = 20)
ggsave("cor_plots_pre_filter_wraped.png", cor_plots_pre_filter_wraped, width = 20, height = 20)
ggsave("violin_plots_threshold_filter_only_wraped.png", violin_plots_threshold_filter_only_wraped, width = 20, height = 20)
ggsave("cor_plots_threshold_filter_only_wraped.png", cor_plots_threshold_filter_only_wraped, width = 20, height = 20)
ggsave("violin_plots_outlier_only_wraped.png", violin_plots_outlier_only_wraped, width = 20, height = 20)
ggsave("cor_plots_outlier_only_wraped.png", cor_plots_outlier_only_wraped, width = 20, height = 20)
ggsave("violin_plots_threshold_then_outlier_wraped.png", violin_plots_threshold_then_outlier_wraped, width = 20, height = 20)
ggsave("cor_plots_threshold_then_outlier_wraped.png", cor_plots_threshold_then_outlier_wraped, width = 20, height = 20)
ggsave("violin_plots_outlier_then_threshold_wraped.png", violin_plots_outlier_then_threshold_wraped, width = 20, height = 20)
ggsave("cor_plots_outlier_then_threshold_wraped.png", cor_plots_outlier_then_threshold_wraped, width = 20, height = 20)



# by analizing the results , the method which first we find outliers and then filter the data based on thresholds, is the best method

##first we need to remove the outliers and the filter the data based on thresholds , this is because the outliers can affect the thresholds and the results
##also we do this before merging the data, since the physical properties of the cells are different and the way that experiment has been done is different, so we need to do this for each sample separately


# save the results
saveRDS(outlier_then_threshold, file = "outlier_then_threshold.rds")

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])

samples_seurats <- readRDS(file = "outlier_then_threshold.rds")


for (i in 1:length(samples_seurats)) {
  samples_seurats[[i]] <- work_sample_seurat(samples_seurats[[i]])
}


saveRDS(samples_seurats, file = "samples_seurats_list_qc_finished.rds")

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])

samples_seurats <- readRDS(file = "samples_seurats_list_qc_finished.rds")


merged_samples_uncorrected <- merge(
  x = samples_seurats[[1]], 
  y = c(samples_seurats[[2]], samples_seurats[[3]], samples_seurats[[4]]), 
  add.cell.ids = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1"), 
  project = "MergedSamplesUnCorrected"
)

merged_samples_uncorrected <- work_sample_seurat(merged_samples_uncorrected, "MergedSamplesUnCorrected",doublet = F)
d_unc_umap <- DimPlot(merged_samples_uncorrected, reduction = "umap",group.by =c("orig.ident", "seurat_clusters"))
d_unc_tsne <- DimPlot(merged_samples_uncorrected, reduction = "tsne",group.by =c("orig.ident", "seurat_clusters"))
d_unc_pca <- DimPlot(merged_samples_uncorrected, reduction = "pca",group.by =c("orig.ident", "seurat_clusters"))

features <- SelectIntegrationFeatures(object.list = samples_seurats)
anchors <- FindIntegrationAnchors(object.list = samples_seurats)


merged_samples_corrected[["RNA"]] <-JoinLayers(merged_samples_corrected[["RNA"]])

##
merged_samples_corrected <- IntegrateData(anchorset = anchors)



DefaultAssay(merged_samples_corrected) <- "integrated"
merged_samples_corrected <- work_sample_seurat(merged_samples_corrected, "MergedSamplesCorrected",doublet = F,Scale = T,Feature = F,Normal =F)


d_cor_umap <- DimPlot(merged_samples_corrected, reduction = "umap",group.by =c("orig.ident", "seurat_clusters"))
d_cor_tsne <- DimPlot(merged_samples_corrected, reduction = "tsne",group.by =c("orig.ident", "seurat_clusters"))
d_cor_pca <- DimPlot(merged_samples_corrected, reduction = "pca",group.by =c("orig.ident", "seurat_clusters"))

wrap_umap_unc_cor = wrap_plots(list(d_unc_umap, d_cor_umap), ncol = 1)

wrap_tsne_unc_cor = wrap_plots(list(d_unc_tsne, d_cor_tsne), ncol = 1)
wrap_pca_unc_cor = wrap_plots(list(d_unc_pca, d_cor_pca), ncol = 1)

ggsave("wrap_umap_unc_cor.png", wrap_umap_unc_cor, width = 20, height = 20)
ggsave("wrap_tsne_unc_cor.png", wrap_tsne_unc_cor, width = 20, height = 20)
ggsave("wrap_pca_unc_cor.png", wrap_pca_unc_cor, width = 20, height = 20)



saveRDS(merged_samples_uncorrected, file = "merged_samples_uncorrected.rds")
saveRDS(merged_samples_corrected, file = "merged_samples_corrected.rds")

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])

###
#Yes, batch correction is necessary when merging single-cell RNA sequencing datasets from different samples. Technical variations and batch effects—such as differences in cell lysis, PCR amplification efficiency, reagent lots, and sequencing equipment—can introduce unwanted variability that confounds true biological differences. These technical factors can overshadow the biological signals in the data, making it challenging to accurately address research questions and leading to misleading results in downstream analyses.

#To mitigate these issues, computational batch correction methods like Seurat's integration are employed. Key parameters used in Seurat's integration method include selecting integration features (number of variable genes), identifying integration anchors with specified principal components (e.g., dims = 1:20), and integrating the data based on these anchors. By applying batch correction, we remove technical variations, ensuring that the data reflects true biological differences. This necessity is evident when comparing UMAP plots before and after batch correction: without correction, cells cluster by sample origin due to batch effects; after correction, cells cluster based on biological similarities, demonstrating that batch correction is essential for accurate data interpretation.

corrected <- readRDS(file = "merged_samples_corrected.rds")

cell_dex_dataset <- HumanPrimaryCellAtlasData()

sce = as.SingleCellExperiment(corrected,assay = "integrated")

# singler
predicted_main <- SingleR(test = sce, ref = cell_dex_dataset, labels = cell_dex_dataset$label.main)

predicted_fine <- SingleR(test = sce, ref = cell_dex_dataset, labels = cell_dex_dataset$label.main)

saveRDS(predicted_main, file = "predicted_main_singler.rds")

saveRDS(predicted_fine, file = "predicted_fine_singler.rds")



predicted_main <- readRDS("predicted_main_singler.rds")

predicted_fine <- readRDS("predicted_fine_singler.rds")


prediction_combination <- combineCommonResults(
  list("Broad"=predicted_main,"Fine"=predicted_fine)
)

corrected$predicted_singler <- prediction_combination$pruned.labels


table_df <- as.data.frame.matrix(table(corrected$predicted_singler, corrected$seurat_clusters))
assigned_cell_types <- vector("character", length = ncol(table_df))
names(assigned_cell_types) <- colnames(table_df)

# Assign the cell type with the highest count to each cluster
for (cluster in colnames(table_df)) {
  # Get the counts for all cell types in the current cluster
  counts <- table_df[, cluster]
  
  # Assign names to counts
  names(counts) <- rownames(table_df)
  
  # Identify the cell type(s) with the maximum count
  max_count <- max(counts)
  dominant_types <- names(counts)[counts == max_count]
  
  # If there's a tie, concatenate the cell types or handle as needed
  if (length(dominant_types) == 1) {
    assigned_cell_types[cluster] <- dominant_types
  } else {
    # Example: Concatenate with a separator
    assigned_cell_types[cluster] <- paste(dominant_types, collapse = "/")
    # Alternatively, leave as "Ambiguous" or handle based on additional criteria
  }
}


corrected@meta.data$Single_r_predict <- assigned_cell_types[as.character(corrected$seurat_clusters)]


corrected.markers <- FindAllMarkers(corrected, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.3)


manual_markers = list(
  "Stem/Progenitor Cells:Hematopoeitic Stem Cells(HSC)" = c("CD34","CD38","Sca1","Kit"),
  "Stem/Progenitor Cells:Lymphoid-primed multipotent progenitors (LMPP)" = c("CD38","CD52","CSF3R","ca1","Kit","CD34","Flk2"),
  "Stem/Progenitor Cells:Common Lymphoid Progenitor (CLP)" = c("IL7R"),
  "Stem/Progenitor Cells:Granulocyte-Monocyte Progenitor (GMP)/Neutrophiles" = c("ELANE"),
  "Stem/Progenitor Cells:Common Myeloid Progenitor" = c("IL3","GM-CSF","M-CSF"),
  "B-cells"=c("CD19"),
  "B-cells:B Cells (B)" = c("CD19","CD20","CD38"),
  "B-cells:Pre B-cell progenitors(Pre B)" = c("CD19","CD34"),
  "B-cells:Plasma" = c("SDC1","IGHA1","IGLC1","MZB1","JCHAIN"),
  "T-cells" = c("CD3D"),
  "T-cells:CD8+ T Cells (CD8)" = c("CD3D","CD3E","CD8A","CD8B"),
  "T-cells:CD4+ T Cells (CD4)" = c("CD3D","CD3E","CD4"),
  "NK cells:Natural Killer Cells (NK)" = c("FCGR3A","NCAM1","NKG7","KLRB1"),
  "Myeloid cells:Erythrocytes"=c("GATA1","HBB","HBA1","HBA2"),
  "Myeloid cells:pDC"=c("IRF8","IRF4","IRF7"),
  "Myeloid cells:cDC"=c("CD1C","CD207","ITGAM","NOTCH2","SIRPA"),
  "Myeloid cells:CD14+ Monocytes (CD14)"=c("CD14","CCL3","CCL4","IL1B"),
  "Myeloid cells:CD16+ Monocytes (CD16)"=c("FCGR3A","CD68","S100A12"),
  "Myeloid cells:Basophils"=c("GATA2")
)

manual_markers <- lapply(manual_markers, toupper)
corrected_copy <-corrected.markers
corrected_copy$gene <- toupper(corrected_copy$gene)

annotations <- names(manual_markers)
clusters <- unique(corrected.markers$cluster)
score_matrix <- matrix(0, nrow = length(annotations), ncol = length(clusters),
                       dimnames = list(annotations, clusters))

for (annotation in annotations) {
  genes <- manual_markers[[annotation]]
  # Filter markers that are in manual markers
  relevant_markers <- corrected.markers %>%
    filter(gene %in% genes)
  
  # Aggregate scores per cluster
  cluster_scores <- relevant_markers %>%
    group_by(cluster) %>%
    summarise(score = sum(avg_log2FC, na.rm = TRUE)) # Option 2
  
  # Populate the score_matrix
  score_matrix[annotation, cluster_scores$cluster] <- cluster_scores$score
}
score_matrix_scaled <- t(scale(t(score_matrix)))

top_annotations <- apply(score_matrix_scaled, 2, function(x) {
  annotation <- names(x)[which.max(x)]
  return(annotation)
})

names(top_annotations) <- colnames(score_matrix_scaled)

corrected$manual_annotation <- plyr::mapvalues(
  x = Idents(corrected),
  from = names(top_annotations),
  to = top_annotations
)

metadata <- corrected@meta.data

cell_counts <- metadata %>%
  group_by(orig.ident, manual_annotation) %>%
  summarise(count = n()) %>%
  ungroup()

cell_props <- cell_counts %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

cell_props$orig.ident <- factor(cell_props$orig.ident, levels = unique(cell_props$orig.ident))

ggplot(cell_props, aes(x = orig.ident, y = proportion, fill = manual_annotation)) +
  geom_bar(stat = "identity") +
  labs(title = "Cell-Type Proportions per Sample",
       x = "orig.ident",
       y = "Proportion (%)",
       fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set3")

unique(metadata$manual_annotation)

# plot umap by manual_annotation and labels of manual_annotation

d_manual_annotation <- DimPlot(corrected, reduction = "umap",group.by =c("manual_annotation", "Single_r_predict"))

# Function to generate Violin and UMAP plots for marker genes
plot_cluster_markers <- function(seurat_obj, meta_col, top_n = 1) {
  # Ensure the column exists in metadata
  if (!meta_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("The column", meta_col, "is not found in the Seurat object metadata."))
  }
  
  # Set the active identity class to the specified metadata column
  Idents(seurat_obj) <- seurat_obj@meta.data[[meta_col]]
  
  # Identify marker genes for each group in the specified metadata column
  marker_genes <- FindAllMarkers(
    object = seurat_obj,
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  
  # Select top marker genes for each group
  top_markers <- marker_genes %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n) %>%
    pull(gene)
  
  # Ensure unique markers and limit to a manageable size
  top_markers <- unique(top_markers)
  if (length(top_markers) > 20) {
    warning("Too many marker genes selected. Limiting to the first 20.")
    top_markers <- top_markers[1:20]
  }
  
  # Violin plot for top markers
  violin_plot <- VlnPlot(
    object = seurat_obj,
    features = top_markers,
    group.by = meta_col,
    pt.size = 0.1,
    combine = TRUE
  ) + ggtitle("Gene Expression Violin Plot by Cluster") + theme(plot.title = element_text(hjust = 0.5))
  
  # UMAP plots for each marker gene
  umap_plots <- lapply(top_markers, function(gene) {
    FeaturePlot(
      object = seurat_obj,
      features = gene,
      reduction = "umap"
    ) + ggtitle(paste("UMAP Plot for", gene)) + theme(plot.title = element_text(hjust = 0.5))
  })
  
  # Combine UMAP plots
  combined_umap <- wrap_plots(umap_plots, ncol = 1)
  
  # Save plots
  ggsave("ViolinPlot_ClusterMarkers.png", violin_plot, width = 12, height = 6)
  ggsave("UMAPPlots_ClusterMarkers.png", combined_umap, width = 12, height = length(top_markers) * 4)
  
  # Return the plots and top markers
  return(list(violin_plot = violin_plot, umap_plots = combined_umap, top_markers = top_markers))
}

plot_cluster_markers(corrected, "manual_annotation")

group1 <- "B-cells"
group2 <- "T-cells"

#corrected <-readRDS("./final_annotated.rds")

# Subset the Seurat object to include only B cells and T cells
b_vs_t_cells <- subset(corrected, subset = manual_annotation %in% c(group1, group2))


unique(corrected@meta.data$manual_annotation)

# Set the identity to 'manual_annotation' for DE analysis
Idents(b_vs_t_cells) <- "manual_annotation"

# Perform DE analysis using the Wilcoxon Rank Sum test
de_b_vs_t <- FindMarkers(
  object = b_vs_t_cells,
  ident.1 = group1,
  ident.2 = group2,
  logfc.threshold = 0.25,       # Adjust threshold as needed
  min.pct = 0.1,                # Adjust minimum percentage as needed
)

# Add gene names as a column
de_b_vs_t$gene <- rownames(de_b_vs_t)

create_volcano_plot <- function(de_results, title, p_cutoff = 0.05, fc_cutoff = 0.25) {
  de_results <- de_results %>%
    mutate(
      neg_log10_pval = -log10(p_val_adj),
      significance = case_when(
        p_val_adj < p_cutoff & avg_log2FC > fc_cutoff ~ "Upregulated",
        p_val_adj < p_cutoff & avg_log2FC < -fc_cutoff ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Define color palette
  colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")
  
  # Create the volcano plot
  p <- ggplot(de_results, aes(x = avg_log2FC, y = neg_log10_pval)) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(
      title = title,
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~Adjusted~P-value),
      color = "Significance"
    ) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black")
  
  return(p)
}
title_b_vs_t <- "Differential Expression: B-cells vs. T-cells"

create_volcano_plot(
  de_results = de_b_vs_t,
  title = title_b_vs_t,
  p_cutoff = 0.05,
  fc_cutoff = 0.25
)


# Define the groups
group1_t_vs_mon <- "T-cells"
group2_t_vs_mon <- "Monocytes"  # Aggregated Monocytes

# Create a new metadata column to aggregate Monocytes
corrected$manual_annotation_agg <- as.character(corrected$manual_annotation)

# Assign "Monocytes" to the specified Monocyte subtypes
monocyte_subtypes <- c("Myeloid cells:CD14+ Monocytes (CD14)", 
                       "Myeloid cells:CD16+ Monocytes (CD16)")
corrected$manual_annotation_agg[corrected$manual_annotation_agg %in% monocyte_subtypes] <- "Monocytes"

corrected$manual_annotation_agg <- as.factor(corrected$manual_annotation_agg)


# Subset the Seurat object to include only T-cells and Monocytes
t_vs_monocytes <- subset(corrected, subset = manual_annotation_agg %in% c(group1_t_vs_mon, group2_t_vs_mon))

# Set the identity to 'manual_annotation_agg' for DE analysis
Idents(t_vs_monocytes) <- "manual_annotation_agg"

# Perform DE analysis using the  Rank Sum test
de_t_vs_monocytes <- FindMarkers(
  object = t_vs_monocytes,
  ident.1 = group1_t_vs_mon,
  ident.2 = group2_t_vs_mon,
  logfc.threshold = 0.25,
  min.pct = 0.1,
)

# Add gene names as a column
de_t_vs_monocytes$gene <- rownames(de_t_vs_monocytes)


title_t_vs_monocytes <- "Differential Expression: T-cells vs. Monocyts"

create_volcano_plot(
  de_results = de_t_vs_monocytes,
  title = title_t_vs_monocytes,
  p_cutoff = 0.05,
  fc_cutoff = 0.25
)

top5_b_vs_t <- de_b_vs_t %>%
  arrange(-abs(avg_log2FC)) %>%        # Sort by adjusted p-value
  slice_head(n = 5) %>%         # Select top 5
  mutate(Comparison = "B vs T") # Ad

top5_t_vs_mon <- de_t_vs_monocytes %>%
  arrange(-abs(avg_log2FC)) %>%                # Sort by adjusted p-value
  slice_head(n = 5) %>%                 # Select top 5
  mutate(Comparison = "T vs Monocytes") #

top_genes_combined <- bind_rows(top5_b_vs_t, top5_t_vs_mon)
top_genes_combined <- top_genes_combined %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj))

volcano_dot_plot <- ggplot(top_genes_combined, aes(x = Comparison, y = gene)) +
  geom_point(aes(color = avg_log2FC, size = -log10(p_val_adj))) +
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0,
    name = expression(Log[2]~Fold~Change)
  ) +
  scale_size_continuous(
    name = "Significance (-Log10 Adjusted P-value)",
    range = c(3, 10)  # Adjust the range as needed
  ) +
  theme_minimal() +
  labs(
    title = "Top 5 Differentially Expressed Genes per Comparison",
    x = "Comparison",
    y = "Gene"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Display the plot
print(volcano_dot_plot)

unique(corrected@meta.data$orig.ident)
# Create a new metadata column 'group' based on 'orig.ident'
corrected$group <- ifelse(grepl("^BMMC", corrected$orig.ident), "BMMC",
                           ifelse(grepl("^CD34", corrected$orig.ident), "CD34", NA))

# Verify the new 'group' column
table(corrected$group)

# Set active identity to 'group'
Idents(corrected) <- "group"

de_results_all <- FindMarkers(corrected,
                              ident.1 = "BMMC",
                              ident.2 = "CD34",
                              test.use = "wilcox",           # Wilcoxon Rank Sum Test
                              logfc.threshold = 0.25,        # Minimum log2 fold change
                              min.pct = 0.1,                 # Minimum percentage of cells expressing the gene
                              only.pos = FALSE)              # Return both up and downregulated genes

# View the top DEGs
head(de_results_all)


top_5_orig_vs_orig_logFC <- de_results_all %>%
  arrange(-abs(avg_log2FC)) %>%        # Sort by adjusted p-value
  slice_head(n = 5)

top_5_orig_vs_orig_adjPval<- de_results_all %>%
  arrange(p_val_adj) %>%        # Sort by adjusted p-value
  slice_head(n = 5)

top_genes_combined_all <- bind_rows(top_5_orig_vs_orig_logFC, top_5_orig_vs_orig_adjPval)


monoc_seurat <- subset(corrected, subset = manual_annotation_agg == "Monocytes")

de_results_all_monocyts  <-FindMarkers(corrected,
                                       ident.1 = "BMMC",
                                       ident.2 = "CD34",
                                       logfc.threshold = 0.25,        # Minimum log2 fold change
                                       min.pct = 0.1,                 # Minimum percentage of cells expressing the gene
                                       only.pos = FALSE) 
de_results_all_monocyts

top_5_monocyts_orig_vs_orig_logFC <- de_results_all_monocyts %>%
  arrange(-abs(avg_log2FC)) %>%        # Sort by adjusted p-value
  slice_head(n = 5)

top_5_monocyts_orig_vs_orig_adjPval<- de_results_all_monocyts %>%
  arrange(p_val_adj) %>%        # Sort by adjusted p-value
  slice_head(n = 5)
top_genes_combined_all_monocyts <- bind_rows(top_5_monocyts_orig_vs_orig_logFC, top_5_monocyts_orig_vs_orig_adjPval)


# get all databases for enrich r
all_databases <- listEnrichrDbs()

dataBases_fetch = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")

# Run Enrichr for all genes in de_results_all
enrichr_results_all <- enrichr(rownames(de_results_all), dataBases_fetch)

# use Deenrichrplot
DEENRICH_GO_BIOLOGICAL_F <- DEenrichRPlot(corrected,ident.1 = "BMMC", ident.2 = "CD34",assay="integrated",enrich.database = dataBases_fetch[1],max.genes =Inf)
DEENRICH_GO_CELLULAR_F <- DEenrichRPlot(corrected,ident.1 = "BMMC", ident.2 = "CD34",assay="integrated",enrich.database = dataBases_fetch[2],max.genes =Inf)
DEENRICH_GO_MOLECULAR_F <- DEenrichRPlot(corrected,ident.1 = "BMMC", ident.2 = "CD34",assay="integrated",enrich.database = dataBases_fetch[3],max.genes =Inf)

DEENRICH_GO_BIOLOGICAL_R <- DEenrichRPlot(corrected,ident.1 = "CD34", ident.2 = "BMMC",assay="integrated",enrich.database = dataBases_fetch[1],max.genes =Inf)
DEENRICH_GO_CELLULAR_R <- DEenrichRPlot(corrected,ident.1 = "CD34", ident.2 = "BMMC",assay="integrated",enrich.database = dataBases_fetch[2],max.genes =Inf)
DEENRICH_GO_MOLECULAR_R <- DEenrichRPlot(corrected,ident.1 = "CD34", ident.2 = "BMMC",assay="integrated",enrich.database = dataBases_fetch[3],max.genes =Inf)


GO.Biological = enrichr_results_all$GO_Biological_Process_2023
GO.Celular = enrichr_results_all$GO_Cellular_Component_2023
GO.Molecular = enrichr_results_all$GO_Molecular_Function_2023

# sort by  p value in GO biological
GO.Biological <- GO.Biological[order(GO.Biological$P.value),]
GO.Celular <- GO.Celular[order(GO.Celular$P.value),]
GO.Molecular <- GO.Molecular[order(GO.Molecular$P.value),]

firstBiological <- GO.Biological[1,]
firstCelular <- GO.Celular[1,]
firstMolecular <- GO.Molecular[1,]

# save corrected
saveRDS(corrected, file = "final_annotated.rds")

# save all DEENRICH_GO_BIOLOGICAL_F
saveRDS(DEENRICH_GO_BIOLOGICAL_F, file = "DEENRICH_GO_BIOLOGICAL_F.rds")
saveRDS(DEENRICH_GO_CELLULAR_F, file = "DEENRICH_GO_CELLULAR_F.rds")
saveRDS(DEENRICH_GO_MOLECULAR_F, file = "DEENRICH_GO_MOLECULAR_F.rds")

saveRDS(DEENRICH_GO_BIOLOGICAL_R, file = "DEENRICH_GO_BIOLOGICAL_R.rds")
saveRDS(DEENRICH_GO_CELLULAR_R, file = "DEENRICH_GO_CELLULAR_R.rds")
saveRDS(DEENRICH_GO_MOLECULAR_R, file = "DEENRICH_GO_MOLECULAR_R.rds")

saveRDS(GO.Biological, file = "GO.Biological.rds")
saveRDS(GO.Celular, file = "GO.Celular.rds")
saveRDS(GO.Molecular, file = "GO.Molecular.rds")



DEENRICH_GO_BIOLOGICAL_F <- readRDS("./DEENRICH_GO_BIOLOGICAL_F.rds")
DEENRICH_GO_CELLULAR_F <- readRDS("./DEENRICH_GO_CELLULAR_F.rds")
DEENRICH_GO_MOLECULAR_F <- readRDS("./DEENRICH_GO_MOLECULAR_F.rds")

DEENRICH_GO_BIOLOGICAL_R <- readRDS("./DEENRICH_GO_BIOLOGICAL_R.rds")
DEENRICH_GO_CELLULAR_R <- readRDS("./DEENRICH_GO_CELLULAR_R.rds")
DEENRICH_GO_MOLECULAR_R <- readRDS("./DEENRICH_GO_MOLECULAR_R.rds")

GO.Biological <- readRDS("./GO.Biological.rds")
GO.Celular <- readRDS("./GO.Celular.rds")
GO.Molecular <- readRDS("./GO.Molecular.rds")



rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])


# select a group of cells  for trajectory

corrected <- readRDS(file = "final_annotated.rds")

unique(corrected@meta.data$manual_annotation) 

lymphoid_annotations <- c(
  "Stem/Progenitor Cells:Common Lymphoid Progenitor (CLP)",                  
  "B-cells",                                                                 
  "T-cells",                                                                 
  "Stem/Progenitor Cells:Granulocyte-Monocyte Progenitor (GMP)/Neutrophiles",
  "NK cells:Natural Killer Cells (NK)",                                      
  "Stem/Progenitor Cells:Hematopoeitic Stem Cells(HSC)"    
)



# Subset the Seurat object to include only Lymphoid lineage cells
lymphoid_seurat <- subset(corrected, subset = manual_annotation %in% lymphoid_annotations)


lymphoid_seurat[["RNA"]] <-JoinLayers(lymphoid_seurat[["RNA"]])

counts_rna <- GetAssayData(lymphoid_seurat, assay = "RNA", slot = "counts")


gene_metadata <- data.frame(gene_short_name = rownames(counts_rna))

rownames(gene_metadata) <- rownames(counts_rna)


# Convert the Seurat object to Monocle 3 CellDataSet
lymphoid_cds <- new_cell_data_set(
  expression_data = counts_rna,
  cell_metadata = lymphoid_seurat@meta.data,
  gene_metadata = gene_metadata
)


# add my manual anotation meta data to lymphoid_cds
lymphoid_cds$manual_annotation <- lymphoid_seurat$manual_annotation


lymphoid_cds<- estimate_size_factors(lymphoid_cds)

lymphoid_cds <- preprocess_cds(lymphoid_cds,num_dim = 100, method = "PCA")

lymphoid_cds <- align_cds(lymphoid_cds,"PCA",alignment_group = "manual_annotation")

lymphoid_cds <- reduce_dimension(lymphoid_cds,max_components = 2,reduction_method = "UMAP",preprocess_method ="PCA" )

lymphoid_cds<-cluster_cells(lymphoid_cds,reduction_method = "UMAP")

# plot umap the clusters

lymphoid_cds <- learn_graph(lymphoid_cds)

plot_cells(
  lymphoid_cds,
  color_cells_by = "manual_annotation", 
) + ggtitle("Trajectory Analysis of Lymphoid Lineage Cells")

lymphoid_cds <- order_cells(lymphoid_cds,root_cells = c("Stem/Progenitor Cells:Hematopoeitic Stem Cells(HSC)"))

plot_cells(lymphoid_cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE)




