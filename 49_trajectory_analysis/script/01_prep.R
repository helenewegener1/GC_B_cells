library(glue)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(grDevices)
library(RColorBrewer)

source("10_broad_annotation/script/color_palette.R")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")
meta <- readRDS("45_immcantation/out/rds/meta_4_Gina_list.rds")

# Add meta data
seurat_integrated <- AddMetaData(seurat_integrated, meta)

# Wrangle
seurat_integrated[[]] <- seurat_integrated[[]] %>% 
  mutate(
    L1_annotation = ifelse(L1_annotation == "GC_Bcells", "GC_B_cells", L1_annotation),
    sample_clean_plot = paste0(sample_clean, "_", manual_ADT_ID),
    sample_clean_plot = str_remove_all(sample_clean_plot, "_NA")
  )

# Check distribution 
seurat_integrated$L1_annotation %>% table()
seurat_integrated$patient %>% table()
seurat_integrated$sample_clean_plot %>% table()

# ------------------------------------------------------------------------------
# Warm-up plots
# ------------------------------------------------------------------------------

DimPlot(seurat_integrated, group.by = "L1_annotation", split.by = "patient", reduction = "RNA_umap.harmony") + 
  scale_color_manual(values = L1_colors)

# ------------------------------------------------------------------------------
# Investigate largest CRC clone
# ------------------------------------------------------------------------------

HH <- "HH117"
clone <- "1849_1"

# clone_nr <- 3
# # Get clone
# clone <- seurat_integrated[[]] %>% 
#   filter(
#     patient == HH, 
#     !is.na(clone_subgroup_id)
#   ) %>% 
#   dplyr::count(clone_subgroup_id, sort = TRUE) %>% 
#   dplyr::slice(clone_nr) %>% 
#   pull(clone_subgroup_id)

# Subset seurat object 
sub <- subset(seurat_integrated, patient == HH & clone_subgroup_id == clone)
ncol(sub)

# Re-normalize 
sub <- NormalizeData(sub)
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)

# Re-integrate 
# No need to integrate since we are working with 1 patient 

Reductions(sub)

# Re-cluster on new embedding
sub <- FindNeighbors(sub, reduction = "pca")
sub <- FindClusters(sub, res = 0.2)
sub <- RunUMAP(sub, reduction = "pca", dims = 1:30)

# Plot
DimPlot(sub, group.by = "sample_clean_plot", reduction = "umap")

DimPlot(sub, group.by = "L1_annotation", reduction = "umap") + 
  scale_color_manual(values = L1_colors)

DimPlot(sub, group.by = "c_call", reduction = "umap") + 
  scale_color_manual(values = isotype_colors_custom)

DimPlot(sub, group.by = "seurat_clusters", reduction = "umap")

# ------------------------------------------------------------------------------
# slingshot 
# ------------------------------------------------------------------------------

# Transform Seurat object to single cell experiment 
sce <- SingleCellExperiment(
  assays = list(
    counts = GetAssayData(sub, layer = "counts"),
    logcounts = GetAssayData(sub, layer = "data")
  ),
  colData = sub@meta.data,
  reducedDims = list(
    PCA = Embeddings(sub, "pca"),
    UMAP = Embeddings(sub, "umap")
  )
)

# N cells 
ncol(sce)

# N unique clones 
sce$junction %>% unique() %>% length()

# Run slingshot
# sce <- slingshot(sce, clusterLabels = 'junction', reducedDim = 'PCA')
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')

# Plot
summary(sce$slingPseudotime_1)

# Trajectory
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

png(glue("49_trajectory_analysis/plot/{HH}_{clone}_slingshot.png"), res = 1000, width = 4000, height = 3000)

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
lines(SlingshotDataSet(sce), lwd=2, col='black')

dev.off()

# Isotype
colors <- isotype_colors_custom
plotcol <- colors[sce$c_call]

png(glue("49_trajectory_analysis/plot/{HH}_{clone}_isotype.png"), res = 1000, width = 4000, height = 3000)
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
lines(SlingshotDataSet(sce), lwd=2, col='black')

present <- unique(sce$c_call)
legend("topright",
       legend = names(colors[present]),
       col = colors[present],
       pch = 16,
       bty = "n",
       cex = 0.8)

dev.off()

# Cell type
colors <- L1_colors
plotcol <- colors[sce$L1_annotation]

png(glue("49_trajectory_analysis/plot/{HH}_{clone}_L1.png"), res = 1000, width = 4000, height = 3000)

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
lines(SlingshotDataSet(sce), lwd=2, col='black')

present <- unique(sce$L1_annotation)
legend("topright",
       legend = names(colors[present]),
       col = colors[present],
       pch = 16,
       bty = "n",
       cex = 0.8)

dev.off()












