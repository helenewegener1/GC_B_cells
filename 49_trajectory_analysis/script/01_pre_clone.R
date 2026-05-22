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

# HH <- "HH117"
# clone <- "1849_1" # 3rd clone

# HH <- "HH119"
# clone <- "28075_1" # 3rd clone

all_clones <- list(
  "HH117" = c(
    "4221_1", "2628_1", "1849_1", "3709_1"#, "2301_1",
    # "1320_1", "5941_1", "6115_1", "1910_1", "2169_1"
  ),
  
  "HH119" = c(
    "28075_1", "12120_1", "15287_1", "23124_1",
    "8372_1", "3869_1", "25158_1", "7913_1"
  )
)

for (HH in names(all_clones)){
  
  HH_clones <- all_clones[[HH]]
  
  for (clone in HH_clones){
    
    # Subset seurat object 
    sub <- subset(seurat_integrated, patient == HH & clone_subgroup_id == clone)
    ncol(sub)
    
    # Re-normalize 
    sub <- NormalizeData(sub)
    sub <- FindVariableFeatures(sub)
    sub <- ScaleData(sub)
    sub <- RunPCA(sub, npcs = 10)
    
    # Re-integrate 
    # No need to integrate since we are working with 1 patient 
    
    Reductions(sub)
    
    # Re-cluster on new embedding
    sub <- FindNeighbors(sub, reduction = "pca")
    sub <- FindClusters(sub, res = 0.3)
    sub <- RunUMAP(sub, reduction = "pca", dims = 1:10)
    
    # Plot
    DimPlot(sub, group.by = "sample_clean_plot", reduction = "umap")
    
    DimPlot(sub, group.by = "L1_annotation", reduction = "umap") + 
      scale_color_manual(values = L1_colors)
    
    DimPlot(sub, group.by = "L1_annotation", reduction = "pca") + 
      scale_color_manual(values = L1_colors)
    
    DimPlot(sub, group.by = "c_call", reduction = "umap") + 
      scale_color_manual(values = isotype_colors_custom)
    
    DimPlot(sub, group.by = "seurat_clusters", reduction = "umap")
    
    # ------------------------------------------------------------------------------
    # slingshot - prep
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
    
    # Define clusterLabel 
    counts <- table(sce$junction)
    clusterLabel <- ifelse(counts[sce$junction] == 1, -1, sce$junction)
    table(clusterLabel)
    
    # prep outdir 
    outdir <- glue("49_trajectory_analysis/plot/{HH}_{clone}/")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    # ------------------------------------------------------------------------------
    # slingshot 
    # ------------------------------------------------------------------------------
    
    # Run slingshot
    # sce <- slingshot(sce, clusterLabels = clusterLabel, reducedDim = 'PCA')
    sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')
    
    # Plot
    summary(sce$slingPseudotime_1)
    
    # Trajectory
    png(glue("{outdir}/{HH}_{clone}_slingshot.png"), res = 1000, width = 9000, height = 6000)
    
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
    plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
    lines(SlingshotDataSet(sce), lwd=2, col='black')
    
    dev.off()
    
    # # Seurat clusters 
    # colors <- c("red", "blue")
    # plotcol <- colors[sce$seurat_clusters]
    # plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
    # lines(SlingshotDataSet(sce), lwd=2, col='black')
    
    # Isotype
    png(glue("{outdir}/{HH}_{clone}_isotype.png"), res = 1000, width = 9000, height = 6000)
    
    colors <- isotype_colors_custom
    plotcol <- colors[sce$c_call]
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
    
    png(glue("{outdir}/{HH}_{clone}_L1.png"), res = 1000, width = 9000, height = 6000)
    
    colors <- L1_colors
    plotcol <- colors[sce$L1_annotation]
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
    
    # sample
    png(glue("{outdir}/{HH}_{clone}_sample.png"), res = 1000, width = 9000, height = 6000)
    
    colors <- sample_clean_plot_colors
    plotcol <- colors[sce$sample_clean_plot]
    plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
    lines(SlingshotDataSet(sce), lwd=2, col='black')
    
    present <- unique(sce$sample_clean_plot)
    legend("topright",
           legend = names(colors[present]),
           col = colors[present],
           pch = 16,
           bty = "n",
           cex = 0.8)
    
    dev.off()
    
  }
  
}




# Combine plots 


imgs <- image_read(list.files("49_trajectory_analysis/plot/HH117_1849_1/", full.names = TRUE))

# e.g. 2 columns
rows <- split(imgs, ceiling(seq_along(imgs) / 2))
rows <- lapply(rows, function(r) image_append(image_join(r)))
combined <- image_append(image_join(rows), stack = TRUE)

image_write(combined, "combined_grid.png")


# ------------------------------------------------------------------------------
# slingshot::getLineages
# ------------------------------------------------------------------------------

# lin1 <- getLineages(sce, clusterLabel)
lin1 <- getLineages(sce, 'seurat_clusters')

# Isotype
colors <- isotype_colors_custom
plotcol <- colors[sce$c_call]
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
lines(SlingshotDataSet(lin1), lwd=2, col='black')

present <- unique(sce$c_call)
legend("topright",
       legend = names(colors[present]),
       col = colors[present],
       pch = 16,
       bty = "n",
       cex = 0.8)

# Cell type
colors <- L1_colors
plotcol <- colors[sce$L1_annotation]
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
lines(SlingshotDataSet(lin1), lwd=2, col='black')

present <- unique(sce$L1_annotation)
legend("topright",
       legend = names(colors[present]),
       col = colors[present],
       pch = 16,
       bty = "n",
       cex = 0.8)






