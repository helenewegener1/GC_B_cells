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

for (HH in names(all_clones)){
  
  # HH <- "HH117"
  
  fols <- seurat_integrated[[]] %>% filter(patient == HH) %>% pull(manual_ADT_ID) %>% unique()
  fols <- fols %>% na.omit()
  fols <- fols[!(fols %in% c("Negative", "Doublet"))]
  
  for (fol in fols){
    
    # fol <- "Fol-1"
    # Subset seurat object 
    # sub <- subset(seurat_integrated, patient == HH & manual_ADT_ID == fol & L1_annotation != "Tfh_cells")
    sub <- subset(seurat_integrated, patient == HH & manual_ADT_ID == fol & L1_annotation %in% c("GC_B_cells", "Naive_Bcells"))
    ncol(sub)
    
    # Re-normalize 
    sub <- NormalizeData(sub)
    sub <- FindVariableFeatures(sub)
    sub <- ScaleData(sub)
    sub <- RunPCA(sub, npcs = 30)
    
    # Re-integrate 
    # No need to integrate since we are working with 1 patient 
    
    Reductions(sub)
    
    # Re-cluster on new embedding
    sub <- FindNeighbors(sub, reduction = "pca")
    sub <- FindClusters(sub, res = 0.3)
    sub <- RunUMAP(sub, reduction = "pca", dims = 1:30)
    
    # Plot
    DimPlot(sub, group.by = "clone_subgroup_id", reduction = "umap") + NoLegend()
    
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
    
    # prep outdir 
    outdir <- glue("49_trajectory_analysis/plot/{HH}_{fol}/")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    # ------------------------------------------------------------------------------
    # slingshot 
    # ------------------------------------------------------------------------------
    
    # Run slingshot
    # sce <- slingshot(sce, clusterLabels = clusterLabel, reducedDim = 'PCA')
    # sce <- slingshot(sce, clusterLabels = 'L1_annotation', reducedDim = 'PCA', start.clus = "Naive_Bcells")
    sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA', start.clus = "3")
    
    # Plot
    summary(sce$slingPseudotime_1)
    
    # Trajectory
    png(glue("{outdir}/{HH}_{fol}_slingshot.png"), res = 1000, width = 9000, height = 6000)
    
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
    plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: {fol}"))
    lines(SlingshotDataSet(sce), lwd=2, col='black')
    
    dev.off()
    
    # # Seurat clusters 
    # colors <- c("red", "blue")
    # plotcol <- colors[sce$seurat_clusters]
    # plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: clone {clone}"))
    # lines(SlingshotDataSet(sce), lwd=2, col='black')
    
    # Isotype
    png(glue("{outdir}/{HH}_{fol}_isotype.png"), res = 1000, width = 9000, height = 6000)
    
    colors <- isotype_colors_custom
    plotcol <- colors[sce$c_call]
    plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: {fol}"))
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
    png(glue("{outdir}/{HH}_{fol}_L1.png"), res = 1000, width = 9000, height = 6000)
    
    colors <- L1_colors
    plotcol <- colors[sce$L1_annotation]
    plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1, main = glue("{HH}: {fol}"))
    lines(SlingshotDataSet(sce), lwd=2, col='black')
    
    present <- unique(sce$L1_annotation)
    legend("topright",
           legend = names(colors[present]),
           col = colors[present],
           pch = 16,
           bty = "n",
           cex = 0.8)
    
    dev.off()
    
    # Clone
    png(glue("{outdir}/{HH}_{fol}_clone.png"), res = 1000, width = 9000, height = 6000)
    
    present <- unique(sce$clone_subgroup_id)
    
    # Count cells per clone and get top 20
    clone_counts <- sort(table(sce$clone_subgroup_id), decreasing = TRUE)
    top20 <- names(clone_counts)[1:min(20, length(clone_counts))]
    
    # Make colors: distinct colors for top 20, grey for the rest
    top20_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(top20)), top20)
    all_colors <- setNames(rep("grey80", length(present)), present)
    all_colors[top20] <- top20_colors[top20]
    
    colors <- all_colors
    plotcol <- colors[sce$clone_subgroup_id]
    plot(reducedDims(sce)$PCA, col = plotcol, pch = 16, asp = 1, main = glue("{HH}: {fol}"))
    lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
    
    # Only show top 20 in legend (grey would be too many entries)
    legend("topright",
           legend = c(names(top20_colors), "Other"),
           col = c(top20_colors, "grey80"),
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






