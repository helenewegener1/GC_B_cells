library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
# seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")
seurat_obj_singlets_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")

# TODO: Remember to clean up clusters in this flow OR at an earlier point.
# Like removing contaminating clusters. 

# Initialize filtered list
seurat_obj_singlets_annotated_list <- rep(0, length(seurat_obj_singlets_list)) %>% as.list()
names(seurat_obj_singlets_annotated_list) <- names(seurat_obj_singlets_list)


sample_names <- names(seurat_obj_singlets_list)

################################################################################
sample_name <- "HH119-SI-PP-CD19-Pool1" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Memory_B_cells", 
  "1" = "Naive_B_cells",
  "2" = "Plasma_cells",
  "3" = "GC_B_cells",
  "4" = "Low_quality",
  "5" = "GC_B_cells",
  "6" = "Plasma_cells",
  "7" = "T_cells"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# See DimPlot
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE) + NoLegend() 

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-CD19-Pool2" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naive_B_cells",
  "1" = "Memory_B_cells",
  "2" = "Plasma_cells",
  "3" = "Low_quality",
  "4" = "GC_B_cells",
  "5" = "GC_B_cells",
  "6" = "Plasma_cells"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# See DimPlot
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE) + NoLegend() 

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "GC_B_cells",
  "1" = "TFH",
  "2" = "GC_B_cells",
  "3" = "Plasma_blasts",
  "4" = "GC_B_cells",
  "5" = "T_cells",
  "6" = "Contamination_stromal?"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# See DimPlot
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE) + NoLegend() 

################################################################################
sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "GC_B_cells",
  "1" = "GC_B_cells",
  "2" = "TFH",
  "3" = "GC_B_cells",
  "4" = "GC_B_cells",
  "5" = "Plasma_blasts",
  "6" = "Plasma_cells",
  "7" = "TFH"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# See DimPlot
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE) + NoLegend() 

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Memory_B_cells", 
  "1" = "T_cells",
  "2" = "GC_B_cells",
  "3" = "GC_B_cells",
  "4" = "Plasma_cells",
  "5" = "B_cells_ribosomal",
  "6" = "DCs",
  "7" = "B_cells_ambient",
  "8" = "GC_B_cells",
  "9" = "DCs"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# See DimPlot
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE) + NoLegend() 

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################

# sample_names
# 
# sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
# 
# seurat_obj <- seurat_obj_singlets_list[[sample_name]]
# 
# # Map clusters to cell types
# cluster_to_celltype <- c(
#   "0" = "Memory_B_cells", 
#   "1" = "T_cells",
#   "2" = "GC_B_cells",
#   "3" = "GC_B_cells",
#   "4" = "Plasma_cells",
#   "5" = "B_cells_ribosomal",
#   "6" = "DCs",
#   "7" = "B_cells_ambient",
#   "8" = "GC_B_cells",
#   "9" = "DCs"
# )
# 
# celltype_broad <- cluster_to_celltype[
#   as.character(seurat_obj$seurat_clusters)
# ] %>% as.data.frame()
# 
# rownames(celltype_broad) <- colnames(seurat_obj)
# 
# # Add broad cell type annotations
# seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")
# 
# # Check that mapping when correctly
# table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)
# 
# # See DimPlot
# DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE) + NoLegend() 

################################################################################

# Export seurat object list 
saveRDS(seurat_obj_singlets_annotated_list, "11_broad_annotation/out/seurat_obj_singlets_annotated_list.rds")


