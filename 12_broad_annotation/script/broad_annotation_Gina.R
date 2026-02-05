library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
# seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")
seurat_obj_singlets_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")
seurat_obj_singlets_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")


# TODO: Remember to clean up clusters in this flow OR at an earlier point.
# Like removing contaminating clusters. 

# Initialize filtered list
seurat_obj_singlets_annotated_list <- rep(0, length(seurat_obj_singlets_list)) %>% as.list()
names(seurat_obj_singlets_annotated_list) <- names(seurat_obj_singlets_list)

sample_names <- names(seurat_obj_singlets_list)

# Defining color scheme for each cell type for streamlined plotting
source("12_broad_annotation/script/color_palette.R")

################################################################################
sample_name <- "HH119-SI-PP-CD19-Pool1" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "Naïve_memory_B_cells",
  "2" = "GC_like_B_cells",
  "3" = "GC_like_B_cells",
  "4" = "Naïve_memory_B_cells",
  "5" = "GC_like_B_cells",
  "6" = "PCs_PBs",
  "7" = "Contamination_ambiguous"#(T_cell-MNP signature) 
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-CD19-Pool2" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells",
  "1" = "Naïve_memory_B_cells",
  "2" = "GC_like_B_cells",
  "3" = "Naïve_memory_B_cells",
  "4" = "GC_like_B_cells",
  "5" = "PCs_PBs"
  #"6" = "Plasma_cells"/ I don't have a 6th cluster 
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "GC_like_B_cells",
  "1" = "Tfh_like_cells",
  "2" = "GC_like_B_cells",
  "3" = "PCs_PBs",
  "4" = "GC_like_B_cells",
  "5" = "Tfh_like_cells",
  "6" = "GC_like_B_cells" # yes, there might be some contamination here but we keep for now
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "GC_like_B_cells",
  "1" = "GC_like_B_cells",
  "2" = "Tfh_like_cells",
  "3" = "GC_like_B_cells",
  "4" = "GC_like_B_cells",
  "5" = "Naïve_memory_B_cells",
  "6" = "PCs_PBs",
  "7" = "Tfh_like_cells"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "Naïve_memory_B_cells",
  "2" = "Naïve_memory_B_cells",
  "3" = "Naïve_memory_B_cells",
  "4" = "PCs_PBs",
  "5" = "Tfh_like_cells",
  "6" = "Naïve_memory_B_cells"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-COLP-PC" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs",
  "2" = "Contamination_stroma"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "Naïve_memory_B_cells",
  "2" = "Tfh_like_cells",
  "3" = "GC_like_B_cells",
  "4" = "GC_like_B_cells",
  "5" = "PCs_PBs",
  "6" = "Naïve_memory_B_cells", 
  "7" = "Tfh_like_cells"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SILP-PC" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs",
  "2" = "PCs_PBs",
  "3" = "GC_like_B_cells",
  "4" = "Contamination_mast_cells",
  "5" = "Contamination_MNPs",
  "6" = "Contamination_γδT_cell"
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
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "Tfh_like_cells",
  "2" = "GC_B_cells",
  "3" = "GC_B_cells",
  "4" = "PCs_PBs",
  "5" = "Naïve_memory_B_cells",
  "6" = "DCs_MNPs",
  "7" = "Naïve_memory_B_cells",
  "8" = "GC_B_cells",
  "9" = "DCs_MNPs"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SI-MILF-INF-HLADR-AND-CD19"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "DCs_MNPs",
  "2" = "Naïve_memory_B_cells",
  "3" = "PCs_PBs",
  "4" = "Contamination_stroma",
  "5" = "PCs_PBs"
  
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SI-MILF-nonINF-HLADR-AND-CD19"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "DCs_MNPs",
  "2" = "PCs_PBs",
  "3" = "DCs_MNPs",
  "4" = "DCs_MNPs",
  "5" = "DCs_MNPs",
  "6" = "GC_B_cells",
  "7" = "Contamination_stroma"
  
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SILP-INF-PC"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs",
  "2" = "PCs_PBs",
  "3" = "PCs_PBs",
  "4" = "PCs_PBs",
  "5" = "Contamination_mast_cells"
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SILP-nonINF-PC"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs",
  "2" = "PCs_PBs",
  "3" = "PCs_PBs" #cluster 3 possibly contamination or early PBs
 
)

celltype_broad <- cluster_to_celltype[
  as.character(seurat_obj$seurat_clusters)
] %>% as.data.frame()

rownames(celltype_broad) <- colnames(seurat_obj)

# Add broad cell type annotations
seurat_obj <- AddMetaData(seurat_obj, celltype_broad, "celltype_broad")

# Check that mapping when correctly
table(seurat_obj$seurat_clusters, seurat_obj$celltype_broad)

# # Subset: remove contamination_* and DCs_MNPs
# seurat_obj <- subset(
#  seurat_obj,
#  subset = !grepl("^Contamination_|^DCs_MNPs", celltype_broad))
# # Check removal 
# table(seurat_obj$celltype_broad)

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("12_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

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
saveRDS(seurat_obj_singlets_annotated_list, "12_broad_annotation/out/seurat_obj_singlets_annotated_list.rds")


