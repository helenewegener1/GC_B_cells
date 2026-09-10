library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
# seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")
seurat_obj_singlets_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")

# Initialize filtered list
seurat_obj_singlets_annotated_list <- rep(0, length(seurat_obj_singlets_list)) %>% as.list()
names(seurat_obj_singlets_annotated_list) <- names(seurat_obj_singlets_list)

sample_names <- names(seurat_obj_singlets_list)

# Defining color scheme for each cell type for streamlined plotting
source("10_broad_annotation/script/color_palette.R")

# ANNOTATIONS ARE BASED ON RESLOLUTION 0.3

################################################################################
sample_name <- "HH119-SI-PP-CD19-Pool1" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "Naïve_memory_B_cells",
  "2" = "GC_B_cells",
  "3" = "GC_B_cells",
  "4" = "GC_B_cells",
  "5" = "PCs_PBs",
  "6" = "Naïve_memory_B_cells",
  "7" = "Naïve_memory_B_cells", 
  "8" = "Contamination_ambiguous", # T cell & MNP signature
  "9" = "Naïve_memory_B_cells"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-CD19-Pool2" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells",
  "1" = "Naïve_memory_B_cells",
  "2" = "GC_B_cells",
  "3" = "GC_B_cells",
  "4" = "GC_B_cells",
  "5" = "Naïve_memory_B_cells", 
  "6" = "Naïve_memory_B_cells", 
  "7" = "GC_B_cells", 
  "8" = "PCs_PBs"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "GC_B_cells",
  "1" = "Tfh_like_cells",
  "2" = "GC_B_cells",
  "3" = "GC_B_cells",
  "4" = "PCs_PBs",
  "5" = "GC_B_cells",
  "6" = "Tfh_like_cells", 
  "7" = "Naïve_memory_B_cells", # contaminating and not completely sure
  "8" = "GC_B_cells", 
  "9" = "PCs_PBs", 
  "10" = "PCs_PBs"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "GC_B_cells",
  "1" = "GC_B_cells",
  "2" = "Tfh_like_cells",
  "3" = "GC_B_cells",
  "4" = "GC_B_cells",
  "5" = "Naïve_memory_B_cells",
  "6" = "GC_B_cells", # a bit unsure about this cluster but likely GCB
  "7" = "GC_B_cells", 
  "8" = "PCs_PBs", 
  "9" = "Tfh_like_cells", 
  "10" = "Tfh_like_cells"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

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
  "4" = "Naïve_memory_B_cells",
  "5" = "Naïve_memory_B_cells",
  "6" = "Naïve_memory_B_cells", 
  "7" = "PCs_PBs", 
  "8" = "Naïve_memory_B_cells", 
  "9" = "Tfh_like_cells", 
  "10" = "PCs_PBs"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH119-COLP-PC" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs",
  "2" = "PCs_PBs", 
  "3" = "PCs_PBs",
  "4" = "PCs_PBs",
  "5" = "PCs_PBs",
  "6" = "Contamination_stroma", 
  "7" = "Contamination_mast_cells", 
  "8" = "PCs_PBs" 
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

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
  "3" = "GC_B_cells",
  "4" = "PCs_PBs",
  "5" = "GC_B_cells",
  "6" = "GC_B_cells", 
  "7" = "Naïve_memory_B_cells", 
  "8" = "Tfh_like_cells",
  "9" = "PCs_PBs",
  "10" = "Naïve_memory_B_cells"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

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
  "3" = "PCs_PBs",
  "4" = "PCs_PBs",
  "5" = "Naïve_memory_B_cells_&_GC_B_cells",
  "6" = "Contamination_mast_cells", 
  "7" = "Contamination_MNPs", # DCs or macrophage contamination
  "8" = "Contamination_T_cells" 
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "GC_B_cells",
  "2" = "Tfh_like_cells",
  "3" = "GC_B_cells",
  "4" = "Naïve_memory_B_cells",
  "5" = "PCs_PBs",
  "6" = "DCs_MNPs",
  "7" = "Naïve_memory_B_cells",
  "8" = "Tfh_like_cells",
  "9" = "GC_B_cells",
  "10" = "PCs_PBs", 
  "11" = "DCs_MNPs"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SI-MILF-INF-HLADR-AND-CD19"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs",
  "2" = "DCs_MNPs",
  "3" = "Naïve_memory_B_cells",
  "4" = "DCs_MNPs",
  "5" = "PCs_PBs", 
  "6" = "Contamination_ambiguous", # stroma & epithelial signature
  "7" = "PCs_PBs", 
  "8" = "PCs_PBs", 
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH117-SI-MILF-nonINF-HLADR-AND-CD19"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "Naïve_memory_B_cells",
  "2" = "DCs_MNPs",
  "3" = "PCs_PBs",
  "4" = "DCs_MNPs",
  "5" = "DCs_MNPs",
  "6" = "DCs_MNPs",
  "7" = "DCs_MNPs", 
  "8" = "GC_B_cells", 
  "9" = "Contamination_stroma"
  
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

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
  "5" = "PCs_PBs", 
  "6" = "Contamination_mast_cells", 
  "7" = "PCs_PBs"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

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
  "3" = "PCs_PBs" , 
  "4" = "PCs_PBs", 
  "5" = "PCs_PBs", 
  "6" = "PCs_PBs", # maybe contaminating B cells

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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

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


################################################################################
sample_name <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Blue"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "2" = "Naïve_memory_B_cells", 
  "3" = "Tfh_like_cells_&_GC_B_cells",   # To few cells for the cluster to split, hence mixing       
  "4" = "Naïve_memory_B_cells",        
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

# See DimPlot
n_cells <- ncol(seurat_obj)
DimPlot(seurat_obj, group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
  labs(subtitle = sample_name, 
       caption = glue("N cells: {n_cells}"))
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Green"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Tfh_like_cells", 
  "1" = "Tfh_like_cells",       
  "2" = "Naïve_memory_B_cells", 
  "3" = "Naïve_memory_B_cells",           
  "4" = "GC_B_cells", 
  "5" = "Naïve_memory_B_cells",               
  "6" = "GC_B_cells", 
  "7" = "Naïve_memory_B_cells", 
  "8" = "PCs_PBs", 
  "9" = "Contamination_ambiguous" # Probably both Myeloid/DCs and mast cells
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Red"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", # somewhere there might be some GCB cells but so so few
  "1" = "Naïve_memory_B_cells", 
  "2" = "Naïve_memory_B_cells", 
  "3" = "PCs_PBs",              
  "4" = "Tfh_like_cells"        
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Yellow"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "Naïve_memory_B_cells", 
  "2" = "Tfh_like_cells",       
  "3" = "Naïve_memory_B_cells",           
  "4" = "Naïve_memory_B_cells", 
  "5" = "Tfh_like_cells",               
  "6" = "GC_B_cells", 
  "7" = "PCs_PBs", 
  "8" = "Contamination_MNPs" # DCs or macrophage contamination
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH151-SILP-INF-PC"

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs",              
  "1" = "PCs_PBs",              
  "2" = "PCs_PBs",             
  "3" = "PCs_PBs", 
  "4" = "PCs_PBs",    
  "5" = "Contamination_MNPs", # DCs or macrophage contamination 
  "6" = "Contamination_mast_cells", 
  "7" = "PCs_PBs", 
  "8" = "Contamination_T_cells"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH151-SILP-nonINF-PC" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs", 
  "2" = "PCs_PBs", 
  "3" = "PCs_PBs", 
  "4" = "Contamination_MNPs", # DCs or macrophage contamination 
  "5" = "Contamination_mast_cells", 
  "6" = "PCs_PBs", 
  "7" = "Contamination_T_cells",
  "8" = "PCs_PBs", 
  "9" = "Contamination_stroma"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH153-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB-Pool1" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "Naïve_memory_B_cells", 
  "1" = "GC_B_cells", 
  "2" = "Tfh_like_cells", 
  "3" = "GC_B_cells", 
  "4" = "GC_B_cells", 
  "5" = "GC_B_cells", 
  "6" = "GC_B_cells", 
  "7" = "GC_B_cells",
  "8" = "Tfh_like_cells", 
  "9" = "GC_B_cells", 
  "10" = "Naïve_memory_B_cells", 
  "11" = "PCs_PBs", 
  "12" = "Contamination_MNPs" # DCs or macrophage contamination 
  
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH153-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB-Pool2" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "GC_B_cells", 
  "1" = "GC_B_cells", 
  "2" = "Naïve_memory_B_cells", 
  "3" = "Tfh_like_cells", 
  "4" = "GC_B_cells", 
  "5" = "GC_B_cells", 
  "6" = "Tfh_like_cells", 
  "7" = "Contamination_MNPs", # DCs or macrophage contamination 
  "8" = "PCs_PBs"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH153-SILP-INF-PC" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs", 
  "2" = "PCs_PBs", 
  "3" = "Contamination_stroma", 
  "4" = "Contamination_epithelial_cells", 
  "5" = "Contamination_ambiguous",  # T cell and MNPs signature
  "6" = "Contamination_ambiguous", # T cell and MNPs signature
  "7" = "Contamination_mast_cells"
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

################################################################################
sample_name <- "HH153-SILP-nonINF-PC" 

seurat_obj <- seurat_obj_singlets_list[[sample_name]]

# Map clusters to cell types
cluster_to_celltype <- c(
  "0" = "PCs_PBs", 
  "1" = "PCs_PBs", 
  "2" = "PCs_PBs", 
  "3" = "Contamination_mast_cells", 
  "4" = "Contamination_stroma", 
  "5" = "Contamination_ambiguous", # B cell and myeloid signature
  "6" = "Contamination_T_cells", 
  "7" = "Contamination_MNPs", # DCs or macrophage contamination 
  "8" = "Contamination_MNPs", # DCs or macrophage contamination 
  "9" = "Contamination_mast_cells", 
  "10" = "Contamination_epithelial_cells", 
  "11" = "PCs_PBs", 
  "12" = "Contamination_glial_cells" # not sure about this one
  
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
ggsave(glue("10_broad_annotation/plot/{sample_name}.png"), width = 8, height = 8)

# Export
seurat_obj_singlets_annotated_list[[sample_name]] <- seurat_obj

# Export seurat object list 
saveRDS(seurat_obj_singlets_annotated_list, "10_broad_annotation/out/seurat_obj_singlets_annotated_list.rds")


