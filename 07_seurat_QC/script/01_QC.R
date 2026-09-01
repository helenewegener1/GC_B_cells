setwd("~/gcb/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(glmGamPoi)
# devtools::install_github("constantAmateur/SoupX", ref='devel')
library(SoupX)
library(multtest)
library(scDblFinder)
library(decontX)
library(readxl)
library(purrr)
library(patchwork)
library(png)
library(grid)

version <- "v9"

# Load data
seurat_obj_list <- readRDS(glue("06_seurat_load/out/seurat_obj_list_{version}.rds")) # cellranger filtered

# Load Gina annotation file 
broad_annot_file <- read_excel("00_data/Gene_markers_GL_HW_new.xlsx", sheet = "Very broad level")
detailed_annot_file <- read_excel("00_data/Gene_markers_GL_HW_new.xlsx", sheet = "More detailed level")

# Extract cell marker genes as lists
broad_markers <- broad_annot_file %>% as.list()
broad_markers <- map(broad_markers, ~ .x[!is.na(.x)])
names(broad_markers) <- c("T_cell", "B_cell", "DC", "plasmablast_plasma_cell")

detailed_markers <- detailed_annot_file %>% as.list() %>% na.omit()
detailed_markers <- map(detailed_markers, ~ .x[!is.na(.x)])
names(detailed_markers) <- c("TFH_cell", "Naive_B_cell", "Memory_B_cell", "GC_B_cell", 
                             "Activation_markers", "Immunoglobulin_subsets", "Other_markers")

# Load sample to get genes that are in the data.  
seurat_obj <- seurat_obj_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

source("09_seurat_QC_clusters/script/functions.R")

# Update marker format for the broad markers
broad_markers <- update_marker_names(broad_markers, seurat_obj)

########################### Define major cell types ############################

# Prep to save clustered seurat objects
seurat_obj_clustered_list <- rep(0, length(seurat_obj_list)) %>% as.list()
names(seurat_obj_clustered_list) <- names(seurat_obj_list)

# N PCs
# n_pcs <- list(
#   "HH117-SILP-INF-PC" = ,                        
#   "HH117-SILP-nonINF-PC" = ,                           
#   "HH117-SI-MILF-INF-HLADR-AND-CD19" = ,                 
#   "HH117-SI-MILF-nonINF-HLADR-AND-CD19" = ,             
#   "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = ,
#   "HH119-COLP-PC" = ,                                   
#   "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH" = ,    
#   "HH119-SILP-PC" = ,                                   
#   "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH" = ,      
#   "HH119-SI-PP-CD19-Pool1" = ,                        
#   "HH119-SI-PP-CD19-Pool2" = ,                        
#   "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = ,                
#   "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = 
# )

# sample_names <- names(seurat_obj_list)
seurat_obj_clustered_list_old <- readRDS(glue("07_seurat_QC/out/seurat_obj_clustered_list_{version}.rds"))

sample_names <- names(seurat_obj_list)[!(names(seurat_obj_list) %in% names(seurat_obj_clustered_list_old))]

res <- 0.1
n_dims <- 10

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-COLP-PC"
  # sample_name <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Green"  
  
  print(glue("--- Processing: {sample_name} ---"))
  
  # Get seurat object 
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # Create directory for plots of specific sample
  out_dir <- glue("07_seurat_QC/plot_{version}/01_clusters/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Seurat workflow so I can UMAP
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  
  ElbowPlot(seurat_obj) + labs(title = sample_name) #to determine dimentions used for following steps in doublet detection. Adjust dims. 
  ggsave(glue("{out_dir}/{sample_name}_elbow.png"), width = 9, height = 5.5)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims, verbose = FALSE)

  n_cells <- ncol(seurat_obj)

  # Plot
  DimPlot(seurat_obj, reduction = 'umap', label = TRUE) + NoLegend() + 
    labs(title = glue("Seurat clusters {version}"),
         subtitle = sample_name, 
         caption = glue("N cells: {n_cells}\nN dim: {n_dims}\nresolution: {res}"))
  ggsave(glue("{out_dir}/{sample_name}_clusters.png"), width = 8, height = 8)
  
  ####################### FeaturePlot with broad_markers ####################### 
  for (markers in names(broad_markers)){
    
    FeaturePlot(seurat_obj, features = broad_markers[[markers]], ncol = 3) + 
      plot_annotation(title = glue("{markers} - {version}"),
                      subtitle = sample_name,
                      caption = glue("N cells: {n_cells}\nN dim: {n_dims}\nresolution: {res}"))
    
    # Adjust heigh of plot to number of markers
    n_markers <- broad_markers[[markers]] %>% length()
    height <- (n_markers/3) * 4
    
    ggsave(glue("{out_dir}/{sample_name}_broad_{markers}.png"), width = 14, height = height)
    
  }
  
  ############################################################################## 
  
  ##################### FeaturePlot with detailed_markers ######################
  for (markers in names(detailed_markers)){
    
    FeaturePlot(seurat_obj, features = detailed_markers[[markers]], ncol = 3) + 
      plot_annotation(title = glue("{markers} - {version}"),
                      subtitle = sample_name,
                      caption = glue("N cells: {n_cells}\nN dim: {n_dims}\nresolution: {res}"))
    
    # Adjust heigh of plot to number of markers
    n_markers <- detailed_markers[[markers]] %>% length()
    height <- (n_markers/3) * 4
    
    ggsave(glue("{out_dir}/{sample_name}_detailed_{markers}.png"), width = 14, height = height)
    
  }
  
  ##############################################################################
  
  # Save object
  seurat_obj_clustered_list[[sample_name]] <- seurat_obj
  
}


# Combine 
seurat_obj_clustered_list_combined <- c(seurat_obj_clustered_list_old, seurat_obj_clustered_list)

saveRDS(seurat_obj_clustered_list_combined, glue("07_seurat_QC/out/seurat_obj_clustered_list_{version}.rds"))

# Get number of clusters
for (sample_name in sample_names){

  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  n_clusters <- seurat_obj_clustered_list[[sample_name]]$seurat_clusters %>% levels() %>% length()

  print(sample_name)
  print(glue("N clusters: {n_clusters}"))
  print("---------------------------------------------")

}

################################################################################
############################ GINA CLUSTER MERGING ##############################
################################################################################

# Gina got the plots from 07_seurat_QC/plot/clusters.
# Gina drew which clusters should be merged based on expression.
# This was to make sure that doublet finder did not merge two too similar cells and called them a doublet.
# Gina's drawing can be seen in 07_seurat_QC/plot/gina_merged_clusters.

