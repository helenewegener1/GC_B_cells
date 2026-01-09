# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

source("11_broad_annotation/script/functions.R")

# Load data
# seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")
seurat_obj_singlets_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")

sample_names <- names(seurat_obj_singlets_list)

# Get top 100 DEGs for each sample
lapply(sample_names, function(sample_name) {
  
  # sample_name <- "HH119-SI-PP-CD19-Pool2"
  # sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"
  seurat_obj <- seurat_obj_singlets_list[[sample_name]]
  
  top_DEGs_to_excel(
    seurat_obj = seurat_obj,
    sample_name = sample_name,
    n_DEGs = 100
  )
  
})


