# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
# seurat_obj_list <- readRDS("11_broad_annotation/out/")
seurat_obj_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list.rds")

# MERGE THE noDC OBJECT AND THE SINGLETS AND HAVE THEM BROADLY ANNOTATED

sample_names <- names(seurat_obj_list) 

# Prep to save seurat objects
seurat_obj_prepped_list <- rep(0, length(seurat_obj_list)) %>% as.list()
names(seurat_obj_prepped_list) <- names(seurat_obj_list)

# Add metadata
for (sample_name in sample_names){
  
  # sample_name <- "HH119-SILP-PC"
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # Define metadata
  sample <- sample_name
  patient <- str_split_i(sample_name, "-", 1)
  inflammed <- str_detect(sample_name, "-INF-")
  # Tissue
  tissue_1 <- str_split_i(sample_name, "-", 2)
  tissue_2 <- ifelse(nchar(tissue_1) < 4, paste0("-", str_split_i(sample_name, "-", 3)), "")
  tissue <- glue("{tissue_1}{tissue_2}")
  
  # print(sample)
  # print(patient)
  # print(inflammed)
  # print(tissue)
  # print("------------------")
  
  # Add metadata to seurat object
  seurat_obj@meta.data$sample <- sample
  seurat_obj@meta.data$patient <- patient
  seurat_obj@meta.data$inflammed <- inflammed
  seurat_obj@meta.data$tissue <- tissue

  seurat_obj_prepped_list[[sample_name]] <- seurat_obj
  
}

saveRDS(seurat_obj_prepped_list, "12_prep_integration/out/seurat_obj_prepped_list.rds")
