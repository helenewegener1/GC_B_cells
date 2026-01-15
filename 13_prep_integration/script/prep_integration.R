# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
seurat_obj_list <- readRDS("12_broad_annotation/out/seurat_obj_singlets_annotated_list.rds")

sample_names <- names(seurat_obj_list)

source("12_broad_annotation/script/color_palette.R")

################################# REMOVE DCs ###################################

n_dim <- 10
res <- 0.1

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH117-SILP-INF-PC"
  
  seurat_obj <- seurat_obj_list[[sample_name]]
  note <- ""
  
  # Remove the DCs if there are any
  if ("DCs_MNPs" %in% seurat_obj$celltype_broad) {
    
    seurat_obj_subset <- subset(seurat_obj, celltype_broad != "DCs_MNPs")
    
    # Seurat workflow after subsetting
    seurat_obj_subset <- NormalizeData(seurat_obj_subset, verbose = FALSE)
    seurat_obj_subset <- FindVariableFeatures(seurat_obj_subset, verbose = FALSE)
    seurat_obj_subset <- ScaleData(seurat_obj_subset, verbose = FALSE)
    seurat_obj_subset <- RunPCA(seurat_obj_subset, verbose = FALSE)
    
    seurat_obj_subset <- FindNeighbors(seurat_obj_subset,  dims = 1:n_dim, verbose = FALSE)
    seurat_obj_subset <- FindClusters(seurat_obj_subset, resolution = res, verbose = FALSE)
    seurat_obj_subset <- RunUMAP(seurat_obj_subset, reduction = "pca", dims = 1:n_dim, verbose = FALSE)
    
    seurat_obj_list[[sample_name]] <- seurat_obj_subset
    
    note <- "DCs removed"
    
  }
  
  n_cells <- ncol(seurat_obj_list[[sample_name]])
  DimPlot(seurat_obj_list[[sample_name]], group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
    labs(subtitle = sample_name, 
         caption = glue("N cells: {n_cells}\n{note}"))
  ggsave(glue("13_prep_integration/plot/{sample_name}.png"), width = 8, height = 8)
  
}

 ################################ ADD META DATA #################################

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
  # condition <- # Define CRC/UC/CD
  # Tissue
  tissue_1 <- str_split_i(sample_name, "-", 2)
  tissue_2 <- ifelse(nchar(tissue_1) < 4, paste0("-", str_split_i(sample_name, "-", 3)), "")
  tissue <- glue("{tissue_1}{tissue_2}")
  # group <- paste(patient, ifelse(inflammed, "INF", "UNINF"), tissue, sep = "_")
  
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
  # seurat_obj@meta.data$group <- group
  seurat_obj@meta.data$group <- sub("-Pool[0-9]+$", "", seurat_obj$sample)

  seurat_obj_prepped_list[[sample_name]] <- seurat_obj
  
}

saveRDS(seurat_obj_prepped_list, "13_prep_integration/out/seurat_obj_prepped_list.rds")


