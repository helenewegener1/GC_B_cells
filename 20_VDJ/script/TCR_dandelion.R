# Following this tutorial:
# https://bioconductor.org/packages/release/bioc/vignettes/dandelionR/inst/doc/dandelionR.html

# devtools::install_github("CRAN/smoother")
# devtools::install_github("tuonglab/dandelionR", dependencies = TRUE)
library(dandelionR) 
library(scRepertoire)
library(scater)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

# contig.list <- loadContigs(input = demo_airr, format = "AIRR")
# 
# # Format to `scRepertoire`'s requirements and some light filtering
# combined.TCR <- combineTCR(contig.list,
#                            removeNA = TRUE,
#                            removeMulti = FALSE,
#                            filterMulti = TRUE
# )

# Load data 
combined.TCR.filtered <- readRDS("20_VDJ/out/combined.TCR.filtered.rds")
seurat_object_list_ADT <- readRDS("10_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds")
seurat_object_list_annotated <- readRDS("12_broad_annotation/out/seurat_obj_singlets_annotated_list.rds")

# Make seurat object in same format as combined.BCR.filtered
sample_names <- names(combined.TCR.filtered) %>% str_split_i("_", 1) %>% unique()

seurat_object_split <- list()
for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  seurat_object <- seurat_object_list[[sample_name]]
  
  # Combine combined.TCR.filtered per sample
  combined.TCR_combined <- bind_rows(combined.TCR.filtered[grep(sample_name, names(combined.TCR.filtered))])

  # Match barcode of seurat obejct to those of combined.TCR
  colnames(seurat_object) <- paste(sample_name, seurat_object[[]]$manual_ADT_ID, colnames(seurat_object), sep = "_")
  
  # Exclude Doublets and Negatives 
  combined.TCR_combined_clean <- combined.TCR_combined[grep("Negative|Doublet", combined.TCR_combined$barcode, invert = TRUE) ,]
  seurat_object_clean <- subset(seurat_object, grep("Negative|Doublet", seurat_object$manual_ADT_ID, invert = TRUE))
  
  # Make single cell experiment
  sub_sce <- as.SingleCellExperiment(seurat_object_clean)
  
  # Merging VDJ data with gene expression data
  sce <- combineExpression(combined.TCR_combined_clean, sub_sce)
  
  # Dandelion workflow 
  sce <- setupVdjPseudobulk(sce,
                            mode_option = "abT",
                            already.productive = TRUE,
                            subsetby = "manual_ADT_ID")
  
  plotUMAP(sce, color_by = "manual_ADT_ID") 
  
  # -----------------------------------
  # Milo object and neighbourhood graph construction
  # -----------------------------------
  
  library(miloR)
  milo_object <- Milo(sce)
  milo_object <- buildGraph(milo_object, k = 30, d = 20, reduced.dim = "X_scvi")
  milo_object <- makeNhoods(milo_object,
                            reduced_dims = "X_scvi", d = 20,
                            prop = 0.3
  )

  
}

