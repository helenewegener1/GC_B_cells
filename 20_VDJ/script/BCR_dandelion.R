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
combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.rds")
seurat_object_list <- readRDS("10_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds")

# Make seurat object in same format as combined.BCR.filtered
sample_names <- names(combined.BCR.filtered) %>% str_split_i("_", 1) %>% unique()

seurat_object_split <- list()
for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  seurat_object <- seurat_object_list[[sample_name]]
  
  if ("manual_ADT_ID" %in% colnames(seurat_object[[]])){
    
    seurat_object$split_by <- paste(sample_name, seurat_object$manual_ADT_ID, sep = "_")
    seurat_object_split <- c(seurat_object_split, SplitObject(seurat_object, split.by = "split_by"))
    # names(seurat_object_split)
  
  } else {
    
    seurat_object_split <- c(seurat_object_split, setNames(list(seurat_object), sample_name))
    
  }
  
}

# Check that names are the same 
names(combined.BCR.filtered) %>% length()
names(seurat_object_split) %>% length()

table(names(combined.BCR.filtered) == names(seurat_object_split))

# Make into single cell experiment object 
seurat_object_merged <- merge(seurat_object_split[[1]], y = seurat_object_split[-1])
sce <- as.SingleCellExperiment(seurat_object_merged)
# sce <-

sce <- combineExpression(combined.TCR, demo_sce)
