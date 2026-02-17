library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(readxl)
library(tidyr)
library(tibble)
library(purrr)
library(scRepertoire)
library(APackOfTheClones)

seurat_obj <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")
combined.BCR.filtered_all <- readRDS("20_VDJ/out/combined.BCR.filtered.clean_all.rds")

# ------------------------------------------------------------------------------
# Check barcodes 
# ------------------------------------------------------------------------------

table(colnames(seurat_obj) %in% combined.BCR.filtered_all$barcode)

# ------------------------------------------------------------------------------
# Combine 
# ------------------------------------------------------------------------------

# Combine combined.BCR.filtered_sample and seurat_obj
seurat_obj_BCR <- combineExpression(
  combined.BCR.filtered_sample,
  seurat_obj,
  cloneCall = "strict",
  proportion = TRUE
)
  

# ------------------------------------------------------------------------------
# Plot 
# ------------------------------------------------------------------------------

Idents(seurat_obj_BCR) <- "celltype_broad"

UMAPPlot(seurat_obj_BCR) + 
  labs(title = "UMAP", subtitle = sample_name)
ggsave(glue("20_VDJ/plot/BCR_PackOfClones/{sample_name}_UMAP.png"), width = 14, height = 10)

# Fix legend
vizAPOTC(seurat_obj_BCR, verbose = FALSE, legend_text_size = 3.5) + 
  labs(title = "BCR PackOfClones", subtitle = sample_name)
ggsave(glue("20_VDJ/plot/BCR_PackOfClones/{sample_name}_BCR_PackOfClones.png"), width = 14, height = 10)

