library(tidyverse)
library(glue)
library(Seurat)
library(SeuratObject)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# df_tcr <- readRDS("70_TCR_analysis/out/df_tcr_filtered.rds")
seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

patients <- seurat_integrated@meta.data$patient %>% unique()

# outdir <- ""
# dir.create()

# ------------------------------------------------------------------------------
# Check expression of ILs across follicles 
# ------------------------------------------------------------------------------

HH <- "HH117"

# Focus on one patient 
seurat_integrated_HH <- subset(seurat_integrated, subset = patient == HH & !is.na(manual_ADT_ID) & !(manual_ADT_ID %in% c("Doublet", "Negative")))

# Reintegrate

seurat_integrated_HH$manual_ADT_ID %>% table()

grep("IL", rownames(seurat_integrated_HH), value = TRUE)
grep("TCF7", rownames(seurat_integrated_HH), value = TRUE)

gene_list <- c()

VlnPlot(seurat_integrated_HH, features = "TCF7", split.by = "manual_ADT_ID", slot = "counts.HH117")











