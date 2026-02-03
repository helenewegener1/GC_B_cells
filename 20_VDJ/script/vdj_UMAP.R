getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)

# Load data
seurat_obj_list <- readRDS("08_seurat_QC/out/seurat_obj_finalQC_list.rds")

# Which samples have TCR and BCR data 
tcr_mask <- lapply(names(seurat_obj_list), function(x) {"tcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()
bcr_mask <- lapply(names(seurat_obj_list), function(x) {"bcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with TCR and BCR data 
tcr_seurat_obj_list <- seurat_obj_list[tcr_mask]
bcr_seurat_obj_list <- seurat_obj_list[bcr_mask]

# Colnums to go through 
gene_letters <- c("c", "v", "j", "d")

# Make DimPlots for TCR 
for (sample_name in names(tcr_seurat_obj_list)){
  
  # Create directory for plots of specific sample
  out_dir <- glue("12_vdj/plot/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE)
  
  # sample_name <- names(tcr_seurat_obj_list)[[1]]
  
  # Define specific seurat object 
  tcr_seurat_obj <- tcr_seurat_obj_list[[sample_name]]
  
  for (gene_letter in gene_letters){

    # Dimplot for each contig
    for (contig_nr in c(1, 2, 3, 4)){
      
      # Define group: gene and contig combination 
      group <- glue("tcr_{gene_letter}_gene_contig_{contig_nr}")
      
      # Check group is in seurat object (not all contigs for all genes are in all samples)
      if (!group %in% colnames(tcr_seurat_obj[[]])){
        next
      } else {

        p <- DimPlot(tcr_seurat_obj, reduction = "umap", group.by = group) +
          labs(title = glue("TCR {gene_letter} gene, contig {contig_nr}"), 
               subtitle = sample_name) + 
          theme(legend.text = element_text(size = 8))
        
        # Save in different format if a lot of unique groups 
        unique_groups <- length(unique(tcr_seurat_obj@meta.data[[group]]))
        if (unique_groups > 20) {
          ggsave(glue("{out_dir}/{sample_name}_DimPlot_tcr_{gene_letter}_gene_contig_{contig_nr}.png"), width = 14, height = 10)
        } else {
          ggsave(glue("{out_dir}/{sample_name}_DimPlot_tcr_{gene_letter}_gene_contig_{contig_nr}.png"), width = 7, height = 6)
        }

      }
      
    }

  }
  
}

# Make DimPlots for BCR (copy of the TCR workflow )
for (sample_name in names(bcr_seurat_obj_list)){
  
  # Create directory for plots of specific sample
  out_dir <- glue("12_vdj/plot/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE)
  
  # sample_name <- names(bcr_seurat_obj_list)[[1]]
  
  # Define specific seurat object 
  bcr_seurat_obj <- bcr_seurat_obj_list[[sample_name]]
  
  for (gene_letter in gene_letters){
    
    # Dimplot for each contig
    for (contig_nr in c(1, 2, 3, 4)){
      
      # Define group: gene and contig combination 
      group <- glue("bcr_{gene_letter}_gene_contig_{contig_nr}")
      
      # Check group is in seurat object (not all contigs for all genes are in all samples)
      if (!group %in% colnames(bcr_seurat_obj[[]])){
        next
      } else {
        
        p <- DimPlot(bcr_seurat_obj, reduction = "umap", group.by = group) +
          labs(title = glue("TCR {gene_letter} gene, contig {contig_nr}"), 
               subtitle = sample_name) + 
          theme(legend.text = element_text(size = 8))
        
        # Save in different format if a lot of unique groups 
        unique_groups <- length(unique(bcr_seurat_obj@meta.data[[group]]))
        if (unique_groups > 20) {
          ggsave(glue("{out_dir}/{sample_name}_DimPlot_bcr_{gene_letter}_gene_contig_{contig_nr}.png"), width = 14, height = 10)
        } else {
          ggsave(glue("{out_dir}/{sample_name}_DimPlot_bcr_{gene_letter}_gene_contig_{contig_nr}.png"), width = 7, height = 6)
        }
        
      }
      
    }
    
  }
  
}



