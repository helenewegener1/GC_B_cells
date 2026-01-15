setwd("~/ciir/people/helweg/projects/GC_B_cells/")

library(tidyverse)
library(Seurat)
library(glue)

# Samples
samples <- list.files("05_run_cellranger/out_v8/") %>% str_split_i("_", 2)

# Prep out list
seurat_obj_list <- list()

for (sample in samples){
  
  # sample <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample <- samples[[1]]
  
  # Define path to cellranger output
  # The files in the per_sample_outs directory have been demultiplexed to single samples.
  # Read more about output of cellrange multi: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-3p-outputs-cellplex
  OUTS_DIR <- glue("05_run_cellranger/out_v8/res_{sample}/outs/per_sample_outs/res_{sample}")
  
  # Read GEX counts (the default method for 10x data)
  gex.data <- Read10X(data.dir = glue("{OUTS_DIR}/count/sample_filtered_feature_bc_matrix"))
  
  # Create the base Seurat object
  
  ######################## If antibody data available ######################## 
  if (length(gex.data) == 2 & "Gene Expression" %in% names(gex.data)){ 
    
    
    # Import Gene Expression (GEX) Data 
    seurat_obj <- CreateSeuratObject(
      counts = gex.data$`Gene Expression`,
      project = "Multi_Analysis", 
      min.cells = 3,
      min.features = 200
    )
    
    # Import Antibody Capture (ADT) Data 
    
    # Extract the ADT count matrix from the list
    adt_counts <- gex.data$`Antibody Capture`
    
    # Get the list of cell barcodes that exist in the GEX-based Seurat object
    gex_cells <- colnames(seurat_obj)
    
    # Subset the ADT matrix to include ONLY those GEX cells
    adt_counts_aligned <- adt_counts[, gex_cells]
    
    # Add ADT data as a new assay (e.g., "ADT")
    seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts_aligned) 
    
    ################### If only gene expression data is available ################### 
  } else { 
    
    # Import Gene Expression (GEX) Data 
    seurat_obj <- CreateSeuratObject(
      counts = gex.data,
      project = "Multi_Analysis", 
      min.cells = 3,
      min.features = 200
    )
    
  }
  
  #################### Import VDJ Annotation Files (TCR and BCR) ##################### 

  # tcr_meta <- NULL
  # bcr_meta <- NULL
  # 
  # # Check if BCR data is available 
  # if ("vdj_t" %in% list.files(OUTS_DIR)){
  #   tcr_meta <- read.csv(glue("{OUTS_DIR}/vdj_t/filtered_contig_annotations.csv"))
  # } 
  # 
  # if ("vdj_b" %in% list.files(OUTS_DIR)){
  #   bcr_meta <- read.csv(glue("{OUTS_DIR}/vdj_b/filtered_contig_annotations.csv"))
  # } 
  # 
  # # Wrangling
  # 
  # # Can we safely remove barcode from contig ids
  # # table(tcr_meta$barcode == str_split_i(tcr_meta$contig_id, "_", 1))
  # # table(bcr_meta$barcode == str_split_i(bcr_meta$contig_id, "_", 1))
  # 
  # # Wrangle TCR data so it can be added to the metadata of the seurat object
  # if (!is.null(tcr_meta)){
  #   
  #   index_barcode <- which(colnames(tcr_meta) == "barcode")
  #   index_contig_id <- which(colnames(tcr_meta) == "contig_id")
  #   
  #   tcr_meta_clean <- tcr_meta %>% 
  #     mutate(contig_id = str_remove_all(contig_id, glue("{barcode}_"))) %>% 
  #     pivot_wider(id_cols = barcode,
  #                 names_from = contig_id, 
  #                 values_from = colnames(tcr_meta)[-c(index_barcode, index_contig_id)]) %>% 
  #     column_to_rownames("barcode")
  #   
  #   # Mutate all "" to NA
  #   tcr_meta_clean[tcr_meta_clean == ""] <- NA
  #   
  #   # # Filter on UMI if 2*TRB?
  #   # tcr_meta_clean %>% select(chain_contig_1, chain_contig_2, umis_contig_1, umis_contig_2) %>% 
  #   #   mutate(same = ifelse(chain_contig_1 == chain_contig_2, TRUE, FALSE)) %>% 
  #   #   filter(same == TRUE)
  #   
  #   # Add prefix to column names as bcr data has the same colnames
  #   colnames(tcr_meta_clean) <- paste0("tcr_", colnames(tcr_meta_clean))
  #   
  #   # Add TCR data to seurat object
  #   seurat_obj <- AddMetaData(
  #     seurat_obj,
  #     metadata = tcr_meta_clean
  #   )
  #   
  # }
  # 
  # # Wrangle BCR data so it can be added to the metadata of the seurat object
  # if (!is.null(bcr_meta)){
  #   
  #   index_barcode <- which(colnames(bcr_meta) == "barcode")
  #   index_contig_id <- which(colnames(bcr_meta) == "contig_id")
  #   
  #   bcr_meta_clean <- bcr_meta %>% 
  #     mutate(contig_id = str_remove_all(contig_id, glue("{barcode}_"))) %>% 
  #     pivot_wider(id_cols = barcode,
  #                 names_from = contig_id, 
  #                 values_from = colnames(bcr_meta)[-c(index_barcode, index_contig_id)]) %>% 
  #     column_to_rownames("barcode")
  #   
  #   # Mutate all "" to NA
  #   bcr_meta_clean[bcr_meta_clean == ""] <- NA
  #   
  #   # Add prefix to column names as tcr data has the same colnames
  #   colnames(bcr_meta_clean) <- paste0("bcr_", colnames(bcr_meta_clean))
  #   
  #   # Add BCR data to seurat object
  #   seurat_obj <- AddMetaData(
  #     seurat_obj,
  #     metadata = bcr_meta_clean
  #   )
  #   
  # }
  

  # CRITICAL STEP: Add clonotype info (the actual T/B cell ID)
  # This usually requires running a specialized function to simplify the contig data.
  # E.g., Use the 'clonotype' and 'raw_clonotype_id' columns from the filtered_contig_annotations.csv
  # and link them to the cell barcodes (BCR-CELL_BARCODE-1, TCR-CELL_BARCODE-1). 
  # Be mindful of the barcode naming convention differences between GEX and VDJ outputs!
  
  # Add seurat object to list 
  seurat_obj_list[[sample]] <- seurat_obj
  
}

######################## Export list of seurat objects ######################### 

saveRDS(seurat_obj_list, "06_seurat_load/out/seurat_obj_list.rds")


