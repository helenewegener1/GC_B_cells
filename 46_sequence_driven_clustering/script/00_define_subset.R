library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

bcr_data <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")

patients <- names(bcr_data)

# ------------------------------------------------------------------------------
# Setsub 1: 2000 PCs per sample and 1000 Navie 
# ------------------------------------------------------------------------------

version <- "subset_1"

subset_1 <- lapply(patients, function(HH){
  
  # HH <- "HH117"
  
  unique(bcr_data[[HH]]$L1_annotation)
  
  # PCs
  ## Extract PC cell IDs
  PC_cell_id <- bcr_data[[HH]] %>% filter(L1_annotation == "PCs") %>% pull(cell_id)
  N_PCs <- length(PC_cell_id)
  
  ## Sample 2000 random PCs
  if (N_PCs > 2000){
    PC_cell_id_sample <- sample(PC_cell_id, 2000)
  } else {
    PC_cell_id_sample <- PC_cell_id
  }
  
  # Naive
  ## Extract PC cell IDs
  naive_cell_id <- bcr_data[[HH]] %>% filter(L1_annotation == "Naive_Bcells") %>% pull(cell_id)
  N_naive <- length(naive_cell_id)
  
  ## Sample 1000 random Naive B cells
  if (N_naive > 1000){
    naive_cell_id_sample <- sample(naive_cell_id, 2000)
  } else {
    naive_cell_id_sample <- naive_cell_id
  }
  
  ## Extract 
  HH_subset <- bcr_data[[HH]] %>% 
    filter(cell_id %in% c(PC_cell_id_sample, naive_cell_id_sample))
  
  # Return 
  return(HH_subset)
  
}) %>% bind_rows()

saveRDS(subset_1, glue("46_sequence_driven_clustering/out/{version}.rds"))

# ------------------------------------------------------------------------------
# Setsub 2: 5000 PCs per sample 
# ------------------------------------------------------------------------------

version <- "subset_2"

n <- 5000

subset_2 <- lapply(patients, function(HH){
  
  # HH <- "HH117"
  
  unique(bcr_data[[HH]]$L1_annotation)
  
  # PCs
  ## Extract PC cell IDs
  PC_cell_id <- bcr_data[[HH]] %>% filter(L1_annotation == "PCs") %>% pull(cell_id)
  N_PCs <- length(PC_cell_id)
  
  # Sample 2000 random PCs
  if (N_PCs > n){
    PC_cell_id_sample <- sample(PC_cell_id, n)
  } else {
    PC_cell_id_sample <- PC_cell_id
  }
  
  ## Extract 
  HH_subset <- bcr_data[[HH]] %>% 
    filter(cell_id %in% PC_cell_id_sample)
  
  # Return 
  return(HH_subset)
  
}) %>% bind_rows()

saveRDS(subset_2, glue("46_sequence_driven_clustering/out/{version}.rds"))



