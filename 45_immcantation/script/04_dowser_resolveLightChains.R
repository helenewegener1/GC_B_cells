setwd("~/gcb/")

library(dowser)
library(alakazam)
library(dplyr)
library(tidyverse)
library(shazam)
library(glue)

# Following: https://dowser.readthedocs.io/en/stable/vignettes/Germlines-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

files <- list.files("45_immcantation/out/")
sample_names <- files[1:length(files)-1]

# # Load light chain corrected
# clone_10x_list <- lapply(sample_names, function(x){
# 
#   # x <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"
# 
#   # These files are created using light_clusters.py via 06_add_light_chain. 
#   # We are not doing this anymore. See 07_summarise_clones.R
#   clone_10x <- read.delim(glue("45_immcantation/out/{x}/{x}_10X_clone-pass.tsv"))
#   
#   clone_10x$sample_id <- x %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
#   clone_10x$cell_id <- paste(clone_10x$sample_id, clone_10x$cell_id, sep = "_")
#   clone_10x$sequence_id <- paste(clone_10x$sequence_id, "Heavy", sep = "_")
# 
#   return(clone_10x)
# 
# }) %>% setNames(sample_names)
# clone_10x_combined <- bind_rows(clone_10x_list)

# clone_10x_combined <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds") %>% bind_rows()
# version <- ""


# data <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_gmm_threshold.rds") 
data <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds") 
clone_10x_combined <- data %>% bind_rows()

# HH <- "HH117"
# clone_10x_combined %>% filter(patient_id == HH) %>% nrow()

light_chain_list <- lapply(sample_names, function(x){
  
  # x <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"
  # x <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Yellow"
  
  light_chain <- read.delim(glue("45_immcantation/out/{x}/{x}_light_germ-pass_QC.tsv"))
  # light_chain <- read.delim(glue("45_immcantation/out/{x}/{x}_light_germ-pass.tsv"))
  
  # light_chain$sample_id <- x %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  light_chain$cell_id <- paste(x, light_chain$cell_id, sep = "_")
  light_chain$sequence_id <- paste(light_chain$sequence_id, "Light", sep = "_")
  
  # Heavy chain cells 
  # meta_x <- clone_10x_combined %>% filter(sample_id == x) %>% 
  #   select(
  #     cell_id
  #   )
  # heavy_cells <- meta_x$cell_id
  
  nrow(light_chain)
  
  light_chain_qc <- light_chain %>% 
    # inner_join(meta_x, by = "cell_id") %>% 
    mutate(
      sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
    ) %>% 
    select(-L1_annotation)
  
  nrow(light_chain_qc)
  
  return(light_chain_qc)
  
}) %>% setNames(sample_names)

light_chain_combined <- bind_rows(light_chain_list)

# Check cell IDs of heavy and light chain
cell_id_heavy <- clone_10x_combined$cell_id %>% str_split_i("_", 1) %>% unique() %>% sort()
cell_id_light <- light_chain_combined$cell_id %>% str_split_i("_", 1) %>% unique() %>% sort()
table(cell_id_heavy == cell_id_light)

table(light_chain_combined$cell_id %in% clone_10x_combined$cell_id)
table(clone_10x_combined$cell_id %in% light_chain_combined$cell_id)

# Combine heavy and light chain in one df
light_chain_combined$barcode_suffix <- as.character(light_chain_combined$barcode_suffix)

all_combined <- bind_rows(clone_10x_combined, light_chain_combined)

all_combined$cell_id %>% str_split_i("_", 2) %>% unique()

patients <- all_combined$patient_id %>% unique()

all_combined_list <- lapply(patients, function(HH){
  all_combined %>% filter(patient_id == HH)
}) %>% setNames(patients)

saveRDS(all_combined_list, "45_immcantation/out/rds/04_bcr_heavy_light.rds")

# spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")
# 
# HH <- "HH119"
# HH_spec_clones_vj <- spec_clones_vj[[HH]]
# 
# nrow(HH_spec_clones_vj)
# 
# # Check that neccessary columns are present
# HH_spec_clones_vj$sequence_alignment %>% head()
# HH_spec_clones_vj$germline_alignment_d_mask %>% head()

# ------------------------------------------------------------------------------
# resolveLightChains
# ------------------------------------------------------------------------------

# Check that sequence_id and cell_id are unique
for (HH in patients){
  # HH <- "HH151"
  print(HH)
  print(all_combined_list[[HH]]$sequence_id %>% length())
  print(all_combined_list[[HH]]$sequence_id %>% unique() %>% length())
  print(clone_10x_combined %>% filter(patient_id == HH) %>% pull(cell_id) %>% length())
  print(clone_10x_combined %>% filter(patient_id == HH) %>% pull(cell_id) %>% unique() %>% length())
  print(light_chain_combined %>% filter(patient_id == HH) %>% pull(cell_id) %>% length())
  print(light_chain_combined %>% filter(patient_id == HH) %>% pull(cell_id) %>% unique() %>% length())
  cat("\n")
}

# # light_chain_combined %>% filter(patient_id == "HH119") %>% count(cell_id, sort = TRUE)
# 
# all_combined_list$HH117$sequence_id %>% length()
# all_combined_list$HH117$sequence_id %>% unique() %>% length()
# clone_10x_combined %>% filter(patient_id == "HH117") %>% pull(cell_id) %>% length()
# clone_10x_combined %>% filter(patient_id == "HH117") %>% pull(cell_id) %>% unique() %>% length()
# light_chain_combined %>% filter(patient_id == "HH117") %>% pull(cell_id) %>% length()
# light_chain_combined %>% filter(patient_id == "HH117") %>% pull(cell_id) %>% unique() %>% length()

# Max 2 chains per cell id: heavy and light chain
for (HH in patients){
  print(HH)
  print(all_combined_list[[HH]] %>% dplyr::count(cell_id, sort = TRUE) %>% pull(n) %>% table())
}
  
HH117_heavy <- clone_10x_combined %>% filter(patient_id == "HH117") %>% pull(cell_id) 
HH117_light <- light_chain_combined %>% filter(patient_id == "HH117") %>% pull(cell_id) 

HH117_heavy %>% sort() %>% head()
HH117_light %>% sort() %>% head()

table(HH117_heavy %in% HH117_light)

HH119_heavy <- clone_10x_combined %>% filter(patient_id == "HH119") %>% pull(cell_id) 
HH119_light <- light_chain_combined %>% filter(patient_id == "HH119") %>% pull(cell_id) 

HH119_heavy %>% sort() %>% head()
HH119_light %>% sort() %>% head()

table(HH119_heavy %in% HH119_light)

# Run resolveLightChains
# patients <- names(all_combined_list)

# resolve_LC_list <- lapply(patients, function(HH){
# 
#   # HH <- "HH117"
#   
#   # all_combined_list[[HH]] %>% filter(locus == "IGH") %>% nrow()
#   
#   # resolve_LC_HH <- resolve_LC_list[[HH]]
#   resolve_LC_HH <- resolveLightChains(all_combined_list[[HH]])
#   
#   # resolve_LC_HH %>% filter(locus == "IGH") %>% nrow()
# 
#   # table(resolve_LC_HH$celltype_broad, useNA = "always")
#   
#   # Get meta data from heavy chains
#   meta <- resolve_LC_HH %>%
#     filter(locus == "IGH", !is.na(celltype_broad)) %>%
#     select(cell_id, celltype_broad)
#   
#   # nrow(meta)
#   # table(resolve_LC_HH$cell_id %in% meta$cell_id)
# 
#   # Clear and re-add 
#   resolve_LC_HH_final <- resolve_LC_HH %>%
#     select(-celltype_broad) %>%
#     left_join(meta, by = "cell_id")
#   
#   # resolve_LC_HH_test %>% filter(locus == "IGH") %>% nrow()
# 
#   return(resolve_LC_HH_final)
# 
# }) %>% setNames(patients)
# 
# saveRDS(resolve_LC_list, "45_immcantation/out/rds/resolve_LC_list.rds")
# 
# ## Test
# HH <- "HH117"
# clone_10x_combined %>% filter(patient_id == HH) %>% nrow()
# resolve_LC_list[[HH]] %>% filter(locus == "IGH") %>% nrow()
# 
# HH <- "HH119"
# clone_10x_combined %>% filter(patient_id == HH) %>% nrow()
# resolve_LC_list[[HH]] %>% filter(locus == "IGH") %>% nrow()




