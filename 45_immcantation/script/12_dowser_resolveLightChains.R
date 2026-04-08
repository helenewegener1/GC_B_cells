setwd("~/ciir/people/helweg/projects/GC_B_cells")

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

clone_10x_combined <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds") %>% bind_rows()

light_chain_list <- lapply(sample_names, function(x){
  
  # x <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"
  
  light_chain <- read.delim(glue("45_immcantation/out/{x}/{x}_light_germ-pass_QC.tsv"))
  
  # light_chain$sample_id <- x %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  light_chain$cell_id <- paste(light_chain$sample_id, light_chain$cell_id, sep = "_")
  light_chain$sequence_id <- paste(light_chain$sequence_id, "Light", sep = "_")
  
  light_chain <- light_chain %>%  mutate(
    sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
  )
  
  return(light_chain)
  
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

both_combined_all <- bind_rows(clone_10x_combined, light_chain_combined)

both_combined_all$cell_id %>% str_split_i("_", 2) %>% unique()

both_combined <- list(
  "HH117" = both_combined_all %>% filter(subject_id == "HH117"),
  "HH119" = both_combined_all %>% filter(subject_id == "HH119")
)

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
both_combined$HH119$sequence_id %>% length()
both_combined$HH119$sequence_id %>% unique() %>% length()
clone_10x_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% length()
clone_10x_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% unique() %>% length()
light_chain_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% length()
light_chain_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% unique() %>% length()
# light_chain_combined %>% filter(subject_id == "HH119") %>% count(cell_id, sort = TRUE)

both_combined$HH117$sequence_id %>% length()
both_combined$HH117$sequence_id %>% unique() %>% length()
clone_10x_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% length()
clone_10x_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% unique() %>% length()
light_chain_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% length()
light_chain_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% unique() %>% length()

# Max 2 chains per cell id: heavy and light chain
both_combined$HH119 %>% count(cell_id, sort = TRUE) %>% pull(n) %>% table()
both_combined$HH117 %>% count(cell_id, sort = TRUE) %>% pull(n) %>% table()

HH117_heavy <- clone_10x_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) 
HH117_light <- light_chain_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) 

HH117_heavy %>% sort() %>% head()
HH117_light %>% sort() %>% head()

table(HH117_heavy %in% HH117_light)

HH119_heavy <- clone_10x_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) 
HH119_light <- light_chain_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) 

HH119_heavy %>% sort() %>% head()
HH119_light %>% sort() %>% head()

table(HH119_heavy %in% HH119_light)

# Run resolveLightChains
patients <- names(both_combined)

resolve_LC_list <- lapply(patients, function(HH){

  # HH <- "HH119"
  
  # both_combined[[HH]] %>% filter(locus == "IGH") %>% nrow()
  
  # resolve_LC_HH <- resolve_LC_list[[HH]]
  resolve_LC_HH <- resolveLightChains(both_combined[[HH]])
  
  # resolve_LC_HH %>% filter(locus == "IGH") %>% nrow()

  table(resolve_LC_HH$celltype_broad, useNA = "always")
  
  # Get meta data from heavy chains
  meta <- resolve_LC_HH %>%
    filter(!is.na(celltype_broad)) %>%
    select(cell_id, celltype_broad) %>%
    distinct()

  # Clear celltype
  resolve_LC_HH$celltype_broad <- NULL

  # Add metadata
  resolve_LC_HH <- resolve_LC_HH %>% left_join(meta, by = "cell_id")

  return(resolve_LC_HH)

}) %>% setNames(patients)

saveRDS(resolve_LC_list, "45_immcantation/out/rds/resolve_LC_list.rds")
