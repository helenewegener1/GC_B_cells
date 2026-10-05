library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

L3_GCB_annotation <- readRDS("00_data/GCB_meta_GL.rds")

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

# Load both patients
df_all <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
    filter(
      locus == "IGH"
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows()

# Define largest CRC clone
large_crc_clone <- df_all %>%
  filter(patient_id == "HH119", L1_annotation == "GC_B_cells") %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  head(1) %>% 
  pull(clone_subgroup_id_90_similarity)

# Create outdit
outdir <- glue("50_BACH2/plot/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Wrangle data
# ------------------------------------------------------------------------------

df_all$cell_id
L3_GCB_annotation$cell_id

L3_GCB_annotation_clean <- L3_GCB_annotation %>% 
  select(cell_id, Initial_annotation_comb, L2_GCB_annotation, L3_GCB_annotation)

df_all <- df_all %>% left_join(L3_GCB_annotation_clean, by = "cell_id")


# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

# Overview
df_all_GC <- df_all %>% filter(L1_annotation == "GC_B_cells")
table(df_all_GC$patient, df_all_GC$L3_GCB_annotation)

table(df_all_GC$c_call, df_all_GC$L3_GCB_annotation, useNA = "always")

# Largest clone
df_all %>% 
  filter(patient_id == "HH119", clone_subgroup_id_90_similarity == large_crc_clone) %>% 
  pull(L3_GCB_annotation) %>% 
  table()


# 03_heavy_bcr_data_qc_annot.rds
data <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds") 

HH <- "HH119"
GC_B_cells <- L3_GCB_annotation_clean %>% filter(!is.na(L3_GCB_annotation), L3_GCB_annotation == "BACH2+")
df_heavy <- data[[HH]] %>% filter(locus == "IGH") 

table(GC_B_cells$cell_id %in% df_heavy$cell_id)

cell_id_heavy <- GC_B_cells$cell_id[GC_B_cells$cell_id %in% df_heavy$cell_id]
length(cell_id_heavy)




# 04_bcr_heavy_light.rds
all_combined_list <- readRDS("45_immcantation/out/rds/04_bcr_heavy_light.rds")

HH <- "HH119"
GC_B_cells <- L3_GCB_annotation_clean %>% filter(!is.na(L3_GCB_annotation), L3_GCB_annotation == "BACH2+")
df_heavy <- all_combined_list[[HH]] %>% filter(locus == "IGH") 
df_light <- all_combined_list[[HH]] %>% filter(locus != "IGH") 

table(GC_B_cells$cell_id %in% df_heavy$cell_id)
table(GC_B_cells$cell_id %in% df_light$cell_id)

cell_id_heavy <- GC_B_cells$cell_id[GC_B_cells$cell_id %in% df_heavy$cell_id]
cell_id_light <- GC_B_cells$cell_id[GC_B_cells$cell_id %in% df_light$cell_id]

length(cell_id_heavy)
intersect(cell_id_heavy, cell_id_light) %>% length()

# df <- all_combined_list$HH119 %>% left_join(L3_GCB_annotation_clean, by = "cell_id")



  