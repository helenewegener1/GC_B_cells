library(glue)
library(tidyverse)
library(alakazam)
library(scatterpie)
library(patchwork)
library(ggbreak)

source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

L3_GCB_annotation <- readRDS("00_data/GCB_meta_GL.rds")

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

# Prep output
outdir = glue("45_immcantation/plot/15_clones_L3_GCB_annotation/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
    filter(
      locus == "IGH"
      # !is.na(manual_ADT_full_ID)
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows()

# Define largest CRC clone
large_crc_clone <- df_both %>%
  filter(patient_id == "HH119", L1_annotation == "GC_B_cells") %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  head(1) %>% 
  pull(clone_subgroup_id_90_similarity)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

df_both$cell_id
L3_GCB_annotation$cell_id

L3_GCB_annotation_clean <- L3_GCB_annotation %>% 
  select(cell_id, Initial_annotation_comb, L2_GCB_annotation, L3_GCB_annotation)

df_both <- df_both %>% left_join(L3_GCB_annotation_clean, by = "cell_id")

n_clones <- 20

for (HH in patients) {
  
  # HH <- "HH117"
  
  top_clones <- df_both %>% 
    filter(
      patient_id == HH, 
      L1_annotation == "GC_B_cells"
    ) %>% 
    dplyr::count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
    head(n_clones) %>% 
    pull(clone_subgroup_id_90_similarity)
  
  df_plot <- df_both %>%
    filter(
      patient_id == HH, 
      L1_annotation == "GC_B_cells",
      clone_subgroup_id_90_similarity %in% top_clones
    ) %>% 
    count(clone_subgroup_id_90_similarity, L3_GCB_annotation) %>% 
    mutate(
      clone_subgroup_id_90_similarity = factor(clone_subgroup_id_90_similarity, top_clones)
    )
  
  df_plot %>% 
    ggplot(aes(x = clone_subgroup_id_90_similarity, y = n, fill = L3_GCB_annotation)) +
    geom_col(position = "fill") + 
    scale_fill_manual(values = L3_GC_pretty_colors) + 
    scale_y_continuous(
      labels = scales::percent,
      breaks = scales::breaks_width(0.1)
    ) +
    theme_bw() + 
    labs(
      y = "Count", 
      x = "Clone", 
      fill = "GCB annotation", 
      title = glue("{HH}: GCB annotation of top {n_clones} clones")
    )
  
  ggsave(glue("{outdir}/{HH}_top{n_clones}_L3_GCB_annotation.png"), width = 13)
  

}



