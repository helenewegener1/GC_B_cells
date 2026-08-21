library(glue)
library(tidyverse)
library(Seurat)
library(scRepertoire)

source("10_broad_annotation/script/color_palette.R")

# ==============================================================================
# NOTES / ASSUMPTIONS
# ==============================================================================
# - Input: 70_TCR_analysis/out/rds/TCR_<CELL_TYPE_FILTER>.rds, written by
#   01_prep_TCR_clone_data.R (one row per Tfh cell, clone id = CTstrict).
# - "Isotype barplot" in the BCR pipeline (60_PC_clones/05_PC_general_figures.R)
#   has no TCR equivalent (T cells don't class-switch). The analog used here
#   is TRBV gene usage per compartment, extracted from CTgene.
# - CTgene's exact string layout depends on your scRepertoire version. Inspect
#   `df_tcr$CTgene %>% head()` once and adjust `extract_trbv()` below if the
#   regex doesn't match.
# ==============================================================================

df_tcr <- readRDS("70_TCR_analysis/out/df_tcr_filtered.rds")

df_tcr$CTgene %>% head()

outdir <- glue("70_TCR_analysis/plot/02_general_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Pull out the dominant TRBV gene call per cell from CTgene
# (expects something like "TRAV12-3.TRAJ23_TRBV20-1.TRBJ2-1"; adjust if needed)
extract_trbv <- function(ctgene) {
  str_extract(ctgene, "TRBV[0-9A-Za-z-]+")
}

df_tcr <- df_tcr %>%
  mutate(v_trb_call = extract_trbv(CTgene))

# ==============================================================================
# N clones per compartment
# ==============================================================================

outdir1 <- glue("{outdir}/N_cells/")
dir.create(outdir1, recursive = TRUE, showWarnings = FALSE)

for (HH in names(patient_names)) {
  
  df_tcr %>%
    filter(patient_id == HH, !is.na(CTstrict)) %>%
    mutate(
      sample_clean_fol_plot = str_remove_all(sample_clean_fol, glue("{HH}-")) %>% str_remove_all("SI-PP_")
    ) %>% 
    count(sample_clean_fol_plot) %>%
    ggplot(aes(x = sample_clean_fol_plot, y = n)) +
    geom_col() +
    geom_text(aes(label = n), size = 3, vjust = -0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = glue("{HH}: N Tfh cells with TCR per compartment"),
      y = "N Tfh cells",
      x = "Compartment"
    ) 
  
  ggsave(glue("{outdir1}/{HH}_N_cells_per_sample.png"), width = 12, height = 6)
  
}


# ==============================================================================
# N clones per compartment
# ==============================================================================

outdir2 <- glue("{outdir}/N_clones/")
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

for (HH in names(patient_names)) {
  
  df_tcr %>%
    filter(patient_id == HH, !is.na(CTstrict)) %>%
    mutate(
      sample_clean_fol_plot = str_remove_all(sample_clean_fol, glue("{HH}-")) %>% str_remove_all("SI-PP_")
    ) %>% 
    select(sample_clean_fol_plot, CTstrict) %>%
    distinct() %>% 
    count(sample_clean_fol_plot) %>%
    ggplot(aes(x = sample_clean_fol_plot, y = n)) +
    geom_col() +
    geom_text(aes(label = n), size = 3, vjust = -0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = glue("{HH}: N clones with TCR per compartment"),
      y = "N clones",
      x = "Compartment"
    ) 
  
  ggsave(glue("{outdir2}/{HH}_N_clones_per_sample.png"), width = 12, height = 6)
  
}


# ==============================================================================
# Clone size frequency plot
# ==============================================================================

outdir3 <- glue("{outdir}/clone_size/")
dir.create(outdir3, recursive = TRUE, showWarnings = FALSE)

df_plot <- df_tcr %>%
  filter(!is.na(CTstrict)) %>%
  count(sample_clean, CTstrict) %>%
  dplyr::rename(clone_size = n) %>%
  mutate(
    clone_size_group = case_when(
      clone_size == 1 ~ "Singleton",
      clone_size == 2 ~ "2",
      clone_size == 3 ~ "3",
      clone_size == 4 ~ "4",
      clone_size == 5 ~ "5",
      clone_size > 5  & clone_size <= 10  ~ "6-10",
      clone_size > 10 & clone_size <= 20  ~ "11-20",
      clone_size > 20 & clone_size <= 50  ~ "21-50",
      clone_size > 50 & clone_size <= 100 ~ "51-100",
      clone_size > 100 ~ "100+"
    ),
    clone_size_group = factor(clone_size_group, levels = c("2", "3", "4", "5", "6-10", "11-20", "21-50", "51-100", "100+"))
  ) %>%
  count(sample_clean, clone_size_group)

df_plot %>%
  filter(!is.na(clone_size_group)) %>%
  ggplot(aes(x = clone_size_group, y = n, color = sample_clean)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_line(aes(group = sample_clean)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Frequency of TCR clone sizes",
    x = "Clone size (N cells)",
    y = "N clones",
    color = "Compartment"
  )

ggsave(glue("{outdir3}/freq_of_clone_size.png"), width = 10, height = 6)

# ==============================================================================
# TRBV gene usage (isotype-barplot analog)
# ==============================================================================

# outdir4 <- glue("{outdir}/gene_usage/")
# dir.create(outdir4, recursive = TRUE, showWarnings = FALSE)
# 

# df_plot <- df_tcr %>%
#   filter(!is.na(CTstrict), !is.na(v_trb_call)) %>%
#   count(sample_clean, v_trb_call)
# 
# df_plot %>%
#   ggplot(aes(x = sample_clean, y = n, fill = v_trb_call)) +
#   geom_col() +
#   scale_fill_viridis_d(option = "turbo") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(
#     title = "TRBV gene usage",
#     x = "Compartment",
#     y = "N cells",
#     fill = "TRBV gene"
#   )
# 
# ggsave(glue("{outdir3}/TRBV_usage_barplot.png"), width = 12, height = 6.5)

# ==============================================================================
# APackOfTheClones - clone-size UMAP overlay
# ==============================================================================
# Requires the APackOfTheClones package and a Seurat object subset to the
# same cells as df_tcr (join on cell_id, same pattern as
# 60_PC_clones/05_PC_general_figures.R). Left commented out since it depends
# on SEURAT_OBJ_PATH from 01_prep_TCR_clone_data.R being correctly wired up
# -- uncomment once that join has been verified.

# library(APackOfTheClones)
#
# outdir4 <- glue("{outdir}/APackOfTheClones/")
# dir.create(outdir4, recursive = TRUE, showWarnings = FALSE)
#
# seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")
#
# for (HH in unique(df_tcr$patient_id)) {
#
#   # HH <- "HH117"
#   df_HH <- df_tcr %>% filter(patient_id == HH)
#
#   seurat_HH <- subset(seurat_integrated, cells = intersect(colnames(seurat_integrated), df_HH$cell_id))
#   seurat_HH[[]] <- seurat_HH[[]] %>%
#     rownames_to_column("cell_id") %>%
#     left_join(df_HH, by = "cell_id") %>%
#     column_to_rownames("cell_id")
#
#   Idents(seurat_HH) <- "sample_clean_fol"
#   vizAPOTC(
#     seurat_HH,
#     reduction_base = "umap.unintegrated",
#     clonecall = "CTstrict",
#     show_labels = TRUE,
#     legend_text_size = 3,
#     label_size = 3
#   ) +
#     plot_annotation(title = glue("{HH}: {CELL_TYPE_FILTER} clone sizes"))
#
#   ggsave(glue("{outdir4}/{HH}_APackOfTheClones.png"), width = 7, height = 6.5)
#
# }
