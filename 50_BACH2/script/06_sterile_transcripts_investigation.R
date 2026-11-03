library(Seurat)
library(tidyverse)
library(glue)
source("10_broad_annotation/script/color_palette.R")

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

seurat_obj <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")
df_sterile_cells <- readRDS("50_BACH2/out/df_sterile_cells.rds")
broad_annotaions <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs_metadata.rds")
GC_annotaions <- readRDS("00_data/GCB_meta_GL_new.rds")

outdir <- "50_BACH2/plot/06_sterile_transcripts"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Prep data
# ------------------------------------------------------------------------------

# Wrange cell id 
## Seurat 
seurat_obj$cell_id <- colnames(seurat_obj) %>% str_remove("_\\d+")
seurat_obj$cell_id <- glue("{seurat_obj$sample}_{seurat_obj$cell_id}")
seurat_obj$cell_id %>% head()

df_sterile_cells$cell_id %>% head()

## Broad annotations  
broad_annotaions$cell_id <- rownames(broad_annotaions) %>% str_remove("_\\d+") 
broad_annotaions$cell_id <- glue("{broad_annotaions$sample}_{broad_annotaions$cell_id}") %>% str_remove("_(Green|Yellow|Blue|Red)")
head(broad_annotaions$cell_id)

## GC annotations 
GC_annotaions$cell_id <- rownames(GC_annotaions) %>% str_remove("_\\d+") 
GC_annotaions$cell_id <- glue("{GC_annotaions$Sample_name}_{GC_annotaions$cell_id}") %>% str_remove("_(Green|Yellow|Blue|Red)")
head(GC_annotaions$cell_id)

GC_annotaions$cell_id %>% str_remove("_[ATCG]+-1") %>% unique()

# Check overlapping cells
table(df_sterile_cells$cell_id %in% seurat_obj$cell_id)
table(broad_annotaions$cell_id %in% seurat_obj$cell_id)
table(GC_annotaions$cell_id %in% seurat_obj$cell_id)

# Clean 
## Seurat 
# df_sterile_cells$cell_id %>% length()
# df_sterile_cells$cell_id %>% unique() %>% length()
df_sterile_cells_clean <- df_sterile_cells %>%
  select(cell_id, sterile_feature, type_switch, n_molecules) %>%
  pivot_wider(
    names_from = c(sterile_feature, type_switch),
    names_sep = "_",
    values_from = n_molecules
  ) %>%
  mutate(
    total_sequential = rowSums(across(ends_with("_sequential")), na.rm = TRUE),
    total_same       = rowSums(across(ends_with("_same")), na.rm = TRUE),
    total_reverse    = rowSums(across(ends_with("_reverse")), na.rm = TRUE)
  ) 
  
df_sterile_cells_clean$cell_id %>% length()
df_sterile_cells_clean$cell_id %>% unique() %>% length()

## Broad annotations
broad_annotaions_clean <- broad_annotaions %>% 
  select(
    cell_id, L1_annotation
  )

## GC annotations
GC_annotaions_clean <- GC_annotaions %>% 
  mutate(
    L_annotation_Helene = case_when(
            !is.na(L3_GCB_annotation) ~ L3_GCB_annotation,
            .default = L1_annotation
          )
  ) %>% 
  select(
    cell_id, L_annotation_Helene, L3_GCB_annotation, L4_GCB_annotation, L6_GCB_annotation
  )

# Add to seurat metadata  
meta_add <- seurat_obj[[]] %>%
  rownames_to_column("cell") %>%
  select(cell, cell_id) %>%
  left_join(broad_annotaions_clean, by = "cell_id") %>%
  left_join(GC_annotaions_clean, by = "cell_id") %>%
  left_join(df_sterile_cells_clean, by = "cell_id") %>%
  column_to_rownames("cell") %>%
  select(-cell_id)

stopifnot(nrow(meta_add) == ncol(seurat_obj))

seurat_obj_2 <- AddMetaData(seurat_obj, meta_add)

# ------------------------------------------------------------------------------
# Check 
# ------------------------------------------------------------------------------

DimPlot(seurat_obj_2, group.by = "L1_annotation", cols = L1_colors)

FeaturePlot(seurat_obj_2, features = "IGHA1-S_sequential", raster=FALSE)
FeaturePlot(seurat_obj_2, features = "IGHA2-S_sequential", raster=FALSE)

switch_types <- c("sequential", "same", "reverse")

switch_description <- c(
  sequential = "the transcript isotype lies downstream of the cell's BCR isotype",
  same       = "the transcript isotype is the cell's BCR isotype",
  reverse    = "the transcript isotype lies upstream of the cell's BCR isotype"
)

annotation_levels <- list(
  L1 = list(col = "L1_annotation", colors = L1_colors, width = 10),
  L3 = list(col = "L3_GCB_annotation", colors = L3_GC_pretty_colors, width = 10)
)

for (switch_type in switch_types) {
  
  features <- grep(glue("-S_{switch_type}$"), colnames(seurat_obj_2@meta.data), value = TRUE)
  
  for (level in names(annotation_levels)) {
    
    annot_col <- annotation_levels[[level]]$col
    
    df_frac <- seurat_obj_2[[]] %>%
      select(all_of(c(annot_col, features))) %>%
      filter(!is.na(.data[[annot_col]])) %>%
      pivot_longer(-all_of(annot_col), names_to = "feature", values_to = "n_molecules") %>%
      filter(!is.na(n_molecules)) %>%                      # cells where this switch type applies
      mutate(feature = str_remove(feature, glue("_{switch_type}$"))) %>%
      group_by(.data[[annot_col]], feature) %>%
      summarise(n_cells = n(), frac = mean(n_molecules > 0), .groups = "drop")
    
    p <- ggplot(df_frac, aes(x = .data[[annot_col]], y = frac, fill = .data[[annot_col]])) +
      geom_col() +
      geom_text(aes(label = n_cells), vjust = -0.4, size = 2.5) +
      scale_fill_manual(values = annotation_levels[[level]]$colors) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
      facet_wrap(~ feature, ncol = 4) +
      labs(
        x = NULL,
        y = "Cells with sterile transcript",
        title = glue("Sterile transcripts: {switch_type}"),
        caption = glue(
          "Bars: percentage of cells with at least one sterile molecule of the given isotype.\n",
          "Numbers: cells per cell type for which {switch_description[[switch_type]]} (denominator)."
        )
      ) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    
    ggsave(glue("{outdir}/{switch_type}_transcripts_{level}.png"), p,
           width = annotation_levels[[level]]$width, height = 5.5)
  }
}

