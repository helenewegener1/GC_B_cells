library(glue)
library(tidyverse)

source("10_broad_annotation/script/color_palette.R")

# ==============================================================================
# NOTES / ASSUMPTIONS - please check before running
# ==============================================================================
# 1. Starts from 20_VDJ/out/combined.TCR.filtered.rds, produced by
#    20_VDJ/script/TCR_make_combined.filtered.R (scRepertoire::combineTCR).
#    That object is a named list, one element per sample/follicle, one row
#    per cell, with CTstrict = the "strict" paired TRA+TRB clonotype ID.
# 2. There is no TCR equivalent of somatic hypermutation / germline
#    reconstruction, so this script does NOT reproduce 45_immcantation /
#    45_changeo / 48_GCtree. Clonal identity comes straight from CTstrict.
# 3. sample_clean / sample_clean_fol / patient_id are re-derived here with the
#    SAME suffix-stripping rules used for BCR in
#    20_VDJ/script/BCR_01_make_combined.filtered.R, so that TCR compartment
#    names line up with the existing sample_clean_plot_colors palette in
#    10_broad_annotation/script/color_palette.R. Re-check the regexes below
#    still cover all of your sample names (new samples/pools would need new
#    patterns added).
# 4. L1_annotation (used to select Tfh_cells) is pulled in by joining to a
#    Seurat object on a constructed `cell_id`. The join key mirrors the
#    pattern used in 60_PC_clones/05_PC_general_figures.R's APackOfTheClones
#    section (`glue("{sample}_{rownames(.)}")`), but that pattern was written
#    for the BCR object -- confirm it matches for the T cell object/barcodes
#    before trusting downstream filtering. If cell_id doesn't line up, the
#    L1_annotation join below will silently produce NAs (see the sanity
#    `table()` check).
# ==============================================================================

seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

combined.TCR.filtered <- readRDS("20_VDJ/out/combined.TCR.filtered.rds")

names(combined.TCR.filtered) %>% head(20)

# ------------------------------------------------------------------------------
# Flatten list -> one row per cell
# ------------------------------------------------------------------------------

# Bind row together and add column with sample_fol_name (name of list)
df_tcr <- imap(combined.TCR.filtered, function(df, sample_fol_name) {
  df %>% mutate(sample_fol_name = sample_fol_name)
}) %>%
  bind_rows()

# Drop HTO Doublet / Negative calls (mirrors BCR pipeline convention)
df_tcr <- df_tcr %>%
  filter(!str_detect(sample_fol_name, "Doublet|Negative"))

nrow(df_tcr)
table(df_tcr$sample_high_level)

# ------------------------------------------------------------------------------
# Clean sample names -> sample_clean / sample_clean_fol / patient_id
# (Same suffix-stripping rules as BCR_01_make_combined.filtered.R, so labels
#  match the existing sample_clean_plot_colors / L1_colors palette)
# ------------------------------------------------------------------------------

df_tcr <- df_tcr %>%
  mutate(
    sample_clean_fol = sample_fol_name %>%
      str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC") %>%
      str_remove_all("-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-CD19-Pool1|-CD19-Pool2"),
    sample_clean = sample_clean_fol %>% str_split_i("_", 1),
    patient_id   = sample_clean %>% str_split_i("-", 1)
  )

df_tcr$sample_clean %>% unique()
df_tcr$sample_clean_fol %>% unique() %>% head(20)
df_tcr$patient_id %>% table()

# ------------------------------------------------------------------------------
# Join Seurat metadata for L1_annotation (Tfh_cells)
# ------------------------------------------------------------------------------

meta <- seurat_integrated@meta.data %>%
  rownames_to_column("row_barcode") %>%
  mutate(cell_id = glue("{sample}_{row_barcode}") %>% str_remove("_\\d+")) %>%  
  select(cell_id, L1_annotation)

# meta$cell_id %>% tail()
meta %>%  
  filter(str_detect(cell_id, "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH") & L1_annotation == "Tfh_cells")

df_tcr <- df_tcr %>%
  mutate(cell_id = glue("{sample_high_level}_{barcode}") %>% str_remove("_\\d$"))  

df_tcr %>% 
  filter(str_detect(cell_id, "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH"))

# df_tcr$cell_id %>% tail()

# Sanity check before trusting the join
# 180 cells have TCR informtaion but is not found in seurat object --> these cells were probably filtered out during GEX QC.
table(df_tcr$cell_id %in% meta$cell_id)

# df_tcr$cell_id[df_tcr$cell_id %in% meta$cell_id] # Checking which cells these are

df_tcr_filtered <- df_tcr %>%
  inner_join(meta, by = "cell_id") %>%  # Inner join to not keep the bad quality cells 
  filter(L1_annotation == "Tfh_cells") # Only keep Tfh cells 
  
# ------------------------------------------------------------------------------
# Filter to the T cell population of interest and export
# ------------------------------------------------------------------------------

nrow(df_tcr_filtered)

saveRDS(df_tcr_filtered, "70_TCR_analysis/out/df_tcr_filtered.rds")

