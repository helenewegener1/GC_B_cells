library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

df_all <- readRDS("50_BACH2/out/df_all.rds")

# Define largest CRC clone
large_crc_clone <- df_all %>%
  filter(patient_id == "HH119", L1_annotation == "GC_B_cells") %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  head(1) %>% 
  pull(clone_subgroup_id_90_similarity)

# Create outdit
outdir <- glue("50_BACH2/plot/03_isotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Overview
# ------------------------------------------------------------------------------

df_all_GC <- df_all %>% filter(L1_annotation == "GC_B_cells")
# table(df_all_GC$patient, df_all_GC$L3_GCB_annotation)

table(df_all_GC$c_call_grouped, df_all_GC$L3_GCB_annotation)

# Plot table
df <- as.data.frame(
  table(
    c_call = df_all_GC$c_call_grouped,
    L3_annotation = df_all_GC$L3_GCB_annotation,
    Patient = df_all_GC$patient
  )
)

ggplot(df, aes(x = c_call, y = L3_annotation, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Isotype",
    y = "GC B cell subset",
    fill = "Cell count"
  ) +
  facet_wrap(vars(Patient)) +
  theme_bw(base_size = 14)

ggsave(glue("{outdir}/GC_subsets_isotype_heatmap.png"), width = 13, height = 7.5)
