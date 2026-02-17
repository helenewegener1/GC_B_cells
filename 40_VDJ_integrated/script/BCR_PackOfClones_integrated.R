library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(readxl)
library(tidyr)
library(tibble)
library(purrr)
library(scRepertoire)
library(APackOfTheClones)

seurat_obj <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")
combined.BCR.filtered_all <- readRDS("20_VDJ/out/combined.BCR.filtered.clean_all.rds")

# ------------------------------------------------------------------------------
# Check barcodes 
# ------------------------------------------------------------------------------

colnames(seurat_obj) %>% head() # "HH117-SILP-INF_AAACCAAAGCCGACAT-1_1" "HH117-SILP-INF_AAACCAAAGGGTAGGC-1_1"
combined.BCR.filtered_all$barcode %>% head() # "HH117-SILP-INF_ACGACAGAGAGTGCGA-1" "HH117-SILP-INF_AGCGAGTCAGTGTTCG-1"

# ------------------------------------------------------------------------------
# Fix barcodes 
# ------------------------------------------------------------------------------

# Check what number corresponds to what sample
# Create a lookup: sample name -> suffix number
sample_suffix_map <- seurat_obj[[]] %>%
  select(sample) %>%
  rownames_to_column("barcode") %>%
  mutate(suffix = str_split_i(barcode, "_", 3)) %>%
  select(sample, suffix) %>%
  distinct() %>% 
  rename(sample_high_level = sample)

sample_suffix_map

# Then fix BCR barcodes to match integrated object
combined.BCR.filtered_all <- combined.BCR.filtered_all %>%
  left_join(sample_suffix_map, by = "sample_high_level") %>%
  mutate(
    barcode_integrated = paste0(
      barcode,
      "_", suffix                            # Add integration suffix
    )
  ) %>% 
  rename(barcode_rm = barcode, 
         barcode = barcode_integrated)

# Check matches
table(colnames(seurat_obj) %in% combined.BCR.filtered_all$barcode) # TRUE = 76175

seurat_obj[[]] %>% filter(!is.na(CTstrict)) %>% rownames() %>% length() # 76175
combined.BCR.filtered_all$barcode %>% length() # 76175

# ------------------------------------------------------------------------------
# Combine 
# ------------------------------------------------------------------------------

# Combine combined.BCR.filtered_sample and seurat_obj
seurat_obj_BCR <- combineExpression(
  combined.BCR.filtered_all,
  seurat_obj,
  cloneCall = "strict",
  proportion = TRUE
)
  
# ------------------------------------------------------------------------------
# Plot 
# ------------------------------------------------------------------------------

# My colors <3
source("10_broad_annotation/script/color_palette.R")

df_celltype_colors <- celltype_colors %>% as.data.frame() %>% rownames_to_column("celltype_broad")
colnames(df_celltype_colors) <- c("celltype_broad", "new_color")

# -----------------------------------
# UMAP
# -----------------------------------

Idents(seurat_obj_BCR) <- "celltype_broad"

UMAPPlot(seurat_obj_BCR, label = TRUE) + 
  labs(title = "UMAP - Integrated") + 
  NoLegend() + 
  scale_color_manual(values = celltype_colors)

ggsave("40_VDJ_integrated/plot/BCR_PackOfClones/UMAP.png", width = 12, height = 10)

# -----------------------------------
# vizAPOTC
# -----------------------------------
# Define initial plot
p <- vizAPOTC(
  seurat_obj_BCR, 
  clonecall = "strict", 
  verbose = FALSE, 
  reduction_base = "RNA_umap.harmony", # integrated
  legend_text_size = 3, 
  # legend_spacing = 5, 
  add_size_legend = FALSE,
  legend_position = "top_right",
  show_labels = TRUE, 
  label_size = 3
)

# Extract celltype_broads without NA in CTstrict
df_colors <- seurat_obj_BCR[[]] %>% filter(!is.na(CTstrict)) %>% group_by(celltype_broad) %>% distinct(CTstrict) %>% summarise(count = n())

# Join with colors I want 
df_colors <- df_colors %>% left_join(df_celltype_colors, b = "celltype_broad")

# Extract colors from original plot
p_colors <- data.frame(
  old_color = rle(p$data$color)$values,
  count = rle(p$data$color)$lengths
)

# Make dataframe for updated colors 
df_colors <- p_colors %>% left_join(df_colors, by = "count")

# Make vector with new colors - one item per clone...
color_vector <- rep(df_colors$new_color, df_colors$count)

# Change colors of initial plot with the color_vector
p$data$color <- color_vector

p <- p + labs(title = "BCR PackOfClones - Integrated")

ggsave("40_VDJ_integrated/plot/BCR_PackOfClones/BCR_PackOfClones.png", p, width = 14, height = 10)
  
# ------------------------------------------------------------------------------
# Plot - Highlight clone
# ------------------------------------------------------------------------------

# Define unique patients
unique_patients <- seurat_obj_BCR[[]]$patient %>% unique()

# N top clones 
n_clones <- 5

# Define most expanded clones in each patient
top_clones <- lapply(
  unique_patients, 
  function(x) {
    seurat_obj_BCR[[]] %>% 
      filter(patient == x & !is.na(CTstrict)) %>% 
      group_by(CTstrict) %>% 
      summarize(CTstrict_abundance = n()) %>% 
      arrange(desc(CTstrict_abundance)) %>%
      head(n_clones) %>% 
      pull(CTstrict)
  }
)

names(top_clones) <- unique_patients

# Combine into one list for showCloneHighlight 
clones_to_highlight <- top_clones %>% unlist()
names(clones_to_highlight) <- NULL 

# Define colors 
blue_shades <- RColorBrewer::brewer.pal(n_clones, "Blues") %>% rev() # dark to light - HH117
green_shades <- RColorBrewer::brewer.pal(n_clones, "Greens") %>% rev() # dark to light - HH119 

clone_shades <- c(blue_shades, green_shades)

# Initiate plot
p_clone <- vizAPOTC(
  seurat_obj_BCR, 
  clonecall = "strict", 
  verbose = FALSE, 
  reduction_base = "RNA_umap.harmony", # integrated
  legend_text_size = 3, 
  # legend_spacing = 5, 
  add_size_legend = FALSE,
  show_labels = TRUE, 
  label_size = 3
) %>%
  showCloneHighlight(clones_to_highlight,
                     color_each = clone_shades, 
                     fill_legend = TRUE) 

# Finish plot
p_clone <- p_clone + labs(title = "BCR PackOfClones", subtitle = glue("Top {n_clones} clones per patient")) 

# Save plot
ggsave("40_VDJ_integrated/plot/BCR_PackOfClones/BCR_PackOfClones_CloneHighlight.png", p_clone, width = 14, height = 10)
