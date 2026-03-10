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
combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")

# ------------------------------------------------------------------------------
# Check barcodes 
# ------------------------------------------------------------------------------

(seurat_obj) %>% head() # "HH117-SILP-INF_AAACCAAAGCCGACAT-1_1" "HH117-SILP-INF_AAACCAAAGGGTAGGC-1_1"
combined.BCR.filtered$HH117$barcode %>% head() # "HH117-SILP-INF_ACGACAGAGAGTGCGA-1" "HH117-SILP-INF_AGCGAGTCAGTGTTCG-1"

# Clean seurat barcodes
seurat_obj$sample_high_level <- seurat_obj$sample %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
seurat_obj[[]] <- seurat_obj[[]] %>% mutate(sample = ifelse(!is.na(manual_ADT_ID), glue("{sample_high_level}_{manual_ADT_ID}"), sample_high_level))
barcode_seurat <- rownames(seurat_obj[[]]) %>% str_extract("(?<=_)[^_]+_[^_]+")
seurat_obj$barcode <- paste(seurat_obj$sample, barcode_seurat, sep = "_")
seurat_obj$barcode <- seurat_obj$barcode %>% str_remove_all("-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-CD19-Pool1|-CD19-Pool2")

# Check - all should be true
table(combined.BCR.filtered$HH117$barcode %in% seurat_obj$barcode)
table(combined.BCR.filtered$HH119$barcode %in% seurat_obj$barcode)

# Assign new barcodes as colnames
seurat_obj$barcode %>% unique() %>% length() == seurat_obj$barcode %>% length()
colnames(seurat_obj) <- seurat_obj$barcode

# ------------------------------------------------------------------------------
# Combine 
# ------------------------------------------------------------------------------

# Combine combined.BCR.filtered_sample and seurat_obj
seurat_obj_BCR <- combineExpression(
  combined.BCR.filtered %>% bind_rows(),
  seurat_obj,
  cloneCall = "manual_cluster",
  proportion = TRUE
)

saveRDS(seurat_obj_BCR, "40_VDJ_integrated/out/seurat_obj_BCR.rds")
  
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

Idents(seurat_obj_BCR) <- "celltype_broad"

# Define initial plot
p <- vizAPOTC(
  seurat_obj_BCR, 
  clonecall = "manual_cluster",
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
df_colors <- seurat_obj_BCR[[]] %>% filter(!is.na(manual_cluster)) %>% group_by(celltype_broad) %>% distinct(manual_cluster) %>% summarise(count = n())

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
      filter(patient == x & !is.na(manual_cluster) & !str_detect(manual_cluster, "NA") & celltype_broad == "GC_B_cells") %>% 
      summarize(CT_abundance = n(), .by = "manual_cluster") %>% 
      arrange(desc(CT_abundance)) %>%
      head(n_clones, n = n_clones) %>% 
      pull(manual_cluster)
  }
)

names(top_clones) <- unique_patients

# Size of clone types (how many cells)
seurat_obj_BCR[[]] %>% filter(manual_cluster == "cluster.4792_cluster.461") %>% nrow()
seurat_obj_BCR[[]] %>% filter(manual_cluster == "cluster.4834_cluster.461") %>% nrow()
seurat_obj_BCR[[]] %>% filter(manual_cluster == "cluster.99_cluster.392") %>% nrow()

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
  clonecall = "manual_cluster", 
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
p_clone <- p_clone + labs(title = "BCR PackOfClones", subtitle = glue("Top {n_clones} GCB clones per patient")) 

# Save plot
ggsave("40_VDJ_integrated/plot/BCR_PackOfClones/BCR_PackOfClones_CloneHighlight.png", p_clone, width = 14, height = 10)

# Inspect genes of clones
seurat_obj_BCR[[]] %>% filter(manual_cluster == "cluster.4792_cluster.461") %>%
  select(CTgene) %>%
  distinct() %>% 
  head()

seurat_obj_BCR[[]] %>% filter(manual_cluster == "cluster.4834_cluster.461") %>%
  select(CTgene) %>%
  distinct() %>% 
  head()

seurat_obj_BCR[[]] %>% filter(manual_cluster == "cluster.381_cluster.97") %>%
  select(CTgene) %>%
  distinct() %>% 
  head()
