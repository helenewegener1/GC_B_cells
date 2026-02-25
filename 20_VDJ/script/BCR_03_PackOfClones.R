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

# seurat_obj_list <- readRDS("13_add_metadata/out/seurat_obj_prepped_list.rds")
# combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")
# 
# # Merge combined.BCR.filtered and update barcode column 
# combined.BCR.filtered_all <- bind_rows(combined.BCR.filtered)
# combined.BCR.filtered_all$barcode_new <- paste0(
#   combined.BCR.filtered_all$sample_clean, "_", 
#   combined.BCR.filtered_all$barcode %>% str_split_i("_", -1)
# )
# 
# # Rename barcode_new to barcode for combineExpression - THIS is what should match with colnames(seurat_obj)
# combined.BCR.filtered_all <- combined.BCR.filtered_all %>%
#   select(!barcode) %>%
#   rename(barcode = barcode_new) 
# 
# combined.BCR.filtered_all$barcode %>% head()
# 
# saveRDS(combined.BCR.filtered_all, "20_VDJ/out/combined.BCR.filtered.clean_all.rds")
# 
# # ------------------------------------------------------------------------------
# # Prep data - Per sample  
# # ------------------------------------------------------------------------------
# 
# seurat_obj_BCR_list <- list()
# 
# sample_names <- names(seurat_obj_list)
# 
# for (sample_name in sample_names){
#   
#   # Get your sample
#   # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
#   # sample_name <- "HH119-SI-PP-CD19-Pool1"
#   
#   seurat_obj <- seurat_obj_list[[sample_name]]
#   sample_name_clean <- unique(seurat_obj$sample_clean)
#   
#   # Add barcode_new to Seurat
#   seurat_obj$barcode_new <- paste0(seurat_obj$sample_clean, "_", colnames(seurat_obj))
#   
#   if (seurat_obj$barcode_new %>% length() != seurat_obj$barcode_new %>% unique() %>% length()){
#     print("barcode_new not unique not in seurat_obj")
#   } else {
#     colnames(seurat_obj) <- seurat_obj$barcode_new
#   }
# 
#   # Filter BCR data
#   combined.BCR.filtered_sample <- combined.BCR.filtered_all %>% 
#     filter(sample_high_level == sample_name)
#   
#   if (combined.BCR.filtered_sample %>% select(barcode) %>% nrow() != combined.BCR.filtered_sample %>% select(barcode) %>% unique() %>% nrow()){
#     print("barcode not unique not in combined.BCR.filtered_sample")
#   }
#   
#   # Check barcodes match
#   combined.BCR.filtered_sample$barcode %>% head()
#   colnames(seurat_obj) %>% head()
#   
#   table(combined.BCR.filtered_sample$barcode %in% colnames(seurat_obj))
#   
#   # Combine combined.BCR.filtered_sample and seurat_obj
#   seurat_obj_BCR <- combineExpression(
#     combined.BCR.filtered_sample,
#     seurat_obj,
#     cloneCall = "strict",
#     proportion = TRUE
#   )
#   
#   seurat_obj_BCR_list[[sample_name]] <- seurat_obj_BCR
# 
# }
# 
# saveRDS(seurat_obj_BCR_list, "20_VDJ/out/seurat_obj_BCR_list.rds")
# 
# 
# # ------------------------------------------------------------------------------
# # Check 
# # ------------------------------------------------------------------------------
# 
# # HH117-SILP-INF-PC
# seurat_obj_BCR_list$`HH117-SILP-INF-PC`[[]] %>% filter(!is.na(CTstrict)) %>% nrow() # 5642
# combined.BCR.filtered_all %>% filter(sample_clean == "HH117-SILP-INF" & !is.na(CTstrict)) %>% nrow() # 5642
# 
# # HH119-SI-PP
# seurat_obj_BCR_list$`HH119-SI-PP-CD19-Pool1`[[]] %>% filter(!is.na(CTstrict)) %>% nrow()
# seurat_obj_BCR_list$`HH119-SI-PP-CD19-Pool2`[[]] %>% filter(!is.na(CTstrict)) %>% nrow()
# seurat_obj_BCR_list$`HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1`[[]] %>% filter(!is.na(CTstrict)) %>% nrow()
# seurat_obj_BCR_list$`HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2`[[]] %>% filter(!is.na(CTstrict)) %>% nrow()
# 
# 7322 + 9617 + 2807 + 4784 # 24530
# 
# combined.BCR.filtered_all %>% filter(sample_clean == "HH119-SI-PP" & !is.na(CTstrict)) %>% nrow() # 24530
# 
# # HH117-SI-PP-nonINF
# seurat_obj_BCR_list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH`[[]] %>% filter(!is.na(CTstrict)) %>% nrow() # 3155
# combined.BCR.filtered_all %>% filter(sample_clean == "HH117-SI-PP-nonINF" & !is.na(CTstrict)) %>% nrow() # 3155

# ------------------------------------------------------------------------------
# Plot - Per sample  
# ------------------------------------------------------------------------------

seurat_obj_BCR_list <- readRDS("20_VDJ/out/seurat_obj_BCR_list.rds")
# combined.BCR.filtered_all <- readRDS("20_VDJ/out/combined.BCR.filtered.clean_all.rds")

# My colors <3
source("10_broad_annotation/script/color_palette.R")

df_celltype_colors <- celltype_colors %>% as.data.frame() %>% rownames_to_column("celltype_broad")
colnames(df_celltype_colors) <- c("celltype_broad", "new_color")

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH"

  seurat_obj_BCR <- seurat_obj_BCR_list[[sample_name]]
  sample_name_clean <- unique(seurat_obj_BCR$sample_clean)
  
  Idents(seurat_obj_BCR) <- "celltype_broad"
  
  # -----------------------------------
  # UMAP
  # -----------------------------------
  
  UMAPPlot(seurat_obj_BCR, label = TRUE) + 
    labs(title = "UMAP", subtitle = sample_name) + 
    NoLegend() + 
    scale_color_manual(values = celltype_colors)
  ggsave(glue("20_VDJ/plot/BCR_PackOfClones/{sample_name}_UMAP.png"), width = 14, height = 10)
  
  
  # -----------------------------------
  # vizAPOTC
  # -----------------------------------
  # Define initial plot
  p <- vizAPOTC(seurat_obj_BCR, clonecall = "strict", verbose = FALSE, 
           legend_text_size = 3, legend_spacing = 0.7, show_labels = TRUE, label_size = 3)
  
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

  p <- p + labs(title = "BCR PackOfClones", subtitle = sample_name)

  ggsave(glue("20_VDJ/plot/BCR_PackOfClones/{sample_name}_BCR_PackOfClones.png"), p, width = 14, height = 10)
  
}

# ------------------------------------------------------------------------------
# Plot - Highlight clone
# ------------------------------------------------------------------------------

# Define sample
sample_name <- "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH"

# Define clone
clone <- "IGH:Cluster.4529.IGHV4-34_IGLC:Cluster.172.IGLV1-40" # can also be a list

# Define colors 
library(wesanderson)
# clone_color <- wes_palette("GrandBudapest1")[2]
clone_color <- wes_palette("Cavalcanti1")[1] # if clone is a list, this should be a list. 

# Get seurat_obj_BCR for given sample
seurat_obj_BCR <- seurat_obj_BCR_list[[sample_name]]
sample_name_clean <- unique(seurat_obj_BCR$sample_clean)

Idents(seurat_obj_BCR) <- "celltype_broad"

# Initiate plot
p_clone <- vizAPOTC(seurat_obj_BCR, clonecall = "strict", verbose = FALSE,
                    legend_text_size = 2, legend_spacing = 0.7, show_labels = TRUE, label_size = 3) %>%
  showCloneHighlight(clone) 

# Update clone color 
p_clone$data$color <- str_replace_all(p_clone$data$color, "#F8766D", clone_color)

# Finish plot
p_clone <- p_clone + labs(title = "BCR PackOfClones", subtitle = glue("{sample_name}\n{clone}"))

# Save plot
ggsave(glue("20_VDJ/plot/BCR_PackOfClones/{sample_name}_BCR_PackOfClones_CloneHighlight.png"), p_clone, width = 14, height = 10)
