# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
# seurat_obj_list <- readRDS("12_broad_annotation/out/seurat_obj_singlets_annotated_list.rds")
seurat_obj_list <- read("11_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds")

sample_names <- names(seurat_obj_list)

source("10_broad_annotation/script/color_palette.R")

# ------------------------------------------------------------------------------
# Cell frequency per sample
# ------------------------------------------------------------------------------

# Extract cell counts by celltype_broad for each sample
celltype_counts <- lapply(names(seurat_obj_list), function(sample_name) {
  n_cells <- ncol(seurat_obj_list[[sample_name]])
  seurat_obj_list[[sample_name]][[]] %>%
    count(celltype_broad) %>%
    mutate(sample_name = sample_name,
           n_cells = n_cells)
}) %>%
  bind_rows() 

# Create stacked barplot
celltype_counts %>%
  ggplot(aes(x = sample_name, y = n, fill = celltype_broad)) +
  geom_col() +
  geom_text(
    data = celltype_counts %>% distinct(sample_name, n_cells),
    aes(x = sample_name, y = n_cells, label = n_cells),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(2:5)]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(
    title = "Number of Cells per Sample by Cell Type",
    subtitle = "Filtered seurat object", 
    x = "",
    y = "Number of cells",
    fill = "Cell Type"
  )

ggsave("13_prep_integration/plot/00_N_cells_celltypes.png", width = 12, height = 7)


################################# REMOVE DCs ###################################

n_dim <- 10
res <- 0.1

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH117-SILP-INF-PC"
  
  seurat_obj <- seurat_obj_list[[sample_name]]
  note <- ""

  n_cells <- ncol(seurat_obj_list[[sample_name]])
  DimPlot(seurat_obj_list[[sample_name]], group.by = "celltype_broad", label = TRUE, cols = celltype_colors) + NoLegend() + 
    labs(subtitle = sample_name, 
         caption = glue("N cells: {n_cells}\n{note}"))
  ggsave(glue("13_prep_integration/plot/{sample_name}.png"), width = 8, height = 8)
  
}

 ################################ ADD META DATA #################################

# Prep to save seurat objects
seurat_obj_prepped_list <- rep(0, length(seurat_obj_list)) %>% as.list()
names(seurat_obj_prepped_list) <- names(seurat_obj_list)

# Add metadata
for (sample_name in sample_names){
  
  # sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" 
  # sample_name <- "HH117-SI-MILF-INF-HLADR-AND-CD19" 
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # Define metadata
  sample <- sample_name
  patient <- str_split_i(sample_name, "-", 1)
  inflammed <- str_detect(sample_name, "-INF-")
  condition <- ifelse(patient == "HH117", "Chrons", "Control") # TODO: Update when more samples come
  # Tissue
  tissue_1 <- str_split_i(sample_name, "-", 2)
  tissue_2 <- ifelse(nchar(tissue_1) < 4, paste0("-", str_split_i(sample_name, "-", 3)), "")
  tissue <- glue("{tissue_1}{tissue_2}")
  # group <- paste(patient, ifelse(inflammed, "INF", "UNINF"), tissue, sep = "_")
  
  # print(sample)
  # print(patient)
  # print(inflammed)
  # print(tissue)
  # print("------------------")
  
  # Add metadata to seurat object
  seurat_obj@meta.data$sample <- sample
  seurat_obj@meta.data$patient <- patient
  seurat_obj@meta.data$condition <- condition
  seurat_obj@meta.data$inflammed <- inflammed
  seurat_obj@meta.data$tissue <- tissue
  sample_clean_inflammed <- ifelse(condition == "Chrons", ifelse(inflammed, '-INF', '-nonINF'), "")
  seurat_obj@meta.data$sample_clean <- glue("{patient}-{tissue}{sample_clean_inflammed}")
  # print(sample_name)
  # print("--------------------")
  # seurat_obj@meta.data$group <- group
  # seurat_obj@meta.data$group <- sub("-Pool[0-9]+$", "", seurat_obj$sample)

  seurat_obj_prepped_list[[sample_name]] <- seurat_obj
  
}

saveRDS(seurat_obj_prepped_list, "13_prep_integration/out/seurat_obj_prepped_list.rds")


