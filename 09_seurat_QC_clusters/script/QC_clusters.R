getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
library(purrr) # map funciton 
library(tidyr)

# Load data
seurat_obj_QC_filtered_singlets_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_singlets_list.rds")

sample_names <- names(seurat_obj_QC_filtered_singlets_list)

source("09_seurat_QC_clusters/script/functions.R")

################################################################################ 
################################## Annotation ##################################
################################################################################ 

# Check that doublets are removed 
# seurat_obj_QC_filtered_singlets_list$`HH117-SILP-INF-PC`$scDblFinder.class %>% table()
# seurat_obj_QC_filtered_singlets_list$`HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1`$scDblFinder.class %>% table()

# Load Gina annotation file 
broad_annot_file <- read_excel("00_data/Gene_markers_GL_HW_new.xlsx", sheet = "Very broad level")
detailed_annot_file <- read_excel("00_data/Gene_markers_GL_HW_new.xlsx", sheet = "More detailed level")

# Extract cell marker genes as lists
broad_markers <- broad_annot_file %>% as.list()
broad_markers <- map(broad_markers, ~ .x[!is.na(.x)])
names(broad_markers) <- c("T_cell", "B_cell", "DC", "plasmablast_plasma_cell")

detailed_markers <- detailed_annot_file %>% as.list() %>% na.omit()
detailed_markers <- map(detailed_markers, ~ .x[!is.na(.x)])
names(detailed_markers) <- c("TFH_cell", "Naive_B_cell", "Memory_B_cell", "GC_B_cell", 
                             "Activation_markers", "Immunoglobulin_subsets", "Other_markers")

#### Make sure the markers are in the same format as in the seurat object ###### 

# grep("IGHA", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

# Load sample to get genes that are in the data.  
seurat_obj <- seurat_obj_QC_filtered_singlets_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

# Update marker format for the broad markers
broad_markers <- update_marker_names(broad_markers, seurat_obj)

# grep("CD21", rownames(seurat_obj), value = TRUE)

# Update marker format for the detailed markers
detailed_markers <- update_marker_names(detailed_markers, seurat_obj)

############################# Get cell cycle score #############################

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes # MKI67 in here 

################################################################################

################################### Singlets ###################################

# Prep to save clustered seurat objects
seurat_obj_clustered_singlets_list <- rep(0, length(seurat_obj_QC_filtered_singlets_list)) %>% as.list()
names(seurat_obj_clustered_singlets_list) <- names(seurat_obj_QC_filtered_singlets_list)

# Define number of dimensions for clustering and resolution of clusters. 
n_dim <- 10
res <- 0.1

# for (method in c("standard", "SCT")){
# method <- "standard"
for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-COLP-PC"
  # sample_name <- "HH117-SILP-INF-PC"
  print(glue("--- Processing: {sample_name} ---"))
  
  # Create directory for plots of specific sample
  out_dir <- glue("09_seurat_QC_clusters/plot/singlets/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define specific seurat object 
  seurat_obj <- seurat_obj_QC_filtered_singlets_list[[sample_name]]
  
  # Workflow
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  DefaultAssay(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:n_dim, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:n_dim, verbose = FALSE)
  
  # DimPlot(seurat_obj, label = TRUE, group.by = "seurat_clusters") + NoLegend() 
  
  # generate plots 
  all_plots(seurat_obj = seurat_obj, sample_name = sample_name, n_dim = n_dim, 
            extra_title = "Singlets", out_dir = out_dir)
  
  # Save clustered seurat object 
  seurat_obj_clustered_singlets_list[[sample_name]] <- seurat_obj
  
}

# Save object 
saveRDS(seurat_obj_clustered_singlets_list, "09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")

rm(seurat_obj_QC_filtered_singlets_list)

################################################################################ 
################################### Doublets ###################################
################################################################################ 

# # Load data
# seurat_obj_QC_filtered_doublets_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_doublets_list.rds")
# 
# for (sample_name in sample_names){
#   
#   # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
#   n_doublets <- seurat_obj_QC_filtered_doublets_list[[sample_name]] %>% ncol()
#   
#   print(sample_name)
#   print(glue("N doublets: {n_doublets}"))
#   print("---------------------------------------------")
#   
# }
# 
# rm(seurat_obj_QC_filtered_doublets_list)

################################################################################ 
##################################### All ######################################
################################################################################ 

# Load data
seurat_obj_QC_filtered_all_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list.rds")

# Prep to save clustered seurat objects
seurat_obj_clustered_all_list <- rep(0, length(seurat_obj_QC_filtered_all_list)) %>% as.list()
names(seurat_obj_clustered_all_list) <- names(seurat_obj_QC_filtered_all_list)

n_dim <- 10
res <- 0.1

for (sample_name in sample_names){

  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH117-SILP-INF-PC"
  
  print(glue("--- Processing: {sample_name} ---"))
  
  # Create directory for plots of specific sample
  out_dir <- glue("09_seurat_QC_clusters/plot/all/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  seurat_obj <- seurat_obj_QC_filtered_all_list[[sample_name]]

  # Calculate and print doublet percentage
  n_doublets <- sum(seurat_obj$scDblFinder.class == "doublet")
  n_cells <- seurat_obj %>% ncol()
  percentage_doublet <- round(n_doublets/n_cells * 100, 1)

  # print(sample_name)
  # print(glue("Doublet percentage: {percentage_doublet}"))
  # print("---------------------------------------------")

  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

  DefaultAssay(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:n_dim, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:n_dim, verbose = FALSE)
  
  # generate plots 
  all_plots(seurat_obj = seurat_obj, sample_name = sample_name, n_dim = n_dim, 
            extra_title = "All", out_dir = out_dir)
  
  n_doublet <- sum(seurat_obj$scDblFinder.class == "doublet")
  
  DimPlot(seurat_obj, group.by = "scDblFinder.class") + 
    labs(
      title = "All: scDblFinder.class", 
      subtitle = sample_name,
      caption = glue("N cells: {n_cells}\nN doublet: {n_doublet}")
    )
  
  ggsave(glue("{out_dir}/{sample_name}_scDblFinder.class.png"), width = 8, height = 7)
  
  # Doublet percentage across cell cycle phases. 
  doublet_N_genes(seurat_obj = seurat_obj, sample_name = sample_name, out_dir = out_dir)
  
  # Save object
  seurat_obj_clustered_all_list[[sample_name]] <- seurat_obj

}

saveRDS(seurat_obj_clustered_all_list, "09_seurat_QC_clusters/out/seurat_obj_clustered_all_list.rds")
# rm(seurat_obj_QC_filtered_all_list)

################################################################################ 
################################## PC loadings #################################
################################################################################ 

# Sheet names in excel can only be 31 chars
sheet_names <- list(
  "HH117-SILP-INF-PC"                                = "HH117_SILP_INF_PC",
  "HH117-SILP-nonINF-PC"                             = "HH117_SILP_nonINF_PC",
  "HH117-SI-MILF-INF-HLADR-AND-CD19"                 = "HH117_MILF_INF_HLA_CD19",
  "HH117-SI-MILF-nonINF-HLADR-AND-CD19"              = "HH117_MILF_nonINF_HLA_CD19",
  "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = "HH117_PP_nonINF_HLA_CD19_GC_TFH",
  "HH119-COLP-PC"                                   = "HH119_COLP_PC",
  "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH"        = "HH119_SMILF_CD19_GC_PB_TFH",
  "HH119-SILP-PC"                                   = "HH119_SILP_PC",
  "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH"         = "HH119_MILF_CD19_GC_PB_TFH",
  "HH119-SI-PP-CD19-Pool1"                           = "HH119_PP_CD19_P1",
  "HH119-SI-PP-CD19-Pool2"                           = "HH119_PP_CD19_P2",
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"              = "HH119_PP_GC_PB_TFH_P1",
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"              = "HH119_PP_GC_PB_TFH_P2"
)

seurat_obj_singlets_clustered_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")

save_pc_loading(seurat_obj_clustered_all_list, extra_title = "all")
save_pc_loading(seurat_obj_singlets_clustered_list, extra_title = "singlets")

rm(seurat_obj_singlets_clustered_list)

# # PC DimHeatmap
# # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
# sample_name <- "HH117-SILP-INF-PC"
# 
# seurat_object <- seurat_obj_QC_filtered_singlets_list[[sample_name]]
# 
# DimHeatmap(seurat_object, dims = 1:10, cells = 500, balanced = TRUE)
# 
# seurat_object <- JackStraw(seurat_object, num.replicate = 100)
# seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
# JackStrawPlot(seurat_object, dims = 1:20)

################################################################################ 
################################ DC Annotation #################################
################################################################################ 

# Venla annotate after human atlas projection of v.8 object?
# But for now we do this. 

# Define DC cluster numbers
# Look at 09_seurat_QC_clusters/plot/clusters/

# All cells
# seurat_obj_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_all_list.rds")
# out_dir <- glue("09_seurat_QC_clusters/plot/all_DCs_removed/{sample_name}")
# dc_clusters_all <- list(
#   "HH117-SILP-INF-PC"                                = NULL,
#   "HH117-SILP-nonINF-PC"                             = NULL,
#   "HH117-SI-MILF-INF-HLADR-AND-CD19"                 = c("1"),
#   "HH117-SI-MILF-nonINF-HLADR-AND-CD19"              = c("1", "3", "5", "6 "),
#   "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = c("5", "8", "9"),
#   "HH119-COLP-PC"                                    = NULL,
#   "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH"        = NULL,
#   "HH119-SILP-PC"                                    = NULL,
#   "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH"         = NULL,
#   "HH119-SI-PP-CD19-Pool1"                           = NULL,
#   "HH119-SI-PP-CD19-Pool2"                           = NULL,
#   "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"              = NULL,
#   "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"              = NULL
# )

# Singlets
seurat_obj_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_clustered_list_singlets.rds")

dc_clusters <- list(
  "HH117-SILP-INF-PC"                                = NULL,
  "HH117-SILP-nonINF-PC"                             = NULL,
  "HH117-SI-MILF-INF-HLADR-AND-CD19"                 = c("1"),
  "HH117-SI-MILF-nonINF-HLADR-AND-CD19"              = c("1", "3", "5", "6"),
  "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = c("6"),
  "HH119-COLP-PC"                                    = NULL,
  "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH"        = NULL,
  "HH119-SILP-PC"                                    = NULL,
  "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH"         = NULL,
  "HH119-SI-PP-CD19-Pool1"                           = NULL,
  "HH119-SI-PP-CD19-Pool2"                           = NULL,
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"              = NULL,
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"              = NULL
)


# Prep to save seurat objects without DCs
seurat_obj_nonDC_list <- rep(0, length(seurat_obj_list)) %>% as.list()
names(seurat_obj_nonDC_list) <- names(seurat_obj_list)

sample_names <- names(seurat_obj_list)

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH117-SILP-INF-PC"
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  if (is.null(dc_clusters[[sample_name]])){
    seurat_obj_nonDC_list[[sample_name]] <- seurat_obj
    next
  }
  
  # Create a new logical column
  seurat_obj$DC_bool <- ifelse(seurat_obj$seurat_clusters %in% dc_clusters[[sample_name]], TRUE, FALSE)
  table(seurat_obj$DC_bool)
  n_cells <- seurat_obj %>% ncol()
  n_DCs <- sum(seurat_obj$DC_bool)
  
  # Plot
  # DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) + NoLegend() + 
  #   labs(title = "DCs not removed",
  #        subtitle = sample_name,
  #        caption = glue("N cells: {n_cells}")) 
  # ggsave(glue("{out_dir}/{sample_name}_clusters.png"), width = 8, height = 8)
  
  out_dir <- glue("09_seurat_QC_clusters/plot/singlets_DCs_removed/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  DimPlot(seurat_obj, group.by = "DC_bool", label = TRUE, cols = c("grey", "orange")) + 
    NoLegend() + 
    labs(caption = glue("N cells: {n_cells}\nN DCs: {n_DCs}"))
  ggsave(glue("{out_dir}/{sample_name}_DC_bool.png"), width = 8, height = 8)
  
  # # Split object in DCs and non-DCs
  # # First, DCs
  # seurat_obj_DC_list <- list()
  # 
  # seurat_obj_DC <- subset(seurat_obj, subset = DC_bool == TRUE)
  # seurat_obj_DC[["ADT"]] <- NULL
  # 
  # seurat_obj_DC_list[[sample_name]] <- seurat_obj_DC
  
  # Then, other
  seurat_obj_nonDC <- subset(seurat_obj, subset = DC_bool == FALSE)
  seurat_obj_nonDC <- NormalizeData(seurat_obj_nonDC, verbose = FALSE)
  seurat_obj_nonDC <- FindVariableFeatures(seurat_obj_nonDC, verbose = FALSE)
  seurat_obj_nonDC <- ScaleData(seurat_obj_nonDC, verbose = FALSE)
  
  DefaultAssay(seurat_obj_nonDC)
  seurat_obj_nonDC <- RunPCA(seurat_obj_nonDC, verbose = FALSE)
  ElbowPlot(seurat_obj_nonDC)
  
  n_dim <- 10
  res = 0.1
  
  seurat_obj_nonDC <- FindNeighbors(seurat_obj_nonDC,  dims = 1:n_dim, verbose = FALSE)
  seurat_obj_nonDC <- FindClusters(seurat_obj_nonDC, resolution = res, verbose = FALSE)
  seurat_obj_nonDC <- RunUMAP(seurat_obj_nonDC, reduction = "pca", dims = 1:n_dim, verbose = FALSE)
  
  # Plot
  n_cells <- ncol(seurat_obj_nonDC)
  DimPlot(seurat_obj_nonDC, group.by = "seurat_clusters", label = TRUE) + NoLegend() + 
    labs(title = glue("DCs removed"),
         subtitle = sample_name, 
         caption = glue("N cells: {n_cells}\nN dim: {n_dim}\nResolution: {res}"))
  ggsave(glue("{out_dir}/{sample_name}_clusters_DCs_removed.png"), width = 8, height = 8)
  
  # Save seurat object without DCs in list
  seurat_obj_nonDC_list[[sample_name]] <- seurat_obj_nonDC

}

names(seurat_obj_nonDC_list)

# Save lists of DC and nonDC seurat objects 
# saveRDS(seurat_obj_DC_list, "09_seurat_QC_clusters/out/seurat_obj_DC_list.rds")
saveRDS(seurat_obj_nonDC_list, "09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")

# seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")

################################## Annotation ##################################

