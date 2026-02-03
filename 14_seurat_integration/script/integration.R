# Load libraries 

# library(reticulate)
# use_condaenv("~/miniconda3/envs/scvi-conda/", required = TRUE)
# py_config()
# sc <- import('scanpy', convert = FALSE)
# scvi <- import("scvi", convert = FALSE)

library(harmony)
library(SeuratObject)
library(Seurat)
library(SeuratData)
library(glmGamPoi)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
# remotes::install_github('satijalab/azimuth', ref = 'master')
library(Azimuth)
library(clustree)
options(future.globals.maxSize = 1e9)

integration_var <- "patient"

# # Load data
# seurat_obj_list <- readRDS("13_prep_integration/out/seurat_obj_prepped_list.rds")
# 
# ############################# RNA integration prep #############################
# 
# # Merge objects in list
# seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])
# 
# # Layers are joined based on sample
# Layers(seurat_merged)
# 
# # Join layers so we can split them on the wanted variable
# seurat_merged <- JoinLayers(seurat_merged)
# 
# # Split object on patient variable
# seurat_merged[["RNA"]] <- split(seurat_merged[["RNA"]], f = seurat_merged[[integration_var]] %>% pull())
# 
# Layers(seurat_merged[["RNA"]]) # Check that object is split
# DefaultAssay(seurat_merged) <- "RNA"
# 
# ############################### Seurat workflow ################################
# 
# res <- 0.2
# 
# seurat_merged <- NormalizeData(seurat_merged, verbose = FALSE)
# seurat_merged <- FindVariableFeatures(seurat_merged, verbose = FALSE)
# seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
# seurat_merged <- RunPCA(seurat_merged, verbose = FALSE)
# 
# 
# # ElbowPlot(seurat_merged)
# ElbowPlot(seurat_merged, ndims = 50)
# # DimHeatmap(seurat_merged, dims = 1:30, cells = 500)
# 
# n_dim <- 10
# 
# seurat_merged <- FindNeighbors(seurat_merged, dims = 1:n_dim, verbose = FALSE)
# seurat_merged <- FindClusters(seurat_merged, resolution = res, verbose = FALSE, cluster.name = "unintegrated_clusters")
# seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:n_dim, verbose = FALSE, reduction.name = "umap.unintegrated")
# 
# DefaultAssay(seurat_merged)
# 
# Reductions(seurat_merged)
# 
# # Export merged dataset
# saveRDS(seurat_merged, "14_seurat_integration/out/seurat_obj_merged_list.rds")
# seurat_merged <- readRDS("14_seurat_integration/out/seurat_obj_merged_list.rds")
# 
# # Pre integration UMAPs
# ## Cell type and sample
# DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = "celltype_broad", split.by = "sample", label = TRUE, ncol = 3, label.size = 2) +
#   NoLegend() + 
#   labs(title = "UMAP - RNA - pre integration") +
#   theme(
#     plot.title = element_text(size = 10),     # title font
#     strip.text = element_text(size = 6)      # facet labels (the sample names)
#   )
# 
# ggsave("14_seurat_integration/plot/UMAP_PRE_integration_split_sample.png",width = 10, height = 10)
# 
# ## Cell type and integration var (patient)
# DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = "celltype_broad", split.by = integration_var, label = TRUE, ncol = 3, label.size = 2) +
#   NoLegend() + 
#   labs(title = "UMAP - RNA - pre integration") 
# 
# ggsave(glue("14_seurat_integration/plot/UMAP_PRE_integration_{integration_var}_split.png"), width = 12, height = 7)
# 
# ## More UMAPs... 
# for (var in c(integration_var, "sample", "celltype_broad", "tissue", "inflammed")){
#   
#   # var <- "patient"
#   
#   DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = var) +
#     labs(title = "UMAP - RNA - pre integration",
#          subtitle = var) + theme(legend.text = element_text(size = 6))
#   
#   ggsave(glue("14_seurat_integration/plot/UMAP_PRE_integration_{var}.png"),
#          width = 8,
#          height = 7)
#   
# }


################################# Integration ##################################
n_dim <- 10
seurat_merged <- readRDS("14_seurat_integration/out/seurat_obj_merged_list.rds")

seurat_integrated <- seurat_merged

DefaultAssay(seurat_integrated) <- "RNA"

Layers(seurat_integrated[["RNA"]])

# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = CCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "RNA_integrated.cca",
#   assay = "RNA",
#   verbose = FALSE
# )

seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "RNA_harmony",
  assay = "RNA", 
  verbose = FALSE
)

# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated, 
#   method = scVIIntegration,
#   new.reduction = "RNA_scvi",
#   conda_env = "~/miniconda3/envs/scvi-env", 
#   verbose = FALSE
# )

# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = RPCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "RNA_integrated.rpca",
#   assay = "RNA",
#   verbose = FALSE
# )
# 
# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = FastMNNIntegration,
#   orig.reduction = "pca",
#   new.reduction = "RNA_integrated.mnn",
#   assay = "RNA",
#   verbose = FALSE
# )

################### Export list of integrated Seurat objects ################### 
 
Reductions(seurat_integrated)

####################### Run UMAP using Harmony embedding #######################

reductions <- list(

  # c("RNA_integrated.cca", "RNA_umap.cca", "RNA_cca_clusters"),
  c("RNA_harmony", "RNA_umap.harmony", "RNA_harmony_clusters")
  # c("RNA_integrated.mnn", "RNA_umap.mnn", "RNA_mnn_clusters"),
  # c("RNA_integrated.rpca", "RNA_umap.rpca", "RNA_rpca_clusters")
  
)

# n_dim <- 20

for (red in reductions){
  
  red <- c("RNA_harmony", "RNA_umap.harmony", "RNA_harmony_clusters")
  
  reduction <- red[[1]]
  umap_reduction.name <- red[[2]]
  cluster.name <- red[[3]]
  
  # Either RNA or SCT
  assay <- str_split_i(reduction, "_", 1)
  
  # Set default assay 
  DefaultAssay(seurat_integrated) <- assay
  
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:n_dim)
  seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.2, cluster.name = cluster.name)
  seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction, dims = 1:n_dim, reduction.name = umap_reduction.name)
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:n_dim)
  
  Idents(seurat_integrated) <- integration_var
  
  # Visualize with UMAP stratified by dataset - post integration 
  vars <- c("sample", "patient", "inflammed", "tissue", "celltype_broad")
  
  for (var in vars){
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = var) +
      labs(title = glue("UMAP - post {reduction}"), 
           subtitle = var) + 
      theme(legend.text = element_text(size = 8))
    
    ggsave(glue("14_seurat_integration/plot/UMAP_{reduction}_{var}.png"), 
           width = 8, 
           height = 7)
  }
  
  DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = integration_var, group.by = "celltype_broad", ncol = 3, label.size = 2, label = TRUE) +
    NoLegend() + 
    labs(title = glue("UMAP - post {reduction}")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("14_seurat_integration/plot/UMAP_{reduction}_{integration_var}_split.png"), 
         width = 12, 
         height = 7)
  
    # Visualize with UMAP stratified by seurat clusters - post integration 
  res_list <- seq(0.1, 0.5, by = 0.1)
  
  for (res in res_list){
    
    # res <- 0.2
    
    seurat_integrated <- FindClusters(seurat_integrated, resolution = res, cluster.name = glue("{cluster.name}_res.{res}"))
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{cluster.name}_res.{res}"), label = TRUE) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{cluster.name}_res.{res}"))
    
    ggsave(glue(glue("14_seurat_integration/plot/clusters/UMAP_{reduction}_{assay}_snn_res_{res}.png")), 
           width = 8, 
           height = 7)
    
  }

}

# clustree

# cluster.name <- "RNA_cca_clusters"
# pdf(file = glue("04_SILP_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 12)
# clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()

cluster.name <- "RNA_harmony_clusters"
pdf(file = glue("14_seurat_integration/plot/clusters/clustree/clustree_{cluster.name}.pdf"), width = 12, height = 12)
clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
dev.off()

# cluster.name <- "RNA_rpca_clusters"
# pdf(file = glue("04_SILP_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 12)
# clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()
# 
# cluster.name <- "RNA_mnn_clusters"
# pdf(file = glue("04_SILP_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 12)
# clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()

######################### Save as h5ad file for python ######################### 

Reductions(seurat_integrated)

saveRDS(seurat_integrated, "14_seurat_integration/out/seurat_integrated_v5_RNA.rds")
# seurat_integrated <- readRDS("11_seurat_integration/out/seurat_integrated_v5_RNA.rds")

############################ Investigate clusters  #############################

VlnPlot(
  seurat_integrated,
  features = "rna_CD8A", group.by = "unintegrated_clusters", raster=FALSE
) + NoLegend() + ggtitle("CD8A - Unintegrated Clusters")

vars <- c("percent.ribo", "percent.mt", "G2M.Score", "S.Score")
var <- "percent.mt"
VlnPlot(
  seurat_integrated,
  features = var, group.by = "unintegrated_clusters", raster=FALSE
) + NoLegend() + ggtitle(glue("{var} - Unintegrated Clusters"))

p2 <- VlnPlot(
  seurat_integrated, "rna_CD8A",
  group.by = "RNA_harmony_clusters"
) + NoLegend() + ggtitle("CD8A - CCA Clusters")

p1 | p2 

######################## Save as h5ad file for python #########################

library(SeuratDisk)
library(rhdf5)

# All in one file 
obj_tmp <- seurat_integrated 

DefaultAssay(obj_tmp)

Layers(obj_tmp[["RNA"]])

# IMPORTANT: Join layers before export (h5ad doesn't support split layers)
obj_tmp[["RNA"]] <- JoinLayers(obj_tmp[["RNA"]])

obj_tmp[["RNA3"]] <- as(object = obj_tmp[["RNA"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "RNA3"
obj_tmp[["RNA"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, RNA3 = 'RNA')

# Check new names
reductions <- Reductions(obj_tmp)[str_detect(Reductions(obj_tmp), "integrated")]

# 1. Extract the embeddings
embeddings <- lapply(reductions, function(x) Embeddings(obj_tmp, reduction = x))
names(embeddings) <- reductions

# 2. Add the clean matrix as a new reduction
for(x in reductions) {
  
  new_name <- str_replace(x, "\\.", "_")
  
  new_key <- str_replace(x, "_integrated\\.", "") %>% 
    str_to_upper() %>% 
    paste0("_")
  
  obj_tmp[[new_name]] <- CreateDimReducObject(
    embeddings = embeddings[[x]], 
    key = new_key,
    assay = DefaultAssay(obj_tmp),
    global = TRUE
  )
  
  # Remove the problematic reduction to keep things clean
  obj_tmp[[x]] <- NULL
  
}

# Handle PCA 
obj_tmp[["PCA"]] <- CreateDimReducObject(
  embeddings = Embeddings(obj_tmp, reduction = "pca"), 
  key = "PCA_",
  assay = DefaultAssay(obj_tmp),
  global = TRUE
)

# Remove the problematic reduction to keep things clean
obj_tmp[["pca"]] <- NULL

umap_reductions <- Reductions(obj_tmp)[str_detect(Reductions(obj_tmp), "umap")]
for (x in umap_reductions) {
  obj_tmp[[x]] <- NULL
}

# Verify the new reduction name (optional)
Reductions(obj_tmp)

filename <- glue("11_seurat_integration/out/mydata_RNA_v5.h5Seurat")
SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)














