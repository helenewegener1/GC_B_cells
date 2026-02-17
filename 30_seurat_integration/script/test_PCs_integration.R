# Load libraries 
library(harmony)
library(SeuratObject)
library(Seurat)
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

integration_var <- "patient"

# Load data
seurat_obj_list <- readRDS("13_prep_integration/out/seurat_obj_prepped_list.rds")

############################# RNA integration prep #############################

# Merge objects in list
seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])

# Layers are joined based on sample
Layers(seurat_merged)

# Join layers so we can split them on the wanted variable
seurat_merged <- JoinLayers(seurat_merged)

# Split object on patient variable
seurat_merged[["RNA"]] <- split(seurat_merged[["RNA"]], f = seurat_merged[[integration_var]] %>% pull())

Layers(seurat_merged[["RNA"]]) # Check that object is split
DefaultAssay(seurat_merged) <- "RNA"

############################### Seurat workflow ################################

res <- 0.2

seurat_merged <- NormalizeData(seurat_merged, verbose = FALSE)
seurat_merged <- FindVariableFeatures(seurat_merged, verbose = FALSE)
seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
seurat_merged <- RunPCA(seurat_merged, verbose = FALSE)


# ElbowPlot(seurat_merged)
ElbowPlot(seurat_merged, ndims = 50) + theme_classic()
ggsave("14_seurat_integration/PCs_explanied/ElbowPlot.png", width = 10, height = 7)

pdf("14_seurat_integration/PCs_explanied/DimHeatmap.pdf", width = 12, height = 14)
DimHeatmap(seurat_merged, dims = 1:10, cells = 500)
DimHeatmap(seurat_merged, dims = 11:20, cells = 500)
DimHeatmap(seurat_merged, dims = 21:30, cells = 500)
dev.off()

n_dim <- 30

seurat_merged <- FindNeighbors(seurat_merged, dims = 1:n_dim, verbose = FALSE)
seurat_merged <- FindClusters(seurat_merged, resolution = res, verbose = FALSE, cluster.name = "unintegrated_clusters")
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:n_dim, verbose = FALSE, reduction.name = "umap.unintegrated")

DefaultAssay(seurat_merged)

Reductions(seurat_merged)

################################ PCs explained #################################

library(writexl)

pct <- seurat_merged[["pca"]]@stdev / sum(seurat_merged[["pca"]]@stdev) * 100
pca_df <- data.frame(
  PC = 1:length(pct),
  Variance_Percent = pct,
  Cumulative_Percent = cumsum(pct)
)

write_xlsx(pca_df, "14_seurat_integration/PCs_explanied/PCA_variance.xlsx")


# Alternative: All PCs in one sheet (wide format)
# This creates columns: PC1_gene, PC2_gene, etc.
all_pcs <- data.frame(Rank = 1:100)

n_pcs <- ncol(seurat_merged[["pca"]]@feature.loadings)

for(i in 1:n_pcs) {
  top_genes <- TopFeatures(seurat_merged[["pca"]], dim = i, nfeatures = 100,
                           projected = FALSE, 
                           balanced = FALSE) # Returns only top positive genes 
  all_pcs[[paste0("PC", i)]] <- top_genes
}

write_xlsx(all_pcs, "14_seurat_integration/PCs_explanied/PC_top100_all_in_one.xlsx")


################################################################################

# Export merged dataset
# saveRDS(seurat_merged, "14_seurat_integration/out/seurat_obj_merged_list.rds")
# seurat_merged <- readRDS("14_seurat_integration/out/seurat_obj_merged_list.rds")

# Pre integration UMAPs
## Cell type and sample
# DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = "celltype_broad", split.by = "sample", label = TRUE, ncol = 3, label.size = 2) +
#   NoLegend() + 
#   labs(title = "UMAP - RNA - pre integration") +
#   theme(
#     plot.title = element_text(size = 10),     # title font
#     strip.text = element_text(size = 6)      # facet labels (the sample names)
#   )
# 
# ggsave("14_seurat_integration/plot/UMAP_PRE_integration_split_sample.png",width = 10, height = 10)

# ## Cell type and integration var (patient)
# DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = "celltype_broad", split.by = integration_var, label = TRUE, ncol = 3, label.size = 2) +
#   NoLegend() + 
#   labs(title = "UMAP - RNA - pre integration") 
# 
# ggsave(glue("14_seurat_integration/plot/UMAP_PRE_integration_split_{integration_var}.png"), width = 12, height = 7)

## More UMAPs... 
# for (var in c(integration_var, "sample", "celltype_broad", "tissue", "inflammed")){
for (var in c("patient", "celltype_broad")){
    
  # var <- "patient"
  
  DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = var) +
    labs(title = "UMAP - RNA - pre integration",
         subtitle = var,
         caption = glue("{n_dim}PCs")) + theme(legend.text = element_text(size = 6))
  
  ggsave(glue("14_seurat_integration/plot/test_PCs/UMAP_PRE_integration_{n_dim}PCs_{var}.png"),
         width = 8,
         height = 7)

}

################################# Integration ##################################

seurat_integrated <- seurat_merged

DefaultAssay(seurat_integrated) <- "RNA"

Layers(seurat_integrated[["RNA"]])

# options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
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

# options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
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

# for (red in reductions){
  
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
  
  # # Visualize with UMAP stratified by dataset - post integration 
  # vars <- c("sample", "patient", "inflammed", "tissue", "celltype_broad")
  vars <- c("patient", "celltype_broad")
  
  for (var in vars){
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = var) +
      labs(title = glue("UMAP - post {reduction}"), 
           subtitle = var,
           caption = glue("{n_dim}PCs")) + 
      theme(legend.text = element_text(size = 8))
    
    ggsave(glue("14_seurat_integration/plot/test_PCs/UMAP_{reduction}_{n_dim}PCs_{var}.png"), 
           width = 8, 
           height = 7)
  }
  
  # DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = integration_var, group.by = "celltype_broad", ncol = 3, label.size = 2, label = TRUE) +
  #   NoLegend() + 
  #   labs(title = glue("UMAP - post {reduction}")) + 
  #   theme(legend.text = element_text(size = 8))
  # 
  # ggsave(glue("14_seurat_integration/plot/UMAP_{reduction}_{integration_var}_split.png"), 
  #        width = 12, 
  #        height = 7)
  # 
  #   # Visualize with UMAP stratified by seurat clusters - post integration 
  # res_list <- seq(0.1, 0.5, by = 0.1)
  # 
  # for (res in res_list){
  #   
  #   # res <- 0.2
  #   
  #   seurat_integrated <- FindClusters(seurat_integrated, resolution = res, cluster.name = glue("{cluster.name}_res.{res}"))
  #   
  #   DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{cluster.name}_res.{res}"), label = TRUE) +
  #     labs(title = glue("UMAP - post {reduction}"),
  #          subtitle = glue("{cluster.name}_res.{res}"))
  #   
  #   ggsave(glue(glue("14_seurat_integration/plot/clusters/UMAP_{reduction}_{assay}_snn_res_{res}.png")), 
  #          width = 8, 
  #          height = 7)
  #   
  # }

# }
