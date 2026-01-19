setwd("~/ciir/people/helweg/projects/GC_B_cells/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(glmGamPoi)
# devtools::install_github("constantAmateur/SoupX", ref='devel')
library(SoupX)
library(multtest)
library(scDblFinder)
library(decontX)
library(readxl)
library(purrr)
library(patchwork)
library(png)
library(grid)

source("07_seurat_QC/script/merged_clusters.R")

version <- "v8"

# Load data
seurat_obj_list <- readRDS(glue("06_seurat_load/out/seurat_obj_list_{version}.rds")) # cellranger filtered

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

# Load sample to get genes that are in the data.  
seurat_obj <- seurat_obj_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

source("09_seurat_QC_clusters/script/functions.R")

# Update marker format for the broad markers
broad_markers <- update_marker_names(broad_markers, seurat_obj)

########################### Define major cell types ############################

# Prep to save clustered seurat objects
seurat_obj_clustered_list <- rep(0, length(seurat_obj_list)) %>% as.list()
names(seurat_obj_clustered_list) <- names(seurat_obj_list)

# N PCs
# n_pcs <- list(
#   "HH117-SILP-INF-PC" = ,                        
#   "HH117-SILP-nonINF-PC" = ,                           
#   "HH117-SI-MILF-INF-HLADR-AND-CD19" = ,                 
#   "HH117-SI-MILF-nonINF-HLADR-AND-CD19" = ,             
#   "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = ,
#   "HH119-COLP-PC" = ,                                   
#   "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH" = ,    
#   "HH119-SILP-PC" = ,                                   
#   "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH" = ,      
#   "HH119-SI-PP-CD19-Pool1" = ,                        
#   "HH119-SI-PP-CD19-Pool2" = ,                        
#   "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = ,                
#   "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = 
# )

sample_names <- names(seurat_obj_list)

res <- 0.1
n_dims <- 10

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-COLP-PC"
  
  print(glue("--- Processing: {sample_name} ---"))
  
  # Get seurat object 
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # Create directory for plots of specific sample
  out_dir <- glue("07_seurat_QC/plot_{version}/01_clusters/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Seurat workflow so I can UMAP
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  
  ElbowPlot(seurat_obj) + labs(title = sample_name) #to determine dimentions used for following steps in doublet detection. Adjust dims. 
  ggsave(glue("{out_dir}/{sample_name}_elbow.png"), width = 9, height = 5.5)
  
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  
  n_cells <- ncol(seurat_obj)
  
  # Plot
  DimPlot(seurat_obj, reduction = 'umap', label = TRUE) + NoLegend() + 
    labs(title = glue("Seurat clusters {version}"),
         subtitle = sample_name, 
         caption = glue("N cells: {n_cells}\nN dim: {n_dims}\nresolution: {res}"))
  ggsave(glue("{out_dir}/{sample_name}_clusters.png"), width = 8, height = 8)
  
  ####################### FeaturePlot with broad_markers ####################### 
  for (markers in names(broad_markers)){
    
    FeaturePlot(seurat_obj, features = broad_markers[[markers]], ncol = 3) + 
      plot_annotation(title = glue("{markers} - {version}"),
                      subtitle = sample_name,
                      caption = glue("N cells: {n_cells}\nN dim: {n_dims}\nresolution: {res}"))
    
    # Adjust heigh of plot to number of markers
    n_markers <- broad_markers[[markers]] %>% length()
    height <- (n_markers/3) * 4
    
    ggsave(glue("{out_dir}/{sample_name}_broad_{markers}.png"), width = 14, height = height)
    
  }
  
  ############################################################################## 
  
  ##################### FeaturePlot with detailed_markers ######################
  for (markers in names(detailed_markers)){
    
    FeaturePlot(seurat_obj, features = detailed_markers[[markers]], ncol = 3) + 
      plot_annotation(title = glue("{markers} - {version}"),
                      subtitle = sample_name,
                      caption = glue("N cells: {n_cells}\nN dim: {n_dims}\nresolution: {res}"))
    
    # Adjust heigh of plot to number of markers
    n_markers <- detailed_markers[[markers]] %>% length()
    height <- (n_markers/3) * 4
    
    ggsave(glue("{out_dir}/{sample_name}_detailed_{markers}.png"), width = 14, height = height)
    
  }
  
  ##############################################################################
  
  # Save object
  seurat_obj_clustered_list[[sample_name]] <- seurat_obj
  
}

# saveRDS(seurat_obj_clustered_list, glue("07_seurat_QC/out/seurat_obj_clustered_list_{version}.rds"))

# Get number of clusters
for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  n_clusters <- seurat_obj_clustered_list[[sample_name]]$seurat_clusters %>% levels() %>% length()
  
  print(sample_name)
  print(glue("N clusters: {n_clusters}"))
  print("---------------------------------------------")
  
}

################################################################################
############################ GINA CLUSTER MERGING ##############################
################################################################################

# Gina got the plots from 07_seurat_QC/plot/clusters.
# Gina drew which clusters should be merged based on expression.
# This was to make sure that doublet finder did not merge two too similar cells and called them a doublet.
# Gina's drawing can be seen in 07_seurat_QC/plot/gina_merged_clusters.

merged_clusters <- merged_clusters_all[[version]]

# sample_names <- names(seurat_obj_clustered_list)
for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  seurat_obj <- seurat_obj_clustered_list[[sample_name]]
  n_cells <- ncol(seurat_obj)
  
  # Create directory for plots of specific sample
  out_dir <- glue("07_seurat_QC/plot_{version}/03_clusters_for_scDblFinder/")
  dir.create(out_dir, showWarnings = FALSE)
  
  # Extract old clusters
  old_clusters <- as.character(seurat_obj$seurat_clusters)
  names(old_clusters) <- colnames(seurat_obj)
  
  # Extract mapping of new and old clusterss
  mapping_list <- merged_clusters[[sample_name]]
  
  # Create named vector: old cluster -> new cluster
  cluster_mapping <- unlist(lapply(names(mapping_list), function(new_cluster) {
    old_clusters <- mapping_list[[new_cluster]]
    setNames(rep(new_cluster, length(old_clusters)), old_clusters)
  }))
  
  # Map merged clusters
  new_clusters <- cluster_mapping[old_clusters]
  names(new_clusters) <- colnames(seurat_obj)
  
  # Add merged clusters to metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = as.data.frame(new_clusters), col.name = "merged_clusters")
  seurat_obj_clustered_list[[sample_name]] <- seurat_obj
  
  # Check
  table("old" = seurat_obj$seurat_clusters, "new" = seurat_obj$merged_clusters)
  
  # Plot
  img <- readPNG(glue("07_seurat_QC/plot_{version}/02_gina_merged_clusters/{sample_name}_clusters.png"))
  p_gina_draw <- ggplot() +
    annotation_custom(
      rasterGrob(img)
    ) +
    theme_void()
  
  p_new <- DimPlot(seurat_obj, reduction = 'umap', group.by = "merged_clusters", label = TRUE) + NoLegend() +
    labs(title = "Seurat clusters merged",
         subtitle = sample_name,
         caption = glue("N cells: {n_cells}\nN dim: {n_dims}\nresolution: {res}"))
  
  p_final <- p_gina_draw + p_new + plot_annotation(title = version)
  
  ggsave(glue("{out_dir}/{sample_name}_merged_clusters.png"), p_final, width = 12, height = 7)
  
}

# saveRDS(seurat_obj_clustered_list, glue("07_seurat_QC/out/seurat_obj_clustered_list_{version}.rds"))
seurat_obj_clustered_list <- readRDS(glue("07_seurat_QC/out/seurat_obj_clustered_list_{version}.rds"))


################################################################################
################################################################################

# # Investigate need for removal of empty droplets
# # https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html
# print("---------------------------------------------------------------")
# print("Investigate need for removal of empty droplets")
# for (sample_name in names(seurat_obj_list)){
#
#   print("----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----")
#   print(sample_name)
#
#   ######################## CHECK EMPTY DROPLETS IN RAW #########################
#
#   # Load raw (totally non-filtered) counts from cellranger
#   raw_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/multi/count/raw_feature_bc_matrix"))
#   seurat_obj_raw <- CreateSeuratObject(counts = raw_counts)
#
#   # Transform to SingleCellExperiment object
#   if (length(Layers(seurat_obj_raw))){
#     sc_exp_raw <- seurat_obj_raw %>% JoinLayers() %>% as.SingleCellExperiment()
#   } else {
#     sc_exp_raw <- as.SingleCellExperiment(seurat_obj_raw)
#   }
#
#   # Run emptyDrops
#   set.seed(100)
#   e.out <- emptyDrops(counts(sc_exp_raw))
#
#   # See ?emptyDrops for an explanation of why there are NA values.
#   is.cell <- table(e.out$FDR <= 0.01, useNA = "ifany")
#
#   # Computing barcode ranks
#   br.out <- barcodeRanks(counts(sc_exp_raw))
#
#   # Making Barcode Rank Plot.
#   pdf(glue("07_seurat_QC/plot/emptyDrops/emptyDrops_{sample_name}_raw.pdf"), width = 8, height = 6)
#
#   plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",
#        main=glue("{sample_name} raw (prefilter) Barcode Rank Plot"),
#        sub = glue("Cells: {is.cell[['TRUE']]}, Not cells: {is.cell[['FALSE']]}, NAs: {is.cell[[3]]}"))
#   o <- order(br.out$rank)
#   lines(br.out$rank[o], br.out$fitted[o], col="red")
#
#   abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
#   abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
#   legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
#          legend=c("knee", "inflection"))
#
#   dev.off()
#
#   ############################ Ambiant RNA with SoupX ############################
#
#   # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
#
#   raw_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/multi/count/raw_feature_bc_matrix"))
#   cell_counts <- Read10X(data.dir = glue("05_run_cellranger/out/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/count/sample_filtered_feature_bc_matrix"))
#
#   if (is.null(names(raw_counts))){
#     sc <- SoupChannel(raw_counts, cell_counts, calcSoupProfile=FALSE)
#   } else if (length(names(raw_counts)) > 1) {
#     sc <- SoupChannel(raw_counts$`Gene Expression`, cell_counts$`Gene Expression`, calcSoupProfile=FALSE)
#   }
#
#   sc <- estimateSoup(sc)
#
#   # Perform quick clustering steps
#   seurat_temp <- CreateSeuratObject(counts = sc$toc)
#   seurat_temp <- NormalizeData(seurat_temp, verbose = FALSE)
#   seurat_temp <- FindVariableFeatures(seurat_temp, verbose = FALSE)
#   seurat_temp <- ScaleData(seurat_temp, verbose = FALSE)
#   seurat_temp <- RunPCA(seurat_temp, npcs = 20, verbose = FALSE)
#   seurat_temp <- FindNeighbors(seurat_temp, dims = 1:20, verbose = FALSE)
#   seurat_temp <- FindClusters(seurat_temp, resolution = 0.7, verbose = FALSE)
#   seurat_temp <- RunUMAP(seurat_temp, dims = 1:20)
#
#   DimPlot(seurat_temp, group.by = "seurat_clusters") + labs(title = "Pre-SoupX")
#   ggsave(glue("07_seurat_QC/plot/SoupX/pre_SoupX_{sample_name}.png"), width = 7, height = 6)
#
#   # Extract and set the clusters
#   clusters <- seurat_temp@meta.data$seurat_clusters
#   sc <- setClusters(sc, clusters)
#   sc <- autoEstCont(sc)
#   out <- adjustCounts(sc)
#
#   # Get the UMAP coordinates and add them to the SoupChannel object
#   # umap_coords <- seurat_temp@reductions$umap@cell.embeddings
#   # sc <- setDR(sc, umap_coords)
#   # plotMarkerMap(sc, "MKI67")
#
#   # Use the 'out' matrix to create a Seurat object
#   final_seurat_obj <- CreateSeuratObject(counts = out, project = "SoupX_Corrected_scRNA")
#   # Add the per-cell rho to the final Seurat object
#   final_seurat_obj$soupX_rho_per_cell <- sc$metaData$rho
#
#   # Visualize the contamination fraction across your UMAP
#   final_seurat_obj <- NormalizeData(final_seurat_obj, verbose = FALSE)
#   final_seurat_obj <- FindVariableFeatures(final_seurat_obj, verbose = FALSE)
#   final_seurat_obj <- ScaleData(final_seurat_obj, verbose = FALSE)
#   final_seurat_obj <- RunPCA(final_seurat_obj, npcs = 20, verbose = FALSE)
#   final_seurat_obj <- FindNeighbors(final_seurat_obj, dims = 1:20, verbose = FALSE)
#   final_seurat_obj <- FindClusters(final_seurat_obj, resolution = 0.7, verbose = FALSE)
#   final_seurat_obj <- RunUMAP(final_seurat_obj, dims = 1:20)
#
#   DimPlot(final_seurat_obj, group.by = "seurat_clusters") + labs(title = "Post-SoupX")
#   ggsave(glue("07_seurat_QC/plot/SoupX/post_SoupX_{sample_name}.png"), width = 7, height = 6)
#
# }

############################# Get cell cycle score #############################

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes # MKI67 in here

################################ scDblFinder on cellranger filtered ################################

# seurat_obj_clustered_list <- readRDS(glue("07_seurat_QC/out/seurat_obj_clustered_list_{version}.rds"))
n_dims <- 10

# Initialize final QC list
seurat_obj_QC <- list()

sample_names <- names(seurat_obj_clustered_list)

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-SILP-PC"
  # print("---------------------------------------------------------------")
  # print(sample_name)
  #
  # ############################ Ambiant RNA with decontX ############################
  # print("----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----")
  # print("decontX")
  
  raw_counts <- Read10X(data.dir = glue("05_run_cellranger/out_{version}/res_{sample_name}/outs/multi/count/raw_feature_bc_matrix"))
  cell_counts <- Read10X(data.dir = glue("05_run_cellranger/out_{version}/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/count/sample_filtered_feature_bc_matrix"))
  
  if (is.null(names(raw_counts))){
    
    sce <- decontX(cell_counts, background = raw_counts)
    
  } else if (length(names(raw_counts)) > 1) {
    
    sce <- decontX(cell_counts$`Gene Expression`, background = raw_counts$`Gene Expression`)
    
  }
  
  # Get seurat object
  seurat_obj <- seurat_obj_clustered_list[[sample_name]]
  
  # # Add decontX to contamination "score" to metadata
  seurat_obj <- AddMetaData(seurat_obj, sce$contamination, "sce_contamination")
  
  # seurat_obj@meta.data$sce_contamination
  
  # print("----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----")
  # print("scDblFinder")
  
  # Create directory for plots of specific sample
  out_dir <- glue("07_seurat_QC/plot_{version}/04_scDblFinder/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get count matrix
  # counts <- seurat_obj@assays$RNA$counts
  
  # n_dims <- 10
  # # Seurat workflow so I can UMAP
  # seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  # seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  # seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  # # seurat_obj <- SCTransform(seurat_obj)
  # #Doublet detection
  # seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  # ElbowPlot(seurat_obj) #to determine dimentions used for following steps in doublet detection. Adjust dims.
  # seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  # seurat_obj <- FindClusters(seurat_obj, resolution = 0.05, verbose = FALSE)
  
  # n_clusters - no relevant if we set the number of clusters in scDblFinder
  # n_clusters <- seurat_obj$seurat_clusters %>% unique()
  # n_clusters
  # table(seurat_obj$seurat_clusters)
  
  # Create SingelSellExperiment object
  sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
  
  # Run scDblFinder
  # clusters can also be a number - like N sorted cell types
  # set dbr or not? - if 10x, fine to leave undefined.
  sce <- scDblFinder(sce, clusters = colData(sce)$merged_clusters, verbose = FALSE)
  
  # Clustered works best
  n_cells <- ncol(seurat_obj)
  n_doublets <- table(sce$scDblFinder.class)[["doublet"]]
  percentage_doublet <- round((n_doublets/n_cells) * 100, 1)
  print(glue("N cells: {n_cells}"))
  print(glue("N doublet: {n_doublets}"))
  print(glue("Doublet percentage: {percentage_doublet}"))
  # print("---------------------------------")
  
  # sce <- scDblFinder(sce, verbose = FALSE)
  #
  # n_cells <- ncol(seurat_obj)
  # n_doublets <- table(sce$scDblFinder.class)[["doublet"]]
  # percentage_doublet <- round((n_doublets/n_cells) * 100, 1)
  # print(glue("without clusters: {percentage_doublet}"))
  # print("---------------------------------")
  
  
  table(sce$scDblFinder.class)
  table(sce$scDblFinder.cluster)
  
  # Access doublets and make metadata
  doublet_metadata <- data.frame(scDblFinder.class = sce$scDblFinder.class,
                                 scDblFinder.score = sce$scDblFinder.score,
                                 scDblFinder.cluster = sce$scDblFinder.cluster,
                                 row.names = colnames(sce)
  )
  
  # Add doublet analysis to metadata
  seurat_obj <- AddMetaData(seurat_obj, doublet_metadata)
  
  # Number of singlet and doublet - Add to plot
  result <- table(seurat_obj@meta.data$scDblFinder.class, useNA = "ifany")
  
  # Plot
  DimPlot(seurat_obj, reduction = 'umap', group.by = "merged_clusters") +
    labs(
      title = glue("scDblFinder - {version}"),
      subtitle = glue("N doublets: {result[[2]]}, N singlets: {result[[1]]}"),
      caption = glue("Doublet percentage: {percentage_doublet}")
    )
  ggsave(glue("{out_dir}/{sample_name}_scDblFinder_merged_clusters.png"), width = 7, height = 6)
  
  DimPlot(seurat_obj, reduction = 'umap', group.by = "scDblFinder.class", order = TRUE) +
    labs(
      title = glue("scDblFinder - {version}"),
      subtitle = glue("N doublets: {result[[2]]}, N singlets: {result[[1]]}"),
      caption = glue("Doublet percentage: {percentage_doublet}")
    )
  ggsave(glue("{out_dir}/{sample_name}_scDblFinder.png"), width = 7, height = 6)
  
  # Cell cycle score
  seurat_obj <- CellCycleScoring(seurat_obj,
                                 s.features = s.genes,
                                 g2m.features = g2m.genes)
  
  DimPlot(seurat_obj, reduction = 'umap', group.by = "Phase", order = TRUE) +
    labs(title = glue("Phase - {version}"),
         subtitle = glue("N doublets: {result[[2]]}, N singlets: {result[[1]]}"),
         caption = glue("Doublet percentage: {percentage_doublet}"))
  ggsave(glue("{out_dir}/{sample_name}_scDblFinder_CellCyclePhase.png"), width = 7, height = 6)
  
  table(seurat_obj$Phase, seurat_obj$scDblFinder.class)
  
  seurat_obj_QC[[sample_name]] <- seurat_obj
  
}

#################### Export list of Seurat objects with QC metrices in metadata ####################

saveRDS(seurat_obj_QC, glue("07_seurat_QC/out/seurat_obj_QC_{version}.rds"))

# seurat_obj_QC$`HH117-SILP-INF-PC`$merged_clusters
# seurat_obj$merged_clusters

# doublet check
# Doublet class VS nFeature_RNA and nCount_RNA

# for (sample_name in names(seurat_obj_QC)){
#
#   # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
#   seurat_obj <- seurat_obj_QC[[sample_name]]
#
#   n_singlets <- seurat_obj$scDblFinder.class %>% str_count("singlet") %>% sum()
#   n_doublet <- seurat_obj$scDblFinder.class %>% str_count("doublet") %>% sum()
#
#   p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", group.by = "scDblFinder.class", layer = "counts") + NoLegend() +
#     labs(caption = glue("singlets: {n_singlets}\ndoublets: {n_doublet}"))
#   p2 <- VlnPlot(seurat_obj, features = "nCount_RNA", group.by = "scDblFinder.class", layer = "counts") + NoLegend()
#
#   p <- p1+p2
#
#   ggsave(glue("07_seurat_QC/plot/doublets/{sample_name}_doublets_vs_nGENES.png"), width = 8, height = 12)
#
#   FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by = "scDblFinder.class", log = TRUE)
#   ggsave(glue("07_seurat_QC/plot/doublets/{sample_name}_nFeature_vs_nCount.png"), width = 10, height = 8)
#
#
# }
