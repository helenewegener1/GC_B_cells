############################ GINA TEST REGRESS ############################ 
library(scater)

seurat_obj <- seurat_obj_QC_filtered_singlets_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]
# prior, do Normalize + FindVariableFeatures

sce <- SingleCellExperiment(assays = list(
  logcounts = as.matrix(GetAssayData(seurat_obj, slot = "data")))
)
# Add metadata
colData(sce) <- DataFrame(seurat_obj@meta.data)

# The percentage of variance that is explained by cell cycle phase
vars <- c("Phase")
var_explained_matrix <- getVarianceExplained(sce, variables = vars)

# Look at the genes where > 5 % of the variance is explained by phase
var_explained_matrix %>% as.data.frame() %>% rownames_to_column("Gene") %>%
  arrange(desc(Phase)) %>% filter(Phase > 5) %>% nrow()

# Extract them
genes_rm <- var_explained_matrix %>% as.data.frame() %>% rownames_to_column("Gene") %>%
  arrange(desc(Phase)) %>% filter(Phase > 5) %>% pull(Gene)

# How many of these genes are in the variable features
VariableFeatures(seurat_obj) %in% genes_rm %>% table()

# Remove them from the variable features
VariableFeatures(seurat_obj) <- setdiff(VariableFeatures(seurat_obj), genes_rm)

seurat_obj <- ScaleData(seurat_obj, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score"))

DefaultAssay(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj,  dims = 1:n_dim, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:n_dim, verbose = FALSE)
