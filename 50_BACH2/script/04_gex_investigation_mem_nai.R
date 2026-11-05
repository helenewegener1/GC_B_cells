library(tidyverse)
library(Seurat)
library(glue)
library(qs2)
library(GeneNMF)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
source("10_broad_annotation/script/color_palette.R")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# THIS IS ONLY OUR SAMPLES. GET NAIVE AND MEMORY FROM ATLAS IF WE WANT TO TEST WITH ALL THE CELLS.
# --> BALANCE CELL TYPES --> WE DON'T WANT THE NAIVE AND MEMORY TO BE 90% OF THE CELLS 

full_obj <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs_annotated.rds")
gc_obj <- qs_read("00_data/GCB_cells_comb_wo_119pop_Oct_2026_annotated.qs2")

gc_obj$L3_GCB_annotation %>% unique()
gc_obj$L6_GCB_annotation %>% unique()

# Add GC annotations to full_obj
gc_meta <- gc_obj[[]] %>% select(L3_GCB_annotation, L6_GCB_annotation)
table(rownames(gc_meta) %in% colnames(full_obj))
full_obj <- AddMetaData(full_obj, gc_meta)

# Add new combined annotation columns
full_obj[[]] <- full_obj[[]] %>% 
  mutate(
    L1_L3_annotations = case_when(
      !is.na(L3_GCB_annotation) ~ L3_GCB_annotation, 
      .default = L1_annotation 
    ),
    L1_L6_annotations = case_when(
      !is.na(L6_GCB_annotation) ~ L6_GCB_annotation, 
      .default = L1_annotation 
    )
  )

# Create outdit
outdir <- glue("50_BACH2/plot/04_gex_mem_nai")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Quick checks 
# ------------------------------------------------------------------------------

L1_L3_colors <- c(L1_colors, L3_GC_pretty_colors)

DimPlot(full_obj)

DimPlot(full_obj, group.by = "L1_L3_annotations", cols = L1_L3_colors)
ggsave(glue("{outdir}/Dimplot_full_L1_L3.png"), width = 9, height = 6)

# ------------------------------------------------------------------------------
# Subset to GC_B_cell subsets (excluding large clone labeled as GC B cell in L1 and NA in L3) + memory + naive 
# ------------------------------------------------------------------------------

full_obj$L1_L3_annotations %>% table()

cell_subset <- c(
  levels(full_obj$L3_GCB_annotation), 
  "Memory_B_cells", "Naive_B_cells", "Unconventional_B_cells"
)

subset_obj <- subset(full_obj, subset = L1_L3_annotations %in% cell_subset)

# Reintegrate 
# Seurat workflow 
res <- 0.3
subset_obj <- NormalizeData(subset_obj, verbose = FALSE)
subset_obj <- FindVariableFeatures(subset_obj, verbose = FALSE)
subset_obj <- ScaleData(subset_obj, verbose = FALSE)
subset_obj <- RunPCA(subset_obj, verbose = FALSE)

DefaultAssay(subset_obj)
Reductions(subset_obj)
subset_obj_integrated <- subset_obj

DefaultAssay(subset_obj_integrated) <- "RNA"
Layers(subset_obj_integrated[["RNA"]])
Idents(subset_obj_integrated) %>% unique()

subset_obj_integrated <- IntegrateLayers(
  object = subset_obj_integrated, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "subset_RNA_harmony",
  # group.by.vars = "patient",
  assay = "RNA", 
  verbose = FALSE
)

subset_obj_integrated <- FindNeighbors(subset_obj_integrated, reduction = "subset_RNA_harmony", dims = 1:30)
subset_obj_integrated <- FindClusters(subset_obj_integrated, resolution = 0.3)
subset_obj_integrated <- RunUMAP(subset_obj_integrated, reduction = "pca", dims = 1:10, verbose = FALSE, reduction.name = "subset_umap.harmony")

Reductions(subset_obj_integrated)

DimPlot(subset_obj_integrated, group.by = "L1_L3_annotations", cols = L1_L3_colors, reduction = "subset_umap.harmony")
ggsave(glue("{outdir}/Dimplot_subset_L1_L3.png"), width = 9, height = 6)

DimPlot(subset_obj_integrated, group.by = "L1_L3_annotations", cols = L1_L3_colors, split.by = "patient", reduction = "subset_umap.harmony")
ggsave(glue("{outdir}/Dimplot_subset_L1_L3_split.png"), width = 15, height = 4.5)

# DotPlot(gc_obj %>% subset(L6_GCB_annotation %in% names(L3_GC_pretty_colors)), )

 # Plot features 
features <- c(
  "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", 
  "BACH2", "PRDM1", "MKI67", "AIDA", "PAX5", "BCL6", "XBP1", "SDC1"
)

grep("PRDM1", rownames(subset_obj_integrated), value = TRUE)
# Bach2 expression: Repressing Blimp1 and the plasma cell fate.
# CSR Factor Expression: AID, Pax5, Bcl6. 
# Absence of Plasma Cell Markers: Low expression of Blimp1, CD138, Xbp1.

VlnPlot(subset_obj_integrated, features = features, pt.size = 0, group.by = "L1_L3_annotations", cols = L1_L3_colors)
ggsave(glue("{outdir}/VlnPlot_subset.png"), width = 12, height = 8)

# ------------------------------------------------------------------------------
# Investigate mitochondrial content 
# ------------------------------------------------------------------------------

# # Find MT genes
# mt_genes <- grep("^MT-", rownames(gc_obj), value = TRUE)
# VlnPlot(gc_obj, features = mt_genes, pt.size = 0)
# 
# # Long data frame: one row per cell × gene, grouped by active identities
# mt_df <- FetchData(gc_obj, vars = mt_genes, layer = "data", assay = "RNA") %>%
#   mutate(ident = Idents(gc_obj)) %>%
#   pivot_longer(-ident, names_to = "gene", values_to = "expr") %>%
#   mutate(gene = factor(gene, levels = mt_genes))
# 
# ggplot(mt_df, aes(x = ident, y = expr, fill = ident)) +
#   geom_violin(scale = "width", alpha = 0.7, colour = NA) +
#   geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
#   facet_wrap(~ gene, ncol = 4) +
#   labs(x = NULL, y = NULL, title = "Mitochondrial genes") +
#   theme_classic() +
#   theme(
#     plot.title       = element_text(face = "bold", size = 20, hjust = 0.5),
#     strip.background = element_blank(),
#     strip.text       = element_text(face = "bold", hjust = 0.5),
#     axis.text.x      = element_text(angle = 45, hjust = 1),
#     legend.position  = "none"
#   )
# ggsave(glue("{outdir}/MT_genes.png"), width = 14, height = 7.5)
# 
# # Add MT module score
# mt_genes_flat <- str_flatten(mt_genes, collapse = ", ")
# gc_obj <- AddModuleScore(gc_obj, features = list(mt_genes), name = "MT_score") # stored as gc_obj$MT_score1
# FeaturePlot(gc_obj, features = "MT_score1")
# VlnPlot(gc_obj, features = c(features, "MT_score1"), pt.size = 0, ncol = 5) + plot_annotation(title = glue("MT module score ({mt_genes_flat})"))
# ggsave(glue("{outdir}/MT_module_score.png"), width = 13, height = 8)
# 
# plot_vars <- c(features, "MT_score1")
# 
# vln_df <- FetchData(gc_obj, vars = plot_vars, layer = "data", assay = "RNA") %>%
#   mutate(ident = Idents(gc_obj)) %>%
#   pivot_longer(-ident, names_to = "feature", values_to = "expr") %>%
#   mutate(feature = factor(feature, levels = plot_vars))
# 
# ggplot(vln_df, aes(x = ident, y = expr, fill = ident)) +
#   geom_violin(scale = "width", alpha = 0.7, colour = NA) +
#   geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
#   facet_wrap(~ feature, ncol = 5, scales = "free_y") +
#   labs(x = NULL, y = NULL) +
#   theme_classic() +
#   theme(
#     plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
#     strip.background = element_blank(),
#     strip.text       = element_text(face = "bold"),
#     axis.text.x      = element_text(angle = 45, hjust = 1),
#     legend.position  = "none"
#   )
# ggsave(glue("{outdir}/MT_module_score_pretty.png"), width = 16, height = 7)
# 


# ------------------------------------------------------------------------------
# NMF gene programs for BACH2
# ------------------------------------------------------------------------------

# Subset to population and reintegrate 
# subset_obj <- subset(gc_obj, L6_GCB_annotation %in% c("BACH2+", "BACH2_cycling"))
# subset_obj <- subset(gc_obj, L6_GCB_annotation %in% c("BACH2+"))

subset_obj_integrated$L1_L3_annotations %>% table()

# version <- "non_cycling"
# non_cycling_subset <- c("BACH2+", "DZ", "LZ", "LZ_transitioning", "Pre_MBC", "Pre_PB", "Memory_B_cells", "Naive_B_cells", "Unconventional_B_cells")
# subset_obj <- subset(subset_obj_integrated, L1_L3_annotations %in% non_cycling_subset)
version <- "cycling"
cycling_subset <- c("BACH2_cycling", "DZ_cycling", "LZ_cycling", "Pre_MBC_cycling", "Memory_B_cells", "Naive_B_cells", "Unconventional_B_cells")
subset_obj <- subset(subset_obj_integrated, L1_L3_annotations %in% cycling_subset)

subset_obj$L1_L3_annotations %>% table()

# Downsample to have more balanced populations 
set.seed(123)
max_cells <- 3000   # cap per cell type, about the size of LZ

cells_keep <- subset_obj@meta.data %>% 
  rownames_to_column("cell") %>% 
  group_by(L1_L3_annotations) %>% 
  slice_sample(n = max_cells) %>%          # takes all cells if a group has fewer than max_cells
  ungroup() %>% 
  pull(cell)

subset_obj_balanced <- subset(subset_obj, cells = cells_keep)
subset_obj_balanced$L1_L3_annotations %>% table()

# Split object on patient variable
# subset_obj[["RNA"]] <- split(subset_obj[["RNA"]], f = subset_obj[["patient"]] %>% pull())

Layers(subset_obj_balanced[["RNA"]]) # Check that object is split
DefaultAssay(subset_obj_balanced) <- "RNA"

# Seurat workflow 
res <- 0.3
subset_obj_balanced <- NormalizeData(subset_obj_balanced, verbose = FALSE)
subset_obj_balanced <- FindVariableFeatures(subset_obj_balanced, verbose = FALSE)
subset_obj_balanced <- ScaleData(subset_obj_balanced, verbose = FALSE)
subset_obj_balanced <- RunPCA(subset_obj_balanced, verbose = FALSE)

DefaultAssay(subset_obj_balanced)
Reductions(subset_obj_balanced)
subset_obj_int <- subset_obj_balanced

DefaultAssay(subset_obj_int) <- "RNA"
Layers(subset_obj_int[["RNA"]])
Idents(subset_obj_int) %>% unique()

subset_obj_int <- IntegrateLayers(
  object = subset_obj_int, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = glue("{version}_RNA_harmony"),
  # group.by.vars = "patient",
  assay = "RNA", 
  verbose = FALSE
)

subset_obj_int <- FindNeighbors(subset_obj_int, reduction = glue("{version}_RNA_harmony"), dims = 1:30)
subset_obj_int <- FindClusters(subset_obj_int, resolution = 0.3)
subset_obj_int <- RunUMAP(subset_obj_int, reduction = "pca", dims = 1:10, verbose = FALSE, reduction.name = glue("{version}_umap.harmony"))

# Plot
DimPlot(
  subset_obj_int, reduction = glue("{version}_umap.harmony"), 
  group.by = "L1_L3_annotations", split.by = "patient", 
  ncol = 4, cols = L1_L3_colors
) + 
  plot_annotation(
    caption = glue("Cell types downsampled to at most {max_cells} cells each")
  )
ggsave(glue("{outdir}/Dimplot_{version}_L1_L3_split.png"), width = 15, height = 4.5) 
# DimPlot(subset_obj_int, reduction = "bach2_umap.harmony", split.by = "Patient", ncol = 5)

# Split 
gc_obj.list <- SplitObject(subset_obj_int, split.by = "patient")

# Run NMF across a range of k
geneNMF.programs <- multiNMF(gc_obj.list, assay="RNA", slot="data", k=2:7, nfeatures = 1000)

# Combine gene programs identified over multiple samples and numbers of k into 
# metaprograms (MPs) / consensus programs that are robustly identified across NMF runs.
nMP <- 8

geneNMF.metaprograms <- getMetaPrograms(
  geneNMF.programs,
  nMP=nMP,
  weight.explained = 0.7,
  max.genes=100
)

# Heatmap
anno_colors <- brewer.pal(n=nMP, name="Paired")
names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

ph <- plotMetaPrograms(geneNMF.metaprograms, annotation_colors = anno_colors)
ggsave(glue("{outdir}/{version}_heatmap.png"), ph, width = 8, height = 8)

# Signature scores for gene programs
mp.genes <- geneNMF.metaprograms$metaprograms.genes
subset_obj_int <- UCell::AddModuleScore_UCell(subset_obj_int, features = mp.genes, assay="RNA", ncores=4, name = "")

VlnPlot(
  subset_obj_int, features=names(mp.genes),
  pt.size = 0, ncol=4, group.by = "L1_L3_annotations", cols = L1_L3_colors
)

for (category in c("H", "C5")){
  
  # category <- "C5"
  
  # Intepretation of gene programs by overrepresentation analysis
  top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
    runGSEA(program, universe=rownames(subset_obj_int), category = category)
  })
  
  top_p
  
  top_p_to_MP <- lapply(names(top_p), function(MP){
    
    # MP <- "MP2"
    
    top_p[[MP]] %>% 
      select(pathway) %>% 
      head(3) %>% 
      pull() %>% 
      paste0(collapse = "\n")
    
  }) %>% 
    setNames(names(top_p))
  
  df_top_p_to_MP <- data.frame(
    MP = names(top_p_to_MP), 
    pathway = top_p_to_MP %>% unlist()
  ) %>% 
    mutate(
      MP_top_p = glue("{MP}\n{pathway}")
    )
  
  # Boxplot of population and MPs
  df_plot <- subset_obj_int@meta.data %>% select(names(mp.genes), L1_L3_annotations) %>% 
    pivot_longer(cols = starts_with("MP"), names_to = "MP") %>% 
    left_join(df_top_p_to_MP, by = "MP") 
  
  # Boxplot of population and MPs
  df_plot %>% 
    ggplot(aes(x = L1_L3_annotations, y = value, color = L1_L3_annotations)) + 
    geom_boxplot() + 
    theme_bw() + 
    facet_wrap(vars(MP_top_p)) +
    scale_color_manual(values = L1_L3_colors) + 
    labs(
      title = glue("{version}: ORA category: {category}"),
      y = "MP module score",
      subtitle = glue("Cell types downsampled to at most {max_cells} cells each")
      # x = "Module"
    ) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) 
  
  # ggsave(glue("{outdir}/BACH2_BACH2_cycling_MP{nMP}_ORA_{category}.png"), width = 13, height = 6)
  ggsave(glue("{outdir}/{version}_MP{nMP}_ORA_{category}.png"), width = 13, height = 8)
  
}
 
