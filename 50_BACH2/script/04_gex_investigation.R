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

gc_obj <- qs_read("00_data/GCB_cells_comb_wo_119pop_Oct_2026_annotated.qs2")

gc_obj$BACH2_bool <- Idents(gc_obj) == "BACH2+"
gc_obj$BACH2cyc_bool <- Idents(gc_obj) == "BACH2_cycling"

# Create outdit
outdir <- glue("50_BACH2/plot/04_gex")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Quick checks 
# ------------------------------------------------------------------------------

DimPlot(gc_obj)

DimPlot(gc_obj, group.by = "BACH2_bool", cols = c("grey", L3_GC_pretty_colors[["BACH2+"]]))
DimPlot(gc_obj, group.by = "BACH2cyc_bool", cols = c("grey", L3_GC_pretty_colors[["BACH2_cycling"]]))

DimPlot(gc_obj %>% subset(L6_GCB_annotation %in% names(L3_GC_pretty_colors)), group.by = "L6_GCB_annotation", split.by = "Patient", ncol = 7, cols = L3_GC_pretty_colors)
ggsave(glue("{outdir}/Dimplot_L6_split.png"), width = 13, height = 5)

DimPlot(gc_obj, group.by = "L6_GCB_annotation", cols = L3_GC_pretty_colors)
ggsave(glue("{outdir}/Dimplot_L6.png"), width = 10, height = 6)

# DotPlot(gc_obj %>% subset(L6_GCB_annotation %in% names(L3_GC_pretty_colors)), )

 # Plot features 
features <- c(
  "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", 
  "BACH2", "PRDM1", "MKI67", "AIDA", "PAX5", "BCL6", "XBP1", "SDC1"
)

grep("PRDM1", rownames(gc_obj), value = TRUE)
# Bach2 expression: Repressing Blimp1 and the plasma cell fate.
# CSR Factor Expression: AID, Pax5, Bcl6. 
# Absence of Plasma Cell Markers: Low expression of Blimp1, CD138, Xbp1.

VlnPlot(gc_obj, features = features, pt.size = 0)

# ------------------------------------------------------------------------------
# Quick checks 
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Investigate mitochondrial content 
# ------------------------------------------------------------------------------

# Find MT genes
mt_genes <- grep("^MT-", rownames(gc_obj), value = TRUE)
VlnPlot(gc_obj, features = mt_genes, pt.size = 0)

# Long data frame: one row per cell × gene, grouped by active identities
mt_df <- FetchData(gc_obj, vars = mt_genes, layer = "data", assay = "RNA") %>%
  mutate(ident = Idents(gc_obj)) %>%
  pivot_longer(-ident, names_to = "gene", values_to = "expr") %>%
  mutate(gene = factor(gene, levels = mt_genes))

ggplot(mt_df, aes(x = ident, y = expr, fill = ident)) +
  geom_violin(scale = "width", alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  facet_wrap(~ gene, ncol = 4) +
  labs(x = NULL, y = NULL, title = "Mitochondrial genes") +
  theme_classic() +
  theme(
    plot.title       = element_text(face = "bold", size = 20, hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", hjust = 0.5),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "none"
  )
ggsave(glue("{outdir}/MT_genes.png"), width = 14, height = 7.5)

# Add MT module score
mt_genes_flat <- str_flatten(mt_genes, collapse = ", ")
gc_obj <- AddModuleScore(gc_obj, features = list(mt_genes), name = "MT_score") # stored as gc_obj$MT_score1
FeaturePlot(gc_obj, features = "MT_score1")
VlnPlot(gc_obj, features = c(features, "MT_score1"), pt.size = 0, ncol = 5) + plot_annotation(title = glue("MT module score ({mt_genes_flat})"))
ggsave(glue("{outdir}/MT_module_score.png"), width = 13, height = 8)

plot_vars <- c(features, "MT_score1")

vln_df <- FetchData(gc_obj, vars = plot_vars, layer = "data", assay = "RNA") %>%
  mutate(ident = Idents(gc_obj)) %>%
  pivot_longer(-ident, names_to = "feature", values_to = "expr") %>%
  mutate(feature = factor(feature, levels = plot_vars))

ggplot(vln_df, aes(x = ident, y = expr, fill = ident)) +
  geom_violin(scale = "width", alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  facet_wrap(~ feature, ncol = 5, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "none"
  )
ggsave(glue("{outdir}/MT_module_score_pretty.png"), width = 16, height = 7)



# ------------------------------------------------------------------------------
# NMF gene programs for BACH2
# ------------------------------------------------------------------------------

# Subset to population and reintegrate 
# gc_obj_bach2 <- subset(gc_obj, L6_GCB_annotation %in% c("BACH2+", "BACH2_cycling"))
# gc_obj_bach2 <- subset(gc_obj, L6_GCB_annotation %in% c("BACH2+"))

gc_obj$L6_GCB_annotation %>% table()

version <- "BACH2"
gc_obj_bach2 <- subset(gc_obj, L6_GCB_annotation %in% c("BACH2+", "DZ", "LZ", "Pre_MBC", "PB"))

# version <- "BACH2_cycling"
# gc_obj_bach2 <- subset(gc_obj, L6_GCB_annotation %in% c("BACH2_cycling", "DZ_cycling", "LZ_cycling", "Pre_MBC_cycling"))

gc_obj_bach2$L2_GCB_annotation %>% table()

# Split object on patient variable
gc_obj_bach2[["RNA"]] <- split(gc_obj_bach2[["RNA"]], f = gc_obj_bach2[["Patient"]] %>% pull())

Layers(gc_obj_bach2[["RNA"]]) # Check that object is split
DefaultAssay(gc_obj_bach2) <- "RNA"

# Seurat workflow 
res <- 0.3
gc_obj_bach2 <- NormalizeData(gc_obj_bach2, verbose = FALSE)
gc_obj_bach2 <- FindVariableFeatures(gc_obj_bach2, verbose = FALSE)
gc_obj_bach2 <- ScaleData(gc_obj_bach2, verbose = FALSE)
gc_obj_bach2 <- RunPCA(gc_obj_bach2, verbose = FALSE)
gc_obj_bach2 <- FindNeighbors(gc_obj_bach2, dims = 1:10, verbose = FALSE)
gc_obj_bach2 <- FindClusters(gc_obj_bach2, resolution = res, verbose = FALSE, cluster.name = "bach2_unintegrated_clusters")
gc_obj_bach2 <- RunUMAP(gc_obj_bach2, reduction = "pca", dims = 1:10, verbose = FALSE, reduction.name = "bach2_umap.unintegrated")

DefaultAssay(gc_obj_bach2)
Reductions(gc_obj_bach2)
gc_obj_bach2_integrated <- gc_obj_bach2

DefaultAssay(gc_obj_bach2_integrated) <- "RNA"
Layers(gc_obj_bach2_integrated[["RNA"]])
Idents(gc_obj_bach2_integrated) %>% unique()

gc_obj_bach2_integrated <- IntegrateLayers(
  object = gc_obj_bach2_integrated, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "bach2_RNA_harmony",
  # group.by.vars = "patient",
  assay = "RNA", 
  verbose = FALSE
)

gc_obj_bach2_integrated <- FindNeighbors(gc_obj_bach2_integrated, reduction = "bach2_RNA_harmony", dims = 1:30)
gc_obj_bach2_integrated <- FindClusters(gc_obj_bach2_integrated, resolution = 0.3)
gc_obj_bach2_integrated <- RunUMAP(gc_obj_bach2_integrated, reduction = "pca", dims = 1:10, verbose = FALSE, reduction.name = "bach2_umap.harmony")

# Plot
DimPlot(gc_obj_bach2_integrated, reduction = "bach2_umap.harmony", group.by = "L6_GCB_annotation", split.by = "Patient", ncol = 5, cols = L3_GC_pretty_colors)
# DimPlot(gc_obj_bach2_integrated, reduction = "bach2_umap.harmony", split.by = "Patient", ncol = 5)

# Split 
gc_obj.list <- SplitObject(gc_obj_bach2_integrated, split.by = "Patient")

# Run NMF across a range of k
geneNMF.programs <- multiNMF(gc_obj.list, assay="RNA", slot="data", k=2:7, nfeatures = 1000)

# Combine gene programs identified over multiple samples and numbers of k into 
# metaprograms (MPs) / consensus programs that are robustly identified across NMF runs.
nMP <- 6

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
ph
ggsave(glue("{outdir}/{version}_heatmap.png"), ph, width = 8, height = 8)

# Signature scores for gene programs
mp.genes <- geneNMF.metaprograms$metaprograms.genes
gc_obj_bach2_integrated <- UCell::AddModuleScore_UCell(gc_obj_bach2_integrated, features = mp.genes, assay="RNA", ncores=4, name = "")

VlnPlot(gc_obj_bach2_integrated, features=names(mp.genes),
        pt.size = 0, ncol=5, group.by = "L6_GCB_annotation")


for (category in c("H", "C5")){
  
  # category <- "C5"
  
  # Intepretation of gene programs by overrepresentation analysis
  top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
    runGSEA(program, universe=rownames(gc_obj_bach2_integrated), category = category)
  })
  
  top_p
  
  top_p_to_MP <- lapply(names(top_p), function(MP){
    
    # MP <- "MP2"
    
    top_p[[MP]] %>% 
      select(pathway) %>% 
      head(5) %>% 
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
  df_plot <- gc_obj_bach2_integrated@meta.data %>% select(names(mp.genes), L6_GCB_annotation) %>% 
    pivot_longer(cols = starts_with("MP"), names_to = "MP") %>% 
    left_join(df_top_p_to_MP, by = "MP") 
  
  # Boxplot of population and MPs
  df_plot %>% 
    ggplot(aes(x = L6_GCB_annotation, y = value)) + 
    geom_boxplot() + 
    theme_bw() + 
    facet_wrap(vars(MP_top_p)) +
    labs(
      title = glue("ORA category: {category}"),
      y = "MP module score"
      # subtitle = "Only focus on MP1"
      # x = "Module"
    ) 
  
  # ggsave(glue("{outdir}/BACH2_BACH2_cycling_MP{nMP}_ORA_{category}.png"), width = 13, height = 6)
  ggsave(glue("{outdir}/{version}_MP{nMP}_ORA_{category}.png"), width = 13, height = 6)
  
}

