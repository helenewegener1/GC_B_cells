library(tidyverse)
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

# Plot features 
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "BACH2")
VlnPlot(gc_obj, features = features, pt.size = 0)

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
VlnPlot(gc_obj, features = c(features, "MT_score1"), pt.size = 0) + plot_annotation(title = glue("MT module score ({mt_genes_flat})"))
ggsave(glue("{outdir}/MT_module_score.png"), width = 13, height = 7.5)

# ------------------------------------------------------------------------------
# NMF gene programs for BACH2
# ------------------------------------------------------------------------------

# Subset to population and reintegrate 
gc_obj_bach2 <- subset(gc_obj, L6_GCB_annotation == "BACH2+")

# Split object on patient variable
gc_obj_bach2[["RNA"]] <- split(gc_obj_bach2[["RNA"]], f = gc_obj_bach2[["patient_id"]] %>% pull())

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

# DO THIS!!!
seurat_obj_integrated <- FindNeighbors(seurat_obj_integrated, reduction = "subset_RNA_harmony", dims = 1:30)
seurat_obj_integrated <- FindClusters(seurat_obj_integrated, resolution = 0.5)

Reductions(gc_obj_bach2_integrated)
DefaultDimReduc(gc_obj_bach2_integrated)
DimPlot(gc_obj_bach2_integrated, reduction = "bach2_umap.unintegrated", group.by = "patient_id")
DimPlot(gc_obj_bach2_integrated, reduction = "umap_harmony", group.by = "patient_id")

# Split 
gc_obj.list <- SplitObject(gc_obj_bach2_integrated, split.by = "patient_id")

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
ph

# Signature scores for gene programs
mp.genes <- geneNMF.metaprograms$metaprograms.genes
gc_obj_bach2_integrated <- AddModuleScore_UCell(gc_obj_bach2_integrated, features = mp.genes, assay="RNA", ncores=4, name = "")

VlnPlot(gc_obj_bach2_integrated, features=names(mp.genes),
        pt.size = 0, ncol=5)

# Intepretation of gene programs by overrepresentation analysis
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(gc_obj_bach2_integrated), category = "C5")
})

head(top_p$MP4)

top_p_to_MP <- top_p$MP4 %>% 
  select(pathway) %>% 
  rownames_to_column("MP") %>% 
  mutate(
    MP = glue("MP{MP}") %>% as.character(), 
    MP_top_p = glue("{MP}: {pathway}")
  ) 

# Boxplot of population and MPs
df_plot <- gc_obj_bach2_integrated@meta.data %>% select(names(mp.genes), L6_GCB_annotation) %>% 
  pivot_longer(cols = starts_with("MP"), names_to = "MP") %>% 
  left_join(top_p_to_MP, by = "MP")

# Boxplot of population and MPs
df_plot %>% 
  ggplot(aes(x = MP_top_p, y = value)) + 
  geom_boxplot() + 
  theme_bw() + 
  labs(
    # title = glue("{population}"),
    y = "MP module score",
    x = "Module"
  ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



