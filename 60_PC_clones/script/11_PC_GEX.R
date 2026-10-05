library(qs2)
library(Seurat)
library(glue)
library(UCell)

# ==============================================================================
# Load data 
# ==============================================================================

PC_obj <- qs_read(file = "00_data/PC_obj_sept_2026_noIgs_new_samples.qs2")

outdir <- "60_PC_clones/plot/11_PC_GEX"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Check markers 
# ==============================================================================

# Short-lived profile: CD19⁺, MKI67⁺, CD95 high, Bcl2 low.
# Long-lived profile: CD19⁻CD45⁻, Ki67⁻, CD95 lower, Bcl2 higher, CD56⁺ and CD28⁺.

markers <- c(
  CD19 = "CD19",     # BCR co-receptor, B cell marker
  CD45 = "PTPRC",    # Pan-leukocyte phosphatase, immune cell marker 
  Ki67 = "MKI67",    # Proliferation marker
  CD95 = "FAS",      # Death receptor, apoptosis sensitivity
  Bcl2 = "BCL2",     # Anti-apoptotic survival factor
  CD56 = "NCAM1",    # Niche adhesion, long-lived PCs
  CD28 = "CD28"      # Niche interaction, long-lived PCs
)

# Check missing markers 
markers[!markers %in% rownames(PC_obj)] 

Reductions(PC_obj)

DimPlot(PC_obj, reduction = "umap_harmony")

FeaturePlot(PC_obj, features = "CD19", reduction = "umap_harmony", pt.size = 2, order = T, raster=FALSE)

FeaturePlot(PC_obj, features = markers, reduction = "umap_harmony", raster=FALSE)
ggsave(glue("{outdir}/FeaturePlot.png"), width = 11, height = 7.5, dpi = 900)

FeaturePlot(PC_obj, features = markers, reduction = "umap_harmony", order = T, raster=FALSE)
ggsave(glue("{outdir}/FeaturePlot_order.png"), width = 11, height = 7.5, dpi = 900)

PC_obj$PC_clusters <- Idents(PC_obj)

VlnPlot(PC_obj, features = c("CD19", "PTPRC"), group.by = "PC_clusters", raster=FALSE)

# ==============================================================================
# Module score (with UCell)
# ==============================================================================

sigs <- list(
  shortlived = c("CD19+", "MKI67+", "FAS+", "BCL2-"),
  longlived  = c("CD19-", "PTPRC-", "MKI67-", "FAS-", "BCL2+", "NCAM1+", "CD28+")
)

subtitle <- paste0(
  "Short-lived markers: ", str_flatten(sigs$shortlived, ", "), 
  "\n", 
  "Long-lived markers: ", str_flatten(sigs$longlived, ", ")
)


PC_obj <- AddModuleScore_UCell(PC_obj, features = sigs)

VlnPlot(PC_obj, features = c("shortlived_UCell", "longlived_UCell"),
        group.by = "PC_clusters", pt.size = 0) + 
  patchwork::plot_annotation(
    subtitle = subtitle
  )

ggsave(glue("{outdir}/UCell_ModuleScore_Vln.png"), width = 9.5, height = 5)






