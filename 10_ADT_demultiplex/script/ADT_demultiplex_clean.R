getwd()

library(SeuratObject)
library(DropletUtils)
library(Seurat)
library(tidyverse)
library(glue)
library(patchwork)
library(readxl)
# library(ztable)
library(pheatmap)
library(ggsci)
library(gtools)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")

sample_names <- names(seurat_obj_nonDC_list)
sample_names

# Check samples that have ADT data available 
for (sample_name in sample_names){
  
  seurat_obj <- seurat_obj_nonDC_list[[sample_name]]
  
  # Check sample has an ADT assay
  if ("ADT" %in% names(seurat_obj@assays)) {
    print(sample_name)
  }
}

# seurat_obj_nonDC_list$`HH119-SI-PP-CD19-Pool1`@assays$ADT$counts
# seurat_obj_nonDC_list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH`@assays$ADT$counts

# ---------------------------------------------------------------------------
# Plot log_10 counts for each ADT for each sample to define zero points
# ---------------------------------------------------------------------------

seurat_obj_ADT <- list()

for (sample_name in sample_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH117-SILP-INF-PC"
  seurat_obj <- seurat_obj_nonDC_list[[sample_name]]
  
  # Check sample has an ADT assay
  if (!("ADT" %in% names(seurat_obj@assays))) {
    next
  }
  
  # DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) + NoLegend()
  # seurat_obj$scDblFinder.class %>% table()
  # seurat_obj$DC_bool %>% table()
  
  # N cells 
  n_cells <- seurat_obj %>% ncol()
  
  # Follicols 
  fols <- rownames(seurat_obj[["ADT"]]$data)

  # CLR normalization of ADT
  # seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
  
  # Save object 
  seurat_obj_ADT[[sample_name]] <- seurat_obj
  
  # Initialize directory 
  output_dir <- glue("10_ADT_demultiplex/plot/{sample_name}/log10_counts_per_fol")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Plot distribution of each ADT 
  for (fol in fols){
    
    # fol <- "Fol-3"
    seurat_obj[["ADT"]]$counts %>% t() %>% as.data.frame() %>%
      pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "counts") %>% 
      filter(ADT == fol) %>%
      ggplot(aes(x = log10(counts))) + # log counts
      geom_density() +
      # geom_histogram(alpha = 0.5, binwidth = 0.05) +
      scale_x_continuous(
        breaks = seq(0, 10, by = 0.5),        # labeled ticks
        minor_breaks = seq(0, 10, by = 0.1)   # unlabeled splits
      ) +
      theme_bw() + 
      theme(legend.position = "none") + 
      labs(
        title = fol, 
        subtitle = sample_name, 
        x = "log10 counts"
      )
    
    ggsave(glue("{output_dir}/{fol}_ADT_log10_count_values_density.png"), width = 14, height = 6)
    
  }
  
  # All together
  seurat_obj[["ADT"]]$counts %>% t() %>% as.data.frame() %>% 
    pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "counts") %>% 
    ggplot(aes(x = log(counts), fill = ADT)) + 
    geom_density(alpha = 0.5) + 
    facet_wrap(vars(ADT)) + 
    theme_bw() + 
    theme(legend.position = "none") + 
    labs(subtitle = sample_name)
  
  ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/ADT_log_counts_density.png"), width = 10, height = 7)
  
}

names(seurat_obj_ADT)

saveRDS(seurat_obj_ADT, "10_ADT_demultiplex/out/seurat_obj_ADT.rds")

# ---------------------------------------------------------------------------
# Define zero points for each sample in the sourced script
# ---------------------------------------------------------------------------

source("10_ADT_demultiplex/script/ADT_zero_points.R")

# ---------------------------------------------------------------------------
# Demultiplexing
# ---------------------------------------------------------------------------

# seurat_obj_ADT <- readRDS("10_ADT_demultiplex/out/seurat_obj_ADT.rds")
sample_names_ADT <- names(seurat_obj_ADT)

seurat_obj_ADT_demultiplexed <- list()

source("10_ADT_demultiplex/script/functions.R")

for (sample_name in sample_names_ADT){
  
  # Define object of sample
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  seurat_obj <- seurat_obj_ADT[[sample_name]]
  
  # Get raw ADT counts 
  ADT_counts <- seurat_obj@assays$ADT$counts 
  
  # Data wrangle
  ADT_counts_t <- ADT_counts %>% as.matrix() %>% t()
  
  # Divide by zero-points
  ## Keeps values positive (log works)
  ## 1 → exactly at the noise threshold
  ## >1 → above noise (signal)
  ## <1 → below noise (background)
  zero_point <- ADT_zero_point[[sample_name]]
  ADT_counts_t_corrected <- sweep(ADT_counts_t, 2, zero_point[colnames(ADT_counts_t)], FUN = "/")
  
  # Make long format and calculate log ratios
  ADT_demultiplexed <- ADT_counts_t_corrected %>%
    as.data.frame() %>% 
    rownames_to_column("Cell") %>%
    pivot_longer(
      cols = starts_with("Fol"),
      names_to = "Fol",
      values_to = "Count"
    ) %>% 
    group_by(Cell) %>%
    arrange(desc(Count), .by_group = TRUE) %>% 
    # ratio 
    mutate(
      Fol = factor(Fol, levels = Fol),
      next_Count = lead(Count),
      ratio = next_Count / Count,
      label_y = pmin(Count, next_Count) * 1.05,
      label_x = row_number() + 0.5,
      rank = row_number(),
      n_above_1 = sum(Count > 1)
    ) %>%
    # Classification
    mutate(
      top_signal = Count[rank == 1],
      ratio_12 = ratio[rank == 1],
      
      # Demultiplexing 
      manual_ADT_class = case_when(
        top_signal < 1                  ~ "Negative",  # Below the manual "zero point" --> no signal 
        n_above_1 == 1                  ~ "Singlet",   # A single ADT above "zero point" --> singlet 
        n_above_1 > 1 & ratio_12 < 0.5  ~ "Singlet",   # Multiple ADTs above "zero point" + if highest signal is twice as large as the second --> singlet 
        TRUE                            ~ "Doublet"    # Else --> doublet 
      ), 
      
      manual_ADT_ID = case_when(
        
        top_signal < 1                  ~ "Negative",  
        n_above_1 == 1                  ~ Fol[1],
        n_above_1 > 1 & ratio_12 < 0.5  ~ Fol[1],
        TRUE                            ~ "Doublet"
        
      ),
      
      manual_ADT_full_ID = case_when(
        
        top_signal < 1                  ~ "Negative",
        n_above_1 == 1                  ~ Fol[1],
        n_above_1 > 1 & ratio_12 < 0.5  ~ Fol[1],
        TRUE                            ~ paste0(sort(c(Fol[1], Fol[2]))[1],  # first ADT in canonical order
                                                 "-",
                                                 sort(c(Fol[1], Fol[2]))[2])  # second ADT in canonical order
                                                 
      )
      
    ) %>% 
    ungroup()
  
  # Stats
  ADT_demultiplexed %>% 
    select(Fol, Count) %>% 
    mutate(above_zero_point = ifelse(Count >= 1, TRUE, FALSE)) %>% 
    group_by(Fol) %>% 
    summarize(N_above_zero_point = matrixStats::count(above_zero_point)) %>% 
    ggplot(aes(x = reorder(Fol, -N_above_zero_point), y = N_above_zero_point)) + 
    geom_col() + 
    labs(
      subtitle = sample_name, 
      x = "",
      y = "N cells above zero-point"
    ) + 
    theme_bw()
  
  ggsave(
    glue("10_ADT_demultiplex/plot/{sample_name}/N_above_zero_point_{sample_name}.png"),
    width = 10, 
    height = 6
  )
  
  outdir_plot <- glue("10_ADT_demultiplex/plot/{sample_name}/individual_cell_annotation")
  dir.create(outdir_plot, recursive = TRUE, showWarnings = FALSE)
  
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 1)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 2)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 3)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 4)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 5)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 6)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 112)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 8) 
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 13) 
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 1224)
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 18) 
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 9) 
  # plot_ADT_cells(df = ADT_demultiplexed, cell_nr = 26) 
  
  for (cell_nr in sample(1:8000, size = 10)){
    plot_ADT_cells(df = ADT_demultiplexed, cell_nr = cell_nr)
  }
  
  # Clean up to prep for adding demultiplexing to seurat object 
  metadata_ADT_demultiplexed <- ADT_demultiplexed %>% 
    select(Cell, manual_ADT_class, manual_ADT_ID, manual_ADT_full_ID) %>% 
    distinct() %>% 
    column_to_rownames("Cell")
  
  # Add to seurat object
  seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_ADT_demultiplexed)
  
  seurat_obj_ADT_demultiplexed[[sample_name]] <- seurat_obj
  
}

# ---------------------------------------------------------------------------
# Demultiplexing stats 
# ---------------------------------------------------------------------------

source("10_ADT_demultiplex/script/functions.R")

sample_names_ADT <- names(seurat_obj_ADT_demultiplexed)

for (sample_name in sample_names_ADT) {
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"
  seurat_obj <- seurat_obj_ADT_demultiplexed[[sample_name]]
  
  # ------------------------------
  # Demultiplexing stats bar plots
  # ------------------------------
  plot_barplot(data.frame(manual_ADT_class = seurat_obj$manual_ADT_class), "manual_ADT_class")
  ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/barplot_ADT_class.png"), width = 10, height = 6)
  plot_barplot(data.frame(manual_ADT_ID = seurat_obj$manual_ADT_ID), "manual_ADT_ID")
  ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/barplot_ADT_ID.png"), width = 10, height = 6)
  
  # ------------------------------
  # Pattern in doublets
  # ------------------------------
  heatmap_data <- seurat_obj$manual_ADT_full_ID[str_detect(seurat_obj$manual_ADT_full_ID, "-Fol")] %>% 
    table() %>% 
    as.data.frame() %>% 
    separate(., 1, into = c("ADT_1", "ADT_2"), sep = "-(?=Fol)", remove = FALSE) %>% 
    select(!.) %>% 
    pivot_wider(names_from = ADT_1, values_from = Freq, values_fill = 0) %>% 
    column_to_rownames("ADT_2") %>% 
    as.matrix() %>% 
    .[rev(mixedsort(rownames(.))), rev(mixedsort(colnames(.)))]

  # Create heatmap
  pdf(glue("10_ADT_demultiplex/plot/{sample_name}/doublet_pattern_heatmap.pdf"), width = 8, height = 7, useDingbats = FALSE)
  pheatmap(heatmap_data, cluster_rows = F, cluster_cols = F, display_numbers = T, main = sample_name)
  dev.off()
  
  # ------------------------------
  # UMAPs
  # ------------------------------
  
  outdir_umap <- glue("10_ADT_demultiplex/plot/{sample_name}/UMAP_plots")
  dir.create(outdir_umap, recursive = TRUE, showWarnings = FALSE)

  # Dim plots
  dimplot_vars <- c("manual_ADT_class", "Phase")
  for (var in dimplot_vars){
    DimPlot(seurat_obj, group.by = var) + labs(subtitle = sample_name) 
    ggsave(glue("{outdir_umap}/UMAP_{var}.png"), width = 8, height = 8)
  }
  
  # Feature plots
  featureplot_vars <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "sce_contamination")
  for (var in featureplot_vars){
    FeaturePlot(seurat_obj, features = var) + labs(subtitle = sample_name) 
    ggsave(glue("{outdir_umap}/UMAP_{var}.png"), width = 8, height = 8)
  }
  
  # FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA") + theme_classic()
  # ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/FeatureScatter.png"), width = 8, height = 6)

  VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
          group.by = "manual_ADT_class", ncol = 2)
  ggsave(glue("{outdir_umap}/ADT_VlnPlot.png"), width = 12, height = 10)
  
  # ------------------------------
  # UMAPs and tSNEs
  # ------------------------------
  
  Idents(seurat_obj) <- "manual_ADT_class"
  
  # First, we will remove negative cells from the object
  seurat_obj.subset <- subset(seurat_obj, idents = "Negative", invert = TRUE)
  
  # Calculate a tSNE embedding of the HTO data
  DefaultAssay(seurat_obj.subset) <- "ADT"
  seurat_obj.subset <- ScaleData(seurat_obj.subset, features = rownames(seurat_obj.subset), verbose = FALSE)
  seurat_obj.subset <- RunPCA(seurat_obj.subset, features = rownames(seurat_obj.subset), approx = FALSE)
  seurat_obj.subset <- RunUMAP(seurat_obj.subset, dims = 1:8, perplexity = 100)
  seurat_obj.subset <- RunTSNE(seurat_obj.subset, dims = 1:8, perplexity = 100)
  
  DimPlot(seurat_obj.subset, reduction = "umap", group.by = "manual_ADT_ID", label = TRUE) + NoLegend() + labs(subtitle = sample_name)
  ggsave(glue("{outdir_umap}/UMAP_ADT_ID.png"), width = 9, height = 7)
  DimPlot(seurat_obj.subset, reduction = "umap", group.by = "manual_ADT_class") + labs(subtitle = sample_name)
  ggsave(glue("{outdir_umap}/UMAP_ADT_class.png"), width = 9, height = 7)
  
  DimPlot(seurat_obj.subset, reduction = "tsne", group.by = "manual_ADT_ID", label = TRUE) + NoLegend() + labs(subtitle = sample_name)
  ggsave(glue("{outdir_umap}/tSNE_ADT_ID.png"), width = 9, height = 7)
  DimPlot(seurat_obj.subset, reduction = "tsne", group.by = "manual_ADT_class") + labs(subtitle = sample_name)
  ggsave(glue("{outdir_umap}/tSNE_ADT_class.png"), width = 9, height = 7)
  
  # ------------------------------
  # Rigde plot 
  # ------------------------------
  
  Idents(seurat_obj) <- "manual_ADT_ID"
  RidgePlot(seurat_obj, assay = "ADT", features = rownames(seurat_obj[["ADT"]]))
  ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/RidgePlot.png"), width = 16, height = 20)
  
  # ------------------------------
  # Compare to HTODemux
  # ------------------------------

  # Demultiplexing with HTODemux
  seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
  seurat_obj <- HTODemux(seurat_obj, assay = "ADT", positive.quantile = 0.999)
  seurat_obj$HTODemux_ADT_ID <- lapply(seurat_obj@meta.data$ADT_classification, function(x) {ifelse(str_detect(x, "_"), "Doublet", x)}) %>% unlist()
  # seurat_obj$HTODemux_ADT_ID %>% table()
  
  # Compare demultiplexing in heatmap
  seurat_obj[[]] %>% 
    count(manual_ADT_ID, HTODemux_ADT_ID) %>% 
    ggplot(aes(x = HTODemux_ADT_ID, y = manual_ADT_ID, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n), color = "black") +
    scale_fill_gradient(low = "lightgrey", high = "red") +
    theme_minimal() + 
    labs(
      subtitle = sample_name
    )
  
  ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/compare_manual_HTODemux.png"), width = 10, height = 6.5)
  
}

# ---------------------------------------------------------------------------
# Make and export final list of seurat objects 
# ---------------------------------------------------------------------------

not_ADT_samples <- names(seurat_obj_nonDC_list)[!names(seurat_obj_nonDC_list) %in% names(seurat_obj_ADT_demultiplexed)]

seurat_obj_ADT_demultiplexed_all <- c(seurat_obj_nonDC_list[not_ADT_samples], seurat_obj_ADT_demultiplexed)

saveRDS(seurat_obj_ADT_demultiplexed_all, "10_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds")
