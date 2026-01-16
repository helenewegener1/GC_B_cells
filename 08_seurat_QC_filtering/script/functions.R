
compute_QC_metrics <- function(seurat_obj){
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HBA|^HBB")
  
  return(seurat_obj)
  
}

plot_qc <- function(seurat_obj, sample_name, version = "raw", filtering = ""){
  
  out_dir <- glue("08_seurat_QC_filtering/plot_{cellRversion}/{sample_name}")
  
  n_cells <- ncol(seurat_obj)
  
  p_ribo <- VlnPlot(seurat_obj, features = "percent.ribo", layer = "counts") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  p_mt <- VlnPlot(seurat_obj, features = "percent.mt", layer = "counts") + geom_hline(yintercept = 20, color = "black") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  p_feature <- VlnPlot(seurat_obj, features = "nFeature_RNA", layer = "counts") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  p_count <- VlnPlot(seurat_obj, features = "nCount_RNA", layer = "counts") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')
  
  p_final <- (p_ribo | p_mt) / (p_feature | p_count) + 
    plot_annotation(title = glue("{tools::toTitleCase(version)} QC plots of sample {sample_name}"),
                    subtitle = cellRversion, 
                    caption = glue("N cells: {n_cells}"))
  
  if (filtering != ""){
    p_final <- p_final + 
      plot_annotation(subtitle = filtering)
  }
  
  ggsave(plot = p_final,
         filename = glue("{out_dir}/{sample_name}_{version}_QC_plot.png"), 
         width = 9, 
         height = 8)
  
  return(p_final)
  
}

pre_filter_pipeline <- function(seurat_obj){

  # Remove doublets 
  # seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")
  
  # Create directory for plots of specific sample
  out_dir <- glue("08_seurat_QC_filtering/plot_{cellRversion}/{sample_name}")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Compute QC metrics
  seurat_obj <- compute_QC_metrics(seurat_obj)
  
  # Plot QC metrics in violin plots
  plot_qc(seurat_obj = seurat_obj, 
          sample_name = sample_name, 
          version = "raw", 
          filtering = "")
  
  FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt") 
  ggsave(glue("{out_dir}/{sample_name}_nFeature_vs_MT_plot.png"), width = 10, height = 8)
  
  return(seurat_obj)
  
}
