getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
library(scRepertoire)

# Load data
combined.TCR.filtered <- readRDS("20_VDJ/out/combined.TCR.filtered.rds")

# Quantifying Unique Clones
# clonalQuant() returns the total or relative numbers of unique clones
for (group in c("sample_high_level", "patient")){
  
  # group <- "sample_high_level"
  
  clonalQuant(
    combined.TCR.filtered, 
    cloneCall="strict",
    chain = "both", 
    group.by = group, 
    scale = TRUE # relative percentage of unique clones scaled by the total repertoire size
  ) + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 80, hjust = 1)) + 
    labs(title = "TCR: Relative numbers of unique clones", 
         x = "")
  
  ggsave(glue("20_VDJ/plot/TCR_scRepertoire/clonalQuant_{group}.png"), width = 8, height = 10)
  
}

# Compare clones within each patient 
res <- lapply(
  
  c("HH117", "HH119"), function(x) {
    
    # x <- "HH119"
    
    HH_mask <- grep(x, names(combined.TCR.filtered))
    combined.TCR.filtered_HH <- combined.TCR.filtered[HH_mask]
    # names(combined.TCR.filtered_HH119)
    
    clonalCompare(combined.TCR.filtered_HH, 
                  top.clones = 10, 
                  group.by = "sample_high_level",
                  cloneCall="aa", 
                  graph = "alluvial") + 
      labs(title = glue("TCR: Compare {x}")) + 
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(glue("20_VDJ/plot/TCR_scRepertoire/clonalCompare_{x}.png"), width = 13, height = 7.5)
    
    return(combined.TCR.filtered_HH)
    
  }
)

combined.TCR.filtered_HH117 <- res[[1]]
names(combined.TCR.filtered_HH117)

combined.TCR.filtered_HH119 <- res[[2]]
names(combined.TCR.filtered_HH119)

sample_names <- names(combined.TCR.filtered) %>% str_split_i("_", 1) %>% unique()

# Split by samples since too many sample in one plot with each follicle. 
for (sample_name in sample_names) {
  
  # sample_name <- "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH"
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  ############################## clonalAbundance ############################### 
  # total number of clones at specific frequencies within a sample or group
  # General: There's one clone that has the highest abundance within each sample. 
  
  combined.TCR_subset <- combined.TCR.filtered[startsWith(names(combined.TCR.filtered), sample_name)]
  
  outdir <- glue("20_VDJ/plot/TCR_scRepertoire/{sample_name}")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (length(names(combined.TCR_subset)) > 1){
    
    # Exclude Doublet and Negative subsets
    names_fols_mask <- grep("Doublet|Negative", names(combined.TCR_subset), invert = TRUE)
    combined.TCR_subset <- combined.TCR_subset[names_fols_mask]
    # names(combined.TCR_subset)
    
    # Sort naturally by the numeric part
    sorted_x <- mixedsort(names(combined.TCR_subset), decreasing = TRUE)
    combined.TCR_subset <- combined.TCR_subset[sorted_x]
    names(combined.TCR_subset)
    
    n_cells <- lapply(combined.TCR_subset, function(x) nrow(x))
    
    labels <- combined.TCR_subset[startsWith(names(combined.TCR_subset), sample_name)] %>% names() %>% str_split_i("_", 2)
    
    labels_final <- paste(labels, "-", n_cells, "cells")
    
  } else if (length(names(combined.TCR_subset)) == 1) {
    
    labels <- sample_name
    
    n_cells <- nrow(combined.TCR_subset)
    
    labels_final <- paste(labels, "-", n_cells, "cells")
    
  }
  
  # Not scaled
  clonalAbundance(combined.TCR_subset, 
                  cloneCall = "gene", 
                  scale = FALSE) + 
    labs(title = "TCR: Clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    # scale_color_discrete(labels = labels) +
    scale_color_discrete(labels = labels_final) +
    theme_dark()
  
  ggsave(glue("{outdir}/{sample_name}_clonalAbundance.png"), width = 10, height = 6.5)
  
  # Scaled 
  clonalAbundance(combined.TCR_subset, 
                  cloneCall = "gene", 
                  scale = TRUE) + 
    labs(title = "TCR: Scaled clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    # scale_color_discrete(labels = labels) + 
    scale_color_discrete(labels = labels_final) + 
    theme_dark()
  
  ggsave(glue("{outdir}/{sample_name}_clonalAbundance_scaled.png"), width = 10, height = 6.5)
  
  # Distribution of Sequence Lengths: looking at the length distribution of the CDR3 sequences
  for (chain in c("TRA", "TRB")){
    # chain <- "TRA"
    clonalLength(combined.TCR_subset, 
                 cloneCall="aa", 
                 chain = chain, 
                 scale = TRUE) + 
      labs(title = glue("TCR: Sequence Lengths of {chain}"), subtitle = sample_name) +
      scale_fill_discrete(labels = labels_final)
    
    ggsave(glue("{outdir}/{sample_name}_clonalLength_{chain}.png"), width = 10, height = 6.5)
    
  }
  
  # Compare clones
  clonalCompare(combined.TCR_subset, 
                top.clones = 10, 
                samples = c(
                  names(combined.TCR_subset)
                ), 
                cloneCall="aa", 
                graph = "alluvial") + 
    scale_x_discrete(labels = labels_final) +
    labs(title = glue("TCR: Compare {sample_name}")) + 
    theme(
      legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(glue("{outdir}/{sample_name}_clonalCompare.png"), width = 13, height = 7.5)
  
  # Clonal scatter 
  # clonalScatter(combined.TCR_subset, 
  #               cloneCall ="gene", 
  #               x.axis = "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_Fol-18", 
  #               y.axis = "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_Fol-19",
  #               dot.size = "total",
  #               graph = "proportion")
  
  clonalHomeostasis(combined.TCR_subset, 
                    cloneCall = "gene") + 
    scale_x_discrete(labels = labels_final) + 
    labs(title = glue("TCR: Clonal homeostasis {sample_name}")) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) 
  ggsave(glue("{outdir}/{sample_name}_clonalHomeostasis.png"), width = 13, height = 7.5)
  
}
