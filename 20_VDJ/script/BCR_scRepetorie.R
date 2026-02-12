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
library(gtools)

combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.rds")

# Quantifying Unique Clones
# clonalQuant() returns the total or relative numbers of unique clones
for (group in c("sample_high_level", "patient")){
  
  # group <- "sample_high_level"
  
  # Manual calculation example:
  # CTaa is a combination of the heavy and light chain --> paste(cdr3_aa1, cdr3_aa2, sep = "_")
  # CTstrict is the most stringent clonotype definition in scRepertoire. It combines both the V(D)J gene usage AND the nucleotide CDR3 sequences.
  # FIGURE OUT WHICH DEFINITION OF CLONOTYPES THERE ARE AND WHICH YOU WANT TO USE 
  # total_clones <- combined.BCR.filtered$`HH117-SILP-INF-PC`$CTstrict %>% length()
  # unique_clones <- combined.BCR.filtered$`HH117-SILP-INF-PC`$CTstrict %>% unique() %>% length() 
  # unique_clones/total_clones * 100
  # 
  # total_clones <- combined.BCR.filtered$`HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH`$CTstrict %>% length()
  # unique_clones <- combined.BCR.filtered$`HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH`$CTstrict %>% unique() %>% length() 
  # unique_clones/total_clones * 100
  
  clonalQuant(
    combined.BCR.filtered, 
    cloneCall="strict",
    chain = "both", 
    group.by = group, 
    scale = TRUE # relative percentage of unique clones scaled by the total repertoire size
  ) + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 80, hjust = 1)) + 
    labs(title = "BCR: Relative numbers of unique clones", 
         x = "")
  
  ggsave(glue("20_VDJ/plot/BCR_scRepertoire/clonalQuant_{group}.png"), width = 15, height = 10)
  
}

# Compare clones within each patient 
res <- lapply(
  c("HH117", "HH119"), function(x) {
    
    HH_mask <- grep(x, names(combined.BCR.filtered))
    combined.BCR.filtered_HH <- combined.BCR.filtered[HH_mask]
    # names(combined.BCR.filtered_HH119)
    
    clonalCompare(combined.BCR.filtered_HH, 
                  top.clones = 10, 
                  group.by = "sample_high_level",
                  cloneCall="aa", 
                  graph = "alluvial") + 
      labs(title = glue("BCR: Compare {x}")) + 
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(glue("20_VDJ/plot/BCR_scRepertoire/clonalCompare_{x}.png"), width = 13, height = 7.5)
    
    return(combined.BCR.filtered_HH)
    
  }
)

combined.BCR.filtered_HH117 <- res[[1]]
names(combined.BCR.filtered_HH117)

combined.BCR.filtered_HH119 <- res[[2]]
names(combined.BCR.filtered_HH119)

sample_names <- names(combined.BCR.filtered) %>% str_split_i("_", 1) %>% unique()

# Split by samples since too many sample in one plot with each follicle. 
for (sample_name in sample_names) {
  
  # sample_name <- "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH"
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  ############################## clonalAbundance ############################### 
  # total number of clones at specific frequencies within a sample or group
  # General: There's one clone that has the highest abundance within each sample. 
  
  combined.BCR_subset <- combined.BCR.filtered[startsWith(names(combined.BCR.filtered), sample_name)]
  
  outdir <- glue("20_VDJ/plot/BCR_scRepertoire/{sample_name}")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (length(names(combined.BCR_subset)) > 1){
    
    # Exclude Doublet and Negative subsets
    names_fols_mask <- grep("Doublet|Negative", names(combined.BCR_subset), invert = TRUE)
    combined.BCR_subset <- combined.BCR_subset[names_fols_mask]
    # names(combined.BCR_subset)

    # Sort naturally by the numeric part
    sorted_x <- mixedsort(names(combined.BCR_subset), decreasing = TRUE)
    combined.BCR_subset <- combined.BCR_subset[sorted_x]
    names(combined.BCR_subset)
    
    n_cells <- lapply(combined.BCR_subset, function(x) nrow(x))
    
    labels <- combined.BCR_subset[startsWith(names(combined.BCR_subset), sample_name)] %>% names() %>% str_split_i("_", 2)
    
    labels_final <- paste(labels, "-", n_cells, "cells")
    
  } else if (length(names(combined.BCR_subset)) == 1) {
    
    labels <- sample_name
    
    n_cells <- nrow(combined.BCR_subset)
    
    labels_final <- paste(labels, "-", n_cells, "cells")
  
  }
  
  # Not scaled
  clonalAbundance(combined.BCR_subset, 
                  cloneCall = "gene", 
                  scale = FALSE) + 
    labs(title = "BCR: Clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    # scale_color_discrete(labels = labels) +
    scale_color_discrete(labels = labels_final) +
    theme_dark()
  
  ggsave(glue("{outdir}/{sample_name}_clonalAbundance.png"), width = 10, height = 6.5)
  
  # Scaled 
  clonalAbundance(combined.BCR_subset, 
                  cloneCall = "gene", 
                  scale = TRUE) + 
    labs(title = "BCR: Scaled clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    # scale_color_discrete(labels = labels) + 
    scale_color_discrete(labels = labels_final) + 
    theme_dark()
  
  ggsave(glue("{outdir}/{sample_name}_clonalAbundance_scaled.png"), width = 10, height = 6.5)
  
  # Distribution of Sequence Lengths: looking at the length distribution of the CDR3 sequences
  for (chain in c("IGH", "IGL")){
    # chain <- "IGL"
    clonalLength(combined.BCR_subset, 
                 cloneCall="aa", 
                 chain = chain, 
                 scale = TRUE) + 
      labs(title = glue("BCR: Sequence Lengths of {chain}"), subtitle = sample_name) +
      scale_fill_discrete(labels = labels_final)
    
    ggsave(glue("{outdir}/{sample_name}_clonalLength_{chain}.png"), width = 10, height = 6.5)
    
  }
  
  # Compare clones
  clonalCompare(combined.BCR_subset, 
                top.clones = 10, 
                samples = c(
                  names(combined.BCR_subset)
                ), 
                cloneCall="aa", 
                graph = "alluvial") + 
    scale_x_discrete(labels = labels_final) +
    labs(title = glue("BCR: Compare {sample_name}")) + 
    theme(
      legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(glue("{outdir}/{sample_name}_clonalCompare.png"), width = 13, height = 7.5)
  
  # Clonal scatter 
  # clonalScatter(combined.BCR_subset, 
  #               cloneCall ="gene", 
  #               x.axis = "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_Fol-18", 
  #               y.axis = "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_Fol-19",
  #               dot.size = "total",
  #               graph = "proportion")
  
  clonalHomeostasis(combined.BCR_subset, 
                    cloneCall = "gene") + 
    scale_x_discrete(labels = labels_final) + 
    labs(title = glue("BCR: Clonal homeostasis {sample_name}")) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) 
  ggsave(glue("{outdir}/{sample_name}_clonalHomeostasis.png"), width = 13, height = 7.5)
  
}
