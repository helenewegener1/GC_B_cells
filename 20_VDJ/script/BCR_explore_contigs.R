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

# Load data
seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds") # Broad annotation
seurat_obj_ADT_demulti <- readRDS("10_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds") # Follicol information 

# Investigate non-unique chains across contigs which should be filtered on umi. 

# Subset cells with BCR respectively. 
bcr_mask <- lapply(names(seurat_obj_ADT_demulti), function(x) {"bcr_v_gene_contig_1" %in% colnames(seurat_obj_ADT_demulti[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with BCR data 
bcr_seurat_obj_list <- seurat_obj_ADT_demulti[bcr_mask]

# Loading Data into scRepertoire
b_contigs.list <- list()
for (sample_name in names(bcr_seurat_obj_list)){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  # Load contig annotation file for sample 
  b_contigs <- read.csv(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/filtered_contig_annotations.csv"))
  
  # Load Seurat object
  seurat_obj <- seurat_obj_ADT_demulti[[sample_name]]
  
  # If sample is muliplexed, split contigs into list of contig files. 
  if ("manual_ADT_ID" %in% colnames(seurat_obj@meta.data)){
    
    b_contigs <- createHTOContigList(b_contigs, 
                                     seurat_obj, 
                                     group.by = "manual_ADT_ID")
    
    # Have name include sample_name
    names(b_contigs) <- paste(sample_name, names(b_contigs), sep = "_")
    
    # Merge list with b_contigs.list
    b_contigs.list <- c(b_contigs.list, b_contigs)
    
  } else {
    
    # Append contiguous annotation of non-multiplexed sample to b_contigs.list
    b_contigs.list[[sample_name]] <- b_contigs
    
  }
  
}

# Check names of contig list 
names(b_contigs.list)

# Here we have one row per contig
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 536 - N contigs
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$contig_id %>% str_split_i("_", 1) %>% unique() %>% length() # 295 - N cells
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3 %>% unique() %>% length() # 486 - amino acid sequence of CDR3
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt %>% unique() %>% length() # 500 - nucleotide sequence of CDR3 (More unique because silent mutations)

# ------------------------------------------------------------------------------
# Understanding start
# ------------------------------------------------------------------------------

# df <- b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`
# 
# # Count chains per barcode
# chain_counts <- df %>%
#   group_by(barcode, chain) %>%
#   summarise(n_chains = n(), .groups = "drop") %>%
#   pivot_wider(names_from = chain, 
#               values_from = n_chains, 
#               values_fill = 0)
# 
# # If you have IGH, IGK, IGL columns, create a summary
# # Adjust column names based on what you actually have
# chain_summary <- chain_counts %>%
#   mutate(
#     total_heavy = if("IGH" %in% names(.)) IGH else 0,
#     total_light = rowSums(select(., matches("IGK|IGL")), na.rm = TRUE),
#     chain_config = paste0("IGH:", total_heavy, " Light:", total_light)
#   ) %>% 
#   count(chain_config) %>%
#   arrange(desc(n))
# 
# ggplot(chain_summary, aes(x = reorder(chain_config, n), y = n)) +
#   geom_col(fill = "steelblue") +
#   geom_text(aes(label = n), hjust = -0.2) +
#   coord_flip() +
#   labs(
#     title = "Distribution of Chain Configurations",
#     x = "Chain Configuration",
#     y = "Number of Cells"
#   ) +
#   theme_minimal(base_size = 10) +
#   theme(plot.title = element_text(face = "bold"))


# What combineBCR does
## Keeps all 295 cells (even those with only heavy or only light) - from plot --> 208 + 62 + 10= 280
## For the 5 cells with 2+ heavy chains: selects one (highest UMI if filterMulti=TRUE)
## For the 13 cells with 2+ light chains: selects one (highest UMI if filterMulti=TRUE)
## Assigns heavy to cdr3_aa1/nt1 and light to cdr3_aa2/nt2
## Cells missing a chain will have NA in the corresponding columns

# ------------------------------------------------------------------------------
# Understanding end
# ------------------------------------------------------------------------------

# combineBCR - one line per cell
# Combine using the default similarity clustering
combined.BCR.NOTfiltered <- combineBCR(b_contigs.list,
                                       samples = names(b_contigs.list), 
                                       # filterNonproductive = FALSE, 
                                       filterMulti = FALSE 
) 

combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 295 - N cells 
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa1 %>% unique() %>% length() # 218 - heavy chain amino acid sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt1 %>% unique() %>% length() # 219 - heavy chain nucleotide sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa2 %>% unique() %>% length() # 256 - ligth chain amino acid sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt2 %>% unique() %>% length() # 267 - ligth chain nucleotide sequence

# How many receptors with NA in any chains 
is.na(combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGH) %>% table() # 62
is.na(combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGLC) %>% table() # 10 

is.na(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGH) %>% table() # 0 
is.na(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGLC) %>% table() # 0 

head(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`, n = 5)

# Combine using the default similarity clustering
combined.BCR.filtered <- combineBCR(b_contigs.list,
                                    samples = names(b_contigs.list), 
                                    removeNA = TRUE,
                                    threshold = 0.85, # Default is 0.85. Oliver used default. 
                                    # filterNonproductive = TRUE, # Default. Removes non-productive contigs , keeping only functional receptor chains. 
                                    filterMulti = TRUE # Default. For cells with more than one heavy or light chain detected, this automatically selects the chain with the highest UMI count and discards the others. 
) 

combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 223 - N cells 
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa1 %>% unique() %>% length() # 209 - heavy chain amino acid sequence
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt1 %>% unique() %>% length() # 210 - heavy chain nucleotide sequence
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa2 %>% unique() %>% length() # 200 - ligth chain amino acid sequence
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt2 %>% unique() %>% length() # 208 - ligth chain nucleotide sequence

# Adding Variables for Plotting: sample_high_level
combined.BCR.filtered <- addVariable(combined.BCR.filtered,
                                     variable.name = "sample_high_level",
                                     variables = str_split_i(names(combined.BCR.filtered), "_", 1))

# Adding Variables for Plotting: Inflamed / non-inflamed
patient_number <- str_split_i(names(b_contigs.list), "-", 1)
patient_number <- ifelse(patient_number == "HH117", "HH117_Crohns", "HH117_Control")

combined.BCR.filtered <- addVariable(combined.BCR.filtered,
                                     variable.name = "patient",
                                     variables = patient_number)

# The CTstrict column contains cluster IDs (e.g., "cluster.1")
head(combined.BCR.filtered[[1]][, c("barcode", "CTstrict", "IGH", "cdr3_aa1")])

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

# Split by samples since too many sample in one plot with each follicle. 
for (sample_name in names(bcr_seurat_obj_list)) {
  
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
  
}
