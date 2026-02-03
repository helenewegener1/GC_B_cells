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
# seurat_obj_list <- readRDS("07_seurat_QC/out/seurat_obj_QC.rds")
# seurat_obj_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_singlets_list.rds")
seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")

# Investigate non-unique chains across contigs which should be filtered on umi. 

# ------------------------------------------------------------------------------
#                                      BCR
# ------------------------------------------------------------------------------

# Subset cells with BCR respectively. 
bcr_mask <- lapply(names(seurat_obj_list), function(x) {"bcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with BCR data 
bcr_seurat_obj_list <- seurat_obj_list[bcr_mask]

# Loading Data into scRepertoire
b_contigs.list <- list()
for (sample_name in names(bcr_seurat_obj_list)){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"
  
  # Load contig annotation file for sample 
  b_contigs <- read.csv(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/filtered_contig_annotations.csv"))
  
  # Load Seurat object
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # If sample is muliplexed, split contigs into list of contig files. 
  if ("ADT_classification" %in% colnames(seurat_obj@meta.data)){
    
    b_contigs <- createHTOContigList(b_contigs, 
                                     seurat_obj, 
                                     group.by = "ADT_maxID")
    
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

# combineBCR - one line per cell
# Combine using the default similarity clustering
combined.BCR <- combineBCR(b_contigs.list, 
                           samples = names(b_contigs.list), 
                           filterNonproductive = TRUE, # Default. Removes non-productive contigs , keeping only functional receptor chains. 
                           filterMulti = TRUE # Default. For cells with more than one heavy or light chain detected, this automatically selects the chain with the highest UMI count and discards the others. 
) 

# Adding Variables for Plotting: sample_high_level
combined.BCR <- addVariable(combined.BCR,
                            variable.name = "sample_high_level",
                            variables = str_split_i(names(combined.BCR), "_", 1))

# Adding Variables for Plotting: Inflamed / non-inflamed
patient_number <- str_split_i(names(b_contigs.list), "-", 1)
patient_number <- ifelse(patient_number == "HH117", "HH117_Crohns", "HH117_Control")

combined.BCR <- addVariable(combined.BCR,
                            variable.name = "patient",
                            variables = patient_number)

# The CTstrict column contains cluster IDs (e.g., "cluster.1")
head(combined.BCR[[1]][, c("barcode", "CTstrict", "IGH", "cdr3_aa1")])

# Quantifying Unique Clones
# clonalQuant() returns the total or relative numbers of unique clones
clonalQuant(
  combined.BCR, 
  cloneCall="strict",
  chain = "both", 
  scale = TRUE # relative percentage of unique clones scaled by the total repertoire size
) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 80, hjust = 1)) + 
  labs(title = "BCR: Relative numbers of unique clones", 
       x = "")

ggsave("20_VDJ/plot/BCR_scRepertoire/clonalQuant.png", width = 15, height = 10)

# Group by sample_high_level - overlook individual follicles.
# clonalQuant(
#   combined.BCR,
#   group.by = "sample_high_level",
#   cloneCall="strict",
#   chain = "both",
#   scale = TRUE # relative percentage of unique clones scaled by the total repertoire size
# ) +
#   labs(title = "BCR: Relative numbers of unique clones",
#        x = "")

# Split by samples since too many sample in one plot with each follicle. 
for (sample_name in names(bcr_seurat_obj_list)) {
  
  # sample_name <- "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH"
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  ############################## clonalAbundance ############################### 
  # total number of clones at specific frequencies within a sample or group
  # General: There's one clone that has the highest abundance within each sample. 
  
  combined.BCR_subset <- combined.BCR[startsWith(names(combined.BCR), sample_name)]
  
  if (length(names(combined.BCR_subset)) > 1){
    labels <- combined.BCR[startsWith(names(combined.BCR), sample_name)] %>% names() %>% str_split_i("_", 2)
  } else if (length(names(combined.BCR_subset)) == 1) {
    labels <- sample_name
  }
  
  # Not scaled
  clonalAbundance(combined.BCR_subset, 
                  cloneCall = "gene", 
                  scale = FALSE) + 
    labs(title = "BCR: Clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    scale_color_discrete(labels = labels) + 
    theme_dark()
  
  ggsave(glue("20_VDJ/plot/BCR_scRepertoire/{sample_name}_clonalAbundance.png"), width = 10, height = 6.5)
  
  # Scaled 
  clonalAbundance(combined.BCR_subset, 
                  cloneCall = "gene", 
                  scale = TRUE) + 
    labs(title = "BCR: Scaled clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    scale_color_discrete(labels = labels) + 
    theme_dark()
  
  ggsave(glue("20_VDJ/plot/BCR_scRepertoire/{sample_name}_clonalAbundance_scaled.png"), width = 10, height = 6.5)
  
  ############################## clonalAbundance ############################### 
  # Distribution of Sequence Lengths: looking at the length distribution of the CDR3 sequences
  for (chain in c("IGH", "IGL")){
    
    clonalLength(combined.BCR_subset, 
                 cloneCall="aa", 
                 chain = chain, 
                 scale = TRUE) + 
      labs(title = glue("BCR: Sequence Lengths of {chain}"), subtitle = sample_name) +
      scale_fill_discrete(labels = labels) 
    
    ggsave(glue("20_VDJ/plot/BCR_scRepertoire/{sample_name}_clonalLength_{chain}.png"), width = 10, height = 6.5)
    
  }
  
}

# Compare  
sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
group_A <- "Fol-1"
group_B <- "Fol-2"

clonalCompare(combined.BCR, 
              top.clones = 10, 
              samples = c(
                glue("{sample_name}_{group_A}"),
                glue("{sample_name}_{group_B}")
              ), 
              cloneCall="aa", 
              graph = "alluvial") + 
  scale_x_discrete(labels = c(group_A, group_B)) +
  labs(title = glue("BCR: Compare {sample_name}"),
       subtitle = glue("{group_A} vs {group_B}")) 
  # theme(legend.position = "none")

# All follicols 
fols <- names(combined.BCR)[names(combined.BCR) %>% startsWith(sample_name)] %>% str_split_i("_", 2)

clonalCompare(combined.BCR, 
              top.clones = 10, 
              samples = c(
                # glue("{sample_name}_{group_A}"), 
                # glue("{sample_name}_{group_B}"),
                # glue("{sample_name}_Fol-8")
                names(combined.BCR)[names(combined.BCR) %>% startsWith(sample_name)]
              ), 
              cloneCall="aa", 
              graph = "alluvial") + 
  # scale_x_discrete(labels = c(group_A, group_B)) + 
  scale_x_discrete(labels = fols) + 
  labs(title = glue("BCR: Compare {sample_name}")) + 
  theme(legend.position = "none")


