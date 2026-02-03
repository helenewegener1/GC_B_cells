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
# seurat_obj_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_doublets_list.rds")
seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")

# Investigate non-unique chains across contigs which should be filtered on umi. 

# Subset cells with TCR respectively. 
# Which samples have TCR data 
tcr_mask <- lapply(names(seurat_obj_list), function(x) {"tcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with TCR data 
tcr_seurat_obj_list <- seurat_obj_list[tcr_mask]

################################# scRepertorie - TCR ################################# 

# Loading Data into scRepertoire
t_contigs.list <- list()
for (sample_name in names(tcr_seurat_obj_list)){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  # sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"
  
  # Load contig annotation file for sample 
  t_contigs <- read.csv(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_t/filtered_contig_annotations.csv"))
  
  # Load Seurat object
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # If sample is muliplexed, split contigs into list of contig files. 
  if ("ADT_classification" %in% colnames(seurat_obj@meta.data)){
    
    t_contigs <- createHTOContigList(t_contigs, 
                                     seurat_obj, 
                                     group.by = "ADT_maxID")
    
    # Have name include sample_name
    names(t_contigs) <- paste(sample_name, names(t_contigs), sep = "_")
    
    # Merge list with t_contigs.list
    t_contigs.list <- c(t_contigs.list, t_contigs)
    
  } else {
    
    # Append contiguous annotation of non-multiplexed sample to t_contigs.list
    t_contigs.list[[sample_name]] <- t_contigs
    
  }
  
}

# Check names of contig list 
names(t_contigs.list) 

# Clonal identification: Group cells by shared CDR3 sequence using combineTCR()
combined.TCR <- combineTCR(t_contigs.list,
                           # samples = names(t_contigs.list), # does not work for some reason 
                           removeNA = FALSE,
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

names(combined.TCR) <- names(t_contigs.list) 

# Output is one cell per line. 
head(combined.TCR[[1]])
names(combined.TCR)
colnames(combined.TCR[[1]])

# Adding Variables for Plotting: sample_high_level
combined.TCR <- addVariable(combined.TCR,
                            variable.name = "sample_high_level",
                            variables = str_split_i(names(combined.TCR), "_", 1))

# Quantifying Unique Clones
# clonalQuant() returns the total or relative numbers of unique clones
clonalQuant(
  combined.TCR, 
  cloneCall="strict", 
  chain = "both", 
  scale = TRUE # relative percentage of unique clones scaled by the total repertoire size
) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 80, hjust = 1)) + 
  labs(title = "TCR: Relative numbers of unique clones", 
       x = "")

ggsave("20_VDJ/plot/TCR_scRepertoire/clonalQuant.png", width = 15, height = 10)

# Group by sample_high_level - overlook individual follicles. 
# clonalQuant(
#   combined.TCR,
#   group.by = "sample_high_level",
#   cloneCall="strict",
#   chain = "both",
#   scale = TRUE # relative percentage of unique clones scaled by the total repertoire size
# ) +
#   labs(title = "TCR: Relative numbers of unique clones",
#        x = "")

# Combined plots 

clonalAbundance(combined.TCR, 
                cloneCall = "gene", 
                scale = FALSE)
ggsave("20_VDJ/plot/TCR_scRepertoire/clonalAbundance.png", width = 10, height = 6.5)

# Split by samples since too many sample in one plot with each follicle. 
for (sample_name in names(tcr_seurat_obj_list)) {
  
  # sample_name <- "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH"
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  ############################## clonalAbundance ############################### 
  # total number of clones at specific frequencies within a sample or group
  # General: There's one clone that has the highest abundance within each sample. 
  
  combined.TCR_subset <- combined.TCR[startsWith(names(combined.TCR), sample_name)]
  
  if (length(names(combined.TCR_subset)) > 1){
    labels <- combined.TCR[startsWith(names(combined.TCR), sample_name)] %>% names() %>% str_split_i("_", 2)
  } else if (length(names(combined.TCR_subset)) == 1) {
    labels <- sample_name
  }

  # Not scaled
  clonalAbundance(combined.TCR_subset, 
                  cloneCall = "gene", 
                  scale = FALSE) + 
    labs(title = "TCR: Clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    scale_color_discrete(labels = labels) + 
    theme_dark()
  
  ggsave(glue("20_VDJ/plot/TCR_scRepertoire/clonalAbundance_{sample_name}.png"), width = 10, height = 6.5)
  
  # Scaled 
  clonalAbundance(combined.TCR_subset, 
                  cloneCall = "gene", 
                  scale = TRUE) + 
    labs(title = "TCR: Scaled clonal abundance - The total number of clones at specific frequencies",
         subtitle = sample_name) + 
    scale_color_discrete(labels = labels) + 
    theme_dark()
  
  ggsave(glue("20_VDJ/plot/TCR_scRepertoire/clonalAbundance_scaled_{sample_name}.png"), width = 10, height = 6.5)
  
  ############################## clonalAbundance ############################### 
  # Distribution of Sequence Lengths: looking at the length distribution of the CDR3 sequences
  clonalLength(combined.TCR_subset, 
               cloneCall="aa", 
               chain = "TRA", 
               scale = TRUE) + 
    labs(title = "TCR: Sequence Lengths", subtitle = sample_name) +
    scale_fill_discrete(labels = labels) 
  
  ggsave(glue("20_VDJ/plot/TCR_scRepertoire/clonalLength_scaled_{sample_name}.png"), width = 10, height = 6.5)
  
}

# Compare samples
group_A <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
group_B <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"

clonalCompare(combined.TCR, 
              top.clones = 10, 
              samples = c(
                group_A,
                group_B
              ), 
              cloneCall="aa", 
              graph = "alluvial") + 
  scale_x_discrete(labels = c(group_A, group_B)) +
  labs(title = glue("TCR: Compare {sample_name}"),
       subtitle = glue("{group_A} vs {group_B}"))

group_A <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"
group_B <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"

clonalCompare(combined.TCR, 
              top.clones = 10, 
              samples = c(
                group_A,
                group_B
              ), 
              cloneCall="aa", 
              graph = "alluvial") + 
  scale_x_discrete(labels = c(group_A, group_B)) +
  labs(title = glue("TCR: Compare {sample_name}"),
       subtitle = glue("{group_A} vs {group_B}"))

# Compare  
sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
group_A <- "Fol-1"
group_B <- "Fol-5"

clonalCompare(combined.TCR, 
              top.clones = 10, 
              samples = c(
                glue("{sample_name}_{group_A}"),
                glue("{sample_name}_{group_B}")
              ), 
              cloneCall="aa", 
              graph = "alluvial") + 
  scale_x_discrete(labels = c(group_A, group_B)) +
  labs(title = glue("TCR: Compare {sample_name}"),
       subtitle = glue("{group_A} vs {group_B}")) + 
  theme(legend.position = "none")

# All follicles 
fols <- names(combined.TCR)[names(combined.TCR) %>% startsWith(sample_name)] %>% str_split_i("_", 2)

clonalCompare(combined.TCR, 
              top.clones = 10, 
              samples = c(
                # glue("{sample_name}_{group_A}"), 
                # glue("{sample_name}_{group_B}"),
                # glue("{sample_name}_Fol-8")
                names(combined.TCR)[names(combined.TCR) %>% startsWith(sample_name)]
              ), 
              cloneCall="aa", 
              graph = "alluvial") + 
  # scale_x_discrete(labels = c(group_A, group_B)) + 
  scale_x_discrete(labels = fols) + 
  labs(title = glue("TCR: Compare {sample_name}")) + 
  theme(legend.position = "none")

