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

