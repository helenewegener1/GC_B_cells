library(alakazam)
library(scoper)
library(dplyr)
library(shazam)
library(dowser)
library(tidyverse)
library(glue)

# Newest versions
packageVersion("scoper")
packageVersion("alakazam")

# load the built-in human SHM targeting model
data(HH_S5F)  # or S5F, RS5NF depending on your preference

# Following this SCOPer tutorial: 
# https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

subset_nr <- "1"

subset <- readRDS(glue("46_sequence_driven_clustering/out/subset_{subset_nr}.rds"))

subset <- subset %>%
  mutate(
    sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
  )

# ------------------------------------------------------------------------------
# vj method: Groups clones based on junction sequences and SHM in V and J sequences
# ------------------------------------------------------------------------------

subset_clones_vj <- spectralClones(
  
  subset,
  method="vj",
  # threshold = list_thresholds[[HH]]$density,
  germline = "germline_alignment_d_mask",
  cell_id = "cell_id",
  junction = "junction",
  first = FALSE,
  targeting_model = HH_S5F
  
)


# ------------------------------------------------------------------------------
# Define top clones vj
# ------------------------------------------------------------------------------

top_clones <- subset_clones_vj %>%
  count(clone_id, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(clone_id)

# ------------------------------------------------------------------------------
# Wrangle data
# ------------------------------------------------------------------------------

subset_clones_vj <- subset_clones_vj %>% 
  mutate(
    v_call_subgroup = v_call %>% str_remove_all("\\*0[0-9]+"),
    v_call_subgroup = sapply(strsplit(v_call_subgroup, ","), function(x) unique(trimws(x))) %>%
      sapply(paste, collapse = ","),
    
    j_call_subgroup = j_call %>% str_remove_all("\\*0[0-9]+"),
    j_call_subgroup = sapply(strsplit(j_call_subgroup, ","), function(x) unique(trimws(x))) %>%
      sapply(paste, collapse = ",")
  ) %>% 
  group_by(clone_id) %>% 
  mutate(
    v_call_subgroup_clone = paste(unique(unlist(strsplit(v_call_subgroup, ","))), collapse = ","),
    j_call_subgroup_clone = paste(unique(unlist(strsplit(j_call_subgroup, ","))), collapse = ",")
  )

# Save 
saveRDS(subset_clones_vj, glue("46_sequence_driven_clustering/out/subset_{subset_nr}_clones_vj.rds"))


# ------------------------------------------------------------------------------
# Visualize top clones vj
# ------------------------------------------------------------------------------

source("10_broad_annotation/script/color_palette.R")

n_clones <- length(top_clones)

for (clone_nr in 1:length(top_clones)){
  
  # clone_nr <- 1
  clone <- top_clones[clone_nr]
  
  plot_df <- subset_clones_vj %>% 
    filter(clone_id == clone) %>% 
    mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
    add_count(sample_clean_fol, name = "Count")
  
  # Meta data
  n_cells <- plot_df %>% nrow()
  v_gene <- plot_df$v_call_subgroup_clone %>% unique()
  j_gene <- plot_df$j_call_subgroup_clone %>% unique()
  
  # Color by cell type
  plot_df %>%   
    ggplot(aes(y = sample_clean_fol, fill = L1_annotation)) +
    geom_bar() + 
    scale_fill_manual(values = L1_colors) + 
    geom_text(aes(x = Count, label = Count), 
              hjust = -0.2, color = "black",
              stat = "unique") +
    labs(
      title = glue("Subset {subset_nr}: Top {clone_nr} clone"),
      subtitle = glue("Clone ID: {clone}"),
      caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
      y = ""
    ) +
    theme_bw()
  
  ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}_clone_nr_{clone_nr}.png"), width = 15, height = 8.5)
  
}
