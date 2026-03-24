library(alakazam)
library(scoper)
library(dplyr)
library(shazam)

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

bcr_data <- readRDS("45_immcantation/out/rds/bcr_data_qc_annot.rds")

patients <- names(bcr_data)

bcr_data <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  bcr_data_HH <- bcr_data[[HH]]
  
  bcr_data_HH <- bcr_data_HH %>%
    mutate(
      sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
    )
  
  return(bcr_data_HH)
  
}) %>% setNames(patients)


# ------------------------------------------------------------------------------
# Identifying clones by sequence identity
# ------------------------------------------------------------------------------
# 
# # Clonal assignment using identical nucleotide sequences
# seq_clones <- lapply(patients, function(HH){
#   
#   identicalClones(
#     bcr_data[[HH]], 
#     method="nt", 
#     cell_id = "cell_id", 
#     junction = "junction", 
#     first = TRUE
#   )
#   
# }) %>% setNames(patients)
# 
# saveRDS(seq_clones, "45_immcantation/out/rds/seq_clones.rds")
# 
# # -------------------
# # HH119
# # -------------------
# 
# HH <- "HH119"
# 
# seq_clones[[HH]] %>% 
#   count(clone_id, sort = TRUE)
#   # count(clone_id, v_call, j_call, sort = TRUE)
#   
# # -------------------
# # HH117
# # -------------------
# 
# HH <- "HH117"
# 
# seq_clones[[HH]] %>% 
# count(clone_id, sort = TRUE)
# # count(clone_id, v_call, j_call, sort = TRUE)

# ------------------------------------------------------------------------------
# Identifying clones by spectral clustering
# novj method: Groups clones only based on junction sequences
# ------------------------------------------------------------------------------

spec_clones_novj <- lapply(patients, function(HH){
  
  spectralClones(
    bcr_data[[HH]], 
    method="novj", 
    junction="junction",# In the paper they said that using cdr3 instead of junction improved performance. 
    cell_id = "cell_id",
    first = TRUE # Taking only the first gene call per cell
  )
  
}) %>% setNames(patients)

saveRDS(spec_clones_novj, "45_immcantation/out/rds/spec_clones_novj.rds")

# -------------------
# HH119 - novj method
# -------------------

HH <- "HH119"
spec_clones_novj[[HH]] %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
plot(spec_clones_novj[[HH]], binwidth=0.02)

top_clone <- spec_clones_novj[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_novj[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_novj[[HH]] %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# -------------------
# HH117 - novj method
# -------------------

# HH117 has larger clones than with hierical clustering.

HH <- "HH117"
spec_clones_novj[[HH]] %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
plot(spec_clones_novj[[HH]], binwidth=0.02)

top_clone <- spec_clones_novj[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_novj[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_novj[[HH]] %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# -------------------
# Define top clones novj
# -------------------

top_GC_clones_novj <- lapply(patients, function(HH) {
  
  # find clones that contain at least 1 GC cell
  GC_clones <- spec_clones_novj[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    pull(clone_id) %>%
    unique()
  
  # rank those clones by total size (all cell types) and take top 10
  spec_clones_novj[[HH]] %>%
    filter(clone_id %in% GC_clones) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_id)
  
}) %>% setNames(patients)

# -------------------
# Look at top clones novj
# -------------------

lapply(patients, function(HH){
  
  # HH <- "HH119"
  # HH <- "HH117"
  
  spec_clones_novj[[HH]] %>% 
    filter(clone_id %in% top_GC_clones_novj[[HH]]) %>% 
    mutate(
      v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
    ) %>% 
    mutate(
      j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
    ) %>% 
    count(clone_id, v_gene, j_gene, sort = TRUE) 
  
}) %>% setNames(patients)

# -------------------
# Visualize top clones novj
# -------------------

source("10_broad_annotation/script/color_palette.R")

for (HH in patients){
  
  # HH <- "HH117"
  HH_top_clones <- top_GC_clones_novj[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    plot_df <- spec_clones_novj[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    j_gene <- plot_df$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    
    plot_df %>%   
      ggplot(aes(y = sample_clean_fol, fill = celltype_broad)) +
      geom_bar() + 
      scale_fill_manual(values = celltype_colors) + 
      geom_text(aes(x = Count, label = Count), 
                hjust = -0.2, color = "black",
                stat = "unique") +
      labs(
        title = glue("{HH}: Top {clone_nr} GCB clone spec_novj"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
        y = ""
      ) +
      theme_bw()
    
    ggsave(glue("45_immcantation/plot/GC_clones_spec_novj/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5)
    
  }
}

# ------------------------------------------------------------------------------
# vj method: Groups clones based on junction sequences and SHM in V and J sequences
# ------------------------------------------------------------------------------

spec_clones_vj <- lapply(patients, function(HH){
  
  spectralClones(
    bcr_data[[HH]], 
    method="vj", 
    cell_id = "cell_id", 
    junction = "junction",
    first = TRUE
  )
  
}) %>% setNames(patients)

saveRDS(spec_clones_vj, "45_immcantation/out/rds/spec_clones_vj.rds")

# -------------------
# HH119 - vj method
# -------------------

HH <- "HH119"
spec_clones_vj[[HH]] %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
# plot(spec_clones_vj[[HH]], binwidth=0.02)

top_clone <- spec_clones_vj[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_vj[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_vj[[HH]] %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# -------------------
# HH117 - vj method
# -------------------

# HH117 has larger clones than with hierical clustering.

HH <- "HH117"
spec_clones_vj[[HH]] %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
# plot(spec_clones_vj[[HH]], binwidth=0.02)

top_clone <- spec_clones_vj[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_vj[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_vj[[HH]] %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# -------------------
# Define top clones vj
# -------------------

top_GC_clones_vj <- lapply(patients, function(HH) {
  
  # find clones that contain at least 1 GC cell
  GC_clones <- spec_clones_vj[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    pull(clone_id) %>%
    unique()
  
  # rank those clones by total size (all cell types) and take top 10
  spec_clones_vj[[HH]] %>%
    filter(clone_id %in% GC_clones) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_id)
  
}) %>% setNames(patients)

# -------------------
# Look at top clones vj
# -------------------

lapply(patients, function(HH){

  # HH <- "HH119"
  # HH <- "HH117"
  
  spec_clones_vj[[HH]] %>% 
    filter(clone_id %in% top_GC_clones_vj[[HH]]) %>% 
    mutate(
      v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
    ) %>% 
    mutate(
      j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
    ) %>% 
    count(clone_id, v_gene, j_gene, sort = TRUE) 
  
}) %>% setNames(patients)


# -------------------
# Visualize top clones vj
# -------------------

source("10_broad_annotation/script/color_palette.R")

for (HH in patients){
  
  # HH <- "HH117"
  HH_top_clones <- top_GC_clones_vj[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    plot_df <- spec_clones_vj[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    j_gene <- plot_df$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    
    plot_df %>%   
      ggplot(aes(y = sample_clean_fol, fill = celltype_broad)) +
      geom_bar() + 
      scale_fill_manual(values = celltype_colors) + 
      geom_text(aes(x = Count, label = Count), 
                hjust = -0.2, color = "black",
                stat = "unique") +
      labs(
        title = glue("{HH}: Top {clone_nr} GCB clone spec_vj"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
        y = ""
      ) +
      theme_bw()
    
    ggsave(glue("45_immcantation/plot/GC_clones_spec_vj/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5)
    
  }
}


# ------------------------------------------------------------------------------
# Compare methods: novj VS vj 
# ------------------------------------------------------------------------------

# -------------------
# HH117 
# -------------------

HH <- "HH119"

HH_spec_clones_vj <- spec_clones_vj[[HH]]
HH_spec_clones_novj <- spec_clones_novj[[HH]]

# Not in the same order...
table(HH_spec_clones_vj$cell_id %in% HH_spec_clones_novj$cell_id) 
table(HH_spec_clones_vj$cell_id == HH_spec_clones_novj$cell_id) 
length(HH_spec_clones_vj$cell_id)
length(HH_spec_clones_novj$cell_id)

# reorder novj to match the order of vj
HH_spec_clones_novj <- HH_spec_clones_novj[match(HH_spec_clones_vj$cell_id_seurat, HH_spec_clones_novj$cell_id_seurat), ]

# verify they now match
table(HH_spec_clones_vj$cell_id_seurat == HH_spec_clones_novj$cell_id_seurat)

# Investigate clone overlap clone_ids

# get top N clones from each method
top_n <- 50

top_vj <- HH_spec_clones_vj %>% count(clone_id, sort = TRUE) %>% slice_head(n = top_n) %>% pull(clone_id)
top_novj <- HH_spec_clones_novj %>% count(clone_id, sort = TRUE) %>% slice_head(n = top_n) %>% pull(clone_id)

# build overlap table for top clones only
overlap <- HH_spec_clones_vj %>%
  select(cell_id_seurat, clone_vj = clone_id) %>%
  left_join(HH_spec_clones_novj %>% select(cell_id_seurat, clone_novj = clone_id),
            by = "cell_id_seurat") %>%
  filter(clone_vj %in% top_vj, clone_novj %in% top_novj) %>%
  count(clone_vj, clone_novj)

# heatmap
ggplot(overlap, aes(x = factor(clone_novj), y = factor(clone_vj), fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "noVJ clone", y = "VJ clone", fill = "# cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
