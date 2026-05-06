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

list_thresholds <- readRDS("45_immcantation/out/rds/04_list_thresholds.rds")
bcr_data <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")

bcr_data$HH117$j_call %>% table()

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
#     first = FALSE
#   )
# 
# }) %>% setNames(patients)
# 
# saveRDS(seq_clones, "45_immcantation/out/rds/seq_clones.rds")
# 
# 
# 
# # ------------------------------------------------------------------------------
# # novj method: Groups clones based on junction sequences and SHM in V and J sequences
# # ------------------------------------------------------------------------------
# 
# spec_clones_novj_no_threshold <- lapply(patients, function(HH){
#   
#   # HH <- "HH119"
#   spectralClones(
#     bcr_data[[HH]],
#     method="novj",
#     # threshold = list_thresholds[[HH]]$density,
#     germline = "germline_alignment_d_mask",
#     cell_id = "cell_id",
#     junction = "junction",
#     first = FALSE,
#     targeting_model = HH_S5F, 
#     summarize_clones = TRUE
#   )
#   
# }) %>% setNames(patients)
# 
# saveRDS(spec_clones_novj_no_threshold, "45_immcantation/out/rds/05_spec_clones_novj_no_threshold.rds")
# 
# spec_clones_novj_density_threshold <- lapply(patients, function(HH){
#   
#   # HH <- "HH119"
#   spectralClones(
#     bcr_data[[HH]],
#     method="novj",
#     threshold = list_thresholds[[HH]]$density,
#     germline = "germline_alignment_d_mask",
#     cell_id = "cell_id",
#     junction = "junction",
#     first = FALSE,
#     targeting_model = HH_S5F, 
#     summarize_clones = TRUE
#   )
#   
# }) %>% setNames(patients)
# 
# saveRDS(spec_clones_novj_density_threshold, "45_immcantation/out/rds/05_spec_clones_novj_density_threshold.rds")
# 
# spec_clones_novj_gmm_threshold <- lapply(patients, function(HH){
#   
#   # HH <- "HH119"
#   spectralClones(
#     bcr_data[[HH]],
#     method="novj",
#     threshold = list_thresholds[[HH]]$gmm,
#     germline = "germline_alignment_d_mask",
#     cell_id = "cell_id",
#     junction = "junction",
#     first = FALSE,
#     targeting_model = HH_S5F, 
#     summarize_clones = TRUE
#   )
#   
# }) %>% setNames(patients)
# 
# saveRDS(spec_clones_novj_gmm_threshold, "45_immcantation/out/rds/05_spec_clones_novj_gmm_threshold.rds")
# 
# # -------------------
# # Plot distributions
# # -------------------
# 
# 
# lapply(patients, function(HH){
#   
#   # HH <- "HH117"
#   
#   # No threshold
#   res <- spec_clones_vj_no_threshold[[HH]]
#   # Plot a histogram of inter and intra clonal distances
#   plot(res, binwidth=0.02) 
#   ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_scoper_vj_no_threshold.png"))
#   
#   # density threshold
#   res <- spec_clones_vj_density_threshold[[HH]]
#   # Plot a histogram of inter and intra clonal distances
#   plot(res, binwidth=0.02) + labs(subtitle = glue("Density threshold by findThreshold: {list_thresholds[[HH]]$density}"))
#   ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_scoper_vj_density_threshold.png"))
#   
#   # GMM threshold
#   res <- spec_clones_vj_gmm_threshold[[HH]]
#   # Plot a histogram of inter and intra clonal distances
#   plot(res, binwidth=0.02) + labs(subtitle = glue("GMM threshold by findThreshold: {list_thresholds[[HH]]$gmm}"))
#   ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_scoper_vj_gmm_threshold.png"))
#   
#   
# })

# ------------------------------------------------------------------------------
# vj method: Groups clones based on junction sequences and SHM in V and J sequences
# ------------------------------------------------------------------------------

spec_clones_vj_no_threshold <- lapply(patients, function(HH){

  # HH <- "HH119"
  spectralClones(
    bcr_data[[HH]],
    method="vj",
    # threshold = list_thresholds[[HH]]$density,
    germline = "germline_alignment_d_mask",
    cell_id = "cell_id",
    junction = "junction",
    first = FALSE,
    targeting_model = HH_S5F, 
    summarize_clones = TRUE
  )

}) %>% setNames(patients)

saveRDS(spec_clones_vj_no_threshold, "45_immcantation/out/rds/05_spec_clones_vj_no_threshold.rds")

spec_clones_vj_density_threshold <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  spectralClones(
    bcr_data[[HH]],
    method="vj",
    threshold = list_thresholds[[HH]]$density,
    germline = "germline_alignment_d_mask",
    cell_id = "cell_id",
    junction = "junction",
    first = FALSE,
    targeting_model = HH_S5F, 
    summarize_clones = TRUE
  )
  
}) %>% setNames(patients)

saveRDS(spec_clones_vj_density_threshold, "45_immcantation/out/rds/05_spec_clones_vj_density_threshold.rds")

spec_clones_vj_gmm_threshold <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  spectralClones(
    bcr_data[[HH]],
    method="vj",
    threshold = list_thresholds[[HH]]$gmm,
    germline = "germline_alignment_d_mask",
    cell_id = "cell_id",
    junction = "junction",
    first = FALSE,
    targeting_model = HH_S5F, 
    summarize_clones = TRUE
  )
  
}) %>% setNames(patients)

saveRDS(spec_clones_vj_gmm_threshold, "45_immcantation/out/rds/05_spec_clones_vj_gmm_threshold.rds")

# -------------------
# Plot distributions
# -------------------


spec_clones_vj_no_threshold <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_no_threshold.rds")
spec_clones_vj_density_threshold <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_density_threshold.rds")
spec_clones_vj_gmm_threshold <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_gmm_threshold.rds")



lapply(patients, function(HH){
  
  # HH <- "HH117"
  
  # No threshold
  res <- spec_clones_vj_no_threshold[[HH]]
  # Plot a histogram of inter and intra clonal distances
  plot(res, binwidth=0.02) 
  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_scoper_vj_no_threshold.png"))
  
  # density threshold
  res <- spec_clones_vj_density_threshold[[HH]]
  # Plot a histogram of inter and intra clonal distances
  plot(res, binwidth=0.02) + labs(subtitle = glue("Density threshold by findThreshold: {list_thresholds[[HH]]$density}"))
  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_scoper_vj_density_threshold.png"))
  
  # GMM threshold
  res <- spec_clones_vj_gmm_threshold[[HH]]
  # Plot a histogram of inter and intra clonal distances
  plot(res, binwidth=0.02) + labs(subtitle = glue("GMM threshold by findThreshold: {list_thresholds[[HH]]$gmm}"))
  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_scoper_vj_gmm_threshold.png"))

  
})

# -------------------
# Plot distributions
# -------------------

spec_clones_vj$HH117@db

# -------------------
# Detect main V and J gene for each clone
# -------------------

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

spec_clones_vj <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  spec_clones_vj[[HH]] %>%
    group_by(clone_id) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)

saveRDS(spec_clones_vj, "45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")

# -------------------
# Save clone_ids per sample
# -------------------

# Load data
spec_clones_vj <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")

patients <- names(spec_clones_vj)

# CLONES ARE SPLIT ON SAMPLE LEVEL AND NOT PATIENT LEVEL
# CLONES WILL ONLY BE PER SAMPLE AND NOT ACROSS - so fix this if we choose to go with this. 
# ligth_cluster.py removes cells without light chains which might be too strict.
# resolveLightChains is softer
# for (HH in patients) {
# 
#   # HH <- "HH119"
#   # get the cloned data for this patient
#   spec_clones_vj_HH <- spec_clones_vj[[HH]]
# 
#   # split back by sample and write one file per sample
#   for (s in unique(spec_clones_vj_HH$sample_id)) {
# 
#     # s <- "HH119-SILP-PC"
# 
#     sample_db <- spec_clones_vj_HH %>% filter(sample_id == s)
#     sample_db$cell_id <- sample_db$cell_id %>% str_split_i("_", 2)
# 
#     outfile <- glue("45_immcantation/out/{s}/{s}_heavy_germ-pass_clone-pass.tsv")
# 
#     airr::write_rearrangement(sample_db, outfile)
# 
#   }
# }

# heavy <- read.delim("45_immcantation/out/HH117-SILP-INF-PC/HH117-SILP-INF-PC_heavy_germ-pass_clone-pass.tsv")
# light <- read.delim("45_immcantation/out/HH117-SILP-INF-PC/HH117-SILP-INF-PC_light_parse-select.tsv")
# table(light$cell_id %in% heavy$cell_id)

# clone_10x <- read.delim("45_immcantation/out/HH117-SILP-INF-PC/HH117-SILP-INF-PC_10X_clone-pass.tsv")
# nrow(clone_10x)
# spec_clones_vj[["HH117"]] %>% filter(sample_id == "HH117-SILP-INF-PC")

# -------------------
# Define top clones vj
# -------------------

top_GC_clones_vj <- lapply(patients, function(HH) {
  
  # # find clones that contain at least 1 GC cell
  # GC_clones <- spec_clones_vj[[HH]] %>%
  #   filter(celltype_broad == "GC_B_cells") %>%
  #   pull(clone_id) %>%
  #   unique()
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- spec_clones_vj[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    group_by(clone_id) %>%
    summarise(n_samples = n_distinct(sample_clean_fol)) %>%
    filter(n_samples >= 2) %>%
    pull(clone_id)
  
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
    # count(clone_id, v_gene, j_gene, sort = TRUE)
    count(clone_id, v_call_majority, j_call_majority, sort = TRUE) 
  
}) %>% setNames(patients)

# -------------------
# Visualize top clones vj
# -------------------

source("10_broad_annotation/script/color_palette.R")

for (HH in patients){
  
  # HH <- "HH119"
  HH_top_clones <- top_GC_clones_vj[[HH]]
  
  n_clones <- length(HH_top_clones)
  
  # Heatmap of sample_clean_fol and c_call (Isotype)
  spec_clones_vj[[HH]] %>% 
    filter(clone_id %in% HH_top_clones) %>% 
    mutate(c_call = ifelse(c_call == "", NA, c_call)) %>%
    group_by(clone_id, c_call) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(clone_id, c_call, fill = list(Count = 0)) %>% 
    group_by(clone_id) %>%
    mutate(Pct = Count / sum(Count) * 100) %>%
    ungroup() %>%
    mutate(clone_id = factor(clone_id, HH_top_clones)) %>% 
    ggplot(aes(x = clone_id, y = c_call, fill = Pct)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(Pct > 0, sprintf("%.1f%%", Pct), "")),
              color = "white", size = 2.5) +
    scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                        limits = c(0, 100)) +
    labs(
      x = "Clone ID",
      y = "Ig Class",
      fill = "% per clone",
      title = glue("{HH}: Top {n_clones} GCB clone spec_vj")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "none"
    )
  
  ggsave(glue("45_immcantation/plot/GC_clones_spec_vj/{HH}_heatmap_across_clones_and_isotypes.png"), width = 6, height = 4.5, dpi = 1000)
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    plot_df <- spec_clones_vj[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call_majority %>% unique()
    j_gene <- plot_df$j_call_majority %>% unique()
    
    # Color by cell type
    plot_df %>%   
      ggplot(aes(y = sample_clean_fol, fill = L1_annotation)) +
      geom_bar() + 
      scale_fill_manual(values = L1_colors) + 
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
    
    # Heatmap of sample_clean_fol and c_call (Isotype) and celltype
    plot_df %>%
      mutate(c_call = ifelse(c_call == "", NA, c_call)) %>%
      group_by(sample_clean_fol, c_call, L1_annotation) %>%
      summarise(Count = n(), .groups = "drop") %>%
      complete(sample_clean_fol, c_call, L1_annotation, 
               fill = list(Count = 0)) %>%
      ggplot(aes(x = sample_clean_fol, y = c_call, fill = Count)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(Count > 0, Count, "")),
                color = "white", size = 2.5) +
      scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                          limits = c(0, NA)) +
      facet_wrap(~ L1_annotation, nrow = 2) +
      labs(
        x = NULL,
        y = "Ig Class",
        fill = "Count",
        title = glue("{HH}: Top {clone_nr} GCB clone spec_vj"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "right"
      )
    
    n_samples <- plot_df$sample_clean_fol %>% unique() %>% length()
    
    if (n_samples > 8){
      ggsave(glue("45_immcantation/plot/GC_clones_spec_vj/{HH}_clone_nr_{clone_nr}_heatmap_across_samples_and_isotypes_and_celltype.png"), width = 20, height = 8, dpi = 1000)
    } else {
      ggsave(glue("45_immcantation/plot/GC_clones_spec_vj/{HH}_clone_nr_{clone_nr}_heatmap_across_samples_and_isotypes_and_celltype.png"), width = 8, height = 6, dpi = 1000)
    }
    
  }
}

# ------------------------------------------------------------------------------
# Clones (Not GC specific)
# ------------------------------------------------------------------------------

# -------------------
# Define top clones vj
# -------------------

top_clones_vj <- lapply(patients, function(HH) {
  
  # find clones that have GC cells in at least 2 different sample_ids
  clones <- spec_clones_vj[[HH]] %>%
    count(clone_id, sort = TRUE) %>%
    pull(clone_id)
  
  # rank those clones by total size (all cell types) and take top 10
  spec_clones_vj[[HH]] %>%
    filter(clone_id %in% clones) %>%
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
    filter(clone_id %in% top_clones_vj[[HH]]) %>% 
    # count(clone_id, v_gene, j_gene, sort = TRUE)
    count(clone_id, v_call_majority, j_call_majority, sort = TRUE) 
  
}) %>% setNames(patients)

# -------------------
# Visualize top clones vj
# -------------------

source("10_broad_annotation/script/color_palette.R")

for (HH in patients){
  
  # HH <- "HH117"
  HH_top_clones <- top_clones_vj[[HH]]
  
  n_clones <- length(HH_top_clones)
  
  # Heatmap of sample_clean_fol and c_call (Isotype)
  spec_clones_vj[[HH]] %>% 
    filter(clone_id %in% HH_top_clones) %>% 
    mutate(c_call = ifelse(c_call == "", NA, c_call)) %>%
    group_by(clone_id, c_call) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(clone_id, c_call, fill = list(Count = 0)) %>% 
    mutate(clone_id = factor(clone_id, HH_top_clones)) %>% 
    ggplot(aes(x = clone_id, y = c_call, fill = Count)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(Count > 0, Count, "")),
              color = "white", size = 2.5) +
    scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                        limits = c(0, NA)) +
    labs(
      x = "Clone ID",
      y = "Ig Class",
      fill = "Count",
      title = glue("{HH}: Top {n_clones} clone spec_vj")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(glue("45_immcantation/plot/clones_spec_vj/{HH}_heatmap_across_clones_and_isotypes.png"), width = 8, height = 4.5, dpi = 1000)
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    plot_df <- spec_clones_vj[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call_majority %>% unique()
    j_gene <- plot_df$j_call_majority %>% unique()
    
    # Color by cell type
    plot_df %>%   
      ggplot(aes(y = sample_clean_fol, fill = L1_annotation)) +
      geom_bar() + 
      scale_fill_manual(values = L1_colors) + 
      geom_text(aes(x = Count, label = Count), 
                hjust = -0.2, color = "black",
                stat = "unique") +
      labs(
        title = glue("{HH}: Top {clone_nr} clone spec_vj"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
        y = ""
      ) +
      theme_bw()
    
    ggsave(glue("45_immcantation/plot/clones_spec_vj/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5)
    
    # Heatmap of sample_clean_fol and c_call (Isotype) and celltype
    plot_df %>%
      mutate(c_call = ifelse(c_call == "", NA, c_call)) %>%
      group_by(sample_clean_fol, c_call, L1_annotation) %>%
      summarise(Count = n(), .groups = "drop") %>%
      complete(sample_clean_fol, c_call, L1_annotation, 
               fill = list(Count = 0)) %>%
      ggplot(aes(x = sample_clean_fol, y = c_call, fill = Count)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(Count > 0, Count, "")),
                color = "white", size = 2.5) +
      scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                          limits = c(0, NA)) +
      facet_wrap(~ L1_annotation, nrow = 2) +
      labs(
        x = NULL,
        y = "Ig Class",
        fill = "Count",
        title = glue("{HH}: Top {clone_nr} clone spec_vj"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "right"
      )
    
    n_samples <- plot_df$sample_clean_fol %>% unique() %>% length()
    
    if (n_samples > 8){
      ggsave(glue("45_immcantation/plot/clones_spec_vj/{HH}_clone_nr_{clone_nr}_heatmap_across_samples_and_isotypes_and_celltype.png"), width = 20, height = 8, dpi = 1000)
    } else {
      ggsave(glue("45_immcantation/plot/clones_spec_vj/{HH}_clone_nr_{clone_nr}_heatmap_across_samples_and_isotypes_and_celltype.png"), width = 8, height = 6, dpi = 1000)
    }
    
  }
}


# ------------------------------------------------------------------------------
# Compare methods: novj VS vj 
# ------------------------------------------------------------------------------

# HH <- "HH119"
# 
# HH_spec_clones_vj <- spec_clones_vj[[HH]]
# HH_spec_clones_novj <- spec_clones_novj[[HH]]
# 
# # Not in the same order...
# table(HH_spec_clones_vj$cell_id %in% HH_spec_clones_novj$cell_id) 
# table(HH_spec_clones_vj$cell_id == HH_spec_clones_novj$cell_id) 
# length(HH_spec_clones_vj$cell_id)
# length(HH_spec_clones_novj$cell_id)
# 
# # reorder novj to match the order of vj
# HH_spec_clones_novj <- HH_spec_clones_novj[match(HH_spec_clones_vj$cell_id_seurat, HH_spec_clones_novj$cell_id_seurat), ]
# 
# # verify they now match
# table(HH_spec_clones_vj$cell_id_seurat == HH_spec_clones_novj$cell_id_seurat)
# 
# # Investigate clone overlap clone_ids
# 
# # get top N clones from each method
# top_n <- 50
# 
# top_vj <- HH_spec_clones_vj %>% count(clone_id, sort = TRUE) %>% slice_head(n = top_n) %>% pull(clone_id)
# top_novj <- HH_spec_clones_novj %>% count(clone_id, sort = TRUE) %>% slice_head(n = top_n) %>% pull(clone_id)
# 
# # build overlap table for top clones only
# overlap <- HH_spec_clones_vj %>%
#   select(cell_id_seurat, clone_vj = clone_id) %>%
#   left_join(HH_spec_clones_novj %>% select(cell_id_seurat, clone_novj = clone_id),
#             by = "cell_id_seurat") %>%
#   filter(clone_vj %in% top_vj, clone_novj %in% top_novj) %>%
#   count(clone_vj, clone_novj)
# 
# # heatmap
# ggplot(overlap, aes(x = factor(clone_novj), y = factor(clone_vj), fill = n)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   labs(x = "noVJ clone", y = "VJ clone", fill = "# cells") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
