library(tidyverse) 
library(glue)

# ------------------------------------------------------------------------------
# Having a look 
# ------------------------------------------------------------------------------

heavy_clones <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")
light_clones <- readRDS("45_immcantation/out/rds/05_spec_clones_novj_light.rds")

HH <- "HH117"

heavy_clones_HH <- heavy_clones[[HH]]
light_clones_HH <- light_clones[[HH]]

nrow(heavy_clones_HH)
nrow(light_clones_HH)

# Check cell_id match 
table(heavy_clones_HH$cell_id %in% light_clones_HH$cell_id) 
table(light_clones_HH$cell_id %in% heavy_clones_HH$cell_id)

# N clones
heavy_clones_HH$clone_id %>% unique() %>% length()
light_clones_HH$clone_id %>% unique() %>% length()

# Combine
# Remove light chains without heavy chain match and keep heavy chains without light chain match. 
clones_HH <- heavy_clones_HH %>% 
  left_join(light_clones_HH, by = "cell_id", suffix = c("_heavy", "_light")) %>% 
  mutate(
    clone_id_combine = paste(clone_id_heavy, clone_id_light, sep = "_")
  )

# N clones
clones_HH$clone_id_heavy %>% unique() %>% length()
clones_HH$clone_id_light %>% unique() %>% length()
clones_HH$clone_id_combine %>% unique() %>% length()

# Inspect clones 
clones_HH %>% count(clone_id_heavy, sort = TRUE)

top_clone <- clones_HH %>% count(clone_id_heavy, sort = TRUE) %>% head(n = 1) %>% pull(clone_id_heavy)

clones_HH %>% filter(clone_id_heavy == top_clone) %>% count(clone_id_combine, sort = TRUE)

# ------------------------------------------------------------------------------
# Combine heavy and light chain data
# ------------------------------------------------------------------------------

patients <- names(heavy_clones)

clones_combined <- lapply(patients, function(HH) {
  
  # HH <- "HH117"
  
  heavy_clones_HH <- heavy_clones[[HH]]
  light_clones_HH <- light_clones[[HH]]
  
  # Build a lookup: cell_id -> clone_id for each locus
  heavy_lookup <- heavy_clones_HH %>% select(cell_id, clone_id_heavy = clone_id)
  light_lookup <- light_clones_HH %>% select(cell_id, clone_id_light = clone_id)
  
  # Stack into long format
  rbind(heavy_clones_HH, light_clones_HH) %>%
    # Join both lookups so every row knows both clone IDs
    left_join(heavy_lookup, by = "cell_id") %>%
    left_join(light_lookup, by = "cell_id") %>%
    mutate(
      clone_id_combine = paste(clone_id_heavy, clone_id_light, sep = "_")
    )
  
}) %>% setNames(patients)


# ------------------------------------------------------------------------------
# Have a look at clone_id_combine
# ------------------------------------------------------------------------------

HH <- "HH117"

top_clone <- clones_combined[[HH]] %>% filter(locus == "IGH") %>% count(clone_id, sort = TRUE) %>% head(n = 1) %>% pull(clone_id)
clones_combined[[HH]] %>% filter(locus == "IGH") %>% filter(clone_id == top_clone) %>% count(clone_id_combine, sort = TRUE)

# ------------------------------------------------------------------------------
# Define subclones based on light chain data
# ------------------------------------------------------------------------------

library(stringdist)

clones_NA_resolved <- lapply(patients, function(HH) {
  
  # HH <- "HH117"
  clones_HH <- clones_combined[[HH]]
  
  # Work on cell level using only IGH rows for the imputation logic
  igh_rows <- clones_HH %>% filter(locus == "IGH")
  
  igh_has_light <- igh_rows %>% filter(!is.na(clone_id_light))
  igh_no_light  <- igh_rows %>% filter(is.na(clone_id_light))
  
  if (nrow(igh_no_light) == 0) return(clones_HH %>% mutate(light_imputed = FALSE))
  
  # Impute on IGH rows only
  igh_imputed <- igh_no_light %>%
    rowwise() %>%
    mutate(
      clone_id_light = {
        current_heavy_clone <- pick(clone_id_heavy)$clone_id_heavy
        current_seq         <- pick(sequence_alignment)$sequence_alignment
        
        candidates <- igh_has_light %>%
          filter(clone_id_heavy == current_heavy_clone)
        
        if (nrow(candidates) == 0) {
          NA_character_
        } else {
          mean_dists <- candidates %>%
            group_by(clone_id_light) %>%
            summarise(
              mean_dist = mean(stringdist(
                current_seq,
                sequence_alignment,
                method = "hamming"
              ), na.rm = TRUE)
            )
          mean_dists$clone_id_light[which.min(mean_dists$mean_dist)]
        }
      },
      clone_id_combine = paste(clone_id_heavy, clone_id_light, sep = "_"),
      light_imputed = TRUE
    ) %>%
    ungroup()
  
  # Lookup table: only the columns we want to patch, cleanly named
  imputed_lookup <- igh_imputed %>%
    select(cell_id, 
           new_clone_id_light   = clone_id_light, 
           new_clone_id_combine = clone_id_combine,
           new_light_imputed    = light_imputed)
  
  # Apply imputed values back to ALL rows of those cells
  clones_HH %>%
    left_join(imputed_lookup, by = "cell_id") %>%
    mutate(
      clone_id_light   = coalesce(new_clone_id_light,   clone_id_light),
      clone_id_combine = coalesce(new_clone_id_combine, clone_id_combine),
      light_imputed    = coalesce(new_light_imputed,    FALSE)
    ) %>%
    select(-new_clone_id_light, -new_clone_id_combine, -new_light_imputed)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Have a look at clone_id_combine before and after NA resolvement 
# ------------------------------------------------------------------------------

HH <- "HH117"
# 
# Get top heavy chain clone
top_clone <- clones_combined[[HH]] %>% filter(locus == "IGH") %>% count(clone_id, sort = TRUE) %>% head(n = 1) %>% pull(clone_id)

# Before
clones_combined[[HH]] %>% 
  filter(locus == "IGH") %>% 
  filter(clone_id == top_clone) %>% 
  count(clone_id_combine, sort = TRUE)

# After
clones_NA_resolved[[HH]] %>% 
  filter(locus == "IGH") %>%
  filter(clone_id_heavy == top_clone) %>% 
  count(clone_id_combine, sort = TRUE)

saveRDS(clones_NA_resolved, "45_immcantation/out/rds/07_clones_combined_NA_resolved.rds")
# clones_NA_resolved <- readRDS("45_immcantation/out/rds/07_clones_combined_NA_resolved.rds")

# ------------------------------------------------------------------------------
# Define top clones
# ------------------------------------------------------------------------------

patients <- names(clones_NA_resolved)

top_GC_clones <- lapply(patients, function(HH) {
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- clones_NA_resolved[[HH]] %>%
    filter(locus == "IGH" & celltype_broad == "GC_B_cells") %>%
    group_by(clone_id) %>%
    summarise(n_samples = n_distinct(sample_clean_fol)) %>%
    filter(n_samples >= 2) %>%
    pull(clone_id)
  
  # rank those clones by total size (all cell types) and take top 10
  clones_NA_resolved[[HH]] %>%
    filter(locus == "IGH" & clone_id %in% GC_clones) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_id)
  
}) %>% setNames(patients)


# ------------------------------------------------------------------------------
# Visualize top clones vj
# ------------------------------------------------------------------------------

for (HH in patients){
  
  # HH <- "HH117"
  HH_top_clones <- top_GC_clones[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    df <- clones_NA_resolved[[HH]] %>%
      filter(str_starts(clone_id_combine, paste0(clone, "_"))) %>%
      mutate(
        sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev(),
        clone_subgroup_genes = as.character(clone_id) %>% paste(v_call_majority, j_call_majority, junction_length, sep = "_")
      )
    
    # Meta data
    n_cells <- df %>% filter(locus == "IGH") %>% nrow()
    v_gene <- df %>% filter(locus == "IGH") %>% pull(v_call_majority) %>% unique()
    j_gene <- df %>% filter(locus == "IGH") %>% pull(j_call_majority) %>% unique()
    
    # Plot data
    plot_df <- df %>% 
      filter(locus != "IGH") %>% 
      count(clone_subgroup_genes, celltype_broad, sort = TRUE)
    
    plot_df %>%
      mutate(clone_subgroup_genes = fct_reorder(
        clone_subgroup_genes,
        as.numeric(str_extract(clone_subgroup_genes, "^\\d+"))
      )) %>% 
      ggplot(aes(y = clone_subgroup_genes, x = celltype_broad, fill = n)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(n > 0, n, "")),
                color = "white", size = 2.5) +
      scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                          limits = c(0, NA)) +
      labs(
        title = glue("{HH}: Top {clone_nr} GCB clone - Light chain subgroups after resolveLightChains"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells heavy chain: {n_cells}\nV gene heavy chain: {v_gene}\n J gene heavy chain: {j_gene}"),
        y = ""
      ) +
      theme_classic() + 
      theme(legend.position = "none")
    
    # ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_resolveLightChains/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 8, height = 8.5)
    ggsave(glue("45_immcantation/plot/GC_clones_HW_combined/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 8, height = 8.5)
    
  }
}


# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# clone_10x <- readRDS("45_immcantation/out/rds/clone_10x.rds")
# patients <- names(clone_10x)

# # Load light chain corrected 
# files <- list.files("45_immcantation/out/")
# sample_names <- files[1:length(files)-1]
# 
# clone_10x_list <- lapply(sample_names, function(x){
#   
#   clone_10x <- read.delim(glue("45_immcantation/out/{x}/{x}_10X_clone-pass.tsv"))
#   
#   return(clone_10x)
#   
# }) %>% setNames(sample_names)
# 
# clone_10x_combined <- bind_rows(clone_10x_list)



# ------------------------------------------------------------------------------
# Alter format
# ------------------------------------------------------------------------------

# clone_10x <- list(
#   "HH117" = clone_10x_combined %>% filter(subject_id == "HH117"),
#   "HH119" = clone_10x_combined %>% filter(subject_id == "HH119")
# )

# ------------------------------------------------------------------------------
# Define top clones
# ------------------------------------------------------------------------------

patients <- names(clone_10x)

top_GC_clones_vj <- lapply(patients, function(HH) {
  
  # # find clones that contain at least 1 GC cell
  # GC_clones <- clone_10x[[HH]] %>%
  #   filter(celltype_broad == "GC_B_cells") %>%
  #   pull(clone_id) %>%
  #   unique()
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- clone_10x[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    group_by(clone_id) %>%
    summarise(n_samples = n_distinct(sample_clean_fol)) %>%
    filter(n_samples >= 2) %>%
    pull(clone_id)
  
  # rank those clones by total size (all cell types) and take top 10
  clone_10x[[HH]] %>%
    filter(clone_id %in% GC_clones) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_id)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Look at top clones
# ------------------------------------------------------------------------------

lapply(patients, function(HH){
  
  # HH <- "HH119"
  # HH <- "HH117"
  
  clone_10x[[HH]] %>% 
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
  
  # HH <- "HH119"
  HH_top_clones <- top_GC_clones_vj[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 2
    clone <- HH_top_clones[clone_nr] %>% str_split_i("_", 1)
    
    plot_df <- clone_10x[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    j_gene <- plot_df$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    
    # Color by cell type
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
    
    ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_ligth_chain_corrected/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5)
    
    # Heatmap of sample_clean_fol and c_call (Isotype)
    plot_df %>% 
      mutate(c_call = ifelse(c_call == "", NA, c_call)) %>% 
      group_by(sample_clean_fol, c_call) %>%
      summarise(Count = n(), .groups = "drop") %>%
      complete(sample_clean_fol, c_call, fill = list(Count = 0)) %>%
      ggplot(aes(x = sample_clean_fol, y = c_call, fill = Count)) + 
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(Count > 0, Count, "0")), 
                color = "white", size = 2.5) +
      scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                          limits = c(0, NA)) +
      labs(
        x = NULL, 
        y = "Ig Class", 
        fill = "Count",  
        title = glue("{HH}: Top {clone_nr} GCB clone spec_vj"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}"),
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_ligth_chain_corrected/{HH}_clone_nr_{clone_nr}_heatmap_across_samples_and_isotypes.png"), width = 15, height = 8.5)
    
    # Heatmap of sample_clean_fol and c_call (Isotype) and celltype
    plot_df %>%
      mutate(c_call = ifelse(c_call == "", NA, c_call)) %>%
      group_by(sample_clean_fol, c_call, celltype_broad) %>%
      summarise(Count = n(), .groups = "drop") %>%
      complete(sample_clean_fol, c_call, celltype_broad, 
               fill = list(Count = 0)) %>%
      ggplot(aes(x = sample_clean_fol, y = c_call, fill = Count)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(Count > 0, Count, "")),
                color = "white", size = 2.5) +
      scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                          limits = c(0, NA)) +
      facet_wrap(~ celltype_broad, nrow = 2) +
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
      ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_ligth_chain_corrected/{HH}_clone_nr_{clone_nr}_heatmap_across_samples_and_isotypes_and_celltype.png"), width = 17, height = 8, dpi = 1000)
    } else {
      ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_ligth_chain_corrected/{HH}_clone_nr_{clone_nr}_heatmap_across_samples_and_isotypes_and_celltype.png"), width = 8, height = 8, dpi = 1000)
    }
    
  }
}

