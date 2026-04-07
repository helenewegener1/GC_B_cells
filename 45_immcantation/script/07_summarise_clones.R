library(tidyverse) 
library(glue)

# NOT DOING light_cluster.py SINCE IT'S TOO STRICT (removes cells without light chains)

# ------------------------------------------------------------------------------
# Having a look 
# ------------------------------------------------------------------------------

# x <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"
# # x <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
# # 
# light <- read.delim(glue("45_immcantation/out/{x}/{x}_light_germ-pass.tsv"))
# light_qc <- read.delim(glue("45_immcantation/out/{x}/{x}_light_germ-pass_QC.tsv"))
# heavy <- read.delim(glue("45_immcantation/out/{x}/{x}_heavy_germ-pass.tsv"))
# vj <- read.delim(glue("45_immcantation/out/{x}/{x}_heavy_germ-pass_clone-pass.tsv"))
# clone_10x <- read.delim(glue("45_immcantation/out/{x}/{x}_10X_clone-pass.tsv"))
#
# nrow(light)
# nrow(light_qc)
# nrow(heavy)
# nrow(vj)
# nrow(clone_10x)
# 
# clone_10x %>% count(clone_id, sort = T) %>% head()
# 
# table(light$cell_id %in% heavy$cell_id)
#
# table(light$cell_id %in% vj$cell_id)
# table(heavy$cell_id %in% vj$cell_id)
#
# table(light$cell_id %in% clone_10x$cell_id)
# table(heavy$cell_id %in% clone_10x$cell_id)
# 
# table(vj$cell_id %in% clone_10x$cell_id)

# Compare rds objects 

# spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")
# resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list.rds")
# 
# spec_clones_vj$HH119 %>% count(clone_id, sort = TRUE) %>% head(10)
# 
# light_qc %>% nrow()
# spec_clones_vj$HH119 %>% filter(sample_id == x) %>% nrow()
# resolve_LC_list$HH119 %>% filter(sample_id == x) %>% nrow()
# 
# spec_clones_vj$HH119 %>% filter(sample_id == x) %>% count(clone_id, sort = TRUE) %>% head(10)
# resolve_LC_list$HH119 %>% filter(sample_id == x) %>% count(clone_id, locus, sort = TRUE) %>% head(10)


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

