library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
source("10_broad_annotation/script/color_palette.R")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

# resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")
# plot_version <- "plot"

resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list_gmm_threshold_germlined.rds")
plot_version <- "plot_gmm_threshold"

patients <- names(resolve_LC_list)

ncol(seurat_integrated)

# df <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")
# 
# df$HH119 %>% count(clone_id, sort = TRUE)
# 
# resolve_LC_list$HH117 %>% 
#   filter(locus == "IGH") %>% 
#   count(clone_id, sort = TRUE)

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

# resolve_LC_list$HH117 %>% filter(clone_subgroup_id == "578_1" & locus == "IGH") %>% 
#   count(j_call, v_call, junction_length, sort = TRUE) 
# 
# 
# resolve_LC_list$HH117 %>% filter(clone_subgroup_id == "578_1" & locus != "IGH") %>% 
#   count(j_call, v_call, junction_length, sort = TRUE)

# Get sequences 

# resolve_LC_list$HH117 %>% filter(clone_subgroup_id == "578_1" & locus == "IGH") %>% pull(sequence)

# ==============================================================================
# Export BCR meta data to Gina 
# ==============================================================================

# meta_4_Gina_list <- lapply(patients, function(HH){
# 
#   # HH <- "HH119"
# 
#   seurat_obj <- subset(seurat_integrated, patient == HH)
#   resolve_LC_HH <- resolve_LC_list[[HH]] %>% filter(locus == "IGH")
# 
#   # Check IDs
#   seurat_obj %>% colnames() %>% head()
#   resolve_LC_HH$cell_id %>% head()
# 
#   seurat_obj %>% colnames() %>% length()
#   resolve_LC_HH$cell_id %>% length()
# 
#   (seurat_obj %>% colnames() %>% length()) == (seurat_obj %>% colnames() %>% unique() %>% length())
#   (resolve_LC_HH$cell_id %>% length()) == (resolve_LC_HH$cell_id %>% unique() %>% length())
# 
# 
#   # Wrangle IDs
#   seurat_ids <- seurat_obj %>% colnames()
#   LC_ids <- resolve_LC_HH$cell_id_seurat %>% str_remove(".*?_")
# 
#   # # IDs test 
#   # seurat_ids_sub <- seurat_ids %>% str_split_i("_", 2)
#   # table(seurat_ids_sub, seurat_obj$sample_clean)
#   # 
#   # LC_ids_sub <- LC_ids %>% str_split_i("_", 2)
#   # table(LC_ids_sub, resolve_LC_HH$sample_clean)
#   # # End
# 
#   (seurat_ids %>% length()) == (seurat_ids %>% unique() %>% length())
#   (LC_ids %>% length()) == (LC_ids %>% unique() %>% length())
# 
#   table(LC_ids %in% seurat_ids)
# 
#   # Prep for merge
#   resolve_LC_HH_meta <- resolve_LC_HH %>%
#     mutate(cell_id_seurat_clean = str_remove(cell_id_seurat, ".*?_")) %>%
#     select(cell_id_seurat_clean, c_call, clone_subgroup_id)
# 
#   # Merge and create final meta data for Gina
#   meta_4_Gina <- seurat_obj[[]] %>%
#     select(manual_ADT_class, manual_ADT_ID, manual_ADT_full_ID, sample) %>%
#     rownames_to_column("cell_id_seurat_clean") %>%
#     left_join(resolve_LC_HH_meta, by = "cell_id_seurat_clean") %>%
#     column_to_rownames("cell_id_seurat_clean")
#   
#   
#   # table(meta_4_Gina$L1_annotation, meta_4_Gina$c_call, useNA = "always")
# 
#   # Check
#   meta_4_Gina %>% count(clone_subgroup_id, sort = TRUE)  %>% head()
# 
#   return(meta_4_Gina)
# 
# }) %>% setNames(patients)
# 
# 
# saveRDS(meta_4_Gina_list %>% bind_rows(), "45_immcantation/out/rds/meta_4_Gina_list.rds")

# meta <- readRDS("45_immcantation/out/rds/meta_4_Gina_list.rds")
# 
# meta %>% filter(patient_id == "HH117")
# 
# resolve_LC_list$HH117$L1_annotation %>% table()


# ==============================================================================
# All cells: Summary of cell types in follicles.
# ==============================================================================

outdir_1 <- glue("45_immcantation/{plot_version}/15_poster_figures/Follicle_cell_types")
dir.create(outdir_1, recursive = TRUE)

lapply(patients, function(HH){
  
  # HH <- "HH119"
  p <- patient_names[[HH]]
  
  seurat_obj <- subset(seurat_integrated, patient == HH)
  
  # Define LP samples
  LP_samples <- grep("LP", seurat_obj[[]]$sample_clean, value = TRUE) %>% unique()
  
  # How many Tfh cells with BCR?
  seurat_obj[[]] %>% filter((!is.na(bcr_productive_contig_1) & !is.na(bcr_productive_contig_2) & L1_annotation == "Tfh_cells")) 
  
  # Clean meta data and prep for plotting 
  seurat_meta_clean <- seurat_obj[[]] %>%  
    mutate(L1_annotation = ifelse(L1_annotation == "GC_Bcells", "GC_B_cells", L1_annotation)) %>% 
    filter(
      (str_detect(L1_annotation, "Contamination", negate = TRUE)), # Remove contamination
      !(!is.na(bcr_productive_contig_1) & !is.na(bcr_productive_contig_2) & L1_annotation == "Tfh_cells"), # Remove Tfh cells with BCR
      !(L1_annotation == "GC_B_cells" & sample_clean %in% LP_samples), # Remove GC B cells in LP samples
    ) %>% mutate(
      sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-")),
      sample_clean_plot = fct_infreq(sample_clean_plot) #%>% fct_rev()
    ) %>% 
    add_count(sample_clean_plot, name = "Count") 
  
  # Across follicles 
  HH_fol_sample_clean <- seurat_meta_clean %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))

  # Count
  if (HH == "HH117"){
    width <- 12 
  } else if (HH == "HH119"){
    width <- 15
  }
  png(glue("{outdir_1}/{HH}_N_cells_across_follicles.png"), width = width, height = 7, res = 1000, units = "in")
  
  print(
    seurat_meta_clean %>% 
      filter(!is.na(manual_ADT_ID)) %>% 
      mutate(
        manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
      ) %>% 
      ggplot(aes(x = manual_ADT_ID_plot, fill = L1_annotation)) +
      geom_bar() + 
      scale_fill_manual(
        values = L1_colors, 
        labels = cell_type_names
      ) + 
      scale_x_continuous(
        breaks = function(x) seq(1, ceiling(max(x)), by = 1),
        limits = c(0.5, NA),
        expand = c(0, 0.5)
      ) + 
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Count", 
        title = glue ("{p}: N B cells across {HH_fol_sample_clean} follicles"),
        fill = "Cell type"
      ) + 
      theme(
        plot.title = element_text(face = "bold", size = 26),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)
      )
  )
  
  dev.off()
  
})

# ==============================================================================
# G B cells: Summary of follicles and isotypes
# ==============================================================================

outdir_2 <- glue("45_immcantation/{plot_version}/15_poster_figures/Follicle_GC_B_cells_isotypes")
dir.create(outdir_2, recursive = TRUE)

lapply(patients, function(HH){
  
  # HH <- "HH119"
  p <- patient_names[[HH]]
  
  plot_df <- resolve_LC_list[[HH]] %>% 
    filter(locus == "IGH") %>% 
    mutate(
      sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-")),
      sample_clean_plot = fct_infreq(sample_clean_plot) #%>% fct_rev()
    ) %>%
    add_count(sample_clean_plot, name = "Count") 
  

  # Isotype
  
  ## Freq
  if (HH == "HH117"){
    width <- 12 
  } else if (HH == "HH119"){
    width <- 15
  }
  
  png(glue("{outdir_2}/{HH}_Isotype_freq_across_follicles.png"), width = width, height = 7, res = 1000, units = "in")
  
  print(
    plot_df %>% 
      filter(
        !is.na(manual_ADT_ID), 
        L1_annotation == "GC_B_cells",
        !is.na(c_call)
      ) %>%
      mutate(
        manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
      ) %>% 
      ggplot(aes(x = manual_ADT_ID_plot, fill = c_call)) + 
      geom_bar(position = "fill") +            # <-- changed
      scale_fill_manual(values = isotype_colors_custom) +
      scale_y_continuous(labels = scales::percent) +   # <-- optional: show % instead of 0-1
      scale_x_continuous(
        breaks = function(x) seq(1, ceiling(max(x)), by = 1),
        limits = c(0.5, NA),
        expand = c(0, 0.5)
      ) + 
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Frequency", 
        title = glue("{p}: Isotypes across GC B cells in {HH_fol_sample_clean} follicles"),
        fill = "Isotype"
      ) + 
      theme(
        plot.title = element_text(face = "bold", size = 26),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)
      )
  )
  
  dev.off()
  
})

# ==============================================================================
# CLONE BAR PLOTS
# ==============================================================================

# ------------------------------------------------------------------------------
# Plot data - GC clones in at least two samples
# ------------------------------------------------------------------------------

outdir_3 <- glue("45_immcantation/{plot_version}/15_poster_figures/bar_plots_top_GC_clones")
dir.create(outdir_3, recursive = TRUE)

# Define top GC clone
top_GC_clones <- lapply(patients, function(HH) {
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- resolve_LC_list[[HH]] %>%
    filter(locus == "IGH") %>% 
    filter(L1_annotation == "GC_B_cells") %>%
    group_by(clone_subgroup_id) %>%
    summarise(n_samples = n_distinct(sample_clean_fol)) %>%
    filter(n_samples >= 2) %>%
    pull(clone_subgroup_id)
  
  # rank those clones by total size (all cell types) and take top 10
  resolve_LC_list[[HH]] %>%
    filter(locus == "IGH") %>% 
    filter(clone_subgroup_id %in% GC_clones) %>%
    count(clone_subgroup_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_subgroup_id)
  
}) %>% setNames(patients)


# Look at top clones
lapply(patients, function(HH){
  
  resolve_LC_list[[HH]] %>% 
    filter(locus == "IGH") %>% 
    filter(clone_subgroup_id %in% top_GC_clones[[HH]]) %>% 
    count(clone_subgroup_id, sort = TRUE) 
  
}) %>% setNames(patients)


# Visualize top clones
for (HH in patients){
  
  # HH <- "HH117"
  HH_top_clones <- top_GC_clones[[HH]]
  n_clones <- length(HH_top_clones)
  
  p <- patient_names[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    # Title
    if (clone_nr == 1){
      clone_definition <- "Largest clone with GC B cells across samples"
    } else if (clone_nr == 2) {
      clone_definition <- "Second largest clone with GC B cells across samples"
    } else {
      clone_definition <- glue("{clone_nr}. largest clone with GC B cells across samples")
    }
    
    plot_df <- resolve_LC_list[[HH]] %>% 
      filter(locus == "IGH") %>% 
      filter(clone_subgroup_id == clone) %>%
      mutate(
        sample_clean_fol_plot = sample_clean_fol %>% str_remove_all(glue("{HH}-")),
        sample_clean_fol_plot = ifelse(str_detect(sample_clean_fol_plot, "Fol"), 
                                       str_split_i(sample_clean_fol_plot, "_", 2) %>% str_replace("Fol", "Follicle"),
                                       sample_clean_fol_plot),
        sample_clean_fol_plot = fct_infreq(sample_clean_fol_plot) %>% fct_rev()
      ) %>%
      add_count(sample_clean_fol_plot, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call_majority %>% unique() %>% str_split_i(",", 1)
    j_gene <- plot_df$j_call_majority %>% unique() %>% str_split_i(",", 1)
    
    png(glue("{outdir_3}/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5, res = 1000, units = "in")
    
    # Color by cell type
    print(
      plot_df %>%
        ggplot(aes(y = sample_clean_fol_plot, fill = L1_annotation)) +
        geom_bar() +
        scale_fill_manual(
          values = L1_colors,
          labels = cell_type_names
        ) +
        geom_text(aes(x = Count, label = Count),
                  hjust = -0.2, color = "black",
                  stat = "unique") +
        labs(
          title = glue("{p}: {clone_definition}"),
          # subtitle = glue("Clone ID: {clone}"),
          caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
          y = "Samples",
          fill = "Cell type"
        ) +
        # theme_minimal()
        theme_classic() + 
        theme(plot.title = element_text(face = "bold", size = 16))
    )

    dev.off()


    # Isotype heatmap
    
    if (HH == "HH119" & clone_nr == 1){
      width <- 20
      height <- 8
    } else {
      width <- 12
      height <- 4
    }
    
    png(glue("{outdir_3}/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types_and_isotype.png"), width = width, height = height, units = "in", res = 1000)
    
    print(
      plot_df %>%
        filter(!c_call %in% c("IGHE", "") & !is.na(c_call)) %>%
        group_by(sample_clean_fol_plot, c_call, L1_annotation) %>%
        summarise(Count = n(), .groups = "drop") %>%
        complete(sample_clean_fol_plot, c_call, L1_annotation,
                 fill = list(Count = 0)) %>%
        ggplot(aes(x = L1_annotation, y = c_call, fill = Count)) +
        geom_tile(color = "white", linewidth = 0.3) +
        geom_text(aes(label = ifelse(Count > 0, Count, "")),
                  color = "white", size = 2.5) +
        scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                            limits = c(0, NA)) +
        facet_grid(~ sample_clean_fol_plot) +
        labs(
          x = "Cell type",
          y = "Ig Class",
          title = glue("{p}: {clone_definition}"),
          # subtitle = glue("Clone ID: {clone}"),
          caption = glue("N cells: {n_cells}\nV gene: {v_gene}\nJ gene: {j_gene}")
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 9),
          panel.grid = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 16)
        )
    )
  
    dev.off()
    
    # Upset plot 
    plot_df_upset <- plot_df %>%
      filter(!c_call %in% c("IGHE", "") & !is.na(c_call))
    
    # Assuming plot_df has one row per sequence with c_call and L1_annotation columns
    upset_input_c <- plot_df_upset %>%
      mutate(present = 1L) %>%
      pivot_wider(
        id_cols = sequence_id,
        names_from = c_call,      
        values_from = present,
        values_fill = 0L
      ) %>%
      as.data.frame() %>% 
      select(where(~ any(. != 0)))
    
    upset_input_L1 <- plot_df_upset %>%
      mutate(present = 1L) %>%
      pivot_wider(
        id_cols = sequence_id,
        names_from = L1_annotation, 
        values_from = present,
        values_fill = 0L
      ) %>%
      as.data.frame() %>% 
      select(where(~ any(. != 0)))
    
    upset_input_sample <- plot_df_upset %>%
      mutate(present = 1L) %>%
      pivot_wider(
        id_cols = sequence_id,
        names_from = sample_clean_fol_plot, 
        values_from = present,
        values_fill = 0L
      ) %>%
      as.data.frame() %>% 
      select(where(~ any(. != 0)))
    
    upset_input_final <- full_join(upset_input_c, upset_input_L1, by = "sequence_id") %>% full_join(upset_input_sample, by = "sequence_id")

    # Define order
    c_order <- upset_input_c %>% column_to_rownames("sequence_id") %>% colSums() %>% sort() %>% names()
    L1_order <- upset_input_L1 %>% column_to_rownames("sequence_id") %>% colSums() %>% sort() %>% names()
    sample_order <- upset_input_sample %>% column_to_rownames("sequence_id") %>% colSums() %>% sort() %>% names()
    
    set_order <- c(L1_order, c_order, sample_order)
    
    set_colors <- c(
      rep("blue", length(L1_order)),
      rep("darkred", length(c_order)),
      rep("darkgreen", length(sample_order))
    )

    if (HH == "HH119" & clone_nr == 1){
      width <- 8
      height <- 9.5
      nintersects <- 27
      more_than_N_cells <- 10
    } else if (HH == "HH117" & clone_nr == 1){
      width <- 5.5
      height <- 5.5
      more_than_N_cells <- 2
      nintersects <- 8
    } else {
      width <- 11
      height <- 4
      more_than_N_cells <- 1
    }
    
    # Prep saving plot 
    png(glue("{outdir_3}/{HH}_clone_nr_{clone_nr}_upsetplot.png"), width = width, height = height, units = "in", res = 1000)
    
    # Plot 
    print(UpSetR::upset(
      upset_input_final,
      sets = set_order,
      nintersects = nintersects,
      order.by = "freq",
      keep.order = TRUE, 
      sets.bar.color = set_colors,
      mb.ratio = c(0.4, 0.6),
      line.size = 0.3
    ))
    
    # White background rectangle behind title
    grid.rect(
      x = 0.50, y = 0.98,
      width = 0.7, height = 0.05,   # adjust size to fit your text
      gp = gpar(fill = "white", col = NA)  # col = NA removes border
    )
    
    # Title
    grid.text(
      glue("{p}: {clone_definition}"),
      x = 0.50, y = 0.97,          # adjust position as needed
      gp = gpar(fontsize = 12, fontface = "bold")
    )
    
    # Caption 
    grid.text(
      glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}\nCombinations with >= {more_than_N_cells} cells are shown."),
      just = "right", 
      x = 0.99, y = 0.03,          # adjust position as needed
      gp = gpar(fontsize = 6)
    )
    
    # Legend
    legend_labels <- c("Sample", "Isotype", "Cell type")
    legend_colors <- c("darkgreen", "darkred", "blue")
    
    x_start <- 0.07
    y_start <- 0.80
    gap     <- 0.03
    
    for (i in seq_along(legend_labels)) {
      grid.rect(
        x = x_start, y = y_start - (i - 1) * gap,
        width = 0.02, height = 0.020,
        just = "left",
        gp = gpar(fill = legend_colors[i], col = NA)
      )
      grid.text(
        legend_labels[i],
        x = x_start + 0.03, y = y_start - (i - 1) * gap,
        just = "left",
        gp = gpar(fontsize = 9)
      )
    }
    
    dev.off()

    
  }
}

# ------------------------------------------------------------------------------
# Plot data - All clones (no restrictions)
# ------------------------------------------------------------------------------

outdir_4 <- glue("45_immcantation/{plot_version}/15_poster_figures/bar_plots_top_all_clones")
dir.create(outdir_4, recursive = TRUE)

# Define top GC clone
top_clones <- lapply(patients, function(HH) {
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- resolve_LC_list[[HH]] %>%
    filter(locus == "IGH") %>% 
    group_by(clone_subgroup_id) %>%
    pull(clone_subgroup_id)
  
  # rank those clones by total size (all cell types) and take top 10
  resolve_LC_list[[HH]] %>%
    filter(locus == "IGH") %>% 
    filter(clone_subgroup_id %in% GC_clones) %>%
    count(clone_subgroup_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_subgroup_id)
  
}) %>% setNames(patients)


# Look at top clones
lapply(patients, function(HH){
  
  resolve_LC_list[[HH]] %>% 
    filter(locus == "IGH") %>% 
    filter(clone_subgroup_id %in% top_clones[[HH]]) %>% 
    count(clone_subgroup_id, sort = TRUE) 
  
}) %>% setNames(patients)

# Visualize top clones
for (HH in patients){
  
  # HH <- "HH119"
  HH_top_clones <- top_clones[[HH]]
  n_clones <- length(HH_top_clones)
  
  p <- patient_names[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    # Title
    if (clone_nr == 1){
      clone_definition <- "Largest clone"
    } else if (clone_nr == 2) {
      clone_definition <- "Second largest clone"
    } else {
      clone_definition <- glue("{clone_nr}. largest clone")
    }
    
    plot_df <- resolve_LC_list[[HH]] %>% 
      filter(locus == "IGH") %>% 
      filter(clone_subgroup_id == clone) %>%
      mutate(
        sample_clean_fol_plot = sample_clean_fol %>% str_remove_all(glue("{HH}-")),
        sample_clean_fol_plot = ifelse(str_detect(sample_clean_fol_plot, "Fol"), 
                                       str_split_i(sample_clean_fol_plot, "_", 2) %>% str_replace("Fol", "Follicle"),
                                       sample_clean_fol_plot),
        sample_clean_fol_plot = fct_infreq(sample_clean_fol_plot) %>% fct_rev()
      ) %>%
      add_count(sample_clean_fol_plot, name = "Count") 
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call_majority %>% unique() %>% str_split_i(",", 1)
    j_gene <- plot_df$j_call_majority %>% unique() %>% str_split_i(",", 1)
    
    # Color by cell type
    
    png(glue("{outdir_4}/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5, units = "in", res = 1000)
    
    print(
      plot_df %>%   
        ggplot(aes(y = sample_clean_fol_plot, fill = L1_annotation)) +
        geom_bar() + 
        scale_fill_manual(
          values = L1_colors, 
          labels = cell_type_names
        ) + 
        geom_text(aes(x = Count, label = Count), 
                  hjust = -0.2, color = "black",
                  stat = "unique") +
        labs(
          title = glue("{p}: {clone_definition}"),
          # subtitle = glue("Clone ID: {clone}"),
          caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
          y = "Samples",
          fill = "Cell type"
        ) +
        # theme_minimal()
        theme_classic() + 
        theme(plot.title = element_text(face = "bold", size = 16))
    )
    
    dev.off()
    
    # ISOTYPE HEATMAP
    
    if (HH == "HH119" & clone_nr == 1){
      width <- 20
      height <- 8
    } else {
      width <- 12
      height <- 4
    }
    
    png(glue("{outdir_4}/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types_and_isotype.png"), width = width, height = height, units = "in", res = 1000)
    
    print(
      plot_df %>%
        filter(!c_call %in% c("IGHE", "") & !is.na(c_call)) %>% 
        group_by(sample_clean_fol_plot, c_call, L1_annotation) %>%
        summarise(Count = n(), .groups = "drop") %>%
        complete(sample_clean_fol_plot, c_call, L1_annotation, 
                 fill = list(Count = 0)) %>%
        ggplot(aes(x = L1_annotation, y = c_call, fill = Count)) +
        geom_tile(color = "white", linewidth = 0.3) +
        geom_text(aes(label = ifelse(Count > 0, Count, "")),
                  color = "white", size = 2.5) +
        scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                            limits = c(0, NA)) +
        facet_grid(~ sample_clean_fol_plot) +
        labs(
          x = "Cell type",
          y = "Ig Class",
          title = glue("{p}: {clone_definition}"),
          # subtitle = glue("Clone ID: {clone}"),
          caption = glue("N cells: {n_cells}\nV gene: {v_gene}\nJ gene: {j_gene}")
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 9),
          panel.grid = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 16)
        )
    )
    
    dev.off()
    
    # if (HH == "HH119" & clone_nr == 1){
    #   width <- 20
    #   height <- 8
    # } else {
    #   width <- 12
    #   height <- 4
    # }
    # 
    # ggsave(glue("{outdir_4}/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types_and_isotype.png"), width = width, height = height, dpi = 1000)
    
  }
}

# ==============================================================================
# ISOTYPES 
# ==============================================================================

outdir_5 <- glue("45_immcantation/{plot_version}/15_poster_figures/isotypes")
dir.create(outdir_5, recursive = TRUE)

lapply(patients, function(HH){
  
  # HH <- "HH119"
  p <- patient_names[[HH]]
  
  plot_df <- resolve_LC_list[[HH]] %>% 
    filter(locus == "IGH") %>% 
    mutate(
      sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-"))
    ) %>%
    add_count(sample_clean_plot, name = "Count") 
  
  png(glue("{outdir_5}/{HH}_isotypes_samples_cell_types.png"), width = 9.5, height = 5.5, units = "in", res = 1000)
  
  # Heatmap of sample_clean_fol and c_call (Isotype) and celltype
  # print(
  #   plot_df %>%
  #     filter(!c_call %in% c("IGHE", "") & !is.na(c_call)) %>%
  #     # filter(!is.na(c_call)) %>% 
  #     group_by(sample_clean_plot, c_call, L1_annotation) %>%
  #     summarise(Count = n(), .groups = "drop") %>%
  #     complete(sample_clean_plot, c_call, L1_annotation, 
  #              fill = list(Count = 0)) %>%
  #     ggplot(aes(x = L1_annotation, y = c_call, fill = Count)) +
  #     geom_tile(color = "white", linewidth = 0.3) +
  #     geom_text(aes(label = ifelse(Count > 0, Count, "")),
  #               color = "white", size = 2.5) +
  #     scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
  #                         limits = c(0, NA)) +
  #     facet_wrap(vars(sample_clean_plot), ncol = 3) +
  #     labs(
  #       x = "Cell type",
  #       y = "Ig Class",
  #       title = glue("{p}: Samples across samples, isotypes and cell types")
  #     ) +
  #     theme_minimal() +
  #     theme(
  #       axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
  #       axis.text.y = element_text(size = 9),
  #       panel.grid = element_blank(),
  #       strip.text = element_text(size = 9, face = "bold"),
  #       legend.position = "none",
  #       plot.title = element_text(face = "bold", size = 16)
  #     )
  # )
  
  print(
    plot_df %>%
      filter(!c_call %in% c("IGHE", "") & !is.na(c_call)) %>%
      filter(L1_annotation == "PCs") %>% 
      ggplot(aes(x = sample_clean_plot, fill = c_call)) +
      geom_bar() + 
      scale_fill_manual(values = isotype_colors_custom) + 
      theme_minimal()
      
      # 
      # geom_tile(color = "white", linewidth = 0.3) +
      # geom_text(aes(label = ifelse(Count > 0, Count, "")),
      #           color = "white", size = 2.5) +
      # scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
      #                     limits = c(0, NA)) +
      # facet_wrap(vars(sample_clean_plot), ncol = 3) +
      # labs(
      #   x = "Cell type",
      #   y = "Ig Class",
      #   title = glue("{p}: Samples across samples, isotypes and cell types")
      # ) +
      # theme_minimal() +
      # theme(
      #   axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      #   axis.text.y = element_text(size = 9),
      #   panel.grid = element_blank(),
      #   strip.text = element_text(size = 9, face = "bold"),
      #   legend.position = "none",
      #   plot.title = element_text(face = "bold", size = 16)
      # )
  )
  
  dev.off()
  
  # ggsave(glue("{outdir_5}/{HH}_isotypes_samples_cell_types.png"), width = 9.5, height = 5.5, dpi = 1000)
  
})

# ==============================================================================
# ISOTYPES - Upset plot
# ==============================================================================

# upset
lapply(patients, function(HH){
  
  # HH <- "HH119"
  p <- patient_names[[HH]]
  
  plot_df <- resolve_LC_list[[HH]] %>%
    filter(locus == "IGH") %>%
    mutate(
      sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-"))
    ) %>%
    filter(!c_call %in% c("IGHE", "") & !is.na(c_call))
  
  # Assuming plot_df has one row per sequence with c_call and L1_annotation columns
  upset_input_c <- plot_df %>%
    mutate(present = 1L) %>%
    pivot_wider(
      id_cols = sequence_id,
      names_from = c_call,
      values_from = present,
      values_fill = 0L
    ) %>%
    as.data.frame()
  
  upset_input_L1 <- plot_df %>%
    mutate(present = 1L) %>%
    pivot_wider(
      id_cols = sequence_id,
      names_from = L1_annotation,
      values_from = present,
      values_fill = 0L
    ) %>%
    as.data.frame()
  
  upset_input_sample <- plot_df %>%
    mutate(present = 1L) %>%
    pivot_wider(
      id_cols = sequence_id,
      names_from = sample_clean_plot, 
      values_from = present,
      values_fill = 0L
    ) %>%
    as.data.frame()
  
  upset_input_final <- full_join(upset_input_c, upset_input_L1, by = "sequence_id") %>% full_join(upset_input_sample, by = "sequence_id")
  
  # Define order
  c_order <- upset_input_c %>% column_to_rownames("sequence_id") %>% colSums() %>% sort() %>% names()
  L1_order <- upset_input_L1 %>% column_to_rownames("sequence_id") %>% colSums() %>% sort() %>% names()
  sample_order <- upset_input_sample %>% column_to_rownames("sequence_id") %>% colSums() %>% sort() %>% names()
  
  set_order <- c(L1_order, c_order, sample_order)
  
  set_colors <- c(
    rep("darkgreen", length(sample_order)),
    rep("darkred", length(c_order)),
    rep("blue", length(L1_order))
  )
  
  
  # Prep saving plot 
  png(glue("{outdir_5}/{HH}_upsetplot.png"), width = 13, height = 6, units = "in", res = 1000)
  
  # Plot 
  print(
    UpSetR::upset(
      upset_input_final,
      sets     = set_order,
      order.by = "freq",
      keep.order = TRUE, 
      sets.bar.color = set_colors
    )
  )
  
  # Title
  grid.text(
    glue("{p}: Sample, isotype and cell type intersections"),
    x = 0.65, y = 0.97,          # adjust position as needed
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  dev.off()
  
})

# ==============================================================================
# Frequency of top clone per follicle 
# ==============================================================================

outdir_6 <- glue("45_immcantation/{plot_version}/15_poster_figures/Freq_of_clones_across_follicles")
dir.create(outdir_6, recursive = TRUE)

n_clones <- 10

lapply(patients, function(HH){
  
  # HH <- "HH117"
  p <- patient_names[[HH]]
  
  # Subset to 5 clones
  top_GC_clones_subset <- top_GC_clones[[HH]][c(1:n_clones)]
  
  plot_df <- resolve_LC_list[[HH]] %>% 
    filter(locus == "IGH") %>% 
    mutate(
      sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-")),
      sample_clean_plot = fct_infreq(sample_clean_plot), #%>% fct_rev()
      clone_subgroup_id_plot = ifelse(clone_subgroup_id %in% top_GC_clones_subset, clone_subgroup_id, "other"),
      clone_subgroup_id_plot = factor(clone_subgroup_id_plot, levels = c(top_GC_clones_subset, "other"))
    ) %>%
    add_count(sample_clean_plot, name = "Count") 
  
  # Across follicles 
  HH_fol_sample_clean <- plot_df %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))
  
  # Define clone colors 
  clone_colors <- list("#E05C8A", "#66CC55", "#5588DD", "#EE9944", "#AA3377",
                       "#44BBAA", "#CC6644", "#4499CC", "#AACC33", "#9955BB",
                       "grey") %>% setNames(c(top_GC_clones_subset, "other"))
  
  # Define clone names
  clone_names <- c(paste("Clone", 1:n_clones), "Other") %>% as.list() %>% setNames(c(top_GC_clones_subset, "other"))
  
  # N clones with >= 3 cells
  N_clones_per_fol <- plot_df %>%
    filter(!is.na(manual_ADT_ID)) %>%
    mutate(
      manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
    ) %>% 
    group_by(manual_ADT_ID_plot) %>% 
    count(clone_subgroup_id) %>% 
    filter(n >= 3) %>% 
    count(manual_ADT_ID_plot) %>% 
    ungroup() %>%
    complete(
      manual_ADT_ID_plot = seq(min(manual_ADT_ID_plot), max(manual_ADT_ID_plot)),
      fill = list(n = 0)
    )

  png(glue("{outdir_6}/{HH}_N_{n_clones}.png"), width = 14, height = 7, units = "in", res = 1000)
  
  print(
    plot_df %>%
      filter(!is.na(manual_ADT_ID)) %>%
      mutate(
        manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
      ) %>% 
      ggplot(aes(x = manual_ADT_ID_plot)) + 
      geom_bar(aes(fill = clone_subgroup_id_plot), position = "fill") + 
      geom_text(
        data = N_clones_per_fol, 
        aes(x = manual_ADT_ID_plot, y = 1.02, label = n)
      ) +
      scale_fill_manual(
        values = clone_colors, 
        labels = clone_names
      ) + 
      scale_x_continuous(
        breaks = function(x) seq(1, ceiling(max(x)), by = 1),
        limits = c(0.5, NA),
        expand = c(0, 0.5)
      ) + 
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Frequency", 
        title = glue("{p}: Distribution of clones across {HH_fol_sample_clean} follicles"),
        subtitle = glue("Top {n_clones} clones highlighted and number of clones with >= 3 cells in each follicle is stated."),
        fill = "Clone"
      ) + 
      theme(plot.title = element_text(face = "bold", size = 16))
  )
  
  dev.off()
  
})


# ==============================================================================
# Frequency of top clone per follicle 
# ==============================================================================

outdir_7 <- glue("45_immcantation/{plot_version}/15_poster_figures/Freq_of_clones_across_samples")
dir.create(outdir_7, recursive = TRUE)

n_clones <- 10

lapply(patients, function(HH){
  
  # HH <- "HH117"
  p <- patient_names[[HH]]
  
  # Subset to 5 clones
  top_GC_clones_subset <- top_GC_clones[[HH]][c(1:n_clones)]
  
  plot_df <- resolve_LC_list[[HH]] %>% 
    filter(locus == "IGH") %>% 
    mutate(
      sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-")),
      sample_clean_plot = fct_infreq(sample_clean_plot), #%>% fct_rev()
      clone_subgroup_id_plot = ifelse(clone_subgroup_id %in% top_GC_clones_subset, clone_subgroup_id, "other"),
      clone_subgroup_id_plot = factor(clone_subgroup_id_plot, levels = c(top_GC_clones_subset, "other"))
    ) %>%
    add_count(sample_clean_plot, name = "Count") 
  
  # Across follicles 
  HH_fol_sample_clean <- plot_df %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))
  
  # Define clone colors 
  clone_colors <- list("#E05C8A", "#66CC55", "#5588DD", "#EE9944", "#AA3377",
                       "#44BBAA", "#CC6644", "#4499CC", "#AACC33", "#9955BB",
                       "grey") %>% setNames(c(top_GC_clones_subset, "other"))
  
  # Define clone names
  clone_names <- c(paste("Clone", 1:n_clones), "Other") %>% as.list() %>% setNames(c(top_GC_clones_subset, "other"))
  
  # N clones with >= 3 cells
  N_clones_per_fol <- plot_df %>% 
    group_by(sample_clean_plot) %>% 
    count(clone_subgroup_id) %>% 
    filter(n >= 3) %>% 
    count(sample_clean_plot) %>% 
    ungroup()
  
  png(glue("{outdir_7}/{HH}_N_{n_clones}.png"), width = 14, height = 7, units = "in", res = 1000)
  
  print(
    plot_df %>%
      ggplot(aes(x = sample_clean_plot)) + 
      geom_bar(aes(fill = clone_subgroup_id_plot), position = "fill") + 
      scale_fill_manual(
        values = clone_colors, 
        labels = clone_names
      ) + 
      geom_text(
        data = N_clones_per_fol, 
        aes(x = sample_clean_plot, y = 1.02, label = n)
      ) +
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Frequency", 
        title = glue("{p}: Distribution of clones across samples"),
        subtitle = glue("Top {n_clones} clones highlighted and number of clones with >= 3 cells in each sample is stated."),
        fill = "Clone"
      ) +
      theme(plot.title = element_text(face = "bold", size = 26))
  )
  
  dev.off()
  
})




