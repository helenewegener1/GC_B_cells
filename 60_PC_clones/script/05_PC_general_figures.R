library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)
patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

# Load data and filter for LP PCs
resolve_LC_list <- lapply(resolve_LC_files, function(x){
  readRDS(glue("45_immcantation/out/rds/{x}")) %>% 
    filter(locus == "IGH" & L1_annotation == "PCs" & str_detect(sample_clean, "LP"))
}) %>% 
  setNames(patients)

# Look at clone IDs
grep("clone", colnames(resolve_LC_list$HH117), value = TRUE)

# clone_subgroup_id_90_similarity

# resolve_LC <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds"))
# table(resolve_LC$locus)
# 
# df_heavy <- resolve_LC %>% filter(locus == "IGH")
# 
# nrow(df_heavy)

# Load seurat object
seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

outdir <- glue("60_PC_clones/plot/05_general_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

resolve_LC_list <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  resolve_LC_list[[HH]] %>%
    group_by(clone_subgroup_id_90_similarity, locus) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)

df_both <- bind_rows(resolve_LC_list)

# ==============================================================================
# N clones
# ==============================================================================

outdir4 <- glue("{outdir}/N_clones/")
dir.create(outdir4, recursive = TRUE, showWarnings = FALSE)

# N GC B cells in PPs
df_plot <- df_both %>% 
  filter(
    !is.na(clone_subgroup_id_90_similarity)
  ) 

df_plot %>% 
  select(sample_clean, clone_subgroup_id_90_similarity) %>% 
  distinct() %>% 
  count(sample_clean) %>% 
  ggplot(aes(x = sample_clean, y = n)) + 
  geom_col() + 
  geom_text(
    aes(label = n), size = 3, vjust = -0.5
  ) + 
  theme_bw() + 
  labs(
    title = "N clones for LP PC samples", 
    y = "N clones",
    x = "Compartment"
  )

ggsave(glue("{outdir4}/N_clones_per_sample.png"))

# ==============================================================================
# Clone size graph (Freds plot)
# ==============================================================================

outdir5 <- glue("{outdir}/clone_size/")
dir.create(outdir5, recursive = TRUE, showWarnings = FALSE)

df_plot <- df_both %>% 
  filter(
    !is.na(clone_subgroup_id_90_similarity)
  ) %>% 
  count(sample_clean, clone_subgroup_id_90_similarity) %>% 
  dplyr::rename(clone_size = n) %>% 
  mutate(
    clone_size_group = case_when(
      clone_size == 1 ~ "Singleton",
      # clone_size > 1 & clone_size <= 5 ~ "2-5", 
      clone_size == 2 ~ "2",
      clone_size == 3 ~ "3",
      clone_size == 4 ~ "4",
      clone_size == 5 ~ "5",
      clone_size > 5 & clone_size <= 10 ~ "6-10",
      clone_size > 10 & clone_size <= 20 ~ "11-20",
      clone_size > 20 & clone_size <= 50 ~ "21-50",
      clone_size > 50 & clone_size <= 100 ~ "51-100",
      clone_size > 100 ~ "100+"
    ),
    # clone_size_group = factor(clone_size_group, levels = c("Singleton", "2-5", "6-10", "11-20", "21-50", "51-100", "100+"))
    clone_size_group = factor(clone_size_group, levels = c("2", "3", "4", "5", "6-10", "11-20", "21-50", "51-100", "100+"))
  ) %>% 
  count(sample_clean, clone_size_group)

df_plot %>% 
  filter(clone_size_group != "Singleton") %>%
  ggplot(aes(x = clone_size_group, y = n, color = sample_clean)) + 
  geom_point(alpha = 0.5, size = 2) + 
  geom_line(aes(group = sample_clean)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(
    title = "Frequency of clone sizes of LP PCs",
    x = "Clone size (N cells)", 
    y = "N clones", 
    color = "Compartment"
  )

ggsave(glue("{outdir5}/freq_of_clones_size.png"), width = 10)

# ==============================================================================
# Isotypes barplots
# ==============================================================================

outdir6 <- glue("{outdir}/isotypes/")
dir.create(outdir6, recursive = TRUE, showWarnings = FALSE)

df_plot <- df_both %>% 
  filter(
    !is.na(clone_subgroup_id_90_similarity) & !is.na(c_call_grouped)
  ) %>% 
  count(sample_clean, c_call_grouped) 

df_plot %>% 
  ggplot(aes(x = sample_clean, y = n, fill = c_call_grouped)) + 
  geom_col() + 
  scale_fill_manual(values = isotype_grouped_colors_custom) + 
  theme_bw() + 
  labs(
    title = "Isotypes of LP PCs clones",
    x = "Samples", 
    y = "N cells", 
    fill = "Isotype"
  )
  
ggsave(glue("{outdir6}/isotype_barplot.png"))

# ==============================================================================
# APackOfTheClones
# ==============================================================================

library(Seurat)
library(scRepertoire)
library(APackOfTheClones)

outdir7 <- glue("{outdir}/APackOfTheClones/")
dir.create(outdir7, recursive = TRUE, showWarnings = FALSE)

lapply(patients, function(HH){
  
  # HH <- "HH117"
  
  # Subset seurat object to LP PCs
  seurat_integrated_PC <- subset(seurat_integrated, patient == HH & str_detect(sample, "LP") & L1_annotation == "PCs")
  df_HH <- df_both %>% filter(patient_id == HH)
  
  rownames(seurat_integrated_PC[[]]) %>% head()
  df_HH$cell_id %>% head()
  
  # Update cell_id
  seurat_integrated_PC[[]] <- seurat_integrated_PC[[]] %>% 
    mutate(
      cell_id = glue("{sample}_{rownames(.)}") %>% str_remove("_\\d")
    )
  
  seurat_integrated_PC$cell_id %>% head()
  df_HH$cell_id %>% head()
  
  table(seurat_integrated_PC$cell_id %in% df_HH$cell_id)
  
  # Subset to cells in df_both 
  seurat_integrated_PC <- subset(seurat_integrated_PC, cell_id %in% df_HH$cell_id)
  
  table(seurat_integrated_PC$cell_id %in% df_HH$cell_id)
  
  # Add clone information to seurat object 
  seurat_integrated_PC[[]] <- seurat_integrated_PC[[]] %>% left_join(df_HH, by = "cell_id")
  
  # APackOfTheClones
  Idents(seurat_integrated_PC) <- "sample"
  vizAPOTC(
    seurat_integrated_PC,
    reduction_base = "umap.unintegrated",
    clonecall = "clone_subgroup_id_90_similarity", 
    show_labels = TRUE, 
    legend_text_size = 3,
    label_size = 3
  ) + 
    plot_annotation(title = glue("{HH}: LP PCs clone sizes"))
  
  ggsave(glue("{outdir7}/{HH}_APackOfTheClones.png"), width = 7, height = 6.5)
  
})


# ------------------------------------------------------------------------------
# Circle-packing: one blob per follicle, circles = clones, sized by clone frequency
# ------------------------------------------------------------------------------

outdir7 <- glue("{outdir}/APackOfTheClones/")
dir.create(outdir7, recursive = TRUE, showWarnings = FALSE)

library(packcircles)

lapply(patients, function(HH){
  
  # HH <- "HH117"
  
  # Subset data to patient, PPs and GC B cells
  df_HH <- df_both %>% filter(patient_id == HH)
  
  # ------------------------------------------------------------------------------
  # Identify top 15 shared clones (present in >1 follicle), by total cell count
  # ------------------------------------------------------------------------------
  
  clone_counts <- df_HH %>% 
    count(sample_clean, clone_subgroup_id_90_similarity, name = "clone_size") %>% 
    arrange(sample_clean, desc(clone_size))
  
  top_shared_clones <- clone_counts %>% 
    group_by(clone_subgroup_id_90_similarity) %>% 
    summarise(n_follicles = n_distinct(sample_clean), total_cells = sum(clone_size), .groups = "drop") %>% 
    filter(n_follicles > 1) %>% 
    arrange(desc(total_cells)) %>% 
    slice_head(n = 15) %>% 
    pull(clone_subgroup_id_90_similarity)
  
  # color: one distinct color per top shared clone, grey for everything else
  top_clone_colors <- set_names(
    scales::hue_pal()(length(top_shared_clones)),
    top_shared_clones
  )
  clone_colors <- c(top_clone_colors, "Other" = "grey80")
  
  clone_counts <- clone_counts %>% 
    mutate(
      clone_color_group = if_else(
        clone_subgroup_id_90_similarity %in% top_shared_clones, 
        clone_subgroup_id_90_similarity, 
        "Other"
      )
    )
  
  # ------------------------------------------------------------------------------
  # Circle packing (same as before, now carrying clone_color_group through)
  # ------------------------------------------------------------------------------
  
  all_circles <- data.frame()
  samples <- unique(clone_counts$sample_clean)
  
  for (sample in samples) {
    
    # sample <- samples[[1]]
    
    df_sample <- clone_counts %>% filter(sample_clean == sample)
    
    inner_layout <- circleProgressiveLayout(sqrt(df_sample$clone_size), sizetype = "radius") %>% 
      mutate(
        sample_clean = sample,
        clone_subgroup_id_90_similarity = df_sample$clone_subgroup_id_90_similarity,
        clone_size = df_sample$clone_size,
        clone_color_group = df_sample$clone_color_group,
        radius = sqrt(df_sample$clone_size)
      )
    
    bounding_radius <- max(sqrt(inner_layout$x^2 + inner_layout$y^2) + inner_layout$radius)
    inner_layout$bounding_radius <- bounding_radius
    
    all_circles <- bind_rows(all_circles, inner_layout)
    
  }
  
  fol_radii <- all_circles %>% distinct(sample_clean, bounding_radius)
  outer_layout <- circleProgressiveLayout(fol_radii$bounding_radius, sizetype = "radius") %>% 
    mutate(sample_clean = fol_radii$sample_clean) %>% 
    select(sample_clean, x_offset = x, y_offset = y, blob_radius = radius)
  
  plot_data <- all_circles %>% 
    left_join(outer_layout %>% select(sample_clean, x_offset, y_offset), by = "sample_clean") %>% 
    mutate(x_final = x + x_offset, y_final = y + y_offset)
  
  plot_circles_df <- circleLayoutVertices(
    data.frame(x = plot_data$x_final, y = plot_data$y_final, radius = plot_data$radius),
    npoints = 50
  ) %>% 
    mutate(
      sample_clean = plot_data$sample_clean[id],
      clone_size = plot_data$clone_size[id],
      clone_color_group = plot_data$clone_color_group[id]
    )
  
  label_data <- outer_layout %>% 
    mutate(label_x = x_offset, label_y = y_offset)
  
  # ------------------------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------------------------
  
  ggplot() + 
    geom_polygon(
      data = plot_circles_df, 
      aes(x = x, y = y, group = id, fill = clone_color_group), 
      color = "white", linewidth = 0.2
    ) + 
    geom_text(
      data = label_data,
      aes(x = label_x, y = label_y, label = sample_clean),
      size = 4, fontface = "bold", color = "black"
    ) + 
    scale_fill_manual(values = clone_colors) + 
    coord_equal() + 
    labs(
      title = glue("{HH}: LP PCs clones"), 
      subtitle = "Top 15 clones shared across tissues, colored; all other clones in grey",
      fill = "Clone"
    ) + 
    theme_void() + 
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(glue("{outdir7}/{HH}_circle_packing_shared_clones.png"), width = 10, height = 8)
  
})



