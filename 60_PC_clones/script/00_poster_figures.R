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
seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs_annotated.rds")

outdir <- glue("60_PC_clones/plot/00_poster_figures")
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

# Condition 
patient_to_condition <- seurat_integrated[[]] %>% 
  select(patient, condition) %>% 
  distinct() %>% 
  dplyr::rename(patient_id = patient) %>% 
  mutate(
    condition = ifelse(condition == "Crohn's", "CD", "CRC") %>% as_factor() %>% fct_relevel("CRC", "CD")
  )
rownames(patient_to_condition) <- NULL

df_both <- df_both %>% left_join(patient_to_condition, by = "patient_id")


# ==============================================================================
# Renaming 
# ==============================================================================

df_both$sample_clean %>% unique()

df_both <- df_both %>% 
  mutate(
    sample_clean = sample_clean %>% 
      str_replace("-SILP", "\nIleum") %>% 
      str_replace("-COLP", "\nColon") %>% 
      str_replace("-INF", "\nINF") %>% 
      str_replace("-nonINF", "\nnon-INF") 
  )

df_both$sample_clean %>% unique()

# ==============================================================================
# N cells
# ==============================================================================

outdir1 <- glue("{outdir}/N_cells/")
dir.create(outdir1, recursive = TRUE, showWarnings = FALSE)

seurat_meta <- seurat_integrated[[]] %>% 
  filter(
    L1_annotation == "PCs", 
    str_detect(sample_clean, "LP")
  ) %>% 
  mutate(
    sample_clean = sample_clean %>% 
      str_replace("-SILP", "\nIleum") %>% 
      str_replace("-COLP", "\nColon") %>% 
      str_replace("-INF", "\nINF") %>% 
      str_replace("-nonINF", "\nnon-INF"), 
    condition = ifelse(condition == "Crohn's", "CD", "CRC") %>% 
      as_factor() %>% 
      fct_relevel("CRC", "CD")
  )

seurat_meta %>% 
  count(condition, sample_clean, L1_annotation) %>% 
  ggplot(aes(x = sample_clean, y = n, fill = L1_annotation)) + 
  geom_col() + 
  scale_fill_manual(values = L1_colors) + 
  scale_y_continuous(
    breaks = scales::breaks_width(2000)
    # minor_breaks = scales::breaks_width(1000)
  ) + 
  theme_bw() + 
  labs(
    title = "Cell count across LP PC compartments", 
    y = "N cells",
    x = "Compartment"
  ) + 
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir1}/N_cells_per_sample.png"), width = 10, height = 6)

# ==============================================================================
# N cells with BCR available
# ==============================================================================

seurat_meta$cell_id <- glue("{seurat_meta$sample}_{rownames(seurat_meta)}") %>% str_remove("_\\d+")

table(df_both$cell_id %in% seurat_meta$cell_id)

seurat_meta <- seurat_meta %>% 
  left_join(df_both %>% select(cell_id, clone_subgroup_id_90_similarity, c_call_grouped), by = "cell_id") %>% 
  mutate(
    has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE)
  )

seurat_meta %>% 
  count(condition, sample_clean, has_bcr) %>% 
  ggplot(aes(x = sample_clean, y = n, fill = has_bcr)) + 
  geom_col() + 
  scale_fill_manual(values = list("grey", "purple")) + 
  scale_y_continuous(
    breaks = scales::breaks_width(2000)
    # minor_breaks = scales::breaks_width(1000)
  ) + 
  theme_bw() + 
  labs(
    title = "Cell count with BCR data available across LP PC compartments", 
    y = "N cells",
    x = "Compartment", 
    fill = "Has BCR data"
  ) + 
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir1}/N_cells_w_bcr_per_sample.png"), width = 10, height = 6)




# ==============================================================================
# N clones
# ==============================================================================

outdir4 <- glue("{outdir}/N_clones/")
dir.create(outdir4, recursive = TRUE, showWarnings = FALSE)

# N LP PCs
df_plot <- df_both %>% 
  filter(
    !is.na(clone_subgroup_id_90_similarity)
  ) 

df_plot %>% 
  select(condition, patient_id, sample_clean, clone_subgroup_id_90_similarity) %>% 
  distinct() %>% 
  count(condition, patient_id, sample_clean) %>% 
  ggplot(aes(x = sample_clean, y = n, fill = patient_id)) + 
  geom_col() + 
  geom_text(
    aes(label = n), size = 4, vjust = -0.5
  ) + 
  scale_fill_manual(values = patient_color_values) + 
  theme_bw() + 
  labs(
    title = "Clone count for LP PC samples", 
    y = "N clones",
    x = "Compartment", 
    fill = "Patient"
  ) + 
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir4}/N_clones_per_sample.png"), width = 10, height = 6)


# ==============================================================================
# Isotypes barplots
# ==============================================================================

outdir6 <- glue("{outdir}/isotypes/")
dir.create(outdir6, recursive = TRUE, showWarnings = FALSE)

df_plot <- df_both %>% 
  filter(
    !is.na(clone_subgroup_id_90_similarity) & !is.na(c_call_grouped)
  ) %>% 
  count(condition, sample_clean, c_call_grouped) 

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
  ) + 
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir6}/isotype_barplot.png"), width = 10, height = 6)

# ==============================================================================
# Clusters per clones
# ==============================================================================

outdir5 <- glue("{outdir}/PC_clusters/")
dir.create(outdir5, recursive = TRUE, showWarnings = FALSE)

PC_meta_GL <- readRDS("00_data/PC_meta_GL.rds")
PC_meta_GL$cell_id <- glue("{PC_meta_GL$sample}_{rownames(PC_meta_GL)}") %>% str_remove("_\\d+")

table(df_both$cell_id %in% PC_meta_GL$cell_id)

# top 20 clones per patient
top_clones <- df_both %>% 
  count(patient_id, clone_subgroup_id_90_similarity, name = "clone_size") %>% 
  group_by(patient_id) %>% 
  slice_max(clone_size, n = 20, with_ties = FALSE) %>% 
  mutate(clone_rank = row_number()) %>%      # 1 = largest
  ungroup()

df_plot <- df_both %>% 
  left_join(PC_meta_GL %>% select(cell_id, RNA_snn_res.0.4.merged), by = "cell_id") %>% 
  inner_join(top_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  filter(!is.na(RNA_snn_res.0.4.merged)) %>% 
  count(condition, patient_id, clone_rank, clone_size, RNA_snn_res.0.4.merged) %>% 
  mutate(patient_id_plot = glue("{condition}: {patient_id}"))

ggplot(df_plot, aes(x = factor(clone_rank), y = n, fill = RNA_snn_res.0.4.merged)) +
  geom_col(position = "fill") +
  facet_wrap(vars(patient_id_plot), scales = "free_x", ncol = 1) +
  scale_fill_manual(values = PC_clusters_colors) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  labs(
    x = "Clone (ranked by size within patient)",
    y = "Percentage of cells",
    fill = "Cluster",
    title = "Cluster composition of the 20 largest PC clones per patient"
  ) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir5}/PC_clusters_across_clones.png"), width = 7.5, height = 8)

 # ==============================================================================
# Isotypes barplots
# ==============================================================================

outdir5 <- glue("{outdir}/Isotype_per_clone/")
dir.create(outdir5, recursive = TRUE, showWarnings = FALSE)

# top 20 clones per patient
top_clones <- df_both %>% 
  count(patient_id, clone_subgroup_id_90_similarity, name = "clone_size") %>% 
  group_by(patient_id) %>% 
  slice_max(clone_size, n = 20, with_ties = FALSE) %>% 
  mutate(clone_rank = row_number()) %>%      # 1 = largest
  ungroup()

df_plot <- df_both %>% 
  filter(!is.na(c_call_grouped)) %>% 
  inner_join(top_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  count(condition, patient_id, clone_rank, clone_size, c_call_grouped) %>% 
  mutate(patient_id_plot = glue("{condition}: {patient_id}")) 

ggplot(df_plot, aes(x = factor(clone_rank), y = n, fill = c_call_grouped)) +
  geom_col(position = "fill") +
  facet_wrap(vars(patient_id_plot), scales = "free_x", ncol = 1) +
  scale_fill_manual(values = isotype_grouped_colors_custom) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  labs(
    x = "Clone (ranked by size within patient)",
    y = "Number of cells",
    fill = "Isotype",
    title = "Isotype composition of the 20 largest PC clones per patient"
  ) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir5}/Isotype_across_clones.png"), width = 7.5, height = 8)


# ------------------------------------------------------------------------------
# Circle-packing: one blob per follicle, circles = clones, sized by clone frequency
# ------------------------------------------------------------------------------

outdir7 <- glue("{outdir}/APackOfTheClones/")
dir.create(outdir7, recursive = TRUE, showWarnings = FALSE)

n_clones <- 15

library(packcircles)

lapply(patients, function(HH){
  
  # HH <- "HH151"
  
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
    slice_head(n = n_clones) %>% 
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
      subtitle = glue("Top {n_clones} clones shared across tissues, colored; all other clones in grey"),
      fill = "Clone"
    ) + 
    theme_void() + 
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 16, face = "bold")
    ) 
  
  ggsave(glue("{outdir7}/{HH}_circle_packing_shared_clones.png"), width = 10, height = 8)
  
})

# ==============================================================================
# DIVERSITY ANALYSIS
# ==============================================================================

library(alakazam)

outdir_diversity <- glue("{outdir}/diversity/")
dir.create(outdir_diversity, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Gini - both patients
# ------------------------------------------------------------------------------

gini_coeff <- function(clone_counts) {
  x <- sort(clone_counts)
  n <- length(x)
  numerator <- sum((2 * seq_along(x) - n - 1) * x)
  denominator <- (n - 1) * sum(x)
  numerator / denominator
}

gini_combined <- df_both %>%
  filter(!is.na(clone_subgroup_id_90_similarity)) %>%
  count(patient_id, sample_clean, clone_subgroup_id_90_similarity, name = "n_cells") %>%
  group_by(patient_id, sample_clean) %>%
  summarise(
    gini = gini_coeff(n_cells),
    total_clones = n_distinct(clone_subgroup_id_90_similarity),
    .groups = "drop"
  ) %>% 
  left_join(patient_to_condition, by = "patient_id")

ggplot(gini_combined, aes(x = sample_clean, y = gini, color = patient_id)) +
  # geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
  # geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
  geom_point(size = 4) + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = patient_color_values) +
  theme_bw() +
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
  labs(
    x = "Sample",
    y = "Gini coefficient",
    title = "Gini coefficient per LP PC compartment",
    color = "Patient"
    # subtitle = "Each point represents one compartment"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir_diversity}/gini_coef_combined.png"), width = 10, height = 6)


# ------------------------------------------------------------------------------
# Shannon - both patients
# ------------------------------------------------------------------------------

# Run alphaDiversity per patient and combine
# lapply(versions, function(version) {

for (min_n in c(30, 100)) {
  
  # min_n <- 30
  shannon_combined <- lapply(patients, function(HH) {
    
    df_sha <- df_both %>%
      filter(patient_id == HH)
    
    set.seed(123)
    curve <- alphaDiversity(
      df_sha,
      group = "sample_clean",
      clone = "clone_subgroup_id_90_similarity",
      min_q = 0, max_q = 4, step_q = 0.1,
      ci = 0.95, nboot = 100, min_n = min_n
    )
    
    curve@diversity %>%
      filter(q == 1) %>%
      mutate(patient_id = HH)
    
  }) %>% bind_rows()
  
  shannon_combined %>% 
    left_join(patient_to_condition, by = "patient_id") %>% 
    ggplot(aes(x = sample_clean, y = d, color = patient_id)) +
    # geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
    # geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
    geom_point(size = 4) +
    scale_color_manual(values = patient_color_values) +
    theme_bw() +
    facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
    labs(
      x = "Sample",
      y = "Shannon diversity (q=1)",
      title = glue("Shannon diversity per compartment - LP PCs"),
      subtitle = glue("Compartments with <{min_n} PCs excluded."),
      caption = "Bootstrapped estimates (100 resamples, 95% CI)", 
      color = "Patient"
    ) +
    theme(
      # legend.position = "none", 
      plot.title = element_text(size = 16, face = "bold"),
      axis.title   = element_text(size = 14),
      axis.text    = element_text(size = 12),
      strip.text   = element_text(size = 12, face = "bold")
    )
  
  ggsave(glue("{outdir_diversity}/shannon_diversity_min_n_{min_n}_combined.png"), width = 10, height = 6)
  
}

# ------------------------------------------------------------------------------
# Dxx plot - PC cell clones
# ------------------------------------------------------------------------------

compute_Dxx <- function(clone_counts, xx = 0.5) {
  sorted <- sort(clone_counts, decreasing = TRUE)
  cumulative <- cumsum(sorted) / sum(sorted)
  which(cumulative >= xx)[1]
}


d50_per_follicle <- df_both %>%
  filter(
    !is.na(clone_subgroup_id_90_similarity)
  ) %>%
  count(sample_clean, clone_subgroup_id_90_similarity, name = "n_cells") %>%
  group_by(sample_clean) %>%
  summarise(
    D50 = compute_Dxx(n_cells, 0.50),
    D20 = compute_Dxx(n_cells, 0.20),
    total_clones = n_distinct(clone_subgroup_id_90_similarity),
    .groups = "drop"
  ) %>%
  arrange(D50) %>% 
  mutate(
    patient_id = sample_clean %>% str_split_i(":", 1)
  ) %>% 
  left_join(patient_to_condition, by = "patient_id")

d50_per_follicle_long <- d50_per_follicle %>%
  mutate(
    D20_segment = D20,
    D50_extra = D50 - D20
  ) %>%
  pivot_longer(cols = c(D20_segment, D50_extra), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("D50_extra", "D20_segment")),
    patient_id = sample_clean %>% str_split_i("-", 1)
  )


d50_per_follicle_long %>% 
  filter(metric == "D20_segment") %>% 
  ggplot(aes(x = sample_clean, y = value, fill = metric)) +
  geom_col() +
  geom_text(
    data = d50_per_follicle,
    aes(x = sample_clean, y = D20, label = total_clones),
    inherit.aes = FALSE,
    vjust = -0.5, size = 4
  ) +
  scale_fill_manual(
    values = c("D20_segment" = "steelblue"),
    labels = c("D20_segment" = "D20")
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(20)
    # minor_breaks = scales::breaks_width(minor_breaks_width),
    # limits = c(0.5, NA)
  ) + 
  theme_bw() +
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
  labs(
    x = "Compartment",
    y = "Number of clones",
    fill = NULL,
    title = glue("Clonal dominance per LP PC compartment"),
    caption = "Numbers above bars indicate total LP PCs clone count per compartment"
  ) + 
  theme(
    # legend.position = "none", 
    plot.title = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold")
  )

ggsave(glue("{outdir_diversity}/clonal_D20_plot.png"), width = 10, height = 6)

# ==============================================================================
# Clonal sharing
# ==============================================================================


outdir_clonal_sharing <- glue("{outdir}/clonal_sharing/")
dir.create(outdir_clonal_sharing, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# For each clone, count presents in compartents per patient
# ------------------------------------------------------------------------------

for (n_min_cells in c(1, 2)){
  
  # Min number of cells in follicle to be sure of presents 
  # n_min_cells <- 2
  
  # Clones are present in both compartments  
  df_n_compartment <- df_both %>% 
    group_by(patient_id, clone_subgroup_id_90_similarity) %>% 
    count(sample_clean) %>% 
    filter(n >= n_min_cells) %>%
    count(sample_clean) %>% 
    count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
    dplyr::rename(n_compartments = n) %>% 
    ungroup()
  
  df_n_compartment %>% nrow()
  
  # N cells per clone 
  df_N_cells <- df_both %>% 
    right_join(
      df_n_compartment, by = c("patient_id", "clone_subgroup_id_90_similarity")
    ) %>% 
    count(clone_subgroup_id_90_similarity, patient_id) %>% 
    dplyr::rename(clone_size = n)
  
  # combine
  nrow(df_n_compartment) == nrow(df_N_cells)
  df_plot <- df_n_compartment %>% 
    left_join(df_N_cells, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
    mutate(
      clone_size_group = case_when(
        clone_size > 1 & clone_size <= 5 ~ "2-5", 
        clone_size > 5 & clone_size <= 10 ~ "6-10",
        clone_size > 10 & clone_size <= 20 ~ "11-20",
        clone_size > 20 & clone_size <= 50 ~ "21-50",
        clone_size > 50 & clone_size <= 100 ~ "51-100",
        clone_size > 100 ~ "100+"
      ),
      clone_size_group = factor(clone_size_group, levels = c("2-5", "6-10", "11-20", "21-50", "51-100", "100+"))#,
      # is_largest_clone = ifelse(patient_id == "HH119" & clone_subgroup_id_90_similarity == large_crc_clone, TRUE, FALSE) 
    ) %>% 
    group_by(n_compartments, patient_id) %>% 
    mutate(
      mean_clone_size = mean(clone_size) %>% round(1),
      median_clone_size = median(clone_size)
    ) %>% 
    ungroup() %>% 
    left_join(patient_to_condition, by = "patient_id")
  
  # N 
  df_N <- df_plot %>%
    select(-clone_subgroup_id_90_similarity) %>% 
    count(condition, patient_id, n_compartments, mean_clone_size, median_clone_size, sort = TRUE) 
  
  # Jitter plot split by patient 
  df_plot %>% 
    # mutate(
    #   patient_id_plot = case_when(
    #     patient_id == "HH119" ~ "HH119 (CRC)",
    #     patient_id == "HH117" ~ "HH117 (CD)",
    #     patient_id == "HH151" ~ "HH151 (CD)",
    #     patient_id == "HH153" ~ "HH153 (CD)"
      # )
    # ) %>% 
    ggplot(aes(x = patient_id, y = n_compartments)) + 
    geom_jitter(
      aes(color = clone_size_group, size = clone_size_group), 
      alpha = 0.5, width = 0.30, height = 0.2#, size = 2.5
    ) + 
    scale_color_viridis_d(option = "plasma", direction = -1) +
    # geom_text(
    #   data = df_N,
    #   aes(x = patient_id, y = n_compartments, label = glue("{n} clones ({mean_clone_size}; {median_clone_size})")),
    #   size = 3,
    #   position = position_nudge(x = 0.45)
    # ) +
    theme_bw() + 
    scale_y_continuous(
      breaks = scales::breaks_width(1),
      minor_breaks = scales::breaks_width(1)
    ) + 
    facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
    labs(
      x = "Patient",
      y = "N compartments",
      title = "Clonal sharing across LP PC compartments",
      subtitle = glue("For each clone, how many compartments is it present in (at least {n_min_cells} cells)"),
      caption = "N clones (mean clone size; median clone size)", 
      size = "Clone size",
      color = "Clone size"
    ) + 
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title   = element_text(size = 14),
      axis.text    = element_text(size = 12),
      strip.text   = element_text(size = 12, face = "bold")
    )
  
  
  ggsave(glue("{outdir_clonal_sharing}/clones_in_n_compartments_{n_min_cells}.png"), width = 9, height = 5)
  
  # ------------------------------------------------------------------------------
  # Clonal sharing proportions 
  # ------------------------------------------------------------------------------
  
  # These plots (above and this one under) only include clones/cells that are present 
  # in the abundance of n_min_cells. This is both the number and the plots. 
  
  df_totals <- df_plot %>% 
    group_by(condition, patient_id) %>% 
    summarise(
      n_clones = n(),
      total_cells = sum(clone_size),
      .groups = "drop"
    )
  
  df_plot %>% 
    mutate(n_follicle_fct = n_compartments %>% as.factor()) %>% 
    ggplot(aes(x = patient_id, fill = n_follicle_fct)) + 
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    geom_bar(position = "fill") + 
    geom_text(
      data = df_totals,
      aes(x = patient_id, y = 0.9, label = glue("{n_clones} clones\n({total_cells} cells)")),
      inherit.aes = FALSE
      # size = 3
    ) +
    scale_y_continuous(labels = scales::percent) +
    facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") + 
    theme_bw() +
    labs(
      fill = "N compartments",
      x = "Patient ID",
      y = "% of clones",
      title = "Percentage of clonal sharing across LP PC compartments",
      subtitle = glue("For each clone, how many compartments is it present in (at least {n_min_cells} cells)")
    ) + 
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title   = element_text(size = 14),
      axis.text    = element_text(size = 12),
      strip.text   = element_text(size = 12, face = "bold")
    )
  
  
  ggsave(glue("{outdir_clonal_sharing}/clones_in_n_compartments_proportions_{n_min_cells}.png"), width = 10, height = 6.5)
  
}
 
