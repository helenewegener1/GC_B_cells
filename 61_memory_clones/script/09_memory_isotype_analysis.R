library(glue)
library(tidyverse)
library(alakazam)
library(scatterpie)
library(patchwork)
library(ggbreak)

source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

# Prep output
outdir = glue("61_memory_clones/plot/09_isotype_analysis/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
    filter(
      locus == "IGH" & L1_annotation == "Memory_Bcells"
    )
}) %>% bind_rows()

# ==============================================================================
# Are different isotypes of Plasma cells in LP clonally overlapping or not.  
# % of smaller set and MHI. Also differences in clonal diversity and size between isotypes? 
# ==============================================================================

# ------------------------------------------------------------------------------
# N isotypes per clone per patient / tissue 
# ------------------------------------------------------------------------------

for (n_min_cells in c(1, 2)){
  
  # Min number of cells in follicle to be sure of presents 
  n_min_cells <- 2
  
  # PER PATIENT
  
  # Define N isotypes per clone 
  df_n_isotypes <- df_both %>% 
    count(patient_id, clone_subgroup_id_90_similarity, c_call_grouped) %>% 
    filter(n >= n_min_cells) %>%
    count(patient_id, clone_subgroup_id_90_similarity) %>% 
    dplyr::rename(n_isotypes = n) %>% 
    mutate(
      n_isotypes = as.factor(n_isotypes)
    )
    
  # Plot
  df_n_isotypes %>% 
    ggplot(aes(x = patient_id, group = n_isotypes, fill = n_isotypes)) + 
    geom_bar() + 
    # geom_text(
    #   aes(label = n_isotypes, y = n_isotypes+100)
    # ) + 
    theme_bw() + 
    scale_fill_viridis_d(option = "magma", direction = -1, na.value = "grey90") +  
    labs(
      title = "N isotypes per clone per patient - Memory B cells",
      subtitle = glue("For each clone, how many isotypes are found (at least {n_min_cells} cells)"),
      x = "Patient ID", 
      y = "N clones", 
      fill = "N isotypes"
    )
    
  ggsave(glue("{outdir}/isotype_per_clone_per_patient_{n_min_cells}.png"), width = 10, height = 6.5)
  
  # PER TISSUE
  
  # Define N isotypes per clone 
  df_n_isotypes <- df_both %>% 
    count(sample_clean, clone_subgroup_id_90_similarity, c_call_grouped) %>% 
    filter(n >= n_min_cells) %>%
    count(sample_clean, clone_subgroup_id_90_similarity) %>% 
    dplyr::rename(n_isotypes = n) %>% 
    mutate(
      n_isotypes = as.factor(n_isotypes)
    )
  
  # Plot
  df_n_isotypes %>% 
    ggplot(aes(x = sample_clean, group = n_isotypes, fill = n_isotypes)) + 
    geom_bar() + 
    # geom_text(
    #   aes(label = n_isotypes, y = n_isotypes+100)
    # ) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_viridis_d(option = "magma", direction = -1, na.value = "grey90") +  
    labs(
      title = "N isotypes per clone per sample_clean - Memory B cells",
      subtitle = glue("For each clone, how many isotypes are found (at least {n_min_cells} cells)"),
      x = "Sample", 
      y = "N clones", 
      fill = "N isotypes"
    )
  
  ggsave(glue("{outdir}/isotype_per_clone_per_sample_clean_{n_min_cells}.png"), width = 10, height = 6.5)
    
}  
  
# ------------------------------------------------------------------------------
# Clonal overlap between isotypes: raw shared clones, % of smaller set, and MHI
# ------------------------------------------------------------------------------

n_min_cells <- 2  # min cells for an isotype to count as present in a clone

clone_isotype_counts <- df_both %>% 
  filter(!is.na(c_call_grouped)) %>% 
  count(patient_id, c_call_grouped, clone_subgroup_id_90_similarity, name = "n_cells") %>% 
  filter(n_cells >= n_min_cells)

df_overlap <- data.frame()

for (HH in unique(clone_isotype_counts$patient_id)) {
  
  # HH <- "HH119"
  
  df_HH <- clone_isotype_counts %>% filter(patient_id == HH)
  isotypes_HH <- unique(df_HH$c_call_grouped)
  isotype_pairs <- combn(isotypes_HH, 2, simplify = FALSE)
  
  for (pair in isotype_pairs) {
    
    # pair <- isotype_pairs[[1]]
    
    iso_a <- pair[1]
    iso_b <- pair[2]
    
    df_a <- df_HH %>% filter(c_call_grouped == iso_a)
    df_b <- df_HH %>% filter(c_call_grouped == iso_b)
    
    clones_a <- df_a$clone_subgroup_id_90_similarity
    clones_b <- df_b$clone_subgroup_id_90_similarity
    
    n_a <- set_names(df_a$n_cells, clones_a)
    n_b <- set_names(df_b$n_cells, clones_b)
    
    # raw N clones shared between the two isotypes
    shared_clones <- intersect(clones_a, clones_b)
    n_shared_clones <- length(shared_clones)
    
    # % overlap relative to the smaller set (size-normalized, presence/absence)
    pct_overlap <- 100 * n_shared_clones / min(length(clones_a), length(clones_b))
    
    # Morisita-Horn index (abundance-based, via proportions)
    all_clones <- union(clones_a, clones_b)
    p <- n_a[all_clones]; p[is.na(p)] <- 0; p <- p / sum(p)
    q <- n_b[all_clones]; q[is.na(q)] <- 0; q <- q / sum(q)
    mhi <- 2 * sum(p * q) / (sum(p^2) + sum(q^2))
    
    # Prep row for current isotype pair
    df_overlap_pair <- data.frame(
      patient_id = HH,
      isotype_a = iso_a,
      isotype_b = iso_b,
      n_clones_a = length(clones_a),
      n_clones_b = length(clones_b),
      n_shared_clones = n_shared_clones,
      pct_overlap_smaller_set = pct_overlap,
      morisita_horn = mhi
    )
    
    # Append to dataframe
    df_overlap <- bind_rows(df_overlap, df_overlap_pair)
    
  }
}

df_overlap

# fix isotype_a/isotype_b to a consistent factor order so the upper-triangle logic below works
df_overlap_plot <- df_overlap %>% 
  mutate(
    isotype_a_plot = glue("{isotype_a}\n({n_clones_a})"),
    isotype_b_plot = glue("{isotype_b}\n({n_clones_b})")
  )

isotype_order_plot <- sort(unique(c(as.character(df_overlap_plot$isotype_a_plot), as.character(df_overlap_plot$isotype_b_plot))))

df_overlap_plot <- df_overlap_plot %>% 
  mutate(
    isotype_a_plot = factor(isotype_a_plot, levels = isotype_order_plot),
    isotype_b_plot = factor(isotype_b_plot, levels = isotype_order_plot)
  ) %>% 
  filter()

# heatmap: raw N shared clones between isotype pairs 
df_overlap_plot %>% 
  ggplot(aes(x = isotype_a_plot, y = isotype_b_plot, fill = n_shared_clones)) + 
  geom_tile(color = "white") + 
  geom_text(aes(label = n_shared_clones), size = 3.5) + 
  scale_fill_gradient(low = "white", high = "darkorange") + 
  scale_x_discrete(limits = isotype_order_plot) + 
  scale_y_discrete(limits = rev(isotype_order_plot)) + 
  facet_wrap(~patient_id) + 
  labs(
    x = NULL, y = NULL, fill = "N shared\nclones",
    title = "Clonal overlap between isotypes of memory B cell clones",
    subtitle = glue("Number of clones containing both isotypes (at least {n_min_cells} cells each)"),
    caption = "Number in parenthesis specifies the number of clones with given isotype"
  ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(glue("{outdir}/mem_isotype_overlap_N.png"), width = 13, height = 7)

# heatmap: % overlap of smaller set
df_overlap_plot %>%
  ggplot(aes(x = isotype_a_plot, y = isotype_b_plot, fill = pct_overlap_smaller_set)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(pct_overlap_smaller_set, 1)), size = 3.5) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) +
  scale_x_discrete(limits = isotype_order_plot) +
  scale_y_discrete(limits = rev(isotype_order_plot)) +
  facet_wrap(~patient_id) +
  labs(
    x = NULL, y = NULL, fill = "% overlap\n(smaller set)",
    title = "Clonal overlap between isotypes of memory B cell clones",
    subtitle = glue("Shared clones as a percentage of the smaller isotype's clone count (at least {n_min_cells} cells each)"),
    caption = "Number in parenthesis specifies the number of clones with given isotype"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(glue("{outdir}/mem_isotype_overlap_pct_smaller_set.png"), width = 13, height = 7)

# heatmap: Morisita-Horn index (upper triangle only)
df_overlap_plot %>% 
  ggplot(aes(x = isotype_a_plot, y = isotype_b_plot, fill = morisita_horn)) + 
  geom_tile(color = "white") + 
  geom_text(aes(label = round(morisita_horn, 2)), size = 3.5) + 
  scale_fill_gradient(low = "white", high = "forestgreen", limits = c(0, 1)) + 
  scale_x_discrete(limits = isotype_order_plot) + 
  scale_y_discrete(limits = rev(isotype_order_plot)) + 
  facet_wrap(~patient_id) + 
  labs(
    x = NULL, y = NULL, fill = "MHI",
    title = "Clonal overlap between PC isotypes of memory B cell clones",
    subtitle = glue("Morisita-Horn index, weighting shared clones by relative cell abundance (at least {n_min_cells} cells each)"),
    caption = "Number in parenthesis specifies the number of clones with given isotype"
  ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(glue("{outdir}/mem_isotype_overlap_MHI.png"), width = 13, height = 7)

# ------------------------------------------------------------------------------
# Clonal diversity per isotype (Hill numbers, alakazam)
# ------------------------------------------------------------------------------

# diversity_curve <- alphaDiversity(
#   df_both, 
#   group = "c_call_grouped", 
#   clone = "clone_subgroup_id_90_similarity",
#   min_n = 30,     # rarefies all isotypes to the smallest group's cell count for fair comparison
#   step_q = 0.25, 
#   nboot = 200
# )
# 
# plot(diversity_curve, colors = isotype_colors_custom) + 
#   labs(
#     title = "PC clonal diversity by isotype in LP",
#     subtitle = "Hill diversity curve, rarefied to the smallest isotype group's cell count"
#   )
# 
# ggsave(glue("{outdir}/pc_isotype_diversity_curve.png"), width = 8, height = 6)
# 
# # per patient too
# diversity_curve_patient <- alphaDiversity(
#   df_both,
#   group = c("patient_id", "c_call_grouped"),
#   clone = "clone_subgroup_id_90_similarity",
#   min_n = 30,
#   step_q = 0.25,
#   nboot = 200
# )
# 
# plot(diversity_curve_patient, colors = isotype_colors_custom) + 
#   facet_wrap(~patient_id) +
#   labs(
#     title = "PC clonal diversity by isotype in LP, per patient",
#     subtitle = "Hill diversity curve, rarefied to the smallest isotype group's cell count"
#   )
# 
# ggsave(glue("{outdir}/pc_isotype_diversity_curve_per_patient.png"), width = 10, height = 6)
# 
# # extract Hill numbers at q = 0 (richness), q = 1 (Shannon exp), q = 2 (Simpson-like) as a table
# diversity_curve@diversity %>% 
#   filter(round(q, 2) %in% c(0, 1, 2)) %>% 
#   select(c_call_grouped, q, d, d_ci_lower = d_lower, d_ci_upper = d_upper)


# ------------------------------------------------------------------------------
# Clone size differences between isotypes
# ------------------------------------------------------------------------------

clone_sizes <- df_both %>% 
  filter(!is.na(c_call_grouped)) %>% 
  count(patient_id, c_call_grouped, clone_subgroup_id_90_similarity, name = "clone_size")

clone_sizes %>% 
  ggplot(aes(x = c_call_grouped, y = clone_size, fill = c_call_grouped)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.7) + 
  scale_y_log10() + 
  facet_wrap(~patient_id) + 
  scale_fill_manual(values = isotype_colors_custom) + 
  labs(
    x = "Isotype", y = "Clone size (N cells, log scale)",
    title = "Clone size by isotype of memory B cell clones",
    subtitle = "Each point is one clone."
  ) + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(glue("{outdir}/mem_isotype_clone_size.png"), width = 10, height = 6)

# Kruskal-Wallis test per patient (non-parametric, since clone sizes are heavily right-skewed)
clone_sizes %>% 
  group_by(patient_id) %>% 
  summarise(kruskal_p = kruskal.test(clone_size ~ c_call_grouped)$p.value)


# # ==============================================================================
# # Within same isotype what is clonal overlap between inflamed and uninflamed (CD). 
# # For CRC ileum vs colon. % of smaller set and MHI. 
# # ==============================================================================
# 
# # ------------------------------------------------------------------------------
# # LP PCs: within-isotype clonal overlap between conditions
# # (HH119 CD: inflamed vs uninflamed; HH117 CRC: ileum vs colon)
# # ------------------------------------------------------------------------------
# 
# table(df_both$patient_id, df_both$sample_clean)
# 
# n_min_cells <- 2  # min cells for an isotype-condition combo to count as present in a clone
# 
# lapply(patients, function(HH){
#   
#   # HH <- "HH117"
#   
#   # Add condition 
#   df_condition <- df_both %>% 
#     filter(
#       patient_id == HH & !is.na(c_call_grouped)
#     ) %>% 
#     mutate(
#       condition = sample_clean %>% str_remove(glue("{HH}-"))
#     )
#   
#   df_condition$condition %>% table()
#   
#   clone_condition_counts <- df_condition %>% 
#     count(patient_id, c_call_grouped, condition, clone_subgroup_id_90_similarity, name = "n_cells") %>% 
#     filter(n_cells >= n_min_cells)
#   
#   conditions_HH <- unique(clone_condition_counts$condition)
#   isotypes_HH <- unique(clone_condition_counts$c_call_grouped)
#     
#   # skip patients without exactly 2 conditions defined
#   # if (length(conditions_HH) != 2) next
#   
#   cond_a <- conditions_HH[1]
#   cond_b <- conditions_HH[2]
#   
#   df_overlap <- data.frame()
#   
#   for (iso in isotypes_HH) {
#     
#     # iso <- isotypes_HH[1]
#     
#     # Zoom in on one isotype
#     df_iso <- clone_condition_counts %>% filter(c_call_grouped == iso)
#     
#     df_a <- df_iso %>% filter(condition == cond_a)
#     df_b <- df_iso %>% filter(condition == cond_b)
#     
#     clones_a <- df_a$clone_subgroup_id_90_similarity
#     clones_b <- df_b$clone_subgroup_id_90_similarity
#     
#     # skip if either condition has no clones for this isotype
#     if (length(clones_a) == 0 | length(clones_b) == 0) next
#     
#     # Clone sizes
#     n_a <- set_names(df_a$n_cells, clones_a)
#     n_b <- set_names(df_b$n_cells, clones_b)
#     
#     # raw N clones shared between the two conditions
#     shared_clones <- intersect(clones_a, clones_b)
#     n_shared_clones <- length(shared_clones)
#     
#     # % overlap relative to the smaller set
#     pct_overlap <- 100 * n_shared_clones / min(length(clones_a), length(clones_b))
#     
#     # Morisita-Horn index
#     all_clones <- union(clones_a, clones_b)
#     p <- n_a[all_clones]; p[is.na(p)] <- 0; p <- p / sum(p)
#     q <- n_b[all_clones]; q[is.na(q)] <- 0; q <- q / sum(q)
#     mhi <- 2 * sum(p * q) / (sum(p^2) + sum(q^2))
#     
#     # Prep row for current isotype
#     df_overlap_iso <- data.frame(
#       patient_id = HH,
#       isotype = iso,
#       condition_a = cond_a,
#       condition_b = cond_b,
#       n_clones_a = length(clones_a),
#       n_clones_b = length(clones_b),
#       n_shared_clones = n_shared_clones,
#       pct_overlap_smaller_set = pct_overlap,
#       morisita_horn = mhi
#     )
#     
#     # Append to dataframe
#     df_overlap <- bind_rows(df_overlap, df_overlap_iso)
#     
#   }
#   
#   df_overlap
#   
#   # ------------------------------------------------------------------------------
#   # Plots
#   # ------------------------------------------------------------------------------
#   
#   
#   condition_title <- list(
#     "HH117" = "inflamed vs uninflamed",
#     "HH119" = "ileum vs colon"
#   )
#   
#   
#   # N shared clones per isotype
#   df_overlap %>% 
#     ggplot(aes(x = isotype, y = n_shared_clones, fill = isotype)) + 
#     geom_col() + 
#     geom_text(aes(label = n_shared_clones), vjust = -0.3, size = 3.5) + 
#     scale_fill_manual(values = isotype_colors_custom) + 
#     labs(
#       x = "Isotype", y = "N shared clones",
#       title = glue("{HH}: Clonal sharing between conditions, within isotype (LP PCs)"),
#       subtitle = glue("N clones present in both conditions ({condition_title[[HH]]})")
#     ) + 
#     theme_bw() + 
#     theme(legend.position = "none")
#   
#   ggsave(glue("{outdir}/{HH}_pc_condition_overlap_n_shared_clones.png"), width = 8, height = 6)
#   
#   # % overlap of smaller set per isotype
#   df_overlap %>% 
#     ggplot(aes(x = isotype, y = pct_overlap_smaller_set, fill = isotype)) + 
#     geom_col() + 
#     geom_text(aes(label = round(pct_overlap_smaller_set, 1)), vjust = -0.3, size = 3.5) + 
#     scale_fill_manual(values = isotype_colors_custom) + 
#     ylim(0, 100) + 
#     labs(
#       x = "Isotype", y = "% overlap (smaller set)",
#       title = glue("{HH}: Clonal overlap between conditions within isotype (LP PCs)"),
#       subtitle = glue("Condition: {condition_title[[HH]]}")
#     ) + 
#     theme_bw() + 
#     theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   ggsave(glue("{outdir}/{HH}_pc_condition_overlap_pct_smaller_set.png"), width = 8, height = 6)
#   
#   # Morisita-Horn index per isotype
#   df_overlap %>% 
#     ggplot(aes(x = isotype, y = morisita_horn, fill = isotype)) + 
#     geom_col() + 
#     geom_text(aes(label = round(morisita_horn, 2)), vjust = -0.3, size = 3.5) + 
#     scale_fill_manual(values = isotype_colors_custom) + 
#     ylim(0, 1) + 
#     labs(
#       x = "Isotype", y = "Morisita-Horn index",
#       title = glue("{HH}: Clonal overlap between conditions within isotype (LP PCs)"),
#       subtitle = glue("Condition: {condition_title[[HH]]}")
#     ) + 
#     theme_bw() + 
#     theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   ggsave(glue("{outdir}/{HH}_pc_condition_overlap_morisita_horn.png"), width = 8, height = 6)
# 
# })
# 
# 
# # ==============================================================================
# # Zoom in on IgA clones
# # ==============================================================================
# 
# lapply(patients, function(HH){
#   
#   # HH <- "HH119"
#   
#   # Add condition and filter for IgA cells
#   df_iga <- df_both %>% 
#     filter(
#       patient_id == HH & c_call_grouped %in% c("IGHA1", "IGHA2")
#     ) %>% 
#     mutate(
#       condition = sample_clean %>% str_remove(glue("{HH}-"))
#     )
#   
#   df_iga$condition %>% table()
#   
#   # Calculate percentage of IgA1/IgA2
#   df_plot <- df_iga %>% 
#     group_by(clone_subgroup_id_90_similarity) %>% 
#     count(c_call_grouped, name = "n_cells") %>% 
#     mutate(
#       n_cells_total = sum(n_cells),
#       IGA1_percentage = ifelse(c_call_grouped == "IGHA1", (n_cells/n_cells_total) * 100, NA)
#       # IGA2_percentage = ifelse(c_call_grouped == "IGHA2", (n_cells/n_cells_total) * 100, NA)
#     ) %>% 
#     ungroup() %>% 
#     filter(!is.na(IGA1_percentage)) %>%  
#     count(IGA1_percentage)
#   
#   max_ys <- sort(df_plot$n, decreasing = TRUE)[1:2]
#   
#   df_plot %>% 
#     ggplot(aes(x = IGA1_percentage, y = n)) + 
#     geom_point(color = "darkblue", alpha = 0.6) + 
#     geom_line(color = "darkblue", alpha = 0.6) + 
#     scale_y_break(c(max_ys[2]+5, max_ys[1]-10), scales = 1/5) +
#     scale_y_continuous(
#       breaks = scales::breaks_width(20)
#     ) +
#     scale_x_continuous(
#       breaks = scales::breaks_width(10)
#     ) +
#     theme_bw() + 
#     labs(
#       title = glue("{HH}: IgA clones - ratio of IgA1 to IgA2"), 
#       x = "% IgA1/(IgA1+IgA2)", 
#       y = "N clones"
#     )
#   
#   ggsave(glue("{outdir}/{HH}_IgA_ratio.png"), width = 15, height = 6)
# 
# })
# 
