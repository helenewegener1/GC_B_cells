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
resolve_LC_files <- grep("resolve_LC_3_definitions", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
patients

# Prep output
outdir = glue("45_immcantation/plot/18_90_similarity/09_clonal_sharing_largest_clone_removed/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Distance between follicles 
fol_distances_list <- lapply(
  patients, 
  function(HH) read_csv(glue("00_data/{HH}_Claude_follicle_distance_matrix.csv"))
) %>% 
  setNames(patients)

# ------------------------------------------------------------------------------
# For each clone, count presents in follicles (GC B cells)
# ------------------------------------------------------------------------------

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>%
    filter(
      locus == "IGH"
      # !is.na(manual_ADT_full_ID)
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows()

# Define largest CRC clone
large_crc_clone <- df_both %>%
  filter(patient_id == "HH119", L1_annotation == "GC_B_cells") %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  head(1) %>% 
  pull(clone_subgroup_id_90_similarity)

# Remove largest clone from patient HH119
df_both <- df_both %>% 
  filter(
    !(patient_id == "HH119" & clone_subgroup_id_90_similarity == large_crc_clone)
    # patient_id == "HH119"
  ) 

# patients <- df_both$patient_id %>% unique()

for (n_min_cells in c(1, 2)){
  
  # Min number of cells in follicle to be sure of presents 
  # n_min_cells <- 2
  
  # Clones are present in N follicles  
  df_n_fol <- df_both %>% 
    filter(
      L1_annotation == "GC_B_cells",
      !is.na(manual_ADT_full_ID)
    ) %>% # should it only be across GC B cells?
    group_by(clone_subgroup_id_90_similarity, patient_id) %>% 
    count(manual_ADT_full_ID) %>% 
    filter(n >= n_min_cells) %>%
    count(manual_ADT_full_ID) %>% 
    count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
    dplyr::rename(n_follicles = n) %>% 
    ungroup()
  
  # N GC B cells per clone 
  df_N_cells <- df_both %>% 
    filter(
      L1_annotation == "GC_B_cells",
      !is.na(manual_ADT_full_ID)
    ) %>% 
    filter(clone_subgroup_id_90_similarity %in% df_n_fol$clone_subgroup_id_90_similarity) %>% 
    count(clone_subgroup_id_90_similarity) %>% 
    dplyr::rename(clone_size = n)
  
  # combine
  nrow(df_n_fol) == nrow(df_N_cells)
  df_plot <- df_n_fol %>% 
    left_join(df_N_cells, by = "clone_subgroup_id_90_similarity") %>% 
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
    group_by(n_follicles, patient_id) %>% 
    mutate(
      mean_clone_size = mean(clone_size) %>% round(1),
      median_clone_size = median(clone_size)
    ) %>% 
    ungroup()
  
  # N 
  df_N <- df_plot %>%
    select(-clone_subgroup_id_90_similarity) %>% 
    count(patient_id, n_follicles, mean_clone_size, median_clone_size, sort = TRUE) 
  
  # Jitter plot split by patient 
  df_plot %>% 
    ggplot(aes(x = patient_id, y = n_follicles)) + 
    geom_jitter(
      aes(color = clone_size_group, size = clone_size_group), 
      alpha = 0.5, width = 0.30, height = 0#, size = 2.5
    ) + 
    scale_color_viridis_d(option = "plasma", direction = -1) +
    geom_text(
      data = df_N,
      aes(x = patient_id, y = n_follicles, label = glue("{n} clones ({mean_clone_size}; {median_clone_size})")),
      size = 3,
      position = position_nudge(x = 0.45)
    ) +
    theme_bw() + 
    scale_y_continuous(
      breaks = scales::breaks_width(5),
      minor_breaks = scales::breaks_width(1)
    ) + 
    labs(
      x = "Patient ID",
      y = "N follicles",
      title = "Clonal sharing across follicles",
      subtitle = glue("For each clone (GC B cells), how many follicles is it present in (at least {n_min_cells} cells)"),
      caption = "N clones (mean clone size; median clone size)"
    )
  
  ggsave(glue("{outdir}/clones_in_n_follicles_{n_min_cells}.png"), width = 12, height = 6)
  
  # Mean clonal size VS presents in N follicles
  range_fun <- function(x) {
    data.frame(y = mean(x), ymin = min(x), ymax = max(x))
  }
  
  df_plot %>%
    select(patient_id, mean_clone_size, n_follicles) %>%
    distinct() %>% 
    ggplot(aes(x = n_follicles, y = mean_clone_size, color = patient_id)) +
    geom_point(size = 2.5, alpha = 0.7) +
    stat_summary(
      data = df_plot,
      aes(y = clone_size),
      fun.data = range_fun,
      geom = "errorbar",
      width = 0.4
    ) + 
    theme_bw() +
    # scale_y_break(c(100, 7050), scales = 1/4) +
    # facet_wrap(vars(patient_id)) +
    scale_x_continuous(
      breaks = 1:10,
      # limits = c(min(follicle_breaks) - 0.5, max(follicle_breaks) + 0.5),
      expand = c(0.05, 0),
      minor_breaks = scales::breaks_width(1)
    ) + 
    labs(
      color = "Patient ID",
      x = "N follicles",
      y = "Mean clone size",
      title = "Mean clonal size VS presents in N follicles",
      subtitle = glue("For each clone (GC B cells), how many follicles is it present in (at least {n_min_cells} cells)")
    )
  
  ggsave(glue("{outdir}/clones_in_n_follicles_errorbar_{n_min_cells}.png"), width = 10, height = 6.5)
  
  # ------------------------------------------------------------------------------
  # Clonal sharing proportions 
  # ------------------------------------------------------------------------------
  
  # These plots (above and this one under) only include clones/cells that are present 
  # in the abundance of n_min_cells. This is both the number and the plots. 
  
  df_totals <- df_plot %>% 
    group_by(patient_id) %>% 
    summarise(
      n_clones = n(),
      total_cells = sum(clone_size),
      .groups = "drop"
    )
  
  df_plot %>% 
    mutate(n_follicle_fct = n_follicles %>% as.factor()) %>% 
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
    theme_bw() +
    labs(
      fill = "N follciles",
      x = "Patient ID",
      y = "% of clones",
      title = "Percentage of clonal sharing",
      subtitle = glue("For each clone (GC B cells), how many follicles is it present in (at least {n_min_cells} cells)")
    )
  
  ggsave(glue("{outdir}/clones_in_n_follicles_proportions_{n_min_cells}.png"), width = 10, height = 6.5)
  
}

# ------------------------------------------------------------------------------
# Different isotypes in different follicles? (GC B cells)
# ------------------------------------------------------------------------------

df_shared_clones <- df_n_fol %>% filter(n_follicles > 1)

# ------------------------------------------------------------------------------
# Build the pie-chart data for a given cell population (all cells, GC B cells, ...)
# ------------------------------------------------------------------------------
build_isotype_pie_data <- function(df_both, df_shared_clones, cell_filter = TRUE) {
  
  df_filtered <- df_both %>% filter({{ cell_filter }})
  
  # ---- ordered site lookup per patient, plus a trailing "Total" column ----
  site_lookup <- df_filtered %>% 
    distinct(patient_id, sample_clean_fol) %>% 
    mutate(
      follicle_num = str_extract(sample_clean_fol, "(?<=Fol-)\\d+") %>% as.integer(),
      is_follicle = !is.na(follicle_num),
      sample_clean_fol_plot = sample_clean_fol %>% 
        str_remove("HH11\\d+-") %>% 
        str_remove("SI-PP_|SI-PP-nonINF_")
    ) %>% 
    arrange(patient_id, is_follicle, follicle_num, sample_clean_fol) %>% 
    group_by(patient_id) %>% 
    mutate(x = row_number()) %>% 
    ungroup() %>% 
    select(patient_id, sample_clean_fol_plot, x)
  
  total_lookup <- site_lookup %>% 
    group_by(patient_id) %>% 
    summarise(x = max(x) + 1, .groups = "drop") %>% 
    mutate(sample_clean_fol_plot = "Combined")
  
  site_lookup <- bind_rows(site_lookup, total_lookup)
  
  # ---- counts per clone x site x isotype ----
  counts_by_site <- df_filtered %>% 
    filter(!is.na(c_call)) %>% 
    inner_join(df_shared_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
    mutate(
      sample_clean_fol_plot = sample_clean_fol %>% 
        str_remove("HH11\\d+-") %>% 
        str_remove("SI-PP_|SI-PP-nonINF_")
    ) %>% 
    count(patient_id, clone_subgroup_id_90_similarity, sample_clean_fol_plot, c_call)
  
  # ---- counts per clone summed across ALL sites -> the "Total" column ----
  counts_total <- counts_by_site %>% 
    group_by(patient_id, clone_subgroup_id_90_similarity, c_call) %>% 
    summarise(n = sum(n), .groups = "drop") %>% 
    mutate(sample_clean_fol_plot = "Combined")
  
  plot_data <- bind_rows(counts_by_site, counts_total) %>% 
    left_join(site_lookup, by = c("patient_id", "sample_clean_fol_plot"))
  
  clone_id_lookup <- plot_data %>% 
    distinct(patient_id, clone_subgroup_id_90_similarity) %>% 
    arrange(patient_id, clone_subgroup_id_90_similarity) %>% 
    group_by(patient_id) %>% 
    mutate(y = row_number()) %>% 
    ungroup()
  
  isotype_cols <- plot_data$c_call %>% unique() %>% na.omit() %>% as.character()
  
  pie_data <- plot_data %>% 
    left_join(clone_id_lookup, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
    pivot_wider(
      id_cols = c(patient_id, clone_subgroup_id_90_similarity, x, y),
      names_from = c_call, 
      values_from = n, 
      values_fill = 0
    ) %>% 
    mutate(total_n = rowSums(across(all_of(isotype_cols))))
  
  list(pie_data = pie_data, clone_id_lookup = clone_id_lookup, site_lookup = site_lookup, isotype_cols = isotype_cols)
}

# ------------------------------------------------------------------------------
# Plot one patient's pie grid, with a dashed divider before the "Total" column
# ------------------------------------------------------------------------------
make_pie_plot <- function(patient_name, dat) {
  
  data_sub     <- dat$pie_data        %>% filter(patient_id == patient_name)
  lookup_sub   <- dat$clone_id_lookup %>% filter(patient_id == patient_name)
  x_lookup_sub <- dat$site_lookup     %>% filter(patient_id == patient_name)
  
  last_real_x <- x_lookup_sub %>% filter(sample_clean_fol_plot != "Combined") %>% pull(x) %>% max()
  
  ggplot() + 
    geom_vline(xintercept = last_real_x + 0.5, color = "grey30") +
    geom_scatterpie(
      data = data_sub, 
      aes(x = x, y = y, r = 0.35),
      cols = dat$isotype_cols, 
      color = NA
    ) + 
    geom_text(
      data = data_sub,
      aes(x = x, y = y, label = total_n),
      size = 2.5,
      nudge_y = -0.5
    ) +
    scale_fill_manual(
      values = isotype_colors_custom,
      breaks = c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE")
    ) + 
    scale_x_continuous(
      breaks = x_lookup_sub$x, 
      labels = x_lookup_sub$sample_clean_fol_plot,
      minor_breaks = scales::breaks_width(1)
    ) + 
    scale_y_continuous(
      breaks = lookup_sub$y, 
      labels = lookup_sub$clone_subgroup_id_90_similarity,
      minor_breaks = scales::breaks_width(1)
    ) + 
    coord_equal() + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Site", y = "Clone", fill = "Isotype")
}

# ------------------------------------------------------------------------------
# Build + save plots for all patients
# ------------------------------------------------------------------------------
plot_isotype_pies <- function(dat, filename_suffix, title_suffix) {
  
  width <- 12
  height <- 10
  
  for (HH in names(patient_names)) {
    
    # dat <- dat_all
    # HH <- "HH117"
    # title_suffix <- "bla"
    
    p <- make_pie_plot(HH, dat) + 
      labs(title = glue("{HH}: Isotype usage in follicle-shared clones - {title_suffix}"))
    
    plot_width <- if (HH == "HH119") width + 4 else width
    plot_height <- if (HH == "HH119") height + 4 else height
    
    ggsave(
      glue("{outdir}/{HH}_isotype_follicle_shared_clones_{filename_suffix}.png"),
      plot = p, width = plot_width, height = plot_height
    )
  }
}

# ------------------------------------------------------------------------------
# Run for both populations
# ------------------------------------------------------------------------------
dat_all <- build_isotype_pie_data(df_both, df_shared_clones, cell_filter = TRUE)
plot_isotype_pies(dat_all, "all_cells", "all cells")

dat_gc <- build_isotype_pie_data(df_both, df_shared_clones, cell_filter = L1_annotation == "GC_B_cells")
plot_isotype_pies(dat_gc, "GC_B_cells", "GC B cells")


# ==============================================================================
# ==============================================================================
# ==============================================================================


# ------------------------------------------------------------------------------
# Isotype diversity per clone (from the Combined row = summed across follicles)
# ------------------------------------------------------------------------------
isotype_diversity <- function(dat) {
  
  combined_x <- dat$site_lookup %>% 
    filter(sample_clean_fol_plot == "Combined") %>% 
    select(patient_id, x)
  
  dat$pie_data %>% 
    inner_join(combined_x, by = c("patient_id", "x")) %>% 
    mutate(n_isotypes = rowSums(across(all_of(dat$isotype_cols), ~ . > 0))) %>% 
    select(patient_id, clone_subgroup_id_90_similarity, n_isotypes, total_n)
}

df_isotype_diversity <- isotype_diversity(dat_gc)  # or dat_all

# ------------------------------------------------------------------------------
# Plot: % of clones by number of distinct isotypes used
# ------------------------------------------------------------------------------

df_isotype_diversity %>% 
  count(patient_id, n_isotypes) %>% 
  group_by(patient_id) %>% 
  mutate(pct = n / sum(n) * 100) %>% 
  ungroup() %>% 
  ggplot(aes(x = factor(n_isotypes), y = pct)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.4, size = 3) +
  facet_wrap(~patient_id) +
  theme_bw() +
  labs(
    x = "N distinct isotypes used (combined across follicles)",
    y = "% of clones",
    title = "Isotype diversity within follicle-shared clones",
    subtitle = "How often a clone sticks to one isotype vs. diversifies across several"
  )

ggsave(glue("{outdir}/N_isotypes_barplot.png"))

# ==============================================================================
# ==============================================================================
# ==============================================================================

# ------------------------------------------------------------------------------
# N shared clones for each pair of follicles 
# ------------------------------------------------------------------------------

df_fols <- df_both %>% 
  filter(!is.na(manual_ADT_ID)) %>% 
  select(patient_id, clone_subgroup_id_90_similarity, manual_ADT_ID) 


for (HH in patients){
  
  # HH <- "HH119"
  # HH <- "HH117"
  
  # clone x follicle presence/absence matrix
  incidence <- df_shared_clones %>% 
    left_join(df_fols, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
    filter(patient_id == HH) %>% 
    select(clone_subgroup_id_90_similarity, manual_ADT_ID) %>% 
    distinct(clone_subgroup_id_90_similarity, manual_ADT_ID) %>% 
    mutate(present = 1) %>% 
    pivot_wider(names_from = manual_ADT_ID, values_from = present, values_fill = 0)
    
  mat <- as.matrix(incidence %>% select(-clone_subgroup_id_90_similarity))
  shared_mat <- t(mat) %*% mat
  
  # ------------------------------------------------------------------------------
  # N shared follicles per follicle VS total cells in that follicle
  # ------------------------------------------------------------------------------
  
  # N follicles each follicle shares >=1 clone with (excluding itself)
  shared_mat_off_diag <- shared_mat
  diag(shared_mat_off_diag) <- 0
  n_shared_follicles <- rowSums(shared_mat_off_diag > 0)
  
  df_degree <- tibble(
    manual_ADT_ID = names(n_shared_follicles),
    n_shared_follicles = as.integer(n_shared_follicles)
  )

  # Total GC B cells per follicle (same filtering as the rest of the pipeline)
  df_follicle_cells <- df_both %>% 
    filter(
      patient_id == HH, L1_annotation == "GC_B_cells", !is.na(manual_ADT_ID)
    ) %>% 
    count(manual_ADT_ID, name = "total_cells")
  
  df_degree_plot <- df_degree %>% 
    left_join(df_follicle_cells, by = "manual_ADT_ID")
  
  ggplot(df_degree_plot, aes(x = total_cells, y = n_shared_follicles)) +
    geom_point(size = 2.5, alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = manual_ADT_ID), size = 3) +
    theme_bw() +
    scale_y_continuous(
      breaks = scales::breaks_width(1),
      minor_breaks = scales::breaks_width(1)
    ) +
    scale_x_continuous(
      breaks = scales::breaks_width(20),
      minor_breaks = scales::breaks_width(10)
    ) + 
    labs(
      title = glue("{HH}: Clonal sharing degree VS follicle size - largest CRC clone removed"),
      subtitle = "At least 2 cells per follicle",
      x = "Total N GC B cells in follicle",
      y = "N follicles sharing \u2265 1 clone with this follicle"
    )
  
  ggsave(glue("{outdir}/{HH}_n_shared_follicles_VS_total_cells.png"), width = 8, height = 6.5)
  
  
  # ------------------------------------------------------------------------------
  # Pairwise clonal sharing between follicles
  # ------------------------------------------------------------------------------
  
  # long format, one row per unordered follicle pair
  shared_full <- shared_mat %>% 
    as.data.frame() %>% 
    rownames_to_column("follicle_1") %>% 
    pivot_longer(-follicle_1, names_to = "follicle_2", values_to = "n_shared_clones") %>% 
    mutate(
      follicle_1 = str_split_i(follicle_1, "-", 2) %>% as.double(),
      follicle_2 = str_split_i(follicle_2, "-", 2) %>% as.double()
    ) %>% 
    filter(follicle_1 < follicle_2)
  
  follicle_breaks <- sort(unique(c(shared_full$follicle_1, shared_full$follicle_2)))
  
  ggplot(shared_full, aes(x = follicle_1, y = follicle_2, fill = n_shared_clones)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_shared_clones), size = 2.5) +
    scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey90") +
    coord_equal() +
    theme_bw() +
    scale_x_continuous(
      breaks = follicle_breaks,
      limits = c(min(follicle_breaks) - 0.5, max(follicle_breaks) + 0.5),
      expand = c(0, 0),
      minor_breaks = scales::breaks_width(1)
    ) + 
    scale_y_continuous(
      breaks = follicle_breaks,
      limits = c(min(follicle_breaks) - 0.5, max(follicle_breaks) + 0.5),
      expand = c(0, 0),
      minor_breaks = scales::breaks_width(1)
    ) + 
    labs(
      x = "Follicle", y = "Follicle", fill = "N shared\nclones",
      title = glue("{HH}: Pairwise clonal sharing between follicles")
    )
  
  ggsave(glue("{outdir}/{HH}_pairwise_clonal_sharing_follicles.png"), width = 8, height = 8)
  
  # ------------------------------------------------------------------------------
  # N shared clones for each pair of follicles - Include distance
  # ------------------------------------------------------------------------------
  
  # Get pairs in same order
  shared_full_reorder <- shared_full %>% 
    mutate(
      f_lo = pmin(follicle_1, follicle_2),
      f_hi = pmax(follicle_1, follicle_2)
    )
  
  fol_distances_reorder <- fol_distances_list[[HH]] %>% 
    filter(!(str_detect(follicle_1, "_R|_L") | str_detect(follicle_2, "_R|_L"))) %>% # TODO figure out what is follicle 13 and what is follicle 15 
    mutate(
      f_lo = pmin(follicle_1, follicle_2) %>% as.double(),
      f_hi = pmax(follicle_1, follicle_2) %>% as.double()
    ) %>% 
    select(!c(follicle_1, follicle_2)) 
  
  # Join dataframes of N shared clones and distance
  shared_full_dist <- shared_full_reorder %>% 
    left_join(fol_distances_reorder, by = c("f_lo", "f_hi")) %>% 
    select(!c(follicle_1, follicle_2)) %>% 
    mutate(
      n_shared_clones_fct = n_shared_clones %>% as.factor()
    )
  
  # For fragments individually
  if (!"fragment" %in% colnames(shared_full_dist)) {
    shared_full_dist <- shared_full_dist %>% mutate(fragment = "All")
  }
  
  shared_full_dist_plot <- shared_full_dist %>% filter(!is.na(fragment))
  
  # Plot
  shared_full_dist_plot %>% 
    ggplot(aes(x = distance_px, y = n_shared_clones_fct, color = fragment)) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75), alpha = 0.5) + 
    theme_bw() + 
    labs(
      title = glue("{HH}: Pairwise distance and N shared clones between follicles"), 
      x = "Distance between pairs of follicles",
      y = "N shared clones between pairs of follciles"
    )
  
  ggsave(glue("{outdir}/{HH}_pairwise_clonal_sharing_follicles_VS_distance.png"), width = 8, height = 8)

}


