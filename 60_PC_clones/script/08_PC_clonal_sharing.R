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
outdir = glue("60_PC_clones/plot/08_clonal_sharing/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
    filter(
      locus == "IGH" & L1_annotation == "PCs" & str_detect(sample_clean, "LP")
    )
}) %>% bind_rows()

# ------------------------------------------------------------------------------
# Define top clones
# ------------------------------------------------------------------------------


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
    ungroup()
  
  # N 
  df_N <- df_plot %>%
    select(-clone_subgroup_id_90_similarity) %>% 
    count(patient_id, n_compartments, mean_clone_size, median_clone_size, sort = TRUE) 
  
  # Jitter plot split by patient 
  df_plot %>% 
    ggplot(aes(x = patient_id, y = n_compartments)) + 
    geom_jitter(
      aes(color = clone_size_group, size = clone_size_group), 
      alpha = 0.5, width = 0.30, height = 0.2#, size = 2.5
    ) + 
    scale_color_viridis_d(option = "plasma", direction = -1) +
    geom_text(
      data = df_N,
      aes(x = patient_id, y = n_compartments, label = glue("{n} clones ({mean_clone_size}; {median_clone_size})")),
      size = 3,
      position = position_nudge(x = 0.45)
    ) +
    theme_bw() + 
    scale_y_continuous(
      breaks = scales::breaks_width(1),
      minor_breaks = scales::breaks_width(1)
    ) + 
    labs(
      x = "Patient ID",
      y = "N compartments",
      title = "Clonal sharing across LP PC compartments",
      subtitle = glue("For each clone, how many compartments is it present in (at least {n_min_cells} cells)"),
      caption = "N clones (mean clone size; median clone size)"
    )
  
  ggsave(glue("{outdir}/clones_in_n_compartments_{n_min_cells}.png"), width = 12, height = 6)
  
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
    theme_bw() +
    labs(
      fill = "N compartments",
      x = "Patient ID",
      y = "% of clones",
      title = "Percentage of clonal sharing across LP PC compartments",
      subtitle = glue("For each clone, how many compartments is it present in (at least {n_min_cells} cells)")
    )
  
  ggsave(glue("{outdir}/clones_in_n_compartments_proportions_{n_min_cells}.png"), width = 10, height = 6.5)
  
}

# ------------------------------------------------------------------------------
# Different isotypes in different compartments?
# ------------------------------------------------------------------------------

n_min_cells <- 2

# Clones are present in both compartments  
df_n_compartment <- df_both %>% 
  group_by(patient_id, clone_subgroup_id_90_similarity) %>% 
  count(sample_clean) %>% 
  filter(n >= n_min_cells) %>%
  count(sample_clean) %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  dplyr::rename(n_compartments = n) %>% 
  ungroup()

# View top 20 clones per patient
df_shared_clones <- df_both %>%
  left_join(df_n_compartment, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  filter(n_compartments > 1) %>% 
  count(patient_id, clone_subgroup_id_90_similarity, sort = TRUE) %>%
  group_by(patient_id) %>%
  slice_max(order_by = n, n = 20) %>%
  ungroup()

df_shared_clones


# ------------------------------------------------------------------------------
# Build the pie-chart data for a given cell population 
# ------------------------------------------------------------------------------
build_isotype_pie_data <- function(df_both, df_shared_clones) {

  # ---- ordered site lookup per patient, plus a trailing "Total" column ----
  site_lookup <- df_both %>% 
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
  counts_by_site <- df_both %>% 
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
plot_isotype_pies <- function(dat) {

  for (HH in names(patient_names)) {
    
    # dat <- dat_all
    # HH <- "HH117"
    # title_suffix <- "bla"
    
    p <- make_pie_plot(HH, dat) + 
      labs(title = glue("{HH}: Isotype usage in LP PC compartment-shared clones"))
    
    ggsave(
      glue("{outdir}/{HH}_isotype_compartment_shared_clones.png"),
      plot = p, width = 6, height = 10
    )
  }
}

# ------------------------------------------------------------------------------
# Run for both populations
# ------------------------------------------------------------------------------
dat_all <- build_isotype_pie_data(df_both, df_shared_clones)
plot_isotype_pies(dat_all)

