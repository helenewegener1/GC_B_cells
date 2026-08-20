library(glue)
library(tidyverse)
library(scatterpie)
library(patchwork)

source("10_broad_annotation/script/color_palette.R")

# Adapted from 60_PC_clones/08_PC_clonal_sharing.R.
#
# ==============================================================================
# NOTES / ASSUMPTIONS
# ==============================================================================
# - clone = "CTstrict" instead of clone_subgroup_id_90_similarity.
# - SITE_PATTERN restricts to compartments compared for sharing -- defaults to
#   follicles within SI-PP (Tfh/GC dynamics within Peyer's patch), unlike
#   60_PC_clones which looks at "LP" (lamina propria dissemination). Change if
#   you want a different comparison.
# - The isotype pie-chart grid in 60_PC_clones/08 (c_call, IGH-only) has no
#   direct TCR equivalent -- replaced here with a TRBV gene usage pie grid
#   per compartment-shared clone. Spot-check `df_both$CTgene %>% head()`
#   before trusting the extract_trbv() regex.
# ==============================================================================

CELL_TYPE_FILTER <- "Tfh_cells"
SITE_PATTERN <- "SI-PP"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("70_TCR_analysis/out/rds")
patient_files <- grep(glue("_TCR_{CELL_TYPE_FILTER}\\.rds$"), rds_files, value = TRUE) %>%
  grep("^HH", ., value = TRUE)
patients <- lapply(patient_files, function(x) str_split_i(x, "_", 1)) %>% unlist()

outdir <- glue("70_TCR_analysis/plot/04_clonal_sharing/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

df_both <- lapply(patients, function(HH) {
  readRDS(glue("70_TCR_analysis/out/rds/{HH}_TCR_{CELL_TYPE_FILTER}.rds")) %>%
    filter(str_detect(sample_clean, SITE_PATTERN))
}) %>% bind_rows()

extract_trbv <- function(ctgene) str_extract(ctgene, "TRBV[0-9A-Za-z-]+")

df_both <- df_both %>% mutate(v_trb_call = extract_trbv(CTgene))

# ------------------------------------------------------------------------------
# For each clone, count presence across compartments per patient
# ------------------------------------------------------------------------------

for (n_min_cells in c(1, 2)) {

  # n_min_cells <- 2

  df_n_compartment <- df_both %>%
    group_by(patient_id, CTstrict) %>%
    count(sample_clean) %>%
    filter(n >= n_min_cells) %>%
    count(sample_clean) %>%
    count(CTstrict, sort = TRUE) %>%
    dplyr::rename(n_compartments = n) %>%
    ungroup()

  df_N_cells <- df_both %>%
    right_join(df_n_compartment, by = c("patient_id", "CTstrict")) %>%
    count(CTstrict, patient_id) %>%
    dplyr::rename(clone_size = n)

  df_plot <- df_n_compartment %>%
    left_join(df_N_cells, by = c("patient_id", "CTstrict")) %>%
    mutate(
      clone_size_group = case_when(
        clone_size > 1  & clone_size <= 5   ~ "2-5",
        clone_size > 5  & clone_size <= 10  ~ "6-10",
        clone_size > 10 & clone_size <= 20  ~ "11-20",
        clone_size > 20 & clone_size <= 50  ~ "21-50",
        clone_size > 50 & clone_size <= 100 ~ "51-100",
        clone_size > 100 ~ "100+"
      ),
      clone_size_group = factor(clone_size_group, levels = c("2-5", "6-10", "11-20", "21-50", "51-100", "100+"))
    ) %>%
    group_by(n_compartments, patient_id) %>%
    mutate(
      mean_clone_size = mean(clone_size) %>% round(1),
      median_clone_size = median(clone_size)
    ) %>%
    ungroup()

  df_N <- df_plot %>%
    select(-CTstrict) %>%
    count(patient_id, n_compartments, mean_clone_size, median_clone_size, sort = TRUE)

  df_plot %>%
    ggplot(aes(x = patient_id, y = n_compartments)) +
    geom_jitter(
      aes(color = clone_size_group, size = clone_size_group),
      alpha = 0.5, width = 0.30, height = 0.2
    ) +
    scale_color_viridis_d(option = "plasma", direction = -1) +
    geom_text(
      data = df_N,
      aes(x = patient_id, y = n_compartments, label = glue("{n} clones ({mean_clone_size}; {median_clone_size})")),
      size = 3, position = position_nudge(x = 0.45)
    ) +
    theme_bw() +
    scale_y_continuous(
      breaks = scales::breaks_width(1),
      minor_breaks = scales::breaks_width(1)
    ) +
    labs(
      x = "Patient ID",
      y = "N compartments",
      title = glue("Clonal sharing across {CELL_TYPE_FILTER} compartments"),
      subtitle = glue("For each clone, how many compartments is it present in (at least {n_min_cells} cells)"),
      caption = "N clones (mean clone size; median clone size)"
    )

  ggsave(glue("{outdir}/clones_in_n_compartments_{n_min_cells}.png"), width = 12, height = 6)

  # Proportions
  df_totals <- df_plot %>%
    group_by(patient_id) %>%
    summarise(n_clones = n(), total_cells = sum(clone_size), .groups = "drop")

  df_plot %>%
    mutate(n_compartment_fct = n_compartments %>% as.factor()) %>%
    ggplot(aes(x = patient_id, fill = n_compartment_fct)) +
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    geom_bar(position = "fill") +
    geom_text(
      data = df_totals,
      aes(x = patient_id, y = 0.9, label = glue("{n_clones} clones\n({total_cells} cells)")),
      inherit.aes = FALSE
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    labs(
      fill = "N compartments",
      x = "Patient ID",
      y = "% of clones",
      title = glue("Percentage of clonal sharing across {CELL_TYPE_FILTER} compartments"),
      subtitle = glue("For each clone, how many compartments is it present in (at least {n_min_cells} cells)")
    )

  ggsave(glue("{outdir}/clones_in_n_compartments_proportions_{n_min_cells}.png"), width = 10, height = 6.5)

}

# ------------------------------------------------------------------------------
# TRBV usage pie grid for compartment-shared clones (isotype-pie analog)
# ------------------------------------------------------------------------------

n_min_cells <- 2

df_n_compartment <- df_both %>%
  group_by(patient_id, CTstrict) %>%
  count(sample_clean) %>%
  filter(n >= n_min_cells) %>%
  count(sample_clean) %>%
  count(CTstrict, sort = TRUE) %>%
  dplyr::rename(n_compartments = n) %>%
  ungroup()

df_shared_clones <- df_both %>%
  left_join(df_n_compartment, by = c("patient_id", "CTstrict")) %>%
  filter(n_compartments > 1) %>%
  count(patient_id, CTstrict, sort = TRUE) %>%
  group_by(patient_id) %>%
  slice_max(order_by = n, n = 20) %>%
  ungroup()

build_trbv_pie_data <- function(df_both, df_shared_clones) {

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

  counts_by_site <- df_both %>%
    filter(!is.na(v_trb_call)) %>%
    inner_join(df_shared_clones, by = c("patient_id", "CTstrict")) %>%
    mutate(
      sample_clean_fol_plot = sample_clean_fol %>%
        str_remove("HH11\\d+-") %>%
        str_remove("SI-PP_|SI-PP-nonINF_")
    ) %>%
    count(patient_id, CTstrict, sample_clean_fol_plot, v_trb_call)

  counts_total <- counts_by_site %>%
    group_by(patient_id, CTstrict, v_trb_call) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(sample_clean_fol_plot = "Combined")

  plot_data <- bind_rows(counts_by_site, counts_total) %>%
    left_join(site_lookup, by = c("patient_id", "sample_clean_fol_plot"))

  clone_id_lookup <- plot_data %>%
    distinct(patient_id, CTstrict) %>%
    arrange(patient_id, CTstrict) %>%
    group_by(patient_id) %>%
    mutate(y = row_number()) %>%
    ungroup()

  v_cols <- plot_data$v_trb_call %>% unique() %>% na.omit() %>% as.character()

  pie_data <- plot_data %>%
    left_join(clone_id_lookup, by = c("patient_id", "CTstrict")) %>%
    pivot_wider(
      id_cols = c(patient_id, CTstrict, x, y),
      names_from = v_trb_call, values_from = n, values_fill = 0
    ) %>%
    mutate(total_n = rowSums(across(all_of(v_cols))))

  list(pie_data = pie_data, clone_id_lookup = clone_id_lookup, site_lookup = site_lookup, v_cols = v_cols)

}

make_trbv_pie_plot <- function(patient_name, dat) {

  data_sub     <- dat$pie_data        %>% filter(patient_id == patient_name)
  lookup_sub   <- dat$clone_id_lookup %>% filter(patient_id == patient_name)
  x_lookup_sub <- dat$site_lookup     %>% filter(patient_id == patient_name)

  last_real_x <- x_lookup_sub %>% filter(sample_clean_fol_plot != "Combined") %>% pull(x) %>% max()

  ggplot() +
    geom_vline(xintercept = last_real_x + 0.5, color = "grey30") +
    geom_scatterpie(
      data = data_sub, aes(x = x, y = y, r = 0.35), cols = dat$v_cols, color = NA
    ) +
    geom_text(
      data = data_sub, aes(x = x, y = y, label = total_n), size = 2.5, nudge_y = -0.5
    ) +
    scale_fill_viridis_d(option = "turbo") +
    scale_x_continuous(breaks = x_lookup_sub$x, labels = x_lookup_sub$sample_clean_fol_plot, minor_breaks = scales::breaks_width(1)) +
    scale_y_continuous(breaks = lookup_sub$y, labels = lookup_sub$CTstrict, minor_breaks = scales::breaks_width(1)) +
    coord_equal() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Site", y = "Clone", fill = "TRBV gene")

}

dat_all <- build_trbv_pie_data(df_both, df_shared_clones)

for (HH in unique(dat_all$pie_data$patient_id)) {

  # HH <- "HH117"
  p <- make_trbv_pie_plot(HH, dat_all) +
    labs(title = glue("{HH}: TRBV usage in {CELL_TYPE_FILTER} compartment-shared clones"))

  ggsave(glue("{outdir}/{HH}_TRBV_compartment_shared_clones.png"), plot = p, width = 6, height = 10)

}
