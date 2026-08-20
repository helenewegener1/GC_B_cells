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
outdir = glue("45_immcantation/plot/13_clonal_sharing/")
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
  readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
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
    scale_y_break(c(10, 33), scales = 1/4) +
    scale_y_continuous(
      breaks = scales::breaks_width(1),
      minor_breaks = scales::breaks_width(1)
    ) + 
    theme(
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank()
    ) + 
    labs(
      x = "Patient ID",
      y = "N follicles",
      title = "Clonal sharing across follicles",
      subtitle = glue("For each clone (GC B cells), how many follicles is it present in (at least {n_min_cells} cells)"),
      caption = "N clones (mean clone size; median clone size)"
    )
  
  ggsave(glue("{outdir}/clones_in_n_follicles_{n_min_cells}.png"), width = 12, height = 6)

  # ------------------------------------------------------------------------------
  # Permutations
  # ------------------------------------------------------------------------------
  # Null model: shuffle the follicle label across cells (within each patient),
  # keeping the clone label fixed, then recompute "N follicles per clone"
  # exactly as above. Shuffling breaks any real association between clone
  # identity and follicle location while preserving both clone sizes and
  # follicle sizes (a permutation preserves the marginal counts of each
  # column), so any difference between the observed pattern and the permuted
  # pattern is attributable to genuine spatial/clonal structure rather than
  # sampling artifacts (e.g. follicle size differences).
  #
  # Only one of the two labels needs to be shuffled: composing two independent
  # random permutations (one on clone, one on follicle) is itself just a
  # single random permutation of the clone-follicle pairing, so shuffling both
  # produces the exact same null distribution as shuffling one while holding
  # the other fixed -- it was pure redundancy, not extra rigor.

  n_perm <- 1000
  set.seed(123)

  # Same base filter as df_n_fol above
  df_base <- df_both %>%
    filter(
      L1_annotation == "GC_B_cells",
      !is.na(manual_ADT_full_ID)
    ) %>%
    select(patient_id, clone_subgroup_id_90_similarity, manual_ADT_full_ID)

  # Same logic as the df_n_fol computation above, factored out so it can be
  # reused on both the real and the permuted labels
  compute_n_follicles_per_clone <- function(df, n_min_cells) {
    df %>%
      group_by(clone_subgroup_id_90_similarity, patient_id) %>%
      count(manual_ADT_full_ID) %>%
      filter(n >= n_min_cells) %>%
      count(manual_ADT_full_ID) %>%
      count(clone_subgroup_id_90_similarity, sort = TRUE) %>%
      dplyr::rename(n_follicles = n) %>%
      ungroup()
  }

  # Run once, per-clone (not pre-aggregated), so the exact same permutation
  # draws can answer both "how many clones share beyond chance" (Summary 1/2)
  # and "does that depend on clone size" (Summary 3) without re-permuting.
  perm_results_clone <- map_dfr(1:n_perm, function(perm_i) {

    df_base %>%
      group_by(patient_id) %>%
      mutate(
        manual_ADT_full_ID = sample(manual_ADT_full_ID)
      ) %>%
      ungroup() %>%
      compute_n_follicles_per_clone(n_min_cells) %>%
      mutate(perm_id = perm_i)

  })

  # A clone that, under a given permutation, never reaches n_min_cells in any
  # single follicle simply drops out of compute_n_follicles_per_clone()'s
  # output -- make that an explicit n_follicles = 0 rather than a missing row,
  # otherwise the null distribution below is biased upward.
  perm_results_clone_full <- df_plot %>%
    distinct(patient_id, clone_subgroup_id_90_similarity) %>%
    tidyr::crossing(perm_id = 1:n_perm) %>%
    left_join(
      perm_results_clone %>% select(patient_id, clone_subgroup_id_90_similarity, perm_id, n_follicles),
      by = c("patient_id", "clone_subgroup_id_90_similarity", "perm_id")
    ) %>%
    mutate(n_follicles = replace_na(n_follicles, 0))

  # perm_table_dir <- glue("45_immcantation/out/permutation_tests")
  # dir.create(perm_table_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Summary 1: observed vs. permuted % of clones per N-follicles category ---
  # (shows the whole distribution, not just one number, so you can see e.g.
  # whether it's specifically the "shared in 2-3 follicles" bucket that's
  # inflated, or the tail of highly-shared clones)

  df_obs_summary <- df_plot %>%
    count(patient_id, n_follicles, name = "n_clones") %>%
    group_by(patient_id) %>%
    mutate(pct_clones = n_clones / sum(n_clones) * 100) %>%
    ungroup()

  # For each permutation, restrict to clones that reached n_min_cells in >=1
  # follicle (n_follicles >= 1) before computing percentages -- this matches
  # how the observed bars are defined (df_plot only ever contains clones that
  # qualify somewhere), so bars and points are percentages of the same kind
  # of universe.
  perm_counts <- perm_results_clone_full %>%
    filter(n_follicles >= 1) %>%
    count(patient_id, perm_id, n_follicles, name = "n_clones")

  # Calculate % of clones in N follicles (How many % if clones are present in 1 follicle, 2 follicles...)
  perm_counts <- perm_counts %>%
    # group_by(patient_id) %>%
    # # complete() makes every category explicit (0 clones) in every permutation
    # complete(perm_id = 1:n_perm, n_follicles, fill = list(n_clones = 0)) %>%
    # ungroup() %>%
    group_by(patient_id, perm_id) %>%
    mutate(pct_clones = n_clones / sum(n_clones) * 100) %>%
    ungroup()
  
  # Calculate mean, lower and upper quantile for % of clones for each N follicle across all permutations 
  perm_ci <- perm_counts %>%
    group_by(patient_id, n_follicles) %>%
    summarise(
      perm_mean  = mean(pct_clones),
      perm_lower = quantile(pct_clones, 0.025),
      perm_upper = quantile(pct_clones, 0.975),
      .groups = "drop"
    )

  # Sanity check: perm_mean should print ~100 for every patient
  perm_ci %>% group_by(patient_id) %>% summarise(total_perm_mean = sum(perm_mean)) %>% print()

  # Join observed and permutated (summarized) counts and visualize 
  df_obs_summary %>%
    full_join(perm_ci, by = c("patient_id", "n_follicles")) %>%
    mutate(
      n_clones   = replace_na(n_clones, 0),
      pct_clones = replace_na(pct_clones, 0),
      perm_mean  = replace_na(perm_mean, 0),
      perm_lower = replace_na(perm_lower, 0),
      perm_upper = replace_na(perm_upper, 0)
    ) %>%
    ggplot(aes(x = n_follicles)) +
    geom_col(aes(y = pct_clones), fill = "steelblue", alpha = 0.8) +
    geom_pointrange(aes(y = perm_mean, ymin = perm_lower, ymax = perm_upper), color = "black", size = 0.2) +
    facet_wrap(~patient_id, scales = "free_x", ) +
    scale_x_continuous(breaks = scales::breaks_width(1), minor_breaks = scales::breaks_width(1)) +
    coord_cartesian(xlim = c(1, NA)) + 
    scale_y_continuous(breaks = scales::breaks_width(10), minor_breaks = scales::breaks_width(5)) +
    theme_bw() +
    labs(
      title = "Observed vs. permuted clonal sharing across follicles",
      subtitle = glue("Bars = observed % of clones; points/error bars = permutation null (mean ± 95% range, {n_perm} permutations)"),
      x = "N follicles",
      y = "% of clones"
    )

  ggsave(glue("{outdir}/clones_in_n_follicles_permutation_{n_min_cells}.png"), width = 12, height = 6)
  
  # --- Summary 2: single global statistic + empirical p-value ---
  # (useful as a one-line "is this significant" readout per patient, to quote
  # alongside Summary 1's plot - here: % of clones present in >=2 follicles.
  # Note this does NOT control for clone size -- see Summary 3 for that.)
  
  obs_shared_pct <- df_plot %>%
    group_by(patient_id) %>%
    summarise(pct_shared = mean(n_follicles >= 2) * 100, .groups = "drop")
  
  perm_shared_pct <- perm_results_clone_full %>%
    group_by(patient_id, perm_id) %>%
    summarise(pct_shared = mean(n_follicles >= 2) * 100, .groups = "drop") %>%
    ungroup()
  
  df_perm_test <- obs_shared_pct %>%
    left_join(
      perm_shared_pct %>%
        group_by(patient_id) %>%
        summarise(perm_mean = mean(pct_shared), perm_sd = sd(pct_shared), .groups = "drop"),
      by = "patient_id"
    ) %>%
    rowwise() %>%
    mutate(
      z_score = (pct_shared - perm_mean) / perm_sd,
      # H1: observed sharing is GREATER than expected by chance
      p_greater = mean(perm_shared_pct$pct_shared[perm_shared_pct$patient_id == patient_id] >= pct_shared),
      # H1: observed sharing is LESS than expected by chance
      p_less    = mean(perm_shared_pct$pct_shared[perm_shared_pct$patient_id == patient_id] <= pct_shared),
      # two-sided
      p_two_sided = min(1, 2 * min(p_greater, p_less))
    ) %>%
    ungroup() %>%
    mutate(n_min_cells = n_min_cells, n_perm = n_perm)
  
  print(df_perm_test)
  
  # write_csv(df_perm_test, glue("{perm_table_dir}/clonal_sharing_permutation_test_{n_min_cells}.csv"))

  # Same question, binned + significance rate (easier to read off directly)
  df_clone_test %>%
    mutate(
      clone_size_group = case_when(
        clone_size <= 5   ~ "2-5",
        clone_size <= 10  ~ "6-10",
        clone_size <= 20  ~ "11-20",
        clone_size <= 50  ~ "21-50",
        clone_size <= 100 ~ "51-100",
        TRUE ~ "100+"
      ),
      clone_size_group = factor(clone_size_group, levels = c("2-5", "6-10", "11-20", "21-50", "51-100", "100+"))
    ) %>%
    group_by(patient_id, clone_size_group) %>%
    summarise(n_clones = n(), pct_significant = mean(significant) * 100, .groups = "drop") %>%
    ggplot(aes(x = clone_size_group, y = pct_significant, fill = patient_id)) +
    geom_col(position = position_dodge()) +
    geom_text(aes(label = glue("n={n_clones}")), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    theme_bw() +
    labs(
      title = "% of clones sharing across follicles beyond chance (p<0.05), by clone size",
      subtitle = glue("n_min_cells = {n_min_cells}; p-value is per clone, size-matched (Summary 3)"),
      x = "Clone size (N cells)",
      y = "% significant clones",
      fill = "Patient ID"
    )

  ggsave(glue("{outdir}/pct_significant_by_clone_size_{n_min_cells}.png"), width = 9, height = 6)

  # Mean clonal size VS presents in N follicles
  range_fun <- function(x) {
    data.frame(y = median(x), ymin = min(x), ymax = max(x))
  }
  
  df_plot %>%
    select(patient_id, median_clone_size, n_follicles) %>%
    distinct() %>% 
    ggplot(aes(x = n_follicles, y = median_clone_size, color = patient_id)) +
    geom_point(size = 2.5, alpha = 0.7, position = position_dodge(width = 0.5)) +
    stat_summary(
      data = df_plot,
      aes(y = clone_size),
      fun.data = range_fun,
      geom = "errorbar",
      width = 0.4,
      position = position_dodge(width = 0.5)
    ) + 
    theme_bw() +
    scale_x_continuous(
      breaks = scales::breaks_width(5), 
      minor_breaks = scales::breaks_width(1)
    ) +
    scale_y_break(c(100, 7075), scales = 1/4) +
    scale_y_continuous(
      breaks = scales::breaks_width(10),
      minor_breaks = scales::breaks_width(5)
    ) +
    scale_x_continuous(
      breaks = scales::breaks_width(1),
      minor_breaks = scales::breaks_width(1)
    ) +
    theme(
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank()
    ) + 
    labs(
      color = "Patient ID",
      x = "N follicles",
      y = "Median clone size",
      title = "Median clonal size VS presents in N follicles",
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
    scale_y_continuous(
      labels = scales::percent, 
      breaks = scales::breaks_width(0.1)
    ) +
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

n_min_cells <- 2

# Clones are present in N follicles  
df_n_fol <- df_both %>% 
  filter(
    L1_annotation == "GC_B_cells",
    !is.na(manual_ADT_full_ID),
    !is.na(c_call_grouped)
  ) %>% 
  group_by(clone_subgroup_id_90_similarity, patient_id) %>% 
  count(manual_ADT_full_ID) %>% 
  filter(n >= n_min_cells) %>%
  count(manual_ADT_full_ID) %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  dplyr::rename(n_follicles = n) %>% 
  ungroup()

df_shared_clones <- df_n_fol %>% filter(n_follicles > 1)

# CHECK 
# df_shared_clones %>% filter(patient_id == "HH117" & clone_subgroup_id_90_similarity == "2018_1")

df_both %>% filter(
  L1_annotation == "GC_B_cells",
  !is.na(manual_ADT_full_ID), 
  !is.na(c_call_grouped)
) %>% filter(patient == "HH117" & clone_subgroup_id_90_similarity == "2018_1") %>% count(manual_ADT_ID)


# ------------------------------------------------------------------------------
# Build the pie-chart data for GC B cells
# ------------------------------------------------------------------------------
build_isotype_pie_data <- function(df_both, df_shared_clones, cell_filter = TRUE) {
  
  # df_filtered <- df_both %>% filter({{ cell_filter }})
  df_filtered <- df_both %>% filter(
    L1_annotation == "GC_B_cells",
    !is.na(manual_ADT_full_ID), 
    !is.na(c_call_grouped)
  ) 
  
  # ---- ordered site lookup per patient, plus a trailing "Total" column ----
  site_lookup <- df_filtered %>% 
    distinct(patient_id, manual_ADT_full_ID) %>% 
    mutate(
      # follicle_num = str_extract(sample_clean_fol, "(?<=Fol-)\\d+") %>% as.integer(),
      follicle_num = str_split_i(manual_ADT_full_ID, "-", 2) %>% as.integer(),
      is_follicle = !is.na(follicle_num)
      # sample_clean_fol_plot = sample_clean_fol %>%
      #   str_remove("HH11\\d+-") %>%
      #   str_remove("SI-PP_|SI-PP-nonINF_")
    ) %>% 
    arrange(patient_id, is_follicle, follicle_num, manual_ADT_full_ID) %>% 
    group_by(patient_id) %>% 
    mutate(x = row_number()) %>% 
    ungroup() %>% 
    select(patient_id, manual_ADT_full_ID, x)
  
  total_lookup <- site_lookup %>% 
    group_by(patient_id) %>% 
    summarise(x = max(x) + 1, .groups = "drop") %>% 
    mutate(manual_ADT_full_ID = "Combined")
  
  site_lookup <- bind_rows(site_lookup, total_lookup)
  
  qualifying_follicles <- df_filtered %>% 
    count(patient_id, clone_subgroup_id_90_similarity, manual_ADT_full_ID, name = "n_cells_total") %>% 
    filter(n_cells_total >= n_min_cells)
  
  # ---- counts per clone x site x isotype ----
  counts_by_site <- df_filtered %>% 
    inner_join(df_shared_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
    inner_join(
      qualifying_follicles %>% select(patient_id, clone_subgroup_id_90_similarity, manual_ADT_full_ID), 
      by = c("patient_id", "clone_subgroup_id_90_similarity", "manual_ADT_full_ID")
    ) %>% 
    count(patient_id, clone_subgroup_id_90_similarity, manual_ADT_full_ID, c_call_grouped)
  # counts_by_site <- df_filtered %>% 
  #   inner_join(df_shared_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  #   count(patient_id, clone_subgroup_id_90_similarity, manual_ADT_full_ID, c_call_grouped) 
  
  # ---- counts per clone summed across ALL sites -> the "Total" column ----
  counts_total <- counts_by_site %>% 
    group_by(patient_id, clone_subgroup_id_90_similarity, c_call_grouped) %>% 
    summarise(n = sum(n), .groups = "drop") %>% 
    mutate(manual_ADT_full_ID = "Combined")
  
  plot_data <- bind_rows(counts_by_site, counts_total) %>% 
    left_join(site_lookup, by = c("patient_id", "manual_ADT_full_ID"))
  
  clone_id_lookup <- plot_data %>% 
    distinct(patient_id, clone_subgroup_id_90_similarity) %>% 
    arrange(patient_id, clone_subgroup_id_90_similarity) %>% 
    group_by(patient_id) %>% 
    mutate(y = row_number()) %>% 
    ungroup()
  
  isotype_cols <- plot_data$c_call_grouped %>% unique() %>% na.omit() %>% as.character()
  
  pie_data <- plot_data %>% 
    left_join(clone_id_lookup, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
    pivot_wider(
      id_cols = c(patient_id, clone_subgroup_id_90_similarity, x, y),
      names_from = c_call_grouped, 
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
  
  last_real_x <- x_lookup_sub %>% filter(manual_ADT_full_ID != "Combined") %>% pull(x) %>% max()
  
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
      values = isotype_grouped_colors_custom,
      breaks = c("IGHM/D", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE")
    ) + 
    scale_x_continuous(
      breaks = x_lookup_sub$x, 
      labels = x_lookup_sub$manual_ADT_full_ID,
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
  height <- 8
  
  for (HH in names(patient_names)) {
    
    # dat <- df_both
    # HH <- "HH117"
    # title_suffix <- "bla"
    
    p <- make_pie_plot(HH, dat) + 
      labs(
        title = glue("{HH}: Isotype usage in follicle-shared clones - {title_suffix}"),
        subtitle = glue("Clones shown are present in ≥2 follicles, each with ≥{n_min_cells} GC B cells with an annotated isotype")
      )
    
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
# dat_all <- build_isotype_pie_data(df_both, df_shared_clones, cell_filter = TRUE)
# plot_isotype_pies(dat_all, "all_cells", "all cells")

dat_gc <- build_isotype_pie_data(df_both, df_shared_clones, cell_filter = L1_annotation == "GC_B_cells")
plot_isotype_pies(dat_gc, "GC_B_cells", "GC B cells")


# ------------------------------------------------------------------------------
# Pie grid plot --> Table with % of major isotype
# Rows: Patients
# Columns: % IgA1 in IgA1 dominant clones (mean/medin, sd) (for all isotypes, not sure how to handle 50/50 isotypes.)
# ------------------------------------------------------------------------------

# For each clone, isotype composition combined across all follicles (the
# "Combined" column of the pie grid), find its dominant isotype and what % of
# the clone's cells that dominant isotype makes up. Clones where two or more
# isotypes are tied for the top count (e.g. exactly 50/50) can't be assigned a
# single dominant isotype -- they're labeled "Tied" and reported/excluded
# separately below rather than folded into any one isotype's stats.
df_dominant_isotype <- dat_gc$pie_data %>%
  inner_join(
    dat_gc$site_lookup %>% filter(manual_ADT_full_ID == "Combined") %>% select(patient_id, x),
    by = c("patient_id", "x")
  ) %>%
  pivot_longer(
    cols = all_of(dat_gc$isotype_cols),
    names_to = "c_call_grouped",
    values_to = "n"
  ) %>%
  filter(n > 0) %>%
  group_by(patient_id, clone_subgroup_id_90_similarity) %>%
  mutate(
    pct = n / total_n * 100,
    n_tied = sum(n == max(n))
  ) %>%
  filter(n == max(n)) %>%
  mutate(dominant_isotype = if_else(n_tied == 1, c_call_grouped, "Tied")) %>%
  slice(1) %>%
  ungroup() %>%
  select(patient_id, clone_subgroup_id_90_similarity, total_n, dominant_isotype, dominant_pct = pct)

# Sanity check: how many clones per patient couldn't be assigned a single
# dominant isotype (tied top isotypes)
df_dominant_isotype %>% count(patient_id, is_tied = dominant_isotype == "Tied")

# ------------------------------------------------------------------------------
# Summary table: per patient x dominant isotype, how "pure" that isotype call is
# (i.e. among clones dominated by isotype X, what % of their cells are X)
# ------------------------------------------------------------------------------
table_dominant_isotype <- df_dominant_isotype %>%
  filter(dominant_isotype != "Tied") %>%
  group_by(patient_id, dominant_isotype) %>%
  summarise(
    n_clones   = n(),
    mean_pct   = mean(dominant_pct)   %>% round(1),
    median_pct = median(dominant_pct) %>% round(1),
    sd_pct     = sd(dominant_pct)     %>% round(1),  # NA when n_clones == 1
    .groups = "drop"
  ) %>%
  arrange(patient_id, desc(n_clones))

table_dominant_isotype

table_dominant_isotype_clean <- table_dominant_isotype %>%
  mutate(
    combined_info = glue("{n_clones} ({median_pct}, {mean_pct}, {sd_pct})")
  ) %>% 
  select(patient_id, dominant_isotype, combined_info) %>% 
  pivot_wider(names_from = dominant_isotype, values_from = combined_info) 

# table_dir <- glue("45_immcantation/table/13_clonal_sharing/")
# dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
# write_csv(table_dominant_isotype, glue("{table_dir}/dominant_isotype_purity_by_patient.csv"))


# ==============================================================================
# ==============================================================================
# ==============================================================================


# ------------------------------------------------------------------------------
# Isotype diversity per clone (from the Combined row = summed across follicles)
# ------------------------------------------------------------------------------
isotype_diversity <- function(dat) {
  
  combined_x <- dat$site_lookup %>% 
    filter(manual_ADT_full_ID == "Combined") %>% 
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
      breaks = scales::breaks_width(10),
      minor_breaks = scales::breaks_width(5)
    ) + 
    labs(
      title = glue("{HH}: Clonal sharing degree VS follicle size"),
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
  
  shared_full %>% 
    mutate(
      n_shared_clones = ifelse(n_shared_clones == 0, NA, n_shared_clones)
    ) %>% 
    ggplot(aes(x = follicle_1, y = follicle_2, fill = n_shared_clones)) +
      geom_tile(color = "white") +
      geom_text(aes(label = n_shared_clones), size = 3, color = "blue") +
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
        subtitle = "Grey = 0 shared clones",
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
    # TODO figure out what is follicle 13 and what is follicle 15 
    filter(!(str_detect(follicle_1, "_R|_L") | str_detect(follicle_2, "_R|_L"))) %>% 
    mutate(
      follicle_1 = as.double(follicle_1),
      follicle_2 = as.double(follicle_2),
      f_lo = pmin(follicle_1, follicle_2),
      f_hi = pmax(follicle_1, follicle_2)
    ) %>% 
    select(!c(follicle_1, follicle_2)) 
  
  # Join dataframes of N shared clones and distance
  shared_full_dist <- shared_full_reorder %>% 
    select(!c(follicle_1, follicle_2)) %>% 
    left_join(fol_distances_reorder, by = c("f_lo", "f_hi")) %>% 
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

  # ------------------------------------------------------------------------------
  # Permutation test: does clonal sharing fall off with physical distance?
  # ------------------------------------------------------------------------------
  # Shuffle distances-follicle pairs. 
  # All follciles keep all of their cells and clone in the shuffle.
  
  # Observed information
  obs_df <- shared_full_dist_plot %>%
    filter(!is.na(distance_px)) %>%
    select(f_lo, f_hi, n_shared_clones, distance_px)

  # Distance matrix used for permutations
  dist_df <- obs_df %>% select(f_lo, f_hi, distance_px)

  fols_with_distance <- sort(unique(c(dist_df$f_lo, dist_df$f_hi)))

  # Spearman: n_shared_clones is a skewed, zero-heavy count, and
  # we only care whether the relationship is monotonic, not linear.
  obs_cor <- cor(obs_df$distance_px, obs_df$n_shared_clones, method = "spearman")

  n_perm_dist <- 1000
  set.seed(123)

  perm_cors_dist <- map_dbl(1:n_perm_dist, function(perm_i) {

    relabel <- setNames(sample(fols_with_distance), fols_with_distance)
    
    obs_df %>%
      mutate(
        new_f_lo = relabel[as.character(f_lo)],
        new_f_hi = relabel[as.character(f_hi)],
        f_lo = pmin(new_f_lo, new_f_hi),
        f_hi = pmax(new_f_lo, new_f_hi)
      ) %>%
      select(f_lo, f_hi, n_shared_clones) %>%
      inner_join(dist_df, by = c("f_lo", "f_hi")) %>%
      { cor(.$distance_px, .$n_shared_clones, method = "spearman") }

  })

  df_mantel_test <- tibble(
    patient_id  = HH,
    obs_cor     = obs_cor,
    perm_mean   = mean(perm_cors_dist),
    perm_sd     = sd(perm_cors_dist),
    z_score     = (obs_cor - mean(perm_cors_dist)) / sd(perm_cors_dist),
    # H1: sharing rises with distance (positive correlation)
    p_greater   = mean(perm_cors_dist >= obs_cor),
    # H1: sharing falls with distance (negative correlation) -- the
    # biologically expected direction if nearby follicles cross-seed
    p_less      = mean(perm_cors_dist <= obs_cor),
    p_two_sided = min(1, 2 * min(mean(perm_cors_dist >= obs_cor), mean(perm_cors_dist <= obs_cor))),
    n_pairs     = nrow(obs_df),
    n_perm      = n_perm_dist
  )

  print(df_mantel_test)

  # write_csv(df_mantel_test, glue("{perm_table_dir}/{HH}_distance_sharing_mantel_test.csv"))

  ggplot(tibble(perm_cor = perm_cors_dist), aes(x = perm_cor)) +
    geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = obs_cor, color = "firebrick", linewidth = 1) +
    theme_bw() +
    labs(
      title = glue("{HH}: Distance vs. clonal sharing"),
      subtitle = glue("Observed Spearman r = {round(obs_cor, 3)} (red line); {n_perm_dist} follicle-identity permutations"),
      x = "Spearman correlation (distance vs. N shared clones)",
      y = "N permutations"
    )

  ggsave(glue("{outdir}/{HH}_distance_sharing_mantel_null.png"), width = 8, height = 6)

}



