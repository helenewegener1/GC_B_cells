library(glue)
library(tidyverse)
library(alakazam)

source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette/
# Adapted from 60_PC_clones/07_PC_diversity_analysis.R.
#
# ==============================================================================
# NOTES / ASSUMPTIONS
# ==============================================================================
# - clone = "CTstrict" instead of clone_subgroup_id_90_similarity (no SHM-based
#   clustering for TCR -- clonal identity is the exact paired TRA+TRB CDR3 nt).
# - 60_PC_clones compares LP (lamina propria) compartments, because that stage
#   is specifically about plasma cell dissemination out of the follicle. For
#   Tfh cells the more natural comparison is across follicles WITHIN SI-PP
#   (Peyer's patches), since that's where Tfh/GC B cell interactions happen.
#   Change SITE_PATTERN below if you want a different comparison (e.g. "LP" to
#   mirror the PC analysis, or "" to include every compartment).
# ==============================================================================

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

df_tcr <- readRDS("70_TCR_analysis/out/df_tcr_filtered.rds")

patients <- df_tcr$patient_id %>% unique()

SI_PP_sites <- df_tcr %>% 
  filter(str_detect(sample_fol_name, "Fol")) %>% 
  pull(sample_clean) %>% 
  unique()

names(SI_PP_sites) <- str_split_i(SI_PP_sites, "-", 1)
SI_PP_sites

# Focus on one sample 
for (HH in patients){
  
  # HH <- "HH117"
  
  df_tcr_HH <- df_tcr %>%
    filter(sample_clean == SI_PP_sites[[HH]])
  
  df_tcr_HH$sample_clean_fol %>% table()
  
  outdir <- glue("70_TCR_analysis/plot/03_diversity_analysis/{HH}")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  outdir_combined <- glue("70_TCR_analysis/plot/03_diversity_analysis/combined")
  dir.create(outdir_combined, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------------------------
  # Clonal abundance curve
  # ------------------------------------------------------------------------------
  
  set.seed(123)
  curve <- estimateAbundance(df_tcr_HH, group = "sample_clean_fol", ci = 0.95, nboot = 100, clone = "CTstrict")
  
  plot(curve, legend_title = "Sample") + theme_minimal()
  ggsave(glue("{outdir}/clonal_abundance_curve.png"))
  
  # ------------------------------------------------------------------------------
  # Diversity (Hill number) curve
  # ------------------------------------------------------------------------------
  
  set.seed(123)
  sample_curve <- alphaDiversity(df_tcr_HH, group = "sample_clean_fol", clone = "CTstrict",
                                 min_q = 0, max_q = 4, step_q = 0.1,
                                 ci = 0.95, nboot = 100)
  
  plot(sample_curve, main_title = glue("{SI_PP_sites[[HH]]}: Tfh cell diversity"),
       legend_title = "Compartment", shadow = FALSE) + theme_minimal()
  ggsave(glue("{outdir}/diversity_curve.png"))
  
  # ------------------------------------------------------------------------------
  # Shannon diversity (q = 1) at a few minimum-cell thresholds
  # ------------------------------------------------------------------------------
  
  set.seed(123)
  
  for (min_n in c(10, 20, 30)) {
    
    # min_n <- 10
    sample_curve <- alphaDiversity(df_tcr_HH, group = "sample_clean_fol", clone = "CTstrict",
                                   min_q = 0, max_q = 4, step_q = 0.1,
                                   ci = 0.95, nboot = 100, min_n = min_n)
    
    shannon_vals <- sample_curve@diversity %>% filter(q == 1) %>% 
      mutate(
        sample_clean_fol_plot = str_split_i(sample_clean_fol, "-", -1) %>% as.integer() %>% as.factor()
      )
    
    ggplot(shannon_vals, aes(x = sample_clean_fol_plot, y = d)) +
      geom_pointrange(aes(ymin = d_lower, ymax = d_upper), color = "steelblue", size = 0.5) +
      theme_minimal() +
      labs(
        x = "Follicle",
        y = "Shannon diversity (q=1)",
        title = glue("{SI_PP_sites[[HH]]}: Shannon diversity"),
        subtitle = glue("Follicles with <{min_n} cells are excluded (bootstrapping)"),
        caption = "Error bars show 95% CI from 100 bootstrap resamples"
      )
    
    ggsave(glue("{outdir}/shannon_diversity_plot_min_n_{min_n}.png"), width = 10, height = 6)
    
  }
  
  # ------------------------------------------------------------------------------
  # D50 / D20 clonal dominance
  # ------------------------------------------------------------------------------
  
  compute_Dxx <- function(clone_counts, xx = 0.5) {
    sorted <- sort(clone_counts, decreasing = TRUE)
    cumulative <- cumsum(sorted) / sum(sorted)
    which(cumulative >= xx)[1]
  }
  
  d50_per_compartment <- df_tcr_HH %>%
    filter(!is.na(CTstrict)) %>%
    count(sample_clean_fol, CTstrict, name = "n_cells") %>%
    group_by(sample_clean_fol) %>%
    summarise(
      D50 = compute_Dxx(n_cells, 0.50),
      D20 = compute_Dxx(n_cells, 0.20),
      total_clones = n_distinct(CTstrict),
      .groups = "drop"
    ) %>%
    arrange(D50) %>% 
    mutate(
      sample_clean_fol_plot = str_split_i(sample_clean_fol, "-", -1) %>% as.integer() %>% as.factor()
    )
  
  d50_per_compartment_long <- d50_per_compartment %>%
    mutate(
      D20_segment = D20,
      D50_extra = D50 - D20
    ) %>%
    pivot_longer(cols = c(D20_segment, D50_extra), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("D50_extra", "D20_segment")))
  
  ggplot(d50_per_compartment_long, aes(x = sample_clean_fol_plot, y = value, fill = metric)) +
    geom_col() +
    geom_text(
      data = d50_per_compartment,
      aes(x = sample_clean_fol_plot, y = D50, label = total_clones),
      inherit.aes = FALSE, vjust = -0.5, size = 3
    ) +
    scale_fill_manual(
      values = c("D20_segment" = "#E87722", "D50_extra" = "steelblue"),
      labels = c("D20_segment" = "D20", "D50_extra" = "D20 to D50")
    ) +
    theme_minimal() +
    labs(
      x = "Follicle",
      y = "Number of clones",
      fill = NULL,
      title = glue("{SI_PP_sites[[HH]]}: TCR clonal dominance per follicle"),
      caption = "Numbers above bars indicate total clone count per follicle"
    )
  
  ggsave(glue("{outdir}/clonal_D20_D50_plot.png"), width = 10, height = 6)
  
  # ------------------------------------------------------------------------------
  # Gini coefficient
  # ------------------------------------------------------------------------------
  
  gini_coeff <- function(clone_counts) {
    x <- sort(clone_counts)
    n <- length(x)
    numerator <- sum((2 * seq_along(x) - n - 1) * x)
    denominator <- (n - 1) * sum(x)
    numerator / denominator
  }
  
  gini_per_compartment <- df_tcr_HH %>%
    filter(!is.na(CTstrict)) %>%
    count(sample_clean_fol, CTstrict, name = "n_cells") %>%
    group_by(sample_clean_fol) %>%
    summarise(
      gini = gini_coeff(n_cells),
      total_clones = n_distinct(CTstrict),
      .groups = "drop"
    ) %>% 
    mutate(
      sample_clean_fol_plot = str_split_i(sample_clean_fol, "-", -1) %>% as.integer() %>% as.factor()
    )
  
  ggplot(gini_per_compartment, aes(x = sample_clean_fol_plot, y = gini)) +
    geom_point(size = 3, color = "steelblue") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    labs(
      x = "Follicle",
      y = "Gini coefficient",
      title = glue("{SI_PP_sites[[HH]]}: Gini coefficient per follcile")
    )
  
  ggsave(glue("{outdir}/gini_coef.png"), width = 12, height = 6)
  
}

# ------------------------------------------------------------------------------
# Combined patients - Gini + Shannon
# ------------------------------------------------------------------------------

gini_combined <- df_tcr %>%
  filter(!is.na(CTstrict)) %>%
  count(patient_id, sample_clean, CTstrict, name = "n_cells") %>%
  group_by(patient_id, sample_clean) %>%
  summarise(
    gini = gini_coeff(n_cells),
    total_clones = n_distinct(CTstrict),
    .groups = "drop"
  )

ggplot(gini_combined, aes(x = sample_clean, y = gini, color = patient_id)) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("HH117" = "#4C72B0", "HH119" = "#DD8452")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Compartment",
    y = "Gini coefficient",
    title = glue("Gini coefficient per compartment of Tfh cells"),
    subtitle = "Each point represents one compartment"
  )

ggsave(glue("{outdir_combined}/gini_coef_combined.png"), width = 10, height = 6)

for (min_n in c(30, 100)) {

  # min_n <- 30
  shannon_combined <- lapply(patients, function(HH_i) {

    df_sha <- df_tcr %>% filter(patient_id == HH_i)

    set.seed(123)
    curve <- alphaDiversity(
      df_sha, group = "sample_clean", clone = "CTstrict",
      min_q = 0, max_q = 4, step_q = 0.1,
      ci = 0.95, nboot = 100, min_n = min_n
    )

    curve@diversity %>% filter(q == 1) %>% mutate(patient_id = HH_i)

  }) %>% bind_rows()

  ggplot(shannon_combined, aes(x = sample_clean, y = d, color = patient_id)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("HH117" = "#4C72B0", "HH119" = "#DD8452")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(
      x = "Compartment",
      y = "Shannon diversity (q=1)",
      title = glue("Shannon diversity per compartment with Tfh cells"),
      subtitle = glue("Compartments with <{min_n} cells excluded."),
      caption = "Bootstrapped estimates (100 resamples, 95% CI)"
    )

  ggsave(glue("{outdir_combined}/shannon_diversity_min_n_{min_n}_combined.png"), width = 10, height = 6)

}
