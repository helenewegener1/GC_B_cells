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

CELL_TYPE_FILTER <- "Tfh_cells"
SITE_PATTERN <- "SI-PP"   # restrict to compartments matching this pattern; set to "" for no filter

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("70_TCR_analysis/out/rds")
patient_files <- grep(glue("_TCR_{CELL_TYPE_FILTER}\\.rds$"), rds_files, value = TRUE) %>%
  grep("^HH", ., value = TRUE)

patients <- lapply(patient_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
patients

HH <- patients[1]

df_tcr <- readRDS(glue("70_TCR_analysis/out/rds/{HH}_TCR_{CELL_TYPE_FILTER}.rds")) %>%
  filter(str_detect(sample_clean, SITE_PATTERN))

outdir <- glue("70_TCR_analysis/plot/03_diversity_analysis/{HH}")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outdir_combined <- glue("70_TCR_analysis/plot/03_diversity_analysis/combined")
dir.create(outdir_combined, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Clonal abundance curve
# ------------------------------------------------------------------------------

set.seed(123)
curve <- estimateAbundance(df_tcr, group = "sample_clean", ci = 0.95, nboot = 100, clone = "CTstrict")

plot(curve, legend_title = "Sample") + theme_minimal()
ggsave(glue("{outdir}/clonal_abundance_curve.png"))

# ------------------------------------------------------------------------------
# Diversity (Hill number) curve
# ------------------------------------------------------------------------------

set.seed(123)
sample_curve <- alphaDiversity(df_tcr, group = "sample_clean", clone = "CTstrict",
                               min_q = 0, max_q = 4, step_q = 0.1,
                               ci = 0.95, nboot = 100)

plot(sample_curve, main_title = glue("{HH}: {CELL_TYPE_FILTER} diversity"),
     legend_title = "Compartment", shadow = FALSE) + theme_minimal()
ggsave(glue("{outdir}/diversity_curve.png"))

# ------------------------------------------------------------------------------
# Shannon diversity (q = 1) at a few minimum-cell thresholds
# ------------------------------------------------------------------------------

set.seed(123)

for (min_n in c(20, 30, 100)) {

  # min_n <- 30
  sample_curve <- alphaDiversity(df_tcr, group = "sample_clean", clone = "CTstrict",
                                 min_q = 0, max_q = 4, step_q = 0.1,
                                 ci = 0.95, nboot = 100, min_n = min_n)

  shannon_vals <- sample_curve@diversity %>% filter(q == 1)

  ggplot(shannon_vals, aes(x = sample_clean, y = d)) +
    geom_pointrange(aes(ymin = d_lower, ymax = d_upper), color = "steelblue", size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Compartment",
      y = "Shannon diversity (q=1)",
      title = glue("{HH}: Shannon diversity per {CELL_TYPE_FILTER} compartment"),
      subtitle = glue("Compartments with <{min_n} cells are excluded (bootstrapping)"),
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

d50_per_compartment <- df_tcr %>%
  filter(!is.na(CTstrict)) %>%
  count(sample_clean, CTstrict, name = "n_cells") %>%
  group_by(sample_clean) %>%
  summarise(
    D50 = compute_Dxx(n_cells, 0.50),
    D20 = compute_Dxx(n_cells, 0.20),
    total_clones = n_distinct(CTstrict),
    .groups = "drop"
  ) %>%
  arrange(D50)

d50_per_compartment_long <- d50_per_compartment %>%
  mutate(
    D20_segment = D20,
    D50_extra = D50 - D20
  ) %>%
  pivot_longer(cols = c(D20_segment, D50_extra), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("D50_extra", "D20_segment")))

ggplot(d50_per_compartment_long, aes(x = sample_clean, y = value, fill = metric)) +
  geom_col() +
  geom_text(
    data = d50_per_compartment,
    aes(x = sample_clean, y = D50, label = total_clones),
    inherit.aes = FALSE, vjust = -0.5, size = 3
  ) +
  scale_fill_manual(
    values = c("D20_segment" = "#E87722", "D50_extra" = "steelblue"),
    labels = c("D20_segment" = "D20", "D50_extra" = "D20 to D50")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Compartment",
    y = "Number of clones",
    fill = NULL,
    title = glue("{HH}: Clonal dominance per {CELL_TYPE_FILTER} compartment"),
    caption = "Numbers above bars indicate total clone count per compartment"
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

gini_per_compartment <- df_tcr %>%
  filter(!is.na(CTstrict)) %>%
  count(sample_clean, CTstrict, name = "n_cells") %>%
  group_by(sample_clean) %>%
  summarise(
    gini = gini_coeff(n_cells),
    total_clones = n_distinct(CTstrict),
    .groups = "drop"
  )

ggplot(gini_per_compartment, aes(x = sample_clean, y = gini)) +
  geom_point(size = 3, color = "steelblue") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Compartment",
    y = "Gini coefficient",
    title = glue("{HH}: Gini coefficient per {CELL_TYPE_FILTER} compartment")
  )

ggsave(glue("{outdir}/gini_coef.png"), width = 12, height = 6)

# ------------------------------------------------------------------------------
# Combined patients - Gini + Shannon
# ------------------------------------------------------------------------------

df_both <- lapply(patients, function(x) {
  readRDS(glue("70_TCR_analysis/out/rds/{x}_TCR_{CELL_TYPE_FILTER}.rds")) %>%
    filter(str_detect(sample_clean, SITE_PATTERN))
}) %>% bind_rows()

gini_combined <- df_both %>%
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
    title = glue("Gini coefficient per {CELL_TYPE_FILTER} compartment"),
    subtitle = "Each point represents one compartment"
  )

ggsave(glue("{outdir_combined}/gini_coef_combined.png"), width = 10, height = 6)

for (min_n in c(30, 100)) {

  # min_n <- 30
  shannon_combined <- lapply(patients, function(HH_i) {

    df_sha <- df_both %>% filter(patient_id == HH_i)

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
      title = glue("Shannon diversity per {CELL_TYPE_FILTER} compartment"),
      subtitle = glue("Compartments with <{min_n} cells excluded."),
      caption = "Bootstrapped estimates (100 resamples, 95% CI)"
    )

  ggsave(glue("{outdir_combined}/shannon_diversity_min_n_{min_n}_combined.png"), width = 10, height = 6)

}
