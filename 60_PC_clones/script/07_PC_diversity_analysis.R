library(glue)
library(tidyverse)
library(alakazam)

source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC_3_definitions", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
patients

HH <- "HH117"
extra <- ""

# Read rds
df_heavy <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>% 
  filter(
    locus == "IGH" & L1_annotation == "PCs" & str_detect(sample_clean, "LP")
  )

# df_heavy$clone_subgroup_id_90_similarity
# df_heavy$manual_ADT_full_ID

# Remove largest clone as it "takes all the signal"
# df_heavy <- df_heavy %>% filter(clone_subgroup_id_90_similarity != "20693_1")
# extra <- "_largest_removed"

# Prep output
outdir = glue("60_PC_clones/plot/07_diversity_analysis/{HH}{extra}")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outdir_combined <- glue("60_PC_clones/plot/07_diversity_analysis/combined{extra}")
dir.create(outdir_combined, recursive = TRUE, showWarnings = FALSE)

# # Prep colors 
# HH_samples <- names(sample_clean_plot_colors) %>% str_subset(glue("^{HH}")) %>% str_subset("Fol")
# HH_samples_colors <- sample_clean_plot_colors[HH_samples]
# names(HH_samples_colors) <- names(HH_samples_colors) %>% str_split_i("_", 2)
# HH_samples_colors

# ------------------------------------------------------------------------------
# Generate a clonal abundance curve
# ------------------------------------------------------------------------------

# Partitions the data on the sample column
# Calculates a 95% confidence interval via 100 bootstrap realizations
set.seed(123) # For reproducibility of example bootstrap results
curve <- estimateAbundance(df_heavy, group="sample_clean", ci=0.95, nboot=100, clone="clone_subgroup_id_90_similarity")

df_heavy$sample_clean %>% unique()

# Plots a rank abundance curve of the relative clonal abundances
plot(curve, legend_title="Sample") + theme_minimal()
ggsave(glue("{outdir}/clonal_abundance_curve.png"))

# ------------------------------------------------------------------------------
# Generate a diversity curve
# ------------------------------------------------------------------------------

# Compare diversity curve across values in the "sample" column
# q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 increments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 100 resampling realizations are performed (nboot=100)
set.seed(123) # For reproducibility of example alphaDiversity results
sample_curve <- alphaDiversity(df_heavy, group="sample_clean", clone="clone_subgroup_id_90_similarity",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=100)

# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
plot(sample_curve, main_title="Sample diversity", 
     legend_title="Patient", shadow=FALSE) + theme_minimal()
ggsave(glue("{outdir}/diversity_curve.png"))


p <- plot(sample_curve, main_title="Sample diversity", 
          legend_title="Patient") + theme_minimal()

p$layers <- p$layers[!sapply(p$layers, function(x) inherits(x$geom, "GeomRibbon"))]
p

# ------------------------------------------------------------------------------
# Shannon diversity
# ------------------------------------------------------------------------------

set.seed(123) # For reproducibility of example alphaDiversity results

for (min_n in c(20, 30, 100)){
  
  # min_n <- 100 # default
  sample_curve <- alphaDiversity(df_heavy, group="sample_clean", clone="clone_subgroup_id_90_similarity",
                                 min_q=0, max_q=4, step_q=0.1,
                                 ci=0.95, nboot=100, min_n = min_n)
  
  # Get shannon values
  shannon_vals <- sample_curve@diversity %>%
    filter(q == 1) 

  ggplot(shannon_vals, aes(x = sample_clean, y = d)) +
    geom_pointrange(aes(ymin = d_lower, ymax = d_upper), color = "steelblue", size = 0.5) +
    theme_minimal() +
    labs(
      x = "Follicle",
      y = "Shannon diversity (q=1)",
      title = glue("{HH}: Shannon diversity per LP PC compartment"),
      subtitle = glue("Clones with <{min_n} cells are excluded from this analysis because of bootstrapping"),
      caption = "Error bars show 95% CI from 100 bootstrap resamples"
    ) 
  
  ggsave(glue("{outdir}/shannon_diversity_plot_min_n_{min_n}.png"), width = 10, height = 6)
  
}


# ------------------------------------------------------------------------------
# View diversity tests at a fixed diversity order
# ------------------------------------------------------------------------------

# # Test diversity at q=0, q=1 and q=2 (equivalent to species richness, Shannon entropy,
# # Simpson's index) across values in the sample_id column
# # 100 bootstrap realizations are performed (nboot=100)
# set.seed(123) # For reproducibility of example alphaDiversity results
# isotype_test <- alphaDiversity(resolve_LC_list_c_clean, group="c_call",
#                                min_q=0, max_q=2, step_q=1, nboot=100, clone="clone_subgroup_id_90_similarity")
# 
# # Print P-value table
# print(isotype_test@tests)
# 
# # Plot results at q=0 and q=2
# # Plot the mean and standard deviations at q=0 and q=2
# plot(isotype_test, 0, colors=isotype_colors_custom, main_title=isotype_main,
#      legend_title="Isotype")
# 
# plot(isotype_test, 2, colors=isotype_colors_custom, main_title=isotype_main,
#      legend_title="Isotype")


# ------------------------------------------------------------------------------
# Dxx plot - GC B cell clones
# ------------------------------------------------------------------------------

compute_Dxx <- function(clone_counts, xx = 0.5) {
  sorted <- sort(clone_counts, decreasing = TRUE)
  cumulative <- cumsum(sorted) / sum(sorted)
  which(cumulative >= xx)[1]
}


d50_per_follicle <- df_heavy %>%
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
  arrange(D50)

d50_per_follicle_long <- d50_per_follicle %>%
  mutate(
    D20_segment = D20,
    D50_extra = D50 - D20
  ) %>%
  pivot_longer(cols = c(D20_segment, D50_extra), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("D50_extra", "D20_segment"))
  ) 

ggplot(d50_per_follicle_long, aes(x = sample_clean, y = value, fill = metric)) +
  geom_col() +
  geom_text(
    data = d50_per_follicle,
    aes(x = sample_clean, y = D50, label = total_clones),
    inherit.aes = FALSE,
    vjust = -0.5, size = 3
  ) +
  scale_fill_manual(
    values = c("D20_segment" = "#E87722", "D50_extra" = "steelblue"),
    labels = c("D20_segment" = "D20", "D50_extra" = "D20 to D50")
  ) +
  theme_minimal() +
  labs(
    x = "Compartment",
    y = "Number of clones",
    fill = NULL,
    title = glue("{HH}: Clonal dominance per LP PC compartment"),
    caption = "Numbers above bars indicate total LP PCs clone count per compartment"
  ) 

ggsave(glue("{outdir}/clonal_D20_D50_plot.png"), width = 10, height = 6)


# ------------------------------------------------------------------------------
# Gini
# ------------------------------------------------------------------------------

gini_coeff <- function(clone_counts) {
  x <- sort(clone_counts)
  n <- length(x)
  numerator <- sum((2 * seq_along(x) - n - 1) * x)
  denominator <- (n - 1) * sum(x)
  numerator / denominator
}
  
gini_per_follicle <- df_heavy %>%
  filter(!is.na(clone_subgroup_id_90_similarity)) %>%
  count(sample_clean, clone_subgroup_id_90_similarity, name = "n_cells") %>%
  group_by(sample_clean) %>%
  summarise(
    gini = gini_coeff(n_cells),
    total_clones = n_distinct(clone_subgroup_id_90_similarity),
    .groups = "drop"
  )

ggplot(gini_per_follicle, 
       aes(x = sample_clean, y = gini)) +
  geom_point(size = 3, color = "steelblue") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    x = "Follicle",
    y = "Gini coefficient",
    title = glue("{HH}: Gini coefficient per LP PC compartment")
  ) 

ggsave(glue("{outdir}/gini_coef.png"), width = 12, height = 6)
  

# ------------------------------------------------------------------------------
# Both patients - Gini + Shannon comparison
# ------------------------------------------------------------------------------

patients <- c("HH117", "HH119")

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>%
    filter(
      locus == "IGH" & L1_annotation == "PCs" & str_detect(sample_clean, "LP")
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows()

if (extra == "_largest_removed"){
  df_both <- df_both %>% filter(clone_subgroup_id_90_similarity != "20693_1")
}


# ------------------------------------------------------------------------------
# Gini - both patients
# ------------------------------------------------------------------------------

gini_combined <- df_both %>%
  filter(!is.na(clone_subgroup_id_90_similarity)) %>%
  count(patient, sample_clean, clone_subgroup_id_90_similarity, name = "n_cells") %>%
  group_by(patient, sample_clean) %>%
  summarise(
    gini = gini_coeff(n_cells),
    total_clones = n_distinct(clone_subgroup_id_90_similarity),
    .groups = "drop"
  )

ggplot(gini_combined, aes(x = sample_clean, y = gini, color = patient)) +
  # geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
  # geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
  geom_point(size = 2.5) + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("HH117" = "#4C72B0", "HH119" = "#DD8452")) +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Gini coefficient",
    title = "Gini coefficient per LP PC compartment",
    subtitle = "Each point represents one compartment"
  ) #+
  # theme(legend.position = "none")

ggsave(glue("{outdir_combined}/gini_coef_combined.png"), width = 10, height = 6)


 # ------------------------------------------------------------------------------
# Shannon - both patients
# ------------------------------------------------------------------------------

# Run alphaDiversity per patient and combine
# lapply(versions, function(version) {
  
for (min_n in c(30, 100)) {
  
  # min_n <- 30
  shannon_combined <- lapply(patients, function(HH) {
    
    df_sha <- df_both %>%
      filter(patient == HH)
     
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
      mutate(patient = HH)
    
  }) %>% bind_rows()
  
  ggplot(shannon_combined, aes(x = sample_clean, y = d, color = patient)) +
    # geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
    # geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("HH117" = "#4C72B0", "HH119" = "#DD8452")) +
    theme_minimal() +
    labs(
      x = "Sample",
      y = "Shannon diversity (q=1)",
      title = glue("Shannon diversity per follicle"),
      subtitle = glue("Follicles with <{min_n} GC B cells excluded. Each point represents one follicle."),
      caption = "Bootstrapped estimates (100 resamples, 95% CI)"
    ) +
    theme(legend.position = "none")
  
  ggsave(glue("{outdir_combined}/shannon_diversity_min_n_{min_n}_combined.png"), width = 10, height = 6)
  
}
  
