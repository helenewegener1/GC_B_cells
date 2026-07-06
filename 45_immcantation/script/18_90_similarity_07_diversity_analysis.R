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

HH <- "HH119"
extra <- ""

# Read rds
df_heavy <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>% 
  filter(
    locus == "IGH",
    !is.na(manual_ADT_full_ID)
  )

# df_heavy$clone_subgroup_id_90_similarity
# df_heavy$manual_ADT_full_ID

# Remove largest clone as it "takes all the signal"
# df_heavy <- df_heavy %>% filter(clone_subgroup_id_90_similarity != "20693_1")
# extra <- "_largest_removed"

# Prep output
outdir = glue("45_immcantation/plot/18_90_similarity/07_diversity_analysis/{HH}{extra}")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outdir_combined <- glue("45_immcantation/plot/18_90_similarity/07_diversity_analysis/combined{extra}")
dir.create(outdir_combined, recursive = TRUE, showWarnings = FALSE)


# Prep colors 
HH_samples <- names(sample_clean_plot_colors) %>% str_subset(glue("^{HH}")) %>% str_subset("Fol")
HH_samples_colors <- sample_clean_plot_colors[HH_samples]
names(HH_samples_colors) <- names(HH_samples_colors) %>% str_split_i("_", 2)
HH_samples_colors

# ------------------------------------------------------------------------------
# Generate a clonal abundance curve
# ------------------------------------------------------------------------------

# Partitions the data on the sample column
# Calculates a 95% confidence interval via 100 bootstrap realizations
set.seed(123) # For reproducibility of example bootstrap results
curve <- estimateAbundance(df_heavy, group="manual_ADT_full_ID", ci=0.95, nboot=100, clone="clone_subgroup_id_90_similarity")

df_heavy$manual_ADT_full_ID %>% unique()

# Plots a rank abundance curve of the relative clonal abundances
plot(curve, colors = HH_samples_colors, legend_title="Sample") + theme_minimal()
ggsave(glue("{outdir}/clonal_abundance_curve.png"))

# ------------------------------------------------------------------------------
# Generate a diversity curve
# ------------------------------------------------------------------------------

# Compare diversity curve across values in the "sample" column
# q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 increments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 100 resampling realizations are performed (nboot=100)
set.seed(123) # For reproducibility of example alphaDiversity results
sample_curve <- alphaDiversity(df_heavy, group="manual_ADT_full_ID", clone="clone_subgroup_id_90_similarity",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=100)

# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
plot(sample_curve, colors=HH_samples_colors, main_title="Sample diversity", 
     legend_title="Patient", shadow=FALSE) + theme_minimal()
ggsave(glue("{outdir}/diversity_curve.png"))


p <- plot(sample_curve, colors=HH_samples_colors, main_title="Sample diversity", 
          legend_title="Patient") + theme_minimal()

p$layers <- p$layers[!sapply(p$layers, function(x) inherits(x$geom, "GeomRibbon"))]
p

# ------------------------------------------------------------------------------
# Shannon diversity
# ------------------------------------------------------------------------------

set.seed(123) # For reproducibility of example alphaDiversity results

versions <- c("GC_B_cells", "all_cells")

lapply(versions, function(version) {
  
  # version <- "GC_B_cells"
  df_heavy_sha <- if (version == "GC_B_cells") {
    df_heavy %>% filter(L1_annotation == "GC_B_cells")
  } else {
    df_heavy
  }
  
  # N cells per follicle 
  # df_heavy_sha %>% group_by(manual_ADT_full_ID) %>% count() %>% view()
  
  for (min_n in c(20, 30)){
    
    # min_n <- 30 # default
    sample_curve <- alphaDiversity(df_heavy_sha, group="manual_ADT_full_ID", clone="clone_subgroup_id_90_similarity",
                                   min_q=0, max_q=4, step_q=0.1,
                                   ci=0.95, nboot=100, min_n = min_n)
    
    # Get shannon values
    shannon_vals <- sample_curve@diversity %>%
      filter(q == 1) %>%
      mutate(
        manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer(),
        manual_ADT_full_ID_plot = factor(manual_ADT_full_ID_plot, levels = 1:max(manual_ADT_full_ID_plot))
      )
    
    version_txt <- str_replace_all(version, "_", " ")
    
    ggplot(shannon_vals, aes(x = manual_ADT_full_ID_plot, y = d)) +
      geom_pointrange(aes(ymin = d_lower, ymax = d_upper), color = "steelblue", size = 0.5) +
      theme_minimal() +
      labs(
        x = "Follicle",
        y = "Shannon diversity (q=1)",
        title = glue("{HH}: {version_txt} Shannon diversity per follicle"),
        subtitle = glue("Clones with <{min_n} cells are excluded from this analysis because of bootstrapping"),
        caption = "Error bars show 95% CI from 100 bootstrap resamples"
      ) 
    
    ggsave(glue("{outdir}/shannon_diversity_plot_{version}_min_n_{min_n}.png"), width = 10, height = 6)
    
  }
  
})
  

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
# GC B cells count
# ------------------------------------------------------------------------------

df_heavy %>% 
  filter(L1_annotation == "GC_B_cells") %>% 
  count(manual_ADT_full_ID) %>% 
  mutate(manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer() %>% as.factor()) %>% 
  ggplot(aes(x = manual_ADT_full_ID_plot, y = n)) + 
  geom_col() + 
  geom_text(
    aes(x = manual_ADT_full_ID_plot, y = n, label = n),
    inherit.aes = FALSE,
    vjust = -0.5, size = 3
  ) + 
  theme_minimal() + 
  labs(
    title = glue("{HH}: GC B cell (with BCR) count")
  )

ggsave(glue("{outdir}/GC_B_cell_count.png"), width = 12, height = 6)

# ------------------------------------------------------------------------------
# Dxx plot - GC B cell clones
# ------------------------------------------------------------------------------

compute_Dxx <- function(clone_counts, xx = 0.5) {
  sorted <- sort(clone_counts, decreasing = TRUE)
  cumulative <- cumsum(sorted) / sum(sorted)
  which(cumulative >= xx)[1]
}

versions <- c("GC_B_cells", "all_cells")

lapply(versions, function(version) {
  
  # version <- "GC_B_cells"
  
  df_heavy_DXX <- if (version == "GC_B_cells") {
    df_heavy %>% filter(L1_annotation == "GC_B_cells")
  } else {
    df_heavy
  }
  
  d50_per_follicle <- df_heavy_DXX %>%
    filter(
      !is.na(clone_subgroup_id_90_similarity)
    ) %>%
    count(manual_ADT_full_ID, clone_subgroup_id_90_similarity, name = "n_cells") %>%
    group_by(manual_ADT_full_ID) %>%
    summarise(
      D50 = compute_Dxx(n_cells, 0.50),
      D20 = compute_Dxx(n_cells, 0.20),
      total_clones = n_distinct(clone_subgroup_id_90_similarity),
      .groups = "drop"
    ) %>%
    arrange(D50) %>% 
    mutate(
      manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer() %>% as.factor()
    )
  
  d50_per_follicle_long <- d50_per_follicle %>%
    mutate(
      D20_segment = D20,
      D50_extra = D50 - D20
    ) %>%
    pivot_longer(cols = c(D20_segment, D50_extra), names_to = "metric", values_to = "value") %>%
    mutate(
      metric = factor(metric, levels = c("D50_extra", "D20_segment"))
    ) 
  
  version_txt <- str_replace_all(version, "_", " ")
  
  ggplot(d50_per_follicle_long, aes(x = manual_ADT_full_ID_plot, y = value, fill = metric)) +
    geom_col() +
    geom_text(
      data = d50_per_follicle,
      aes(x = manual_ADT_full_ID_plot, y = D50, label = total_clones),
      inherit.aes = FALSE,
      vjust = -0.5, size = 3
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
      title = glue("{HH}: {version_txt} clonal dominance per follicle"),
      caption = "Numbers above bars indicate total GC B cell clone count per follicle"
    ) 
  
  ggsave(glue("{outdir}/clonal_D20_D50_plot_{version}.png"), width = 10, height = 6)
  
})

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

versions <- c("GC_B_cells", "all_cells")

lapply(versions, function(version) {
  
  # version <- "GC_B_cells"
  
  df_heavy_DXX <- if (version == "GC_B_cells") {
    df_heavy %>% filter(L1_annotation == "GC_B_cells")
  } else {
    df_heavy
  }
  
  gini_per_follicle <- df_heavy_DXX %>%
    filter(!is.na(clone_subgroup_id_90_similarity)) %>%
    count(manual_ADT_full_ID, clone_subgroup_id_90_similarity, name = "n_cells") %>%
    group_by(manual_ADT_full_ID) %>%
    summarise(
      gini = gini_coeff(n_cells),
      total_clones = n_distinct(clone_subgroup_id_90_similarity),
      .groups = "drop"
    )
  
  version_txt <- str_replace_all(version, "_", " ")
  
  ggplot(gini_per_follicle, aes(x = factor(
    str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer(),
    levels = 1:max(str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer())
  ), y = gini)) +
    geom_point(size = 3, color = "steelblue") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    labs(
      x = "Follicle",
      y = "Gini coefficient",
      title = glue("{HH}: {version_txt} Gini coefficient per follicle")
    ) 
  
  ggsave(glue("{outdir}/gini_coef_{version}.png"), width = 12, height = 6)
  
})  

# ------------------------------------------------------------------------------
# Both patients - Gini + Shannon comparison
# ------------------------------------------------------------------------------

patients <- c("HH117", "HH119")

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>%
    filter(
      locus == "IGH",
      !is.na(manual_ADT_full_ID)
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows()

if (extra == "_largest_removed"){
  df_both <- df_both %>% filter(clone_subgroup_id_90_similarity != "20693_1")
}

versions <- c("GC_B_cells", "all_cells")

# ------------------------------------------------------------------------------
# Gini - both patients
# ------------------------------------------------------------------------------

lapply(versions, function(version) {
  
  df_gini <- if (version == "GC_B_cells") {
    df_both %>% filter(L1_annotation == "GC_B_cells")
  } else {
    df_both
  }
  
  gini_combined <- df_gini %>%
    filter(!is.na(clone_subgroup_id_90_similarity)) %>%
    count(patient, manual_ADT_full_ID, clone_subgroup_id_90_similarity, name = "n_cells") %>%
    group_by(patient, manual_ADT_full_ID) %>%
    summarise(
      gini = gini_coeff(n_cells),
      total_clones = n_distinct(clone_subgroup_id_90_similarity),
      .groups = "drop"
    )
  
  version_txt <- str_replace_all(version, "_", " ")
  
  ggplot(gini_combined, aes(x = patient, y = gini)) +
    geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
    geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("HH117" = "#4C72B0", "HH119" = "#DD8452")) +
    theme_minimal() +
    labs(
      x = "Patient",
      y = "Gini coefficient",
      title = glue("Gini coefficient per follicle - {version_txt}"),
      subtitle = "Each point represents one follicle"
    ) +
    theme(legend.position = "none")
  
  ggsave(glue("{outdir_combined}/gini_coef_{version}_combined.png"), width = 8, height = 8)
  
})

# ------------------------------------------------------------------------------
# Shannon - both patients
# ------------------------------------------------------------------------------

# Run alphaDiversity per patient and combine
lapply(versions, function(version) {
  
  for (min_n in c(30)) {
    
    shannon_combined <- lapply(patients, function(HH) {
      
      df_sha <- df_both %>%
        filter(patient == HH) %>%
        { if (version == "GC_B_cells") filter(., L1_annotation == "GC_B_cells") else . }
      
      set.seed(123)
      curve <- alphaDiversity(
        df_sha,
        group = "manual_ADT_full_ID",
        clone = "clone_subgroup_id_90_similarity",
        min_q = 0, max_q = 4, step_q = 0.1,
        ci = 0.95, nboot = 100, min_n = min_n
      )
      
      curve@diversity %>%
        filter(q == 1) %>%
        mutate(patient = HH)
      
    }) %>% bind_rows()
    
    version_txt <- str_replace_all(version, "_", " ")
    
    ggplot(shannon_combined, aes(x = patient, y = d)) +
      geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
      geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
      scale_color_manual(values = c("HH117" = "#4C72B0", "HH119" = "#DD8452")) +
      theme_minimal() +
      labs(
        x = "Patient",
        y = "Shannon diversity (q=1)",
        title = glue("Shannon diversity per follicle - {version_txt}"),
        subtitle = glue("Follicles with <{min_n} GC B cells excluded. Each point represents one follicle."),
        caption = "Bootstrapped estimates (100 resamples, 95% CI)"
      ) +
      theme(legend.position = "none")
    
    ggsave(glue("{outdir_combined}/shannon_diversity_{version}_min_n_{min_n}_combined.png"), width = 8, height = 8)
    
  }
  
})

# ------------------------------------------------------------------------------
# Testing if: Intra-follicular GC B cell clones can exhibit greater sequence homology than their inter-follicular counterparts
# ------------------------------------------------------------------------------

library(tidyverse)
library(stringdist)

# Define clone
clone_nr <- 3
clone <- df_heavy %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  slice(clone_nr) %>% 
  pull(clone_subgroup_id_90_similarity)

clone_of_interest <- df_heavy %>% 
  filter(clone_subgroup_id_90_similarity == clone)

seqs <- clone_of_interest$junction
names(seqs) <- clone_of_interest$cell_id

nchar(seqs) %>% unique

dist_mat <- stringdistmatrix(seqs, seqs, method = "hamming")
rownames(dist_mat) <- names(seqs)
colnames(dist_mat) <- names(seqs)

pairs_df <- as_tibble(as.table(dist_mat), .name_repair = "minimal") %>% 
  set_names(c("cell1", "cell2", "distance")) %>% 
  filter(cell1 != cell2) %>% 
  left_join(clone_of_interest %>% select(cell_id, follicle1 = manual_ADT_ID), by = c("cell1" = "cell_id")) %>% 
  left_join(clone_of_interest %>% select(cell_id, follicle2 = manual_ADT_ID), by = c("cell2" = "cell_id")) %>% 
  mutate(
    pair_key = map2_chr(cell1, cell2, ~ paste(sort(c(.x, .y)), collapse = "_"))
  ) %>% 
  distinct(pair_key, .keep_all = TRUE) %>% 
  mutate(pair_type = if_else(follicle1 == follicle2, "intra", "inter"))

pairs_df %>% 
  group_by(pair_type) %>% 
  summarise(
    mean_dist = mean(distance),
    median_dist = median(distance),
    n_pairs = n()
  )


########

library(tidyverse)
library(stringdist)
library(lmerTest)
library(patchwork)

n_permutations <- 100
set.seed(1)

# ---- 1. nest data per clone, keep clones with enough cells, >=2 follicles, and at least one repeated follicle ----

clone_data <- df_both %>% 
  group_by(patient_id, clone_subgroup_id_90_similarity) %>% 
  filter(n() >= 4, n_distinct(manual_ADT_ID) >= 2) %>% 
  nest() %>% 
  ungroup() %>% 
  rename(clone_id = clone_subgroup_id_90_similarity) %>% 
  mutate(n_cells = map_int(data, nrow)) %>% 
  filter(map_lgl(data, ~ any(duplicated(.x$manual_ADT_ID))))

# ---- 2. per clone: distance matrix, upper-triangle pairs, observed diff, permutation p-value ----

clone_results <- clone_data %>% 
  mutate(
    dist_mat = map(data, ~ {
      seqs <- .x$junction
      names(seqs) <- .x$cell_id
      m <- as.matrix(stringdistmatrix(seqs, seqs, method = "lv"))
      dimnames(m) <- list(.x$cell_id, .x$cell_id)
      m
    }),
    follicle_vec = map(data, ~ .x$manual_ADT_ID),
    ut_idx = map(dist_mat, ~ which(upper.tri(.x), arr.ind = TRUE)),
    pairs_df = pmap(list(dist_mat, follicle_vec, ut_idx), ~ {
      tibble(
        distance = ..1[..3],
        follicle1 = ..2[..3[, "row"]],
        follicle2 = ..2[..3[, "col"]]
      ) %>% 
        mutate(pair_type = if_else(follicle1 == follicle2, "intra", "inter"))
    })
  ) %>% 
  mutate(
    observed_diff = map_dbl(pairs_df, ~ {
      means <- .x %>% group_by(pair_type) %>% summarise(m = mean(distance)) %>% deframe()
      if (!all(c("intra", "inter") %in% names(means))) return(NA_real_)
      means[["inter"]] - means[["intra"]]
    }),
    perm_diffs = pmap(list(dist_mat, follicle_vec, ut_idx), ~ {
      m <- ..1; foll <- ..2; idx <- ..3
      d <- m[idx]
      map_dbl(1:n_permutations, ~ {
        foll_perm <- sample(foll)
        same <- foll_perm[idx[, "row"]] == foll_perm[idx[, "col"]]
        if (sum(same) == 0 || sum(!same) == 0) return(NA_real_)
        mean(d[!same]) - mean(d[same])
      })
    }),
    p_value = map2_dbl(observed_diff, perm_diffs, ~ {
      if (is.na(.x)) return(NA_real_)
      mean(.y >= .x, na.rm = TRUE)
    })
  ) %>% 
  select(-dist_mat, -follicle_vec, -ut_idx)

clone_summary <- clone_results %>% 
  select(patient_id, clone_id, n_cells, observed_diff, p_value) %>% 
  arrange(patient_id, desc(n_cells))

# ---- 3. pool all pairs for an overall mixed-effects estimate ----

all_pairs <- clone_results %>% 
  select(patient_id, clone_id, pairs_df) %>% 
  unnest(pairs_df)

mixed_model <- lmer(distance ~ pair_type + (1 | patient_id/clone_id), data = all_pairs)
summary(mixed_model)

# ---- 4. plot: per-clone effects + pooled distributions ----

p_clone_effects <- clone_summary %>% 
  filter(!is.na(observed_diff)) %>% 
  mutate(clone_id = fct_reorder(clone_id, observed_diff)) %>% 
  ggplot(aes(x = clone_id, y = observed_diff)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n_cells, color = p_value < 0.05)) + 
  coord_flip() + 
  facet_wrap(~patient_id, scales = "free_y") + 
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey60"), name = "p < 0.05") +
  labs(
    x = "Clone",
    y = "Inter minus intra follicular\nmean edit distance",
    size = "N cells",
    title = "Per-clone homology effect"
  ) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_overall <- all_pairs %>% 
  ggplot(aes(x = pair_type, y = distance, fill = pair_type)) + 
  geom_violin(alpha = 0.6) + 
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") + 
  facet_wrap(~patient_id) + 
  scale_fill_manual(values = c(intra = "#1b9e77", inter = "#d95f02")) + 
  labs(
    x = NULL,
    y = "Pairwise edit distance",
    title = "Intra vs inter follicular sequence distance",
    subtitle = "Pooled across all clones per patient"
  ) + 
  theme_bw() + 
  theme(legend.position = "none")

p_overall + p_clone_effects + plot_layout(widths = c(1, 1.4))