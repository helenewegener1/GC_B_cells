library(glue)
library(tidyverse)
library(alakazam)

source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

# patient_to_condition
patient_to_condition <- data.frame(
  patient = c("HH117", "HH151", "HH153", "HH119"), 
  condition = c(rep("Crohn's", 3), rep("Control", 1))
)

# for (HH in patients){
#   
#   # HH <- "HH153"
#   p <- patient_names[[HH]]
#   
#   
#   # Read rds
#   df_heavy <- readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>% 
#     filter(
#       locus == "IGH",
#       !is.na(manual_ADT_full_ID)
#     )
#   
#   # df_heavy$clone_subgroup_id_90_similarity
#   # df_heavy$manual_ADT_full_ID
#   
#   # Remove largest clone as it "takes all the signal"
#   # largest_clone <- df_heavy %>% count(clone_subgroup_id_90_similarity, sort = TRUE) %>% head(1) %>% pull(clone_subgroup_id_90_similarity)
#   # df_heavy <- df_heavy %>% filter(clone_subgroup_id_90_similarity != largest_clone)
#   # extra <- "_largest_removed"
#   
#   # Prep output
#   outdir = glue("45_immcantation/plot/12_diversity_analysis/{HH}{extra}")
#   dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#   
#   outdir_combined <- glue("45_immcantation/plot/12_diversity_analysis/combined{extra}")
#   dir.create(outdir_combined, recursive = TRUE, showWarnings = FALSE)
#   
#   # Prep colors 
#   HH_samples <- names(sample_clean_plot_colors) %>% str_subset(glue("^{HH}")) %>% str_subset("Fol")
#   HH_samples_colors <- sample_clean_plot_colors[HH_samples]
#   names(HH_samples_colors) <- names(HH_samples_colors) %>% str_split_i("_", 2)
#   HH_samples_colors
#   
#   
#   # # ------------------------------------------------------------------------------
#   # # Generate a clonal abundance curve
#   # # ------------------------------------------------------------------------------
#   # 
#   # # Partitions the data on the sample column
#   # # Calculates a 95% confidence interval via 100 bootstrap realizations
#   # set.seed(123) # For reproducibility of example bootstrap results
#   # curve <- estimateAbundance(df_heavy, group="manual_ADT_full_ID", ci=0.95, nboot=100, clone="clone_subgroup_id_90_similarity")
#   # 
#   # df_heavy$manual_ADT_full_ID %>% unique()
#   # 
#   # # Plots a rank abundance curve of the relative clonal abundances
#   # plot(curve, colors = HH_samples_colors, legend_title="Sample") + theme_minimal()
#   # ggsave(glue("{outdir}/clonal_abundance_curve.png"))
#   # 
#   # # ------------------------------------------------------------------------------
#   # # Generate a diversity curve
#   # # ------------------------------------------------------------------------------
#   # 
#   # # Compare diversity curve across values in the "sample" column
#   # # q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 increments (step_q=0.05)
#   # # A 95% confidence interval will be calculated (ci=0.95)
#   # # 100 resampling realizations are performed (nboot=100)
#   # set.seed(123) # For reproducibility of example alphaDiversity results
#   # sample_curve <- alphaDiversity(df_heavy, group="manual_ADT_full_ID", clone="clone_subgroup_id_90_similarity",
#   #                                min_q=0, max_q=4, step_q=0.1,
#   #                                ci=0.95, nboot=100)
#   # 
#   # # Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
#   # # Indicate number of sequences resampled from each group in the title
#   # plot(sample_curve, colors=HH_samples_colors, main_title="Sample diversity", 
#   #      legend_title="Patient", shadow=FALSE) + theme_minimal()
#   # ggsave(glue("{outdir}/diversity_curve.png"))
#   # 
#   # 
#   # p <- plot(sample_curve, colors=HH_samples_colors, main_title="Sample diversity", 
#   #           legend_title="Patient") + theme_minimal()
#   # 
#   # p$layers <- p$layers[!sapply(p$layers, function(x) inherits(x$geom, "GeomRibbon"))]
#   # p
#   # 
#   # # ------------------------------------------------------------------------------
#   # # Shannon diversity
#   # # ------------------------------------------------------------------------------
#   # 
#   # set.seed(123) # For reproducibility of example alphaDiversity results
#   # 
#   # versions <- c("GC_B_cells", "all_cells")
#   # 
#   # lapply(versions, function(version) {
#   #   
#   #   # version <- "GC_B_cells"
#   #   df_heavy_sha <- if (version == "GC_B_cells") {
#   #     df_heavy %>% filter(L1_annotation == "GC_B_cells")
#   #   } else {
#   #     df_heavy
#   #   }
#   #   
#   #   # N cells per follicle 
#   #   # df_heavy_sha %>% group_by(manual_ADT_full_ID) %>% count() %>% view()
#   #   
#   #   for (min_n in c(20, 30)){
#   #     
#   #     # min_n <- 30 # default
#   #     sample_curve <- alphaDiversity(df_heavy_sha, group="manual_ADT_full_ID", clone="clone_subgroup_id_90_similarity",
#   #                                    min_q=0, max_q=4, step_q=0.1,
#   #                                    ci=0.95, nboot=100, min_n = min_n)
#   #     
#   #     # Get shannon values
#   #     shannon_vals <- sample_curve@diversity %>%
#   #       filter(q == 1) %>%
#   #       mutate(
#   #         manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer(),
#   #         manual_ADT_full_ID_plot = factor(manual_ADT_full_ID_plot, levels = 1:max(manual_ADT_full_ID_plot))
#   #       )
#   #     
#   #     version_txt <- str_replace_all(version, "_", " ")
#   #     
#   #     ggplot(shannon_vals, aes(x = manual_ADT_full_ID_plot, y = d)) +
#   #       geom_pointrange(aes(ymin = d_lower, ymax = d_upper), color = "steelblue", size = 0.5) +
#   #       theme_minimal() +
#   #       labs(
#   #         x = "Follicle",
#   #         y = "Shannon diversity (q=1)",
#   #         title = glue("{HH}: {version_txt} Shannon diversity per follicle"),
#   #         subtitle = glue("Clones with <{min_n} cells are excluded from this analysis because of bootstrapping"),
#   #         caption = "Error bars show 95% CI from 100 bootstrap resamples"
#   #       ) 
#   #     
#   #     ggsave(glue("{outdir}/shannon_diversity_plot_{version}_min_n_{min_n}.png"), width = 10, height = 6)
#   #     
#   #   }
#   #   
#   # })
#   #   
#   
#   # ------------------------------------------------------------------------------
#   # View diversity tests at a fixed diversity order
#   # ------------------------------------------------------------------------------
#   
#   # # Test diversity at q=0, q=1 and q=2 (equivalent to species richness, Shannon entropy,
#   # # Simpson's index) across values in the sample_id column
#   # # 100 bootstrap realizations are performed (nboot=100)
#   # set.seed(123) # For reproducibility of example alphaDiversity results
#   # isotype_test <- alphaDiversity(resolve_LC_list_c_clean, group="c_call",
#   #                                min_q=0, max_q=2, step_q=1, nboot=100, clone="clone_subgroup_id_90_similarity")
#   # 
#   # # Print P-value table
#   # print(isotype_test@tests)
#   # 
#   # # Plot results at q=0 and q=2
#   # # Plot the mean and standard deviations at q=0 and q=2
#   # plot(isotype_test, 0, colors=isotype_colors_custom, main_title=isotype_main,
#   #      legend_title="Isotype")
#   # 
#   # plot(isotype_test, 2, colors=isotype_colors_custom, main_title=isotype_main,
#   #      legend_title="Isotype")
#   
#   # ------------------------------------------------------------------------------
#   # GC B cells count
#   # ------------------------------------------------------------------------------
#   
#   df_heavy %>% 
#     filter(L1_annotation == "GC_B_cells") %>% 
#     dplyr::count(manual_ADT_full_ID) %>% 
#     mutate(manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer() %>% as.factor()) %>% 
#     ggplot(aes(x = manual_ADT_full_ID_plot, y = n)) + 
#     geom_col() + 
#     geom_text(
#       aes(x = manual_ADT_full_ID_plot, y = n, label = n),
#       inherit.aes = FALSE,
#       vjust = -0.5, size = 3
#     ) + 
#     theme_minimal() + 
#     labs(
#       title = glue("{HH}: GC B cell (with BCR) count")
#     )
#   
#   ggsave(glue("{outdir}/GC_B_cell_count.png"), width = 12, height = 6)
#   
#   # ------------------------------------------------------------------------------
#   # Dxx plot - GC B cell clones
#   # ------------------------------------------------------------------------------
#   
#   compute_Dxx <- function(clone_counts, xx = 0.5) {
#     sorted <- sort(clone_counts, decreasing = TRUE)
#     cumulative <- cumsum(sorted) / sum(sorted)
#     which(cumulative >= xx)[1]
#   }
#   
#   # versions <- c("GC_B_cells", "all_cells")
#   
#   # lapply(versions, function(version) {
# 
#     version <- "GC_B_cells"
#     
#     df_heavy_DXX <- if (version == "GC_B_cells") {
#       df_heavy %>% filter(L1_annotation == "GC_B_cells")
#     } else {
#       df_heavy
#     }
#     
#     # Remove follicles that have less than 5 GC B cells
#     fol_to_rm <- df_heavy_DXX %>% 
#       count(manual_ADT_ID) %>% 
#       filter(n < 5) %>% 
#       pull(manual_ADT_ID)
#     
#     df_heavy_DXX <- df_heavy_DXX %>% 
#       filter(!(manual_ADT_ID %in% fol_to_rm))
#     
#     d50_per_follicle <- df_heavy_DXX %>%
#       filter(
#         !is.na(clone_subgroup_id_90_similarity)
#       ) %>%
#       dplyr::count(manual_ADT_full_ID, clone_subgroup_id_90_similarity, name = "n_cells") %>%
#       group_by(manual_ADT_full_ID) %>%
#       summarise(
#         D50 = compute_Dxx(n_cells, 0.50),
#         D20 = compute_Dxx(n_cells, 0.20),
#         total_clones = n_distinct(clone_subgroup_id_90_similarity),
#         .groups = "drop"
#       ) %>%
#       arrange(D50) %>% 
#       mutate(
#         manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer() %>% as.factor()
#       )
#     
#     d50_per_follicle_long <- d50_per_follicle %>%
#       mutate(
#         D20_segment = D20,
#         D50_extra = D50 - D20
#       ) %>%
#       pivot_longer(cols = c(D20_segment, D50_extra), names_to = "metric", values_to = "value") %>%
#       mutate(
#         metric = factor(metric, levels = c("D50_extra", "D20_segment"))
#       ) 
#     
#     version_txt <- str_replace_all(version, "_", " ")
#     
#     # D20 + D50
#     ggplot(d50_per_follicle_long, aes(x = manual_ADT_full_ID_plot, y = value, fill = metric)) +
#       geom_col() +
#       geom_text(
#         data = d50_per_follicle,
#         aes(x = manual_ADT_full_ID_plot, y = D50, label = total_clones),
#         inherit.aes = FALSE,
#         vjust = -0.5, size = 3
#       ) +
#       scale_fill_manual(
#         values = c("D20_segment" = L1_colors[["GC_B_cells"]], "D50_extra" = "forestgreen"),
#         labels = c("D20_segment" = "D20", "D50_extra" = "D20 to D50")
#       ) +
#       scale_y_continuous(
#         breaks = scales::breaks_width(1)
#         # minor_breaks = scales::breaks_width(minor_breaks_width),
#       ) + 
#       theme_minimal() +
#       labs(
#         x = "Follicle",
#         y = "Number of clones",
#         fill = NULL,
#         title = glue("{HH}: {version_txt} clonal dominance per follicle"),
#         caption = "Numbers above bars indicate total GC B cell clone count per follicle."
#       ) 
#     
#     ggsave(glue("{outdir}/clonal_D20_D50_plot_{version}.png"), width = 10, height = 6)
#     
#     # D20
#     d50_per_follicle_long %>% 
#       filter(metric == "D20_segment") %>% 
#       ggplot(aes(x = manual_ADT_full_ID_plot, y = value, fill = metric)) +
#       geom_col() +
#       geom_text(
#         data = d50_per_follicle,
#         aes(x = manual_ADT_full_ID_plot, y = D20, label = total_clones),
#         inherit.aes = FALSE,
#         vjust = -0.5, size = 4
#       ) +
#       scale_fill_manual(
#         values = c("D20_segment" = L1_colors[["GC_B_cells"]]), 
#         labels = c("D20_segment" = "D20")
#       ) +
#       scale_y_continuous(
#         breaks = scales::breaks_width(1), 
#         minor_breaks = scales::breaks_width(1),
#         limits = c(0, 4)
#       ) + 
#       theme_minimal() +
#       labs(
#         x = "Follicle",
#         y = "Number of clones",
#         fill = NULL,
#         title = glue("{HH} ({p}): {version_txt} clonal dominance per follicle"),
#         caption = "Numbers above bars indicate total GC B cell clone count per follicle."
#       ) +
#       theme(
#         plot.title   = element_text(size = 20, face = "bold"),
#         axis.title   = element_text(size = 14),
#         axis.text    = element_text(size = 12),
#         axis.text.y  = element_text(size = 14),
#         strip.text   = element_text(size = 12, face = "bold"),
#         legend.text  = element_text(size = 13),
#         legend.title = element_text(size = 14)
#       )
#     
#     ggsave(glue("{outdir}/clonal_D20_plot_{version}.png"), width = 10, height = 6)
#     
#   # })
#     
#     # Largest clone fraction, top 2, top 3, top 5.
#   
# }


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

# versions <- c("GC_B_cells", "all_cells")
# 
# lapply(versions, function(version) {
#   
#   # version <- "GC_B_cells"
#   
#   df_heavy_DXX <- if (version == "GC_B_cells") {
#     df_heavy %>% filter(L1_annotation == "GC_B_cells")
#   } else {
#     df_heavy
#   }
#   
#   gini_per_follicle <- df_heavy_DXX %>%
#     filter(!is.na(clone_subgroup_id_90_similarity)) %>%
#     dplyr::count(manual_ADT_full_ID, clone_subgroup_id_90_similarity, name = "n_cells") %>%
#     group_by(manual_ADT_full_ID) %>%
#     summarise(
#       gini = gini_coeff(n_cells),
#       total_clones = n_distinct(clone_subgroup_id_90_similarity),
#       .groups = "drop"
#     )
#   
#   version_txt <- str_replace_all(version, "_", " ")
#   
#   ggplot(gini_per_follicle, aes(x = factor(
#     str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer(),
#     levels = 1:max(str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer())
#   ), y = gini)) +
#     geom_point(size = 3, color = "steelblue") +
#     scale_y_continuous(limits = c(0, 1)) +
#     theme_minimal() +
#     labs(
#       x = "Follicle",
#       y = "Gini coefficient",
#       title = glue("{HH}: {version_txt} Gini coefficient per follicle")
#     ) 
#   
#   ggsave(glue("{outdir}/gini_coef_{version}.png"), width = 12, height = 6)
#   
# })  

# ------------------------------------------------------------------------------
# All patients combined 
# ------------------------------------------------------------------------------

patients <- c("HH117", "HH119", "HH151", "HH153")

# Load both patients
df_all <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
    filter(
      locus == "IGH",
      !is.na(manual_ADT_full_ID)
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows() %>% 
  mutate(
    clone_uid      = paste(patient_id, clone_subgroup_id_90_similarity, sep = "_")
  ) %>% 
  left_join(patient_to_condition, by = "patient")

# extra <- ""
# Remove largest clone as it "takes all the signal"
largest_clone <- df_all %>% filter(patient == "HH119") %>% count(clone_subgroup_id_90_similarity, sort = TRUE) %>% head(1) %>% pull(clone_subgroup_id_90_similarity)
df_all %>% filter(clone_subgroup_id_90_similarity == largest_clone) %>% pull(patient) %>% unique()
df_all <- df_all %>% filter(clone_subgroup_id_90_similarity != largest_clone)
extra <- "_largest_removed"
extra_subtitle <- "Largest clone removed"

# versions <- c("GC_B_cells", "all_cells")
outdir_combined <- glue("45_immcantation/plot/12_diversity_analysis/combined{extra}_update")
dir.create(outdir_combined, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# GC B cells count
# ------------------------------------------------------------------------------

df_all %>% 
  filter(L1_annotation == "GC_B_cells") %>% 
  dplyr::count(patient, condition, manual_ADT_full_ID) %>% 
  mutate(manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer() %>% as.factor()) %>%
  ggplot(aes(x = manual_ADT_full_ID_plot, y = n)) + 
  geom_col() + 
  geom_text(
    aes(x = manual_ADT_full_ID_plot, y = n, label = n),
    inherit.aes = FALSE,
    vjust = -0.5, size = 3
  ) + 
  facet_wrap(vars(patient), scales = "free_x", ncol = 2) +
  theme_bw() + 
  labs(
    title = glue("GC B cell (with BCR) count"), 
    subtitle = extra_subtitle, 
    x = "Follicle"
  )

ggsave(glue("{outdir_combined}/GC_B_cell_count.png"), width = 12, height = 6)

# ------------------------------------------------------------------------------
# D20 plot - GC B cell clones
# ------------------------------------------------------------------------------

compute_Dxx <- function(clone_counts, xx = 0.5) {
  sorted <- sort(clone_counts, decreasing = TRUE)
  cumulative <- cumsum(sorted) / sum(sorted)
  which(cumulative >= xx)[1]
}

# versions <- c("GC_B_cells", "all_cells")

# lapply(versions, function(version) {
version <- "GC_B_cells"

df_all_DXX <- df_all %>% 
  filter(L1_annotation == "GC_B_cells")

# Remove follicles that have less than 5 GC B cells
fol_to_rm <- df_all_DXX %>% 
  count(sample_clean_fol) %>% 
  filter(n < 5) %>% 
  pull(sample_clean_fol)

df_all_DXX <- df_all_DXX %>% 
  filter(!(sample_clean_fol %in% fol_to_rm))

d50_per_follicle <- df_all_DXX %>%
  filter(
    !is.na(clone_subgroup_id_90_similarity)
  ) %>%
  dplyr::count(patient, condition, sample_clean_fol, clone_subgroup_id_90_similarity, name = "n_cells") %>%
  group_by(patient, sample_clean_fol) %>%
  summarise(
    D50 = compute_Dxx(n_cells, 0.50),
    D20 = compute_Dxx(n_cells, 0.20),
    total_clones = n_distinct(clone_subgroup_id_90_similarity),
    .groups = "drop"
  ) %>%
  arrange(D50) %>% 
  mutate(
    sample_clean_fol_plot = sample_clean_fol %>% str_split_i("_", 2) %>% str_remove("Fol-") %>% as.integer() %>% as.factor()
    # manual_ADT_full_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer() %>% as.factor()
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

# D20 + D50
ggplot(d50_per_follicle_long, aes(x = sample_clean_fol_plot, y = value, fill = metric)) +
  geom_col() +
  geom_text(
    data = d50_per_follicle,
    aes(x = sample_clean_fol_plot, y = D50, label = total_clones),
    inherit.aes = FALSE,
    vjust = -0.5, size = 3
  ) +
  scale_fill_manual(
    values = c("D20_segment" = L1_colors[["GC_B_cells"]], "D50_extra" = "forestgreen"),
    labels = c("D20_segment" = "D20", "D50_extra" = "D20 to D50")
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(1),
    minor_breaks = scales::breaks_width(1)
  ) + 
  theme_bw() +
  facet_wrap(vars(patient), scales = "free_x", ncol = 2) +
  labs(
    x = "Follicle",
    y = "Number of clones",
    fill = NULL,
    title = glue("{version_txt} clonal dominance per follicle"),
    subtitle = extra_subtitle, 
    caption = "Numbers above bars indicate total GC B cell clone count per follicle."
  ) 

ggsave(glue("{outdir_combined}/clonal_D20_D50_plot_{version}.png"), width = 17, height = 10)

# D20
d50_per_follicle_long %>% 
  filter(metric == "D20_segment") %>% 
  ggplot(aes(x = sample_clean_fol_plot, y = value, fill = metric)) +
  geom_col() +
  geom_text(
    data = d50_per_follicle,
    aes(x = sample_clean_fol_plot, y = D20, label = total_clones),
    inherit.aes = FALSE,
    vjust = -0.5, size = 4
  ) +
  scale_fill_manual(
    values = c("D20_segment" = L1_colors[["GC_B_cells"]]), 
    labels = c("D20_segment" = "D20")
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(1), 
    minor_breaks = scales::breaks_width(1),
    limits = c(0, 4)
  ) + 
  theme_bw() +
  facet_wrap(vars(patient), scales = "free_x", ncol = 2) +
  labs(
    x = "Follicle",
    y = "Number of clones",
    fill = NULL,
    title = glue("{version_txt} clonal dominance per follicle"),
    caption = "Numbers above bars indicate total GC B cell clone count per follicle."
  ) +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    axis.text.y  = element_text(size = 14),
    strip.text   = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14)
  )

 ggsave(glue("{outdir_combined}/clonal_D20_plot_{version}.png"), width = 17, height = 10)

# })


# ------------------------------------------------------------------------------
# Largest clone fraction, top 2, top 3, top 5.
# ------------------------------------------------------------------------------

compute_topN_fraction <- function(clone_counts, n) {
 sorted <- sort(clone_counts, decreasing = TRUE)
 n <- min(n, length(sorted))   # cap at however many clones this follicle actually has
 sum(sorted[1:n]) / sum(sorted)
}

topN_per_follicle <- df_all_DXX %>%
 filter(
   !is.na(clone_subgroup_id_90_similarity)
 ) %>%
 dplyr::count(patient, condition, sample_clean_fol, clone_subgroup_id_90_similarity, name = "n_cells") %>%
 group_by(patient, sample_clean_fol) %>%
 summarise(
   Top1 = compute_topN_fraction(n_cells, 1),
   Top2 = compute_topN_fraction(n_cells, 2),
   Top5 = compute_topN_fraction(n_cells, 5),
   total_clones = n_distinct(clone_subgroup_id_90_similarity),
   total_cells  = sum(n_cells),
   .groups = "drop"
 ) %>%
 mutate(
   sample_clean_fol_plot = sample_clean_fol %>% str_split_i("_", 2) %>% str_remove("Fol-") %>% as.integer() %>% as.factor()
 )

topN_per_follicle_long <- topN_per_follicle %>%
 pivot_longer(
   cols = c(Top1, Top2, Top5),
   names_to = "topN",
   values_to = "fraction"
 ) %>%
 mutate(
   topN = factor(topN, levels = c("Top1", "Top2", "Top5"))
 )

topN_colors <- c(
 "Top1" = "firebrick",
 "Top2" = "darkorange",
 "Top5" = "goldenrod2"
)

# Per-follicle view, faceted by patient -- same layout as the D20/D50 plots above
for (topN_i in names(topN_colors)){
  
  # topN_i <- "Top1"
  N <- topN_i %>% str_remove("Top")
  topN_per_follicle_long %>% 
    filter(topN == topN_i) %>% 
    ggplot(aes(x = sample_clean_fol_plot, y = fraction, fill = topN, group = topN)) +
    geom_col() +
    scale_fill_manual(values = topN_colors, name = glue("Top {N} clones")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_bw() +
    scale_y_continuous(
      labels = scales::percent,
      breaks = scales::breaks_width(0.2),
      minor_breaks = scales::breaks_width(0.1),
      limits = c(0, 1)
    ) +
    facet_wrap(vars(patient), scales = "free_x", ncol = 2) +
    labs(
      x = "Follicle",
      y = "Fraction of GC B cells",
      title = glue("{version_txt} clonal dominance: fraction of cells in top {N} clones per follicle"),
      subtitle = extra_subtitle
    )
  
  ggsave(glue("{outdir_combined}/clonal_top{N}_fraction_per_follicle_{version}.png"), width = 17, height = 10)
  
}

# ------------------------------------------------------------------------------
# Gini - both patients
# ------------------------------------------------------------------------------

min_n <- 30

# lapply(versions, function(version) {
  
version <- "GC_B_cells"

df_gini <- if (version == "GC_B_cells") {
  df_all %>% filter(L1_annotation == "GC_B_cells")
} else {
  df_all
}

depth <- df_gini %>% 
  count(patient, sample_clean_fol) %>%
  filter(n >= min_n) %>%
  pull(n) %>%
  min()

set.seed(123)
gini_rare <- df_gini %>%
  group_by(patient, sample_clean_fol) %>%
  filter(n() >= min_n) %>%
  summarise(
    boot = list(replicate(100, {
      clones <- sample(clone_subgroup_id_90_similarity, depth, replace = FALSE)
      gini_coeff(as.vector(table(clones)))
    })),
    .groups = "drop"
  ) %>%
  mutate(
    gini       = map_dbl(boot, mean),
    gini_lower = map_dbl(boot, ~ quantile(.x, 0.025)),
    gini_upper = map_dbl(boot, ~ quantile(.x, 0.975))
  ) %>%
  select(-boot) %>%
  left_join(patient_to_condition, by = "patient")

# gini_combined <- df_gini %>%
#   filter(!is.na(clone_subgroup_id_90_similarity)) %>%
#   dplyr::count(patient, manual_ADT_full_ID, clone_subgroup_id_90_similarity, name = "n_cells") %>%
#   group_by(patient, manual_ADT_full_ID) %>%
#   summarise(
#     gini = gini_coeff(n_cells),
#     total_clones = n_distinct(clone_subgroup_id_90_similarity),
#     .groups = "drop"
#   ) %>% 
#   left_join(patient_to_condition, by = "patient")

version_txt <- str_replace_all(version, "_", " ")

# shared jitter position so the error bars move with their point, not independently
pos_jitter <- position_jitter(width = 0.4, height = 0, seed = 123)

ggplot(gini_rare, aes(x = patient, y = gini)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
  geom_errorbar(
    aes(ymin = gini_lower, ymax = gini_upper, color = patient),
    width = 0.05, position = pos_jitter, alpha = 0.6
  ) +
  geom_point(aes(color = patient), position = pos_jitter, size = 2.5, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = patient_color_values) +
  theme_bw() +
  labs(
    x = "Patient",
    y = "Gini coefficient",
    title = glue("Gini coefficient per follicle - {version_txt}"),
    subtitle = glue("Each point represents one follicle. All samples subsampled to {depth} cells without replacement."),
    caption = "Bootstrapped estimates (100 resamples, 95% CI)"
  ) +
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
  theme(legend.position = "none")

ggsave(glue("{outdir_combined}/gini_coef_{version}_combined.png"), width = 10, height = 8)

# })


# ------------------------------------------------------------------------------
# Gini diagnositic 
# ------------------------------------------------------------------------------

gini_coeff <- function(x) {
  x <- sort(x); n <- length(x)
  if (n < 2) return(NA_real_)
  sum((2 * seq_along(x) - n - 1) * x) / ((n - 1) * sum(x))
}

min_n   <- 30
depths  <- c(30, 50, 100, 200, 500, 1000)
n_rep   <- 50


version_txt <- str_replace_all(version, "_", " ")

# ---- Data: one row per cell with a clone, unique follicle ID ----
df_v <- df_all %>%
  { if (version == "GC_B_cells") filter(., L1_annotation == "GC_B_cells") else . } %>%
  filter(!is.na(clone_subgroup_id_90_similarity)) %>%
  distinct(cell_id, .keep_all = TRUE) %>%
  mutate(fol_id = paste(patient, sample_clean_fol, sep = "_"))

# ---- Full-data Gini per follicle ----
gini_full <- df_v %>%
  count(patient, fol_id, clone_subgroup_id_90_similarity, name = "n_cells") %>%
  group_by(patient, fol_id) %>%
  summarise(gini = gini_coeff(n_cells),
            total_cells = sum(n_cells),
            n_clones = n(), .groups = "drop") %>%
  filter(total_cells >= min_n) %>%
  left_join(patient_to_condition, by = "patient")

# ============================================================================
# Diagnostic 1: is full-data Gini driven by follicle size?
# ============================================================================
ct  <- cor.test(log10(gini_full$total_cells), gini_full$gini,
                method = "spearman", exact = FALSE)
fit <- lm(gini ~ log10(total_cells) + patient, data = gini_full)   # size effect within patients

ggplot(gini_full, aes(total_cells, gini, color = patient)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() +
  scale_color_manual(values = patient_color_values) +
  theme_bw() +
  labs(x = "Cells per follicle (log scale)", y = "Gini (full data)",
       title = glue("Diagnostic 1: Gini vs follicle size - {version_txt}"),
       subtitle = glue("Spearman rho = {round(ct$estimate, 2)}, p = {signif(ct$p.value, 2)}; n = {nrow(gini_full)} follicles"))

ggsave(glue("{outdir_combined}/diag1_size_{version}.png"), p1, width = 7, height = 5)

# ============================================================================
# Diagnostic 2: is the follicle ranking stable across depths?
# ============================================================================
set.seed(123)
gini_curve <- df_v %>%
  filter(fol_id %in% gini_full$fol_id) %>%
  group_by(patient, fol_id) %>%
  group_modify(~ {
    cl <- .x$clone_subgroup_id_90_similarity
    map_dfr(depths[depths <= length(cl)], function(d) {
      tibble(depth = d,
             gini  = mean(replicate(n_rep, gini_coeff(as.vector(table(sample(cl, d)))))))
    })
  }) %>%
  ungroup() %>%
  left_join(patient_to_condition, by = "patient")

rank_agree <- gini_curve %>%
  left_join(gini_full %>% select(fol_id, gini_full = gini), by = "fol_id") %>%
  group_by(depth) %>%
  summarise(n_follicles = n(),
            spearman_vs_full = cor(gini, gini_full, method = "spearman"),
            .groups = "drop")

ggplot(gini_curve, aes(depth, gini, group = fol_id, color = patient)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5) +
  scale_x_log10(
    breaks = depths
  ) +
  scale_color_manual(values = patient_color_values) +
  facet_grid(cols = vars(patient)) +
  theme_bw() +
  labs(
    x = "Subsampling depth (cells, log scale)", y = "Gini",
    title = glue("Gini across depths - {version_txt}"),
    # subtitle = "Follicles smaller than a given depth drop out at that depth"
    subtitle = "Gini increases with subsamping depth."
  )

ggsave(glue("{outdir_combined}/diag2_depth_{version}.png"), p2, width = 9, height = 5)

# ============================================================================
# Diagnostic 3: do the compared groups differ in follicle size?
# ============================================================================
compare_col <- "patient"   # change to the grouping you actually compare

size_test <- if (n_distinct(gini_full[[compare_col]]) == 2) {
  wilcox.test(as.formula(glue("total_cells ~ {compare_col}")), data = gini_full)
} else {
  kruskal.test(as.formula(glue("total_cells ~ {compare_col}")), data = gini_full)
}

p3 <- ggplot(gini_full, aes(.data[[compare_col]], total_cells)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
  geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = patient_color_values) +
  theme_bw() +
  labs(x = compare_col, y = "Cells per follicle (log scale)",
       title = glue("Diagnostic 3: follicle size by group - {version_txt}"),
       subtitle = glue("p = {signif(size_test$p.value, 2)}"))

ggsave(glue("{outdir_combined}/diag3_size_by_group_{version}.png"), p3, width = 6, height = 5)

# ---- Printed summary ----
message(glue("\n===== {version} ====="))
print(ct)
print(summary(fit)$coefficients["log10(total_cells)", , drop = FALSE])
print(rank_agree)
print(size_test)

invisible(list(full = gini_full, curve = gini_curve, rank_agree = rank_agree))


















# min_n <- 30   # consider raising, see point 3
# 
# # lapply(versions, function(version) {
# 
# version <- "GC_B_cells" 
# 
# df_v <- df_all %>%
#   { if (version == "GC_B_cells") filter(., L1_annotation == "GC_B_cells") else . } %>%
#   filter(!is.na(clone_subgroup_id_90_similarity), !is.na(clone_uid)) %>%
#   mutate(fol_id = paste(patient, sample_clean_fol, sep = "_"))
# 
# depth <- df_v %>%
#   count(fol_id) %>%
#   filter(n >= min_n) %>%
#   pull(n) %>%
#   min()
# 
# version_txt <- str_replace_all(version, "_", " ")
# message(glue("{version}: depth = {depth}"))
# 
# # ---- Gini (without replacement) ----
# set.seed(123)
# gini_rare <- df_v %>%
#   group_by(patient, sample_clean_fol, fol_id) %>%
#   filter(n() >= min_n) %>%
#   summarise(
#     boot = list(replicate(100, {
#       clones <- sample(clone_subgroup_id_90_similarity, depth, replace = FALSE)
#       gini_coeff(as.vector(table(clones)))
#     })),
#     .groups = "drop"
#   ) %>%
#   mutate(gini = map_dbl(boot, mean)) %>%
#   select(-boot) %>%
#   left_join(patient_to_condition, by = "patient")
# 
# ggplot(gini_rare, aes(x = patient, y = gini)) +
#   geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
#   geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_color_manual(values = patient_color_values) +
#   facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(
#     x = "Patient", y = "Gini coefficient",
#     title = glue("Gini coefficient per follicle - {version_txt}"),
#     subtitle = glue("Each point is one follicle (mean of 100 subsamples of {depth} cells, without replacement)."),
#     caption = glue("Follicles with < {min_n} cells excluded.")
#   )
# 
# # ggsave(glue("{outdir_combined}/gini_coef_{version}_combined.png"), p_gini, width = 8, height = 8)
# 
# ---- Shannon (alphaDiversity, with replacement) ----
# set.seed(123)
# curve_all <- alphaDiversity(
#   df_v,
#   group = "fol_id",
#   clone = "clone_uid",
#   min_q = 0, max_q = 4, step_q = 0.1,
#   ci = 0.95, nboot = 100,
#   min_n = min_n, max_n = depth
# )
# 
# shannon_df <- curve_all@diversity %>%
#   filter(q == 1) %>%
#   left_join(df_v %>% distinct(fol_id, patient, sample_clean_fol), by = "fol_id") %>%
#   left_join(patient_to_condition, by = "patient")
# 
# ggplot(shannon_df, aes(x = patient, y = d)) +
#   geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
#   geom_jitter(aes(color = patient), width = 0.1, size = 2.5, alpha = 0.8) +
#   scale_color_manual(values = patient_color_values) +
#   facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(
#     x = "Patient", y = "Shannon diversity (Hill number, q = 1)",
#     title = glue("Shannon diversity per follicle - {version_txt}"),
#     subtitle = glue("Each point is one follicle (mean of 100 resamples of {depth} cells, with replacement)."),
#     caption = glue("Follicles with < {min_n} cells excluded.")
#   )
#   
#   # ggsave(glue("{outdir_combined}/shannon_diversity_{version}_min_n_{min_n}_combined.png"), p_shannon, width = 8, height = 8)
#   
#   # invisible(NULL)
# # })



# ------------------------------------------------------------------------------
# Shannon - both patients
# ------------------------------------------------------------------------------

# Run alphaDiversity per patient and combine
# lapply(versions, function(version) {
  
  # for (min_n in c(30)) {
    
min_n <- 30

# Common rarefaction depth: the smallest sample that passes min_n
depth <- df_all %>%
  { if (version == "GC_B_cells") filter(., L1_annotation == "GC_B_cells") else . } %>% 
  count(sample_clean_fol) %>%
  filter(n >= min_n) %>%
  pull(n) %>%
  min()

set.seed(123)
curve_all <- alphaDiversity(
  df_all %>% { if (version == "GC_B_cells") filter(., L1_annotation == "GC_B_cells") else . },
  group  = "sample_clean_fol",
  clone  = "clone_uid",
  min_q  = 0, max_q = 4, step_q = 0.1,
  ci     = 0.95, nboot = 100,
  min_n  = min_n, max_n = depth
)

# Map the sample key back to patient, compartment and condition
sample_meta <- df_all %>%
  distinct(patient, sample_clean_fol) %>%
  left_join(patient_to_condition, by = "patient")

shannon_df <- curve_all@diversity %>%
  filter(q == 1) %>%
  left_join(sample_meta, by = "sample_clean_fol")

version_txt <- str_replace_all(version, "_", " ")

# shared jitter position so the error bars move with their point, not independently
pos_jitter <- position_jitter(width = 0.4, height = 0, seed = 123)

shannon_df %>% 
  ggplot(aes(x = patient, y = d)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = "grey90") +
  geom_errorbar(
    aes(ymin = d_lower, ymax = d_upper, color = patient),
    width = 0.05, position = pos_jitter, alpha = 0.6
  ) +
  geom_point(aes(color = patient), position = pos_jitter, size = 2.5, alpha = 0.8) +
  scale_color_manual(values = patient_color_values) +
  theme_bw() +
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
  labs(
    x = "Patient",
    y = "Shannon diversity (q=1)",
    title = glue("Shannon diversity per follicle - {version_txt}"),
    subtitle = glue("All samples rarefied to {depth} cells with replacement."),
    caption = "Bootstrapped estimates (100 resamples, 95% CI)"
  ) +
  theme(legend.position = "none")


ggsave(glue("{outdir_combined}/shannon_diversity_{version}_min_n_{min_n}_combined.png"), width = 10, height = 8)
    
  # }
  
# })


