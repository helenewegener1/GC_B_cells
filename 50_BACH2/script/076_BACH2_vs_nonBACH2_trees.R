library(tidyverse)
library(glue)
library(ape)

# ==============================================================================
# Compare GCtrees of clones with BACH2 cells (50_BACH2) to GCtrees of clones
# without BACH2 cells (48_GCtree)
#
# Non-BACH2 trees are checked so that they contain no BACH2 cells:
#   - clones that are also in the BACH2 set (same patient + clone id) are removed
#   - clones with >= 1 BACH2 cell (any IGH cell of the clone, not only the cells
#     in the tree) are removed
# All exclusions are written to tree_exclusions.tsv with the reason.
#
# Tree-level metrics are compared between the groups with
#   - Wilcoxon rank-sum test (unadjusted)
#   - linear model: metric ~ group + log10(n_cells) + patient, because BACH2
#     clones may be larger, and many metrics depend on tree size
# ==============================================================================

set.seed(1)

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

resolve_LC_list_germlined <- readRDS("50_BACH2/out/resolve_LC_list_germlined_L1_L3.rds")

tree_sets <- list(
  BACH2 = list(
    tree_dir = "50_BACH2/plot",
    idmap_dir = "50_BACH2/out",
    seq_dir = readRDS("50_BACH2/out/rds/seq_dir.rds")
  ),
  nonBACH2 = list(
    tree_dir = "48_GCtree/plot",
    idmap_dir = "48_GCtree/out",
    seq_dir = readRDS("48_GCtree/out/rds/seq_dir.rds")
  )
)

outdir <- "50_BACH2/plot/076_BACH2_vs_nonBACH2_trees"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

celltype_col <- "L1_L3_annotation"
BACH2_labels <- c("BACH2+", "BACH2_cycling")

group_labels <- c(BACH2 = "BACH2 trees", nonBACH2 = "Non-BACH2 trees")
group_colors <- c("BACH2 trees" = "#D55E00", "Non-BACH2 trees" = "grey40")

min_cells_tree <- 3         # trees with fewer cells in the tree are not compared

# Clonal burst, as in 075_BACH2_tree_position.R
burst_min_children <- 3
burst_min_parent_cells <- 10

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Number of descendant nodes for every node in the tree
count_descendants <- function(tree){
  n_desc <- integer(Ntip(tree) + Nnode(tree))
  edges <- reorder(tree, "postorder")$edge
  for (i in seq_len(nrow(edges))){
    n_desc[edges[i, 1]] <- n_desc[edges[i, 1]] + 1L + n_desc[edges[i, 2]]
  }
  n_desc
}

# Shannon diversity of a vector of counts
shannon <- function(x){
  p <- x[x > 0] / sum(x)
  -sum(p * log(p))
}

# ------------------------------------------------------------------------------
# Per tree: cells and tree-level metrics
# ------------------------------------------------------------------------------

analyse_tree <- function(clone_full_name, set_name){

  # clone_full_name <- "HH117_clone_nr_1_clone_1881_1"; set_name <- "BACH2"

  set <- tree_sets[[set_name]]

  parts <- str_match(clone_full_name, "^(HH[0-9]+)_clone_nr_([0-9]+)_clone_(.+)$")
  HH <- parts[, 2]
  clone <- parts[, 4]

  tree_ids <- tibble(
    set = set_name,
    group = group_labels[[set_name]],
    patient = HH,
    clone = clone,
    clone_nr = as.integer(parts[, 3]),
    clone_full_name = clone_full_name
  )

  # All IGH cells of the clone (used for the BACH2 check)
  df_clone <- resolve_LC_list_germlined[[HH]] %>%
    filter(
      locus == "IGH",
      clone_subgroup_id_90_similarity == clone
    ) %>%
    mutate(cell_type = .data[[celltype_col]])

  tree_ids <- tree_ids %>%
    mutate(
      n_cells_clone = nrow(df_clone),
      n_BACH2_clone = sum(df_clone$cell_type %in% BACH2_labels)
    )

  tree_file <- glue("{set$tree_dir}/{clone_full_name}/{clone_full_name}.inference.1.nk")
  idmap_file <- glue("{set$idmap_dir}/{clone_full_name}/idmap.txt")

  if (!file.exists(tree_file) || !file.exists(idmap_file)){
    return(list(tree_metrics = tree_ids %>% mutate(has_tree = FALSE)))
  }

  # ---- Tree ----

  tree <- read.tree(tree_file)
  node_labels <- c(tree$tip.label, tree$node.label)
  node_labels[Ntip(tree) + 1] <- "GL"

  depth <- node.depth.edgelength(tree)
  n_edges_from_GL <- node.depth.edgelength(compute.brlen(tree, 1))
  n_children <- tabulate(tree$edge[, 1], nbins = length(node_labels))
  n_desc <- count_descendants(tree)

  node_df <- tibble(
    node = node_labels,
    mut_from_GL = depth,
    n_edges_from_GL = n_edges_from_GL,
    is_tip = seq_along(node_labels) <= Ntip(tree),
    n_children = n_children,
    n_desc_nodes = n_desc
  )

  # ---- Cells on the tree (cells identical to GL sit on the GL node) ----

  idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id_gctree"), show_col_types = FALSE) %>%
    separate_longer_delim(sequence_id_gctree, delim = ":") %>%
    filter(sequence_id_gctree != "GL")

  sequence_ids <- set$seq_dir[[clone_full_name]]

  cells <- df_clone %>%
    mutate(sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)]) %>%
    select(sequence_id, sequence_id_gctree, cell_type) %>%
    inner_join(idmap, by = "sequence_id_gctree", relationship = "many-to-one") %>%
    inner_join(node_df %>% filter(str_detect(node, "^seq[0-9]+$") | node == "GL"), by = "node") %>%
    add_count(node, name = "node_n_cells")

  cells <- bind_cols(tree_ids[rep(1, nrow(cells)), ] %>% select(set, group, patient, clone, clone_full_name), cells)

  # ---- Tree-level metrics ----

  observed_nodes <- unique(cells$node)
  node_abundance <- cells %>% filter(node != "GL") %>% dplyr::count(node)

  # Clonal bursts: observed node (not GL) with >= burst_min_parent_cells cells and >= burst_min_children children
  burst_idx <- which(
    n_children >= burst_min_children &
    node_labels %in% node_abundance$node[node_abundance$n >= burst_min_parent_cells]
  )
  burst_nodes <- node_labels[burst_idx]
  burst_offspring_nodes <- node_labels[tree$edge[tree$edge[, 1] %in% burst_idx, 2]]

  internal_idx <- which(n_children > 0)

  tree_metrics <- tree_ids %>%
    mutate(
      has_tree = TRUE,

      # Size
      n_cells = nrow(cells),
      n_BACH2_tree = sum(cells$cell_type %in% BACH2_labels),
      n_nodes = length(node_labels),
      n_genotypes = sum(observed_nodes != "GL"),
      n_unobserved_nodes = sum(!node_labels %in% observed_nodes & node_labels != "GL"),

      # Mutation load and depth
      max_mut_from_GL = max(depth),
      median_cell_mut_from_GL = median(cells$mut_from_GL),
      mean_tip_edges_from_GL = mean(n_edges_from_GL[seq_len(Ntip(tree))]),
      frac_cells_on_GL = mean(cells$node == "GL"),

      # Topology
      frac_cells_internal = mean(!cells$is_tip),
      mean_children_internal = mean(n_children[internal_idx]),
      max_children = max(n_children),

      # Expansion and diversity
      genotypes_per_cell = n_genotypes / n_cells,
      largest_genotype_frac = if (nrow(node_abundance) > 0) max(node_abundance$n) / n_cells else NA_real_,
      genotype_shannon = if (nrow(node_abundance) > 0) shannon(node_abundance$n) else NA_real_,

      # Clonal bursts
      n_bursts = length(burst_idx),
      has_burst = n_bursts > 0,
      frac_cells_in_bursts = mean(cells$node %in% c(burst_nodes, burst_offspring_nodes))
    )

  list(tree_metrics = tree_metrics, cells = cells)
}

# ------------------------------------------------------------------------------
# Run all trees of both sets
# ------------------------------------------------------------------------------

results <- imap(tree_sets, function(set, set_name){

  samples <- list.dirs(set$tree_dir, full.names = FALSE, recursive = FALSE) %>%
    str_subset("^HH[0-9]+_clone_nr_[0-9]+_clone_")

  message(glue("{set_name}: {length(samples)} trees"))

  map(samples, ~ analyse_tree(.x, set_name))
}) %>%
  list_flatten()

tree_metrics_all <- map_dfr(results, "tree_metrics")
cells_all <- map_dfr(results, "cells")

# ------------------------------------------------------------------------------
# Exclusions: make sure non-BACH2 trees contain no BACH2 cells
# ------------------------------------------------------------------------------

BACH2_clones <- tree_metrics_all %>%
  filter(set == "BACH2") %>%
  distinct(patient, clone)

tree_metrics_all <- tree_metrics_all %>%
  mutate(
    in_BACH2_set = paste(patient, clone) %in% paste(BACH2_clones$patient, BACH2_clones$clone),
    exclusion_reason = case_when(
      !has_tree ~ "No GCtree output (tree or idmap missing)",
      set == "nonBACH2" & in_BACH2_set ~ "Same clone as a BACH2 tree",
      set == "nonBACH2" & n_BACH2_clone > 0 ~ "Clone contains BACH2 cells",
      set == "BACH2" & n_BACH2_tree == 0 ~ "No BACH2 cells in the tree",
      n_cells < min_cells_tree ~ glue("< {min_cells_tree} cells in the tree"),
      TRUE ~ NA_character_
    )
  )

exclusions <- tree_metrics_all %>%
  filter(!is.na(exclusion_reason)) %>%
  select(set, patient, clone, clone_nr, clone_full_name, n_cells_clone, n_BACH2_clone, n_BACH2_tree, n_cells, exclusion_reason)

write_delim(exclusions, glue("{outdir}/tree_exclusions.tsv"), delim = "\t")

# Overview of included / excluded trees
tree_metrics_all %>%
  mutate(status = replace_na(exclusion_reason, "Included")) %>%
  dplyr::count(set, patient, status) %>%
  pivot_wider(names_from = patient, values_from = n, values_fill = 0) %>%
  print(n = Inf)

tree_metrics <- tree_metrics_all %>%
  filter(is.na(exclusion_reason)) %>%
  mutate(group = factor(group, levels = group_labels))

cells <- cells_all %>%
  semi_join(tree_metrics, by = c("set", "clone_full_name")) %>%
  mutate(group = factor(group, levels = group_labels))

# Final check: no BACH2 cells in the non-BACH2 trees
stopifnot(sum(tree_metrics$n_BACH2_clone[tree_metrics$set == "nonBACH2"]) == 0)
stopifnot(!any(cells$cell_type[cells$set == "nonBACH2"] %in% BACH2_labels))

write_delim(tree_metrics, glue("{outdir}/tree_metrics.tsv"), delim = "\t")

# ==============================================================================
# 1) Tree-level metrics
# ==============================================================================

metric_labels <- c(
  n_cells = "N cells in tree",
  n_genotypes = "N observed genotypes (nodes)",
  max_mut_from_GL = "Max mutations from GL",
  median_cell_mut_from_GL = "Median mutations from GL per cell",
  mean_tip_edges_from_GL = "Mean N edges GL to tip",
  frac_cells_on_GL = "Fraction of cells identical to GL",
  frac_cells_internal = "Fraction of cells on internal nodes",
  mean_children_internal = "Mean N children per internal node",
  max_children = "Max N children of a node",
  genotypes_per_cell = "Observed genotypes per cell",
  largest_genotype_frac = "Fraction of cells in largest genotype",
  genotype_shannon = "Shannon diversity of genotypes",
  n_bursts = glue("N clonal bursts (>= {burst_min_parent_cells} cells, >= {burst_min_children} children)"),
  frac_cells_in_bursts = "Fraction of cells in clonal bursts"
)

df_metrics <- tree_metrics %>%
  mutate(n_cells_size = n_cells) %>%
  select(group, patient, clone_full_name, n_cells_size, all_of(names(metric_labels))) %>%
  pivot_longer(cols = all_of(names(metric_labels)), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = names(metric_labels))) %>%
  filter(!is.na(value))

# ---- Tests ----

# Wilcoxon rank-sum: BACH2 vs non-BACH2 trees (all patients, and per patient)
run_wilcox <- function(df, ...){
  df %>%
    group_by(..., metric) %>%
    summarise(
      n_BACH2 = sum(group == group_labels[["BACH2"]]),
      n_nonBACH2 = sum(group == group_labels[["nonBACH2"]]),
      median_BACH2 = median(value[group == group_labels[["BACH2"]]]),
      median_nonBACH2 = median(value[group == group_labels[["nonBACH2"]]]),
      p_value = if (n_BACH2 >= 2 && n_nonBACH2 >= 2) wilcox.test(value[group == group_labels[["BACH2"]]], value[group == group_labels[["nonBACH2"]]], exact = FALSE)$p.value else NA_real_,
      .groups = "drop"
    ) %>%
    group_by(...) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

wilcox_all <- run_wilcox(df_metrics)
wilcox_patient <- run_wilcox(df_metrics, patient)

# Linear model adjusting for tree size and patient: metric ~ group + log10(n_cells) + patient
lm_results <- df_metrics %>%
  filter(metric != "n_cells") %>%
  group_by(metric) %>%
  group_modify(function(df, key){
    fit <- lm(value ~ group + log10(n_cells_size) + patient, data = df)
    coefs <- summary(fit)$coefficients
    row <- grep("^group", rownames(coefs))
    tibble(
      estimate_BACH2_vs_nonBACH2 = -coefs[row, "Estimate"],   # factor reference is BACH2 trees
      std_error = coefs[row, "Std. Error"],
      p_value_size_adjusted = coefs[row, "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj_size_adjusted = p.adjust(p_value_size_adjusted, method = "BH"))

test_results <- wilcox_all %>%
  left_join(lm_results, by = "metric")

print(test_results, n = Inf, width = Inf)
write_delim(test_results, glue("{outdir}/tree_metrics_tests.tsv"), delim = "\t")
write_delim(wilcox_patient, glue("{outdir}/tree_metrics_wilcox_per_patient.tsv"), delim = "\t")

# ---- Plots ----

format_p <- function(p){
  case_when(is.na(p) ~ "= NA", p < 0.001 ~ "< 0.001", TRUE ~ paste("=", signif(p, 2)))
}

label_df <- df_metrics %>%
  group_by(metric) %>%
  summarise(y = max(value, na.rm = TRUE) * 1.12, .groups = "drop") %>%
  left_join(test_results, by = "metric") %>%
  mutate(p_label = glue("Wilcoxon p {format_p(p_adj)}\nSize-adj. p {format_p(p_adj_size_adjusted)}"))

# Boxplots per metric, all patients
df_metrics %>%
  ggplot(aes(x = group, y = value, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(shape = patient), width = 0.15, height = 0, alpha = 0.5) +
  geom_text(data = label_df, aes(x = 1.5, y = y, label = p_label), inherit.aes = FALSE, size = 2.6) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4, labeller = as_labeller(metric_labels)) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(
    title = "GCtrees of clones with vs. without BACH2 cells",
    subtitle = glue(
      "One point per tree. BACH2 trees: n = {sum(tree_metrics$set == 'BACH2')}, non-BACH2 trees: n = {sum(tree_metrics$set == 'nonBACH2')}.\n",
      "P-values (BH adjusted): Wilcoxon rank-sum test, and linear model adjusted for log10(N cells) and patient"
    ),
    x = NULL,
    y = NULL,
    shape = "Patient"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 8), legend.position = "bottom")

ggsave(glue("{outdir}/tree_metrics_boxplot.png"), width = 13, height = 13)

# Same, per patient
df_metrics %>%
  ggplot(aes(x = patient, y = value, color = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4, labeller = as_labeller(metric_labels)) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "GCtrees of clones with vs. without BACH2 cells, per patient",
    subtitle = "One point per tree. Wilcoxon tests per patient in tree_metrics_wilcox_per_patient.tsv",
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 8), legend.position = "bottom")

ggsave(glue("{outdir}/tree_metrics_boxplot_per_patient.png"), width = 14, height = 13)

# Metrics against tree size: shows whether a difference is explained by size
df_metrics %>%
  filter(metric != "n_cells") %>%
  ggplot(aes(x = n_cells_size, y = value, color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
  scale_x_log10() +
  facet_wrap(~ metric, scales = "free_y", ncol = 4, labeller = as_labeller(metric_labels)) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Tree metrics against tree size",
    subtitle = "One point per tree. Lines: linear fit per group on log10(N cells)",
    x = "N cells in tree (log10)",
    y = NULL,
    color = NULL
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 8), legend.position = "bottom")

ggsave(glue("{outdir}/tree_metrics_vs_size.png"), width = 13, height = 12)

# ==============================================================================
# 2) Subtype composition of the trees
# ==============================================================================
#
# Fractions are computed among the non-BACH2 cells of each tree, so that the
# composition of the remaining cells can be compared between the groups.

composition <- cells %>%
  filter(!cell_type %in% BACH2_labels, !is.na(cell_type)) %>%
  dplyr::count(group, patient, clone_full_name, cell_type) %>%
  complete(nesting(group, patient, clone_full_name), cell_type, fill = list(n = 0)) %>%
  group_by(clone_full_name) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

composition_order <- composition %>%
  group_by(cell_type) %>%
  summarise(n = sum(n)) %>%
  arrange(desc(n)) %>%
  pull(cell_type)

composition_tests <- composition %>%
  group_by(cell_type) %>%
  summarise(
    median_frac_BACH2 = median(frac[group == group_labels[["BACH2"]]]),
    median_frac_nonBACH2 = median(frac[group == group_labels[["nonBACH2"]]]),
    p_value = wilcox.test(frac[group == group_labels[["BACH2"]]], frac[group == group_labels[["nonBACH2"]]], exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    cell_type = factor(cell_type, levels = composition_order)
  ) %>%
  arrange(cell_type)

print(composition_tests)
write_delim(composition, glue("{outdir}/composition_per_tree.tsv"), delim = "\t")
write_delim(composition_tests, glue("{outdir}/composition_tests.tsv"), delim = "\t")

# Mean composition per group and patient
composition %>%
  group_by(group, patient, cell_type) %>%
  summarise(mean_frac = mean(frac), n_trees = n_distinct(clone_full_name), .groups = "drop") %>%
  mutate(
    cell_type = factor(cell_type, levels = composition_order),
    group_plot = glue("{group}\n(n = {n_trees})")
  ) %>%
  ggplot(aes(x = group_plot, y = mean_frac, fill = cell_type)) +
  geom_col() +
  facet_wrap(~ patient, scales = "free_x", nrow = 1) +
  labs(
    title = "Subtype composition of the trees (BACH2 cells excluded)",
    subtitle = "Mean over trees of the fraction of each subtype among the non-BACH2 cells of the tree",
    x = NULL,
    y = "Mean fraction of cells",
    fill = "Subtype"
  ) +
  theme_bw()

ggsave(glue("{outdir}/composition_barplot.png"), width = 11, height = 6)

# Fraction of each subtype per tree
comp_labels <- composition %>%
  group_by(cell_type) %>%
  summarise(y = max(frac) * 1.1, .groups = "drop") %>%
  left_join(composition_tests %>% mutate(cell_type = as.character(cell_type)), by = "cell_type") %>%
  mutate(
    p_label = if_else(p_adj < 0.001, "p < 0.001", glue("p = {signif(p_adj, 2)}")),
    cell_type = factor(cell_type, levels = composition_order)
  )

composition %>%
  mutate(cell_type = factor(cell_type, levels = composition_order)) %>%
  ggplot(aes(x = cell_type, y = frac, color = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(aes(shape = patient), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.5) +
  geom_text(data = comp_labels, aes(x = cell_type, y = y, label = p_label), inherit.aes = FALSE, size = 3) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Subtype fractions per tree (BACH2 cells excluded)",
    subtitle = "One point per tree. P-values: Wilcoxon rank-sum test, BH adjusted",
    x = NULL,
    y = "Fraction of non-BACH2 cells in the tree",
    color = NULL,
    shape = "Patient"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave(glue("{outdir}/composition_boxplot.png"), width = max(8, 1.1 * length(composition_order)), height = 6)
