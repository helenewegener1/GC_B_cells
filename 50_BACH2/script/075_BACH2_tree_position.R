library(tidyverse)
library(glue)
library(ape)

# ==============================================================================
# Where in the GCtrees are the BACH2+ cells?
#
# For each cell, the position of its GCtree node is described by:
#   mut_from_GL   : mutations from the germline root (sum of branch lengths)
#   rel_depth     : mut_from_GL / deepest node in the tree (0 = GL, 1 = deepest)
#   is_internal   : node is an observed ancestor (internal node) rather than a tip
#   node_n_cells  : number of cells sharing the node genotype (expansion)
#   n_desc_nodes  : number of nodes below the node (size of the subtree it gives rise to)
#
# BACH2+ cells are compared with the other GC B cells of the same tree:
#   1) per tree: label permutation test (labels shuffled, tree fixed) for
#        - median rel_depth (BACH2+ closer to GL?)
#        - mean pairwise patristic distance among BACH2+ cells (clustered?)
#   2) across trees: paired Wilcoxon test per patient on the per-tree summaries
# ==============================================================================

set.seed(1)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/06_resolve_LC_germlined.rds")

# Named list: clone_full_name -> named vector (GCtree sequence id -> sequence_alignment),
# same structure as 48_GCtree/out/rds/seq_dir.rds, but for the BACH2 clones
seq_dir <- readRDS("50_BACH2/out/rds/seq_dir.rds")

tree_dir <- "50_BACH2/plot"
idmap_dir <- "50_BACH2/out"

outdir <- "50_BACH2/plot_07_tree_position"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

BACH2_labels <- c("BACH2+", "BACH2_cycling")
n_perm <- 1000
min_cells_group <- 2   # min BACH2+ and other GC B cells in a tree to test it

# All clones with a GCtree output folder
samples <- list.dirs(tree_dir, full.names = FALSE, recursive = FALSE) %>%
  str_subset("^HH[0-9]+_clone_nr_[0-9]+_clone_")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Number of descendant nodes (tips + internal) for every node in the tree
count_descendants <- function(tree){

  n_nodes <- Ntip(tree) + Nnode(tree)
  n_desc <- integer(n_nodes)

  # Postorder: children are always visited before their parent
  edges <- reorder(tree, "postorder")$edge

  for (i in seq_len(nrow(edges))){
    parent <- edges[i, 1]
    child <- edges[i, 2]
    n_desc[parent] <- n_desc[parent] + 1L + n_desc[child]
  }

  n_desc
}

# Permutation p-value (lower tail: observed is smaller than expected)
perm_p_lower <- function(observed, null){
  null <- null[!is.na(null)]
  if (is.na(observed) || length(null) == 0) return(NA_real_)
  (sum(null <= observed) + 1) / (length(null) + 1)
}

# ------------------------------------------------------------------------------
# Per tree: cell positions + permutation tests
# ------------------------------------------------------------------------------

analyse_tree <- function(clone_full_name){

  # clone_full_name <- samples[1]

  parts <- str_match(clone_full_name, "^(HH[0-9]+)_clone_nr_([0-9]+)_clone_(.+)$")
  HH <- parts[, 2]
  clone_nr <- as.integer(parts[, 3])
  clone <- parts[, 4]

  tree_file <- glue("{tree_dir}/{clone_full_name}/{clone_full_name}.inference.1.nk")
  idmap_file <- glue("{idmap_dir}/{clone_full_name}/idmap.txt")

  if (!file.exists(tree_file) || !file.exists(idmap_file)){
    message(glue("Skipping {clone_full_name}: no tree or idmap"))
    return(NULL)
  }

  # Load gctree
  tree <- read.tree(tree_file)
  node_labels <- c(tree$tip.label, tree$node.label)

  # Position of all nodes. The root of the GCtree is the germline (GL)
  depth <- node.depth.edgelength(tree)
  max_depth <- max(depth)

  node_df <- tibble(
    node = node_labels,
    mut_from_GL = depth,
    rel_depth = if (max_depth > 0) depth / max_depth else NA_real_,
    is_internal = seq_along(node_labels) > Ntip(tree),
    n_desc_nodes = count_descendants(tree)
  ) %>%
    # Keep observed-sequence nodes only (drop GL + unobserved ancestors)
    filter(str_detect(node, "^seq[0-9]+$"))

  # Expand collapsed (identical-sequence) nodes back to individual GCtree sequence ids
  idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id_gctree"), show_col_types = FALSE) %>%
    filter(node != "GL") %>%
    separate_longer_delim(sequence_id_gctree, delim = ":")

  # Cells of the clone, matched to their GCtree sequence id
  sequence_ids <- seq_dir[[clone_full_name]]

  df_clone <- resolve_LC_list_germlined[[HH]] %>%
    filter(
      locus == "IGH",
      clone_subgroup_id_90_similarity == clone
    ) %>%
    mutate(
      sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)]
    )

  # One row per cell with the position of its node
  cells <- df_clone %>%
    select(sequence_id, sequence_id_gctree, L1_annotation, L1_L3_annotation) %>%
    inner_join(idmap, by = "sequence_id_gctree", relationship = "many-to-one") %>%
    inner_join(node_df, by = "node") %>%
    add_count(node, name = "node_n_cells") %>%
    mutate(
      group = case_when(
        L1_L3_annotation %in% BACH2_labels ~ "BACH2+",
        L1_annotation == "GC_B_cells" ~ "Other GC B cells",
        TRUE ~ "Non-GC"
      ),
      patient = HH,
      clone = clone,
      clone_nr = clone_nr,
      clone_full_name = clone_full_name,
      .before = 1
    )

  # ---- Permutation tests: BACH2+ vs other GC B cells of the same tree ----

  test_cells <- cells %>% filter(group %in% c("BACH2+", "Other GC B cells"))
  is_BACH2 <- test_cells$group == "BACH2+"

  tree_stats <- tibble(
    patient = HH,
    clone = clone,
    clone_nr = clone_nr,
    clone_full_name = clone_full_name,
    n_nodes = length(node_labels),
    n_BACH2 = sum(is_BACH2),
    n_other_GC = sum(!is_BACH2)
  )

  if (sum(is_BACH2) >= min_cells_group && sum(!is_BACH2) >= min_cells_group && !is.na(max_depth) && max_depth > 0){

    node_dist <- dist.nodes(tree)
    rownames(node_dist) <- colnames(node_dist) <- node_labels

    rel_depth <- test_cells$rel_depth
    cell_nodes <- test_cells$node

    get_stats <- function(lbl){
      c(
        # < 0: BACH2+ cells closer to GL than other GC B cells
        delta_median_rel_depth = median(rel_depth[lbl]) - median(rel_depth[!lbl]),
        # Mean patristic distance between BACH2+ cells (cell level, so same-node pairs = 0)
        mpd_BACH2 = mean(as.dist(node_dist[cell_nodes[lbl], cell_nodes[lbl], drop = FALSE]))
      )
    }

    observed <- get_stats(is_BACH2)
    null <- replicate(n_perm, get_stats(sample(is_BACH2)))

    tree_stats <- tree_stats %>%
      mutate(
        delta_median_rel_depth = observed[["delta_median_rel_depth"]],
        p_closer_to_GL = perm_p_lower(observed[["delta_median_rel_depth"]], null["delta_median_rel_depth", ]),
        mpd_BACH2 = observed[["mpd_BACH2"]],
        mpd_null_mean = mean(null["mpd_BACH2", ]),
        p_clustered = perm_p_lower(observed[["mpd_BACH2"]], null["mpd_BACH2", ])
      )
  }

  list(cells = cells, tree_stats = tree_stats)
}

# ------------------------------------------------------------------------------
# Run all trees
# ------------------------------------------------------------------------------

results <- map(samples, analyse_tree) %>% compact()

cells_all <- map_dfr(results, "cells")
tree_stats_all <- map_dfr(results, "tree_stats")

write_delim(cells_all, glue("{outdir}/BACH2_tree_position_cells.tsv"), delim = "\t")
write_delim(tree_stats_all, glue("{outdir}/BACH2_tree_position_per_tree.tsv"), delim = "\t")

# Number of trees, and how many show BACH2+ closer to GL / clustered (p < 0.05)
tree_stats_all %>%
  group_by(patient) %>%
  summarise(
    n_trees = n(),
    n_tested = sum(!is.na(p_closer_to_GL)),
    n_closer_to_GL = sum(p_closer_to_GL < 0.05, na.rm = TRUE),
    n_clustered = sum(p_clustered < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Across trees: per-tree summaries, BACH2+ vs other GC B cells
# ------------------------------------------------------------------------------

metric_labels <- c(
  median_rel_depth = "Median relative depth\n(0 = GL, 1 = deepest node)",
  median_mut_from_GL = "Median mutations from GL",
  frac_internal = "Fraction of cells on internal nodes",
  mean_n_desc_nodes = "Mean N descendant nodes"
)

df_per_tree <- cells_all %>%
  filter(group %in% c("BACH2+", "Other GC B cells")) %>%
  group_by(patient, clone, clone_nr, group) %>%
  summarise(
    n_cells = n(),
    median_rel_depth = median(rel_depth),
    median_mut_from_GL = median(mut_from_GL),
    frac_internal = mean(is_internal),
    mean_n_desc_nodes = mean(n_desc_nodes),
    .groups = "drop"
  ) %>%
  # Keep trees with enough cells in both groups
  group_by(patient, clone) %>%
  filter(n() == 2, all(n_cells >= min_cells_group)) %>%
  ungroup() %>%
  pivot_longer(cols = all_of(names(metric_labels)), names_to = "metric", values_to = "value")

# Paired Wilcoxon test per patient and metric
# Paired because BACH2+ and other GC B cells come from the same tree
test_results <- df_per_tree %>%
  select(patient, clone, metric, group, value) %>%
  pivot_wider(names_from = group, values_from = value) %>%
  filter(!is.na(`BACH2+`), !is.na(`Other GC B cells`)) %>%
  group_by(patient, metric) %>%
  summarise(
    n_pairs = n(),
    p_value = if (n_pairs >= 3) wilcox.test(`BACH2+`, `Other GC B cells`, paired = TRUE, exact = FALSE)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    is.na(p_value) ~ glue("n = {n_pairs}"),
    p_value < 0.001 ~ glue("p < 0.001\nn = {n_pairs}"),
    TRUE ~ glue("p = {signif(p_value, 2)}\nn = {n_pairs}")
  ))

print(test_results)
write_delim(test_results, glue("{outdir}/BACH2_tree_position_wilcox.tsv"), delim = "\t")

# Prep label with stats
label_pos <- df_per_tree %>%
  group_by(metric) %>%
  mutate(y = max(value, na.rm = TRUE) * 1.12) %>%
  group_by(patient, metric) %>%
  summarise(y = first(y), .groups = "drop") %>%
  left_join(test_results, by = c("patient", "metric"))

# Plot
df_per_tree %>%
  ggplot(aes(x = patient, y = value, color = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.5) +
  geom_text(data = label_pos, aes(x = patient, y = y, label = p_label), inherit.aes = FALSE, size = 3) +
  facet_wrap(~ metric, scales = "free_y", labeller = as_labeller(metric_labels)) +
  scale_color_manual(values = c("BACH2+" = "#D55E00", "Other GC B cells" = "grey40")) +
  labs(
    title = "Position of BACH2+ cells in GCtrees",
    subtitle = "One point per tree. P-values: paired Wilcoxon test, BACH2+ vs. other GC B cells",
    x = "Patient",
    y = NULL,
    color = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(glue("{outdir}/BACH2_tree_position_boxplot.png"), width = 10, height = 8)

# ------------------------------------------------------------------------------
# Per-tree permutation results
# ------------------------------------------------------------------------------

tree_stats_all %>%
  filter(!is.na(p_closer_to_GL)) %>%
  mutate(significant = if_else(p_closer_to_GL < 0.05, "p < 0.05", "n.s.")) %>%
  ggplot(aes(x = patient, y = delta_median_rel_depth, color = significant, size = n_BACH2)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
  scale_color_manual(values = c("p < 0.05" = "#D55E00", "n.s." = "grey50")) +
  labs(
    title = "BACH2+ vs. other GC B cells: relative depth per tree",
    subtitle = glue("Below 0: BACH2+ cells closer to GL. P-value: one-sided label permutation test ({n_perm} permutations)"),
    x = "Patient",
    y = "Median relative depth, BACH2+ minus other GC B cells",
    color = "Permutation test",
    size = "N BACH2+ cells"
  ) +
  theme_bw()

ggsave(glue("{outdir}/BACH2_tree_position_permutation_depth.png"), width = 8, height = 6)

tree_stats_all %>%
  filter(!is.na(p_clustered)) %>%
  mutate(
    mpd_ratio = mpd_BACH2 / mpd_null_mean,
    significant = if_else(p_clustered < 0.05, "p < 0.05", "n.s.")
  ) %>%
  ggplot(aes(x = patient, y = mpd_ratio, color = significant, size = n_BACH2)) +
  geom_hline(yintercept = 1, color = "grey") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
  scale_color_manual(values = c("p < 0.05" = "#D55E00", "n.s." = "grey50")) +
  labs(
    title = "Clustering of BACH2+ cells per tree",
    subtitle = glue("Below 1: BACH2+ cells closer to each other than random cells. P-value: one-sided label permutation test ({n_perm} permutations)"),
    x = "Patient",
    y = "Mean pairwise distance BACH2+ / permutation mean",
    color = "Permutation test",
    size = "N BACH2+ cells"
  ) +
  theme_bw()

ggsave(glue("{outdir}/BACH2_tree_position_permutation_clustering.png"), width = 8, height = 6)
