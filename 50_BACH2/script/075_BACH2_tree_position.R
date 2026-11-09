library(tidyverse)
library(glue)
library(ape)

# ==============================================================================
# Where in the GCtrees are BACH2+ and BACH2_cycling cells, compared with the
# other cell subtypes?
#
# 1) Tree position per cell subtype (one subtype vs. all other cells of the tree)
#      mut_from_GL   : mutations from the germline root (sum of branch lengths)
#      rel_depth     : mut_from_GL / deepest node in the tree (0 = GL, 1 = deepest)
#      is_internal   : node is an observed ancestor (internal node) rather than a tip
#      node_n_cells  : number of cells sharing the node genotype
#      n_desc_nodes  : number of nodes below the node
#
# 2) Upstream / downstream: which subtypes sit on the observed ancestor and
#    observed descendant nodes of each subtype
#      direct : nearest observed ancestor -> node (as in 06_analyse_trees.R)
#      any    : any observed ancestor    -> node
#
# Cells identical to the germline sit on the GL root node (deduplicate --root GL).
# They are kept (mut_from_GL = 0), so unmutated cells such as naive B cells are included.
#
# 3) Clonal bursts: observed nodes with a large parental population
#    (>= burst_min_parent_cells cells) and >= burst_min_children child nodes,
#    following Pae et al. 2025 (Nature). Counts how often each subtype is the
#    burst parent (cells on the bursting node) or burst offspring (cells on its
#    direct child nodes)
# ==============================================================================

set.seed(1)

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

# resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/06_resolve_LC_germlined.rds")
resolve_LC_list_germlined <- readRDS("50_BACH2/out/resolve_LC_list_germlined_L1_L3.rds")
resolve_LC_list_germlined <- resolve_LC_list_germlined[names(resolve_LC_list_germlined) != "HH151"]

# Named list: clone_full_name -> named vector (GCtree sequence id -> sequence_alignment)
seq_dir <- readRDS("50_BACH2/out/rds/seq_dir.rds")

tree_dir <- "50_BACH2/plot"
idmap_dir <- "50_BACH2/out"

outdir <- "50_BACH2/plot/075_tree_position"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Cell subtype column. Every value is analysed as its own subtype
celltype_col <- "L1_L3_annotation"

# Subtypes shown first in plots
focus_types <- c("BACH2+", "BACH2_cycling")

n_perm <- 1000
min_cells_type <- 1        # min cells of a subtype (and of the rest) in a tree to include it
min_trees_type <- 2        # min trees a subtype must be present in to be shown in the position plots
# Clonal burst: >= burst_min_parent_cells cells on the node AND >= burst_min_children child nodes
burst_min_children <- 3
burst_min_parent_cells_values <- c(10, 20)   # each value is run separately, own output folder
burst_exclude_GL <- TRUE                     # the GL node holds unmutated cells, not a GC expansion

# All clones with a GCtree output folder
samples <- list.dirs(tree_dir, full.names = FALSE, recursive = FALSE) %>%
  str_subset("^HH[0-9]+_clone_nr_[0-9]+_clone_")

samples <- samples[str_detect(samples, "HH151", negate = TRUE)]

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

# All observed ancestors of every observed node, with the number of observed
# steps between them (step 1 = nearest observed ancestor)
get_observed_ancestors <- function(tree, node_labels, is_observed){

  parent_idx <- integer(length(node_labels))
  parent_idx[tree$edge[, 2]] <- tree$edge[, 1]

  map_dfr(which(is_observed), function(i){

    ancestors <- integer(0)
    current <- parent_idx[i]

    while (current != 0){
      if (is_observed[current]) ancestors <- c(ancestors, current)
      current <- parent_idx[current]
    }

    if (length(ancestors) == 0) return(NULL)

    tibble(
      descendant = node_labels[i],
      ancestor = node_labels[ancestors],
      step = seq_along(ancestors)
    )
  })
}

# Permutation p-values
perm_p <- function(observed, null, tail = c("lower", "upper")){
  tail <- match.arg(tail)
  null <- null[!is.na(null)]
  if (is.na(observed) || length(null) == 0) return(NA_real_)
  k <- if (tail == "lower") sum(null <= observed) else sum(null >= observed)
  (k + 1) / (length(null) + 1)
}

# ------------------------------------------------------------------------------
# Per tree
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

  tree_ids <- tibble(patient = HH, clone = clone, clone_nr = clone_nr, clone_full_name = clone_full_name)

  # ---- Tree and node positions ----

  tree <- read.tree(tree_file)
  node_labels <- c(tree$tip.label, tree$node.label)

  # The root is the germline. Cells identical to GL are mapped to it by idmap
  node_labels[Ntip(tree) + 1] <- "GL"

  depth <- node.depth.edgelength(tree)
  max_depth <- max(depth)
  n_children <- tabulate(tree$edge[, 1], nbins = length(node_labels))

  node_df <- tibble(
    node = node_labels,
    mut_from_GL = depth,
    rel_depth = if (max_depth > 0) depth / max_depth else NA_real_,
    is_internal = seq_along(node_labels) > Ntip(tree),
    n_children = n_children,
    n_desc_nodes = count_descendants(tree)
  )

  # ---- Cells ----

  # Keep the GL node: it holds the cells identical to the germline. Drop the GL sequence itself
  idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id_gctree"), show_col_types = FALSE) %>%
    separate_longer_delim(sequence_id_gctree, delim = ":") %>%
    filter(sequence_id_gctree != "GL")

  sequence_ids <- seq_dir[[clone_full_name]]

  df_clone <- resolve_LC_list_germlined[[HH]] %>%
    filter(
      locus == "IGH",
      clone_subgroup_id_90_similarity == clone
    ) %>%
    mutate(
      sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)],
      cell_type = .data[[celltype_col]]
    )

  cells <- df_clone %>%
    filter(!is.na(cell_type)) %>%
    select(sequence_id, sequence_id_gctree, cell_type) %>%
    inner_join(idmap, by = "sequence_id_gctree", relationship = "many-to-one") %>%
    inner_join(node_df %>% filter(str_detect(node, "^seq[0-9]+$") | node == "GL"), by = "node") %>%
    add_count(node, name = "node_n_cells")

  # How many cells of each subtype in the clone end up in the tree (and on GL)
  type_mapping <- df_clone %>%
    dplyr::count(cell_type, name = "n_cells_clone") %>%
    left_join(cells %>% dplyr::count(cell_type, name = "n_cells_tree"), by = "cell_type") %>%
    left_join(cells %>% filter(node == "GL") %>% dplyr::count(cell_type, name = "n_cells_on_GL"), by = "cell_type") %>%
    mutate(across(c(n_cells_tree, n_cells_on_GL), ~ replace_na(.x, 0L)))
  type_mapping <- bind_cols(tree_ids[rep(1, nrow(type_mapping)), ], type_mapping)

  if (nrow(cells) == 0) return(list(type_mapping = type_mapping))

  # Observed nodes: nodes holding >= 1 cell (seq nodes and, if it has cells, GL)
  is_observed <- node_labels %in% cells$node

  cells <- bind_cols(tree_ids, cells)

  # Cell subtypes present on each node
  node_types <- cells %>%
    distinct(node, cell_type)

  # ---- 1) Position: each subtype vs. all other cells of the tree ----

  node_dist <- dist.nodes(tree)
  rownames(node_dist) <- colnames(node_dist) <- node_labels

  type_counts <- table(cells$cell_type)
  types_to_test <- names(type_counts)[type_counts >= min_cells_type & (nrow(cells) - type_counts) >= min_cells_type]

  position_tests <- map_dfr(types_to_test, function(this_type){

    is_type <- cells$cell_type == this_type
    rel_depth <- cells$rel_depth
    cell_nodes <- cells$node

    get_stats <- function(lbl){
      c(
        # < 0: subtype closer to GL than the other cells of the tree
        delta_median_rel_depth = median(rel_depth[lbl]) - median(rel_depth[!lbl]),
        # Mean patristic distance between cells of the subtype (same-node pairs = 0)
        mpd = mean(as.dist(node_dist[cell_nodes[lbl], cell_nodes[lbl], drop = FALSE]))
      )
    }

    observed <- get_stats(is_type)
    null <- replicate(n_perm, get_stats(sample(is_type)))

    tibble(
      cell_type = this_type,
      n_type = sum(is_type),
      n_rest = sum(!is_type),
      delta_median_rel_depth = observed[["delta_median_rel_depth"]],
      p_closer_to_GL = perm_p(observed[["delta_median_rel_depth"]], null["delta_median_rel_depth", ], "lower"),
      p_further_from_GL = perm_p(observed[["delta_median_rel_depth"]], null["delta_median_rel_depth", ], "upper"),
      mpd = observed[["mpd"]],
      mpd_null_mean = mean(null["mpd", ], na.rm = TRUE),
      p_clustered = perm_p(observed[["mpd"]], null["mpd", ], "lower")
    )
  })

  if (nrow(position_tests) > 0) position_tests <- bind_cols(tree_ids, position_tests)

  # ---- 2) Upstream / downstream relations between observed nodes ----

  relations <- get_observed_ancestors(tree, node_labels, is_observed)

  if (nrow(relations) > 0){
    relations <- relations %>%
      inner_join(node_types %>% dplyr::rename(ancestor = node, ancestor_type = cell_type),
                 by = "ancestor", relationship = "many-to-many") %>%
      inner_join(node_types %>% dplyr::rename(descendant = node, descendant_type = cell_type),
                 by = "descendant", relationship = "many-to-many")
    relations <- bind_cols(tree_ids[rep(1, nrow(relations)), ], relations)
  }

  # ---- 3) Clonal burst candidates: nodes with >= burst_min_children children ----
  # The parent size threshold is applied later, so several thresholds can be compared

  burst_idx <- which(n_children >= burst_min_children)

  burst_cells <- map_dfr(burst_idx, function(i){

    child_nodes <- node_labels[tree$edge[tree$edge[, 1] == i, 2]]

    bind_rows(
      cells %>% filter(node == node_labels[i]) %>% mutate(burst_role = "Parent"),
      cells %>% filter(node %in% child_nodes) %>% mutate(burst_role = "Offspring")
    ) %>%
      mutate(
        burst_node = node_labels[i],
        burst_n_children = n_children[i],
        burst_parent_n_cells = sum(cells$node == node_labels[i])
      )
  })

  if (nrow(burst_cells) == 0){
    burst_cells <- cells[0, ] %>%
      mutate(burst_role = character(), burst_node = character(), burst_n_children = integer(), burst_parent_n_cells = integer())
  }

  bursts <- tibble(
    burst_node = node_labels[burst_idx],
    burst_n_children = n_children[burst_idx],
    burst_parent_n_cells = map_int(node_labels[burst_idx], ~ sum(cells$node == .x)),
    burst_mut_from_GL = depth[burst_idx]
  ) %>%
    left_join(
      burst_cells %>%
        group_by(burst_node) %>%
        summarise(
          n_cells_parent = sum(burst_role == "Parent"),
          n_cells_offspring = sum(burst_role == "Offspring"),
          offspring_types = paste(sort(unique(cell_type[burst_role == "Offspring"])), collapse = ", "),
          .groups = "drop"
        ),
      by = "burst_node"
    ) %>%
    mutate(across(c(n_cells_parent, n_cells_offspring), ~ replace_na(.x, 0L)))

  if (nrow(bursts) > 0) bursts <- bind_cols(tree_ids[rep(1, nrow(bursts)), ], bursts)

  tree_summary <- tree_ids %>%
    mutate(
      n_cells = nrow(cells),
      n_nodes = length(node_labels),
      n_observed_nodes = sum(is_observed),
      n_cells_on_GL = sum(cells$node == "GL"),
      max_mut_from_GL = max_depth,
      n_burst_candidates = length(burst_idx)
    )

  list(
    cells = cells,
    position_tests = position_tests,
    relations = relations,
    burst_cells = burst_cells,
    bursts = bursts,
    tree_summary = tree_summary,
    type_mapping = type_mapping
  )
}

# ------------------------------------------------------------------------------
# Run all trees
# ------------------------------------------------------------------------------

results <- map(samples, analyse_tree) %>% compact()

cells_all <- map_dfr(results, "cells")
position_tests_all <- map_dfr(results, "position_tests")
relations_all <- map_dfr(results, "relations")
burst_cells_all <- map_dfr(results, "burst_cells")
bursts_all <- map_dfr(results, "bursts")
tree_summary_all <- map_dfr(results, "tree_summary")
type_mapping_all <- map_dfr(results, "type_mapping")

write_delim(cells_all, glue("{outdir}/cells.tsv"), delim = "\t")
write_delim(tree_summary_all, glue("{outdir}/tree_summary.tsv"), delim = "\t")

# Subtype order: focus subtypes first, then by number of cells
celltype_order <- cells_all %>%
  dplyr::count(cell_type, sort = TRUE) %>%
  pull(cell_type) %>%
  { c(intersect(focus_types, .), setdiff(., focus_types)) }

# How many cells of each subtype are in the clones, in the trees, and on the GL node.
# A subtype with n_cells_clone > 0 but n_cells_tree = 0 is missing from the FASTA / seq_dir
type_mapping_overview <- type_mapping_all %>%
  group_by(cell_type) %>%
  summarise(
    n_trees = sum(n_cells_tree > 0),
    across(c(n_cells_clone, n_cells_tree, n_cells_on_GL), sum),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells_clone))

print(type_mapping_overview, n = Inf)
write_delim(type_mapping_overview, glue("{outdir}/celltype_cells_in_trees.tsv"), delim = "\t")

# Overview of subtypes across trees
cells_all %>%
  group_by(cell_type) %>%
  summarise(n_cells = n(), n_trees = n_distinct(clone_full_name), .groups = "drop") %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order)) %>%
  arrange(cell_type)

# ==============================================================================
# 1) Tree position per subtype
# ==============================================================================

outdir_position <- glue("{outdir}/position")
dir.create(outdir_position, recursive = TRUE, showWarnings = FALSE)

write_delim(position_tests_all, glue("{outdir_position}/position_permutation_per_tree.tsv"), delim = "\t")

# Number of trees where each subtype is closer to / further from GL, or clustered (p < 0.05)
position_overview <- position_tests_all %>%
  group_by(cell_type) %>%
  summarise(
    n_trees_tested = n(),
    n_closer_to_GL = sum(p_closer_to_GL < 0.05),
    n_further_from_GL = sum(p_further_from_GL < 0.05),
    n_clustered = sum(p_clustered < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order)) %>%
  arrange(cell_type)

print(position_overview)
write_delim(position_overview, glue("{outdir_position}/position_permutation_overview.tsv"), delim = "\t")

# ---- Per-tree summaries: subtype vs. rest of the tree ----

metric_labels <- c(
  median_rel_depth = "Median relative depth\n(0 = GL, 1 = deepest node)",
  median_mut_from_GL = "Median mutations from GL",
  frac_internal = "Fraction of cells on internal nodes",
  mean_n_desc_nodes = "Mean N descendant nodes"
)

summarise_metrics <- function(df){
  df %>%
    summarise(
      n_cells = n(),
      median_rel_depth = median(rel_depth),
      median_mut_from_GL = median(mut_from_GL),
      frac_internal = mean(is_internal),
      mean_n_desc_nodes = mean(n_desc_nodes),
      .groups = "drop"
    )
}

df_per_tree <- map_dfr(celltype_order, function(this_type){
  cells_all %>%
    mutate(
      cell_type_focal = this_type,
      group = if_else(cell_type == this_type, "Subtype", "Rest of tree")
    ) %>%
    group_by(patient, clone_full_name, cell_type_focal, group) %>%
    summarise_metrics()
}) %>%
  # Keep trees with enough cells of the subtype and of the rest
  group_by(clone_full_name, cell_type_focal) %>%
  filter(n() == 2, all(n_cells >= min_cells_type)) %>%
  group_by(cell_type_focal) %>%
  filter(n_distinct(clone_full_name) >= min_trees_type) %>%
  ungroup() %>%
  pivot_longer(cols = all_of(names(metric_labels)), names_to = "metric", values_to = "value") %>%
  mutate(
    cell_type_focal = factor(cell_type_focal, levels = celltype_order),
    group = factor(group, levels = c("Subtype", "Rest of tree"))
  )

# Paired Wilcoxon test: subtype vs. rest of the same tree, across trees
# All patients together, and per patient
run_wilcox <- function(df, ...){
  df %>%
    select(clone_full_name, ..., metric, group, value) %>%
    pivot_wider(names_from = group, values_from = value) %>%
    group_by(..., metric) %>%
    summarise(
      n_trees = n(),
      median_diff = median(Subtype - `Rest of tree`),
      p_value = if (n_trees >= 3) wilcox.test(Subtype, `Rest of tree`, paired = TRUE, exact = FALSE)$p.value else NA_real_,
      .groups = "drop"
    ) %>%
    group_by(metric) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

test_results <- run_wilcox(df_per_tree, cell_type_focal)
test_results_patient <- run_wilcox(df_per_tree, patient, cell_type_focal)

print(test_results)
write_delim(test_results, glue("{outdir_position}/position_wilcox_all_patients.tsv"), delim = "\t")
write_delim(test_results_patient, glue("{outdir_position}/position_wilcox_per_patient.tsv"), delim = "\t")

# Plot
# P-value labels per patient and metric
label_pos_patient <- df_per_tree %>%
  group_by(metric) %>%
  summarise(y = max(value, na.rm = TRUE) * 1.1, .groups = "drop") %>%
  right_join(test_results_patient, by = "metric") %>%
  mutate(p_label = case_when(
    is.na(p_adj) ~ "",
    p_adj < 0.001 ~ "p < 0.001",
    TRUE ~ glue("p = {signif(p_adj, 2)}")
  )) %>% 
  filter(cell_type_focal != "GC_B_cells")

df_per_tree %>%
  filter(cell_type_focal != "GC_B_cells") %>% 
  ggplot(aes(x = cell_type_focal, y = value, color = group)) +
  # geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  # geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.5) +
  stat_summary(
    fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", width = 0.5, linewidth = 0.4,
    position = position_dodge(width = 0.75)
  ) +
  geom_text(data = label_pos_patient, aes(x = cell_type_focal, y = y, label = p_label), inherit.aes = FALSE, size = 2.5) +
  facet_grid(
    rows = vars(metric),
    cols = vars(patient),
    scales = "free_y",
    labeller = labeller(metric = as_labeller(metric_labels))
  ) +
  scale_color_manual(values = c("Subtype" = "#D55E00", "Rest of tree" = "grey40")) +
  labs(
    title = "Position of cell subtypes in GCtrees, per patient",
    subtitle = "One point per tree.\nP-values: paired Wilcoxon test per patient, subtype vs. rest of the same tree (BH adjusted per metric)",
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(size = 7),
    legend.position = "bottom"
  )

ggsave(
  glue("{outdir_position}/position_boxplot_per_patient.png"),
  width = max(12, 0.9 * length(celltype_order) * n_distinct(df_per_tree$patient)),
  height = 14
)


# ---- Relative depth per tree: subtype vs. rest of the tree ----
min_cells_perm_plot <- 3   # min cells of the subtype in a tree to show its permutation result

position_tests_all %>%
  filter(
    cell_type %in% unique(df_per_tree$cell_type_focal),
    n_type >= min_cells_perm_plot
  ) %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    result = case_when(
      p_closer_to_GL < 0.05 ~ "Closer to GL (p < 0.05)",
      p_further_from_GL < 0.05 ~ "Further from GL (p < 0.05)",
      TRUE ~ "n.s."
    )
  ) %>%
  ggplot(aes(x = cell_type, y = delta_median_rel_depth)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(aes(color = result, size = n_type), width = 0.15, height = 0, alpha = 0.7) +
  stat_summary(
    fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", width = 0.5, linewidth = 0.4
  ) +
  facet_grid(rows = vars(patient)) +
  scale_color_manual(values = c("Closer to GL (p < 0.05)" = "#0072B2", "Further from GL (p < 0.05)" = "#D55E00", "n.s." = "grey60")) +
  labs(
    title = "Relative depth per tree: subtype vs. rest of the tree",
    subtitle = glue(
      "Below 0: subtype closer to GL. Line: median over trees. Trees with >= {min_cells_perm_plot} cells of the subtype.\n",
      "One-sided label permutation tests ({n_perm} permutations)"
    ),
    x = NULL,
    y = "Median relative depth, subtype minus rest",
    color = "Permutation test",
    size = "N cells of subtype"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  glue("{outdir_position}/position_permutation_depth.png"),
  width = max(8, 1.2 * length(celltype_order)),
  height = 3 * n_distinct(position_tests_all$patient) + 2
)

# ---- Clustering of each subtype within its tree ----

position_tests_all %>%
  filter(
    cell_type %in% unique(df_per_tree$cell_type_focal),
    n_type >= min_cells_perm_plot,
    !is.na(p_clustered)
  ) %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    mpd_ratio = mpd / mpd_null_mean,
    result = if_else(p_clustered < 0.05, "Clustered (p < 0.05)", "n.s.")
  ) %>%
  ggplot(aes(x = cell_type, y = mpd_ratio)) +
  geom_hline(yintercept = 1, color = "grey") +
  geom_jitter(aes(color = result, size = n_type), width = 0.15, height = 0, alpha = 0.7) +
  stat_summary(
    fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", width = 0.5, linewidth = 0.4
  ) +
  facet_grid(rows = vars(patient)) +
  scale_color_manual(values = c("Clustered (p < 0.05)" = "#D55E00", "n.s." = "grey60")) +
  labs(
    title = "Clustering of each subtype within its tree",
    subtitle = glue(
      "Below 1: cells of the subtype closer to each other than random cells. Line: median over trees.\n",
      "Trees with >= {min_cells_perm_plot} cells of the subtype. One-sided label permutation test ({n_perm} permutations)"
    ),
    x = NULL,
    y = "Mean pairwise distance / permutation mean",
    color = "Permutation test",
    size = "N cells of subtype"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  glue("{outdir_position}/position_permutation_clustering.png"),
  width = max(8, 1.2 * length(celltype_order)),
  height = 3 * n_distinct(position_tests_all$patient) + 2
)

# ==============================================================================
# 2) Upstream / downstream subtypes
# ==============================================================================
#
# Counted per pair of observed nodes (ancestor node, descendant node), using the
# unique subtypes present on each node (as in 06_analyse_trees.R).
#   Upstream of X   : subtypes on the ancestor nodes of nodes that contain X
#   Downstream of X : subtypes on the descendant nodes of nodes that contain X
#
# Permutation test: within each tree, the subtype composition of the observed
# nodes is shuffled between node positions (each node keeps its own set of
# subtypes; the GL node is kept fixed), and the pair counts are recomputed.
# Expected = mean over permutations. This keeps the tree structure and the
# subtypes present in each clone, and asks whether the subtypes on ancestor and
# descendant nodes are related. Shuffling single cells instead would break up
# nodes whose cells share a subtype and inflate the expected counts.
#
# Subtypes in relation_exclude_types are kept on the tree (so they still occupy
# nodes, in the observed data and in the permutations) but removed from the
# pair counts and plots.

outdir_relations <- glue("{outdir}/upstream_downstream")
dir.create(outdir_relations, recursive = TRUE, showWarnings = FALSE)

relation_exclude_types <- c("GC_B_cells")
n_perm_relations <- 1000

celltype_order_rel <- setdiff(celltype_order, relation_exclude_types)

# Tree structure: all (ancestor node, descendant node) pairs between observed nodes
node_pairs <- relations_all %>%
  distinct(clone_full_name, ancestor, descendant, step)

# Unique subtypes per observed node
node_types <- cells_all %>%
  distinct(clone_full_name, node, cell_type)

# Observed nodes that are shuffled (GL stays at the root)
nodes_to_shuffle <- node_types %>%
  distinct(clone_full_name, node) %>%
  filter(node != "GL")

# Count subtype pairs, for nearest observed ancestor ("direct") and any observed ancestor ("any")
count_relation_pairs <- function(node_types){
  
  pairs <- node_pairs %>%
    inner_join(
      node_types %>% dplyr::rename(ancestor = node, ancestor_type = cell_type),
      by = c("clone_full_name", "ancestor"), relationship = "many-to-many"
    ) %>%
    inner_join(
      node_types %>% dplyr::rename(descendant = node, descendant_type = cell_type),
      by = c("clone_full_name", "descendant"), relationship = "many-to-many"
    ) %>%
    filter(
      !ancestor_type %in% relation_exclude_types,
      !descendant_type %in% relation_exclude_types
    )
  
  bind_rows(
    pairs %>% filter(step == 1) %>% mutate(relation = "direct"),
    pairs %>% mutate(relation = "any")
  ) %>%
    dplyr::count(relation, ancestor_type, descendant_type, name = "n")
}

relation_grid <- expand_grid(
  relation = c("direct", "any"),
  ancestor_type = celltype_order_rel,
  descendant_type = celltype_order_rel
)

# Observed
observed_pairs <- count_relation_pairs(node_types) %>%
  right_join(relation_grid, by = c("relation", "ancestor_type", "descendant_type")) %>%
  mutate(n = replace_na(n, 0L))

# Permutations: move each node's set of subtypes to a random observed node position of the same tree
null_pairs <- map_dfr(seq_len(n_perm_relations), function(i){
  
  if (i %% 100 == 0) message(glue("Permutation {i}/{n_perm_relations}"))
  
  new_positions <- nodes_to_shuffle %>%
    group_by(clone_full_name) %>%
    mutate(new_node = if (n() > 1) sample(node) else node) %>%
    ungroup()
  
  node_types %>%
    left_join(new_positions, by = c("clone_full_name", "node")) %>%
    mutate(node = coalesce(new_node, node)) %>%
    select(-new_node) %>%
    count_relation_pairs() %>%
    mutate(perm = i)
})

null_pairs <- null_pairs %>%
  right_join(expand_grid(relation_grid, perm = seq_len(n_perm_relations)),
             by = c("relation", "ancestor_type", "descendant_type", "perm")) %>%
  mutate(n = replace_na(n, 0L)) %>%
  dplyr::rename(n_null = n)

# Permutation p-values per pair of subtypes
perm_results <- null_pairs %>%
  left_join(observed_pairs, by = c("relation", "ancestor_type", "descendant_type")) %>%
  group_by(relation, ancestor_type, descendant_type) %>%
  summarise(
    n = first(n),
    expected_perm = mean(n_null),
    p_enriched = (sum(n_null >= n) + 1) / (n_perm_relations + 1),
    p_depleted = (sum(n_null <= n) + 1) / (n_perm_relations + 1),
    .groups = "drop"
  ) %>%
  mutate(
    log2_OE_perm = log2((n + 0.5) / (expected_perm + 0.5)),
    p_two_sided = pmin(1, 2 * pmin(p_enriched, p_depleted))
  ) %>%
  group_by(relation) %>%
  mutate(p_adj = p.adjust(p_two_sided, method = "BH")) %>%
  ungroup()

write_delim(perm_results, glue("{outdir_relations}/pair_permutation_results.tsv"), delim = "\t")

relation_titles <- c(
  direct = "nearest observed ancestor",
  any = "any observed ancestor"
)

plot_relation_heatmap <- function(df, fill_col, fill_label, title, subtitle, diverging = FALSE){
  
  p <- df %>%
    ggplot(aes(x = col_type, y = row_type, fill = .data[[fill_col]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 3) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(limits = rev, drop = FALSE) +
    labs(title = title, subtitle = subtitle, fill = fill_label) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank())
  
  if (diverging){
    lim <- max(abs(df[[fill_col]]), na.rm = TRUE)
    p + scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-lim, lim))
  } else {
    p + scale_fill_gradient(low = "white", high = "steelblue")
  }
}

n_size <- 0.8 * length(celltype_order_rel) + 3

for (relation_name in names(relation_titles)){
  
  # relation_name <- "direct"
  
  pair_counts <- perm_results %>%
    filter(relation == relation_name) %>%
    mutate(
      expected_marginal = sum(n) *
        (ave(n, ancestor_type, FUN = sum) / sum(n)) *
        (ave(n, descendant_type, FUN = sum) / sum(n)),
      log2_OE_marginal = log2((n + 0.5) / (expected_marginal + 0.5))
    )
  
  # Number of trees each pair is seen in
  pair_trees <- relations_all %>%
    filter(
      relation_name == "any" | step == 1,
      !ancestor_type %in% relation_exclude_types,
      !descendant_type %in% relation_exclude_types
    ) %>%
    distinct(clone_full_name, ancestor_type, descendant_type) %>%
    dplyr::count(ancestor_type, descendant_type, name = "n_trees")
  
  pair_counts <- pair_counts %>%
    left_join(pair_trees, by = c("ancestor_type", "descendant_type")) %>%
    mutate(n_trees = replace_na(n_trees, 0L))
  
  write_delim(pair_counts, glue("{outdir_relations}/pair_counts_{relation_name}.tsv"), delim = "\t")
  
  # Upstream: rows = subtype, columns = subtypes found upstream (fraction per row)
  upstream <- pair_counts %>%
    group_by(descendant_type) %>%
    mutate(frac = if (sum(n) > 0) n / sum(n) else 0) %>%
    ungroup() %>%
    transmute(
      row_type = factor(descendant_type, levels = celltype_order_rel),
      col_type = factor(ancestor_type, levels = celltype_order_rel),
      frac, n,
      label = if_else(n > 0, glue("{round(100 * frac)}%\n({n})"), "")
    )
  
  # Downstream: rows = subtype, columns = subtypes found downstream (fraction per row)
  downstream <- pair_counts %>%
    group_by(ancestor_type) %>%
    mutate(frac = if (sum(n) > 0) n / sum(n) else 0) %>%
    ungroup() %>%
    transmute(
      row_type = factor(ancestor_type, levels = celltype_order_rel),
      col_type = factor(descendant_type, levels = celltype_order_rel),
      frac, n, log2_OE_marginal, log2_OE_perm, p_adj,
      label = if_else(n > 0, glue("{round(100 * frac)}%\n({n})"), "")
    )
  
  plot_relation_heatmap(
    upstream, "frac", "Fraction of row",
    title = glue("Subtypes upstream of each subtype ({relation_titles[[relation_name]]})"),
    subtitle = "Rows: subtype on the descendant node. Columns: subtype on the ancestor node.\nLabel: % of row (N node pairs)"
  ) +
    labs(x = "Upstream subtype", y = "Subtype")
  
  ggsave(glue("{outdir_relations}/upstream_{relation_name}.png"), width = n_size + 2, height = n_size)
  
  plot_relation_heatmap(
    downstream, "frac", "Fraction of row",
    title = glue("Subtypes downstream of each subtype ({relation_titles[[relation_name]]})"),
    subtitle = "Rows: subtype on the ancestor node. Columns: subtype on the descendant node.\nLabel: % of row (N node pairs)"
  ) +
    labs(x = "Downstream subtype", y = "Subtype")
  
  ggsave(glue("{outdir_relations}/downstream_{relation_name}.png"), width = n_size + 2, height = n_size)
  
  # Observed / expected from the subtype totals (no permutation)
  downstream %>%
    mutate(label = if_else(n > 0, as.character(n), "")) %>%
    plot_relation_heatmap(
      "log2_OE_marginal", "log2(O/E)",
      title = glue("Ancestor -> descendant subtype pairs vs. subtype totals ({relation_titles[[relation_name]]})"),
      subtitle = "Expected from the row and column totals.\nRed: pair seen more often than expected. Label: N node pairs",
      diverging = TRUE
    ) +
    labs(x = "Descendant subtype", y = "Ancestor subtype")
  
  ggsave(glue("{outdir_relations}/ancestor_descendant_OE_{relation_name}.png"), width = n_size + 2, height = n_size)
  
  # Observed / expected from the within-tree permutations, with significance
  downstream %>%
    mutate(label = case_when(
      n == 0 ~ "",
      p_adj < 0.001 ~ glue("{n}\n***"),
      p_adj < 0.01 ~ glue("{n}\n**"),
      p_adj < 0.05 ~ glue("{n}\n*"),
      TRUE ~ as.character(n)
    )) %>%
    plot_relation_heatmap(
      "log2_OE_perm", "log2(O/E)",
      title = glue("Ancestor -> descendant subtype pairs vs. permutation ({relation_titles[[relation_name]]})"),
      subtitle = glue(
        "Expected: mean of {n_perm_relations} permutations of node subtype composition within each tree (GL fixed).\n",
        "Red: pair seen more often than expected. Label: N node pairs; * BH-adjusted two-sided p < 0.05, ** < 0.01, *** < 0.001"
      ),
      diverging = TRUE
    ) +
    labs(x = "Descendant subtype", y = "Ancestor subtype")
  
  ggsave(glue("{outdir_relations}/ancestor_descendant_OE_permutation_{relation_name}.png"), width = n_size + 2.5, height = n_size)
}

# ------------------------------------------------------------------------------
# 2b) Upstream / downstream subtypes per patient, with permutation test
# ------------------------------------------------------------------------------
# Same counting and node-level permutation as above, but pairs are counted and
# tested separately for each patient. BH adjustment within patient and relation.

tree_patient <- cells_all %>%
  distinct(clone_full_name, patient)

patients_rel <- sort(unique(tree_patient$patient))

count_relation_pairs_patient <- function(node_types){
  
  pairs <- node_pairs %>%
    inner_join(
      node_types %>% dplyr::rename(ancestor = node, ancestor_type = cell_type),
      by = c("clone_full_name", "ancestor"), relationship = "many-to-many"
    ) %>%
    inner_join(
      node_types %>% dplyr::rename(descendant = node, descendant_type = cell_type),
      by = c("clone_full_name", "descendant"), relationship = "many-to-many"
    ) %>%
    filter(
      !ancestor_type %in% relation_exclude_types,
      !descendant_type %in% relation_exclude_types
    ) %>%
    left_join(tree_patient, by = "clone_full_name")
  
  bind_rows(
    pairs %>% filter(step == 1) %>% mutate(relation = "direct"),
    pairs %>% mutate(relation = "any")
  ) %>%
    dplyr::count(patient, relation, ancestor_type, descendant_type, name = "n")
}

relation_grid_patient <- expand_grid(patient = patients_rel, relation_grid)

# Observed
observed_pairs_patient <- count_relation_pairs_patient(node_types) %>%
  right_join(relation_grid_patient, by = c("patient", "relation", "ancestor_type", "descendant_type")) %>%
  mutate(n = replace_na(n, 0L))

# Permutations: move each node's set of subtypes to a random observed node position of the same tree
null_pairs_patient <- map_dfr(seq_len(n_perm_relations), function(i){
  
  if (i %% 100 == 0) message(glue("Per-patient permutation {i}/{n_perm_relations}"))
  
  new_positions <- nodes_to_shuffle %>%
    group_by(clone_full_name) %>%
    mutate(new_node = if (n() > 1) sample(node) else node) %>%
    ungroup()
  
  node_types %>%
    left_join(new_positions, by = c("clone_full_name", "node")) %>%
    mutate(node = coalesce(new_node, node)) %>%
    select(-new_node) %>%
    count_relation_pairs_patient() %>%
    mutate(perm = i)
})

null_pairs_patient <- null_pairs_patient %>%
  right_join(
    expand_grid(relation_grid_patient, perm = seq_len(n_perm_relations)),
    by = c("patient", "relation", "ancestor_type", "descendant_type", "perm")
  ) %>%
  mutate(n = replace_na(n, 0L)) %>%
  dplyr::rename(n_null = n)

# Permutation p-values per patient and pair of subtypes
perm_results_patient <- null_pairs_patient %>%
  left_join(observed_pairs_patient, by = c("patient", "relation", "ancestor_type", "descendant_type")) %>%
  group_by(patient, relation, ancestor_type, descendant_type) %>%
  summarise(
    n = first(n),
    expected_perm = mean(n_null),
    p_enriched = (sum(n_null >= n) + 1) / (n_perm_relations + 1),
    p_depleted = (sum(n_null <= n) + 1) / (n_perm_relations + 1),
    .groups = "drop"
  ) %>%
  mutate(
    log2_OE_perm = log2((n + 0.5) / (expected_perm + 0.5)),
    p_two_sided = pmin(1, 2 * pmin(p_enriched, p_depleted))
  ) %>%
  group_by(patient, relation) %>%
  mutate(p_adj = p.adjust(p_two_sided, method = "BH")) %>%
  ungroup()

write_delim(perm_results_patient, glue("{outdir_relations}/pair_permutation_results_per_patient.tsv"), delim = "\t")

# Plots: one panel per patient
for (relation_name in names(relation_titles)){
  
  # relation_name <- "direct"
  
  perm_results_patient %>%
    filter(relation == relation_name) %>%
    mutate(
      row_type = factor(ancestor_type, levels = celltype_order_rel),
      col_type = factor(descendant_type, levels = celltype_order_rel),
      label = case_when(
        n == 0 ~ "",
        p_adj < 0.001 ~ glue("{n}\n***"),
        p_adj < 0.01 ~ glue("{n}\n**"),
        p_adj < 0.05 ~ glue("{n}\n*"),
        TRUE ~ as.character(n)
      )
    ) %>%
    plot_relation_heatmap(
      "log2_OE_perm", "log2(O/E)",
      title = glue("Ancestor -> descendant subtype pairs vs. permutation, per patient ({relation_titles[[relation_name]]})"),
      subtitle = glue(
        "Expected: mean of {n_perm_relations} permutations of node subtype composition within each tree (GL fixed).\n",
        "Red: pair seen more often than expected. Label: N node pairs; * BH-adjusted two-sided p < 0.05, ** < 0.01, *** < 0.001"
      ),
      diverging = TRUE
    ) +
    facet_wrap(~ patient, nrow = 1) +
    labs(x = "Descendant subtype", y = "Ancestor subtype")
  
  ggsave(
    glue("{outdir_relations}/ancestor_descendant_OE_permutation_per_patient_{relation_name}.png"),
    width = n_size * length(patients_rel) + 2,
    height = n_size
  )
}

# ==============================================================================
# 3) Clonal bursts
# ==============================================================================
#
# Definition, following the phylogenetic signature of clonal bursts described by
# Pae et al. 2025 (Nature): a large parental node of identical cells from which
# a star of descendants arises.
#   Burst     : an observed node with >= burst_min_parent_cells cells
#               AND >= burst_min_children child nodes
#   Parent    : cells on the bursting node
#   Offspring : cells on its direct child nodes
# The whole section is run once per value in burst_min_parent_cells_values,
# each in its own folder.

outdir_bursts_main <- glue("{outdir}/clonal_bursts")
dir.create(outdir_bursts_main, recursive = TRUE, showWarnings = FALSE)

# All candidate nodes (>= burst_min_children children), before the parent size threshold
burst_candidates <- bursts_all %>%
  filter(!(burst_exclude_GL & burst_node == "GL"))

write_delim(burst_candidates, glue("{outdir_bursts_main}/burst_candidates.tsv"), delim = "\t")

# ---- Choosing the threshold: parent size vs. number of children ----

burst_candidates %>%
  ggplot(aes(x = burst_parent_n_cells, y = burst_n_children)) +
  geom_jitter(width = 0.15, height = 0.15, alpha = 0.6) +
  # geom_vline(xintercept = burst_min_parent_cells_values, linetype = "dashed", color = "#D55E00") + #Dashed lines: tested parent size thresholds
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500)) +
  facet_wrap(~ patient) +
  labs(
    title = "Candidate burst nodes",
    subtitle = glue("All nodes with >= {burst_min_children} child nodes."),
    x = "N cells on the parent node. (0 = internal node)",
    y = "N child nodes"
  ) +
  theme_bw()

ggsave(glue("{outdir_bursts_main}/burst_candidates_parent_size_vs_children.png"), width = 10, height = 7)

# Number of bursts over a range of parent size thresholds
burst_threshold_overview <- map_dfr(c(1, 5, 10, 15, 20, 30, 50), function(thr){
  burst_candidates %>%
    filter(burst_parent_n_cells >= thr) %>%
    group_by(patient) %>%
    summarise(n_bursts = n(), n_trees_with_burst = n_distinct(clone_full_name), .groups = "drop") %>%
    complete(patient = unique(tree_summary_all$patient), fill = list(n_bursts = 0L, n_trees_with_burst = 0L)) %>%
    mutate(min_parent_cells = thr)
})

print(burst_threshold_overview %>% pivot_wider(id_cols = patient, names_from = min_parent_cells, values_from = n_bursts, names_prefix = ">="))
write_delim(burst_threshold_overview, glue("{outdir_bursts_main}/burst_threshold_overview.tsv"), delim = "\t")

burst_threshold_overview %>%
  ggplot(aes(x = min_parent_cells, y = n_bursts, color = patient)) +
  geom_line() +
  geom_point() +
  # geom_vline(xintercept = burst_min_parent_cells_values, linetype = "dashed", color = "grey50") +
  labs(
    title = "Number of clonal bursts by parent size threshold",
    subtitle = glue("Burst: node with >= threshold cells and >= {burst_min_children} child nodes"),
    x = "Min N cells on the parent node",
    y = "N bursts",
    color = "Patient"
  ) +
  theme_bw()

ggsave(glue("{outdir_bursts_main}/burst_threshold_overview.png"), width = 8, height = 5)

# ---- Analysis per threshold ----
n_perm_bursts <- 1000
for (min_parent_cells in burst_min_parent_cells_values){

  # min_parent_cells <- 10

  burst_definition <- glue(">= {min_parent_cells} cells on the parent node and >= {burst_min_children} child nodes")

  bursts <- burst_candidates %>%
    filter(burst_parent_n_cells >= min_parent_cells)

  if (nrow(bursts) == 0){
    message(glue("No bursts with {burst_definition}"))
    next
  }

  outdir_bursts <- glue("{outdir_bursts_main}/min_parent_cells_{min_parent_cells}")
  dir.create(outdir_bursts, recursive = TRUE, showWarnings = FALSE)

  burst_cells <- burst_cells_all %>%
    semi_join(bursts, by = c("clone_full_name", "burst_node"))

  write_delim(bursts, glue("{outdir_bursts}/bursts.tsv"), delim = "\t")
  write_delim(burst_cells, glue("{outdir_bursts}/burst_cells.tsv"), delim = "\t")

  # Overview per patient
  bursts %>%
    group_by(patient) %>%
    summarise(
      n_trees_with_burst = n_distinct(clone_full_name),
      n_bursts = n(),
      median_parent_cells = median(burst_parent_n_cells),
      median_n_children = median(burst_n_children),
      .groups = "drop"
    ) %>%
    left_join(tree_summary_all %>% dplyr::count(patient, name = "n_trees"), by = "patient") %>%
    print()

  # Burst size: parent cells and number of children
  bursts %>%
    ggplot(aes(x = burst_n_children, y = burst_parent_n_cells, size = n_cells_offspring)) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6, color = "#D55E00") +
    facet_wrap(~ patient) +
    labs(
      title = "Clonal burst sizes",
      subtitle = glue("Burst: {burst_definition}"),
      x = "N child nodes",
      y = "N cells on the parent node",
      size = "N offspring cells"
    ) +
    theme_bw()

  ggsave(glue("{outdir_bursts}/burst_size_distribution.png"), width = 9, height = 6)

  # ---- Subtypes in bursts ----
  # N bursts: number of bursts in which the subtype has >= 1 cell (per role)
  # N cells : number of cells of the subtype in that role

  burst_type_counts <- burst_cells %>%
    group_by(cell_type, burst_role) %>%
    summarise(
      n_bursts = n_distinct(clone_full_name, burst_node),
      n_trees = n_distinct(clone_full_name),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    complete(cell_type = celltype_order, burst_role = c("Parent", "Offspring"), fill = list(n_bursts = 0, n_trees = 0, n_cells = 0))

  # Fraction of each subtype's cells that are burst parent / offspring, vs. all cells.
  # Only trees with >= 1 burst, so trees without bursts do not dilute the comparison
  cells_in_burst_trees <- cells_all %>%
    filter(clone_full_name %in% bursts$clone_full_name)

  burst_enrichment <- map_dfr(c("Parent", "Offspring"), function(role){

    role_ids <- burst_cells %>%
      filter(burst_role == role) %>%
      distinct(clone_full_name, sequence_id) %>%
      mutate(in_role = TRUE)

    cells_in_burst_trees %>%
      left_join(role_ids, by = c("clone_full_name", "sequence_id")) %>%
      mutate(in_role = replace_na(in_role, FALSE)) %>%
      group_by(cell_type) %>%
      summarise(n_cells = n(), n_in_role = sum(in_role), .groups = "drop") %>%
      mutate(
        burst_role = role,
        frac_in_role = n_in_role / n_cells,
        frac_in_role_all = sum(n_in_role) / sum(n_cells),
        enrichment = frac_in_role / frac_in_role_all,
        # Fisher test: subtype vs. all other cells, in role vs. not
        p_value = map2_dbl(n_in_role, n_cells, function(k, n){
          tot_k <- sum(n_in_role); tot_n <- sum(n_cells)
          fisher.test(matrix(c(k, n - k, tot_k - k, (tot_n - n) - (tot_k - k)), nrow = 2))$p.value
        }),
        p_adj = p.adjust(p_value, method = "BH")
      )
  }) %>%
    mutate(
      cell_type = factor(cell_type, levels = celltype_order),
      burst_role = factor(burst_role, levels = c("Parent", "Offspring"))
    ) %>%
    arrange(burst_role, cell_type)

  print(burst_enrichment)
  write_delim(burst_type_counts, glue("{outdir_bursts}/burst_celltype_counts.tsv"), delim = "\t")
  write_delim(burst_enrichment, glue("{outdir_bursts}/burst_enrichment.tsv"), delim = "\t")

  # Number of bursts each subtype takes part in
  burst_type_counts %>%
    mutate(
      cell_type = factor(cell_type, levels = celltype_order),
      burst_role = factor(burst_role, levels = c("Parent", "Offspring"))
    ) %>%
    ggplot(aes(x = cell_type, y = n_bursts, fill = burst_role)) +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.7) +
    geom_text(aes(label = n_bursts), position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Parent" = "#0072B2", "Offspring" = "#E69F00")) +
    labs(
      title = "Clonal bursts per subtype",
      subtitle = glue("N bursts with >= 1 cell of the subtype on the parent node or its child nodes (offspring).\nBurst: {burst_definition}. Total bursts: {nrow(bursts)}"),
      x = NULL,
      y = "N bursts",
      fill = "Role in burst"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(glue("{outdir_bursts}/burst_celltype_counts.png"), width = max(8, 0.9 * length(celltype_order)), height = 6)

  # Composition of each burst: parent and offspring subtypes per burst (largest bursts)
  burst_cells %>%
    mutate(burst_id = glue("{clone_full_name} | {burst_node}")) %>%
    add_count(burst_id, name = "n_burst_cells") %>%
    mutate(
      burst_id = fct_reorder(burst_id, n_burst_cells, .desc = TRUE),
      cell_type = factor(cell_type, levels = celltype_order),
      burst_role = factor(burst_role, levels = c("Parent", "Offspring"))
    ) %>%
    filter(as.integer(burst_id) <= 40) %>%
    ggplot(aes(x = burst_id, fill = cell_type)) +
    geom_bar() +
    facet_wrap(~ burst_role, ncol = 1, scales = "free_y") +
    labs(
      title = "Subtype composition of the 40 largest bursts",
      subtitle = glue("Burst: {burst_definition}"),
      x = "Burst (clone | node)",
      y = "N cells",
      fill = "Subtype"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

  ggsave(glue("{outdir_bursts}/burst_composition.png"), width = 14, height = 9)
  
  # ---- Permutation test ----
  # Subtype labels are shuffled among the cells of each tree with >= 1 burst
  # (bursts and their cells fixed). Statistic: N cells of each subtype that are
  # burst parent / offspring. Computed per patient, and pooled over patients.
  
  role_cells <- cells_in_burst_trees %>%
    select(patient, clone_full_name, sequence_id, cell_type) %>%
    left_join(
      burst_cells %>% filter(burst_role == "Parent") %>% distinct(clone_full_name, sequence_id) %>% mutate(is_parent = TRUE),
      by = c("clone_full_name", "sequence_id")
    ) %>%
    left_join(
      burst_cells %>% filter(burst_role == "Offspring") %>% distinct(clone_full_name, sequence_id) %>% mutate(is_offspring = TRUE),
      by = c("clone_full_name", "sequence_id")
    ) %>%
    mutate(across(c(is_parent, is_offspring), ~ replace_na(.x, FALSE)))
  
  count_roles <- function(df){
    df %>%
      group_by(patient, cell_type) %>%
      summarise(
        n_cells = n(),
        Parent = sum(is_parent),
        Offspring = sum(is_offspring),
        .groups = "drop"
      ) %>%
      pivot_longer(c(Parent, Offspring), names_to = "burst_role", values_to = "n_in_role")
  }
  
  observed_roles <- count_roles(role_cells)
  
  null_roles <- map_dfr(seq_len(n_perm_bursts), function(i){
    role_cells %>%
      group_by(clone_full_name) %>%
      mutate(cell_type = sample(cell_type)) %>%
      ungroup() %>%
      count_roles() %>%
      mutate(perm = i)
  })
  
  summarise_perm <- function(observed, null, ...){
    null %>%
      dplyr::rename(n_null = n_in_role) %>%
      select(..., burst_role, perm, n_null) %>%
      left_join(observed, by = c(sapply(rlang::ensyms(...), rlang::as_string), "burst_role")) %>%
      group_by(..., burst_role) %>%
      summarise(
        n_cells = first(n_cells),
        n_in_role = first(n_in_role),
        expected_perm = mean(n_null),
        p_enriched = (sum(n_null >= n_in_role) + 1) / (n_perm_bursts + 1),
        p_depleted = (sum(n_null <= n_in_role) + 1) / (n_perm_bursts + 1),
        .groups = "drop"
      ) %>%
      mutate(
        frac_in_role = n_in_role / n_cells,
        frac_expected = expected_perm / n_cells,
        log2_OE_perm = log2((n_in_role + 0.5) / (expected_perm + 0.5)),
        p_two_sided = pmin(1, 2 * pmin(p_enriched, p_depleted))
      )
  }
  
  # Per patient
  burst_perm_patient <- summarise_perm(observed_roles, null_roles, patient, cell_type) %>%
    group_by(patient, burst_role) %>%
    mutate(p_adj = p.adjust(p_two_sided, method = "BH")) %>%
    ungroup()
  
  # Pooled over patients (sum of the per-patient counts, per permutation)
  observed_pooled <- observed_roles %>%
    group_by(cell_type, burst_role) %>%
    summarise(n_cells = sum(n_cells), n_in_role = sum(n_in_role), .groups = "drop")
  
  null_pooled <- null_roles %>%
    group_by(cell_type, burst_role, perm) %>%
    summarise(n_in_role = sum(n_in_role), .groups = "drop")
  
  burst_perm <- summarise_perm(observed_pooled, null_pooled, cell_type) %>%
    group_by(burst_role) %>%
    mutate(p_adj = p.adjust(p_two_sided, method = "BH")) %>%
    ungroup()
  
  write_delim(burst_perm, glue("{outdir_bursts}/burst_enrichment_permutation.tsv"), delim = "\t")
  write_delim(burst_perm_patient, glue("{outdir_bursts}/burst_enrichment_permutation_per_patient.tsv"), delim = "\t")
  
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  
  # Trees and bursts per patient, for captions and facet labels
  burst_tree_counts <- tree_summary_all %>%
    dplyr::count(patient, name = "n_trees") %>%
    left_join(
      bursts %>%
        group_by(patient) %>%
        summarise(n_trees_burst = n_distinct(clone_full_name), n_bursts = n(), .groups = "drop"),
      by = "patient"
    ) %>%
    mutate(
      across(c(n_trees_burst, n_bursts), ~ replace_na(.x, 0L)),
      label_short = glue("{patient}: {n_trees_burst}/{n_trees} trees with bursts ({n_bursts} bursts)"),
      facet_label = glue("{patient}\n{n_trees_burst}/{n_trees} trees with bursts, {n_bursts} bursts")
    )
  
  caption_trees <- burst_tree_counts %>%
    filter(n_trees_burst > 0) %>%
    pull(label_short) %>%
    paste(collapse = "\n")
  patient_labeller <- setNames(as.character(burst_tree_counts$facet_label), burst_tree_counts$patient)
  
  # Enrichment of subtypes as burst parent and offspring: observed vs. permutation, pooled
  burst_perm %>%
    mutate(
      cell_type = factor(cell_type, levels = celltype_order),
      burst_role = factor(burst_role, levels = c("Parent", "Offspring")),
      label = glue("{n_in_role}/{n_cells}{stars(p_adj)}")
    ) %>%
    ggplot(aes(x = cell_type)) +
    geom_col(aes(y = frac_in_role, fill = burst_role), width = 0.7) +
    geom_point(aes(y = frac_expected), shape = 23, size = 2.5, fill = "white", color = "black") +
    geom_text(aes(y = pmax(frac_in_role, frac_expected), label = label), vjust = -0.6, size = 3) +
    facet_wrap(~ burst_role, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("Parent" = "#0072B2", "Offspring" = "#E69F00"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title = "Fraction of cells that are burst parent / offspring, per subtype",
      subtitle = glue(
        "Burst: {burst_definition}. Trees with >= 1 burst. Diamond: expected from {n_perm_bursts} permutations of subtype labels within each tree.\n",
        "Label: N cells in role / N cells (observed); * BH-adjusted two-sided permutation p < 0.05, ** < 0.01, *** < 0.001"
      ),
      caption = glue("Trees with clonal bursts per patient:\n{caption_trees}"),
      x = NULL,
      y = "Fraction of cells"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(glue("{outdir_bursts}/burst_enrichment.png"), width = max(9, 0.9 * length(celltype_order)), height = 9)
  
  # Same, per patient
  burst_perm_patient %>%
    mutate(
      cell_type = factor(cell_type, levels = celltype_order),
      burst_role = factor(burst_role, levels = c("Parent", "Offspring")),
      label = glue("{n_in_role}/{n_cells}{stars(p_adj)}")
    ) %>%
    ggplot(aes(x = cell_type)) +
    geom_col(aes(y = frac_in_role, fill = burst_role), width = 0.7) +
    geom_point(aes(y = frac_expected), shape = 23, size = 2, fill = "white", color = "black") +
    geom_text(aes(y = pmax(frac_in_role, frac_expected), label = label), vjust = -0.6, size = 2.3) +
    # facet_grid(rows = vars(burst_role), cols = vars(patient), scales = "free_y") +
    facet_grid(rows = vars(burst_role), cols = vars(patient), scales = "free_y",
               labeller = labeller(patient = patient_labeller)) +
    scale_fill_manual(values = c("Parent" = "#0072B2", "Offspring" = "#E69F00"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = "Fraction of cells that are burst parent / offspring, per subtype and patient",
      subtitle = glue(
        "Burst: {burst_definition}. Diamond: expected from {n_perm_bursts} permutations of subtype labels within each tree.\n",
        "Label: N cells in role / N cells (observed); * BH-adjusted two-sided permutation p < 0.05, ** < 0.01, *** < 0.001"
      ),
      caption = glue("Trees with clonal bursts per patient:\n{caption_trees}"),
      x = NULL,
      y = "Fraction of cells"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    glue("{outdir_bursts}/burst_enrichment_per_patient.png"),
    width = max(12, 0.6 * length(celltype_order) * n_distinct(burst_perm_patient$patient)),
    height = 8
  )
  
}
