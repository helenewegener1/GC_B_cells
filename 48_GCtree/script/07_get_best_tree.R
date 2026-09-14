library(tidyverse)
library(glue)
library(ape)

# ------------------------------------------------------------------------------
# Purpose
# ------------------------------------------------------------------------------
# GCtree can return several equally (or near-equally) parsimonious trees per
# clone (plot/{clone}/{clone}.inference.1.nk, .2.nk, ...). This script uses
# isotype switching biology to narrow down which of those candidate trees is
# most plausible: class-switch recombination (CSR) is a one-directional,
# stepwise process, so a tree in which a descendant's isotype is *upstream*
# of its ancestor's isotype implies a biologically impossible reverse switch.
# For each clone we score every candidate tree by counting these "reverse"
# violations and pick the best one (fewest violations, ties broken by
# GCtree's own likelihood rank).
#
# This only needs files that already exist locally under 48_GCtree/:
#   - gctree_meta/{clone}_gctree_meta.txt  (seq_unique -> isotype mapping,
#     built once by 03_prep_metadata.R from out/{clone}/idmap.txt)
#   - plot/{clone}/{clone}.inference.*.nk  (GCtree's candidate trees)
# It does NOT need out/{clone}/idmap.txt itself, or the big
# 45_immcantation/.../06_resolve_LC_germlined.rds -- both of which only
# exist on the HPC where GCtree was actually run, not in the local project
# folder. Clones are discovered directly from the gctree_meta files present,
# so nothing needs to be hardcoded about which patients/clones exist.

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

meta_dir <- "48_GCtree/gctree_meta"
plot_dir <- "48_GCtree/plot"
outdir <- "48_GCtree/out/best_tree"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Discover clones from the local gctree_meta files
# ------------------------------------------------------------------------------

meta_files <- list.files(meta_dir, pattern = "_gctree_meta\\.txt$")

clones <- tibble(meta_file = meta_files) %>%
  mutate(
    clone_full_name = str_remove(meta_file, "_gctree_meta\\.txt$"),
    HH = str_extract(clone_full_name, "^HH\\d+"),
    clone_nr = str_match(clone_full_name, "_clone_nr_(\\d+)_clone_")[, 2] %>% as.integer(),
    clone = str_match(clone_full_name, "_clone_nr_\\d+_clone_(.+)$")[, 2]
  )

# ------------------------------------------------------------------------------
# Isotype switch order (human IgH class-switch recombination order)
# Non-switched (M/D) -> G3 -> G1 -> A1 -> G2 -> G4 -> E -> A2
# A descendant isotype that ranks *before* its ancestor's isotype in this
# vector is a "reverse" switch and is not biologically plausible.
# ------------------------------------------------------------------------------

isotype_switch_order <- c("IGHM/D", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")

# ------------------------------------------------------------------------------
# Helper functions (same logic as 06_analyse_trees.R)
# ------------------------------------------------------------------------------

# For a given observed node, walk up the tree until hitting another observed
# ("seq*") node. Internal nodes with plain numeric labels ("1", "2", ...) are
# unobserved ancestral genotypes and are skipped over.
get_nearest_observed_ancestor <- function(node, parent_of) {
  current <- parent_of[node]
  while (!is.na(current) && !str_detect(current, "^seq[0-9]+$")) {
    current <- parent_of[current]
  }
  if (is.na(current)) NA_character_ else current
}

# Resolves each observed node's isotype top-down: a mixed node (a collapsed
# GCtree node that maps back to sequences of more than one isotype) prefers
# continuity with its resolved ancestor's isotype (assumes no switch) over
# blindly picking the earliest isotype in the switch order.
resolve_node_isotypes <- function(observed_nodes, ancestor_map, isotype_lookup, isotype_switch_order) {

  resolved <- list()
  remaining <- observed_nodes

  while (length(remaining) > 0) {

    ready <- remaining[map_lgl(remaining, ~ {
      anc <- ancestor_map[[.x]]
      is.na(anc) || anc %in% names(resolved)
    })]

    if (length(ready) == 0) break

    for (node in ready) {

      node_isotypes <- isotype_lookup[[node]]

      if (length(node_isotypes) == 0) {
        resolved[[node]] <- NA_character_
      } else if (length(node_isotypes) == 1) {
        resolved[[node]] <- node_isotypes[[1]]
      } else {

        anc <- ancestor_map[[node]]
        anc_isotype <- if (!is.na(anc) && anc %in% names(resolved)) resolved[[anc]] else NA_character_

        if (!is.na(anc_isotype) && anc_isotype %in% node_isotypes) {
          resolved[[node]] <- anc_isotype
        } else {
          resolved[[node]] <- node_isotypes[[which.min(match(node_isotypes, isotype_switch_order))]]
        }
      }
    }

    remaining <- setdiff(remaining, ready)
  }

  resolved
}

# ------------------------------------------------------------------------------
# Build a seq_unique -> isotype(s) lookup for one clone straight from its
# already-generated gctree_meta text file (no idmap.txt / germline RDS needed)
# ------------------------------------------------------------------------------

build_isotype_lookup <- function(clone_full_name, meta_dir) {

  meta <- read_csv(glue("{meta_dir}/{clone_full_name}_gctree_meta.txt"), show_col_types = FALSE)

  # c_call_grouped is a ":"-joined list of the unique isotypes underlying a
  # (possibly collapsed) seq_unique node, e.g. "IGHA1:IGHA2", or "NA" when no
  # isotype was known. readr converts a cell that is exactly "NA" to a real
  # NA on read, so replace_na("") normalises that back to an empty string
  # before splitting -- either way, "NA" tokens are dropped, same as the
  # na.omit() used when this mapping was built the first time.
  meta <- meta %>%
    mutate(
      isotypes = c_call_grouped %>%
        replace_na("") %>%
        str_split(":") %>%
        map(~ .x[.x != "" & .x != "NA"] %>% unique())
    )

  set_names(meta$isotypes, meta$seq_unique)
}

# ------------------------------------------------------------------------------
# Score a single candidate tree for isotype-switch-order violations
# ------------------------------------------------------------------------------

score_tree_isotype_switching <- function(tree_file, isotype_lookup, isotype_switch_order) {

  tree <- read.tree(tree_file)

  node_labels <- c(tree$tip.label, tree$node.label)
  edge_df <- as_tibble(tree$edge) %>%
    set_names(c("parent_idx", "child_idx")) %>%
    mutate(parent = node_labels[parent_idx], child = node_labels[child_idx])
  parent_of <- set_names(edge_df$parent, edge_df$child)

  observed_nodes <- node_labels[str_detect(node_labels, "^seq[0-9]+$")]

  ancestor_map <- set_names(
    map_chr(observed_nodes, ~ get_nearest_observed_ancestor(.x, parent_of)),
    observed_nodes
  )

  resolved_isotype <- resolve_node_isotypes(observed_nodes, ancestor_map, isotype_lookup, isotype_switch_order)
  resolved_isotype_vec <- unlist(resolved_isotype)

  switch_df <- tibble(node = observed_nodes) %>%
    mutate(
      ancestor = ancestor_map[node],
      node_isotypes = isotype_lookup[node],
      ancestor_isotypes = isotype_lookup[ancestor],
      node_isotype = resolved_isotype_vec[node],
      ancestor_isotype = resolved_isotype_vec[ancestor]
    ) %>%
    filter(!is.na(ancestor), !is.na(node_isotype), !is.na(ancestor_isotype)) %>%
    rowwise() %>%
    mutate(
      ancestor_rank = match(ancestor_isotype, isotype_switch_order),
      node_rank = match(node_isotype, isotype_switch_order),
      switch_type = case_when(
        node_isotype %in% ancestor_isotypes ~ "none",   # descendant's isotype was already present at the ancestor -> no real switch
        node_rank > ancestor_rank ~ "sequential",
        node_rank < ancestor_rank ~ "reverse",
        TRUE ~ "unknown"
      )
    ) %>%
    ungroup()

  list(
    switch_df = switch_df,
    summary = tibble(
      n_edges_scored = nrow(switch_df),
      n_none = sum(switch_df$switch_type == "none"),
      n_sequential = sum(switch_df$switch_type == "sequential"),
      n_reverse = sum(switch_df$switch_type == "reverse")
    )
  )
}

# ------------------------------------------------------------------------------
# For one clone: score every candidate tree GCtree produced and pick the best
# ------------------------------------------------------------------------------

get_best_tree <- function(HH, clone_nr, clone, clone_full_name, isotype_switch_order) {

  # HH <- "HH117"; clone_nr <- 1; clone <- "4252_1"
  # clone_full_name <- glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")

  isotype_lookup <- build_isotype_lookup(clone_full_name, meta_dir)

  # Discover every ranked tree GCtree produced for this clone (not just #1)
  tree_files <- Sys.glob(glue("{plot_dir}/{clone_full_name}/{clone_full_name}.inference.*.nk"))

  if (length(tree_files) == 0) {
    warning(glue("No GCtree trees found for {clone_full_name}"))
    return(NULL)
  }

  # gctree ranks trees by likelihood: rank 1 = most likely branching/abundance assignment
  gctree_rank <- as.integer(str_extract(tree_files, "(?<=\\.inference\\.)\\d+(?=\\.nk)"))
  ord <- order(gctree_rank)
  tree_files <- tree_files[ord]
  gctree_rank <- gctree_rank[ord]

  clone_scores <- tibble()

  for (i in seq_along(tree_files)) {

    scored <- score_tree_isotype_switching(tree_files[[i]], isotype_lookup, isotype_switch_order)

    clone_scores <- bind_rows(
      clone_scores,
      scored$summary %>%
        mutate(
          HH = HH, clone_nr = clone_nr, clone = clone,
          gctree_rank = gctree_rank[[i]], tree_file = tree_files[[i]],
          .before = 1
        )
    )
  }

  # A tree "passes" if it has zero reverse (biologically implausible) isotype switches.
  # Best tree = fewest reverse switches; ties broken by GCtree's own rank (best likelihood first).
  clone_scores <- clone_scores %>%
    mutate(passes_isotype_filter = n_reverse == 0) %>%
    arrange(n_reverse, gctree_rank)

  best_tree_file <- clone_scores$tree_file[[1]]
  clone_scores <- clone_scores %>% mutate(chosen = tree_file == best_tree_file)

  clone_scores

}

# ------------------------------------------------------------------------------
# Run across all discovered clones
# ------------------------------------------------------------------------------

all_tree_scores <- tibble()

for (i in seq_len(nrow(clones))) {

  row <- clones[i, ]

  clone_scores <- tryCatch(
    get_best_tree(row$HH, row$clone_nr, row$clone, row$clone_full_name, isotype_switch_order),
    error = function(e) {
      message(glue("Skipping {row$clone_full_name}: {conditionMessage(e)}"))
      NULL
    }
  )

  all_tree_scores <- bind_rows(all_tree_scores, clone_scores)

}

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

# Full table: every candidate tree for every clone, with its violation counts
write_delim(all_tree_scores %>% select(-tree_file), delim = "\t", file = glue("{outdir}/all_tree_isotype_scores.tsv"))
saveRDS(all_tree_scores, glue("{outdir}/all_tree_isotype_scores.rds"))

# One row per clone: the chosen "best" tree and its file path
best_trees <- all_tree_scores %>% filter(chosen)
write_delim(best_trees %>% select(-tree_file), delim = "\t", file = glue("{outdir}/best_trees.tsv"))
saveRDS(best_trees, glue("{outdir}/best_trees.rds"))

# ------------------------------------------------------------------------------
# Quick overview
# ------------------------------------------------------------------------------

# How many clones actually had >1 candidate tree, and for how many of those
# did the isotype filter actually discriminate (i.e. not all trees passed)?
overview <- all_tree_scores %>%
  group_by(HH, clone_nr, clone) %>%
  summarise(n_trees = n(), n_passing = sum(passes_isotype_filter), .groups = "drop") %>%
  filter(n_trees > 1)

overview %>% count(discriminating = n_passing < n_trees)

# Clones where GCtree's own top-ranked tree (#1) was NOT the isotype-best tree
# -- i.e. cases where the biological filter actually overrode GCtree's ranking
all_tree_scores %>%
  filter(chosen, gctree_rank != 1) %>%
  select(HH, clone_nr, clone, gctree_rank, n_reverse, n_sequential)

# ------------------------------------------------------------------------------
# Diagnostic plot: reverse isotype switches per candidate tree, per clone
# ------------------------------------------------------------------------------

plot_data <- all_tree_scores %>%
  semi_join(overview, by = c("HH", "clone_nr", "clone")) %>%  # only clones with >1 candidate tree
  mutate(clone_plot = glue("{HH} #{clone_nr}\n({clone})"))

if (nrow(plot_data) > 0) {

  ggplot(plot_data, aes(x = factor(gctree_rank), y = n_reverse, fill = chosen)) +
    geom_col() +
    facet_wrap(~clone_plot, scales = "free_x") +
    scale_fill_manual(values = c(`TRUE` = "forestgreen", `FALSE` = "grey70"), name = "Chosen as\nbest tree") +
    labs(
      title = "Isotype-switch-order violations across GCtree's candidate trees",
      x = "GCtree rank (1 = most likely)",
      y = "N reverse (biologically implausible) isotype switches"
    ) +
    theme_bw()

  ggsave(glue("{outdir}/candidate_trees_reverse_switches.png"), width = 12, height = 8)
}
