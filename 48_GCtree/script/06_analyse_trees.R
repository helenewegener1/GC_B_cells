library(tidyverse)
library(glue)
library(ape)
library(ggnewscale)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/06_resolve_LC_germlined.rds")
seq_dir <- readRDS("48_GCtree/out/rds/seq_dir.rds")

patients <- names(resolve_LC_list_germlined)

HH <- "HH119" 

n_clones_begin_list <- list(
  "HH117" = 1, 
  "HH119" = 2
)

n_clones_begin <- n_clones_begin_list[[HH]]

outdir <- glue("48_GCtree/plot_06_analyse_trees/{HH}")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Get top clones
top_clones <- resolve_LC_list_germlined[[HH]] %>% 
  filter(
    locus == "IGH", 
    L1_annotation == "GC_B_cells"
  ) %>% 
  dplyr::count(clone_subgroup_id_90_similarity, sort = TRUE) 


# ==============================================================================
# Using GCtree nk files
# ==============================================================================

# ------------------------------------------------------------------------------
# Tree-based (patristic) distance matrix from GCtree output, aligned to df_clone row order
# ------------------------------------------------------------------------------

clone_nr <- 3
# clone_nr <- 8

# Get clone
clone <- top_clones %>% 
  dplyr::slice(clone_nr) %>% 
  pull(clone_subgroup_id_90_similarity)

# Get clone files 
clone_full_name <- glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")

sequence_ids <- seq_dir[[clone_full_name]]

# filter for clone 
df_clone <- resolve_LC_list_germlined[[HH]] %>% 
  filter(
    locus == "IGH", 
    clone_subgroup_id_90_similarity == clone
  ) %>% 
  mutate(
    sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)]
  )

df_clone$sequence_alignment == sequence_ids

tree_file <- glue("48_GCtree/plot/{clone_full_name}/{clone_full_name}.inference.1.nk")
idmap_file <- glue("48_GCtree/out/{clone_full_name}/idmap.txt")

# Load gctree
tree <- read.tree(tree_file)

# Distance between ALL nodes (tips + internal) -- some observed genotypes
# are internal nodes in a collapsed GCtree, not just tips
node_labels <- c(tree$tip.label, tree$node.label)
node_dist <- dist.nodes(tree)
rownames(node_dist) <- colnames(node_dist) <- node_labels

# Keep only observed-sequence nodes (drop GL root + unobserved ancestral nodes,
# which have plain numeric labels like "1", "2", "3"...)
seq_nodes <- node_labels[str_detect(node_labels, "^seq[0-9]+$")]
seq_dist  <- node_dist[seq_nodes, seq_nodes]

# Expand collapsed (identical-sequence) nodes back to individual sequence_id's
idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id"), show_col_types = FALSE) %>% 
  filter(node != "GL") %>% 
  separate_longer_delim(sequence_id, delim = ":")

dist_df <- idmap %>% 
  dplyr::rename(node_i = node, sequence_id_i = sequence_id) %>% 
  cross_join(idmap %>% dplyr::rename(node_j = node, sequence_id_j = sequence_id)) %>% 
  mutate(distance = seq_dist[cbind(node_i, node_j)])

dist_mat <- dist_df %>% 
  select(sequence_id_i, sequence_id_j, distance) %>% 
  pivot_wider(names_from = sequence_id_j, values_from = distance) %>% 
  column_to_rownames("sequence_id_i") %>% 
  as.matrix()

# Check seq ids 
length(sequence_ids) == length(rownames(dist_mat))
stopifnot(all(names(sequence_ids) %in% rownames(dist_mat)))
# dist_mat[sequence_ids, sequence_ids]
dist_mat[rownames(dist_mat), rownames(dist_mat)]

# Get average sequence similarity
# dist_mean <- mean(dist_mat)
mean(dist_mat)
median(dist_mat)
dist_median <- median(dist_mat)
dist_median

# Get per follicle distance
# Initiate df_plot
df_plot <- data.frame()

# Define follicles of clone
dist_mat %>% dim()
clone_follicles <- df_clone$sample_clean_fol %>% unique() %>% str_subset("Fol")
all_fols_index <- grep("Fol", df_clone$sample_clean_fol)
all_fols_seq_ids <- df_clone$sequence_id_gctree[all_fols_index]

for (i in 1:length(clone_follicles)){

  # i <- 1

  # N cells in follicle
  n_cells <- sum(df_clone$sample_clean_fol == clone_follicles[[i]])

  # Get index for follicle for distance matrix
  clone_follicle_index <- grep(clone_follicles[[i]], df_clone$sample_clean_fol)
  
  # Get sequence ids 
  clone_follicle_seq_ids <- df_clone$sequence_id_gctree[clone_follicle_index]
  other_fols_seq_ids <- setdiff(all_fols_seq_ids, clone_follicle_seq_ids)

  # distance from follicle cells to follicle cells
  intra_fol_distance <- dist_mat[clone_follicle_seq_ids, clone_follicle_seq_ids] %>% median()

  # distance from follicle cells to all other cells
  inter_fol_distance <- dist_mat[clone_follicle_seq_ids, other_fols_seq_ids] %>% median()

  # Prep row for currect follicle
  df_plot_fol <- data.frame(
    clone = clone,
    dist_median = dist_median,
    follicle = clone_follicles[[i]],
    n_cells = n_cells,
    intra_fol_distance = intra_fol_distance,
    inter_fol_distance = inter_fol_distance
  )

  # Append to dataframe
  df_plot <- bind_rows(df_plot, df_plot_fol)

}

df_plot %>%
  mutate(
    follicle_nr = str_split_i(follicle, "-", -1) %>% as.integer(),
    follicle_plot = glue("{follicle_nr}\n({n_cells})") %>% fct_reorder(follicle_nr)
  ) %>%
  pivot_longer(
    cols = contains("dist"),
    names_to = "distance_type",
    values_to = "distance"
  ) %>%
  ggplot(aes(x = follicle_plot, y = distance, color = distance_type)) +
  geom_point() +
  geom_hline(aes(yintercept = dist_median), color = "grey") +
  # scale_color_manual(values = c())
  scale_color_discrete(labels = c(
    dist_median = "Median clone",
    inter_fol_distance = "Median inter-follicular",
    intra_fol_distance = "Median intra-follicular"
  )) +
  labs(
    title = glue("{HH}: Full sequence, clone nr {clone_nr}, {clone}"),
    y = "GCtree patristic distance (# mutations)",
    x = "Follicle\n(N cells)",
    color = "Distance measure"
  ) +
  theme_bw()

ggsave(glue("{outdir}/{HH}_clone_nr_{clone_nr}_clone_{clone}_distance_median_gctree.png"), width = 10, height = 6)


# ------------------------------------------------------------------------------
# Same as ^^ but with jitter
# ------------------------------------------------------------------------------

df_plot <- data.frame()
clone_follicles <- df_clone$sample_clean_fol %>% unique() %>% str_subset("Fol")
all_fols_index <- grep("Fol", df_clone$sample_clean_fol)
all_fols_seq_ids <- df_clone$sequence_id_gctree[all_fols_index]

for (i in seq_along(clone_follicles)){
  
  n_cells <- sum(df_clone$sample_clean_fol == clone_follicles[[i]])
  clone_follicle_index <- grep(clone_follicles[[i]], df_clone$sample_clean_fol)
  clone_follicle_seq_ids <- df_clone$sequence_id_gctree[clone_follicle_index]
  other_fols_seq_ids <- setdiff(all_fols_seq_ids, clone_follicle_seq_ids)
  
  # individual intra-follicular distances (as.dist drops diagonal + duplicate symmetric pairs)
  intra_fol_distance <- dist_mat[clone_follicle_seq_ids, clone_follicle_seq_ids] %>% as.dist() %>% as.vector()
  
  # individual inter-follicular distances
  inter_fol_distance <- dist_mat[clone_follicle_seq_ids, other_fols_seq_ids] %>% as.vector()
  
  df_plot_fol <- bind_rows(
    tibble(clone = clone, dist_median = dist_median, follicle = clone_follicles[[i]],
           n_cells = n_cells, distance_type = "intra_fol_distance", distance = intra_fol_distance),
    tibble(clone = clone, dist_median = dist_median, follicle = clone_follicles[[i]],
           n_cells = n_cells, distance_type = "inter_fol_distance", distance = inter_fol_distance)
  )
  
  df_plot <- bind_rows(df_plot, df_plot_fol)
}

df_plot %>%
  mutate(
    follicle_nr = str_split_i(follicle, "-", -1) %>% as.integer(),
    follicle_plot = glue("{follicle_nr}\n({n_cells})") %>% fct_reorder(follicle_nr)
  ) %>%
  ggplot(aes(x = follicle_plot, y = distance, color = distance_type)) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_summary(
    aes(fill = distance_type),
    fun = median, geom = "point",
    shape = 23, size = 3, color = "black", stroke = 0.8
  ) +
  scale_fill_discrete(guide = "none") + 
  geom_hline(aes(yintercept = dist_median), color = "grey") +
  scale_color_discrete(labels = c(
    inter_fol_distance = "Inter-follicular",
    intra_fol_distance = "Intra-follicular"
  )) +
  labs(
    title = glue("{HH}: Full sequence, clone nr {clone_nr}, {clone}"),
    y = "GCtree patristic distance (# mutations)",
    x = "Follicle\n(N cells)",
    color = "Distance measure"
  ) +
  theme_bw()

ggsave(glue("{outdir}/{HH}_clone_nr_{clone_nr}_clone_{clone}_distance_jitter_gctree.png"), width = 10, height = 6)

# ------------------------------------------------------------------------------
# Sequence similarity - median boxplots across top 10 clones - GCtree distance
# ------------------------------------------------------------------------------

# Initiate df_plot 
df_plot <- data.frame()
n_clones <- 20

for (HH in patients){
  
  # HH <- "HH117"
  # resolve_LC_list_germlined

  # Get top clones
  these_clones <- resolve_LC_list_germlined[[HH]] %>% 
    filter(
      locus == "IGH", 
      L1_annotation == "GC_B_cells"
    ) %>% 
    dplyr::count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
    head(n_clones) %>% 
    pull(clone_subgroup_id_90_similarity)
  
  # for (clone_nr in 1:n_clones){
  for (clone_nr in n_clones_begin:n_clones){
    
    # clone_nr <- 3
    
    # Get clone
    clone <- these_clones[clone_nr]
    
    # Get clone files 
    clone_full_name <- glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")
    
    sequence_ids <- seq_dir[[clone_full_name]]
    # table(df_clone$sequence_alignment == sequence_ids)
    
    # filter for clone 
    df_clone <- resolve_LC_list_germlined[[HH]] %>% 
      filter(
        locus == "IGH", 
        clone_subgroup_id_90_similarity == clone
      ) %>% 
      mutate(
        sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)]
      )
    
    # Get clone files 
    tree_file <- glue("48_GCtree/plot/{clone_full_name}/{clone_full_name}.inference.1.nk")
    idmap_file <- glue("48_GCtree/out/{clone_full_name}/idmap.txt")
    
    # Load gctree
    tree <- read.tree(tree_file)
    
    # Distance between ALL nodes (tips + internal) -- some observed genotypes
    # are internal nodes in a collapsed GCtree, not just tips
    node_labels <- c(tree$tip.label, tree$node.label)
    node_dist <- dist.nodes(tree)
    rownames(node_dist) <- colnames(node_dist) <- node_labels
    
    # Keep only observed-sequence nodes (drop GL root + unobserved ancestral nodes,
    # which have plain numeric labels like "1", "2", "3"...)
    seq_nodes <- node_labels[str_detect(node_labels, "^seq[0-9]+$")]
    seq_dist  <- node_dist[seq_nodes, seq_nodes]
    
    # Expand collapsed (identical-sequence) nodes back to individual sequence_id's
    idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id"), show_col_types = FALSE) %>% 
      filter(node != "GL") %>% 
      separate_longer_delim(sequence_id, delim = ":")
    
    dist_df <- idmap %>% 
      dplyr::rename(node_i = node, sequence_id_i = sequence_id) %>% 
      cross_join(idmap %>% dplyr::rename(node_j = node, sequence_id_j = sequence_id)) %>% 
      mutate(distance = seq_dist[cbind(node_i, node_j)])
    
    dist_mat <- dist_df %>% 
      select(sequence_id_i, sequence_id_j, distance) %>% 
      pivot_wider(names_from = sequence_id_j, values_from = distance) %>% 
      column_to_rownames("sequence_id_i") %>% 
      as.matrix()
    
    # stopifnot(all(sequence_ids %in% rownames(dist_mat)))
    # dist_mat[sequence_ids, sequence_ids]
    dist_mat[rownames(dist_mat), rownames(dist_mat)]
    
    # Get average sequence similarity
    dist_median <- median(dist_mat)
    dist_median
    
    # Get per follicle distance
    
    # Define follicles of clone
    dist_mat %>% dim()
    clone_follicles <- df_clone$sample_clean_fol %>% unique() %>% str_subset("Fol")
    all_fols_index <- grep("Fol", df_clone$sample_clean_fol)
    all_fols_seq_ids <- df_clone$sequence_id_gctree[all_fols_index]
    
    if (length(clone_follicles) == 0){
      next
    }
    
    for (i in 1:length(clone_follicles)){
      
      # i <- 1
      
      # N cells in follicle
      n_cells <- sum(df_clone$sample_clean_fol == clone_follicles[[i]])
      
      # Get index for follicle for distance matrix
      clone_follicle_index <- grep(clone_follicles[[i]], df_clone$sample_clean_fol)
      
      # Get sequence ids 
      clone_follicle_seq_ids <- df_clone$sequence_id_gctree[clone_follicle_index] %>% na.omit()
      other_fols_seq_ids <- setdiff(all_fols_seq_ids, clone_follicle_seq_ids) %>% na.omit()
      
      # distance from follicle cells to follicle cells
      intra_fol_distance <- dist_mat[clone_follicle_seq_ids, clone_follicle_seq_ids] %>% median()
      
      # distance from follicle cells to all other cells
      inter_fol_distance <- dist_mat[clone_follicle_seq_ids, other_fols_seq_ids] %>% median()
      
      # Prep row for currect follicle 
      df_plot_fol <- data.frame(
        patient = HH, 
        clone = clone,
        clone_nr = clone_nr,
        dist_median = dist_median, 
        follicle = clone_follicles[[i]],
        n_cells = n_cells, 
        intra_fol_distance = intra_fol_distance,
        inter_fol_distance = inter_fol_distance
      ) 
      
      # Append to dataframe 
      df_plot <- bind_rows(df_plot, df_plot_fol)
      
    }
    
  }
  
}

# Statistical test: intra- vs. inter-follicular distance, within each patient
# Paired Wilcoxon test
# Pair because the inter- and intra-follicular distances are paired within each clone
# Wilcoxon because the distribution we are comparing are not normal/Gaussian 
test_results <- df_plot %>%
  filter(!is.na(intra_fol_distance), !is.na(inter_fol_distance)) %>%
  group_by(patient) %>%
  summarise(
    n_pairs = n(),
    p_value = wilcox.test(intra_fol_distance, inter_fol_distance, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_label = ifelse(p_value < 0.001, "p < 0.001", glue("p = {signif(p_value, 2)}")))

print(test_results)

# Prep label with stats 
label_pos <- df_plot %>%
  select(patient, intra_fol_distance, inter_fol_distance) %>%
  pivot_longer(cols = contains("dist"), values_to = "distance") %>%
  group_by(patient) %>%
  summarise(y = max(distance, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  left_join(test_results, by = "patient")

# Plot 
df_plot %>%
  select(!dist_median) %>% 
  pivot_longer(
    cols = contains("dist"),
    names_to = "distance_type",
    values_to = "distance"
  ) %>% 
  ggplot(aes(x = patient, y = distance, color = distance_type)) + 
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75), alpha = 0.5) +
  geom_text(data = label_pos, aes(x = patient, y = y, label = p_label), inherit.aes = FALSE, size = 3.5) + 
  scale_color_discrete(labels = c(
    dist_median = "Median clone",
    inter_fol_distance = "Median inter-follicular",
    intra_fol_distance = "Median intra-follicular"
  )) +
  labs(
    title = glue("Full sequence, top {n_clones} GC B cell clones"),
    subtitle = "P-values: paired Wilcoxon test, intra- vs. inter-follicular distance",
    y = "GCtree patristic distance (# mutations)",
    x = "Clone nr", 
    color = "Distance measure"
  ) +
  theme_bw() 

ggsave(glue("48_GCtree/plot_06_analyse_trees/combined_top_20_distance_median_boxplot_gctree.png"))



# ggsave(glue("{outdir}/{HH}_top_10_distance_median_boxplot_gctree.png"), width = 10, height = 6)


# ------------------------------------------------------------------------------
# Isotype switching along GCtree edges (nearest observed-ancestor comparison)
# ------------------------------------------------------------------------------

outdir_isotype_switching <- glue("{outdir}/isotype_switching")
dir.create(outdir_isotype_switching, recursive = TRUE, showWarnings = FALSE)

isotype_switch_order <- c("IGHM/D", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")

# for a given observed node, walk up the tree until hitting another observed ("seq*") node
get_nearest_observed_ancestor <- function(node, parent_of) {
  current <- parent_of[node]
  while (!is.na(current) && !str_detect(current, "^seq[0-9]+$")) {
    current <- parent_of[current]
  }
  if (is.na(current)) NA_character_ else current
}

# resolves each observed node's isotype top-down: a mixed node prefers continuity
# with its resolved ancestor's isotype (assumes no switch) over blindly picking the
# earliest isotype in the switch order
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

lapply(c(10, 20), function(n_clones){
  
  # n_clones <- 20
  switch_results_all <- tibble()
  
  for (clone_nr in n_clones_begin:n_clones) {
    
    # HH <- "HH117"
    # clone_nr <- 10
    
    clone <- top_clones %>% dplyr::slice(clone_nr) %>% pull(clone_subgroup_id_90_similarity)
    clone_full_name <- glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")
    sequence_ids <- seq_dir[[clone_full_name]]
    
    df_clone <- resolve_LC_list_germlined[[HH]] %>% 
      filter(locus == "IGH", clone_subgroup_id_90_similarity == clone) %>% 
      mutate(sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)])
    
    tree_file <- glue("48_GCtree/plot/{clone_full_name}/{clone_full_name}.inference.1.nk")
    idmap_file <- glue("48_GCtree/out/{clone_full_name}/idmap.txt")
    
    tree <- read.tree(tree_file)
    idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id"), show_col_types = FALSE) %>% 
      filter(node != "GL") %>% 
      separate_longer_delim(sequence_id, delim = ":")
    
    node_isotype_map <- idmap %>% 
      left_join(df_clone %>% select(sequence_id_gctree, c_call_grouped), by = c("sequence_id" = "sequence_id_gctree")) %>% 
      group_by(node) %>% 
      summarise(isotypes = list(unique(na.omit(c_call_grouped))), .groups = "drop")
    isotype_lookup <- set_names(node_isotype_map$isotypes, node_isotype_map$node)
    
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
        node_mixed = map_lgl(node_isotypes, ~ length(.x) > 1),
        ancestor_mixed = map_lgl(ancestor_isotypes, ~ length(.x) > 1),
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
      ungroup() %>% 
      mutate(clone = clone, clone_nr = clone_nr)
    
    switch_results_all <- bind_rows(switch_results_all, switch_df)
    
  }
  
  switch_results_all %>% dplyr::count(clone_nr, switch_type)
  
  # Plot result - Bar plot 
  switch_type_colors <- list("none" = "grey", "sequential" = "forestgreen", "reverse" = "purple")
  
  switch_results_all %>%
    mutate(
      clone_plot = glue("{clone_nr}\n({clone})") %>% fct_reorder(clone_nr), 
      switch_type_fct = factor(switch_type, levels = c("none", "sequential", "reverse"))
    ) %>% 
    dplyr::count(clone_plot, switch_type_fct, .drop = FALSE) %>%  # .drop = FALSE keeps zero-count combos
    ggplot(aes(x = clone_plot, y = n, fill = switch_type_fct)) + 
    geom_col(position = position_dodge2(preserve = "single"), width = 0.7) + 
    theme_bw() + 
    scale_fill_manual(values = switch_type_colors, drop = FALSE) + 
    labs(
      title = glue("{HH}: Isotype switching events per clone for top {n_clones} clones"),
      y = "N transitions",
      x = "Clone nr", 
      fill = "Switch type"
    )
  
  ggsave(glue("{outdir_isotype_switching}/{HH}_isotype_switching_barplot_{n_clones}.png"), width = 14, height = 5)
  
  # Zooming in on "sequential" and "reverse" isotype switching
  switch_transition_clones <- switch_results_all %>% 
    filter(switch_type %in% c("sequential", "reverse")) %>% 
    distinct(clone_nr, clone, ancestor_isotype, node_isotype, switch_type) %>% 
    group_by(ancestor_isotype, node_isotype, switch_type) %>% 
    summarise(
      n_clones = n_distinct(clone_nr),
      clone_nrs = paste(sort(unique(clone_nr)), collapse = ", "),
      clones = paste(sort(unique(clone)), collapse = ", "),
      .groups = "drop"
    ) %>% 
    arrange(desc(n_clones))
  
  switch_transition_clones
  write_delim(
    switch_transition_clones, delim = "\t", 
    file = glue("{outdir_isotype_switching}/{HH}_isotype_switching_heatmap_details_{n_clones}.tsv")
  )
  
  plot_data <- switch_results_all %>% 
    filter(switch_type %in% c("sequential", "reverse")) %>% 
    dplyr::count(node_isotype, ancestor_isotype, switch_type) %>% 
    right_join(
      expand_grid(
        ancestor_isotype = isotype_switch_order,
        node_isotype = isotype_switch_order,
        switch_type = c("sequential", "reverse")
      ),
      by = c("node_isotype", "ancestor_isotype", "switch_type")
    ) %>% 
    mutate(
      n = replace_na(n, 0),
      ancestor_isotype = factor(ancestor_isotype, levels = isotype_switch_order),
      node_isotype = factor(node_isotype, levels = isotype_switch_order)
    )
  
  shared_max <- max(plot_data$n)
  
  ggplot() + 
    geom_tile(
      data = plot_data %>% filter(switch_type == "sequential"), 
      aes(x = node_isotype, y = ancestor_isotype, fill = n), 
      color = "white"
    ) + 
    scale_fill_gradient(
      low = "white", high = "forestgreen", name = "N sequential", 
      limits = c(0, shared_max), breaks = seq(0, shared_max, 3)
    ) + 
    new_scale_fill() + 
    geom_tile(
      data = plot_data %>% filter(switch_type == "reverse", n > 0), 
      aes(x = node_isotype, y = ancestor_isotype, fill = n), 
      color = "white"
    ) + 
    scale_fill_gradient(
      low = "white", high = "purple", name = "N reverse", 
      limits = c(0, shared_max), breaks = seq(0, shared_max, 3)
    ) + 
    geom_text(
      data = plot_data %>% filter(n > 0),
      aes(x = node_isotype, y = ancestor_isotype, label = n), 
      size = 3.5
    ) + 
    scale_x_discrete(position = "top") + 
    labs(
      x = "Descendant isotype", 
      y = "Nearest observed ancestor isotype",
      title = glue("{HH}: Sequential vs reverse isotype transitions (top {n_clones} clones)"), 
      caption = "Isotype switch order: IGHM → IGHD → IGHG3 → IGHG1 → IGHA1 → IGHG2 → IGHG4 → GHE → IGHA2"
    ) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0))
  
  ggsave(glue("{outdir_isotype_switching}/{HH}_isotype_switching_heatmap_{n_clones}.png"), width = 7.5, height = 5)
  
})



# ------------------------------------------------------------------------------
# Do we see PCs or memory cells downstream of GC B cells?
# ------------------------------------------------------------------------------

outdir_celltype_transition <- glue("{outdir}/celltype_transition")
dir.create(outdir_celltype_transition, recursive = TRUE, showWarnings = FALSE)

# check what cell types actually appear across your clones first
df <- resolve_LC_list_germlined[[HH]] %>% 
  filter(locus == "IGH") %>% 
  mutate(
    tissue = sample_clean %>% str_remove_all(glue("{HH}-")),
    L1_annotation_sample = glue("{L1_annotation}_{tissue}")
  )

lapply(c(10, 20), function(n_clones){
  
  celltype_results_all <- tibble()
  
  for (clone_nr in n_clones_begin:n_clones) {
    
    # clone_nr <- 1
    
    clone <- top_clones %>% dplyr::slice(clone_nr) %>% pull(clone_subgroup_id_90_similarity)
    clone_full_name <- glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")
    sequence_ids <- seq_dir[[clone_full_name]]
    
    # NOTE: not filtering on L1_annotation_sample here, so all cell types in this clone are kept
    df_clone <- df %>% 
      filter(clone_subgroup_id_90_similarity == clone) %>% 
      mutate(sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)])
    
    tree_file <- glue("48_GCtree/plot/{clone_full_name}/{clone_full_name}.inference.1.nk")
    idmap_file <- glue("48_GCtree/out/{clone_full_name}/idmap.txt")
    
    tree <- read.tree(tree_file)
    idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id"), show_col_types = FALSE) %>% 
      filter(node != "GL") %>% 
      separate_longer_delim(sequence_id, delim = ":")
    
    node_celltype_map <- idmap %>% 
      left_join(df_clone %>% select(sequence_id_gctree, L1_annotation_sample), by = c("sequence_id" = "sequence_id_gctree")) %>% 
      group_by(node) %>% 
      summarise(celltypes = list(unique(na.omit(L1_annotation_sample))), .groups = "drop")
    celltype_lookup <- set_names(node_celltype_map$celltypes, node_celltype_map$node)
    
    node_labels <- c(tree$tip.label, tree$node.label)
    edge_df <- as_tibble(tree$edge) %>% 
      set_names(c("parent_idx", "child_idx")) %>% 
      mutate(parent = node_labels[parent_idx], child = node_labels[child_idx])
    parent_of <- set_names(edge_df$parent, edge_df$child)
    
    observed_nodes <- node_labels[str_detect(node_labels, "^seq[0-9]+$")]
    
    celltype_df <- tibble(node = observed_nodes) %>% 
      mutate(
        ancestor = map_chr(node, ~ get_nearest_observed_ancestor(.x, parent_of)),
        node_celltypes = celltype_lookup[node],
        ancestor_celltypes = celltype_lookup[ancestor],
        node_mixed = map_lgl(node_celltypes, ~ length(.x) > 1),
        ancestor_mixed = map_lgl(ancestor_celltypes, ~ length(.x) > 1)
      ) %>% 
      filter(!is.na(ancestor)) %>% 
      unnest_longer(ancestor_celltypes, values_to = "ancestor_celltype") %>% 
      unnest_longer(node_celltypes, values_to = "node_celltype") %>% 
      filter(!is.na(node_celltype), !is.na(ancestor_celltype)) %>% 
      mutate(
        transition = glue("{ancestor_celltype} -> {node_celltype}"),
        gc_to_output = ancestor_celltype == "GC_B_cells" & node_celltype != "GC_B_cells",
        clone = clone,
        clone_nr = clone_nr,
        patient_id = HH
      )
    
    celltype_results_all <- bind_rows(celltype_results_all, celltype_df)
    
  }
  
  # ---- summary: which clones show GC_B_cells -> PC/Memory transitions ----
  
  celltype_results_all %>% 
    filter(gc_to_output) %>% 
    dplyr::count(clone_nr, clone, ancestor_celltype, node_celltype, sort = TRUE)
  
  # ---- overall transition counts, any direction ----
  
  celltype_results_all %>% 
    dplyr::count(ancestor_celltype, node_celltype, sort = TRUE)
  
  # Table
  celltype_transition_clones <- celltype_results_all %>% 
    filter(ancestor_celltype != node_celltype) %>%   # only actual transitions, not same-type "none"
    distinct(patient_id, clone_nr, clone, ancestor_celltype, node_celltype) %>% 
    group_by(ancestor_celltype, node_celltype) %>% 
    summarise(
      n_clones = n_distinct(clone_nr),
      clone_nrs = paste(sort(unique(clone_nr)), collapse = ", "),
      clones = paste(sort(unique(clone)), collapse = ", "),
      .groups = "drop"
    ) %>% 
    arrange(desc(n_clones))
  
  celltype_transition_clones
  
  write_delim(
    celltype_transition_clones, delim = "\t", 
    file = glue("{outdir_celltype_transition}/{HH}_celltype_transition_heatmap_details_{n_clones}.tsv")
  )
  
  # Plot
  celltype_order <- df$L1_annotation_sample %>% unique() %>% sort()
  
  celltype_plot_data <- celltype_results_all %>% 
    dplyr::count(ancestor_celltype, node_celltype) %>% 
    right_join(
      expand_grid(
        ancestor_celltype = celltype_order,
        node_celltype = celltype_order
      ),
      by = c("ancestor_celltype", "node_celltype")
    ) %>% 
    mutate(
      n = replace_na(n, 0),
      ancestor_celltype = factor(ancestor_celltype, levels = celltype_order),
      node_celltype = factor(node_celltype, levels = celltype_order)
    )
  
  ggplot(celltype_plot_data, aes(x = node_celltype, y = ancestor_celltype, fill = n)) + 
    geom_tile(color = "white") + 
    geom_text(aes(label = if_else(n > 0, as.character(n), "")), size = 3.5) + 
    scale_fill_gradient(low = "white", high = "steelblue", name = "N transitions") + 
    scale_x_discrete(position = "top") + 
    labs(
      x = "Descendant cell type",
      y = "Nearest observed ancestor cell type",
      title = glue("{HH}: Cell type transitions along top {n_clones} GCtrees")
    ) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0))
  
  ggsave(glue("{outdir_celltype_transition}/{HH}_celltype_transition_heatmap_{n_clones}.png"), width = 11, height = 8)
  
  
})