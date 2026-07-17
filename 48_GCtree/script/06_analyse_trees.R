library(tidyverse)
library(glue)
library(ape)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_90_similarity_germlined.rds")
seq_dir <- readRDS("48_GCtree/out_90_similarity/rds/seq_dir.rds")

patients <- names(resolve_LC_list_germlined)

HH <- "HH119" 

outdir <- glue("48_GCtree/plot_90_similarity_06_analyse_trees/{HH}")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Get top clones
top_clones <- resolve_LC_list_germlined[[HH]] %>% 
  filter(
    locus == "IGH", 
    L1_annotation == "GC_B_cells"
  ) %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) 


# ==============================================================================
# Using GCtree nk files
# ==============================================================================

# ------------------------------------------------------------------------------
# Tree-based (patristic) distance matrix from GCtree output, aligned to df_clone row order
# ------------------------------------------------------------------------------

# clone_nr <- 3
clone_nr <- 8

# Get clone
clone <- top_clones %>% 
  dplyr::slice(clone_nr) %>% 
  pull(clone_subgroup_id_90_similarity)

# Get clone files 
clone_full_name <- glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")

sequence_ids <- seq_dir[[clone_full_name]]
df_clone$sequence_alignment == sequence_ids

# filter for clone 
df_clone <- resolve_LC_list_germlined[[HH]] %>% 
  filter(
    locus == "IGH", 
    clone_subgroup_id_90_similarity == clone
  ) %>% 
  mutate(
    sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)]
  )

tree_file <- glue("48_GCtree/plot_90_similarity/{clone_full_name}/{clone_full_name}.inference.1.nk")
idmap_file <- glue("48_GCtree/out_90_similarity/{clone_full_name}/idmap.txt")

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

# Define clones to plot
n_clones <- 10
these_clones <- top_clones %>% 
  head(n_clones) %>% 
  pull(clone_subgroup_id_90_similarity)

# Initiate df_plot 
df_plot <- data.frame()

# for (clone_nr in 1:n_clones){
for (clone_nr in 2:n_clones){
  
  # clone_nr <- 3
  
  # Get clone
  clone <- top_clones %>% 
    dplyr::slice(clone_nr) %>% 
    pull(clone_subgroup_id_90_similarity)
  
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
  tree_file <- glue("48_GCtree/plot_90_similarity/{clone_full_name}/{clone_full_name}.inference.1.nk")
  idmap_file <- glue("48_GCtree/out_90_similarity/{clone_full_name}/idmap.txt")
  
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

df_plot %>%
  mutate(
    clone_plot = glue("{clone_nr}\n({clone})") %>% fct_reorder(clone_nr)
  ) %>% 
  pivot_longer(
    cols = contains("dist"),
    names_to = "distance_type",
    values_to = "distance"
  ) %>% 
  ggplot(aes(x = clone_plot, y = distance, color = distance_type)) + 
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75), alpha = 0.5) + 
  scale_color_discrete(labels = c(
    dist_median = "Median clone",
    inter_fol_distance = "Median inter-follicular",
    intra_fol_distance = "Median intra-follicular"
  )) +
  labs(
    title = glue("{HH}: Full sequence, top {n_clones} GC B cell clones"),
    y = "GCtree patristic distance (# mutations)",
    x = "Clone nr", 
    color = "Distance measure"
  ) +
  theme_bw() 

ggsave(glue("{outdir}/{HH}_top_10_distance_median_boxplot_gctree.png"), width = 10, height = 6)




