library(tidyverse)
library(stringdist)
library(stats)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(fastcluster)
library(glue)

version <- "HH119_two_largest_clones"
# version <- "HH119_two_largest_clones_junction"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list_gmm_threshold_germlined.rds")
table(resolve_LC_list$locus)

# top 10 GC clones for each patient 
# all_clones <- list(
#   "HH117" = c(
#     "4221_1", "2628_1", "1849_1", "3709_1", "2301_1",
#     "1320_1", "5941_1", "6115_1", "1910_1", "2169_1"
#   ),
#   
#   "HH119" = c(
#     "28075_1", "12120_1", "15287_1", "23124_1",
#     "8372_1", "3869_1", "25158_1", "7913_1"
#   )
# )

# ------------------------------------------------------------------------------
# Subset data
# ------------------------------------------------------------------------------

HH <- str_split_i(version, "_", 1)

clones <- resolve_LC_list[[HH]] %>% 
  filter(locus == "IGH") %>% 
  count(clone_subgroup_id, sort = TRUE) %>% 
  slice(1:2) %>% 
  pull(clone_subgroup_id)

subset <- resolve_LC_list[[HH]] %>% 
  filter(
    locus == "IGH",
    clone_subgroup_id %in% clones
  )

nrow(subset)
subset$clone_subgroup_id %>% table()

# ------------------------------------------------------------------------------
# Extract sequences
# ------------------------------------------------------------------------------

# Add metadata 
subset <- subset %>% mutate(
  v_gene = v_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1),
  j_gene = j_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1), 
  v_j_junction = paste(v_gene, j_gene, junction_length, sep = "_"), 
  v_gene_subgroup = v_gene %>% str_split_i("-", 1),
  j_gene_subgroup = j_gene %>% str_split_i("-", 1),
  # v_gene_clan = coalesce(v_gene_clans[v_gene_subgroup], v_gene_subgroup)
)

# Subset based on clones and genes 
# # Look at top heavy chain clones defined with SCOPer
subset %>%
  dplyr::count(clone_id, v_call, j_call, sort = TRUE)

# Subset heavy chain and plasma cells and trim sequence 
subset <- subset %>% 
  mutate(
    sequence_trimmed = str_sub(sequence, v_sequence_start, j_sequence_end),
    sequence_trimmed_300 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-299, nchar(sequence_trimmed)),
    sequence_trimmed_250 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-249, nchar(sequence_trimmed))
  )

# Have a look 
subset$sequence_trimmed %>% nchar() %>% range()
subset$sequence_trimmed %>% nchar() %>% hist()

subset$sequence_trimmed_300 %>% nchar() %>% unique()
subset$sequence_trimmed_250 %>% nchar() %>% unique()

# Extract sequences
seqs <- subset$sequence_trimmed_300
# seqs <- subset$junction
seq_names <- subset$sequence_id

# Get metadata of sequences
seqs_meta <- subset %>% 
  mutate(
    label = sequence_id, 
    junction_length = as.character(junction_length),
    clone_subgroup = as.character(clone_subgroup)
  ) %>% 
  select(label, everything()) %>% # Move label to front
  as.data.frame()

# ------------------------------------------------------------------------------
# Compute Levenshtein as sequences are not of the same length
# ------------------------------------------------------------------------------

# Get distances 
dist_mat <- stringdistmatrix(seqs, method = "hamming")
# dist_mat <- stringdistmatrix(seqs, method = "lv")

# ------------------------------------------------------------------------------
# Hierarchical Clustering
# ------------------------------------------------------------------------------

fit <- hclust(as.dist(dist_mat), method = "complete")

# saveRDS(fit, glue("46_sequence_driven_clustering/out/fit_{version}_{version}.rds"))

# ==============================================================================
# no-zoom tree
# ==============================================================================

outdir_trees <- glue("45_immcantation/17_look_at_clones/{version}/trees")
dir.create(outdir_trees, showWarnings = FALSE, recursive = TRUE)

# Convert to phylo for nicer plotting
tree <- as.phylo(fit)

# Checks 
length(tree$tip.label)
length(seqs_meta$sequence_id)
tree$tip.label <- seqs_meta$sequence_id


# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------


variables <- c(
  "v_gene_subgroup", "j_gene_subgroup", "v_call_subgroup", "j_call_subgroup",
  "v_call",  "j_call",
  "clone_id_plot", "junction_length"
)
#
for (var in variables){
  
  # var <- "clone_subgroup_id"
  # var <- "v_call"
  # var <- "v_gene_subgroup"
  # var <- "junction_length"
  
  ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta + 
    geom_tippoint(aes(color = !!sym(var)), size=0.5, alpha = 0.5) +
    theme_tree2() + 
    guides(color = guide_legend(override.aes = list(size = 4))) + 
    labs(
      title = glue("{version} colored by {var}")
      # subtitle = glue("Clusters: {clusters_string}")
    )
    
  ggsave(glue("{outdir_trees}/{version}_{var}.png"), width = 15, height = 10, dpi = 1000)

}


# ------------------------------------------------------------------------------
# Cut the tree
# ------------------------------------------------------------------------------

k <- 8

# outdir <- glue("46_sequence_driven_clustering/plot/{version}/{k}_clusters")
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Cut tree
clusters <- cutree(fit, k = k)
table(clusters)

# Add to metadata
seqs_meta$clusters <- clusters %>% as.character()

# See patient split
table(seqs_meta$clusters, seqs_meta$clone_subgroup_id)

# Show clusters 
var <- "clusters"
ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = !!sym(var)), size=0.5, alpha = 0.5) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(
    title = glue("{version} colored by {k}_{var}")
    # subtitle = glue("Clusters: {clusters_string}")
  )

ggsave(glue("{outdir_trees}/{version}_{k}_{var}.png"), width = 15, height = 10, dpi = 1000)

# ==============================================================================
# Zoom 1: Investigate individual cluster
# ==============================================================================

outdir_trees <- glue("{outdir}/trees")
dir.create(outdir_trees, showWarnings = FALSE, recursive = TRUE)

# Zoom in to a cluster
cl_zoom_1 <- c(14)
idx <- which(clusters %in% cl_zoom_1)
seqs_sub <- seqs[idx]

# Get meta data of the cluster 
seqs_meta_zoom_1 <- seqs_meta %>% filter(clusters %in% cl_zoom_1)

# Calculate new tree
dist_sub <- as.dist(stringdistmatrix(seqs_sub, method = "hamming"))
fit_sub  <- hclust(dist_sub, method = "complete")

# Convert to phylo for nicer plotting
tree <- as.phylo(fit_sub)

# Checks 
length(tree$tip.label)
length(seqs_meta_zoom_1$sequence_id)
tree$tip.label <- seqs_meta_zoom_1$sequence_id

# ------------------------------------------------------------------------------
# Define clone colors 
# ------------------------------------------------------------------------------

# Wrangle metadata
# these_clones <- seqs_meta_zoom_1 %>% dplyr::count(clone_id, sort = TRUE) %>% head(50) %>% pull(clone_id)
seqs_meta <- seqs_meta %>% mutate(clone_id_plot = ifelse(clone_id %in% clones, clone_subgroup_id, "other"))

# clone_id_plot colors
clone_colors <- setNames(
  c(
    "#E63946", "#2196F3",
    "grey90"
  ),
  c(clones, "other")
) # NAs will be grey


# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

clusters_string <- paste(c(cl_zoom_1), collapse = ", ")
clusters_string_path <- paste0(c(cl_zoom_1), collapse = "_")

variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "v_call_subgroup", "j_call_subgroup", 
               "junction_length", "patient_id", "clone_id_plot", "L1_annotation")

for (var in variables){
  
  # var <- "clone_subgroup_id"
  # var <- "v_gene_subgroup"
  p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta + 
    geom_tippoint(aes(color = !!sym(var)), size=0.5, alpha = 0.5) +
    theme_tree2() + 
    guides(color = guide_legend(override.aes = list(size = 4))) + 
    labs(
      title = glue("{version} colored by {var} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    )
  
  if (var == "clone_id_plot"){
    p <- p + scale_color_manual(values = clone_colors)
  }
  
  ggsave(glue("{outdir}/trees/{version}_{k}_clusters_TREE_clusters_{clusters_string_path}_{var}.png"), plot = p, width = 15, height = 10, dpi = 1000)
  
}
