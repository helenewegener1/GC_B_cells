library(tidyverse)
library(stringdist)
library(stats)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(fastcluster)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")

# ------------------------------------------------------------------------------
# Extract sequences
# ------------------------------------------------------------------------------

HH <- "HH117"

# Define patient
spec_clones_vj_HH <- spec_clones_vj[[HH]] %>% 
  mutate(
    v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
  ) %>% 
  mutate(
    j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
  )

# Look at top heavy chain clones defined with SCOPer
spec_clones_vj_HH %>% 
  count(clone_id, v_gene, j_gene, sort = TRUE)

# Get top 10 clones
top_clones <- spec_clones_vj_HH %>% count(clone_id, sort = TRUE) %>% head(n = 10) %>% pull(clone_id)

# Check N cells in top 10 clones
spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% nrow()

# Get genes for the clones 
v_genes <- spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% pull(v_gene) %>% unique()
j_genes <- spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% pull(j_gene) %>% unique()

# Get all sequences
seqs <- spec_clones_vj_HH %>% filter(v_gene %in% v_genes, j_gene %in% j_genes) %>% pull(sequence)
seq_names <- spec_clones_vj_HH %>% filter(v_gene %in% v_genes, j_gene %in% j_genes) %>% pull(sequence_id)

seqs %>% length()
seqs %>% unique() %>% length()

# Get metadata of sequences
seqs_meta <- spec_clones_vj_HH %>% 
  filter(v_gene %in% v_genes, j_gene %in% j_genes) %>% 
  mutate(
    label = sequence_id, 
    junction_length = as.character(junction_length)
  ) %>% 
  select(label, everything()) %>% # Move label to front
  as.data.frame()

# ------------------------------------------------------------------------------
# Compute Levenshtein as sequences are not of the same length
# ------------------------------------------------------------------------------

# 'lv' is Levenshtein; 'qgram' or 'jaccard' are faster alternatives
dist_matrix <- stringdistmatrix(seqs, seqs, method = "lv")

# ------------------------------------------------------------------------------
# Hierarchical Clustering
# ------------------------------------------------------------------------------

fit <- hclust(as.dist(dist_matrix), method = "complete")

# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

# Convert your hclust object to a 'phylo' object
tree_phylo <- as.phylo(fit)

length(tree_phylo$tip.label)
length(seq_names)
tree_phylo$tip.label <- seq_names

# Plot the tree using a 'fan' layout (best for high density)
# 'mapping' connects the tree tips to your metadata
# ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
#   geom_tippoint(aes(color = sample_clean_fol), size=0.5, alpha=0.8) +
#   theme_tree2()
# 
# ggsave(glue("46_sequence_driven_clustering/plot/{HH}_sample_clean_fol.png"))

# Color by clone ID

## N clones in this pool of sequences
n_clones <- seqs_meta$clone_id %>% unique() %>% length()

## Wrangle metadata to only color the top 10 clones
seqs_meta <- seqs_meta %>% 
  mutate(
    clone_id_top_10 = ifelse(clone_id %in% top_clones, clone_id, "other")
  )

seqs_meta$clone_id_top_10 %>% table(useNA = "always") 

clone_colors <- setNames(
  c("#E63946",  
    "#2196F3",  
    "#4CAF50",  
    "#FF9800",  
    "#9C27B0",  
    "#00BCD4",  
    "#FFEB3B",  
    "#FF4081",  
    "#795548",  
    "#76FF03", 
    "grey90"), 
  c(top_clones, "other")
) # NAs will be grey

ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = clone_id_top_10), size=0.5, alpha=0.8) +
  theme_tree2() + 
  scale_color_manual(values = clone_colors) + 
  guides(color = guide_legend(override.aes = list(size = 4))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} colored by clone_id (N clones: {n_clones})"))

ggsave(glue("46_sequence_driven_clustering/plot/{HH}_clone_ID_fan.png"), width = 18, height = 16, dpi = 1000)

ggtree(tree_phylo, layout="rectangular", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = clone_id_top_10), size=0.5, alpha=0.8) +
  theme_tree2() + 
  scale_color_manual(values = clone_colors) + 
  guides(color = guide_legend(override.aes = list(size = 4))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} colored by clone_id (N clones: {n_clones})"))

ggsave(glue("46_sequence_driven_clustering/plot/{HH}_clone_ID_rectangular.png"), width = 10, height = 30, dpi = 1000)


# Color by V gene
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = v_gene), size=0.5, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(title = glue("{HH} colored by V gene")) 
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_v_gene.png"), width = 18, height = 16, dpi = 1000)

ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = j_gene), size=0.5, alpha=0.8) +
  theme_tree2() + 
  labs(title = glue("{HH} colored by J gene"))
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_j_gene.png"), width = 18, height = 16, dpi = 1000)

ggtree(tree_phylo, layout="fan", size=0.5) %<+% seqs_meta + 
  geom_tippoint(aes(color = junction_length), size=0.5, alpha=0.8) +
  theme_tree2() + 
  labs(title = glue("{HH} colored by junction length"))
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_junction_length.png"), width = 18, height = 16, dpi = 1000)




# Other layouts
# 'rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight', 'ape'
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = v_gene), size=0.5, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(title = glue("{HH} colored by V gene")) 

ggtree(tree_phylo, layout="rectangular", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = v_gene), size=0.5, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(title = glue("{HH} colored by V gene")) 



# ggsave(glue("46_sequence_driven_clustering/plot/{HH}_v_gene.png"), width = 18, height = 16, dpi = 1000)


# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

# Cut tree
clusters <- cutree(fit, k = 10 + 5)

table(clusters)

table(clusters, seqs_meta$clone_id_top_10)



