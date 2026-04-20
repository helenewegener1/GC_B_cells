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

spec_clones_vj <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")

# ------------------------------------------------------------------------------
# Extract sequences
# ------------------------------------------------------------------------------

HH <- "HH117"

# Define patient
spec_clones_vj_HH <- spec_clones_vj[[HH]]

# Look at top heavy chain clones defined with SCOPer
spec_clones_vj_HH %>% 
  count(clone_id, v_call_majority, j_call_majority, sort = TRUE)

# Get top 10 clones
top_clones <- spec_clones_vj_HH %>% count(clone_id, sort = TRUE) %>% head(n = 5) %>% pull(clone_id)

# Check N cells in top 10 clones
spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% nrow()

# Get genes for the clones 
v_genes <- spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% pull(v_call_majority) %>% unique() %>% str_split(",") %>% unlist() %>% unique()
j_genes <- spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% pull(j_call_majority) %>% unique() %>% str_split(",") %>% unlist() %>% unique()

# Build a regex pattern with escaped asterisks
v_gene_regex <- paste(str_replace_all(v_genes, "\\*", "\\\\*"), collapse = "|")
j_gene_regex <- paste(str_replace_all(j_genes, "\\*", "\\\\*"), collapse = "|")

spec_clones_vj_HH %>% 
  filter(str_detect(v_call_majority, v_gene_regex)) %>% 
  count(v_call_majority, sort = TRUE)

# Get all sequences
seqs <- spec_clones_vj_HH %>% filter(str_detect(v_call_majority, v_gene_regex), str_detect(j_call_majority, j_gene_regex)) %>% pull(sequence)
seq_names <- spec_clones_vj_HH %>% filter(str_detect(v_call_majority, v_gene_regex), str_detect(j_call_majority, j_gene_regex)) %>% pull(sequence_id)

seqs %>% length()
seqs %>% unique() %>% length()

# Get metadata of sequences
seqs_meta <- spec_clones_vj_HH %>% 
  filter(str_detect(v_call_majority, v_gene_regex), str_detect(j_call_majority, j_gene_regex)) %>% 
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

# N j_calls

# N v_calls

v_gene_clans <- c(
  "IGKV1" = "K_clan_1",
  "IGKV2" = "K_clan_2", "IGKV3" = "K_clan_2", "IGKV4" = "K_clan_2", "IGKV6" = "K_clan_2",
  "IGKV5" = "K_clan_3", "IGKV7" = "K_clan_3",
  "IGLV1" = "L_clan_1", "IGLV2" = "L_clan_1", "IGLV6" = "L_clan_1", "IGLV10" = "L_clan_1",
  "IGLV3" = "L_clan_2",
  "IGLV7" = "L_clan_3", "IGLV8" = "L_clan_3",
  "IGLV5" = "L_clan_4", "IGLV11" = "L_clan_4",
  "IGLV4" = "L_clan_5", "IGLV9" = "L_clan_5",
  "IGHV" = 
)


## Wrangle metadata to only color the top 10 clones
seqs_meta <- seqs_meta %>% 
  mutate(
    # clone_id_top = ifelse(clone_subgroup_id %in% top_subclones, paste(clone_subgroup_id, v_gene, j_gene, junction_length, sep = "_"), "other"),
    clone_id_top = ifelse(clone_subgroup_id %in% top_subclones, clone_subgroup_id, "other"),
    v_gene = v_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1),
    j_gene = j_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1), 
    v_j_junction = paste(v_gene, j_gene, junction_length, sep = "_"), 
    v_gene_subgroup = v_gene %>% str_split_i("-", 1),
    j_gene_subgroup = j_gene %>% str_split_i("-", 1), 
    v_gene_clan = coalesce(v_gene_clans[v_gene_subgroup], v_gene_subgroup)
  )

seqs_meta$clone_id_top %>% table(useNA = "always") 
# top_clones_names <- seqs_meta$clone_id_top %>% table() %>% names()

seqs_meta$v_gene_clan %>% table(useNA = "always") 


seqs_meta$clone_id_top_10 %>% table(useNA = "always") 

clone_colors <- setNames(
  c("#E63946",  
    "#2196F3",  
    "#4CAF50",  
    "#FF9800",  
    "#9C27B0",  
    # "#00BCD4",  
    # "#FFEB3B",  
    # "#FF4081",  
    # "#795548",  
    # "#76FF03", 
    "grey90"), 
  c(top_clones, "other")
) # NAs will be grey

ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = clone_id_top_10), size=0.1, alpha=0.8) +
  theme_tree2() + 
  scale_color_manual(values = clone_colors) + 
  guides(color = guide_legend(override.aes = list(size = 4))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} colored by clone_id (N clones: {n_clones})"))

ggsave(glue("46_sequence_driven_clustering/plot/{HH}_clone_ID_fan.png"), width = 18, height = 16, dpi = 1000)

ggtree(tree_phylo, layout="rectangular", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = clone_id_top_10), size=0.1, alpha=0.8) +
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

# Color by V call majority
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = v_call_majority), size=0.5, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(title = glue("{HH} colored by V gene")) 
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_v_call_majority.png"), width = 18, height = 16, dpi = 1000)

# Color by J gene
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = j_gene), size=0.2, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(title = glue("{HH} colored by J gene")) 
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_J_gene.png"), width = 18, height = 16, dpi = 1000)

# Color by J call majority
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = j_call_majority), size=0.5, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(title = glue("{HH} colored by J gene")) 
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_j_call_majority.png"), width = 18, height = 16, dpi = 1000)

# Color by junction length
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = junction_length), size=0.5, alpha=0.8) +
  theme_tree2() + 
  labs(title = glue("{HH} colored by junction length"))
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_junction_length.png"), width = 18, height = 16, dpi = 1000)

# Color by v_j_junction
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = v_j_junction), size=0.5, alpha=0.8) +
  theme_tree2() + 
  labs(title = glue("{HH} colored by v_j_junction"))
ggsave(glue("46_sequence_driven_clustering/plot/{HH}_v_j_junction.png"), width = 22, height = 16, dpi = 1000)



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



