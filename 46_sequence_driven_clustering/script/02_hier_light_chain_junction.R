library(tidyverse)
library(stringdist)
library(stats)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(fastcluster)

version <- "ligth_chain_junction"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

spec_clones_vj <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")
resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list.rds")

resolve_LC_list$HH117$j_call_10x %>% head()
resolve_LC_list$HH117$j_call %>% head()


grep("IGHJ5\\*02", resolve_LC_list$HH117$j_call, value = T) %>% unique()
grep("IGHJ4\\*02", resolve_LC_list$HH117$j_call)

index <- grep("IGHJ1\\*01,IGHJ4\\*02,IGHJ5\\*02", resolve_LC_list$HH117$j_call)

resolve_LC_list$HH117$j_call_10x[index]


grep("\\-", resolve_LC_list$HH117$v_call, value = T) %>% unique()

index <- grep("IGHV1-69\\*01,IGHV1-69\\*18,IGHV1-69D\\*01", resolve_LC_list$HH117$v_call)

resolve_LC_list$HH117$v_call_10x[index] %>% unique()


# ------------------------------------------------------------------------------
# Extract sequences
# ------------------------------------------------------------------------------

HH <- "HH117"

# Define patient
spec_clones_vj_HH <- resolve_LC_list[[HH]]

# Look at top heavy chain clones defined with SCOPer
spec_clones_vj_HH %>% 
  filter(locus == "IGH") %>% 
  count(clone_id, v_call_majority, j_call_majority, sort = TRUE)

# Look at top HL chain clones defined with SCOPer
spec_clones_vj_HH %>% 
  filter(locus != "IGH") %>% 
  count(clone_subgroup_id, sort = TRUE)

# Get top clones
top_clones <- spec_clones_vj_HH %>% count(clone_id, sort = TRUE) %>% slice(2) %>% pull(clone_id)

# Check N cells in top clones
spec_clones_vj_HH %>% filter(clone_id == top_clones & locus == "IGH") %>% nrow()

# Get genes for the clones 
# v_genes <- spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% pull(v_call_majority) %>% unique() %>% str_split(",") %>% unlist() %>% unique()
# j_genes <- spec_clones_vj_HH %>% filter(clone_id %in% top_clones) %>% pull(j_call_majority) %>% unique() %>% str_split(",") %>% unlist() %>% unique()
# 
# Build a regex pattern with escaped asterisks
# v_gene_regex <- paste(str_replace_all(v_genes, "\\*", "\\\\*"), collapse = "|")
# j_gene_regex <- paste(str_replace_all(j_genes, "\\*", "\\\\*"), collapse = "|")
# 
# spec_clones_vj_HH %>% 
#   filter(str_detect(v_call_majority, v_gene_regex)) %>% 
#   count(v_call_majority, sort = TRUE)
# 
# Get all sequences
# seqs <- spec_clones_vj_HH %>% filter(str_detect(v_call_majority, v_gene_regex), str_detect(j_call_majority, j_gene_regex)) %>% pull(sequence)
# seq_names <- spec_clones_vj_HH %>% filter(str_detect(v_call_majority, v_gene_regex), str_detect(j_call_majority, j_gene_regex)) %>% pull(sequence_id)

seqs <- spec_clones_vj_HH %>% filter(clone_id == top_clones & locus != "IGH") %>% pull(junction)
# seq_names <- spec_clones_vj_HH %>% filter(clone_id == top_clones & locus != "IGH") %>% pull(junction)

seqs %>% length()
seqs %>% unique() %>% length()

# Length of sequences 
map(seqs, nchar) %>% unlist() %>% table()

# Get metadata of sequences
seqs_meta <- spec_clones_vj_HH %>% 
  filter(clone_id == top_clones & locus != "IGH") %>% 
  mutate(
    label = junction, 
  ) %>% 
  select(label, everything()) %>% # Move label to front
  as.data.frame()

# Check
seqs_meta$clone_id %>% table()
seqs_meta$clone_subgroup_id %>% table()

# ------------------------------------------------------------------------------
# Compute Levenshtein as sequences are not of the same length
# ------------------------------------------------------------------------------

dist_matrix <- stringdistmatrix(seqs_meta$junction, method = "lv")

# Convert to full matrix first
dist_mat <- as.matrix(dist_matrix)

# Normalize by the longer of the two sequences
lengths <- nchar(seqs_meta$junction)
length_mat <- outer(lengths, lengths, pmax)
dist_mat_norm <- dist_mat / length_mat

# Convert back to dist object for hclust
dist_matrix_norm <- as.dist(dist_mat_norm)

# ------------------------------------------------------------------------------
# Hierarchical Clustering
# ------------------------------------------------------------------------------

fit <- hclust(as.dist(dist_matrix_norm), method = "complete")

# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

# Convert your hclust object to a 'phylo' object
tree_phylo <- as.phylo(fit)

length(tree_phylo$tip.label)
length(seqs)
tree_phylo$tip.label <- seqs

# Plot the tree using a 'fan' layout (best for high density)
# 'mapping' connects the tree tips to your metadata
# ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
#   geom_tippoint(aes(color = sample_clean_fol), size=0.5, alpha=0.8) +
#   theme_tree2()
# 
# ggsave(glue("46_sequence_driven_clustering/plot/{HH}_sample_clean_fol.png"))

# Color by clone ID

## N clones in this pool of sequences
n_clones <- seqs_meta$clone_subgroup_id %>% unique() %>% length()

# Define top subclone
top_subclones <- seqs_meta %>% count(clone_subgroup_id, sort = TRUE) %>% head(5) %>% pull(clone_subgroup_id)

#   IGKV clans: 
# clan I: Homo sapiens IGKV1 subgroup genes
# clan II: Homo sapiens IGKV2 , IGKV3, IGKV4 and IGKV6 subgroup genes
# clan III: Homo sapiens IGKV5 and IGKV7 subgroup genes
#   IGLV clans:
# clan I: Homo sapiens IGLV1, IGLV2, IGLV6 and IGLV10 subgroup genes, and pseudogenes IGLV(I)-20, -38, -42, -56, -63, -68 and -70
# clan II: Homo sapiens IGLV3 subgroup genes;
# clan III: Homo sapiens IGLV7 and IGLV8 subgroup genes
# clan IV: Homo sapiens IGLV5 and IGLV11 subgroup genes, and pseudogenes IGLV(IV)-53, -59, -64, -65 and -66-1
# clan V: Homo sapiens IGLV4 and IGLV9 subgroup genes, and pseudogenes IGLV(V)-58 and -66
# clan VI: Homo sapiens pseudogenes IGLV(VI)-22-1 and -25-1
# clan VII: Homo sapiens pseudogene IGLV(VII)-41-1

v_gene_clans <- c(
  "IGKV1" = "K_clan_1",
  "IGKV2" = "K_clan_2", "IGKV3" = "K_clan_2", "IGKV4" = "K_clan_2", "IGKV6" = "K_clan_2",
  "IGKV5" = "K_clan_3", "IGKV7" = "K_clan_3",
  "IGLV1" = "L_clan_1", "IGLV2" = "L_clan_1", "IGLV6" = "L_clan_1", "IGLV10" = "L_clan_1",
  "IGLV3" = "L_clan_2",
  "IGLV7" = "L_clan_3", "IGLV8" = "L_clan_3",
  "IGLV5" = "L_clan_4", "IGLV11" = "L_clan_4",
  "IGLV4" = "L_clan_5", "IGLV9" = "L_clan_5"
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
  c(top_subclones, "other")
  # top_clones_names[c(1,2,3)]
) # NAs will be grey


# -----------------------
# Plot trees
# -----------------------

ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = clone_id_top), size=1, alpha=0.8) +
  theme_tree2() + 
  scale_color_manual(values = clone_colors) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} H clone {top_clones} colored by clone_id (N light clones: {n_clones})"),
       subtitle = version)

ggsave(glue("46_sequence_driven_clustering/plot/{version}/{HH}_{top_clones}_clone_ID_fan.png"), width = 9, height = 6.5, dpi = 1000)

# Color by v_gene_clan
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = v_gene_clan), size=1, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} H clone {top_clones}"),
       subtitle = version)

ggsave(glue("46_sequence_driven_clustering/plot/{version}/{HH}_{top_clones}_v_gene_clan.png"), width = 9, height = 6.5, dpi = 1000)

# Color by v_gene_subgroup
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = v_gene_subgroup), size=1, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} H clone {top_clones}"),
       subtitle = version)

ggsave(glue("46_sequence_driven_clustering/plot/{version}/{HH}_{top_clones}_v_gene_subgroup.png"), width = 9, height = 6.5, dpi = 1000)


# Color by locus
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = locus), size=1, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} H clone {top_clones}"),
       subtitle = version)

ggsave(glue("46_sequence_driven_clustering/plot/{version}/{HH}_{top_clones}_locus.png"), width = 9, height = 6.5, dpi = 1000)


# Color by J call
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = j_call), size=1, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} H clone {top_clones}"),
       subtitle = version)

ggsave(glue("46_sequence_driven_clustering/plot/{version}/{HH}_{top_clones}_j_call.png"), width = 9, height = 6.5, dpi = 1000)

# Color by J gene
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = j_gene), size=1, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} H clone {top_clones}"),
       subtitle = version)

ggsave(glue("46_sequence_driven_clustering/plot/{version}/{HH}_{top_clones}_j_gene.png"), width = 9, height = 6.5, dpi = 1000)

# Color by j_gene_subgroup
ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
  geom_tippoint(aes(color = j_gene_subgroup), size=1, alpha=0.8) +
  theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  # theme(legend.position = "none") + 
  labs(title = glue("{HH} H clone {top_clones}"),
       subtitle = version)

ggsave(glue("46_sequence_driven_clustering/plot/{version}/{HH}_{top_clones}_j_gene_subgroup.png"), width = 9, height = 6.5, dpi = 1000)




# ggsave(glue("46_sequence_driven_clustering/plot/{HH}_v_gene.png"), width = 18, height = 16, dpi = 1000)

# ------------------------------------------------------------------------------
# Distance to nearest 
# ------------------------------------------------------------------------------

dist_matrix <- as.dist(dist_matrix) %>% as.matrix()
diag(dist_matrix) <- NA

# Get the minimum distance to a neighbor for each sequence (closest neighbor)
dist_to_nearest <- apply(dist_matrix, 1, min, na.rm = TRUE)

hist(dist_to_nearest, breaks = 100)

# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

# Cut tree
clusters <- cutree(fit, k = 10 + 5)

table(clusters)

table(clusters, seqs_meta$clone_id_top_10)



