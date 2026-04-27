library(tidyverse)
library(stringdist)
library(stats)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(fastcluster)
library(glue)

version <- "00_PC_tree"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

spec_clones_vj <- readRDS("45_immcantation/out/rds/resolve_LC_list.rds")

v_gene_clans <- c(
  "IGHV1" = "H_clan_1", "IGHV5" = "H_clan_1", "IGHV7" = "H_clan_1",
  "IGHV2" = "H_clan_2", "IGHV4" = "H_clan_2", "IGHV6" = "H_clan_2",
  "IGHV3" = "H_clan_3",
  "IGKV1" = "K_clan_1",
  "IGKV2" = "K_clan_2", "IGKV3" = "K_clan_2", "IGKV4" = "K_clan_2", "IGKV6" = "K_clan_2",
  "IGKV5" = "K_clan_3", "IGKV7" = "K_clan_3",
  "IGLV1" = "L_clan_1", "IGLV2" = "L_clan_1", "IGLV6" = "L_clan_1", "IGLV10" = "L_clan_1",
  "IGLV3" = "L_clan_2",
  "IGLV7" = "L_clan_3", "IGLV8" = "L_clan_3",
  "IGLV5" = "L_clan_4", "IGLV11" = "L_clan_4",
  "IGLV4" = "L_clan_5", "IGLV9" = "L_clan_5"
)

# ------------------------------------------------------------------------------
# Extract sequences
# ------------------------------------------------------------------------------

HH <- "HH117"

# Define patient
spec_clones_vj_HH <- spec_clones_vj[[HH]]

# Add metadata 
spec_clones_vj_HH <- spec_clones_vj_HH %>% mutate(
  v_gene = v_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1),
  j_gene = j_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1), 
  v_j_junction = paste(v_gene, j_gene, junction_length, sep = "_"), 
  v_gene_subgroup = v_gene %>% str_split_i("-", 1),
  j_gene_subgroup = j_gene %>% str_split_i("-", 1), 
  v_gene_clan = coalesce(v_gene_clans[v_gene_subgroup], v_gene_subgroup)
)

# Subset based on clones and genes 
# # Look at top heavy chain clones defined with SCOPer
spec_clones_vj_HH %>%
  count(clone_id, v_gene_clan, j_gene_subgroup, sort = TRUE)

# Subset heavy chain and plasma cells and trim sequence 
PCs_heavy <- spec_clones_vj_HH %>% 
  filter(locus == "IGH" & L1_annotation == "PCs") %>% 
  mutate(
    sequence_trimmed = str_sub(sequence, v_sequence_start, j_sequence_end),
    sequence_trimmed_300 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-299, nchar(sequence_trimmed)),
    sequence_trimmed_250 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-249, nchar(sequence_trimmed))
  )

# Have a look 
PCs_heavy$sequence_trimmed %>% nchar() %>% range()
PCs_heavy$sequence_trimmed %>% nchar() %>% hist()

PCs_heavy$sequence_trimmed_300 %>% nchar() %>% unique()
PCs_heavy$sequence_trimmed_250 %>% nchar() %>% unique()

# Extract sequences
seqs <- PCs_heavy$sequence_trimmed_300
# seqs <- PCs_heavy$sequence_trimmed_250
seq_names <- PCs_heavy$sequence_id

# Get metadata of sequences
seqs_meta <- PCs_heavy %>% 
  mutate(
    label = sequence_id, 
    junction_length = as.character(junction_length)
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

# saveRDS(fit, glue("46_sequence_driven_clustering/out/fit_{version}_{HH}.rds"))

# # ------------------------------------------------------------------------------
# # Visualize clustering as tree
# # ------------------------------------------------------------------------------
# 
# # Convert your hclust object to a 'phylo' object
# tree_phylo <- as.phylo(fit)
# 
# length(tree_phylo$tip.label)
# length(seq_names)
# tree_phylo$tip.label <- seq_names
# 
# # Plot the tree using a 'fan' layout (best for high density)
# # 'mapping' connects the tree tips to your metadata
# # ggtree(tree_phylo, layout="fan", size=0.2) %<+% seqs_meta + 
# #   geom_tippoint(aes(color = sample_clean_fol), size=0.5, alpha=0.8) +
# #   theme_tree2()
# # 
# # ggsave(glue("46_sequence_driven_clustering/plot/{HH}_sample_clean_fol.png"))
# 
# top_clones <- PCs_heavy %>% count(clone_id, sort = TRUE) %>% head(10) %>% pull(clone_id)
# top_clones
# 
# # Color by clone ID
# 
# ## N clones in this pool of sequences
# n_clones <- seqs_meta$clone_id %>% unique() %>% length()
# 
# # clan I: Homo sapiens IGHV1, IGHV5 and IGHV7 subgroup genes
# # clan II: Homo sapiens IGHV2, IGHV4 and IGHV6 subgroup genes
# # clan III: Homo sapiens IGHV3 subgroup genes
# 
# v_gene_clans <- c(
#   "IGHV1" = "H_clan_1", "IGHV5" = "H_clan_1", "IGHV7" = "H_clan_1",
#   "IGHV2" = "H_clan_2", "IGHV4" = "H_clan_2", "IGHV6" = "H_clan_2",
#   "IGHV3" = "H_clan_3"
# )
# 
# ## Wrangle metadata to only color the top 10 clones
# seqs_meta <- seqs_meta %>% 
#   mutate(
#     # clone_id_top = ifelse(clone_subgroup_id %in% top_subclones, paste(clone_subgroup_id, v_gene, j_gene, junction_length, sep = "_"), "other"),
#     clone_id_top = ifelse(clone_id %in% top_clones, clone_id, "other")
#   )
# 
# seqs_meta$clone_id_top %>% table(useNA = "always") 
# # top_clones_names <- seqs_meta$clone_id_top %>% table() %>% names()
# 
# seqs_meta$v_gene_clan %>% table(useNA = "always") 
# 
# clone_colors <- setNames(
#   c("#E63946",  
#     "#2196F3",  
#     "#4CAF50",  
#     "#FF9800",  
#     "#9C27B0",  
#     "#00BCD4",
#     "#FFEB3B",
#     "#FF4081",
#     "#795548",
#     "#76FF03",
#     "grey90"), 
#   c(top_clones, "other")
# ) # NAs will be grey
# 
# top_clones <- "top_5"

# ------------------------------------------------------------------------------
# Cut the tree
# ------------------------------------------------------------------------------

k <- 20

outdir <- glue("46_sequence_driven_clustering/plot/{version}/{HH}/{k}_clusters")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Cut tree
clusters <- cutree(fit, k = k)
table(clusters)

# Add to metadata
seqs_meta$clusters <- clusters %>% as.character()

# ------------------------------------------------------------------------------
# Investigate clusters
# ------------------------------------------------------------------------------

# table(seqs_meta$clusters, seqs_meta$v_gene_subgroup) %>% 
#   as.data.frame() %>% 
#   rename(Cluster = Var1, Count = Freq) %>% 
#   group_by(Cluster) %>% 
#   mutate(Prop = Count/sum(Count)) %>% 
#   ungroup() %>% 
#   ggplot(aes(x = Var2, y = Cluster, fill = Prop)) +
#   geom_tile(color = "white", linewidth = 0.4) +
#   geom_text(aes(label = round(Prop, 1)), size = 3, color = "black") + 
#   scale_fill_gradient(low = "#EEF3FB", high = "#185FA5") +
#   labs(
#     x     = "V gene subgroup",
#     y     = "Cluster",
#     fill  = "Count",
#     title = glue("Cluster × V gene subgroup ({k} clusters)")
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x  = element_text(angle = 45, hjust = 1),
#     panel.grid   = element_blank(),
#     legend.position = "none"
#   )

table(seqs_meta$clusters, seqs_meta$v_gene_subgroup) %>% 
  as.data.frame() %>% 
  rename(Cluster = Var1, Count = Freq) %>% 
  ggplot(aes(x = Var2, y = Cluster, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Count), size = 3, color = "black") + 
  scale_fill_gradient(low = "#EEF3FB", high = "#185FA5") +
  labs(
    x     = "V gene subgroup",
    y     = "Cluster",
    fill  = "Count",
    title = glue("{HH}: Cluster × V gene subgroup ({k} clusters)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "none"
  )
ggsave(glue("{outdir}/{HH}_{k}_clusters_V_gene.png"))
  
table(seqs_meta$clusters, seqs_meta$j_gene_subgroup) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Freq), size = 3, color = "black") + 
  scale_fill_gradient(low = "#EEF3FB", high = "#185FA5") +
  labs(
    x     = "J gene subgroup",
    y     = "Cluster",
    fill  = "Count",
    title = glue("{HH}: Cluster × J gene subgroup ({k} clusters)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "none"
  )
ggsave(glue("{outdir}/{HH}_{k}_clusters_J_gene.png"))

table(seqs_meta$clusters, seqs_meta$c_call) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Freq), size = 3, color = "black") + 
  scale_fill_gradient(low = "#EEF3FB", high = "#185FA5") +
  labs(
    x     = "Isotype subgroup",
    y     = "Cluster",
    fill  = "Count",
    title = glue("{HH}: Cluster × Isotype ({k} clusters)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "none"
  )
ggsave(glue("{outdir}/{HH}_{k}_clusters_Isotypes.png"))

table(seqs_meta$clusters, seqs_meta$junction_length) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Freq), size = 3, color = "black") + 
  scale_fill_gradient(low = "#EEF3FB", high = "#185FA5") +
  labs(
    x     = "Junction length",
    y     = "Cluster",
    fill  = "Count",
    title = glue("{HH}: Cluster × Junciton length ({k} clusters)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "none"
  )
ggsave(glue("{outdir}/{HH}_{k}_clusters_junction_length.png"), width = 12.5, height = 8)



# ------------------------------------------------------------------------------
# Investigate individual cluster
# ------------------------------------------------------------------------------

outdir_trees <- glue("{outdir}/trees")
dir.create(outdir_trees, showWarnings = FALSE, recursive = TRUE)

cl <- c(6)
idx <- which(clusters %in% cl)
seqs_sub <- seqs[idx]

seqs_meta_cl <- seqs_meta %>% filter(clusters %in% cl)

dist_sub <- as.dist(stringdistmatrix(seqs_sub, method = "hamming"))
fit_sub  <- hclust(dist_sub, method = "complete")

# Convert to phylo for nicer plotting
tree <- as.phylo(fit_sub)

length(tree$tip.label)
length(seqs_meta_cl$sequence_id)
tree$tip.label <- seqs_meta_cl$sequence_id

# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

# ggtree(tree, layout="rectangular", size=0.2) %<+% seqs_meta_cl + 
#   geom_tippoint(aes(color = clone_id_top), size=0.5, alpha=0.8) +
#   theme_tree2() + 
#   guides(color = guide_legend(override.aes = list(size = 4))) + 
#   labs(title = glue("{HH} colored by clone_id_top")) 

clusters_string <- paste(cl, collapse = ", ")
clusters_string_path <- paste0(cl, collapse = "_")

variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "junction_length")
for (var in variables){
  
  # var <- "junction_length"
  ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_cl + 
    geom_tippoint(aes(color = !!sym(var)), size=0.5, alpha=0.8) +
    theme_tree2() + 
    guides(color = guide_legend(override.aes = list(size = 4))) + 
    labs(
      title = glue("{HH} colored by {var} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    )
  
  ggsave(glue("{outdir}/trees/{HH}_{k}_clusters_TREE_clusters_{clusters_string_path}_{var}.png"), width = 9, height = 6.5, dpi = 1000)
  
}



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ZOOMING IN
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

outdir_trees_zoom <- glue("{outdir}/trees_zoom")
dir.create(outdir_trees_zoom, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1 cluster - 1 junction length
# ------------------------------------------------------------------------------

# 1 cluster
cl <- c(6)
idx <- which(clusters %in% cl)
seqs_sub <- seqs[idx]

seqs_meta_cl <- seqs_meta %>% filter(clusters %in% cl)

# 1 junction length
jl <- c(60)
idx <- which(seqs_meta_cl$junction_length %in% jl)
seqs_sub_2 <- seqs_sub[idx]

# Update metadata
seqs_meta_final <- seqs_meta_cl %>% filter(junction_length %in% jl)
table(length(seqs_sub_2) == nrow(seqs_meta_final))

# TREE MAKING
dist_sub <- as.dist(stringdistmatrix(seqs_sub_2, method = "hamming"))
fit_sub  <- hclust(dist_sub, method = "complete")

# Convert to phylo for nicer plotting
tree <- as.phylo(fit_sub)

length(tree$tip.label)
length(seqs_meta_final$sequence_id)

tree$tip.label <- seqs_meta_final$sequence_id

# Plotting
clusters_string <- paste("Cluster: ", cl, "\nJunction length: ", jl)
clusters_string_path <- paste0("cluster", cl, "_", "jl", jl)

outdir_trees_zoom_this_cl_jl <- glue("{outdir_trees_zoom}/{clusters_string_path}")
dir.create(outdir_trees_zoom_this_cl_jl, showWarnings = FALSE, recursive = TRUE)

variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "junction_length")
for (var in variables){
  
  # var <- "v_gene_subgroup"
  # var <- "j_gene_subgroup"
  ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_final +
    geom_tippoint(aes(color = !!sym(var)), size=0.5, alpha=0.8) +
    theme_tree2() + 
    guides(color = guide_legend(override.aes = list(size = 4))) + 
    labs(
      title = glue("{HH} colored by {var} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    )
  
  ggsave(glue("{outdir_trees_zoom_this_cl_jl}/{HH}_{k}_clusters_{clusters_string_path}_{var}.png"), width = 9, height = 6.5, dpi = 1000)
  
}


# ------------------------------------------------------------------------------
# 1 cluster - 1 junction length - 1 v gene
# ------------------------------------------------------------------------------

# 1 cluster
cl <- c(6)
idx <- which(clusters %in% cl)
seqs_sub <- seqs[idx]

seqs_meta_cl <- seqs_meta %>% filter(clusters %in% cl)

# 1 junction length
jl <- c(60)
idx <- which(seqs_meta_cl$junction_length %in% jl)
seqs_sub_2 <- seqs_sub[idx]

seqs_meta_jl <- seqs_meta_cl %>% filter(junction_length %in% jl)

# 1 V gene
v_this_gene <- c("IGHV1")
idx <- which(seqs_meta_jl$v_gene_subgroup %in% v_gene)
seqs_sub_3 <- seqs_sub_2[idx]

# Update metadata
seqs_meta_final <- seqs_meta_jl %>% filter(v_gene_subgroup %in% v_this_gene)
table(length(seqs_sub_3) == nrow(seqs_meta_final))

# TREE MAKING
dist_sub <- as.dist(stringdistmatrix(seqs_sub_3, method = "hamming"))
fit_sub  <- hclust(dist_sub, method = "complete")

# Convert to phylo for nicer plotting
tree <- as.phylo(fit_sub)

length(tree$tip.label)
length(seqs_meta_final$sequence_id)

tree$tip.label <- seqs_meta_final$sequence_id

# Plotting
clusters_string <- paste("Cluster: ", cl, "\nJunction length: ", jl, "\nV gene subgroup: ", v_this_gene)
clusters_string_path <- paste0("cluster", cl, "_", "jl", jl, "_", v_this_gene)

outdir_trees_zoom_this_cl_jl_v <- glue("{outdir_trees_zoom}/{clusters_string_path}")
dir.create(outdir_trees_zoom_this_cl_jl_v, showWarnings = FALSE, recursive = TRUE)

variables <- c("clusters", "v_gene_subgroup", "v_gene", "v_call", "j_gene_subgroup", "j_call", "junction_length")
for (var in variables){
  
  # var <- "j_call"
  # var <- "j_gene_subgroup"
  ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_final +
    geom_tippoint(aes(color = !!sym(var)), size=1, alpha=0.8) +
    theme_tree2() + 
    guides(color = guide_legend(override.aes = list(size = 4))) + 
    labs(
      title = glue("{HH} colored by {var} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    )
  
  ggsave(glue("{outdir_trees_zoom_this_cl_jl_v}/{HH}_{k}_clusters_TREE_clusters_{clusters_string_path}_{var}.png"), width = 9, height = 6.5, dpi = 1000)
  
}
