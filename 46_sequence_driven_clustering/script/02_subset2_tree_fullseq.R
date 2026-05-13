library(tidyverse)
library(stringdist)
library(stats)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(fastcluster)
library(glue)

version <- "subset_2"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

subset_clones_vj_all <- readRDS(glue("46_sequence_driven_clustering/out/{version}_clones_vj.rds"))

table(subset_clones_vj_all$locus)

subset_clones_vj <- subset_clones_vj_all %>% filter(locus == "IGH")

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

# Add metadata 
subset_clones_vj <- subset_clones_vj %>% mutate(
  v_gene = v_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1),
  j_gene = j_call %>% str_split_i(",", 1) %>% str_split_i("\\*", 1), 
  v_j_junction = paste(v_gene, j_gene, junction_length, sep = "_"), 
  v_gene_subgroup = v_gene %>% str_split_i("-", 1),
  j_gene_subgroup = j_gene %>% str_split_i("-", 1),
  v_gene_clan = coalesce(v_gene_clans[v_gene_subgroup], v_gene_subgroup)
)

# Subset based on clones and genes 
# # Look at top heavy chain clones defined with SCOPer
subset_clones_vj %>%
  dplyr::count(clone_id, v_gene_clan, j_gene_subgroup, sort = TRUE)

# Subset heavy chain and plasma cells and trim sequence 
subset_heavy <- subset_clones_vj %>% 
  mutate(
    sequence_trimmed = str_sub(sequence, v_sequence_start, j_sequence_end),
    sequence_trimmed_300 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-299, nchar(sequence_trimmed)),
    sequence_trimmed_250 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-249, nchar(sequence_trimmed))
  )

# Have a look 
subset_heavy$sequence_trimmed %>% nchar() %>% range()
subset_heavy$sequence_trimmed %>% nchar() %>% hist()

subset_heavy$sequence_trimmed_300 %>% nchar() %>% unique()
subset_heavy$sequence_trimmed_250 %>% nchar() %>% unique()

# Extract sequences
seqs <- subset_heavy$sequence_trimmed_300
# seqs <- subset_heavy$sequence_trimmed_250
seq_names <- subset_heavy$sequence_id

# Get metadata of sequences
seqs_meta <- subset_heavy %>% 
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

# ------------------------------------------------------------------------------
# Cut the tree
# ------------------------------------------------------------------------------

k <- 20

outdir <- glue("46_sequence_driven_clustering/plot/{version}/{k}_clusters")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Cut tree
clusters <- cutree(fit, k = k)
table(clusters)

# Add to metadata
seqs_meta$clusters <- clusters %>% as.character()

# See patient split
table(seqs_meta$clusters, seqs_meta$patient_id)

# ------------------------------------------------------------------------------
# Investigate clusters
# ------------------------------------------------------------------------------

vars <- c("patient_id", "L1_annotation", "c_call", "v_gene_subgroup", "j_gene_subgroup", "junction_length")

for (var in vars){
  
  # var <- "junction_length"
  
  table(seqs_meta$clusters, seqs_meta[[var]]) %>% 
    as.data.frame() %>% 
    dplyr::rename(Cluster = Var1, Count = Freq) %>% 
    group_by(Var2) %>% 
    mutate(Pct = Count / sum(Count) * 100) %>% 
    ungroup() %>% 
    ggplot(aes(x = Var2, y = Cluster, fill = Pct)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = glue("{round(Pct, 2)}%")), size = 3, color = "black") + 
    scale_fill_gradient(low = "white", high = "#185FA5") +
    labs(
      x     = var,
      y     = "Cluster",
      fill  = "Percentage per group",
      title = glue("{version}: Cluster × {var} ({k} clusters)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.grid   = element_blank(),
      legend.position = "none"
    )
  
  ggsave(glue("{outdir}/{version}_{k}_clusters_{var}.png"), width = 16, height = 8)
  
}


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
these_clones <- seqs_meta_zoom_1 %>% dplyr::count(clone_id, sort = TRUE) %>% head(50) %>% pull(clone_id)
seqs_meta_zoom_1 <- seqs_meta_zoom_1 %>% mutate(clone_id_plot = ifelse(clone_id %in% these_clones, clone_id, "other"))

# clone_id_plot colors
clone_colors <- setNames(
  c(
    "#E63946", "#2196F3", "#3DAA55", "#FF9800", "#9C27B0", "#00BCD4", "#F5C518", "#FF4081", "#6D4C41", "#76FF03",
    "#1565C0", "#00897B", "#558B2F", "#6200EA", "#37474F", "#FF6F00", "#00E5FF", "#D500F9", "#AEEA00", "#BF360C",
    "#1DE9B6", "#FF1744", "#AA00FF", "#FFD600", "#0091EA", "#F06292", "#43A047", "#E040FB", "#FF6D00", "#26C6DA",
    "#8D6E63", "#C6FF00", "#283593", "#FF3D00", "#00BFA5", "#827717", "#4A148C", "#33691E", "#B71C1C", "#006064",
    "#F48FB1", "#80DEEA", "#CCFF90", "#B39DDB", "#FFCC80", "#EF9A9A", "#80CBC4", "#CE93D8", "#FFF176", "#A5D6A7",
    "grey90"
  ),
  c(these_clones, "other")
) # NAs will be grey


# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

clusters_string <- paste(c(cl_zoom_1), collapse = ", ")
clusters_string_path <- paste0(c(cl_zoom_1), collapse = "_")

variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "v_call_subgroup", "j_call_subgroup", 
               "junction_length", "patient_id", "clone_id_plot", "L1_annotation")

for (var in variables){
  
  # var <- "clone_id"
  # var <- "v_gene_subgroup"
  p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_zoom_1 + 
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

# ==============================================================================
# Zoom 2
# ==============================================================================

outdir_trees_zoom <- glue("{outdir}/trees_zoom")
dir.create(outdir_trees_zoom, showWarnings = FALSE, recursive = TRUE)

k <- 3
 
# Cut tree
clusters <- cutree(fit_sub, k = k)
table(clusters)

# Add to metadata
seqs_meta_zoom_2 <- seqs_meta_zoom_1
seqs_meta_zoom_2$clusters <- clusters %>% as.character()

# See patient split
table(seqs_meta_zoom_2$clusters, seqs_meta_zoom_2$v_gene_subgroup)

# Define cluster to zoom into 
cl_zoom_2 <- c(1)

# Zoom in to a cluster
idx <- which(clusters %in% cl_zoom_2)
seqs_sub_2 <- seqs_sub[idx]

# Get meta data of the cluster 
seqs_meta_zoom_2 <- seqs_meta_zoom_2 %>% filter(clusters %in% cl_zoom_2)

# Calculate new tree
dist_sub <- as.dist(stringdistmatrix(seqs_sub_2, method = "hamming"))
fit_sub  <- hclust(dist_sub, method = "complete")

# Convert to phylo for nicer plotting
tree <- as.phylo(fit_sub)

# Checks 
length(tree$tip.label)
length(seqs_meta_zoom_2$sequence_id)
tree$tip.label <- seqs_meta_zoom_2$sequence_id

# ------------------------------------------------------------------------------
# Define clone colors 
# ------------------------------------------------------------------------------

# Wrangle metadata
these_clones <- seqs_meta_zoom_2 %>% dplyr::count(clone_id, sort = TRUE) %>% head(50) %>% pull(clone_id)
seqs_meta_zoom_2 <- seqs_meta_zoom_2 %>% mutate(clone_id_plot = ifelse(clone_id %in% these_clones, clone_id, "other"))

# clone_id_plot colors
clone_colors <- setNames(
  c(
    "#E63946", "#2196F3", "#3DAA55", "#FF9800", "#9C27B0", "#00BCD4", "#F5C518", "#FF4081", "#6D4C41", "#76FF03",
    "#1565C0", "#00897B", "#558B2F", "#6200EA", "#37474F", "#FF6F00", "#00E5FF", "#D500F9", "#AEEA00", "#BF360C",
    "#1DE9B6", "#FF1744", "#AA00FF", "#FFD600", "#0091EA", "#F06292", "#43A047", "#E040FB", "#FF6D00", "#26C6DA",
    "#8D6E63", "#C6FF00", "#283593", "#FF3D00", "#00BFA5", "#827717", "#4A148C", "#33691E", "#B71C1C", "#006064",
    "#F48FB1", "#80DEEA", "#CCFF90", "#B39DDB", "#FFCC80", "#EF9A9A", "#80CBC4", "#CE93D8", "#FFF176", "#A5D6A7",
    "grey90"
  ),
  c(these_clones, "other")
) # NAs will be grey


# ------------------------------------------------------------------------------
# Visualize clustering as tree
# ------------------------------------------------------------------------------

clusters_string <- paste(c(cl_zoom_1, cl_zoom_2), collapse = ", ")
clusters_string_path <- paste0(c(cl_zoom_1, cl_zoom_2), collapse = "_")

variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "v_call_subgroup", "j_call_subgroup", 
               "junction_length", "patient_id", "clone_id_plot", "L1_annotation")

for (var in variables){
  
  # var <- "clone_id"
  # var <- "j_gene_subgroup"
  p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_zoom_2 + 
    geom_tippoint(aes(color = !!sym(var)), size=0.5) +
    theme_tree2() + 
    guides(
      color = guide_legend(override.aes = list(size = 4))
    ) + 
    labs(
      title = glue("{version} colored by {var} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    )
  
  if (var == "clone_id_plot"){
    p <- p + scale_color_manual(values = clone_colors)
  }
  
  ggsave(glue("{outdir_trees_zoom}/{version}_{k}_clusters_TREE_clusters_{clusters_string_path}_{var}.png"), plot = p, width = 15, height = 10, dpi = 1000)
  
}

# Plot the different clone definitions + shape by light chain (clone_subgroup)
clone_definitions <- colnames(seqs_meta_zoom_2)[str_detect(colnames(seqs_meta_zoom_2), "clone")]

lapply(clone_definitions, function(clone_def){
  
  # clone_def <- "clone_id"
  
  # Wrangle metadata
  these_clones <- seqs_meta_zoom_2 %>% dplyr::count(!!sym(clone_def), sort = TRUE) %>% head(50) %>% pull(!!sym(clone_def))
  seqs_meta_zoom_2 <- seqs_meta_zoom_2 %>% mutate(clone_id_plot = ifelse(!!sym(clone_def) %in% these_clones, !!sym(clone_def), "other"))
  
  # clone_id_plot colors
  clone_colors <- setNames(
    c(
      "#E63946", "#2196F3", "#3DAA55", "#FF9800", "#9C27B0", "#00BCD4", "#F5C518", "#FF4081", "#6D4C41", "#76FF03",
      "#1565C0", "#00897B", "#558B2F", "#6200EA", "#37474F", "#FF6F00", "#00E5FF", "#D500F9", "#AEEA00", "#BF360C",
      "#1DE9B6", "#FF1744", "#AA00FF", "#FFD600", "#0091EA", "#F06292", "#43A047", "#E040FB", "#FF6D00", "#26C6DA",
      "#8D6E63", "#C6FF00", "#283593", "#FF3D00", "#00BFA5", "#827717", "#4A148C", "#33691E", "#B71C1C", "#006064",
      "#F48FB1", "#80DEEA", "#CCFF90", "#B39DDB", "#FFCC80", "#EF9A9A", "#80CBC4", "#CE93D8", "#FFF176", "#A5D6A7",
      "grey90"
    ),
    c(these_clones, "other")
  ) # NAs will be grey
  
  
  ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_zoom_2 +
    geom_tippoint(aes(color = clone_id_plot, shape = clone_subgroup), size=0.8) +
    theme_tree2() + 
    guides(
      color = "none", 
      shape = guide_legend(override.aes = list(size = 4))
    ) + 
    labs(
      title = glue("{version} colored by {clone_def} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    ) + 
    scale_color_manual(values = clone_colors)
  
  
  ggsave(glue("{outdir_trees_zoom}/{version}_{k}_clusters_TREE_clusters_{clusters_string_path}_{clone_def}.png"), width = 9, height = 6.5, dpi = 1000)
  
  
})

# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# # ZOOMING IN - OLD
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# 
# outdir_trees_zoom <- glue("{outdir}/trees_zoom")
# dir.create(outdir_trees_zoom, showWarnings = FALSE, recursive = TRUE)
# 
# # # ------------------------------------------------------------------------------
# # # 1 cluster - 1 junction length
# # # ------------------------------------------------------------------------------
# # 
# # # 1 cluster
# # # # cl <- c(1)
# # # idx <- which(clusters %in% cl)
# # # seqs_sub <- seqs[idx]
# # # 
# # # seqs_meta_cl <- seqs_meta %>% filter(clusters %in% cl)
# # 
# # # 1 junction length
# # jl <- c(48)
# # idx <- which(seqs_meta_cl$junction_length %in% jl)
# # seqs_sub_2 <- seqs_sub[idx]
# # 
# # # Update metadata
# # seqs_meta_jl <- seqs_meta_cl %>% filter(junction_length %in% jl)
# # table(length(seqs_sub_2) == nrow(seqs_meta_jl))
# # 
# # # TREE MAKING
# # dist_sub <- as.dist(stringdistmatrix(seqs_sub_2, method = "hamming"))
# # fit_sub  <- hclust(dist_sub, method = "complete")
# # 
# # # Convert to phylo for nicer plotting
# # tree <- as.phylo(fit_sub)
# # 
# # length(tree$tip.label)
# # length(seqs_meta_jl$sequence_id)
# # 
# # tree$tip.label <- seqs_meta_jl$sequence_id
# # 
# # # Plotting
# clusters_string <- paste("Cluster: ", cl, "\nJunction length: ", jl)
# clusters_string_path <- paste0("cluster", cl, "_", "jl", jl)
# 
# outdir_trees_zoom_this_cl_jl <- glue("{outdir_trees_zoom}/{clusters_string_path}")
# dir.create(outdir_trees_zoom_this_cl_jl, showWarnings = FALSE, recursive = TRUE)
# # 
# # variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "v_call_subgroup", "j_call_subgroup", 
# #                "junction_length", "patient_id", "clone_id_plot", "L1_annotation")
# # for (var in variables){
# #   
# #   # var <- "v_gene_subgroup"
# #   # var <- "junction_length"
# #   # var <- "patient_id"
# #   # var <- "clone_id_plot"
# #   p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_jl +
# #     geom_tippoint(aes(color = !!sym(var)), size=0.5, alpha=0.8) +
# #     theme_tree2() + 
# #     guides(color = guide_legend(override.aes = list(size = 4))) + 
# #     labs(
# #       title = glue("{version} colored by {var} - N clusters: {k}"),
# #       subtitle = glue("Clusters: {clusters_string}")
# #     )
# #   
# #   if (var == "clone_id_plot"){
# #     p <- p + scale_color_manual(values = clone_colors)
# #   }
# #   
# #   ggsave(glue("{outdir_trees_zoom_this_cl_jl}/{version}_{k}_clusters_{clusters_string_path}_{var}.png"), plot = p, width = 9, height = 6.5, dpi = 1000)
# #   
# # }
# 
# 
# # ------------------------------------------------------------------------------
# # 1 cluster - 1 junction length - 1 v gene
# # ------------------------------------------------------------------------------
# 
# # 1 cluster
# # cl <- c(6)
# # idx <- which(clusters %in% cl)
# # seqs_sub <- seqs[idx]
# # 
# # seqs_meta_cl <- seqs_meta %>% filter(clusters %in% cl)
# 
# # 1 junction length
# # jl <- c(60)
# # idx <- which(seqs_meta_cl$junction_length %in% jl)
# # seqs_sub_2 <- seqs_sub[idx]
# # 
# # seqs_meta_jl <- seqs_meta_cl %>% filter(junction_length %in% jl)
# 
# # Cut tree
# clusters <- cutree(fit_sub, k = 3)
# table(clusters)
# 
# # Add to metadata
# seqs_meta_cl$clusters <- clusters %>% as.character()
# 
# # See patient split
# table(seqs_meta_cl$clusters, seqs_meta_cl$v_gene_subgroup)
# 
# # Cluster 1 contain only IGH3
# v_this_gene_cluster <- 1
# 
# # 1 V gene
# v_this_gene <- c("IGHV3")
# # idx <- which(seqs_meta_jl$v_gene_subgroup %in% v_this_gene)
# idx <- which(seqs_meta_cl$cluster %in% v_this_gene_cluster)
# # seqs_sub_3 <- seqs_sub_2[idx]
# seqs_sub_3 <- seqs_sub[idx]
# 
# # Update metadata
# # seqs_meta_final <- seqs_meta_jl %>% filter(v_gene_subgroup %in% v_this_gene)
# # table(length(seqs_sub_3) == nrow(seqs_meta_final))
# seqs_meta_final <- seqs_meta_cl %>% filter(clusters == v_this_gene_cluster)
# table(length(seqs_sub_3) == nrow(seqs_meta_final))
# 
# # Wrangle metadata
# these_clones <- seqs_meta_final %>% dplyr::count(clone_id, sort = TRUE) %>% head(50) %>% pull(clone_id)
# seqs_meta_final <- seqs_meta_final %>% mutate(clone_id_plot = ifelse(clone_id %in% these_clones, clone_id, "other"))
# 
# # clone_id_plot colors
# clone_colors <- setNames(
#   c(
#     "#E63946", "#2196F3", "#3DAA55", "#FF9800", "#9C27B0", "#00BCD4", "#F5C518", "#FF4081", "#6D4C41", "#76FF03",
#     "#1565C0", "#00897B", "#558B2F", "#6200EA", "#37474F", "#FF6F00", "#00E5FF", "#D500F9", "#AEEA00", "#BF360C",
#     "#1DE9B6", "#FF1744", "#AA00FF", "#FFD600", "#0091EA", "#F06292", "#43A047", "#E040FB", "#FF6D00", "#26C6DA",
#     "#8D6E63", "#C6FF00", "#283593", "#FF3D00", "#00BFA5", "#827717", "#4A148C", "#33691E", "#B71C1C", "#006064",
#     "#F48FB1", "#80DEEA", "#CCFF90", "#B39DDB", "#FFCC80", "#EF9A9A", "#80CBC4", "#CE93D8", "#FFF176", "#A5D6A7",
#     "grey90"
#   ),
#   c(these_clones, "other")
# ) # NAs will be grey
# 
# # TREE MAKING
# dist_sub <- as.dist(stringdistmatrix(seqs_sub_3, method = "hamming"))
# fit_sub  <- hclust(dist_sub, method = "complete")
# 
# # Convert to phylo for nicer plotting
# tree <- as.phylo(fit_sub)
# 
# length(tree$tip.label)
# length(seqs_meta_final$sequence_id)
# 
# tree$tip.label <- seqs_meta_final$sequence_id
# 
# # Plotting
# clusters_string <- paste("Cluster: ", cl, "\nJunction length: ", jl, "\nV gene subgroup: ", v_this_gene)
# clusters_string_path <- paste0("cluster", cl, "_", "jl", jl, "_", v_this_gene)
# 
# outdir_trees_zoom_this_cl_jl_v <- glue("{outdir_trees_zoom}/{clusters_string_path}")
# dir.create(outdir_trees_zoom_this_cl_jl_v, showWarnings = FALSE, recursive = TRUE)
# 
# variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "v_call_subgroup", "j_call_subgroup", 
#                "junction_length", "patient_id", "clone_id_plot", "L1_annotation")
# for (var in variables){
#   
#   # var <- "j_call"
#   # var <- "clone_id_plot"
#   # var <- "clone_id_junction"
#   p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_final +
#     geom_tippoint(aes(color = !!sym(var)), size=0.5) +
#     theme_tree2() + 
#     guides(color = guide_legend(override.aes = list(size = 4))) + 
#     labs(
#       title = glue("{version} colored by {var} - N clusters: {k}"),
#       subtitle = glue("Clusters: {clusters_string}")
#     )
#   
#   if (var == "clone_id_plot"){
#     p <- p + scale_color_manual(values = clone_colors)
#   }
#   
#   ggsave(glue("{outdir_trees_zoom_this_cl_jl_v}/{version}_{k}_clusters_TREE_clusters_{clusters_string_path}_{var}.png"), plot = p, width = 9, height = 6.5, dpi = 1000)
#   
# }
# 
# 
# # Plot the different clone definitions + shape by light chain (clone_subgroup)
# # clone_definitions <- c("clone_id", "clone_id_junction", "clone_id_vj_junction", "clone_id_vj_junction_95") 
# clone_definitions <- colnames(seqs_meta_final)[str_detect(colnames(seqs_meta_final), "clone")]
# 
# lapply(clone_definitions, function(clone_def){
#   
#   # clone_def <- "clone_id"
#   
#   # Wrangle metadata
#   these_clones <- seqs_meta_final %>% count(!!sym(clone_def), sort = TRUE) %>% head(50) %>% pull(!!sym(clone_def))
#   seqs_meta_final <- seqs_meta_final %>% mutate(clone_id_plot = ifelse(!!sym(clone_def) %in% these_clones, !!sym(clone_def), "other"))
#   
#   # clone_id_plot colors
#   clone_colors <- setNames(
#     c(
#       "#E63946", "#2196F3", "#3DAA55", "#FF9800", "#9C27B0", "#00BCD4", "#F5C518", "#FF4081", "#6D4C41", "#76FF03",
#       "#1565C0", "#00897B", "#558B2F", "#6200EA", "#37474F", "#FF6F00", "#00E5FF", "#D500F9", "#AEEA00", "#BF360C",
#       "#1DE9B6", "#FF1744", "#AA00FF", "#FFD600", "#0091EA", "#F06292", "#43A047", "#E040FB", "#FF6D00", "#26C6DA",
#       "#8D6E63", "#C6FF00", "#283593", "#FF3D00", "#00BFA5", "#827717", "#4A148C", "#33691E", "#B71C1C", "#006064",
#       "#F48FB1", "#80DEEA", "#CCFF90", "#B39DDB", "#FFCC80", "#EF9A9A", "#80CBC4", "#CE93D8", "#FFF176", "#A5D6A7",
#       "grey90"
#     ),
#     c(these_clones, "other")
#   ) # NAs will be grey
#   
#   
#   p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_final +
#     geom_tippoint(aes(color = clone_id_plot, shape = clone_subgroup), size=0.5) +
#     theme_tree2() + 
#     guides(
#       color = "none", 
#       shape = guide_legend(override.aes = list(size = 4))
#     ) + 
#     labs(
#       title = glue("{version} colored by {clone_def} - N clusters: {k}"),
#       subtitle = glue("Clusters: {clusters_string}")
#     ) + 
#     scale_color_manual(values = clone_colors)
#   
#   
#   ggsave(glue("{outdir_trees_zoom_this_cl_jl_v}/{version}_{k}_clusters_TREE_clusters_{clusters_string_path}_{clone_def}.png"), plot = p, width = 9, height = 6.5, dpi = 1000)
# 
#   
# })


# ------------------------------------------------------------------------------
# Cut the tree
# ------------------------------------------------------------------------------

k <- 50

# outdir <- glue("46_sequence_driven_clustering/plot/{version}/{k}_clusters")
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Cut tree
clusters <- cutree(fit_sub, k = k)
table(clusters)

# Add to metadata
seqs_meta_final$clusters <- clusters %>% as.character()


# meta sub
clones_more_than_one <- seqs_meta_final %>% count(clone_id) %>% filter(n > 1) %>% pull(clone_id)
meta_sub <- seqs_meta_final %>% filter(clone_id %in% clones_more_than_one)

# See patient split
table(meta_sub$clusters, meta_sub$clone_id) %>% 
  as.data.frame() %>% 
  dplyr::rename(cluster = Var1, clone_id = Var2) %>%
  ggplot(aes(x = clone_id, y = cluster, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Freq), size = 3, color = "black") + 
  scale_fill_gradient(low = "white", high = "#185FA5") +
  labs(
    # fill  = "Percentage per group",
    title = glue("{version}: clone_id × ({k} clusters)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "none"
  )
# ggsave(glue("{outdir}/{version}_{k}_clusters_{var}.png"), width = 16, height = 8)
  

library(mclust)  # for ARI
library(aricode) # for NMI

# Loop over cuts
results <- data.frame()
for (k in 2:20) {
  clusters <- cutree(fit_sub, k = k)
  ari <- adjustedRandIndex(clusters, seqs_meta_final$clone_id)
  nmi <- NMI(clusters, seqs_meta_final$clone_id)
  results <- rbind(results, data.frame(k = k, ARI = ari, NMI = nmi))
}

# Plot
ggplot(results, aes(x = k)) +
  geom_line(aes(y = ARI, color = "ARI")) +
  geom_line(aes(y = NMI, color = "NMI")) +
  labs(x = "Number of clusters (k)", y = "Score", color = "Metric")

heights <- seq(0.1, 1.0, by = 0.05)
results <- data.frame()
for (h in heights) {
  clusters <- cutree(fit_sub, h = h)
  ari <- adjustedRandIndex(clusters, seqs_meta_final$clone_id)
  results <- rbind(results, data.frame(h = h, k = max(clusters), ARI = ari))
}

library(dendextend)
dend <- as.dendrogram(fit_sub)
dend <- color_branches(dend, k = 5)
dend <- color_labels(dend, col = as.numeric(as.factor(seqs_meta$clone_id)))
plot(dend)















  

# cluster zoom
cl_zoom <- c(1)
idx <- which(clusters %in% cl_zoom)
seqs_sub_4 <- seqs_sub_3[idx]

# Define metadata
seqs_meta_final_sub_4 <- seqs_meta_final %>% filter(clusters %in% cl_zoom)

# TREE MAKING
dist_sub_4 <- as.dist(stringdistmatrix(seqs_sub_4, method = "hamming"))
fit_sub_4  <- hclust(dist_sub_4, method = "complete")

# Convert to phylo for nicer plotting
tree <- as.phylo(fit_sub_4)

length(tree$tip.label)
length(seqs_meta_final_sub_4$sequence_id)

tree$tip.label <- seqs_meta_final_sub_4$sequence_id

# Plotting
clusters_string <- paste("Cluster: ", cl, "\nJunction length: ", jl, "\nV gene subgroup: ", v_this_gene,
                         "\ncluster: ", cl_zoom)
clusters_string_path <- paste0("cluster", cl, "_", "jl", jl, "_", v_this_gene, "_", cl_zoom)

outdir_trees_zoom_this_cl_jl_v <- glue("{outdir_trees_zoom}/{clusters_string_path}")
dir.create(outdir_trees_zoom_this_cl_jl_v, showWarnings = FALSE, recursive = TRUE)

variables <- c("clusters", "v_gene_subgroup", "j_gene_subgroup", "v_call_subgroup", "j_call_subgroup",
               "junction_length", "patient_id", "clone_id_plot", "L1_annotation")

for (var in variables){
  
  # var <- "j_call"
  # var <- "clone_id_plot"
  # var <- "patient_id"
  p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_final_sub_4 +
    geom_tippoint(aes(color = !!sym(var)), size=0.5) +
    theme_tree2() + 
    guides(color = guide_legend(override.aes = list(size = 4))) + 
    labs(
      title = glue("{version} colored by {var} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    )
  
  if (var == "clone_id_plot"){
    p <- p + scale_color_manual(values = clone_colors)
  }
  
  ggsave(glue("{outdir_trees_zoom_this_cl_jl_v}/{version}_{k}_clusters_TREE_clusters_{clusters_string_path}_{var}.png"), plot = p, width = 9, height = 6.5, dpi = 1000)
  
}

# Plot the different clone definitions 
# clone_definitions <- c("clone_id", "clone_id_junction", "clone_id_vj_junction", "clone_id_vj_junction_95") 
clone_definitions <- colnames(seqs_meta_final)[str_detect(colnames(seqs_meta_final), "clone")]

lapply(clone_definitions, function(clone_def){
  
  # clone_def <- "clone_id_vj_junction_95"
  
  # Wrangle metadata
  these_clones <- seqs_meta_final %>% count(!!sym(clone_def), sort = TRUE) %>% head(15) %>% pull(!!sym(clone_def))
  seqs_meta_final <- seqs_meta_final %>% mutate(clone_id_plot = ifelse(!!sym(clone_def) %in% these_clones, !!sym(clone_def), "other"))
  
  # clone_id_plot colors
  clone_colors <- setNames(
    c(
      "#E63946", "#2196F3", "#3DAA55", "#FF9800",
      "#9C27B0", "#00BCD4", "#F5C518", "#FF4081",
      "#6D4C41", "#76FF03", "#1565C0","#00897B",  
      "#558B2F", "#6200EA","#37474F",
      "grey90"
    ),
    c(these_clones, "other")
  ) # NAs will be grey
  
  
  ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_final +
    geom_tippoint(aes(color = clone_id_plot, shape = clone_subgroup), size=0.5) +
    theme_tree2() + 
    guides(
      color = guide_legend(override.aes = list(size = 4)),
      shape = guide_legend(override.aes = list(size = 4))
    ) + 
    labs(
      title = glue("{version} colored by {clone_def} - N clusters: {k}"),
      subtitle = glue("Clusters: {clusters_string}")
    ) + 
    scale_color_manual(values = clone_colors)
  
  
  ggsave(glue("{outdir_trees_zoom_this_cl_jl_v}/{version}_{k}_clusters_TREE_clusters_{clusters_string_path}_{clone_def}.png"), width = 9, height = 6.5, dpi = 1000)
  
  
})



# Compare sequences
seqs_sub_4
dist_sub_4 

nchar(seqs_sub_4)

seqs_meta_final_sub_4$clone_id

library(pheatmap)
pheatmap(as.matrix(dist_sub_4),
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = viridis::viridis(100))



# plot
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# Your sequences vector (named for labeling)
seqs <- seqs_sub_4
names(seqs) <- paste0("Seq", seq_along(seqs))

# Split each sequence into individual characters
seq_matrix <- do.call(rbind, strsplit(seqs, ""))

# Find conserved positions (all sequences identical at that position)
is_conserved <- apply(seq_matrix, 2, function(col) length(unique(col)) == 1)

# Build long-format data frame for ggplot
df <- as.data.frame(seq_matrix) %>%
  mutate(seq_id = factor(names(seqs), levels = rev(names(seqs)))) %>%
  pivot_longer(-seq_id, names_to = "position", values_to = "base") %>%
  mutate(
    position = as.integer(sub("V", "", position)),
    conservation = ifelse(is_conserved[position], "conserved", "variable")
  ) %>% 
  mutate(fill_val = ifelse(conservation == "conserved", "conserved", base))

# --- alignment plot ---
p_align <- ggplot(df, aes(x = position, y = seq_id, fill = fill_val)) +
  geom_tile() +
  scale_fill_manual(values = c(
    conserved = "#378ADD",
    A = "#4CAF50", T = "#FF9800", G = "#E63946", C = "#9C27B0"
  )) +
  labs(x = "Position", y = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6), panel.grid = element_blank(), legend.position = "top")

p_align

# Example metadata — replace with your actual values
metadata <- data.frame(
  seq_id = factor(paste0("Seq", 1:84), levels = rev(paste0("Seq", 1:84))),
  group   = seqs_meta_final_sub_4$clone_id_plot
)

# --- metadata strip plot ---
p_meta <- ggplot(metadata, aes(x = 1, y = seq_id, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = clone_colors) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "group", y = NULL, fill = "Group") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  )

p_meta + p_align + plot_layout(widths = c(1, 20))



