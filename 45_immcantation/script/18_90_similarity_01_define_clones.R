library(tidyverse)
library(stringdist)
library(stats)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(fastcluster)
library(glue)
library(stringdist)
library(igraph)
library(scoper)
library(pheatmap)

# version <- "HH119_two_largest_clones"
# version <- "HH119_two_largest_clones_junction"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list_gmm_threshold_germlined.rds")
# table(resolve_LC_list$HH117$locus)
# table(resolve_LC_list$HH119$locus)

# ------------------------------------------------------------------------------
# Prep data
# ------------------------------------------------------------------------------

HH <- "HH119"
# df <- resolve_LC_list[[HH]]

# # Clean V and J gene by removing allele information 
# df <- df %>%
#   mutate(
#     v_call_no_allele = gsub("\\*\\d+", "", v_call),
#     v_call_no_allele = sapply(strsplit(v_call_no_allele, ","), function(x) paste(unique(x), collapse = ",")),
#     j_call_no_allele = gsub("\\*\\d+", "", j_call),
#     j_call_no_allele = sapply(strsplit(j_call_no_allele, ","), function(x) paste(unique(x), collapse = ","))
#   ) %>% 
#   rename(
#     "clone_id_og_scoper" = clone_id, 
#     "clone_subgroup_og_scoper" = clone_subgroup,
#     "clone_subgroup_id_og_scoper" = clone_subgroup_id
#   )
# 
# # ------------------------------------------------------------------------------
# # SCOPer no allele 
# # ------------------------------------------------------------------------------
# 
# # Threshold
# list_thresholds <- readRDS("45_immcantation/out/rds/04_list_thresholds.rds")
# list_thresholds[[HH]]$gmm
# 
# # Subset heavy
# df_heavy <- df %>% filter(locus == "IGH")
# 
# # SCOPer with no allele information in V and J gene
# res_scoper_no_allele <- spectralClones(
#   
#   df_heavy,
#   method="vj", 
#   v_call = "v_call_no_allele", 
#   j_call = "j_call_no_allele",
#   threshold = list_thresholds[[HH]]$gmm,
#   germline = "germline_alignment_d_mask",
#   cell_id = "cell_id",
#   clone = "clone_id_no_allele_scoper", # output column 
#   junction = "junction",
#   first = FALSE,
#   targeting_model = HH_S5F,
#   summarize_clones = FALSE
#   
# )
# 
# # Add heavy chain clone to df (includes both heavy and light chain information)
# df_heavy_chain_clones <- res_scoper_no_allele %>% select(cell_id, clone_id_no_allele_scoper)
# df_clones <- df %>% left_join(df_heavy_chain_clones, by = "cell_id")
# 
# # dowser to resolve light chain
# res_scoper_no_allele_resolve_LC <- resolveLightChains(
#   df_clones,
#   clone = "clone_id_no_allele_scoper", 
#   v_call = "v_call_no_allele",
#   j_call = "j_call_no_allele"
# )
# 
# # Update clone column names
# res_scoper_no_allele_resolve_LC <- res_scoper_no_allele_resolve_LC %>% rename(
#   "clone_subgroup_no_allele_scoper" = clone_subgroup,
#   "clone_subgroup_id_no_allele_scoper" = clone_subgroup_id
# )
# 
# # ------------------------------------------------------------------------------
# # Heavy chain clone definition: same V and J gene (ignore allele) and cdr3 length + 90% similarity
# # ------------------------------------------------------------------------------
# 
# df_heavy <- res_scoper_no_allele_resolve_LC %>% filter(locus == "IGH")
# 
# # Do the connected components with 90% similarity. 
# res_90_similarity <- hierarchicalClones(
#   df_heavy,
#   threshold = 0.1,        # 1 - 0.9 = 10% dissimilarity = 90% similarity
#   method = "nt",          # or "aa" for amino acid
#   linkage = "single",     # single linkage = connected components
#   junction = "junction", 
#   v_call = "v_call_no_allele",
#   j_call = "j_call_no_allele",
#   clone = "clone_id_90_similarity", # output column 
#   cell_id = "cell_id", # single-cell mode
#   first = FALSE,          # use all ambiguous gene calls for matching
#   summarize_clones = FALSE
# )
# 
# # Resolving light chain clones with dowser
# 
# # Add heavy chain clone to df (includes both heavy and light chain information)
# df_heavy_chain_clones <- res_90_similarity %>% select(cell_id, clone_id_90_similarity)
# df_clones <- res_scoper_no_allele_resolve_LC %>% left_join(df_heavy_chain_clones, by = "cell_id")
# 
# # dowser to resolve light chain 
# resolve_LC_final <- resolveLightChains(
#   df_clones,
#   clone = "clone_id_90_similarity",
#   v_call = "v_call_no_allele",
#   j_call = "j_call_no_allele"
# )
# 
# # Update clone column names
# resolve_LC_final <- resolve_LC_final %>% rename(
#   "clone_subgroup_90_similarity" = clone_subgroup,
#   "clone_subgroup_id_90_similarity" = clone_subgroup_id
# )
# 
# saveRDS(resolve_LC_final, glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds"))
resolve_LC_final <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds"))

# ------------------------------------------------------------------------------
# Heatmaps to compare clone definitions
# ------------------------------------------------------------------------------

df_heavy <- resolve_LC_final %>% filter(locus == "IGH")

# Top clones subgroups
n_clones <- 20
clones_og_scoper <- df_heavy %>% count(clone_subgroup_id_og_scoper, sort = TRUE) %>% head(n_clones) %>% pull(clone_subgroup_id_og_scoper)
clones_no_allele_scoper <- df_heavy %>% count(clone_subgroup_id_no_allele_scoper, sort = TRUE) %>% head(n_clones) %>% pull(clone_subgroup_id_no_allele_scoper)
clones_90_similarity <- df_heavy %>% count(clone_subgroup_id_90_similarity, sort = TRUE) %>% head(n_clones) %>% pull(clone_subgroup_id_90_similarity)

# Prep for plotting
resolve_LC_compare <- df_heavy %>% 
  filter(
    (clone_subgroup_id_og_scoper %in% clones_og_scoper | 
      clone_subgroup_id_no_allele_scoper %in% clones_no_allele_scoper | 
      clone_subgroup_id_90_similarity %in% clones_90_similarity)
  ) 

# Define clone definitions 
clone_definitions <- list(
  "SCOPer (original)" = "clone_subgroup_id_og_scoper", 
  "SCOPer (no allele information)" = "clone_subgroup_id_no_allele_scoper", 
  "90% similarity + no allele information" = "clone_subgroup_id_90_similarity"
)

pairs <- combn(length(clone_definitions), 2, simplify = FALSE)

for (pair in pairs) {
  
  i <- pair[1]
  j <- pair[2]
  
  col_i <- clone_definitions[[i]]
  col_j <- clone_definitions[[j]]
  path_i <- str_remove(col_i, "clone_subgroup_id_")
  path_j <- str_remove(col_j, "clone_subgroup_id_")
  
  resolve_LC_compare %>% 
    dplyr::select(!!sym(col_i), !!sym(col_j)) %>% 
    table() %>% 
    as.matrix() %>%
    as.data.frame() %>%
    ggplot(aes(x = !!sym(col_i), y = !!sym(col_j), fill = Freq)) +
    geom_tile(color = "grey50") +
    geom_text(aes(label = Freq), size = 3, color = "grey80") +
    scale_fill_gradientn(
      colors = c("grey80", "blue", "darkblue"),
      values = scales::rescale(c(0, 0.001, 1))
    ) +
    labs(
      title = glue("{HH}: {names(clone_definitions)[[i]]} VS {names(clone_definitions)[[j]]}"),
      x = names(clone_definitions)[[i]], 
      y = names(clone_definitions)[[j]], 
      fill = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(glue("45_immcantation/plot/18_90_similarity/{HH}_compare_{path_i}_vs_{path_j}.png"), width = 15, height = 6)
  
}

# ------------------------------------------------------------------------------
# TREES of heavy chain - junction 
# ------------------------------------------------------------------------------

# Subset heavy chain 
df_heavy <- resolve_LC_final %>% filter(locus == "IGH") 

# Top clones subgroups
n_clones <- 16
clone_colors_vec <- c(
  "#E63946", "#2196F3", "#3DAA55", "#FF9800",
  "#9C27B0", "#00BCD4", "#F5C518", "#FF4081",
  "#6D4C41", "#76FF03", "#6200EA", "#37474F",
  "#80CBC4", "#CE93D8", "#FFF176", "#A5D6A7",
  "grey90"
)

clone_definitions <- c(
  "clone_subgroup_id_og_scoper",
  "clone_subgroup_id_no_allele_scoper",
  "clone_subgroup_id_90_similarity"
)

top_clones <- lapply(clone_definitions, function(col) {
  df_heavy %>% count(!!sym(col), sort = TRUE) %>% head(n_clones) %>% pull(1)
})

clone_colors <- setNames(
  lapply(top_clones, function(clones) setNames(clone_colors_vec, c(clones, "other"))),
  clone_definitions
)

# Define subset to make tree of
subset <- df_heavy %>% 
  filter(
    clone_subgroup_id_og_scoper %in% clones_og_scoper | 
      clone_subgroup_id_no_allele_scoper %in% clones_no_allele_scoper | 
      clone_subgroup_id_90_similarity %in% clones_90_similarity
  ) 

# N_cells
nrow(subset)

# Get largest clone
top_1_clones <- lapply(top_clones, function(x) x[[1]] ) %>% setNames(clone_definitions)

# Subset heavy chain and plasma cells and trim sequence 
subset <- subset %>% 
  mutate(
    sequence_trimmed = str_sub(sequence, v_sequence_start, j_sequence_end),
    sequence_trimmed_300 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-299, nchar(sequence_trimmed))
  )

# Have a look 
subset$sequence_trimmed %>% nchar() %>% range()
subset$sequence_trimmed %>% nchar() %>% hist()

subset$sequence_trimmed_300 %>% nchar() %>% unique()

# Extract sequences
sequence_types <- c("junction", "sequence_trimmed_300")

for (sequence_type in sequence_types){
  
  # sequence_type <- "junction"
  
  # seqs <- subset$sequence_trimmed_300
  seqs <- subset[[sequence_type]]
  seq_names <- subset$sequence_id
  
  seqs %>% nchar() %>% range()
  
  # Get metadata of sequences
  # clones_subgroup <- resolve_LC %>% count(clone_subgroup_id, sort = TRUE) %>% head(n_clones) %>% pull(clone_subgroup_id)
  # clones_subgroup_scoper <- resolve_LC %>% count(clone_subgroup_id_scoper, sort = TRUE) %>% head(n_clones) %>% pull(clone_subgroup_id_scoper)
  
  seqs_meta <- subset %>% 
    mutate(
      label = sequence_id, 
      junction_length = as.character(junction_length),
      clone_subgroup_id_og_scoper_plot = ifelse(clone_subgroup_id_og_scoper %in% names(clone_colors$clone_subgroup_id_og_scoper), clone_subgroup_id_og_scoper, "other"),
      clone_subgroup_id_no_allele_scoper_plot = ifelse(clone_subgroup_id_no_allele_scoper %in% names(clone_colors$clone_subgroup_id_no_allele_scoper), clone_subgroup_id_no_allele_scoper, "other"),
      clone_subgroup_id_90_similarity_plot = ifelse(clone_subgroup_id_90_similarity %in% names(clone_colors$clone_subgroup_id_90_similarity), clone_subgroup_id_90_similarity, "other")
    ) %>% 
    select(label, everything()) %>% # Move label to front
    as.data.frame()
  
  # ------------------------------------------------------------------------------
  # Compute Levenshtein as sequences are not of the same length
  # ------------------------------------------------------------------------------
  
  # Get distances 
  if (sequence_type == "junction"){
    dist_mat <- stringdistmatrix(seqs, method = "lv") # different length distance measurement 
  } else if (sequence_type == "sequence_trimmed_300"){
    dist_mat <- stringdistmatrix(seqs, method = "hamming") # same length distance measurement 
  }

  # ------------------------------------------------------------------------------
  # Hierarchical Clustering
  # ------------------------------------------------------------------------------
  
  fit <- hclust(as.dist(dist_mat), method = "complete")
  
  # saveRDS(fit, glue("46_sequence_driven_clustering/out/fit_{version}_{version}.rds"))
  
  # ==============================================================================
  # no-zoom tree
  # ==============================================================================
  
  outdir_trees <- glue("45_immcantation/plot/18_90_similarity/trees_{HH}_{sequence_type}/")
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
    "v_call_no_allele",  "j_call_no_allele", "v_call",  "j_call", "junction_length",
    "clone_subgroup_id_og_scoper_plot", "clone_subgroup_id_no_allele_scoper_plot", "clone_subgroup_id_90_similarity_plot"
  )
  
  for (var in variables){
    
    # var <- "clone_subgroup_id_og_scoper_plot"
    
    p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta + 
      geom_tippoint(aes(color = !!sym(var)), size=0.5) +
      theme_tree2() + 
      guides(color = guide_legend(override.aes = list(size = 4))) + 
      labs(
        title = glue("{HH} {sequence_type} colored by {var}")
        # subtitle = glue("Clusters: {clusters_string}")
      )
    
    
    if (str_detect(var, "clone")){
      var_color <- str_remove(var, "_plot") 
      p <- p + scale_color_manual(values = clone_colors[[var_color]])
    }
    
    ggsave(glue("{outdir_trees}/{HH}_{var}.png"), plot = p, width = 15, height = 10, dpi = 1000)
    
  }
  
  # ------------------------------------------------------------------------------
  # Remove largest clone
  # ------------------------------------------------------------------------------
  
  for (clone_def in clone_definitions){
    
    # clone_def <- "clone_subgroup_id_og_scoper"
    
    # Define meta data
    clone_def_plot <- paste0(clone_def, "_plot")
    top_1_clone <- top_1_clones[[clone_def]]
    seqs_meta_subset <- seqs_meta %>% filter(!!sym(clone_def_plot) != top_1_clone)
    nrow(seqs_meta_subset)
    
    # Get sequences
    seqs <- seqs_meta_subset[[sequence_type]]
    seq_names <- seqs_meta_subset$sequence_id
    length(seqs)
    
    # Get distances 
    if (sequence_type == "junction"){
      dist_mat <- stringdistmatrix(seqs, method = "lv") # different length distance measurement 
    } else if (sequence_type == "sequence_trimmed_300"){
      dist_mat <- stringdistmatrix(seqs, method = "hamming") # same length distance measurement 
    }
    
    # Hierarchical Clustering
    fit <- hclust(as.dist(dist_mat), method = "complete")
    
    # Convert to phylo for nicer plotting
    tree <- as.phylo(fit)
    
    # Checks 
    length(tree$tip.label)
    length(seqs_meta_subset$sequence_id)
    tree$tip.label <- seqs_meta_subset$sequence_id
    
    # variables <- c(
    #   "v_call_no_allele",  "j_call_no_allele", "junction_length", clone_def_plot
    # )
    # 
    # for (var in variables){
    #   
    var <- clone_def_plot
    
    p <- ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta_subset + 
      geom_tippoint(aes(color = !!sym(var)), size=0.5) +
      theme_tree2() + 
      guides(color = guide_legend(override.aes = list(size = 4))) + 
      labs(
        title = glue("{HH} {sequence_type} colored by {var} - largest clone removed")
        # subtitle = glue("Clusters: {clusters_string}")
      )
    
    if (str_detect(var, "clone")){
      var_color <- str_remove(var, "_plot") 
      p <- p + scale_color_manual(values = clone_colors[[var_color]][c(2:length(clone_colors[[var_color]]))])
    }
    
    outdir_trees <- glue("45_immcantation/plot/18_90_similarity/trees_{HH}_{sequence_type}/")
    ggsave(glue("{outdir_trees}/{HH}_rm_largest_{var}.png"), plot = p, width = 15, height = 10, dpi = 1000)
      
    # }
      
  }
  
}


# # ------------------------------------------------------------------------------
# # Cut the tree
# # ------------------------------------------------------------------------------
# 
# k <- 8
# 
# # outdir <- glue("46_sequence_driven_clustering/plot/{version}/{k}_clusters")
# # dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# 
# # Cut tree
# clusters <- cutree(fit, k = k)
# table(clusters)
# 
# # Add to metadata
# seqs_meta$clusters <- clusters %>% as.character()
# 
# # See patient split
# table(seqs_meta$clusters, seqs_meta$clone_subgroup_id)
# 
# # Show clusters 
# var <- "clusters"
# ggtree(tree, layout="fan", size=0.2) %<+% seqs_meta + 
#   geom_tippoint(aes(color = !!sym(var)), size=0.5, alpha = 0.5) +
#   theme_tree2() + 
#   guides(color = guide_legend(override.aes = list(size = 4))) + 
#   labs(
#     title = glue("{version} colored by {k}_{var}")
#     # subtitle = glue("Clusters: {clusters_string}")
#   )
# 
# ggsave(glue("{outdir_trees}/{version}_{k}_{var}.png"), width = 15, height = 10, dpi = 1000)
# 
# 
# 
