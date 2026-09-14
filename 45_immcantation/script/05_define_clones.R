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
library(dowser)
library(patchwork)
library(stringdist)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# resolve_LC_list <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")

all_combined_list <- readRDS("45_immcantation/out/rds/04_bcr_heavy_light.rds")

patients <- names(all_combined_list)

for (HH in patients){
  print(HH)
  print(table(all_combined_list[[HH]]$locus))
  cat("\n")
}

# ------------------------------------------------------------------------------
# Prep data
# ------------------------------------------------------------------------------

lapply(patients, function(HH){

  # HH <- "HH119"
  
  df <- all_combined_list[[HH]]
  
  # Clean V and J gene by removing allele information
  df <- df %>%
    mutate(
      v_call_no_allele = gsub("\\*\\d+", "", v_call),
      v_call_no_allele = sapply(strsplit(v_call_no_allele, ","), function(x) paste(unique(x), collapse = ",")),
      j_call_no_allele = gsub("\\*\\d+", "", j_call),
      j_call_no_allele = sapply(strsplit(j_call_no_allele, ","), function(x) paste(unique(x), collapse = ",")),
      sample_clean_fol = ifelse(!is.na(manual_ADT_ID), glue("{sample_clean}_{manual_ADT_ID}"), sample_clean)
    )
  
  # ------------------------------------------------------------------------------
  # Heavy chain clone definition: same V and J gene (ignore allele) and cdr3 length + 90% similarity
  # ------------------------------------------------------------------------------
  
  df_heavy <- df %>% filter(locus == "IGH")
  
  # Do the connected components with 90% similarity.
  res_90_similarity <- hierarchicalClones(
    df_heavy,
    threshold = 0.1,        # 1 - 0.9 = 10% dissimilarity = 90% similarity
    method = "nt",          # or "aa" for amino acid
    linkage = "single",     # single linkage = connected components
    junction = "junction",
    v_call = "v_call_no_allele",
    j_call = "j_call_no_allele",
    clone = "clone_id_90_similarity", # output column
    cell_id = "cell_id", # single-cell mode
    first = FALSE,          # use all ambiguous gene calls for matching
    summarize_clones = FALSE
  )
  
  # Resolving light chain clones with dowser
  
  # Add heavy chain clone to df (includes both heavy and light chain information)
  df_heavy_chain_clones <- res_90_similarity %>% select(cell_id, clone_id_90_similarity)
  df_clones <- df %>% left_join(df_heavy_chain_clones, by = "cell_id")
  
  # dowser to resolve light chain
  resolve_LC_final <- resolveLightChains(
    df_clones,
    clone = "clone_id_90_similarity",
    v_call = "v_call_no_allele",
    j_call = "j_call_no_allele"
  )
  
  # Update clone column names
  resolve_LC_final <- resolve_LC_final %>% dplyr::rename(
    "clone_subgroup_90_similarity" = clone_subgroup,
    "clone_subgroup_id_90_similarity" = clone_subgroup_id
  )
  
  # ------------------------------------------------------------------------------
  # Export 
  # ------------------------------------------------------------------------------
  
  table(resolve_LC_final$locus, resolve_LC_final$sample_clean_fol)
  
  saveRDS(resolve_LC_final, glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds"))


})

# ------------------------------------------------------------------------------
# Inspect top clones
# ------------------------------------------------------------------------------

for (HH in patients){
  
  # HH <- "HH119"
  df <- readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds"))
  
  print(HH)
  print(df %>% dplyr::count(clone_subgroup_id_90_similarity, sort = TRUE) )
  cat("\n")
  
}

# ------------------------------------------------------------------------------
# Inspect top clones
# ------------------------------------------------------------------------------


n_top <- 15

files <- list.files("45_immcantation/out/rds/")
patients <- grep("05_", files, value = TRUE) %>% str_split_i("_", 2)

for (HH in patients) {
  
  # HH <- "HH153"
  
  df <- readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
    filter(locus == "IGH")
  
  # Top N clones by number of sequences 
  df_meta <- df %>%
    filter(
      L1_annotation == "GC_B_cells"
    ) %>% 
    count(clone_subgroup_id_90_similarity, junction_length, sort = TRUE) %>%
    slice_head(n = n_top) %>%
    dplyr::rename(clone = clone_subgroup_id_90_similarity) %>% 
    mutate(
      label = glue("{clone} (n={n}, jl={junction_length})") %>% as.factor()
    ) 
  
  # All pairwise Levenshtein distances within each clone
  df_seqs <- df %>%
    filter(clone_subgroup_id_90_similarity %in% df_meta$clone) %>% 
    select(clone_subgroup_id_90_similarity, junction) %>% 
    distinct()
  
  dist_list <- list()
  
  for (cl in df_meta$clone) {
    
    # cl <- "3100_1"
    
    seqs <- df_seqs %>%
      filter(clone_subgroup_id_90_similarity == cl) %>%
      pull(junction)
    
    if (length(seqs) < 2) next
    
    dist_list[[cl]] <- tibble(
      clone = cl,
      dist  = as.vector(stringdistmatrix(seqs, method = "lv"))
    )
  }
  
  df_dist <- bind_rows(dist_list) %>%
    left_join(df_meta, by = "clone") %>%
    mutate(label = factor(label, levels = df_meta$label))
  
  df_dist %>%
    ggplot(aes(x = label, y = dist)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", width = 0.6) +
    geom_jitter(width = 0.25, height = 0.15, size = 0.6, alpha = 0.4) +
    scale_y_continuous(breaks = scales::breaks_width(1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Clone",
      y = "Levenshtein distance (N mutations)",
      title = glue("{HH}: pairwise junction distances within top {n_top} GC B clones"),
      subtitle = "Unique junctions only; sequences have whole numbers of mutations"
    )
  
  ggsave(
    glue("45_immcantation/plot/05_define_clones/{HH}_junction_distance.png"),
    width = 16, height = 8, dpi = 300
  )

}
