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

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# resolve_LC_list <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")

resolve_LC_list <- readRDS("45_immcantation/out/rds/04_bcr_heavy_light.rds")

table(resolve_LC_list$HH117$locus)
table(resolve_LC_list$HH119$locus)

patients <- names(resolve_LC_list)

# ------------------------------------------------------------------------------
# Prep data
# ------------------------------------------------------------------------------

# lapply(patients, function(HH){
#   
HH <- "HH119"
#   
# })

df <- resolve_LC_list[[HH]]

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

