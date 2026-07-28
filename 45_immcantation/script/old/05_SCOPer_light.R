library(alakazam)
library(scoper)
library(dplyr)
library(shazam)
library(dowser)
library(tidyverse)
library(glue)

# Newest versions
packageVersion("scoper")
packageVersion("alakazam")

# Following this SCOPer tutorial: 
# https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# list_thresholds <- readRDS("45_immcantation/out/rds/list_thresholds.rds")

bcr_data <- readRDS("45_immcantation/out/rds/03_light_bcr_data_qc_annot.rds")

patients <- names(bcr_data)

bcr_data <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  bcr_data_HH <- bcr_data[[HH]]
  
  bcr_data_HH <- bcr_data_HH %>%
    mutate(
      sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
    )
  
  return(bcr_data_HH)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Identifying clones by sequence identity
# ------------------------------------------------------------------------------
# 
# # Clonal assignment using identical nucleotide sequences
# seq_clones <- lapply(patients, function(HH){
# 
#   identicalClones(
#     bcr_data[[HH]],
#     method="nt",
#     cell_id = "cell_id",
#     junction = "junction",
#     first = FALSE
#   )
# 
# }) %>% setNames(patients)
# 
# saveRDS(seq_clones, "45_immcantation/out/rds/seq_clones.rds")
# 
# # -------------------
# # HH119
# # -------------------
# 
# HH <- "HH119"
# 
# seq_clones[[HH]] %>% 
#   count(clone_id, sort = TRUE)
#   # count(clone_id, v_call, j_call, sort = TRUE)
#   
# # -------------------
# # HH117
# # -------------------
# 
# HH <- "HH117"
# 
# seq_clones[[HH]] %>% 
# count(clone_id, sort = TRUE)
# # count(clone_id, v_call, j_call, sort = TRUE)

# ------------------------------------------------------------------------------
# Identifying clones by spectral clustering
# novj method: Groups clones only based on junction sequences
# ------------------------------------------------------------------------------

# Lambda chain 
chain <- "IGL"
spec_clones_novj_IGL <- lapply(patients, function(HH){

  # HH <- "HH117"
  data <- bcr_data[[HH]] %>% filter(locus == chain)
  data$locus <- "IGH" # Work around to run spectralClones on light chain
  
  df <- spectralClones(
    data,
    method="novj",
    germline = "germline_alignment_d_mask",
    junction="junction", # In the paper they said that using cdr3 instead of junction improved performance. Mats say that the longer the sequence the better
    cell_id = "cell_id",
    first = FALSE
  )
  
  df$locus <- chain
  
  return(df)

}) %>% setNames(patients)


# Kappa chain 
chain <- "IGK"
spec_clones_novj_IGK <- lapply(patients, function(HH){
  
  # HH <- "HH117"
  data <- bcr_data[[HH]] %>% filter(locus == chain)
  data$locus <- "IGH"
  
  df <- spectralClones(
    data,
    method="novj",
    germline = "germline_alignment_d_mask",
    junction="junction", # In the paper they said that using cdr3 instead of junction improved performance. Mats say that the longer the sequence the better
    cell_id = "cell_id",
    first = FALSE
  )
  
  df$locus <- chain
  
  return(df)
  
}) %>% setNames(patients)

# Combine per patient
spec_clones_novj_light <- lapply(patients, function(HH) {
  
  bind_rows(spec_clones_novj_IGK[[HH]], spec_clones_novj_IGL[[HH]])
  
}) %>% setNames(patients)

# -------------------
# Detect main V and J gene for each clone
# -------------------

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

spec_clones_novj_light <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  spec_clones_novj_light[[HH]] %>%
    group_by(clone_id) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)

# spec_clones_novj_light$HH117 %>% count(clone_id, v_call, j_call, v_call_majority, j_call_majority, sort = TRUE)

saveRDS(spec_clones_novj_light, "45_immcantation/out/rds/05_spec_clones_novj_light.rds")

# -------------------
# HH119 - novj method
# -------------------

HH <- "HH119"
spec_clones_novj_light[[HH]] %>%
  count(clone_id, v_call, j_call, sort = TRUE)

top_clone <- spec_clones_novj_light[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>%
  head(1) %>% pull(clone_id)

spec_clones_novj_light[[HH]] %>%
  filter(clone_id == top_clone) %>%
  count(sample_clean, v_call, j_call, sort = TRUE)

# -------------------
# HH117 - novj method
# -------------------

# HH117 has larger clones than with hierical clustering.

HH <- "HH117"
spec_clones_novj_light[[HH]] %>%
  count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
plot(spec_clones_novj[[HH]], binwidth=0.02)

top_clone <- spec_clones_novj[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>%
  head(1) %>% pull(clone_id)

spec_clones_novj[[HH]] %>%
  filter(clone_id == top_clone) %>%
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_novj[[HH]] %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# -------------------
# Define top clones novj
# -------------------

top_GC_clones_novj <- lapply(patients, function(HH) {

  # find clones that contain at least 1 GC cell
  GC_clones <- spec_clones_novj[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    pull(clone_id) %>%
    unique()

  # rank those clones by total size (all cell types) and take top 10
  spec_clones_novj[[HH]] %>%
    filter(clone_id %in% GC_clones) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_id)

}) %>% setNames(patients)

# -------------------
# Look at top clones novj
# -------------------

lapply(patients, function(HH){

  # HH <- "HH119"
  # HH <- "HH117"

  spec_clones_novj[[HH]] %>%
    filter(clone_id %in% top_GC_clones_novj[[HH]]) %>%
    mutate(
      v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
    ) %>%
    mutate(
      j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
    ) %>%
    count(clone_id, v_gene, j_gene, sort = TRUE)

}) %>% setNames(patients)

