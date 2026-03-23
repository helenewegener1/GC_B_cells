library(scoper)
library(dplyr)

# Following this SCOPer tutorial: 
# https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

bcr_data <- readRDS("45_immcantation/out/rds/bcr_data_qc_annot.rds")
patients <- names(bcr_data)

# ------------------------------------------------------------------------------
# Identifying clones by sequence identity
# ------------------------------------------------------------------------------

# Clonal assignment using identical nucleotide sequences
seq_clones <- lapply(patients, function(HH){
  
  identicalClones(
    bcr_data[[HH]], 
    method="nt", 
    cell_id = "cell_id"
  )
  
}) %>% setNames(patients)

saveRDS(seq_clones, "45_immcantation/out/rds/seq_clones.rds")

# -------------------
# HH119
# -------------------

HH <- "HH119"

seq_clones[[HH]] %>% 
  count(clone_id, sort = TRUE)
  # count(clone_id, v_call, j_call, sort = TRUE)
  
# -------------------
# HH117
# -------------------

HH <- "HH117"

seq_clones[[HH]] %>% 
count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# ------------------------------------------------------------------------------
# Identifying clones by hierical clustering (already done in 03_changeO_flow)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Identifying clones by spectral clustering
# ------------------------------------------------------------------------------

# -------------------
# novj method: Groups clones only based on junction sequences
# -------------------

spec_clones_novj <- lapply(patients, function(HH){
  
  spectralClones(
    bcr_data[[HH]], 
    method="novj", 
    junction="junction",
    # junction="cdr3", # In the paper they said that using cdr3 instead of junction improved performance. 
    cell_id = "cell_id"
  )
  
}) %>% setNames(patients)

saveRDS(spec_clones_novj, "45_immcantation/out/rds/spec_clones_novj.rds")

# -------------------
# HH119 - novj method
# -------------------

HH <- "HH119"
spec_clones_novj[[HH]]@db %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
plot(spec_clones_novj[[HH]], binwidth=0.02)

top_clone <- spec_clones_novj[[HH]]@db %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_novj[[HH]]@db %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_novj[[HH]]@db %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# -------------------
# HH117 - novj method
# -------------------

# HH117 has larger clones than with hierical clustering.

HH <- "HH117"
spec_clones_novj[[HH]]@db %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
plot(spec_clones_novj[[HH]], binwidth=0.02)

top_clone <- spec_clones_novj[[HH]]@db %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_novj[[HH]]@db %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_novj[[HH]]@db %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)


# -------------------
# vj method: Groups clones based on junction sequences and SHM in V and J sequences
# -------------------

spec_clones_vj <- lapply(patients, function(HH){
  
  spectralClones(
    bcr_data[[HH]], 
    method="vj", 
    summarize_clones = TRUE, 
    cell_id = "cell_id"
  )
  
}) %>% setNames(patients)

saveRDS(spec_clones_vj, "45_immcantation/out/rds/spec_clones_vj.rds")

# -------------------
# HH119 - vj method
# -------------------

HH <- "HH119"
spec_clones_vj[[HH]]@db %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
plot(spec_clones_vj[[HH]], binwidth=0.02)

top_clone <- spec_clones_vj[[HH]]@db %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_vj[[HH]]@db %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_vj[[HH]]@db %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# -------------------
# HH117 - vj method
# -------------------

# HH117 has larger clones than with hierical clustering.

HH <- "HH117"
spec_clones_vj[[HH]]@db %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
plot(spec_clones_vj[[HH]], binwidth=0.02)

top_clone <- spec_clones_vj[[HH]]@db %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

spec_clones_vj[[HH]]@db %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
spec_clones_vj[[HH]]@db %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# ------------------------------------------------------------------------------
# Compare methods: novj VS vj 
# ------------------------------------------------------------------------------

# -------------------
# HH117 
# -------------------

HH <- "HH119"

HH_spec_clones_vj <- spec_clones_vj[[HH]]@db
HH_spec_clones_novj <- spec_clones_novj[[HH]]@db

# Not in the same order...
table(HH_spec_clones_vj$cell_id %in% HH_spec_clones_novj$cell_id) 
table(HH_spec_clones_vj$cell_id == HH_spec_clones_novj$cell_id) 
length(HH_spec_clones_vj$cell_id)
length(HH_spec_clones_novj$cell_id)

# reorder novj to match the order of vj
HH_spec_clones_novj <- HH_spec_clones_novj[match(HH_spec_clones_vj$cell_id_seurat, HH_spec_clones_novj$cell_id_seurat), ]

# verify they now match
table(HH_spec_clones_vj$cell_id_seurat == HH_spec_clones_novj$cell_id_seurat)

# Investigate clone overlap clone_ids

# get top N clones from each method
top_n <- 50

top_vj <- HH_spec_clones_vj %>% count(clone_id, sort = TRUE) %>% slice_head(n = top_n) %>% pull(clone_id)
top_novj <- HH_spec_clones_novj %>% count(clone_id, sort = TRUE) %>% slice_head(n = top_n) %>% pull(clone_id)

# build overlap table for top clones only
overlap <- HH_spec_clones_vj %>%
  select(cell_id_seurat, clone_vj = clone_id) %>%
  left_join(HH_spec_clones_novj %>% select(cell_id_seurat, clone_novj = clone_id),
            by = "cell_id_seurat") %>%
  filter(clone_vj %in% top_vj, clone_novj %in% top_novj) %>%
  count(clone_vj, clone_novj)

# heatmap
ggplot(overlap, aes(x = factor(clone_novj), y = factor(clone_vj), fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "noVJ clone", y = "VJ clone", fill = "# cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
