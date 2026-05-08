library(alakazam)
library(scoper)
library(dplyr)
library(shazam)
library(dowser)
library(tidyverse)
library(glue)
library(patchwork)

# Newest versions
packageVersion("scoper")
packageVersion("alakazam")

# load the built-in human SHM targeting model
data(HH_S5F)  # or S5F, RS5NF depending on your preference

# Following this SCOPer tutorial: 
# https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

subset_nr <- "2"

subset <- readRDS(glue("46_sequence_driven_clustering/out/subset_{subset_nr}.rds"))

subset <- subset %>%
  mutate(
    sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
  )


# ------------------------------------------------------------------------------
# Distiance to nearest plot
# ------------------------------------------------------------------------------


# Calculate nearest-neighbor Hamming distance distribution
# Per sequence, what is the Hamming distance to the nearest-neighbor (1 distance per sequence)
dist_nearest <- distToNearest(
  subset,
  sequenceColumn = "junction",
  cellIdColumn="cell_id" # Important for invoking in run in single-cell mode
)

# Check how many NAs
dist_nearest$dist_nearest %>% is.na() %>% table()


# -------------------
# Manual
# -------------------

# Nearest-neighbor Hamming distance distribution/histogram
dist_nearest %>% 
  filter(!is.na(dist_nearest)) %>% 
  ggplot(aes(x = dist_nearest)) + 
  geom_histogram(color = "white", binwidth = 0.02) +
  labs(
    title = glue("Subset {subset_nr}: Nearest-neighbor Hamming distance distribution"), 
    x = "Hamming distance", 
    y = "Count"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme_bw()

ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}/nearest_neighbor_distance_histogram.png"), width = 11, height = 6.5)

# -------------------
# Automatic - Density
# -------------------
# find threshold for cloning automatically using density 
# density finds the valley between the two modes directly
threshold_output <- shazam::findThreshold(
  dist_nearest$dist_nearest,
  method = "density" 
)

threshold_density <- threshold_output@threshold
threshold_density 

plot(threshold_output, binwidth = 0.02, silent = TRUE) +
  theme(axis.title = element_text(size = 12)) + 
  plot_annotation(
    title = glue("Subset {subset_nr}: Nearest-neighbor Hamming distance distribution"), 
    subtitle = glue("Automatic density threshold: {threshold_density}")
  )

ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}/nearest_neighbor_distance_histogram_automatic_density_threshold.png"), width = 12, height = 6.5)


# -------------------
# Automatic - GMM 
# -------------------

# find threshold for cloning automatically
threshold_output <- shazam::findThreshold(
  dist_nearest$dist_nearest,
  method = "gmm",
  model = "gamma-norm",
  cutoff = "intersect",
  # spc = 0.995  # specificity
  # spc = 0.99 # slight improve # HH119 - 0.04334775
)

threshold_gmm <- threshold_output@threshold
threshold_gmm

plot(threshold_output, binwidth = 0.02, silent = TRUE) +
  theme(axis.title = element_text(size = 12)) +
  plot_annotation(
    title = glue("Subset {subset_nr}: Nearest-neighbor Hamming distance distribution"),
    # subtitle = glue("GMM, gamma-norm model, specificity = 0.995 threshold: {threshold_gmm}")
    subtitle = glue("GMM, gamma-norm model: {threshold_gmm}")
  )

ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}/nearest_neighbor_distance_histogram_automatic_gmm_threshold.png"), width = 12, height = 6.5)


# ------------------------------------------------------------------------------
# vj method: Groups clones based on junction sequences and SHM in V and J sequences
# ------------------------------------------------------------------------------

subset_clones_vj <- spectralClones(
  
  subset,
  method="vj",
  threshold = threshold_gmm,
  # threshold = threshold_density,
  germline = "germline_alignment_d_mask",
  cell_id = "cell_id",
  junction = "junction",
  first = FALSE,
  targeting_model = HH_S5F,
  summarize_clones = TRUE
  
)


plot(subset_clones_vj)

# ------------------------------------------------------------------------------
# Clones VS patient 
# ------------------------------------------------------------------------------

# data <- subset_clones_novj
data <- subset_clones_vj@db
# version <- "scoper_vj_density_threshold"
version <- "scoper_vj_gmm_threshold"
# version <- "scoper_vj_no_threshold"

# N clones total 
data$clone_id %>% unique() %>% length()

# N clones (> 20 cells )
min_cells <- 10
mask <- table(data$clone_id) >= min_cells
table(mask)
these_clones <- names(which(mask))

subset_clones_vj_mask <- data %>% filter(clone_id %in% these_clones)

subset_clones_vj_mask %>%
  count(clone_id, patient_id) %>%
  ggplot(aes(x = patient_id, y = clone_id, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white", size = 3) +
  # scale_fill_viridis_c() +
  theme_minimal() +
  labs(x = "Patient", y = "Clone", fill = "Count",
       title = glue("SCOPer clones across patients - subset {subset_nr}"),
       subtitle = version,
       caption = glue("min N cells in clone: {min_cells}"))

ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}/clones_vs_patients_{version}.png"), width = 8, height = 8.5)
  

# ------------------------------------------------------------------------------
# Define top clones vj
# ------------------------------------------------------------------------------

subset_clones_vj <- subset_clones_vj@db

top_clones <- subset_clones_vj %>%
  count(clone_id, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(clone_id)

# ------------------------------------------------------------------------------
# Wrangle data
# ------------------------------------------------------------------------------

subset_clones_vj <- subset_clones_vj %>% 
  mutate(
    v_call_subgroup = v_call %>% str_remove_all("\\*[0-9]+"),
    v_call_subgroup = sapply(strsplit(v_call_subgroup, ","), function(x) unique(trimws(x))) %>%
      sapply(paste, collapse = ","),
    
    j_call_subgroup = j_call %>% str_remove_all("\\*[0-9]+"),
    j_call_subgroup = sapply(strsplit(j_call_subgroup, ","), function(x) unique(trimws(x))) %>%
      sapply(paste, collapse = ",")
  ) %>% 
  group_by(clone_id) %>% 
  mutate(
    v_call_subgroup_clone = paste(unique(unlist(strsplit(v_call_subgroup, ","))), collapse = ","),
    j_call_subgroup_clone = paste(unique(unlist(strsplit(j_call_subgroup, ","))), collapse = ",")
  ) %>% 
  ungroup()

# ------------------------------------------------------------------------------
# Visualize top clones vj
# ------------------------------------------------------------------------------

source("10_broad_annotation/script/color_palette.R")

n_clones <- length(top_clones)

for (clone_nr in 1:length(top_clones)){
  
  # clone_nr <- 1
  clone <- top_clones[clone_nr]
  
  plot_df <- subset_clones_vj %>% 
    filter(clone_id == clone) %>% 
    mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
    add_count(sample_clean_fol, name = "Count")
  
  # Meta data
  n_cells <- plot_df %>% nrow()
  v_gene <- plot_df$v_call_subgroup_clone %>% unique()
  j_gene <- plot_df$j_call_subgroup_clone %>% unique()
  
  # Color by cell type
  plot_df %>%   
    ggplot(aes(y = sample_clean_fol, fill = L1_annotation)) +
    geom_bar() + 
    scale_fill_manual(values = L1_colors) + 
    geom_text(aes(x = Count, label = Count), 
              hjust = -0.2, color = "black",
              stat = "unique") +
    labs(
      title = glue("Subset {subset_nr}: Top {clone_nr} clone"),
      subtitle = glue("Clone ID: {clone}"),
      caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
      y = ""
    ) +
    theme_bw()
  
  ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}_clone_nr_{clone_nr}.png"), width = 15, height = 8.5)
  
}

# ------------------------------------------------------------------------------
# Strict clone definition: 100% similar juntion 
# ------------------------------------------------------------------------------

clone_version <- "clone_id_junction"

top_seqs <- subset_clones_vj %>% count(junction, sort = TRUE) %>% head(20) %>% pull(junction)

# calls
subset_clones_vj %>% 
  filter(junction %in% top_seqs) %>% 
  count(junction, patient_id, v_call, j_call)

# calls_subgroup
subset_clones_vj %>% 
  filter(junction %in% top_seqs) %>% 
  count(junction, patient_id, v_call_subgroup, j_call_subgroup)

# Define clone ids from junction 
junction_clones <- subset_clones_vj$junction %>% unique()
names(junction_clones) <- paste0("clone_", 1:length(junction_clones))

# Add strict clone definition 
subset_clones_vj <- subset_clones_vj %>% 
  mutate(
    clone_id_junction = names(junction_clones)[match(junction, junction_clones)]
  )

# N clones (> 20 cells )
min_cells <- 5
mask <- table(subset_clones_vj$clone_id_junction) >= min_cells
table(mask)
these_clones <- names(which(mask))

subset_clones_vj_mask <- subset_clones_vj %>% filter(clone_id_junction %in% these_clones)

subset_clones_vj_mask %>%
  count(clone_id_junction, patient_id) %>%
  ggplot(aes(x = patient_id, y = clone_id_junction, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white", size = 3) +
  # scale_fill_viridis_c() +
  theme_minimal() +
  labs(x = "Patient", y = "Clone", fill = "Count",
       title = glue("SCOPer clones across patients - subset {subset_nr} - {clone_version}"),
       caption = glue("min N cells in clone: {min_cells}"))

ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}/{clone_version}_vs_patients.png"), width = 8, height = 10)

# ------------------------------------------------------------------------------
# Strict clone definition: same V and J and Junction length + 100% similar juntion 
# ------------------------------------------------------------------------------

clone_version <- "clone_id_vj_junction"

# Define clone ids from junction and vj call
vj_junction_clones <- paste0(subset_clones_vj$v_call_subgroup, "_", subset_clones_vj$j_call_subgroup, "_", subset_clones_vj$junction)
names(vj_junction_clones) <- paste0("clone_", 1:length(vj_junction_clones))

# Add strict clone definition 
subset_clones_vj <- subset_clones_vj %>% 
  mutate(
    clone_id_vj_junction = names(junction_clones)[match(junction, junction_clones)]
  )

# N clones (> 20 cells )
min_cells <- 5
mask <- table(subset_clones_vj$clone_id_vj_junction) >= min_cells
table(mask)
these_clones <- names(which(mask))

subset_clones_vj_mask <- subset_clones_vj %>% filter(clone_id_vj_junction %in% these_clones)

subset_clones_vj_mask %>%
  count(clone_id_junction, patient_id) %>%
  ggplot(aes(x = patient_id, y = clone_id_junction, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white", size = 3) +
  # scale_fill_viridis_c() +
  theme_minimal() +
  labs(x = "Patient", y = "Clone", fill = "Count",
       title = glue("SCOPer clones across patients - subset {subset_nr} - {clone_version}"),
       subtitle = glue("min N cells in clone: {min_cells+1}"))

ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}/{clone_version}_vs_patients.png"), width = 8, height = 10)

# ------------------------------------------------------------------------------
# Strict clone definition: same V + J + junction length + 95% similar junction
# ------------------------------------------------------------------------------
library(stringdist)
library(igraph)


# Group by V, J, and junction length first (required before similarity comparison)
subset_clones_vj <- subset_clones_vj %>%
  mutate(
    vj_group = paste0(v_call_subgroup, "_", j_call_subgroup, "_", junction_length)
  )

for (threshold in c(0.5, 0.4, 0.35, 0.3, 0.2, 0.1)){
  
  # For each VJ group, cluster junctions at 95% similarity
  # threshold <- 0.40  # 5% distance = 95% similarity
  similarity <- (1 - threshold) * 100
  clone_version <- glue("clone_id_vj_junction_{similarity}")
  
  vj_groups <- unique(subset_clones_vj$vj_group)
  
  clone_assignments <- lapply(vj_groups, function(grp) {
    rows <- which(subset_clones_vj$vj_group == grp)
    juncs <- subset_clones_vj$junction[rows]
    
    if (length(juncs) == 1) {
      return(data.frame(row = rows, clone_id_vj_junction_95 = paste0(grp, "_c1")))
    }
    
    # Pairwise normalized edit distance
    dmat <- stringdistmatrix(juncs, juncs, method = "hamming") / nchar(juncs[1])
    
    # Build graph: connect sequences within threshold
    adj <- dmat <= threshold
    diag(adj) <- FALSE
    g <- graph_from_adjacency_matrix(adj, mode = "undirected")
    
    # Connected components = clones
    comps <- components(g)$membership
    
    data.frame(
      row = rows,
      clone_id_vj_junction_95 = paste0(grp, "_c", comps)
    )
  }) %>% bind_rows()
  
  # Assign back to dataframe
  subset_clones_vj$clone_id_vj_junction_95 <- NA
  subset_clones_vj$clone_id_vj_junction_95[clone_assignments$row] <- clone_assignments$clone_id_vj_junction_95
  
  # Filter to clones with > min_cells
  min_cells <- 15
  mask <- table(subset_clones_vj$clone_id_vj_junction_95) >= min_cells
  these_clones <- names(which(mask))
  subset_clones_vj_mask <- subset_clones_vj %>% 
    filter(clone_id_vj_junction_95 %in% these_clones)
  
  subset_clones_vj_mask %>%
    count(clone_id_vj_junction_95, patient_id) %>%
    ggplot(aes(x = patient_id, y = clone_id_vj_junction_95, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n), color = "white", size = 3) +
    theme_minimal() +
    labs(
      x = "Patient", y = "Clone", fill = "Count",
      title = glue("SCOPer clones across patients - subset {subset_nr}"),
      subtitle = glue("{clone_version}"),
      caption = glue("min N cells in clone: {min_cells}")
    )
  
  ggsave(glue("46_sequence_driven_clustering/plot/subset_{subset_nr}/{clone_version}_vs_patients.png"), 
         width = 8, height = 10)
  
  
}

# Save 
saveRDS(subset_clones_vj, glue("46_sequence_driven_clustering/out/subset_{subset_nr}_clones_vj.rds"))

