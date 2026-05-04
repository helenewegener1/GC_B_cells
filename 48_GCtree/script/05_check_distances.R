library(glue)
library(tidyverse)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")

patients <- names(resolve_LC_list)


# ------------------------------------------------------------------------------
# Test 
# ------------------------------------------------------------------------------

HH <- "HH117"
clone <- "578_1"


# Heavy chain
resolve_LC_list[[HH]] %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% 
  count(j_call, v_call, junction_length, sort = TRUE) 

# Light chain 
resolve_LC_list$HH117 %>% filter(clone_subgroup_id == "578_1" & locus != "IGH") %>% 
  count(j_call, v_call, junction_length, sort = TRUE)

# Get sequences
seqs <- resolve_LC_list$HH117 %>% filter(clone_subgroup_id == "578_1" & locus == "IGH") %>% pull(sequence_alignment)

length(seqs)
seqs %>% nchar() %>% range()

# ------------------------------------------------------------------------------
# Compute distances
# ------------------------------------------------------------------------------

# Get distances 
dist_mat <- stringdistmatrix(seqs, method = "hamming")
# dist_mat <- stringdistmatrix(seqs, method = "lv")


# Distributions of distances. 
dist_vec <- as.vector(as.matrix(dist_mat))
dist_vec <- dist_vec[dist_vec > 0]  # remove diagonal zeros
summary(dist_vec)
hist(dist_vec, breaks = 100, main = "Pairwise Hamming distances")


# Distance to nearest neighbour
dist_mat_full <- as.matrix(dist_mat)
diag(dist_mat_full) <- NA
dist_to_nearest <- apply(dist_mat_full, 1, min, na.rm = TRUE)

# Plot
ggplot(tibble(d = dist_to_nearest), aes(x = d)) +
  geom_histogram(bins = 50) +
  labs(x = "Distance to nearest neighbour", 
       title = "Hamming distance to nearest")


# ------------------------------------------------------------------------------
# Hierarchical Clustering
# ------------------------------------------------------------------------------

fit <- hclust(as.dist(dist_mat), method = "complete")

# Convert your hclust object to a 'phylo' object
tree_phylo <- as.phylo(fit)

length(tree_phylo$tip.label)
length(seq_names)
tree_phylo$tip.label <- seq_names

# Plot the tree using a 'fan' layout (best for high density)
# 'mapping' connects the tree tips to your metadata
ggtree(tree_phylo, layout="fan", size=0.2)  +
  geom_tippoint(aes(color = sample_clean_fol), size=0.5, alpha=0.8) +
  theme_tree2()











