library(tidyverse)
library(stringdist)
library(stats)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(fastcluster)
library(glue)
library(patchwork)

version <- "subset_2"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

subset_clones_vj <- readRDS(glue("46_sequence_driven_clustering/out/{version}_clones_vj.rds"))

# ------------------------------------------------------------------------------
# Compare sequences - SCOPer
# ------------------------------------------------------------------------------

# Get clone
clone <- subset_clones_vj %>% count(clone_id, sort = TRUE) %>% head(1) %>% pull(clone_id)

subset_clones_vj %>% filter(clone_id == clone) %>% count(clone_id, v_call, j_call, sort = TRUE)

# Pull sequences from that clone
seqs <- subset_clones_vj %>% filter(clone_id == clone) %>% arrange(sequence_alignment) %>% pull(sequence_alignment)
# subset_clones_vj %>% filter(clone_id == clone) %>% pull(germline_alignment_d_mask) %>% unique()

# Your sequences vector (named for labeling)
names(seqs) <- paste0("Seq", seq_along(seqs))

# Get clone meta data
meta_data_clone <- subset_clones_vj %>% filter(sequence_alignment %in% seqs) %>% arrange(sequence_alignment)
meta_data_clone$seq_id <- paste0("Seq", seq_along(seqs))
# Check order
table(meta_data_clone$sequence_alignment == seqs)

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
    conserved = "grey95",
    A = "#4CAF50", T = "#FF9800", G = "#E63946", C = "#9C27B0"
  )) +
  labs(x = "Position", y = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6), panel.grid = element_blank(), legend.position = "top")

p_align

# --- metadata strip plot ---
var <- "v_call"
p_meta <- meta_data_clone %>%
  mutate(seq_id = factor(seq_id, levels = rev(seq_id))) %>%
  ggplot(aes(x = 1, y = seq_id, fill = !!sym(var))) +
  geom_tile() +
  # scale_fill_manual(values = clone_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    x = var, y = NULL, fill = "Group",
    title = glue("SCOPer clone {clone}. N cells: {length(seqs)}")
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "left"
  )

p_meta + p_align + plot_layout(widths = c(1, 20))

ggsave(glue("46_sequence_driven_clustering/plot/{version}/compare_seqs/compare_seqs_clone_{clone}_{var}.png"), width = 18, height = 8)


