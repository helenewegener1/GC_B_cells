library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
library(reshape2)
library(ggtree)
library(treeio)
library(stringdist)
library(igraph)
library(dowser)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC_3_definitions", rds_files, value = TRUE)
patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
patients

resolve_LC_list <- lapply(resolve_LC_files, function(x) readRDS(glue("45_immcantation/out/rds/{x}"))) %>% 
  setNames(patients)

# Look at clone IDs
grep("clone", colnames(resolve_LC_list$HH117), value = TRUE)

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

resolve_LC_list <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  resolve_LC_list[[HH]] %>%
    group_by(clone_subgroup_id_90_similarity, locus) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)


# ------------------------------------------------------------------------------
# Make germline
# ------------------------------------------------------------------------------

# resolve_LC_germline_list <- list()
# 
# for (HH in patients){
#   
#   # HH <- "HH119"
#   
#   resolve_LC_HH <- resolve_LC_list[[HH]]
#   
#   # Read in IMGT-gapped sequences
#   # references <- readIMGT(dir = "../packages/immcantation/scripts/germlines/human/vdj")
#   references <- readIMGT(dir = "00_data/vdj")
#   
#   # remove germline alignment columns for this example
#   db <- select(resolve_LC_HH, -"germline_alignment", -"germline_alignment_d_mask")
#   
#   # Reconstruct germline sequences
#   resolve_LC_HH_germline <- createGermlines(db, references, nproc=1, clone = "clone_subgroup_id_90_similarity")
#   
#   # append to list 
#   resolve_LC_germline_list[[HH]] <- resolve_LC_HH_germline
#   
# }
# 
# saveRDS(resolve_LC_germline_list, "45_immcantation/out/rds/resolve_LC_90_similarity_germlined.rds")
resolve_LC_germline_list <- readRDS("45_immcantation/out/rds/resolve_LC_90_similarity_germlined.rds")

resolve_LC_HH_germline <- resolve_LC_germline_list[[HH]]

# Make outdir
outdir <- glue("45_immcantation/plot/18_90_similarity/06_clone_mutations/{HH}")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Compare sequences - SCOPer
# ------------------------------------------------------------------------------

df_heavy <- resolve_LC_HH_germline %>% filter(locus == "IGH")
df_light <- resolve_LC_HH_germline %>% filter(locus != "IGH")

clone_nr <- 1

# Get clone
clone <- df_heavy %>% 
  filter(L1_annotation == "GC_B_cells") %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  slice(clone_nr) %>% pull(clone_subgroup_id_90_similarity)

# Look at clone
df_heavy %>% filter(clone_subgroup_id_90_similarity == clone) %>% count(clone_subgroup_id_90_similarity, v_call_no_allele, j_call_no_allele, sort = TRUE)

# Pull sequences from that clone
seqs <- df_heavy %>% filter(clone_subgroup_id_90_similarity == clone) %>% arrange(sequence_alignment) %>% pull(sequence_alignment)

length(seqs)

# Get clone meta data
df_clone <- df_heavy %>% filter(sequence_alignment %in% seqs & clone_subgroup_id_90_similarity == clone) %>% arrange(sequence_alignment)
jl <- df_clone$junction_length %>% unique()

# Your sequences vector (named for labeling)
names(seqs) <- paste0("Seq", seq_along(seqs))

df_clone$seq_id <- paste0("Seq", seq_along(seqs))
# Check order
# table(df_clone$sequence_alignment == seqs)

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
ggplot(df, aes(x = position, y = seq_id, fill = fill_val)) +
  geom_tile() +
  scale_fill_manual(values = c(
    conserved = "grey95",
    A = "#4CAF50", T = "#FF9800", G = "#E63946", C = "#9C27B0"
  )) +
  labs(
    title = glue("{HH}: Clone: {clone}"), 
    subtitle = glue("Junction length: {jl}"),
    x = "Position", y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6), panel.grid = element_blank(), legend.position = "top")

ggsave(glue("{outdir}/{clone}_compare_seqs.png"), width = 18, height = 8)


# ------------------------------------------------------------------------------
# Trees
# ------------------------------------------------------------------------------

# Define trimmmed sequences 
df_clone <- df_clone %>% 
  mutate(
    sequence_trimmed = str_sub(sequence, v_sequence_start, j_sequence_end),
    sequence_trimmed_300 = str_sub(sequence_trimmed, nchar(sequence_trimmed)-299, nchar(sequence_trimmed))
  )

# sequence_alignment

# Extract sequences
# sequence_types <- c("junction", "sequence_trimmed_300")

# for (sequence_type in sequence_types){

# sequence_type <- "sequence_trimmed_300"
sequence_type <- "junction"

seqs <- df_clone[[sequence_type]]
seq_names <- df_clone$seq_id

seqs %>% nchar() %>% range()

# Compute Levenshtein as sequences are not of the same length
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
length(df_clone$seq_id)
tree$tip.label <- df_clone$seq_id
df_clone$sequence_id <- df_clone$seq_id


# Visualize clustering as tree
vars <- c("L1_annotation", "c_call", "manual_ADT_ID")

for (var in vars){
  
  ggtree(tree, layout="fan", size=0.2) %<+% df_clone + 
    geom_tippoint(aes(color = !!sym(var)), size=1) +
    geom_tiplab() +
    theme_tree2(linesize = 0.2) + 
    guides(color = guide_legend(override.aes = list(size = 4))) + 
    labs(
      title = glue("{HH}, {clone}, {sequence_type}"),
      subtitle = glue("Junction length: {jl}")
      # subtitle = glue("Clusters: {clusters_string}")
    )
  
  ggsave(glue("{outdir}/{clone}_tree_{sequence_type}_{var}.png"), width = 15, height = 10, dpi = 1000)
  
}
 
# ------------------------------------------------------------------------------
# N mutations
# ------------------------------------------------------------------------------

# Dendogram
attr(dist_mat, "Labels") <- seq_names

png(glue("{outdir}/{clone}_dendogram_{sequence_type}.png"), width = 8000, height = 4000, res = 300)
hclust(dist_mat, method = "complete") %>% 
  plot(
    main = glue("{HH}, {clone}, {sequence_type}, junction length: {jl}"),
    xlab = "", sub = "",
    ylab = "N mutations"
  )
dev.off()


# Heatmap 
mat <- as.matrix(dist_mat)
png(glue("{outdir}/{clone}_heatmap_{sequence_type}.png"), width = 3000, height = 3000, res = 300)
melt(mat) |>
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), size = 2) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(x = "", y = "", fill = "N mutations",
       title = glue("{HH}, {clone}, {sequence_type}, junction length: {jl}"))
dev.off()

# ------------------------------------------------------------------------------
# N mutations across top 10 clones
# ------------------------------------------------------------------------------

df_heavy <- resolve_LC_HH_germline %>% filter(locus == "IGH")
df_light <- resolve_LC_HH_germline %>% filter(locus != "IGH")

# ------------------------------------------------------------------------------
# N mutations across top 15 clones (boxplot)
# ------------------------------------------------------------------------------

df_heavy <- resolve_LC_HH_germline %>% filter(locus == "IGH")

# Identify top 15 clones by cell count (restrict to GC B cells, matching earlier logic)
top_clones <- df_heavy %>%
  filter(L1_annotation == "GC_B_cells") %>%
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>%
  slice_head(n = 15) %>%
  pull(clone_subgroup_id_90_similarity)

# For each clone, compute pairwise Levenshtein distances between junction sequences
clone_dist_list <- lapply(top_clones, function(cl) {
  
  # cl <- "500_1"
  
  df_cl <- df_heavy %>%
    filter(clone_subgroup_id_90_similarity == cl) %>%
    arrange(sequence_alignment)
  
  junction_len <- df_cl$junction_length %>% unique()
  
  seqs <- df_cl$junction
  
  # Need at least 2 sequences to compute a distance
  if (length(seqs) < 2) return(NULL)
  
  dist_mat <- stringdistmatrix(seqs, method = "lv")
  dists <- as.numeric(as.dist(dist_mat))
  
  tibble(
    clone_subgroup_id_90_similarity = cl,
    n_cells = length(seqs),
    distance = dists, 
    junction_len = junction_len
  )
  
}) %>% bind_rows() 

# Order clones by number of cells (largest first), keep as factor for plotting
clone_order <- clone_dist_list %>%
  distinct(clone_subgroup_id_90_similarity, n_cells) %>%
  arrange(desc(n_cells)) %>%
  pull(clone_subgroup_id_90_similarity)

clone_dist_list <- clone_dist_list %>%
  mutate(clone_label = glue("{clone_subgroup_id_90_similarity} (n={n_cells}, jl={junction_len})")) %>%
  mutate(clone_label = factor(clone_label, levels = unique(clone_label[match(clone_order, clone_subgroup_id_90_similarity)])))

ggplot(clone_dist_list, aes(x = clone_label, y = distance)) +
  geom_boxplot(outlier.size = 0, fill = "grey90") +
  geom_jitter(width = 0.15, height = 0.3, size = 0.5, alpha = 0.4) +
  # geom_point(size = 0.5, alpha = 0.4) +
  scale_y_continuous(breaks = seq(0, max(clone_dist_list$distance), by = 2)) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = glue("{HH}: pairwise junction distances within top 15 clones"),
    subtitle = "Jitter plot: sequnces only have whole number of mutations",
    x = "Clone",
    y = "Levenshtein distance (N mutations)"
  )

ggsave(glue("{outdir}/00_top15_clones_distance_boxplot.png"), width = 12, height = 6, dpi = 300)

