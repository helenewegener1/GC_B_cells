library(dowser)
library(alakazam)
library(dplyr)
library(ggplot2)
library(shazam)
library(glue)

# Following: https://dowser.readthedocs.io/en/stable/vignettes/Germlines-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

files <- list.files("45_immcantation/out/")
sample_names <- files[1:length(files)-1]

# Load light chain corrected 
clone_10x_list <- lapply(sample_names, function(x){
  
  clone_10x <- read.delim(glue("45_immcantation/out/{x}/{x}_10X_clone-pass.tsv"))
  clone_10x$sample_id <- x %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  clone_10x$cell_id <- paste(clone_10x$cell_id, clone_10x$sample_id, sep = "_")
  clone_10x$sequence_id <- paste(clone_10x$sequence_id, "Heavy", sep = "_")
  
  # # Filter cells based on multiple heavy chains
  # clone_10x <- clone_10x %>%
  #   group_by(cell_id) %>%
  #   arrange(desc(umi_count), .by_group = TRUE) %>%
  #   mutate(
  #     n_heavy = n(),
  #     dominant = case_when(
  #       n_heavy == 1 ~ TRUE,                          # only one heavy chain, keep it
  #       umi_count[1] >= 2 * umi_count[2] ~ row_number() == 1,  # dominant contig has 2x UMIs
  #       TRUE ~ FALSE                                  # ambiguous, drop all contigs for this cell
  #     )
  #   ) %>%
  #   filter(dominant) %>%
  #   select(-n_heavy, -dominant) %>%
  #   ungroup()
  
  return(clone_10x)
  
}) %>% setNames(sample_names)
clone_10x_combined <- bind_rows(clone_10x_list)

# Load light chains 
light_chain_list <- lapply(sample_names, function(x){
  
  # x <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1"
  
  light_chain <- read.delim(glue("45_immcantation/out/{x}/{x}_light_parse-select.tsv"))
  
  light_chain$sample_id <- x %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  light_chain$sample_clean <- x %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC|-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2")
  light_chain$cell_id <- paste(light_chain$cell_id, light_chain$sample_id, sep = "_")
  light_chain$subject_id <- str_split_i(x, "-", 1)
  light_chain$sequence_id <- paste(light_chain$sample_id, light_chain$sequence_id, "Light", sep = "_")
  
  # Filter cells based on multiple light chains
  light_chain <- light_chain %>%
    group_by(cell_id) %>%
    arrange(desc(umi_count), .by_group = TRUE) %>%
    mutate(
      n_light = n(),
      dominant = case_when(
        n_light == 1 ~ TRUE,                          # only one light chain, keep it
        umi_count[1] >= 2 * umi_count[2] ~ row_number() == 1,  # dominant contig has 2x UMIs
        TRUE ~ FALSE                                  # ambiguous, drop all contigs for this cell
      )
    ) %>%
    filter(dominant) %>%
    select(-n_light, -dominant) %>%
    ungroup()
  
  return(light_chain)
  
}) %>% setNames(sample_names)
light_chain_combined <- bind_rows(light_chain_list)


# Combine heavy and light chain in one df
clone_10x_combined$cell_id %>% str_split_i("_", 2) %>% unique()
light_chain_combined$cell_id %>% str_split_i("_", 2) %>% unique()

both_combined_all <- bind_rows(clone_10x_combined, light_chain_combined)

both_combined_all$cell_id %>% str_split_i("_", 2) %>% unique()

both_combined <- list(
  "HH117" = both_combined_all %>% filter(subject_id == "HH117"),
  "HH119" = both_combined_all %>% filter(subject_id == "HH119")
)

# spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")
# 
# HH <- "HH119"
# HH_spec_clones_vj <- spec_clones_vj[[HH]]
# 
# nrow(HH_spec_clones_vj)
# 
# # Check that neccessary columns are present
# HH_spec_clones_vj$sequence_alignment %>% head()
# HH_spec_clones_vj$germline_alignment_d_mask %>% head()

# ------------------------------------------------------------------------------
# resolveLightChains
# ------------------------------------------------------------------------------

# Check that sequence_id and cell_id are unique
both_combined$HH119$sequence_id %>% length()
both_combined$HH119$sequence_id %>% unique() %>% length()
clone_10x_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% length()
clone_10x_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% unique() %>% length()
light_chain_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% length()
light_chain_combined %>% filter(subject_id == "HH119") %>% pull(cell_id) %>% unique() %>% length()
# light_chain_combined %>% filter(subject_id == "HH119") %>% count(cell_id, sort = TRUE)

both_combined$HH117$sequence_id %>% length()
both_combined$HH117$sequence_id %>% unique() %>% length()
clone_10x_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% length()
clone_10x_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% unique() %>% length()
light_chain_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% length()
light_chain_combined %>% filter(subject_id == "HH117") %>% pull(cell_id) %>% unique() %>% length()


# Run resolveLightChains
patients <- names(both_combined)

resolve_LC_list <- lapply(patients, function(HH){
  
  # HH <- "HH117"
  
  resolve_LC_HH <- resolveLightChains(both_combined[[HH]])
  
  table(resolve_LC_HH$celltype_broad, useNA = "always")
  
  # Get meta data from heavy chains
  meta <- resolve_LC_HH %>% 
    filter(!is.na(celltype_broad)) %>% 
    select(cell_id, celltype_broad)
  
  # Clear celltype 
  resolve_LC_HH$celltype_broad <- NULL 
  
  # Add metadata
  resolve_LC_HH <- resolve_LC_HH %>% left_join(meta, by = "cell_id")
  
  table(resolve_LC_HH$celltype_broad, useNA = "always")
  
  return(resolve_LC_HH)
  
}) %>% setNames(patients)

# saveRDS(resolve_LC_list, "45_immcantation/out/rds/resolve_LC_list.rds")
resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list.rds")

# Check output
resolve_LC_list$HH119$clone_id %>% head()
resolve_LC_list$HH119$clone_subgroup %>% head()
resolve_LC_list$HH119$clone_subgroup_id %>% head()

resolve_LC_list$HH117$clone_subgroup %>% table()
resolve_LC_list$HH119$clone_subgroup %>% table()

clone_10x %>% filter(subject_id == "HH117") %>% nrow()
light_chain_combined %>% filter(subject_id == "HH117") %>% nrow()

resolve_LC_list$HH117$cell_id %>% length()
resolve_LC_list$HH117$cell_id %>% unique() %>% length()

resolve_LC_list$HH119$locus %>% table()

# ------------------------------------------------------------------------------
# Define top clones
# ------------------------------------------------------------------------------

patients <- names(resolve_LC_list)

top_GC_clones_vj <- lapply(patients, function(HH) {

  # # find clones that contain at least 1 GC cell
  # GC_clones <- resolve_LC_list[[HH]] %>%
  #   filter(celltype_broad == "GC_B_cells") %>%
  #   pull(clone_id) %>%
  #   unique()
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- resolve_LC_list[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    group_by(clone_id) %>%
    summarise(n_samples = n_distinct(sample_clean_fol)) %>%
    filter(n_samples >= 2) %>%
    pull(clone_id)
  
  # rank those clones by total size (all cell types) and take top 10
  resolve_LC_list[[HH]] %>%
    filter(clone_id %in% GC_clones) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_id)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Look at top clones
# ------------------------------------------------------------------------------

lapply(patients, function(HH){
  
  # HH <- "HH119"
  # HH <- "HH117"
  
  resolve_LC_list[[HH]] %>% 
    filter(clone_id %in% top_GC_clones_vj[[HH]]) %>% 
    mutate(
      v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
    ) %>% 
    mutate(
      j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
    ) %>% 
    count(clone_id, clone_subgroup_id, v_gene, j_gene, sort = TRUE) %>% 
    head(n = 20)
  
}) %>% setNames(patients)

# -------------------
# Look into subclones of top clone
# -------------------

HH <- "HH119"
clone_nr <- 1

resolve_LC_list[[HH]] %>% 
  filter(clone_id == top_GC_clones_vj[[HH]][[clone_nr]]) %>% 
  mutate(
    v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
  ) %>% 
  mutate(
    j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
  ) %>% 
  count(clone_id, clone_subgroup_id, v_gene, j_gene)

# -------------------
# Visualize top clones vj
# -------------------

source("10_broad_annotation/script/color_palette.R")

for (HH in patients){
  
  # HH <- "HH119"
  HH_top_clones <- top_GC_clones_vj[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 2
    clone <- HH_top_clones[clone_nr]
    
    plot_df <- resolve_LC_list[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    j_gene <- plot_df$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    
    plot_df %>%   
      # ggplot(aes(y = sample_clean_fol, fill = celltype_broad)) +
      ggplot(aes(y = sample_clean_fol, fill = clone_subgroup_id)) +
      geom_bar() + 
      # scale_fill_manual(values = celltype_colors) + 
      geom_text(aes(x = Count, label = Count), 
                hjust = -0.2, color = "black",
                stat = "unique") +
      labs(
        title = glue("{HH}: Top {clone_nr} GCB clone spec_vj"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
        y = ""
      ) +
      theme_bw()
    
    ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_resolveLightChains/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5)
    
  }
}

# ------------------------------------------------------------------------------
# Construct clonal germlines
# ------------------------------------------------------------------------------

library(patchwork)

HH <- "HH117"
HH_spec_clones_vj <- resolve_LC_list[[HH]]

# Download (IMGT) germline reference database
# In terminal
# git clone https://github.com/immcantation
# cd immcantation/scripts
# fetch_imgtdb.sh -o germlines # Run script to obtain IMGT gapped sequences

# Read in IMGT-gapped sequences
# references <- readIMGT(dir = "../packages/immcantation/scripts/germlines/human/vdj")

# remove germline alignment columns for this example
db <- select(HH_spec_clones_vj, -"germline_alignment", -"germline_alignment_d_mask")

# Reconstruct germline sequences
HH_spec_clones_vj <- createGermlines(db, references, nproc=1)

# Check germline of first row
HH_spec_clones_vj$germline_alignment_d_mask[1]

# ------------------------------------------------------------------------------
# Build Lineage Trees...
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Format clones
# ------------------------------------------------------------------------------

# Top clone
clone_nr <- 2
clone <- top_GC_clones_vj[[HH]][[clone_nr]]

# Subset data for this example
HH_spec_clones_vj_clone <- HH_spec_clones_vj[HH_spec_clones_vj$clone_id == clone,]

# Meta data
n_cells <- HH_spec_clones_vj_clone %>% pull(cell_id) %>% unique() %>% length()
v_gene <- HH_spec_clones_vj_clone$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
j_gene <- HH_spec_clones_vj_clone$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")

# Add count for identical clones - used for tipsize of tree
HH_spec_clones_vj_clone <- HH_spec_clones_vj_clone %>%
  group_by(clone_id, sequence_alignment) %>%
  mutate(n_identical = n()) %>%
  ungroup()

# Process example data using default settings
clones <- formatClones(
  HH_spec_clones_vj_clone, 
  text_fields = c("c_call", "celltype_broad", "sample_clean_fol"), 
  num_fields=c("n_identical"),
  chain = "HL", 
  light_traits = TRUE
)

print(clones)

# ------------------------------------------------------------------------------
# Build trees
# ------------------------------------------------------------------------------

# Maximum parsimony trees using phangorn.
# clones <- getTrees(clones, nproc=1)

# Maximum likelihood trees using phangorn
clones <- getTrees(clones, build="pml")

# ------------------------------------------------------------------------------
# Plot tree
# ------------------------------------------------------------------------------

# plotTrees(clones)

plotTrees(
  clones, 
  tips="c_call", 
  tipsize="n_identical",
  title = FALSE
)[[1]] + 
  plot_annotation(
    title = glue("{HH}: Clone number {clone_nr} ({clone})"),
    subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
  )

ggsave(glue("45_immcantation/plot/12_dowser_resolve_LC/{HH}_dowser_tree_clone_{clone_nr}_c_call.png"), width = 15, height = 25, dpi = 1000)

plotTrees(
  clones, 
  tips="celltype_broad", 
  tipsize="n_identical",
  title = FALSE
)[[1]] + 
  plot_annotation(
    title = glue("{HH}: Clone number {clone_nr} ({clone})"),
    subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
  )

ggsave(glue("45_immcantation/plot/12_dowser_resolve_LC/{HH}_dowser_tree_clone_{clone_nr}_celltype.png"), width = 15, height = 25, dpi = 1000)

plotTrees(
  clones, 
  tips="sample_clean_fol", 
  tipsize="n_identical",
  title = FALSE
)[[1]] + 
  plot_annotation(
    title = glue("{HH}: Clone number {clone_nr} ({clone})"), 
    subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
  )

ggsave(glue("45_immcantation/plot/11_dowser/{HH}_dowser_tree_clone_{clone_nr}_sample_clean_fol.png"), width = 25, height = 25, dpi = 1000)

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
