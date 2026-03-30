library(tidyverse) 
library(glue)

# ------------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------------

x <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"
# x <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

light <- read.delim(glue("45_immcantation/out/{x}/{x}_light_parse-select.tsv"))
heavy <- read.delim(glue("45_immcantation/out/{x}/{x}_heavy_germ-pass.tsv"))
vj <- read.delim(glue("45_immcantation/out/{x}/{x}_heavy_germ-pass_clone-pass.tsv"))
clone_10x <- read.delim(glue("45_immcantation/out/{x}/{x}_10X_clone-pass.tsv"))

nrow(light)
nrow(heavy)
nrow(vj)
nrow(clone_10x)

# table(light$cell_id %in% heavy$cell_id)
# 
# table(light$cell_id %in% vj$cell_id)
# table(heavy$cell_id %in% vj$cell_id)
# 
# table(light$cell_id %in% clone_10x$cell_id)
# table(heavy$cell_id %in% clone_10x$cell_id)

table(vj$cell_id %in% clone_10x$cell_id)

dowser::resolveLightChains

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# clone_10x <- readRDS("45_immcantation/out/rds/clone_10x.rds")
# patients <- names(clone_10x)

# Load light chain corrected 
files <- list.files("45_immcantation/out/")
sample_names <- files[1:length(files)-1]

clone_10x_list <- lapply(sample_names, function(x){
  
  clone_10x <- read.delim(glue("45_immcantation/out/{x}/{x}_10X_clone-pass.tsv"))
  
  return(clone_10x)
  
}) %>% setNames(sample_names)

clone_10x_combined <- bind_rows(clone_10x_list)

# ------------------------------------------------------------------------------
# Alter format
# ------------------------------------------------------------------------------

clone_10x <- list(
  "HH117" = clone_10x_combined %>% filter(subject_id == "HH117"),
  "HH119" = clone_10x_combined %>% filter(subject_id == "HH119")
)

# ------------------------------------------------------------------------------
# Inspect HH119 
# ------------------------------------------------------------------------------

HH <- "HH119"
clone_10x[[HH]] %>% 
  count(clone_id, sort = TRUE) %>% head()
# count(clone_id, v_call, j_call, sort = TRUE)

# Plot a histogram of inter and intra clonal distances
# plot(clone_10x[[HH]], binwidth=0.02)

top_clone <- clone_10x[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

clone_10x[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE) %>% 
  head(n = 10 )

# ------------------------------------------------------------------------------
# Inspect HH117
# ------------------------------------------------------------------------------

# HH117 has larger clones than with hierical clustering.

HH <- "HH117"
clone_10x[[HH]] %>% 
  count(clone_id, sort = TRUE) %>% 
# count(clone_id, v_call, j_call, sort = TRUE)
  head()

# Plot a histogram of inter and intra clonal distances
# plot(clone_10x[[HH]], binwidth=0.02)

top_clone <- clone_10x[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

clone_10x[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
clone_10x[[HH]] %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# ------------------------------------------------------------------------------
# Define top clones
# ------------------------------------------------------------------------------

top_GC_clones_vj <- lapply(patients, function(HH) {
  
  # # find clones that contain at least 1 GC cell
  # GC_clones <- clone_10x[[HH]] %>%
  #   filter(celltype_broad == "GC_B_cells") %>%
  #   pull(clone_id) %>%
  #   unique()
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- clone_10x[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    group_by(clone_id) %>%
    summarise(n_samples = n_distinct(sample_clean_fol)) %>%
    filter(n_samples >= 2) %>%
    pull(clone_id)
  
  # rank those clones by total size (all cell types) and take top 10
  clone_10x[[HH]] %>%
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
  
  clone_10x[[HH]] %>% 
    filter(clone_id %in% top_GC_clones_vj[[HH]]) %>% 
    mutate(
      v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
    ) %>% 
    mutate(
      j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
    ) %>% 
    count(clone_id, v_gene, j_gene, sort = TRUE) 
  
}) %>% setNames(patients)

# -------------------
# Visualize top clones vj
# -------------------

source("10_broad_annotation/script/color_palette.R")

for (HH in patients){
  
  # HH <- "HH117"
  HH_top_clones <- top_GC_clones_vj[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    plot_df <- clone_10x[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    j_gene <- plot_df$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    
    plot_df %>%   
      ggplot(aes(y = sample_clean_fol, fill = celltype_broad)) +
      geom_bar() + 
      scale_fill_manual(values = celltype_colors) + 
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
    
    ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_ligth_chain_corrected/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5)
    
  }
}

