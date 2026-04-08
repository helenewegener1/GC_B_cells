library(dowser)
library(alakazam)
library(dplyr)
library(tidyverse)
library(shazam)
library(glue)

spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")
resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list.rds")

spec_clones_vj$HH117 %>% nrow()
resolve_LC_list$HH117 %>% filter(locus == "IGH") %>% nrow()
resolve_LC_list$HH117$celltype_broad %>% table(useNA = "always")

spec_clones_vj$HH119 %>% nrow()
resolve_LC_list$HH119 %>% filter(locus == "IGH") %>% nrow()
resolve_LC_list$HH119$celltype_broad %>% table(useNA = "always")

# Check output
resolve_LC_list$HH119$clone_id %>% head()
resolve_LC_list$HH119$clone_subgroup %>% head()
resolve_LC_list$HH119$clone_subgroup_id %>% head()

resolve_LC_list$HH117$clone_subgroup %>% table()
resolve_LC_list$HH119$clone_subgroup %>% table()

resolve_LC_list$HH117$cell_id %>% length()
resolve_LC_list$HH117$cell_id %>% unique() %>% length()

resolve_LC_list$HH119$cell_id %>% length()
resolve_LC_list$HH119$cell_id %>% unique() %>% length()

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
    count(clone_id, clone_subgroup_id, locus, v_gene, j_gene, sort = TRUE) %>%
    head(n = 20)

}) %>% setNames(patients)

# -------------------
# Look into subclones
# -------------------

HH <- "HH117"
clone_nr <- 1
clone <- top_GC_clones_vj[[HH]][[clone_nr]]

# Clones
resolve_LC_list[[HH]] %>%
  filter(clone_id == clone) %>%
  mutate(
    v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
  ) %>%
  mutate(
    j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
  ) %>%
  count(clone_id, clone_subgroup_id, v_gene, j_gene, locus, sort = TRUE) 

# Light chain
# The subgroups are based on v gene, j gene and junction length of light chain. 
# SHM is not considered. 
resolve_LC_list[[HH]] %>%
  filter(clone_id == clone) %>%
  mutate(
    v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
  ) %>%
  mutate(
    j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
  ) %>%
  filter(locus != "IGH") %>% 
  count(v_gene, j_gene, junction_length, sort = TRUE) %>% 
  count(v_gene, j_gene, sort = TRUE)

resolve_LC_list[[HH]] %>% 
  mutate(
    v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
  ) %>%
  mutate(
    j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
  ) %>%
  filter(clone_id == clone, v_gene == "IGKV1-39", j_gene == "IGKJ1") %>% 
  count(clone_id, clone_subgroup_id, v_gene, j_gene, junction_length)

# -------------------
# Visualize top clones vj
# -------------------

for (HH in patients){

  # HH <- "HH117"
  HH_top_clones <- top_GC_clones_vj[[HH]]

  for (clone_nr in 1:length(HH_top_clones)){

    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]

    df <- resolve_LC_list[[HH]] %>%
      filter(clone_id == clone) %>%
      mutate(
        v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
      ) %>%
      mutate(
        j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
      ) %>% 
      mutate(
        sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev(), 
        clone_subgroup_genes = as.character(clone_subgroup) %>% paste(v_gene, j_gene, junction_length, sep = "_")
      )
    
    # Meta data
    n_cells <- df %>% filter(locus == "IGH") %>% nrow()
    v_gene <- df %>% filter(locus == "IGH") %>% pull(v_gene) %>% unique()
    j_gene <- df %>% filter(locus == "IGH") %>% pull(j_gene) %>% unique()
      
    # Plot data
    plot_df <- df %>% 
      filter(locus != "IGH") %>% 
      count(clone_subgroup_genes, celltype_broad, sort = TRUE)

    plot_df %>%
      mutate(clone_subgroup_genes = fct_reorder(
        clone_subgroup_genes,
        as.numeric(str_extract(clone_subgroup_genes, "^\\d+"))
      )) %>% 
      ggplot(aes(y = clone_subgroup_genes, x = celltype_broad, fill = n)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(n > 0, n, "")),
                color = "white", size = 2.5) +
      scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                          limits = c(0, NA)) +
      labs(
        title = glue("{HH}: Top {clone_nr} GCB clone - Light chain subgroups after resolveLightChains"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells heavy chain: {n_cells}\nV gene heavy chain: {v_gene}\n J gene heavy chain: {j_gene}"),
        y = ""
      ) +
      theme_classic() + 
      theme(legend.position = "none")

    ggsave(glue("45_immcantation/plot/GC_clones_spec_vj_resolveLightChains/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 8, height = 8.5)

  }
}

# ------------------------------------------------------------------------------
# Define top subclones for trees
# ------------------------------------------------------------------------------

patients <- names(resolve_LC_list)

top_GC_clones_vj <- lapply(patients, function(HH) {
  
  # # find clones that contain at least 1 GC cell
  # GC_clones <- resolve_LC_list[[HH]] %>%
  #   filter(celltype_broad == "GC_B_cells") %>%
  #   pull(clone_id) %>%
  #   unique()
  
  # HH <- "HH117"
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- resolve_LC_list[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    group_by(clone_subgroup_id) %>%
    summarise(n_samples = n_distinct(sample_clean_fol)) %>%
    filter(n_samples >= 2) %>%
    pull(clone_subgroup_id)
  
  # rank those clones by total size (all cell types) and take top 10
  resolve_LC_list[[HH]] %>%
    filter(clone_subgroup_id %in% GC_clones) %>%
    count(clone_subgroup_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_subgroup_id)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Construct clonal germlines
# ------------------------------------------------------------------------------
source("10_broad_annotation/script/color_palette.R")
library(patchwork)
library(ggtree)

HH <- "HH119"
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
HH_spec_clones_vj <- createGermlines(db, references, nproc=1, clone = "clone_subgroup_id")

# Check germline of first row
HH_spec_clones_vj$germline_alignment_d_mask[1]

# ------------------------------------------------------------------------------
# Build Lineage Trees...
# ------------------------------------------------------------------------------
# https://dowser.readthedocs.io/en/stable/vignettes/Resolve-Light-Chains-Vignette/

# ------------------------------------------------------------------------------
# Format clones
# ------------------------------------------------------------------------------
source("10_broad_annotation/script/color_palette.R")
library(scatterpie)
library(ggtree)

clone_nrs <- 1:5
# clone_nrs <- 2:5

for (clone_nr in clone_nrs){

  # Top clone
  # clone_nr <- 2
  clone <- top_GC_clones_vj[[HH]][[clone_nr]]

  # Subset data for this example
  HH_spec_clones_vj_clone <- HH_spec_clones_vj %>%
    # filter(clone_subgroup_id == paste0(clone, "_1") & celltype_broad == "GC_B_cells")   # _1 = dominant subgroup
    filter(clone_subgroup_id == clone) 

  # Meta data
  n_cells <- HH_spec_clones_vj_clone %>% pull(cell_id) %>% unique() %>% length()
  v_gene <- HH_spec_clones_vj_clone$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
  j_gene <- HH_spec_clones_vj_clone$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")

  # Add count for identical clones - used for tipsize of tree
  HH_spec_clones_vj_clone <- HH_spec_clones_vj_clone %>%
    group_by(clone_id, locus, sequence_alignment) %>%
    mutate(n_identical = n()) %>%
    ungroup()

  # Process example data using default settings
  clones <- formatClones(
    HH_spec_clones_vj_clone,
    clone = "clone_subgroup_id",
    text_fields = c("c_call", "celltype_broad", "sample_clean_fol"),
    num_fields=c("n_identical"),
    chain = "HL",
    light_traits = TRUE
  )

  # print(clones)

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

  if (HH == "HH119" & clone_nr == 1){
    width <- 15
    height <-  25
  } else {
    width <- 10
    height <- 10
  }

  # plotTrees(clones)

  plotTrees(
    clones,
    tips="c_call",
    tipsize="n_identical",
    title = FALSE
  )[[1]] +
    plot_annotation(
      title = glue("{HH}: Clone number {clone_nr} ({clone}_1)"),
      subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
    )

  ggsave(glue("45_immcantation/plot/13_dowser_resolve_LC/{HH}_dowser_tree_clone_{clone_nr}_c_call.png"), width = width, height = height, dpi = 1000)

  plotTrees(
    clones,
    tips="celltype_broad",
    tipsize="n_identical",
    title = FALSE
  )[[1]] +
    plot_annotation(
      title = glue("{HH}: Clone number {clone_nr} ({clone}_1)"),
      subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
    )

  # # extract node data
  # node_data <- p$data %>%
  #   filter(!is.na(celltype_broad)) %>%
  #   select(node, celltype_broad) %>%
  #   separate_rows(celltype_broad, sep = ",") %>%
  #   count(node, celltype_broad) %>%
  #   pivot_wider(names_from = celltype_broad, values_from = n, values_fill = 0)
  #
  # # add pie charts
  # pies <- nodepie(node_data, cols = 2:ncol(node_data),
  #                 color = celltype_colors[c("GC_B_cells", "PCs_PBs")])
  #
  # library(ggimage)
  # p + node_data +
  #   geom_inset(pies, width = 0.1, height = 0.1)

  ggsave(glue("45_immcantation/plot/13_dowser_resolve_LC/{HH}_dowser_tree_clone_{clone_nr}_celltype.png"), width = width, height = height, dpi = 1000)

  plotTrees(
    clones,
    tips="sample_clean_fol",
    tipsize="n_identical",
    title = FALSE
  )[[1]] +
    plot_annotation(
      title = glue("{HH}: Clone number {clone_nr} ({clone}_1)"),
      subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
    )

  ggsave(glue("45_immcantation/plot/13_dowser_resolve_LC/{HH}_dowser_tree_clone_{clone_nr}_sample_clean_fol.png"), width = width, height = height, dpi = 1000)


  # join with your metadata
  tree <- clones$trees[[1]]
  tree_data <- fortify(tree)

  table(tree_data$label[str_detect(tree_data$label, "Heavy")] %>% na.omit() %in% HH_spec_clones_vj_clone$sequence_id[str_detect(HH_spec_clones_vj_clone$sequence_id, "Heavy")])

  tree_data <- tree_data %>%
    left_join(
      HH_spec_clones_vj_clone %>% select(sequence_id, sample_clean_fol,
                         celltype_broad, n_identical, sequence_alignment),
      by = c("label" = "sequence_id")
    )

  # build plot manually
  tree_data %>%
    ggplot(aes(x = x, y = y)) +
    geom_tree() +
    geom_tippoint(aes(color = sample_clean_fol,
                      shape = celltype_broad,
                      size = n_identical)) +
    scale_size_continuous(
      range = c(2, 8), breaks = scales::breaks_width(1)  # only whole number breaks
    ) +
    theme_tree2() +
    labs(
      color = "Sample", shape = "Cell type", size = "N identical",
      title = glue("{HH}: Clone number {clone_nr} ({clone}_1)"),
      subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
    )

  ggsave(glue("45_immcantation/plot/13_dowser_resolve_LC/{HH}_dowser_tree_clone_{clone_nr}_sample_clean_fol_celltype.png"), width = width, height = height, dpi = 1000)

}

# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------
