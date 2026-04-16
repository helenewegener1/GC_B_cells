library(dowser)
library(alakazam)
library(dplyr)
library(tidyverse)
library(shazam)
library(glue)

spec_clones_vj <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")
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

# Cell types across samples

table(spec_clones_vj$HH117$celltype_broad, spec_clones_vj$HH117$sample_clean)
table(spec_clones_vj$HH119$celltype_broad, spec_clones_vj$HH119$sample_clean)

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
      filter(str_starts(clone_subgroup_id, paste0(clone, "_"))) %>%
      # filter(clone_id == clone) %>%
      mutate(
        sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev(), 
        clone_subgroup_genes = as.character(clone_subgroup) %>% paste(v_call, j_call, junction_length, sep = "_")
      )
    
    # Meta data
    n_cells <- df %>% filter(locus == "IGH") %>% nrow()
    v_gene <- df %>% filter(locus == "IGH") %>% pull(v_call_majority) %>% unique()
    j_gene <- df %>% filter(locus == "IGH") %>% pull(j_call_majority) %>% unique()
      
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

top_GC_subclones_vj <- lapply(patients, function(HH) {
  
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

# Download (IMGT) germline reference database
# In terminal
# git clone https://github.com/immcantation
# cd immcantation/scripts
# fetch_imgtdb.sh -o germlines # Run script to obtain IMGT gapped sequences

# # Read in IMGT-gapped sequences
# references <- readIMGT(dir = "../packages/immcantation/scripts/germlines/human/vdj")
# 
# resolve_LC_list_germlined <- list()
# 
# for (HH in patients){
#   
#   HH_spec_clones_vj <- resolve_LC_list[[HH]]
#   
#   HH_spec_clones_vj <- HH_spec_clones_vj %>% 
#     mutate(
#       v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
#     ) %>%
#     mutate(
#       j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
#     )
#   
#   # remove germline alignment columns for this example
#   db <- select(HH_spec_clones_vj, -"germline_alignment", -"germline_alignment_d_mask")
#   
#   # Reconstruct germline sequences
#   HH_spec_clones_vj <- createGermlines(db, references, nproc=1, clone = "clone_subgroup_id")
#   
#   # Check germline of first row
#   HH_spec_clones_vj$germline_alignment_d_mask[1]
#   
#   resolve_LC_list_germlined[[HH]] <- HH_spec_clones_vj
#   
# }
# 
# saveRDS(resolve_LC_list_germlined, "45_immcantation/out/rds/resolve_LC_list_germlined.rds")
resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")

# ------------------------------------------------------------------------------
# Build Lineage Trees per subclone
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

versions <- c("all cells", "only GCs")

for (HH in patients){
  
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Top clone
    # clone_nr <- 2
    clone <- top_GC_subclones_vj[[HH]][[clone_nr]]
    
    # Subset data for this example
    for (version in versions){
      
      # version <- "all cells"
      # version <- "only GCs"
      
      if (version == "all cells"){
        
        HH_spec_clones_vj_clone <- HH_spec_clones_vj %>%
          filter(clone_subgroup_id == clone)
        
        outdir <- "45_immcantation/plot/13_dowser_resolve_LC_all_cells"
        dir.create(outdir, recursive = TRUE)
        
      } else if (version == "only GCs"){
        
        HH_spec_clones_vj_clone <- HH_spec_clones_vj %>%
          filter(clone_subgroup_id == clone & celltype_broad == "GC_B_cells") 
        
        outdir <- "45_immcantation/plot/13_dowser_resolve_LC_only_GCs"
        dir.create(outdir, recursive = TRUE)
        
      }
      
      # Meta data
      n_cells <- HH_spec_clones_vj_clone %>% pull(cell_id) %>% unique() %>% length()
      v_gene <- HH_spec_clones_vj_clone %>% filter(locus == "IGH") %>% pull(v_gene) %>% unique()
      j_gene <- HH_spec_clones_vj_clone %>% filter(locus == "IGH") %>% pull(j_gene) %>% unique()
      
      # Add count for identical clones - used for tipsize of tree
      HH_spec_clones_vj_clone <- HH_spec_clones_vj_clone %>%
        group_by(clone_subgroup_id, locus, sequence_alignment) %>%
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
        width <- 20
        height <- 30
      } else {
        width <- 10
        height <- 15
      }
      
      # plotTrees(clones)
      
      plotTrees(
        clones,
        tips="c_call",
        tipsize="n_identical",
        title = FALSE
      )[[1]] +
        plot_annotation(
          title = glue("{HH} {version}: Clone number {clone_nr} ({clone})"),
          subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
        )
      
      # TRY TO ADD THIS FOR LARGER COLORS DOTS IN THE LEGEND
      # guides(color = guide_legend(override.aes = list(size = 4))) +
      
      ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_c_call.png"), width = width, height = height, dpi = 1000)
      
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
      
      ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_celltype.png"), width = width, height = height, dpi = 1000)
      
      if (!(HH == "HH119" & clone_nr == 1)){
        
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
        
        ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_sample_clean_fol.png"), width = width, height = height, dpi = 1000)
        
      }
      
      if (version == "all cells" & !(HH == "HH119" & clone_nr == 1)){
        # join with your metadata
        tree <- clones$trees[[1]]
        tree_data <- fortify(tree)
        
        # table(tree_data$label[str_detect(tree_data$label, "Heavy")] %>% na.omit() %in% HH_spec_clones_vj_clone$sequence_id[str_detect(HH_spec_clones_vj_clone$sequence_id, "Heavy")])
        
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
            title = glue("{HH}: Clone number {clone_nr} ({clone})"),
            subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
          )
        
        ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_sample_clean_fol_celltype.png"), width = width, height = height, dpi = 1000)
        
      }
      
    }
    
  }
  
}


# ------------------------------------------------------------------------------
# Clones per follicle
# ------------------------------------------------------------------------------

clone_nrs <- 1:5
# clone_nrs <- 2:5

# Get non-follcile samples
sample_names <- HH_spec_clones_vj$sample_clean_fol %>% unique()
non_fol_samples <- sample_names[!str_detect(sample_names, "Fol")]

# Get follicle sample names 
follicle_sample_names <- sample_names[str_detect(sample_names, "Fol")]

for (HH in patients){
  
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # clone_nr <- 1
    
    for (sample_name in follicle_sample_names){
      
      # Define sample and clone_nr to look at
      # sample_name <- "HH117-SI-PP-nonINF_Fol-2"
      # sample_name <- "HH119-SI-PP_Fol-22"
      
      outdir <- glue("45_immcantation/plot/13_dowser_resolve_LC_per_follicle/{sample_name}")
      dir.create(outdir, recursive = TRUE)
      
      # Get top "clone_nr" clone from given sample
      clone <- HH_spec_clones_vj %>%
        filter(sample_clean_fol == sample_name) %>% 
        count(clone_subgroup_id, sort = TRUE) %>% 
        slice(clone_nr) %>% 
        pull(clone_subgroup_id)
      
      # Inspect clone
      # HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(v_gene, j_gene, sort = TRUE)
      
      # Get clone from specic follicle and the non-fol samples
      HH_spec_clones_vj_clone <- HH_spec_clones_vj %>% 
        filter(clone_subgroup_id == clone & sample_clean_fol %in% c(sample_name, non_fol_samples))
      
      # Meta data
      n_cells <- HH_spec_clones_vj_clone %>% pull(cell_id) %>% unique() %>% length()
      # sub_clones <- HH_spec_clones_vj_clone %>% count(clone_subgroup_id, sort = TRUE)
      v_gene <- HH_spec_clones_vj_clone %>% filter(locus == "IGH") %>% pull(v_gene) %>% unique()
      j_gene <- HH_spec_clones_vj_clone %>% filter(locus == "IGH") %>% pull(j_gene) %>% unique()
      
      # Add count for identical clones - used for tipsize of tree
      HH_spec_clones_vj_clone <- HH_spec_clones_vj_clone %>%
        group_by(clone_subgroup_id, locus, sequence_alignment) %>%
        # group_by(clone_id, locus, sequence_alignment) %>%
        mutate(n_identical = n()) %>%
        ungroup()
      
      if (nrow(HH_spec_clones_vj_clone) < 3){
        next
      }
      
      # Process example data using default settings
      clones <- formatClones(
        HH_spec_clones_vj_clone,
        clone = "clone_subgroup_id",
        # clone = "clone_id",
        text_fields = c("c_call", "celltype_broad", "sample_clean_fol"),
        num_fields=c("n_identical"),
        chain = "HL",
        light_traits = TRUE
      )
      
      if (nrow(clones) == 0){
        next
      }
      
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
        width <- 20
        height <- 30
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
          title = glue("{HH} {sample_name}: Clone number {clone_nr} ({clone})"),
          subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
        )
      
      ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_c_call.png"), width = width, height = height, dpi = 1000)
      
      # plotTrees(
      #   clones,
      #   tips="celltype_broad",
      #   tipsize="n_identical",
      #   title = FALSE
      # )[[1]] +
      #   plot_annotation(
      #     title = glue("{HH} {sample_name}: Clone number {clone_nr} ({clone}_1)"),
      #     subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
      #   )
      # 
      # ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_celltype.png"), width = width, height = height, dpi = 1000)
      # 
      # plotTrees(
      #   clones,
      #   tips="sample_clean_fol",
      #   tipsize="n_identical",
      #   title = FALSE
      # )[[1]] +
      #   plot_annotation(
      #     title = glue("{HH} {sample_name}: Clone number {clone_nr} ({clone})"),
      #     subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
      #   )
      # 
      # ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_sample_clean_fol.png"), width = width, height = height, dpi = 1000)
      # 
      # join with your metadata
      tree <- clones$trees[[1]]
      tree_data <- fortify(tree)
      
      # table(tree_data$label[str_detect(tree_data$label, "Heavy")] %>% na.omit() %in% HH_spec_clones_vj_clone$sequence_id[str_detect(HH_spec_clones_vj_clone$sequence_id, "Heavy")])
      
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
        # scale_size_continuous(
        #   range = c(2, 8), breaks = scales::breaks_width(1)  # only whole number breaks
        # ) +
        theme_tree2() +
        labs(
          color = "Sample", shape = "Cell type", size = "N identical",
          title = glue("{HH} {sample_name}: Clone number {clone_nr} ({clone})"),
          subtitle = glue("N cells: {n_cells}, V gene: {v_gene}, J gene: {j_gene}")
        ) 
      
      ggsave(glue("{outdir}/{HH}_dowser_tree_clone_{clone_nr}_sample_clean_fol_celltype.png"), width = width, height = height, dpi = 1000)
      
    }
    
  }
  
}

# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Compare methods (spectralClones on ligth chain VS resolveLigthChains)
# ------------------------------------------------------------------------------

heavy_clones <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds") # Before subclustering
resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list.rds")
clones_NA_resolved <- readRDS("45_immcantation/out/rds/07_clones_combined_NA_resolved.rds")

HH <- "HH117"

# N cells
heavy_clones[[HH]] %>% nrow()
resolve_LC_list[[HH]] %>% filter(locus == "IGH") %>% nrow() 
clones_NA_resolved[[HH]] %>% filter(locus == "IGH") %>% nrow()

top_clone <- clones_NA_resolved[[HH]] %>% filter(locus == "IGH") %>% count(clone_id, sort = TRUE) %>% head(n = 1) %>% pull(clone_id)
# resolve_LC_list[[HH]] %>% filter(locus == "IGH") %>% count(clone_id, sort = TRUE) %>% head(n = 1) %>% pull(clone_id)

# N cells in top clone
heavy_clones[[HH]] %>% filter(clone_id == top_clone) %>% nrow()
resolve_LC_list[[HH]]  %>% filter(locus == "IGH" & clone_id == top_clone) %>% nrow()
clones_NA_resolved[[HH]] %>% filter(locus == "IGH" & clone_id == top_clone) %>% nrow()

# N subclones in top clone
resolve_LC_list[[HH]] %>% filter(locus != "IGH", str_starts(clone_subgroup_id, paste0(top_clone, "_"))) %>% select(clone_subgroup_id) %>% unique() %>% nrow()
clones_NA_resolved[[HH]] %>% filter(locus != "IGH", str_starts(clone_id_combine, paste0(top_clone, "_"))) %>% select(clone_id_combine) %>% unique() %>% nrow()

# N cells in subclones in top clone
resolve_LC_list[[HH]] %>% filter(locus != "IGH", str_starts(clone_subgroup_id, paste0(top_clone, "_"))) %>% count(clone_subgroup_id, sort = TRUE)
clones_NA_resolved[[HH]] %>% filter(locus != "IGH", str_starts(clone_id_combine, paste0(top_clone, "_"))) %>% count(clone_id_combine, sort = TRUE)

# Which cells end up in which subclones?
df_LC <- resolve_LC_list[[HH]] %>% filter(locus != "IGH", str_starts(clone_subgroup_id, paste0(top_clone, "_"))) %>% select(clone_subgroup_id, v_call, j_call, junction_length)
df_sp <- clones_NA_resolved[[HH]] %>% filter(locus != "IGH", str_starts(clone_id_combine, paste0(top_clone, "_"))) %>% select(clone_id_combine, v_call, j_call, junction_length)

LC <- df_LC %>% mutate(LC_clones = as.character(clone_subgroup_id) %>% paste(v_call, j_call, junction_length, sep = "_")) %>% pull(LC_clones)
sp <- df_sp %>% mutate(sp_clones = as.character(clone_id_combine) %>% paste(v_call, j_call, junction_length, sep = "_")) %>% pull(sp_clones)

# Build co-occurrence table
cooccurrence <- table(LC = LC, sp = sp) %>% as.data.frame()
cooccurrence <- cooccurrence %>% arrange(desc(Freq)) 

LC_levels <- cooccurrence$LC %>% unique()
sp_levels <- cooccurrence$sp %>% unique()

cooccurrence %>% mutate(
  LC = factor(LC, levels = LC_levels),
  sp = factor(sp, levels = sp_levels)
) %>% 
  ggplot(aes(x = LC, y = sp, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white", size = 3) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_gradient(low = "white", high = "#0d2a4e",
                      limits = c(0, NA)) +
  labs(
    x = "LC (resolveLightChains)",
    y = "sp (spectralClones)",
    fill = "N cells",
    title = glue("{HH} clone {top_clone}")
  )

ggsave(glue("45_immcantation/plot/13_compare_methods/{HH}_clone{top_clone}.png"), dpi = 1000, height = 25, width = 25)


# Investigate top subclones
LC_top_clone <- LC %>% table() %>% which.max() %>% names()
sp_top_clone <- sp %>% table() %>% which.max() %>% names()

resolve_LC_list[[HH]] %>% filter(locus != "IGH", clone_subgroup_id == LC_top_clone) %>% count(v_call, j_call)


