# Function to clean and update marker names for a single cell type list
update_marker_names <- function(marker_list, seurat_obj) {
  
  # Use map to iterate over each cell type vector in the list
  updated_list <- map(marker_list, function(marker_vector) {
    
    # Use map_chr to iterate over each individual marker in the vector
    map_chr(marker_vector, function(marker) {
      
      # 1. Find gene name (case-insensitive)
      new_name <- grep(marker, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
      
      # 2. Check and return the appropriate name
      if (length(new_name) == 1) {
        return(new_name) # Found exactly one match, use the official name
      } else if (length(new_name) > 1) {
        # Multiple hits found. The original code's logic here was flawed (if (marker %in% marker)),
        # but the intent seems to be: if we can't be sure, keep the original name to avoid errors.
        # We'll print a warning and keep the original name.
        # warning(glue("Multiple hits found for '{marker}'. Keeping original name."))
        return(marker) 
      } else {
        # No matches found. Print an error message and keep the original name.
        # Keeping the original name allows you to manually fix it later.
        print(glue("'{marker}' not found in Seurat object."))
        return(marker) 
      }
    })
  })
  
  # Assign the corrected list back to the original names
  names(updated_list) <- names(marker_list)
  return(updated_list)
}

# Function to generate plots
all_plots <- function(seurat_obj, sample_name, n_dim, extra_title = "", out_dir) {
  
  # N cells 
  n_cells <- ncol(seurat_obj)
 
  ################## Feature plots with continuous features ####################
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", 
                "scDblFinder.score", "S.Score", "G2M.Score", "sce_contamination")
  for (feature in features){
    
    FeaturePlot(seurat_obj, features = feature) + 
      labs(
        title = glue("{extra_title}: {feature}"), 
        subtitle = sample_name,
        caption = glue("N cells: {n_cells}")
      )
    
    ggsave(glue("{out_dir}/{sample_name}_{feature}.pdf"), width = 8, height = 7)
    
  }
  
  ############################################################################## 
  
  ###################### Dimplots with discrete features #######################
  # groups <- c("scDblFinder.class", "Phase")
  groups <- c("Phase")
  
  for (group in groups){
    
    # group = "Phase"
    DimPlot(seurat_obj, group.by = group) + 
      labs(
        title = glue("{extra_title}: {group}"), 
        subtitle = sample_name,
        caption = glue("N cells: {n_cells}")
      )
    
    ggsave(glue("{out_dir}/{sample_name}_{group}.pdf"), width = 8, height = 7)
    
  }
  
  ############################################################################## 
  
  ####################### FeaturePlot with broad_markers ####################### 
  for (markers in names(broad_markers)){
    
    FeaturePlot(seurat_obj, features = broad_markers[[markers]], ncol = 3) + 
      plot_annotation(title = glue("{extra_title}: {markers}"),
                      subtitle = sample_name,
                      caption = glue("N cells: {n_cells}"))
    
    # Adjust heigh of plot to number of markers
    n_markers <- broad_markers[[markers]] %>% length()
    height <- (n_markers/3) * 4

    ggsave(glue("{out_dir}/{sample_name}_broad_{markers}.pdf"), width = 14, height = height)
    
    }
  
  ############################################################################## 
  
  ##################### FeaturePlot with detailed_markers ######################
  for (markers in names(detailed_markers)){
    
    FeaturePlot(seurat_obj, features = detailed_markers[[markers]], ncol = 3) + 
      plot_annotation(title = glue("{extra_title}: {markers}"), 
                      subtitle = sample_name,
                      caption = glue("N cells: {n_cells}"))
    
    # Adjust heigh of plot to number of markers
    n_markers <- detailed_markers[[markers]] %>% length()
    height <- (n_markers/3) * 4
    
    ggsave(glue("{out_dir}/{sample_name}_detailed_{markers}.pdf"), width = 14, height = height)
    
  }
  
  ############################################################################## 
  
  ################################# Clusters ################################# 
  DimPlot(seurat_obj, label = TRUE, group.by = "seurat_clusters") + NoLegend() + 
    labs(
      title = glue("{extra_title}: Seurat clusters"), 
      subtitle = sample_name, 
      caption = glue("N cells: {n_cells}\nN dim: {n_dim}\nResolution: {res}")
    )
  ggsave(glue("{out_dir}/{sample_name}_clusters.pdf"), width = 8, height = 8)
  
}

# Save loadings
save_pc_loading <- function(seurat_obj_list, extra_title){
  
  # Prep to save clustered seurat objects
  pc_loadings <- rep(0, length(seurat_obj_list)) %>% as.list()
  names(pc_loadings) <- sheet_names
  
  for (sample_name in sample_names){
    
    seurat_obj <- seurat_obj_list[[sample_name]]
    
    loadings <- Loadings(seurat_obj, reduction = "pca")
    
    top_genes_per_pc <- loadings %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      pivot_longer(
        cols = starts_with("PC_"),
        names_to = "PC",
        values_to = "loading"
      ) %>%
      mutate(
        PC = factor(PC, levels = paste0("PC_", 1:20))
      ) %>%
      group_by(PC) %>%
      slice_max(abs(loading), n = 20)
    
    # Export to excel file
    sheet_name <- sheet_names[[sample_name]]
    pc_loadings[[sheet_name]] <- top_genes_per_pc
    
    # top_genes_per_pc
    # View(top_genes_per_pc)
    
  }
  
  # Export xlsx file 
  out_file <- glue("09_seurat_QC_clusters/out/PC_loadings_{extra_title}_cells.xlsx")
  
  # Use openxlsx::write.xlsx, which takes the named list and writes
  # each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
  openxlsx::write.xlsx(
    x = pc_loadings,
    file = out_file,
    overwrite = TRUE # Overwrite the file if it already exists
  )
  
}

# Function to do scatter plot of markers for one lineage (e.g. T cells) VS another (e.g. B cells)
# colored by cell cycle phase and doublet classification (singlet/doublet)
# Idea is to detect cells that are marked as doublets but are in fact proliferating cells that we want to save. 
# Also, we are more sure that it is a doublet if it expresses markers from two different lineages 
# (the droplet might contain a T cell and a B cell)
# It is difficult to detect a droplet containing two T cells... 
doublet_dual_lineages <- function(seurat_obj, sample_name, marker_1, marker_2) {
  
  seurat_obj[[]] %>% 
    ggplot(aes(x = !!sym(marker_1), 
               y = !!sym(marker_2), 
               color = scDblFinder.class,
               shape = Phase)) +
    geom_point(alpha = 0.5) + 
    scale_color_manual(values = c("grey", "red")) + 
    theme_bw()
  
  ggsave(glue("09_seurat_QC_clusters/plot/{sample_name}/doublet/{sample_name}_doublet_Phase_{marker_1}_VS_{marker_2}_red.png"))
  
  seurat_obj[[]] %>% 
    mutate(doublet_Phase = glue("{scDblFinder.class}_{Phase}")) %>% 
    ggplot(aes(x = !!sym(marker_1), 
               y = !!sym(marker_2), 
               color = doublet_Phase,
               shape = scDblFinder.class)) +
    geom_point(alpha = 0.5) + 
    scale_color_manual(values = c("grey", "blue", "green", "grey","grey","grey")) +
    theme_bw()
  
  ggsave(glue("09_seurat_QC_clusters/plot/{sample_name}/doublet/{sample_name}_doublet_Phase_{marker_1}_VS_{marker_2}.png"))
  
}

# Function to do violin plots of nFeature_RNA and nCount_RNA
# split on x axis by doublet/singlet class
# color by cell cycle phase
# Idea is to detect cells that are marked as doublets but are infact proliferating cells that we want to save. 
doublet_N_genes <- function(seurat_obj, sample_name, out_dir) {
  
  # Doublet percentage across cell cycle phases. 
  phase_doublet_percentage <- list()
  for (phase in unique(seurat_obj$Phase)){
    
    n_doublets <- sum(seurat_obj$scDblFinder.class == "doublet" & seurat_obj$Phase == phase)
    n_total <- sum(seurat_obj$Phase == phase)
    n_doublet_percentage <- round((n_doublets/n_total) * 100, 1)
    
    phase_doublet_percentage[[phase]] <- n_doublet_percentage
    
  }
  
  text_box <- glue("Doublet percentage across cell cycle phases
                   G1: {phase_doublet_percentage[['G1']]}%
                   GM2: {phase_doublet_percentage[['G2M']]}%
                   S: {phase_doublet_percentage[['S']]}%
                   ")
  
  # nFeature_RNA
  max_nFeature_RNA <- seurat_obj$nFeature_RNA %>% max()
  
  seurat_obj[[]] %>% 
    ggplot(aes(x = scDblFinder.class, 
               y = nFeature_RNA)) + 
    geom_violin() + 
    geom_jitter(aes(color = Phase), width = 0.3, alpha = 0.3) + 
    scale_color_manual(values = c("grey", "red", "blue")) + 
    theme_bw() + 
    annotate("text", x=1.5, y=max_nFeature_RNA-500, label= text_box) 
  
  ggsave(glue("{out_dir}/00_{sample_name}_doublet_VS_phase_nFeature_RNA.png"))
  
  # nCount_RNA
  max_nCount_RNA <- seurat_obj$nCount_RNA %>% max()
  
  seurat_obj[[]] %>%
    ggplot(aes(x = scDblFinder.class,
               y = nCount_RNA)) +
    geom_violin() +
    geom_jitter(aes(color = Phase), width = 0.3, alpha = 0.3) +
    scale_color_manual(values = c("grey", "red", "blue")) +
    theme_bw() + 
    annotate("text", x=1.5, y=max_nCount_RNA-5000, label= text_box) 

  ggsave(glue("{out_dir}/00_{sample_name}_doublet_VS_phase_nCount_RNA.png"))
  
}
