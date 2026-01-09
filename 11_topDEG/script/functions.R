top_DEGs_to_excel <- function(seurat_obj, sample_name, n_DEGs = 100) {
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) + NoLegend()
  
  cluster.name <- "seurat_clusters"
  # cluster.name <- "merged_clusters"
  
  # Find all markers 
  all_markers <- FindAllMarkers(seurat_obj, 
                                group.by = cluster.name, 
                                only.pos = TRUE)
  
  # Process markers
  top_markers_list <- all_markers %>%
    filter(p_val_adj < 0.05) %>% # Filter for significant markers
    group_by(cluster) %>%
    # Sort by adjusted p-value (most significant) and then avg_log2FC (highest expression)
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = n_DEGs) %>%
    ungroup() %>%
    # The split() function creates the named list needed for openxlsx
    split(., .$cluster)
  
  top_markers_list_print <- all_markers %>%
    filter(p_val_adj < 0.05) %>% # Filter for significant markers
    group_by(cluster) %>%
    # Sort by adjusted p-value (most significant) and then avg_log2FC (highest expression)
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup() %>%
    # The split() function creates the named list needed for openxlsx
    split(., .$cluster)
  
  print(top_markers_list_print)
  
  # Export xlsx file 
  out_file <- glue("11_topDEG/out/{sample_name}_Top_{n_DEGs}_DEGs.xlsx")
  
  # Use openxlsx::write.xlsx, which takes the named list and writes
  # each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
  openxlsx::write.xlsx(
    x = top_markers_list,
    file = out_file,
    overwrite = TRUE # Overwrite the file if it already exists
  )
  
}