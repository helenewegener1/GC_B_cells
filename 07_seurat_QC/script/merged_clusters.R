merged_clusters_all <- list(
  
  "v9" = list(
    "HH117-SILP-INF-PC" = list( 
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "3", "4"),
      "1" = c("1"),
      "2" = c("5")
    ),
    
    "HH117-SILP-nonINF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2"),
      "1" = c("1", "3")
    ),
    
    "HH117-SI-MILF-INF-HLADR-AND-CD19" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "3"),
      "1" = c("1"),
      "2" = c("4"),
      "3" = c("5")
    ),
    
    "HH117-SI-MILF-nonINF-HLADR-AND-CD19" = list(
      # new_cluster = c(old clusters)
      "0" = c("0"),
      "1" = c("1", "3", "5", "7"),
      "2" = c("2"),
      "3" = c("4"),
      "6" = c("6"),
      "5" = c("8")
    ),
    
    "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = list(
      # new_cluster = c(old clusters)
      "0" = c("0"),
      "1" = c("1"),
      "2" = c("2", "3"),
      "3" = c("4"),
      "4" = c("5", "6", "7")
    ),
    
    "HH119-COLP-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0"),
      "1" = c("1"),
      "2" = c("2")
    ),
    
    "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "4"),
      "1" = c("1"),
      "2" = c("3"),
      "3" = c("5")
    ),
    
    "HH119-SILP-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2"),
      "1" = c("3"),
      "2" = c("4"),
      "3" = c("5"),
      "4" = c("6"),
      "5" = c("7")
    ),
    
    "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "6"),
      "1" = c("2", "7"),
      "2" = c("3", "4"),
      "3" = c("5")
    ),
    
    "HH119-SI-PP-CD19-Pool1" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "5"),
      "1" = c("2", "3", "4"),
      "2" = c("6"),
      "3" = c("7")
    ),
    
    "HH119-SI-PP-CD19-Pool2" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "4"),
      "1" = c("2", "3"),
      "2" = c("5"),
      "3" = c("6")
    ),
    
    "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "3", "5"),
      "1" = c("1", "6"),
      "2" = c("4"),
      "3" = c("7")
    ),
    
    "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "3"),
      "1" = c("1", "6"),
      "2" = c("4"),
      "3" = c("5")
    ),
    
    # --- HH151 (OCM) samples ---
    # Based on 07_seurat_QC/plot_v9/01_clusters/<sample>/*_clusters.png and the
    # *_broad_{B_cell,T_cell,DC,plasmablast_plasma_cell}.png feature plots:
    # clusters sharing the same broad-marker signature (CD19/CD79A/CD79B/MS4A1/
    # CD40/CD74 for B cells; CD3D/E/G/CD4/TRBC1/CD2/CD7(/CD8) for T cells;
    # JCHAIN/PRDM1/XBP1/MZB1 for plasmablast/plasma cells; HLA-II+LYZ for
    # myeloid/DC) were merged into one cluster. Clusters negative across all
    # four broad panels were kept separate rather than guessed into a category.
    
    "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Blue" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2"),   # B cell (Memory)
      "1" = c("4"),        # T cell
      "2" = c("3"),        # Plasmablast/PC
      "3" = c("1")         # unclassified (negative on all broad panels)
    ),
    
    "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Green" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "3"),   # B cell
      "1" = c("1"),        # T cell
      "2" = c("5"),        # Plasmablast/PC
      "3" = c("2", "4")    # unclassified (negative on all broad panels)
    ),
    
    "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Red" = list(
      # new_cluster = c(old clusters)
      "0" = c("1", "2"),   # B cell
      "1" = c("0"),        # T cell
      "2" = c("3")         # Plasmablast/PC
    ),
    
    "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Yellow" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "3"),   # B cell
      "1" = c("2"),        # T cell
      "2" = c("5"),        # Plasmablast/PC
      "3" = c("1"),        # unclassified, near B cell cluster
      "4" = c("4")         # unclassified, near T cell cluster
    ),
    
    "HH151-SILP-INF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("2"),        # Myeloid/DC (HLA-II very high, LYZ+, CD4+)
      "2" = c("3"),        # B cell (MS4A1+)
      "3" = c("4"),        # T cell
      "4" = c("5")         # unclassified
    ),
    
    "HH151-SILP-nonINF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("2"),        # Myeloid/DC (HLA-II very high, LYZ+)
      "2" = c("3"),        # unclassified CD4+ (CD3-negative, HLA-II-negative)
      "3" = c("4"),        # T cell (CD8+)
      "4" = c("5"),        # B cell
      "5" = c("6")         # unclassified
    )
  ),
  
  "v8" = list(
      
      "HH117-SILP-INF-PC" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "2", "3", "4"),
        "1" = c("1"),
        "2" = c("5")
      ),
      
      "HH117-SILP-nonINF-PC" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "2"),
        "1" = c("1", "3")
      ),
      
      "HH117-SI-MILF-INF-HLADR-AND-CD19" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "2", "4", "6"),
        "1" = c("1"),
        "2" = c("3"),
        "3" = c("5")
      ),
      
      "HH117-SI-MILF-nonINF-HLADR-AND-CD19" = list(
        # new_cluster = c(old clusters)
        "0" = c("0"),
        "1" = c("1", "3", "5", "6"),
        "2" = c("2"),
        "3" = c("4"),
        "4" = c("7"),
        "5" = c("8")
      ),
      
      "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = list(
        # new_cluster = c(old clusters)
        "0" = c("0"),
        "1" = c("1"),
        "2" = c("2", "3"),
        "3" = c("4"),
        "4" = c("5", "6", "7")
      ),
      
      "HH119-COLP-PC" = list(
        # new_cluster = c(old clusters)
        "0" = c("0"),
        "1" = c("1"),
        "2" = c("2")
      ),
      
      "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "2", "4"),
        "1" = c("1"),
        "2" = c("3"),
        "3" = c("5")
      ),
      
      "HH119-SILP-PC" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "1", "2"),
        "1" = c("3"),
        "2" = c("4"),
        "3" = c("5"),
        "4" = c("6"),
        "5" = c("7")
      ),
      
      "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "1", "6"),
        "1" = c("2", "7"),
        "2" = c("3", "5"),
        "3" = c("4")
      ),
      
      "HH119-SI-PP-CD19-Pool1" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "1", "5"),
        "1" = c("2", "3", "4"),
        "2" = c("6"),
        "3" = c("7")
      ),
      
      "HH119-SI-PP-CD19-Pool2" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "1", "3"),
        "1" = c("2", "4"),
        "2" = c("5"),
        "3" = c("6")
      ),
      
      "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "2", "3", "5"),
        "1" = c("1", "6"),
        "2" = c("4"),
        "3" = c("7")
      ),
      
      "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = list(
        # new_cluster = c(old clusters)
        "0" = c("0", "2", "3"),
        "1" = c("1", "6"),
        "2" = c("4"),
        "3" = c("5")
      )
      
    )
  
  
  
)
