# Resolution 0.1
merged_clusters_all_res01 <- list(
  
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
      "0" = c("0", "1", "2", "4"),
      "1" = c("3"),
      "2" = c("5")
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
    ),
    
    "HH153-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB-Pool1" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "3", "4"),   # B cell
      "1" = c("2"),                  # T cell
      "2" = c("5"),                  # Plasmablast/PC
      "3" = c("6")                   # Myeloid/DC (HLA-II+, LYZ+)
    ),
    
    "HH153-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB-Pool2" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2", "5"),   # B cell
      "1" = c("4"),                  # T cell
      "2" = c("6"),                  # Plasmablast/PC
      "3" = c("3")                   # unclassified (negative on all broad panels)
    ),
    
    "HH153-SILP-INF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("2"),        # Myeloid/DC (HLA-II very high)
      "2" = c("3"),        # T cell
      "3" = c("4")         # unclassified
    ),
    
    "HH153-SILP-nonINF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("3", "7"),   # Myeloid/DC (HLA-II+, LYZ+, ITGAX+)
      "2" = c("4"),        # T cell
      "3" = c("6"),        # B cell (MS4A1+)
      "4" = c("2"),        # unclassified
      "5" = c("5")         # unclassified
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

# Resolution 0.3
merged_clusters_all_res03 <- list(
  
  "v9" = list(
    "HH117-SILP-INF-PC" = list( 
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "3", "4", "5"), # PCs
      "1" = c("1"), # PCs
      "2" = c("6")  # unclassified
    ),
    
    "HH117-SILP-nonINF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "3", "5"), # PCs
      "1" = c("2", "4", "6") # PCs
    ),
    
    "HH117-SI-MILF-INF-HLADR-AND-CD19" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2", "5", "8", "9"), # PCs
      "1" = c("3", "6", "10"), # DCs
      "2" = c("4"), # B mem + Naive
      "3" = c("7") # unclassified 
    ),
    
    "HH117-SI-MILF-nonINF-HLADR-AND-CD19" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1"), # B mem + Naive
      "1" = c("2", "4", "6", "7", "10"), # DCs
      "2" = c("3"), # PCs
      "3" = c("5"), # T cells
      "6" = c("8"), # GCB
      "5" = c("9") # unclassified
    ),
    
    "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "5", "6"), # B mem + Naive
      "1" = c("1", "3", "8"), # GCB
      "2" = c("2", "9"), # T cell
      "3" = c("4", "10"), # PCs 
      "4" = c("7", "11", "12", "13") # DCs cont.
    ),
    
    "HH119-COLP-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2", "3", "5"), # PCs
      "1" = c("4"), # Likely PCs
      "2" = c("6") # unclassified
    ),
    
    "HH119-CO-SMILF-CD19-AND-GC-AND-PB-AND-TFH" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1","2", "3", "4", "5","6", "8"), # B mem + Naive
      "1" = c("7"), # T cell 
      "2" = c("9") # PCs 
    ),
    
    "HH119-SILP-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2", "3", "4"), # PCs
      "1" = c("7", "8"), # B cell cont. 
      "2" = c("5"), # unclassified
      "3" = c("6"), # DC cont.
      "4" = c("9"), # unclassified
      "5" = c("10") # T cell cont.
    ),
    
    "HH119-SI-MILF-CD19-AND-GC-AND-PB-AND-TFH" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "6"), # B mem + Naive
      "1" = c("2", "8"), # T cell
      "2" = c("3", "5", "7"), # GCB
      "3" = c("4", "9") # PCs
    ),
    
    "HH119-SI-PP-CD19-Pool1" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "5", "9"), # B mem + Naive 
      "1" = c("2", "3", "4", "7"), # GCB
      "2" = c("6"), # PCs
      "3" = c("8") # T cell cont.
    ),
    
    "HH119-SI-PP-CD19-Pool2" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "4", "8"), # B mem + Naive 
      "1" = c("2", "3", "5", "6"), # GCB
      "2" = c("7"), # PCs
      "3" = c("9") # T cell cont.
    ),
    
    "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "3", "5", "7"), # GCB
      "1" = c("1", "6"), # T cell
      "2" = c("4", "9"), # PCs
      "3" = c("8 ") # B cell cont.
    ),
    
    "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "3", "4", "7"), # GCB
      "1" = c("2", "8", "9"), # T cell
      "2" = c("5"), # B cell cont.
      "3" = c("6") # PCs
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
      "0" = c("0", "1", "2", "4"),   # Mem B cell + Naive 
      "1" = c("3"),                  # Plasmablast/PC
      "2" = c("5"),                  # Likely GCB but very few cells
      "3" = c("6")                   # T cell
    ),
    
    "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Green" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2"),             # T cell
      "1" = c("1", "5", "7", "8"),   # Mem B cell + Naive 
      "2" = c("3", "4"),             # GCB
      "3" = c("6")                   # Plasmablast/PC
    ),
    
    "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Red" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2"),   # Mem B cell + Naive 
      "1" = c("3"),             # GCB? Uncertain
      "2" = c("4")              # T cell
    ),
    
    "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB_Yellow" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2", "7"),   # Mem B cell + Naive 
      "1" = c("3", "5"),             # T cell
      "2" = c("4"),                  # GCB
      "3" = c("6")                   # Plasmablast/PC
    ),
    
    "HH151-SILP-INF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2", "3", "6", "7"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("4"),                            # Myeloid/DC (HLA-II very high, LYZ+, CD4+) weirdly abundant?
      "2" = c("5"),                            # unclassified
      "3" = c("8")                             # unclassified
    ),
    
    "HH151-SILP-nonINF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2", "3", "6", "7"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("4"),                            # Myeloid/DC (HLA-II very high, LYZ+) weirdly abundant?
      "2" = c("5"),                            # Maybe PCs/ B cells 
      "3" = c("8")                             # T cell (CD8+) cont.
    ),
    
    "HH153-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB-Pool1" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "3", "6", "7"),   # GCB 
      "1" = c("2", "8"),                  # Mem B cell + Naive 
      "2" = c("4", "5"),                  # T cell
      "3" = c("9"),                       # Plasmablast/PC
      "4" = c("10")                       # Myeloid/DC (HLA-II very high, LYZ+) cont.
    ),
    
    "HH153-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB-Pool2" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "2", "3", "5", "6"),   # GCB
      "1" = c("1", "8"),                  # Mem B cell + Naive 
      "2" = c("4", "7"),                  # T cell
      "3" = c("9"),                       # Myeloid/DC (HLA-II very high, LYZ+) cont.
      "4" = c("10")                       # Plasmablast/PC
    ),
    
    "HH153-SILP-INF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "2"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("3"),             # Maybe PCs/ B cells or Myeloid?
      "2" = c("4"),             # unclassified
      "3" = c("5", "6")         # Myeloid/DC (HLA-II very high, LYZ+) cont.
    ),
    
    "HH153-SILP-nonINF-PC" = list(
      # new_cluster = c(old clusters)
      "0" = c("0", "1", "6", "10"),   # Plasmablast/PC (dominant, PC-sorted tissue)
      "1" = c("2", "3"),              # Maybe PCs/B cells
      "2" = c("4"),                   # T cell cont.
      "3" = c("5", "11", "12"),       # unclassified
      "4" = c("7", "8"),              # Myeloid/DC (HLA-II very high, LYZ+) cont.
      "5" = c("9")                    # unclassified
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

