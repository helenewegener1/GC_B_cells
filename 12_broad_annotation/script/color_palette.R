# Defining color scheme for each cell type for streamlined plotting

# celltype_colors <- c(
#   # B-cell continuum (Moonrise + Zissou)
#   "Naïve_memory_B_cells"     = "#3B9AB2",  # Moonrise1
#   "GC_like_B_cells"          = "#1F7F5C",  # Moonrise2
#   "GC_B_cells"               = "#78B7C5",  # Zissou1 light blue
# 
#   # Tfh-like (GrandBudapest + Royal)
#   "Tfh_like_cells"           = "#9A6395",  # GrandBudapest1 purple
#   "DCs_MNPs"                 = "#5B1A75",  # Royal1 deep purple
# 
#   # Plasma & PBs — warm (Darjeeling)
#   "PCs_PBs"                  = "#E07B91",  # Darjeeling1 pink/red
# 
#   # Contaminants — soft neutrals (Cavalcanti)
#   "Contamination_ambiguous"  = "#C9A66B",
#   "Contamination_stroma"     = "#A89F91",
#   "Contamination_mast_cells" = "#90745A",
#   "Contamination_MNPs"       = "#736558",
#   "Contamination_γδT_cell"   = "#BFBFBF"   # neutral grey
# )

# All unique celltypes across all samples. 
celltype_colors <- c(
  # B-cell lineage (light/bright teals)
  "Naïve_memory_B_cells"     = "#56B1C8",  # light bright teal
  "GC_like_B_cells"          = "#2FA67F",  # bright turquoise-green
  "GC_B_cells"               = "#9FD3DD",  # pale aqua
  
  # Tfh-like cells (soft lilac)
  "Tfh_like_cells"           = "#C18ACB",  # light lilac
  "DCs_MNPs"                 = "#8A5BAA",  # medium-bright violet
  
  # Plasma / PBs (brighter coral/pink)
  "PCs_PBs"                  = "#F49CA9",  # coral pink
  
  # Contamination categories (pastel earth tones)
  "Contamination_ambiguous"  = "#D9B678",  # light honey
  "Contamination_stroma"     = "#C6BEB3",  # warm beige-grey
  "Contamination_mast_cells" = "#B99674",  # soft tan
  "Contamination_MNPs"       = "#9F8C7A",  # muted mocha
  "Contamination_γδT_cell"   = "#D1D1D1"   # light grey
)
