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

library(wesanderson)
main_celltypes <- wes_palette("FantasticFox1")
names(main_celltypes) <- c("DCs_MNPs", "Naïve_memory_B_cells", "Tfh_like_cells", "GC_B_cells", "PCs_PBs")

# All unique celltypes across all samples. 
celltype_colors <- c(
  main_celltypes,

  # "Naïve_memory_B_cells"     = "#E58601",  # light bright teal
  # "GC_B_cells"               = "#46ACC8",  # pale aqua
  # 
  # # Tfh-like cells (soft lilac)
  # "Tfh_like_cells"           = "green",  # light lilac
  # "DCs_MNPs"                 = "#DD8D29",  # medium-bright violet
  # 
  # # Plasma / PBs (brighter coral/pink)
  # "PCs_PBs"                  = "#F49CA9",  # coral pink
  
  # Contamination categories (pastel earth tones)
  "Contamination_ambiguous"  = "#D9B678",  # light honey
  "Contamination_stroma"     = "#C6BEB3",  # warm beige-grey
  "Contamination_mast_cells" = "#B99674",  # soft tan
  "Contamination_MNPs"       = "#9F8C7A",  # muted mocha
  "Contamination_γδT_cell"   = "#D1D1D1"   # light grey
)


# L1 annotation 
main_celltypes <- c(wes_palette("FantasticFox1"), "forestgreen")
names(main_celltypes) <- c("Tfh_cells", "Naive_Bcells", "Memory_Bcells", "GC_B_cells", "PCs", "Unconventional_Bcells")

main_celltypes[[1]] <- "hotpink"


# L1_anno$L1_annotation %>% table()

# All unique celltypes across all samples. 
L1_colors <- c(
  
  "Tfh_cells"                  = "#E8608A",
  "Naive_Bcells"               = "#D4C420",
  "Memory_Bcells"              = "#2AAAC8",
  "GC_B_cells"                 = "#E08C20",
  "PCs"                        = "#C42030",
  "Unconventional_Bcells"      = "#8855CC",
  
  # "Tfh_cells"                  = "#F2A0BC",
  # "Naive_Bcells"               = "#E0D875",
  # "Memory_Bcells"              = "#7ACFDF",
  # "GC_B_cells"                 = "#F0B870",
  # "PCs"                        = "#E07880",
  # "Unconventional_Bcells"      = "#B899E8",
    
  "Contamination_mast_cells"  = "#B99674", 
  "Contamination_myeloid_stroma"     = "#D1D1D1"
)


# Name mapping 
patient_names <- list(
  
  "HH117" = "Crohn's Disease", 
  "HH119" = "Colorectal Cancer"
  
)

# Cell type names
cell_type_names <- list(
  
  "Tfh_cells"                  = "Tfh cells",
  "Naive_Bcells"               = "Naive B cells",
  "Memory_Bcells"              = "Memory B cells",
  "GC_B_cells"                 = "GC B cells",
  "PCs"                        = "Plasma cells",
  "Unconventional_Bcells"      = "Unconventional B cells"
  
)
