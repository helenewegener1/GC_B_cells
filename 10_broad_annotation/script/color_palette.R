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

isotype_colors_custom <- c(
  "IGHA1" = "#FF7F00",
  "IGHA2" = "#E31A1C",
  "IGHM" = "#1F78B4",
  "IGHD" = "#00E5FF",
  "IGHG1" = "#3F007D",
  "IGHG2" = "#54278F",
  "IGHG3" = "#756BB1",
  "IGHG4" = "#9E9AC8"
)

isotype_grouped_colors_custom <- c(
  "IGHA1" = "#FF7F00",
  "IGHA2" = "#E31A1C",
  "IGHM/D" = "#1F78B4",
  "IGHG1" = "#3F007D",
  "IGHG2" = "#54278F",
  "IGHG3" = "#756BB1",
  "IGHG4" = "#9E9AC8"
)


# Sample colors 
sample_clean_plot_colors <- c(
    # HH117 samples
    "HH117-SI-MILF-INF"          = "#C42030",
    "HH117-SILP-INF"             = "#E08C20",
    "HH117-SI-MILF-nonINF"       = "#2AAAC8",
    "HH117-SILP-nonINF"          = "#1A6090",
    
    "HH117-SI-PP-nonINF_Fol-1"   = "#8855CC",
    "HH117-SI-PP-nonINF_Fol-2"   = "#A066DD",
    "HH117-SI-PP-nonINF_Fol-3"   = "#55CC77",
    "HH117-SI-PP-nonINF_Fol-4"   = "#44BB66",
    "HH117-SI-PP-nonINF_Fol-5"   = "#33AA55",
    "HH117-SI-PP-nonINF_Fol-6"   = "#D4C420",
    "HH117-SI-PP-nonINF_Fol-7"   = "#C4B410",
    "HH117-SI-PP-nonINF_Fol-8"   = "#E8608A",
    "HH117-SI-PP-nonINF_Fol-9"   = "#D85070",
    "HH117-SI-PP-nonINF_Fol-10"  = "#EF19EC",
    "HH117-SI-PP-nonINF_Fol-11"  = "#CC10CC",
    "HH117-SI-PP-nonINF_Fol-12"  = "#FF8C00",
    "HH117-SI-PP-nonINF_Fol-13"  = "#20B2AA",
    "HH117-SI-PP-nonINF_Fol-14"  = "#10A090",
    "HH117-SI-PP-nonINF_Fol-15"  = "#6A5ACD",
    "HH117-SI-PP-nonINF_Fol-16"  = "#7B68EE",
    "HH117-SI-PP-nonINF_Fol-17"  = "#228B22",
    "HH117-SI-PP-nonINF_Fol-18"  = "#32CD32",
    
    # HH119 samples
    "HH119-CO-SMILF"      = "#C42030",
    "HH119-COLP"          = "#E08C20",
    "HH119-SILP"          = "#1A6090",
    "HH119-SI-MILF"       = "#2AAAC8",
    "HH119-SI-PP_Fol-1"   = "#8855CC",
    "HH119-SI-PP_Fol-2"   = "#A066DD",
    "HH119-SI-PP_Fol-3"   = "#55CC77",
    "HH119-SI-PP_Fol-4"   = "#44BB66",
    "HH119-SI-PP_Fol-5"   = "#33AA55",
    "HH119-SI-PP_Fol-6"   = "#D4C420",
    "HH119-SI-PP_Fol-7"   = "#C4B410",
    "HH119-SI-PP_Fol-8"   = "#E8608A",
    "HH119-SI-PP_Fol-9"   = "#D85070",
    "HH119-SI-PP_Fol-10"  = "#EF19EC",
    "HH119-SI-PP_Fol-11"  = "#CC10CC",
    "HH119-SI-PP_Fol-12"  = "#FF8C00",
    "HH119-SI-PP_Fol-13"  = "#20B2AA",
    "HH119-SI-PP_Fol-14"  = "#10A090",
    "HH119-SI-PP_Fol-15"  = "#6A5ACD",
    "HH119-SI-PP_Fol-16"  = "#7B68EE",
    "HH119-SI-PP_Fol-17"  = "#228B22",
    "HH119-SI-PP_Fol-18"  = "#32CD32",
    "HH119-SI-PP_Fol-19"  = "#C42030",
    "HH119-SI-PP_Fol-20"  = "#E08C20",
    "HH119-SI-PP_Fol-21"  = "#2AAAC8",
    "HH119-SI-PP_Fol-22"  = "#1A6090",
    "HH119-SI-PP_Fol-23"  = "#8855CC",
    "HH119-SI-PP_Fol-24"  = "#55CC77",
    "HH119-SI-PP_Fol-25"  = "#D4C420",
    "HH119-SI-PP_Fol-26"  = "#E8608A",
    "HH119-SI-PP_Fol-27"  = "#EF19EC",
    "HH119-SI-PP_Fol-28"  = "#FF8C00",
    "HH119-SI-PP_Fol-29"  = "#20B2AA",
    "HH119-SI-PP_Fol-30"  = "#6A5ACD",
    "HH119-SI-PP_Fol-31"  = "#228B22",
    "HH119-SI-PP_Fol-32"  = "#A066DD",
    "HH119-SI-PP_Fol-33"  = "#44BB66",
    "HH119-SI-PP_Fol-34"  = "#C4B410"
  )
