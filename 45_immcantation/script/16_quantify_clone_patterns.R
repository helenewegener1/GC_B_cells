library(glue)
library(tidyverse)
library(grid)
source("10_broad_annotation/script/color_palette.R")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

# resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")
# plot_version <- "plot"

resolve_LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list_gmm_threshold_germlined.rds")
plot_version <- "plot_gmm_threshold"

patients <- names(resolve_LC_list)

ncol(seurat_integrated)

# df <- readRDS("45_immcantation/out/rds/05_spec_clones_vj_heavy.rds")
# 
# df$HH119 %>% count(clone_id, sort = TRUE)
# 
# resolve_LC_list$HH117 %>% 
#   filter(locus == "IGH") %>% 
#   count(clone_id, sort = TRUE)

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

resolve_LC_list <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  resolve_LC_list[[HH]] %>%
    group_by(clone_subgroup_id, locus) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)

# top 10 GC clones for each patient 
all_clones <- list(
  "HH117" = c(
    "4221_1", "2628_1", "1849_1", "3709_1", "2301_1",
    "1320_1", "5941_1", "6115_1", "1910_1", "2169_1"
  ),
  
  "HH119" = c(
    "28075_1", "12120_1", "15287_1", "23124_1",
    "8372_1", "3869_1", "25158_1", "7913_1"
  )
)

# ------------------------------------------------------------------------------
# Investigate 1 clone
# ------------------------------------------------------------------------------

HH <- "HH117"
clone_nr <- 3

# Define clone for investigation 
clone <- all_clones[[HH]][[clone_nr]]

# Filter for defined clone and extract relevant meta data 
subset <- resolve_LC_list[[HH]] %>% 
  filter(
    locus == "IGH", 
    clone_subgroup_id == clone
  ) %>% 
  select(sample_clean_fol, L1_annotation, c_call)

# Have a look 
subset %>% count(sample_clean_fol, L1_annotation, c_call)






