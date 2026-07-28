library(dowser)
library(alakazam)
library(dplyr)
library(tidyverse)
library(shazam)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)
patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

resolve_LC_list <- lapply(resolve_LC_files, function(x) readRDS(glue("45_immcantation/out/rds/{x}"))) %>% 
  setNames(patients)

# Look at clone IDs
grep("clone", colnames(resolve_LC_list$HH117), value = TRUE)

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

resolve_LC_list <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  resolve_LC_list[[HH]] %>%
    group_by(clone_subgroup_id_90_similarity, locus) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)
# ------------------------------------------------------------------------------
# Make germline
# ------------------------------------------------------------------------------

resolve_LC_germline_list <- list()

for (HH in patients){

  # HH <- "HH119"

  resolve_LC_HH <- resolve_LC_list[[HH]]

  # Read in IMGT-gapped sequences
  references <- readIMGT(dir = "../packages/immcantation/scripts/germlines/human/vdj") # computerome
  # references <- readIMGT(dir = "00_data/vdj") # local

  # remove germline alignment columns for this example
  db <- select(resolve_LC_HH, -"germline_alignment", -"germline_alignment_d_mask")

  # Reconstruct germline sequences
  resolve_LC_HH_germline <- createGermlines(db, references, nproc=1, clone = "clone_subgroup_id_90_similarity")

  # append to list
  resolve_LC_germline_list[[HH]] <- resolve_LC_HH_germline

}

saveRDS(resolve_LC_germline_list, "45_immcantation/out/rds/06_resolve_LC_germlined.rds")

