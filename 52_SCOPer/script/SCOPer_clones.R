library(tidyverse)
library(airr)
library(scoper)

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

airr_files <- list.files("05_run_cellranger/airr_files")

sample_names <- airr_files %>% str_remove_all("_airr_rearrangement.tsv")


# Load and combine all AIRR files, adding sample and patient info
base_path <- "05_run_cellranger/airr_files/{sample_name}_airr_rearrangement.tsv"

airr_data <- lapply(sample_names, function(s) {
  path <- gsub("\\{sample_name\\}", s, base_path)
  db <- read_rearrangement(path)
  db$sample_id <- s
  db$patient_id <- s %>% str_split_i("-", 1)
  return(db)
}) %>% bind_rows()

# ------------------------------------------------------------------------------
# SCOPer
# ------------------------------------------------------------------------------

# Split by patient and run SCOPer steps per patient
results <- lapply(split(airr_data, airr_data$patient_id), function(patient_db) {
  
  # Step 1: Threshold estimation
  threshold <- findThreshold(patient_db$junction, method = "density")
  
  # Step 2: Hierarchical cloning
  hier <- hierarchicalClones(patient_db,
                             threshold = threshold,
                             junction = "junction",
                             v_call = "v_call",
                             j_call = "j_call")
  
  # Step 3: Spectral cloning
  spectral <- spectralClones(hier,
                             method = "novj",
                             germline = "germline_alignment_d_mask",
                             junction = "junction",
                             v_call = "v_call",
                             j_call = "j_call")
  
  return(spectral)
  
})

# Combine results back together if needed
final_db <- bind_rows(results)
