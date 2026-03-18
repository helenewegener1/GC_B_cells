library(tidyverse)
library(airr)
library(scoper)
library(shazam)

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

airr_files <- list.files("05_run_cellranger/airr_files")

sample_names <- airr_files %>% str_remove_all("_airr_rearrangement.tsv")

test <- read_delim("05_run_cellranger/airr_files/HH117-SI-MILF-INF-HLADR-AND-CD19_airr_rearrangement.tsv")

# Load and combine all AIRR files, adding sample and patient info
base_path <- "05_run_cellranger/airr_files/{sample_name}_airr_rearrangement.tsv"

airr_data <- lapply(sample_names, function(s) {
  path <- gsub("\\{sample_name\\}", s, base_path)
  db <- read_rearrangement(path)
  db$sample_id <- s
  db$patient_id <- s %>% str_split_i("-", 1)
  return(db)
}) %>% bind_rows()

# Add locus information 
airr_data <- airr_data %>%
  mutate(locus = case_when(
    str_detect(v_call, "^IGHV") ~ "IGH",
    str_detect(v_call, "^IGKV") ~ "IGK",
    str_detect(v_call, "^IGLV") ~ "IGL",
    TRUE ~ NA_character_
  ))

# Removing cells with multiple heavy chains from the single cell data
# tmp solution, let's look at UMIs for the chains...
heavy_count <- table(filter(airr_data, locus=="IGH")$cell_id)
multi_heavy_cells <- names(heavy_count)[heavy_count > 1]
airr_data <- filter(airr_data, !cell_id %in% multi_heavy_cells)

# ------------------------------------------------------------------------------
# Identifying clones by hierarchical clustering
# ------------------------------------------------------------------------------

# Split by patient and run SCOPer steps per patient
results <- lapply(split(airr_data, airr_data$patient_id), function(patient_db) {
  
  # patient_db <- split(airr_data, airr_data$patient_id)$HH119

  # Use nucleotide Hamming distance and normalize by junction length
  patient_db <- distToNearest(
    patient_db, 
    sequenceColumn="junction", 
    vCallColumn="v_call",
    jCallColumn="j_call", 
    model="ham", # Hamming
    normalize="len", 
    nproc=1,  
    cellIdColumn = "cell_id" #important for single-cell
  )
  
  # Threshold estimation
  threshold_output <- findThreshold(patient_db$dist_nearest, method = "gmm")
  threshold <- threshold_output@threshold
  
  # Clean 
  patient_db_clean <- patient_db  %>% select(-c(clone_id, cell_id))

  # Hierarchical cloning
  hier <- hierarchicalClones(
    patient_db_clean,
    threshold = threshold,
    junction = "junction",
    v_call = "v_call",
    j_call = "j_call", 
    summarize_clones = TRUE
  )

  ## Plot a histogram of inter and intra clonal distances
  plot(hier, binwidth=0.02)
  
  # Spectral cloning
  spectral <- spectralClones(
    patient_db_clean,
    method = "novj",
    germline = "germline_alignment_d_mask",
    junction = "junction",
    v_call = "v_call",
    j_call = "j_call",
    summarize_clones = TRUE
  )
  
  plot(spectral, binwidth=0.02)
  
  return(spectral)
  
})

# Combine results back together if needed
final_db <- bind_rows(results)

