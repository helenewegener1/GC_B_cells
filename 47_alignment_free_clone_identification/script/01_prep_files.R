library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

bcr_data <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")

patients <- names(bcr_data)

# Wrangle data
bcr_data <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  bcr_data_HH <- bcr_data[[HH]]
  
  bcr_data_HH <- bcr_data_HH %>%
    mutate(
      sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
    )
  
  return(bcr_data_HH)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------


lapply(patients, function(HH){
  
  # Define data
  # HH <- "HH117"
  bcr_data_HH <- bcr_data[[HH]]
  
  # Raw sequnces contain constant region 
  bcr_data_HH$sequence %>% nchar() %>% range()
  
  # Trim raw seqeunce using v_sequence_start and j_sequence_end.
  repertoire <- bcr_data_HH %>% 
    mutate(
      sequence_trimmed = str_sub(sequence, v_sequence_start, j_sequence_end)
    ) %>% 
    filter(!is.na(sequence_trimmed), sequence_trimmed != "") %>% 
    select(
      ID       = cell_id,
      SEQUENCE = sequence_trimmed
    )
  
  head(repertoire)
  
  # Stats
  bcr_data_HH$sequence %>% nchar() %>% range()
  repertoire$SEQUENCE %>% nchar() %>% range()
  
  bcr_data_HH$sequence %>% nchar() %>% hist()
  repertoire$SEQUENCE %>% nchar() %>% hist()
  
  # Make negation sequences: Use sequences from other patients as negation sequences
  table_neg <- bcr_data[names(bcr_data) != HH] %>%
    bind_rows() %>%
    mutate(
      sequence_trimmed = str_sub(sequence, v_sequence_start, j_sequence_end)
    ) %>% 
    filter(!is.na(sequence_trimmed), sequence_trimmed != "") %>% 
    select(
      SEQUENCE = sequence_trimmed
    ) %>% 
    # Sample a reasonable number - doesn't need to be all of them
    slice_sample(n = nrow(repertoire))
  
  write_tsv(repertoire, glue("47_alignment_free_clone_identification/out/{HH}_repertoire.tsv"))
  write_tsv(table_neg,  glue("47_alignment_free_clone_identification/out/{HH}_negation.tsv"))

})



