library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Define samples with OCMs
# ------------------------------------------------------------------------------

samples_w_ocm <- c(
  "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB"
)

# 5'OCM color order 
# https://kb.10xgenomics.com/s/article/35179665979149-Gel-Bead-Color-Codes-and-Ordering-for-GEM-X-Universal-Multiplex-a-k-a-On-Chip-Multiplexing-OCM

color_order <- list(
  "Blue" = "OB1",
  "Red" = "OB2",
  "Green" = "OB3", 
  "Yellow" = "OB4"
)

# ------------------------------------------------------------------------------
# Update config files to handle OCMs
# ------------------------------------------------------------------------------

ocm_lines <- c(
  "[samples]",
  "sample_id,ocm_barcode_ids",
  paste(names(color_order), unlist(color_order), sep = ",")
)

for (sample_name in samples_w_ocm) {
  
  # sample_name <- samples_w_ocm[[1]]
  config_file <- glue("04_prep_config/out/multi_config_{sample_name}.csv")
  
  cat(ocm_lines, file = config_file, sep = "\n", append = TRUE)
  
}




