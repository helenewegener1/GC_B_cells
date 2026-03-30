library(alakazam)
library(dplyr)
library(ggplot2)
library(shazam)
library(glue)

# Following: 

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")

HH <- "HH119"
HH_spec_clones_vj <- spec_clones_vj[[HH]]

nrow(HH_spec_clones_vj)

# Check that neccessary columns are present
HH_spec_clones_vj$sequence_alignment %>% head()
HH_spec_clones_vj$germline_alignment_d_mask %>% head()

# Subset data
db <- subset(HH_spec_clones_vj)# & sample_clean == "")

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
