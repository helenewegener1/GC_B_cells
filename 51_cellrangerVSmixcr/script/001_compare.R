library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Define sample
# ------------------------------------------------------------------------------

# sample_name <- "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" # That sample
sample_name <- "HH119-SILP-PC"
# sample_name <- "HH117-SI-MILF-INF-HLADR-AND-CD19"

# ------------------------------------------------------------------------------
# Load data 
# ------------------------------------------------------------------------------

# Per contig
cellr <- read.delim(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/filtered_contig_annotations.csv"), sep = ",") # Cell Ranger
mixcr <- read.delim(glue("50_MiXCR/mixcr_output/BCR-{sample_name}.clones.tsv")) # MiXCR
colnames(mixcr)

# Per cell - AIRR format 
cellr_airr <- read.delim(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/airr_rearrangement.tsv")) 
mixcr_airr <- read.delim(glue("50_MiXCR/mixcr_output/BCR-{sample_name}.airr.tsv")) 

nrow(mixcr_airr)
nrow(cellr_airr)

head(cellr_airr)
colnames(cellr_airr)
colnames(mixcr_airr)

# ------------------------------------------------------------------------------
# Look at genes
# ------------------------------------------------------------------------------

cellr_airr$v_call %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% head()
mixcr_airr$v_call %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% head()

mixcr_airr %>% view()

