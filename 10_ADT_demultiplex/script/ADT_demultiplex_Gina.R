getwd()

library(SeuratObject)
library(DropletUtils)
library(Seurat)
library(tidyverse)
library(glue)
library(patchwork)
library(readxl)
# library(ztable)
library(pheatmap)

# Load data your data
# seurat_obj <-  


# CLR normalization of ADT
seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")

# Demultiplexing with HTODemux and MULTIseqDemux
seurat_obj <- HTODemux(seurat_obj, assay = "ADT", positive.quantile = 0.99)
seurat_obj <- MULTIseqDemux(seurat_obj, assay = "ADT", autoThresh = TRUE)

# seurat_obj$MULTI_ID
# seurat_obj$ADT_classification 

# Mutate columns so result from HTODemux and MULTIseqDemux become more comparable 
seurat_obj@meta.data$ADT_ID <- lapply(seurat_obj@meta.data$ADT_classification, function(x) {ifelse(str_detect(x, "_"), "Doublet", x)}) %>% unlist()
seurat_obj$MULTI_classification.global <- seurat_obj$MULTI_ID %>% str_replace("Fol-\\d+", "Singlet")

# Look at results
seurat_obj$ADT_ID %>% table()
seurat_obj$MULTI_ID %>% table()

seurat_obj$ADT_classification.global %>% table()
seurat_obj$MULTI_classification.global %>% table()


###################### CONSENSUS OF HTODemux AND HTODemux ######################

# Compare result of methods 
table(seurat_obj$ADT_ID, seurat_obj$MULTI_ID)

# Consensus data used to plot and extract consensus column for final seurat object 
df_consens <- seurat_obj@meta.data %>% 
  select(ADT_classification, ADT_classification.global, ADT_ID, ADT_maxID, ADT_secondID, MULTI_ID, MULTI_classification, MULTI_classification.global) %>% 
  # Extract MULTI_maxID and MULTI_secondID to compare to HTODemux annotations 
  mutate(
    MULTI_maxID = str_split_i(MULTI_classification, "_", 1),
    MULTI_secondID = str_split_i(MULTI_classification, "_", 2)
  ) %>% 
  mutate(
    
    # Everything same, also same exact doublets, else NA
    ADT_consensus_very_hard = ifelse(ADT_classification == MULTI_classification, ADT_classification, NA),
    
    # Same singlet, negative and doublet detection, else NA
    ADT_consensus_hard = ifelse(ADT_ID == MULTI_ID, ADT_ID, NA),
    
    # Since HTODemux detects many doublets and MULTIseqDemux detects more negatives. 
    ADT_consensus_medium = case_when(
      # if HTODemux and MULTIseqDemux agree - keep label
      ADT_ID == MULTI_ID ~ ADT_ID,
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
      ADT_classification.global == "Singlet" & MULTI_classification.global == "Singlet" & ADT_classification != MULTI_maxID ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet IF ADT_maxID == MULTI_ID
      ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet" & ADT_maxID == MULTI_ID ~ MULTI_ID,
      # else if ADT_maxID != MULTI_ID
      ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet" & ADT_maxID != MULTI_ID ~ "Doublet",
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet IF ADT_classification == MULTI_maxID
      # MULTI_maxID will be Negative if MULTI_ID is Negative 
      # ADT_classification.global == "Singlet" & MULTI_ID == "Negative" & ADT_classification == MULTI_maxID ~ ADT_classification,
      # else if ADT_maxID != MULTI_ID
      # ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ "Negative",
      ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ ADT_classification,
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
      ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical but the second ones aren't.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID != MULTI_secondID ~ ADT_maxID,
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet if the two first and the second follicles are identical.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID == MULTI_secondID ~ "Doublet",
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"

    ),
    ADT_consensus_soft = case_when(
      # if HTODemux and MULTIseqDemux agree - keep label
      ADT_ID == MULTI_ID ~ ADT_ID,
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
      ADT_classification.global == "Singlet" & MULTI_classification.global == "Singlet" & ADT_classification != MULTI_maxID ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet
      ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet"  ~ MULTI_ID,
      # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet
      ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ ADT_classification,
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
      ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID  ~ ADT_maxID,
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"
      
    )
  )

head(df_consens) 

# Add to meta data
df_consens$ADT_consensus_hard.global <- df_consens$ADT_consensus_hard %>% str_replace("Fol-\\d+", "Singlet") 
df_consens$ADT_consensus_medium.global <- df_consens$ADT_consensus_medium %>% str_replace("Fol-\\d+", "Singlet") 
df_consens$ADT_consensus_soft.global <- df_consens$ADT_consensus_soft %>% str_replace("Fol-\\d+", "Singlet")

#################### Add consensus colums to seurat object #####################

# Check that same number of rows 
seurat_obj@meta.data %>% dim()
df_consens %>% dim()

seurat_obj@meta.data$ADT_consensus_very_hard <- df_consens$ADT_consensus_very_hard
seurat_obj@meta.data$ADT_consensus_hard <- df_consens$ADT_consensus_hard
seurat_obj@meta.data$ADT_consensus_medium <- df_consens$ADT_consensus_medium
seurat_obj@meta.data$ADT_consensus_soft<- df_consens$ADT_consensus_soft

seurat_obj@meta.data$ADT_consensus_hard.global<- df_consens$ADT_consensus_hard.global
seurat_obj@meta.data$ADT_consensus_medium.global<- df_consens$ADT_consensus_medium.global
seurat_obj@meta.data$ADT_consensus_soft.global <- df_consens$ADT_consensus_soft.global

# Look at results 
table(seurat_obj@meta.data$ADT_consensus_hard.global, useNA = "always")
table(seurat_obj@meta.data$ADT_consensus_medium.global, useNA = "always")
table(seurat_obj@meta.data$ADT_consensus_soft.global, useNA = "always")

# Look at results per cell type (GINA USE THE COLUMN YOU ANNOTATED IN)
table(seurat_obj@meta.data$ADT_consensus_hard.global, seurat_obj@meta.data$GINA_CELLTYPE, useNA = "ifany")
table(seurat_obj@meta.data$ADT_consensus_medium.global, seurat_obj@meta.data$GINA_CELLTYPE, useNA = "ifany")
table(seurat_obj@meta.data$ADT_consensus_soft.global, seurat_obj@meta.data$GINA_CELLTYPE, useNA = "ifany")


