library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(readxl)
library(tidyr)
library(tibble)
library(purrr)

source("10_broad_annotation/script/color_palette.R")
# celltype_colors

combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")
names(combined.BCR.filtered)

names(combined.BCR.filtered) <- names(combined.BCR.filtered) %>%
  str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")

# ------------------------------------------------------------------------------
# Hieu - Making GCtree input file... 
# ------------------------------------------------------------------------------

# containing V J genes, full sequences, CDR3 sequences, abundance

# CDR3 sequence; nucleotide or amino acid sequence?: Both (one in each file)
# V J genes: IGHV1-8 OR Depending on the database (IMGT or 10x reference), we might need to add some additional suffix, e.g. IGHV1-8-1* something, to match the genes in the reference. 
# Full sequence: by full sequence I mean a  sequence of V gene + CDR3 sequence + J gene. This is what we used in our manuscript.  In the single-cell VDJ data, we already have the nucleotide and amino acid sequence for each component. 
# "Since we have both light and heavy chain information, should it be in the same file or in separate?”: I think both are fine, as long as you still keep the “barcode” column, 
# keep the “barcode” column

# Oliver paper: "clonal trees were calculated based on PCs from SI"
# For us that would be: GINA CHECK
# HH117-SILP-INF-PC
# HH117-SILP-nonINF-PC
# HH119-SILP-PC
# HH119-COLP-PC
# And then top 5 clones in each sample?


sample_name <- "HH117-SILP-INF"

# top clone 
top_clones <- combined.BCR.filtered[[sample_name]] %>% 
  group_by(CTstrict) %>% summarise(count = n()) %>% 
  arrange(desc(count)) %>% head(5) %>% pull(CTstrict)

# clone <- "IGH:Cluster.370.IGHV1-8_IGLC:Cluster.1.IGKV1-5"
clone <- top_clones[1]

combined.BCR.clone <- combined.BCR.filtered[[sample_name]] %>% filter(CTstrict == clone)

# Extracting the columns 
# barcode, heavy or light chain, V gene, J gene, CDR3 sequence, full sequences
df_clone <- combined.BCR.clone %>% 
  mutate( 
    IGH_Vgene = IGH %>% str_split_i("\\.", 1),
    IGH_Jgene = IGH %>% str_split_i("\\.", 3),
    IGH_CDR3_nt = cdr3_nt2,
    IGH_CDR3_aa = cdr3_aa2,
    IGLC_Vgene = IGLC %>% str_split_i("\\.", 1),
    IGLC_Jgene = IGLC %>% str_split_i("\\.", 2),
    IGLC_CDR3_nt = cdr3_nt1,
    IGLC_CDR3_aa = cdr3_nt2
  )

# TODO: MAP V and J genes to database (IMGT or 10x reference)

# df_clone_nt: nucleotide seqeunce 
df_clone_nt <- df_clone %>% 
  select(barcode, CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3_nt, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3_nt) %>% 
  rename(IGH_CDR3 = IGH_CDR3_nt, IGLC_CDR3 = IGLC_CDR3_nt) 

# Both chains in same line, adding abundance
df_clone_nt_abundance <- df_clone_nt %>% 
  select(-barcode) %>% 
  group_by(CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3) %>% 
  summarise(count = n()) 

# Heavy chain, adding abundance
df_clone_nt_IGH_abundance <- df_clone_nt %>% 
  select(-barcode) %>% 
  group_by(CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3) %>% 
  summarise(count = n()) 

# Light chain, adding abundance
# df_clone_nt_IGLC_abundance <- df_clone_nt %>% 
#   select(-barcode) %>% 
#   group_by(CTstrict, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3) %>% 
#   summarise(count = n()) 

# Both chains different lines, keeping barcodes (no abundance )
df_clone_nt_long <- df_clone_nt %>% pivot_longer(cols = matches("gene|CDR3")) %>% 
  mutate(chain = str_split_i(name, "_", 1),
         name = str_split_i(name, "_", 2)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  select(-chain)

# Save tables

# ------------------------------------------------------------------------------
# One file per CTstrict and then do top 5 CTstrict for each patient?
# ------------------------------------------------------------------------------

# combined.BCR.filtered_all <- readRDS("20_VDJ/out/combined.BCR.filtered.clean_all.rds")
# # QUESTION FOR GINA 
# # ACROSS ALL SAMPLES OR CHOOSE ONE? - YOU SAID SOMETHING WITH PLASMA CELLS AND THEN MAPPING ON THE FOLLICLES?
# 
# # the abundance of each  sequence in each clone.
# top_clones <- readRDS("40_VDJ_integrated/out/top_5_clones_per_patient.rds")
# 
# top_clones$HH117
# 
# combined.BCR.filtered_all




