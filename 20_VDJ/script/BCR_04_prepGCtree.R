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
library(readr)

source("10_broad_annotation/script/color_palette.R")
# celltype_colors

# combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")
combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.joined.rds")
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
# And then top 5 clones (CTstrict) in each sample?
# https://www.youtube.com/watch?v=ExOCDX1HLA4
# Gina: Did you oversequence? 10*

sample_name <- "HH117-SILP-INF"

# top clone 
top_clones <- combined.BCR.filtered[[sample_name]] %>% 
  group_by(CTstrict) %>% summarise(count = n()) %>% 
  arrange(desc(count)) %>% head(5) %>% pull(CTstrict)

# clone <- "IGH:Cluster.370.IGHV1-8_IGLC:Cluster.1.IGKV1-5"
clone_nr <- 1
clone <- top_clones[clone_nr]

combined.BCR.clone <- combined.BCR.filtered[[sample_name]] %>% filter(CTstrict == clone)

# Extracting the columns 
# barcode, V gene, J gene, CDR3 sequence, full sequences for heavy and light chain
df_clone <- combined.BCR.clone %>% 
  mutate( 
    IGH_Vgene = IGH %>% str_split_i("\\.", 1),
    IGH_Jgene = IGH %>% str_split_i("\\.", 3),
    IGH_CDR3_nt = cdr3_nt1, # Heavy
    IGH_CDR3_aa = cdr3_aa1,
    IGLC_Vgene = IGLC %>% str_split_i("\\.", 1),
    IGLC_Jgene = IGLC %>% str_split_i("\\.", 2),
    IGLC_CDR3_nt = cdr3_nt2, # Light 
    IGLC_CDR3_aa = cdr3_aa2
  )

# TODO: MAP V and J genes to database (IMGT or 10x reference)?

# ---------------------------------------------------
# Nucleotide seqeunce
# ---------------------------------------------------

# df_clone_nt: nucleotide seqeunce 
df_clone_nt <- df_clone %>% 
  select(barcode, CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3_nt, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3_nt, IGH_full_sequence) %>% 
  rename(IGH_CDR3 = IGH_CDR3_nt, IGLC_CDR3 = IGLC_CDR3_nt) 

# # CDR3 is now in full sequence 
# cdr3 <- df_clone_nt$IGH_CDR3[[1]]
# full <- df_clone_nt$IGH_full_sequence[[1]]
# 
# df_clone_nt$IGH_CDR3 %>% length()
# df_clone_nt$IGH_CDR3 %>% unique() %>% length()
# 
# grep(cdr3, full)

# Both chains in same line --> abundance
df_clone_nt_abundance <- df_clone_nt %>% 
  select(-barcode) %>% 
  summarise(abundance = n(), .by = c(CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3))

write_delim(df_clone_nt_abundance, "20_VDJ/table/BCR_clone_nt_both_abundance.tsv")

# Heavy chain --> abundance
df_clone_nt_IGH_abundance <- df_clone_nt %>% 
  select(-barcode) %>% 
  summarise(abundance = n(), .by = c(CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3)) 

write_delim(df_clone_nt_IGH_abundance, "20_VDJ/table/BCR_clone_nt_IGH_abundance.tsv")

# Light chain --> abundance
df_clone_nt_IGLC_abundance <- df_clone_nt %>%
  select(-barcode) %>%
  summarise(abundance = n(), .by = c(CTstrict, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3))

write_delim(df_clone_nt_IGLC_abundance, "20_VDJ/table/BCR_clone_nt_IGLC_abundance.tsv")

# Both chains different lines, keeping barcodes (no abundance)
df_clone_nt_long <- df_clone_nt %>% 
  pivot_longer(cols = matches("gene|CDR3")) %>% 
  mutate(chain = str_split_i(name, "_", 1),
         name = str_split_i(name, "_", 2)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  select(-chain)

write_delim(df_clone_nt_long, "20_VDJ/table/BCR_df_clone_nt_long.tsv")

# ---------------------------------------------------
# Amino acid seqeunce
# ---------------------------------------------------

# df_clone_aa: amino acid seqeunce 
df_clone_aa <- df_clone %>% 
  select(barcode, CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3_aa, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3_aa) %>% 
  rename(IGH_CDR3 = IGH_CDR3_aa, IGLC_CDR3 = IGLC_CDR3_aa) 

# Both chains in same line --> abundance
df_clone_aa_abundance <- df_clone_aa %>% 
  select(-barcode) %>% 
  summarise(abundance = n(), .by = c(CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3))

write_delim(df_clone_nt_abundance, "20_VDJ/table/BCR_clone_aa_both_abundance.tsv")

# Heavy chain --> abundance
df_clone_aa_IGH_abundance <- df_clone_aa %>% 
  select(-barcode) %>% 
  summarise(abundance = n(), .by = c(CTstrict, IGH_Vgene, IGH_Jgene, IGH_CDR3)) 

write_delim(df_clone_aa_IGH_abundance, "20_VDJ/table/BCR_clone_aa_IGH_abundance.tsv")

# Light chain --> abundance
df_clone_aa_IGLC_abundance <- df_clone_aa %>%
  select(-barcode) %>%
  summarise(abundance = n(), .by = c(CTstrict, IGLC_Vgene, IGLC_Jgene, IGLC_CDR3))

write_delim(df_clone_aa_IGLC_abundance, "20_VDJ/table/BCR_clone_aa_IGLC_abundance.tsv")

# Both chains different lines, keeping barcodes (no abundance)
df_clone_aa_long <- df_clone_aa %>% 
  pivot_longer(cols = matches("gene|CDR3")) %>% 
  mutate(chain = str_split_i(name, "_", 1),
         name = str_split_i(name, "_", 2)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  select(-chain)

write_delim(df_clone_aa_long, "20_VDJ/table/BCR_df_clone_aa_long.tsv")


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

# ------------------------------------------------------------------------------
# FASTA file writing
# ------------------------------------------------------------------------------

library(Biostrings) # writing fasta files 
library(msa) # mulitple sequence alignment 

# assuming df has columns: barcode, CTstrict, sequence_vdj
# clones <- split(df, df$CTstrict)
# 
# for (clone_id in names(clones)) {
#   clone_df <- clones[[clone_id]]
#   seqs <- DNAStringSet(clone_df$sequence_vdj)
#   names(seqs) <- clone_df$barcode
#   writeXStringSet(seqs, filepath = paste0(clone_id, ".fasta"))
# }


# Heavy chain: CTscrict, full seqeunce  abundance
clone_df <- df_clone_nt %>% 
  select(-barcode) %>% 
  summarise(abundance = n(), .by = c(CTstrict, IGH_Vgene, IGH_Jgene, IGH_full_sequence, IGH_CDR3)) 

# Get V and J gene of heavy chain 
v_gene <- df_clone_nt$IGH_Vgene %>% unique()
j_gene <- df_clone_nt$IGH_Jgene %>% unique()

df_clone_nt$CTstrict %>% unique()

# Get sequences
seqs <- DNAStringSet(clone_df$IGH_full_sequence)
# names(seqs) <- glue("Sample:{sample_name}|Clone:{clone_df$CTstrict}|CDR3nt:{clone_df$IGH_CDR3}|Abundance:{clone_df$abundance}")
# names(seqs) <- sprintf("Sample:%s|Clone:%s|CDR3nt:%s|Abundance:%s", 
#                        sample_name, 
#                        clone_df$CTstrict,
#                        clone_df$IGH_CDR3,
#                        clone_df$abundance) 

names(seqs) <- clone_df$abundance

# writeXStringSet(seqs, filepath = glue("20_VDJ/fasta/{sample_name}_clone_{clone_nr}.fasta"))

# Adding germline reference
# reference <- readDNAStringSet("05_run_cellranger/out_v9/res_HH117-SILP-INF-PC/outs/vdj_reference/fasta/regions.fa")
reference <- readDNAStringSet("00_data/Human_sequences.fasta")


## V gene
ref_v_gene <- grep(v_gene, names(reference), value = T)

if (length(ref_v_gene) == 1){
  v_gene_germline <- reference[ref_v_gene]
} else{
  v_gene_germline <- reference[ref_v_gene[1]]
}

## J gene
ref_j_gene <- grep(j_gene, names(reference), value = T)

if (length(ref_j_gene) == 1){
  j_gene_germline <- reference[ref_j_gene]
} else{
  j_gene_germline <- reference[ref_j_gene[1]]
}

## Length of CDR3
cdr3_length <- df_clone_nt$IGH_CDR3 %>% nchar() %>% unique()

if (length(cdr3_length) != 1){
  print("CDR3s are not the same length in this clone type")
}

## Create N mask of CDR3 length
cdr3_mask <- paste(rep("N", cdr3_length), collapse = "")

# Combine 

# ## Convert to character
# v_seq <- as.character(v_gene_germline)
# j_seq <- as.character(j_gene_germline)
# 
# ### Remove leader sequence (L-REGION)
# find_fr1_start <- function(imgt_seq, clone_seqs) {
#   # Use first 20 nt of the clone sequence to find FR1 start
#   first_20 <- substr(clone_seqs[1], 1, 20)
#   pos <- regexpr(first_20, imgt_seq)
#   return(as.integer(pos))
# }
# 
# # Then trim
# trim_v_gene <- function(imgt_seq, clone_seqs) {
#   start_pos <- find_fr1_start(imgt_seq, clone_seqs)
#   if (start_pos == -1) {
#     warning("Could not find FR1 start in IMGT sequence")
#     return(NULL)
#   }
#   return(substr(imgt_seq, start_pos, nchar(imgt_seq)))
# }
# Use it
# v_seq_trimmed <- trim_v_gene(v_seq, clone_df$IGH_full_sequence)
# germline_seq <- paste0(v_seq_trimmed, cdr3_mask, j_seq)


## Combine V and J with the Ns
germline_seq <- paste0(v_gene_germline, cdr3_mask, j_gene_germline)

## Convert back to DNAStringSet
germline_dna <- DNAStringSet(germline_seq)
names(germline_dna) <- "GL"

# Add germline to rest 
seqs_final <- c(seqs, germline_dna)

# Align sequence
## Run alignment (ClustalW, ClustalOmega, or MUSCLE available)
alignment <- msa(seqs_final, method = "ClustalOmega")

## Convert to DNAStringSet for downstream use
alignment_dna <- as(alignment, "DNAStringSet")

# Export as FASTA file
writeXStringSet(alignment_dna, filepath = glue("20_VDJ/fasta/{sample_name}_clone_{clone_nr}_aligned_dna.fasta"))










