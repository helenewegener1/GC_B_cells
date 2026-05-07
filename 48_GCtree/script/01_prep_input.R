library(Biostrings) # writing fasta files 
# library(msa) # mulitple sequence alignment 
library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")
resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_list_gmm_threshold_germlined.rds")

version <- "gmm_threshold_GC_clones"

patients <- names(resolve_LC_list_germlined)

# ------------------------------------------------------------------------------
# Make overview
# ------------------------------------------------------------------------------

# Get top "clone_nr" clone from given sample
# clone_nr <- 1
# clone <- HH_spec_clones_vj %>%
#   count(clone_subgroup_id, sort = TRUE) %>%
#   slice(clone_nr) %>%
#   pull(clone_subgroup_id)
# 
# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% nrow()
# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence_alignment, sort = TRUE)
# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence, sort = TRUE)
# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(germline_alignment_d_mask, sort = TRUE)

# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence_alignment, sort = TRUE) %>% pull(sequence_alignment)
# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(germline_alignment_d_mask, sort = TRUE) %>% pull(germline_alignment_d_mask)

# ------------------------------------------------------------------------------
# FASTA file writing
# ------------------------------------------------------------------------------

# Define top GC clone
top_GC_clones <- lapply(patients, function(HH) {
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- resolve_LC_list[[HH]] %>%
    filter(locus == "IGH") %>% 
    filter(L1_annotation == "GC_B_cells") %>%
    count(clone_subgroup_id, sort = TRUE) %>% 
    head(10) %>% 
    pull(clone_subgroup_id)
  
}) %>% setNames(patients)

clone_nrs <- 1:10

for (HH in patients){
  
  # HH <- "HH117"
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 1

    # rank those clones by total size (all cell types) and take top 10
    clone <- top_GC_clones[[HH]][[clone_nr]]
    
    # Extract sequences of clone, each's abundance and the germline sequence (GL)
    seqs <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% pull(sequence_alignment)
    GL <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(germline_alignment_d_mask, sort = TRUE) %>% pull(germline_alignment_d_mask)
    
    # PHYLIP format no longer allows dots in sequence (ValueError when running gctree duplicate)
    seqs <- seqs %>% str_replace_all("\\.", "-")
    GL <- GL %>% str_replace_all("\\.", "-")
    GL <- GL %>% str_replace_all("N", "-")
    
    # Make sure there's only one germline sequence
    length(GL) == 1
    
    # Check length of seqeunce to make sure that they are aligned
    map(seqs, nchar) %>% unique()
    nchar(GL)
    
    # Convert sequence into fasta format and name each sequence by their abundance
    seqs_fasta <- DNAStringSet(seqs)
    names(seqs_fasta) <- paste0("sequence_", 1:length(seqs_fasta))
    
    # Convert germline sequence into fasta format and name the sequence "GL"
    GL_fasta <- DNAStringSet(GL)
    names(GL_fasta) <- "GL"
    
    # Add germline to rest 
    final_fasta <- c(seqs_fasta, GL_fasta)
    
    # Align sequence - Skip this since sequences are already aligned
    ## Run alignment (ClustalW, ClustalOmega, or MUSCLE available)
    # alignment <- msa(final_fasta, method = "ClustalOmega")
    ## Convert to DNAStringSet for downstream use
    # alignment_dna <- as(alignment, "DNAStringSet")
    
    # Export as FASTA file
    outdir <- glue("48_GCtree/fasta/{version}")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    writeXStringSet(final_fasta, filepath = glue("{outdir}/{HH}_clone_nr_{clone_nr}_clone_{clone}.fasta"))
    
  }
  
}


