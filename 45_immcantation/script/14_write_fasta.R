library(Biostrings) # writing fasta files 
library(msa) # mulitple sequence alignment 
library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")

patients <- names(resolve_LC_list_germlined)

# ------------------------------------------------------------------------------
# Make overview
# ------------------------------------------------------------------------------

# Get top "clone_nr" clone from given sample
# clone_nr <- 1
clone <- HH_spec_clones_vj %>%
  count(clone_subgroup_id, sort = TRUE) %>% 
  slice(clone_nr) %>% 
  pull(clone_subgroup_id)

HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% nrow()
HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence_alignment, sort = TRUE)
HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence, sort = TRUE)
HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(germline_alignment_d_mask, sort = TRUE)

# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence_alignment, sort = TRUE) %>% pull(sequence_alignment)
# HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(germline_alignment_d_mask, sort = TRUE) %>% pull(germline_alignment_d_mask)

# ------------------------------------------------------------------------------
# FASTA file (with abundance) writing
# ------------------------------------------------------------------------------

clone_nrs <- 1:5

for (HH in patients){
  
  # HH <- "HH119"
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 1
    clone <- HH_spec_clones_vj %>%
      count(clone_subgroup_id, sort = TRUE) %>% 
      slice(clone_nr) %>% 
      pull(clone_subgroup_id)

    # Extract sequences of clone, each's abundance and the germline sequence (GL)
    seqs <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence_alignment, sort = TRUE) %>% pull(sequence_alignment)
    seq_abundance <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(sequence_alignment, sort = TRUE) %>% pull(n)
    GL <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(germline_alignment_d_mask, sort = TRUE) %>% pull(germline_alignment_d_mask)
    
    # Make sure there's only one germline sequence
    length(GL) == 1
    
    # Check length of seqeunce to make sure that they are aligned
    map(seqs, nchar) %>% unique()
    nchar(GL)
    
    # Convert sequence into fasta format and name each sequence by their abundance
    seqs_fasta <- DNAStringSet(seqs)
    names(seqs_fasta) <- seq_abundance
    
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
    writeXStringSet(final_fasta, filepath = glue("45_immcantation/fasta/abundance/{HH}_clone_nr_{clone_nr}_clone_{clone}.fasta"))
    
  }

}

# ------------------------------------------------------------------------------
# FASTA file writing
# ------------------------------------------------------------------------------

clone_nrs <- 1:5

for (HH in patients){
  
  # HH <- "HH119"
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 1
    clone <- HH_spec_clones_vj %>%
      count(clone_subgroup_id, sort = TRUE) %>% 
      slice(clone_nr) %>% 
      pull(clone_subgroup_id)
    
    # Extract sequences of clone, each's abundance and the germline sequence (GL)
    seqs <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% pull(sequence_alignment)
    GL <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% count(germline_alignment_d_mask, sort = TRUE) %>% pull(germline_alignment_d_mask)
    
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
    writeXStringSet(final_fasta, filepath = glue("45_immcantation/fasta/non_abundance/{HH}_clone_nr_{clone_nr}_clone_{clone}.fasta"))
    
  }
  
}




