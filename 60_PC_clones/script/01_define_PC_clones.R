library(glue)
library(tidyverse)
library(Biostrings) 

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

df_both <- readRDS("45_immcantation/out/rds/06_resolve_LC_germlined.rds")
patients <- names(df_both)

large_clone <- df_both$HH119 %>% 
  dplyr::count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  head(1) %>% 
  pull(clone_subgroup_id_90_similarity)

# ------------------------------------------------------------------------------
# Define top PC clones (excluding clones with no sequence variation)
# ------------------------------------------------------------------------------

PC_clones <- list()

for (HH in patients){
  
  # HH <- "HH119"
  
  df_HH <- df_both[[HH]]
  
  LP_sites <- grep("LP", unique(df_HH$sample_clean), value = TRUE)
  
  for (site in LP_sites){
    
    # site <- "HH119-SILP"
    
    # clones with more than one unique sequence -> can build a tree
    clones_with_variation <- df_HH %>% 
      filter(locus == "IGH", L1_annotation == "PCs", sample_clean == site, clone_subgroup_id_90_similarity != large_clone) %>% 
      group_by(clone_subgroup_id_90_similarity) %>% 
      summarise(n_unique_seqs = n_distinct(sequence_alignment), .groups = "drop") %>% 
      filter(n_unique_seqs > 1) %>% 
      pull(clone_subgroup_id_90_similarity)
    
    site_clones <- df_HH %>% 
      filter(
        locus == "IGH", 
        L1_annotation == "PCs", 
        sample_clean == site,
        clone_subgroup_id_90_similarity %in% clones_with_variation
      ) %>% 
      dplyr::count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
      head(10) %>% 
      pull(clone_subgroup_id_90_similarity)
    
    PC_clones[[site]] <- site_clones
  }
  
}

PC_clones

# ------------------------------------------------------------------------------
# Export top clones in fasta files
# ------------------------------------------------------------------------------

# Prep export of fasta files
outdir <- glue("60_PC_clones/fasta")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clone_nrs <- 1:10

# prep seq name dir 
seq_dir <- list()

for (site in names(PC_clones)){
  
  # site <- "HH119-SILP"
  HH <- site %>% str_split_i("-", 1)
  df_HH <- df_both[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 7
    
    # rank those clones by total size (all cell types) and take top 10
    clone <- PC_clones[[site]][[clone_nr]]
    
    # Extract sequences of clone, each's abundance and the germline sequence (GL)
    seqs <- df_HH %>% filter(clone_subgroup_id_90_similarity == clone & locus == "IGH") %>% pull(sequence_alignment)
    GL <- df_HH %>% filter(clone_subgroup_id_90_similarity == clone & locus == "IGH") %>% dplyr::count(germline_alignment_d_mask, sort = TRUE) %>% pull(germline_alignment_d_mask)
    
    # PHYLIP format no longer allows dots in sequence (ValueError when running gctree duplicate)
    seqs <- seqs %>% str_replace_all("\\.", "-")
    GL <- GL %>% str_replace_all("\\.", "-")
    GL <- GL %>% str_replace_all("N", "-")
    
    # Make sure there's only one germline sequence
    if (length(GL) != 1){
      print(glue("NB! GL len != 1 in {site}"))
    }
   
    # Check length of seqeunce to make sure that they are aligned
    map(seqs, nchar) %>% unique()
    nchar(GL)
    
    # Convert sequence into fasta format and name each sequence by their abundance
    seqs_fasta <- DNAStringSet(seqs)
    names(seqs_fasta) <- paste0("sequence_", 1:length(seqs_fasta))
    
    # Export seq_dir
    seq_dir_tmp <- seqs
    seq_dir_tmp <- seq_dir_tmp %>% str_replace_all("-", "\\.")
    names(seq_dir_tmp) <- names(seqs_fasta)
    seq_dir[[glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")]] <- seq_dir_tmp
    
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
    writeXStringSet(final_fasta, filepath = glue("{outdir}/{site}_clone_nr_{clone_nr}_clone_{clone}.fasta"))
    
  }
  
}

saveRDS(seq_dir, "60_PC_clones/out/rds/seq_dir.rds")

  

  





