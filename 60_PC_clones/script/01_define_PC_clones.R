library(glue)
library(tidyverse)
library(Biostrings) 

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# rds_files <- list.files("45_immcantation/out/rds") 
# resolve_LC_files <- grep("resolve_LC_3_definitions", rds_files, value = TRUE)
# 
# patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
# patients

# ------------------------------------------------------------------------------
# Get data
# ------------------------------------------------------------------------------

# Load both patients
# df_both <- lapply(patients, function(HH) {
#   readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>%
#     filter(
#       locus == "IGH"
#       # !is.na(manual_ADT_full_ID)
#     ) %>%
#     mutate(patient = HH)
# }) %>% setNames(patients)

df_both <- readRDS("45_immcantation/out/rds/resolve_LC_90_similarity_germlined.rds")
patients <- names(df_both)

# df_both$HH117$clone_subgroup_id_90_similarity

# ------------------------------------------------------------------------------
# Define top PC clones 
# ------------------------------------------------------------------------------

PC_clones <- list()

for (HH in patients){
  
  # HH <- "HH117"
  
  df_HH <- df_both[[HH]]
  
  LP_sites <- grep("LP", unique(df_HH$sample_clean_fol), value = TRUE)
  
  for (site in LP_sites){
    
    site_clones <- df_HH %>% 
      filter(
        locus == "IGH", 
        L1_annotation == "PCs", 
        sample_clean_fol == sym(site)
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

for (site in names(PC_clones)){
  
  # site <- "HH117-SILP-INF"
  HH <- site %>% str_split_i("-", 1)
  df_HH <- df_both[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 1
    
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



  

  





