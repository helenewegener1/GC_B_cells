library(Biostrings) # writing fasta files 
# library(msa) # mulitple sequence alignment 
library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/06_resolve_LC_germlined.rds")
# annotation <- readRDS("00_data/GCB_meta_GL_new.rds")
# 
# table(resolve_LC_list_germlined$HH117$locus, resolve_LC_list_germlined$HH117$sample_clean_fol)
# 
# patients <- names(resolve_LC_list_germlined)
# 
# # Prep outdir
outdir_fasta <- glue("50_BACH2/fasta/")
dir.create(outdir_fasta, recursive = TRUE, showWarnings = FALSE)
# 
# # Not germline resolved
# rds_files <- list.files("45_immcantation/out/rds") 
# resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)
# patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
# patients
# 
# # Add finer annotation 
# annotation$cell_id <- rownames(annotation) %>% str_remove("_\\d+")
# annotation$cell_id <- glue("{annotation$Sample_name}_{annotation$cell_id}") %>% str_remove("_(Blue|Red|Grenn|Yellow)")
# annotation_clean <- annotation %>% select(cell_id, L3_GCB_annotation)
# 
# # table(annotation$cell_id %in% resolve_LC_list_germlined[["HH153"]]$cell_id)
# 
# resolve_LC_list_germlined <- lapply(patients, function(HH){
#   
#   clones <- resolve_LC_list_germlined[[HH]] %>% 
#     left_join(annotation_clean, by = "cell_id") %>% 
#     mutate(
#       L1_L3_annotation = case_when(
#         !is.na(L3_GCB_annotation) ~ L3_GCB_annotation, 
#         .default = L1_annotation
#       )
#     )
#   
#   clones
#   
# }) %>% 
#   setNames(patients)
# 
# resolve_LC_list_germlined$HH117 %>% filter(locus == "IGH") %>% pull(L1_L3_annotation) %>% table(useNA = "always")
# 
# saveRDS(resolve_LC_list_germlined, "50_BACH2/out/resolve_LC_list_germlined_L1_L3.rds")
resolve_LC_list_germlined <- readRDS("50_BACH2/out/resolve_LC_list_germlined_L1_L3.rds")
patients <- names(resolve_LC_list_germlined)

# ------------------------------------------------------------------------------
# Clones with BACH2+ and/or BACH2_cycling cells
# ------------------------------------------------------------------------------

BACH2_clones <- lapply(patients, function(HH){
  
  # HH <- "HH153"
  clones <- resolve_LC_list_germlined[[HH]] %>% 
    filter(locus == "IGH") %>% 
    count(clone_subgroup_id_90_similarity, L1_L3_annotation) %>% 
    filter(L1_L3_annotation %in% c("BACH2+", "BACH2_cycling")) %>% 
    pull(clone_subgroup_id_90_similarity) %>% 
    unique()
  
  # Filter "against" singletons 
  clones_clean <- resolve_LC_list_germlined[[HH]] %>% 
    filter(
      locus == "IGH",
      clone_subgroup_id_90_similarity %in% clones
    ) %>% 
    count(clone_subgroup_id_90_similarity, sequence_alignment, sort = TRUE) %>% 
    filter(n>1) %>% 
    pull(clone_subgroup_id_90_similarity) %>% 
    unique()
  
  
  # Exclude large clone 
  if (HH == "HH119"){
    clones_clean <- clones_clean[clones_clean != "21854_1"]
  } 
  
  clones_clean 
  
}) %>% 
  setNames(patients)

lengths(BACH2_clones)

# ------------------------------------------------------------------------------
# FASTA file writing
# ------------------------------------------------------------------------------

# prep seq name dir 
seq_dir <- list()
seq_cell_id_dir <- list()

for (HH in patients){
  
  # HH <- "HH117"
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  clone_nrs <- 1:length(BACH2_clones[[HH]])
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 3

    # rank those clones by total size (all cell types) and take top 10
    clone <- BACH2_clones[[HH]][[clone_nr]]
    
    # Extract sequences of clone, each's abundance and the germline sequence (GL)
    seqs <- HH_spec_clones_vj %>% filter(clone_subgroup_id_90_similarity == clone & locus == "IGH") %>% pull(sequence_alignment)
    seqs_cell_id <- HH_spec_clones_vj %>% filter(clone_subgroup_id_90_similarity == clone & locus == "IGH") %>% pull(cell_id)
    GL <- HH_spec_clones_vj %>% filter(clone_subgroup_id_90_similarity == clone & locus == "IGH") %>% 
      dplyr::count(germline_alignment_d_mask, sort = TRUE) %>% pull(germline_alignment_d_mask)
    
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
    
    # Export seq_dir
    seq_dir_tmp <- seqs
    seq_dir_tmp <- seq_dir_tmp %>% str_replace_all("-", "\\.")
    names(seq_dir_tmp) <- names(seqs_fasta)
    seq_dir[[glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")]] <- seq_dir_tmp
    
    # Export seq name to cell_id map 
    names(seqs_cell_id) <- names(seqs_fasta)
    seq_cell_id_dir[[glue("{HH}_clone_nr_{clone_nr}_clone_{clone}")]] <- seqs_cell_id
    
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
    writeXStringSet(final_fasta, filepath = glue("{outdir_fasta}/{HH}_clone_nr_{clone_nr}_clone_{clone}.fasta"))

    
  }
  
}

# Prep outdir
outdir_rds <- glue("50_BACH2/out/rds/")
dir.create(outdir_rds, recursive = TRUE, showWarnings = FALSE)

saveRDS(seq_dir, glue("{outdir_rds}/seq_dir.rds"))
saveRDS(seq_cell_id_dir, glue("{outdir_rds}/seq_cell_id_dir.rds"))





