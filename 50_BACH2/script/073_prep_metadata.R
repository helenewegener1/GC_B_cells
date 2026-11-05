library(Biostrings) # writing fasta files 
# library(msa) # mulitple sequence alignment 
library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list_germlined <- readRDS("50_BACH2/out/resolve_LC_list_germlined_L1_L3.rds")
dir.create(glue("50_BACH2/gctree_meta"), recursive = TRUE)

patients <- names(resolve_LC_list_germlined)

fasta_path <- glue("50_BACH2/fasta/")
fasta_files <- list.files(fasta_path)

# resolve_LC_list_germlined$HH117$clone_subgroup_id_90_similarity

# Load finer annotations
# L3_GCB_annotation <- readRDS("00_data/GCB_meta_GL.rds")
L3_GCB_annotation <- readRDS("00_data/GCB_meta_GL_new.rds")
L3_GCB_annotation$cell_id <- glue("{L3_GCB_annotation$Sample_name}_{rownames(L3_GCB_annotation)}") %>% str_remove("_\\d+")

# Clean it up 
L3_GCB_annotation_clean <- L3_GCB_annotation %>% select(cell_id, L3_GCB_annotation)

# ------------------------------------------------------------------------------
# Clones with BACH2+ and/or BACH2_cycling cells
# ------------------------------------------------------------------------------

BACH2_clones <- lapply(patients, function(HH){
  
  # HH <- "HH117"
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
# Meta data file writing
# ------------------------------------------------------------------------------

for (HH in patients){
  
  # HH <- "HH117"
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  clone_nrs <- 1:length(BACH2_clones[[HH]])
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 1
    filename <- grep(glue("{HH}_clone_nr_{clone_nr}_"), fasta_files, value = TRUE)
    fasta <- readDNAStringSet(filepath = glue("{fasta_path}/{filename}"))
    clone <- str_extract(filename, "\\d+_\\d+(?=\\.fasta)")
    
    # Extract metadata
    seqs_meta <- HH_spec_clones_vj %>% 
      filter(clone_subgroup_id_90_similarity == clone & locus == "IGH") %>% 
      select(L1_annotation, L1_L3_annotation, c_call_grouped, sample_clean_fol)
    
    # Map seq_names on meta data
    seq_names <- names(fasta)[1:length(fasta)-1]
    length(seq_names) == nrow(seqs_meta)
    seqs_meta$seq_name <- seq_names
    
    # Read idmap.txt file
    sample <- str_split_i(filename, "\\.", 1)
    idmap <- read.csv(glue("50_BACH2/out/{sample}/idmap.txt"), header = FALSE, col.names = c("seq_unique", "seq_name")) 
    
    # Remove GL
    idmap <- idmap %>% filter(seq_unique != "GL")
    
    # Map metadata on seq_unique
    gctree_meta <- idmap %>%
      separate_rows(seq_name, sep = ":") %>%
      left_join(seqs_meta, by = "seq_name") %>% 
      group_by(seq_unique) %>%
      summarise(
        seq_name = paste(seq_name, collapse = ":"),
        L1_annotation = paste(unique(L1_annotation), collapse = ":"),
        L1_L3_annotation = paste(unique(L1_L3_annotation), collapse = ":"),
        # c_call = paste(unique(c_call), collapse = ":"),
        c_call_grouped = paste(unique(c_call_grouped), collapse = ":"),
        sample_clean_fol = paste(unique(sample_clean_fol), collapse = ":")
      ) %>% 
      mutate(
        L1_annotation_int = as.integer(factor(L1_annotation)),
        L1_L3_annotation_int = as.integer(factor(L1_L3_annotation)),
        c_call_grouped_int = as.integer(factor(c_call_grouped)),
        sample_clean_fol_int = as.integer(factor(sample_clean_fol))
      )
    
    # gctree_meta$L1_annotation %>% table()
    # gctree_meta$c_call_grouped %>% table()
    # gctree_meta$sample_clean_fol %>% table()
    
    # gctree_meta$L1_annotation_int %>% table()
    # gctree_meta$c_call_grouped %>% table()
    # gctree_meta$sample_clean_fol %>% table()
  
    
    write.csv(gctree_meta, glue("50_BACH2/gctree_meta/{sample}_gctree_meta.txt"), row.names = FALSE, quote = FALSE)
    
    # L1_counts.csv
    L1_counts <- gctree_meta %>%
      select(seq_unique, seq_name) %>%       
      separate_rows(seq_name, sep = ":") %>% 
      left_join(seqs_meta %>% select(seq_name, L1_annotation), by = "seq_name") %>%  
      group_by(seq_unique, L1_annotation) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = L1_annotation, values_from = count, values_fill = 0)

    write.csv(L1_counts, glue("50_BACH2/gctree_meta/{sample}_L1_counts.csv"), row.names = FALSE)
    
    # L3_counts.csv
    L3_counts <- gctree_meta %>%
      select(seq_unique, seq_name) %>%       
      separate_rows(seq_name, sep = ":") %>% 
      left_join(seqs_meta %>% select(seq_name, L1_L3_annotation), by = "seq_name") %>%  
      group_by(seq_unique, L1_L3_annotation) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = L1_L3_annotation, values_from = count, values_fill = 0)
    
    write.csv(L1_counts, glue("50_BACH2/gctree_meta/{sample}_L1_L3_counts.csv"), row.names = FALSE)
    
    # isotype_counts.csv
    isotype_counts <- gctree_meta %>%
      select(seq_unique, seq_name) %>%       
      separate_rows(seq_name, sep = ":") %>% 
      left_join(seqs_meta %>% select(seq_name, c_call_grouped), by = "seq_name") %>%  
      group_by(seq_unique, c_call_grouped) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = c_call_grouped, values_from = count, values_fill = 0)
    
    write.csv(isotype_counts, glue("50_BACH2/gctree_meta/{sample}_isotype_counts.csv"), row.names = FALSE)
    
    # sample_clean_fol_counts.csv
    sample_clean_fol_counts <- gctree_meta %>%
      select(seq_unique, seq_name) %>%       
      separate_rows(seq_name, sep = ":") %>% 
      left_join(seqs_meta %>% select(seq_name, sample_clean_fol), by = "seq_name") %>%  
      group_by(seq_unique, sample_clean_fol) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = sample_clean_fol, values_from = count, values_fill = 0)

    write.csv(sample_clean_fol_counts, glue("50_BACH2/gctree_meta/{sample}_sample_clean_fol_counts.csv"), row.names = FALSE)

  }
  
}

