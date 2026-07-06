library(Biostrings) # writing fasta files 
library(msa) # mulitple sequence alignment 
library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")
# version <- ""

# resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_list_gmm_threshold_germlined.rds")
resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_90_similarity_germlined.rds")
version <- "90_similarity"
dir.create(glue("48_GCtree/gctree_meta_{version}"), recursive = TRUE)

patients <- names(resolve_LC_list_germlined)

# fasta_path <- "48_GCtree/fasta/GC_clones/"
fasta_path <- glue("48_GCtree/fasta/{version}")
fasta_files <- list.files(fasta_path)

# resolve_LC_list_germlined$HH117$clone_subgroup_id_90_similarity

# ------------------------------------------------------------------------------
# Meta data file writing
# ------------------------------------------------------------------------------

# L1_int_mapping <- c(
#   "Tfh_cells"                  = 1,
#   "Naive_Bcells"               = 2,
#   "Memory_Bcells"              = 3,
#   "GC_B_cells"                 = 4,
#   "PCs"                        = 5,
#   "Unconventional_Bcells"      = 6
# )

clone_nrs <- 2:10

for (HH in patients){
  
  # HH <- "HH117"
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 1
    filename <- grep(glue("{HH}_clone_nr_{clone_nr}_"), fasta_files, value = TRUE)
    fasta <- readDNAStringSet(filepath = glue("{fasta_path}/{filename}"))
    clone <- str_extract(filename, "\\d+_\\d+(?=\\.fasta)")
    
    # Extract metadata
    seqs_meta <- HH_spec_clones_vj %>% filter(clone_subgroup_id_90_similarity == clone & locus == "IGH") %>% select(L1_annotation, c_call, sample_clean_fol)
    
    # Map seq_names on meta data
    seq_names <- names(fasta)[1:length(fasta)-1]
    length(seq_names) == nrow(seqs_meta)
    seqs_meta$seq_name <- seq_names
    
    # Read idmap.txt file
    sample <- str_split_i(filename, "\\.", 1)
    idmap <- read.csv(glue("48_GCtree/out_{version}/{sample}/idmap.txt"), header = FALSE, col.names = c("seq_unique", "seq_name")) 
    
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
        c_call = paste(unique(c_call), collapse = ":"),
        sample_clean_fol = paste(unique(sample_clean_fol), collapse = ":")
      ) %>% 
      # mutate(
      #   L1_annotation_int = sapply(L1_annotation, function(x) {
      #     parts <- str_split(x, ":")[[1]]
      #     paste(L1_int_mapping[parts], collapse = ":")
      #   })
      # )
      mutate(
        L1_annotation_int = as.integer(factor(L1_annotation)),
        c_call_int = as.integer(factor(c_call)),
        sample_clean_fol_int = as.integer(factor(sample_clean_fol))
      )
    
    # gctree_meta$L1_annotation %>% table()
    # gctree_meta$c_call %>% table()
    # gctree_meta$sample_clean_fol %>% table()
    
    # gctree_meta$L1_annotation_int %>% table()
    # gctree_meta$c_call %>% table()
    # gctree_meta$sample_clean_fol %>% table()
  
    
    write.csv(gctree_meta, glue("48_GCtree/gctree_meta_{version}/{sample}_gctree_meta.txt"), row.names = FALSE, quote = FALSE)
    
    # L1_counts.csv
    L1_counts <- gctree_meta %>%
      select(seq_unique, seq_name) %>%       
      separate_rows(seq_name, sep = ":") %>% 
      left_join(seqs_meta %>% select(seq_name, L1_annotation), by = "seq_name") %>%  
      group_by(seq_unique, L1_annotation) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = L1_annotation, values_from = count, values_fill = 0)

    write.csv(L1_counts, glue("48_GCtree/gctree_meta_{version}/{sample}_L1_counts.csv"), row.names = FALSE)
    
    # isotype_counts.csv
    isotype_counts <- gctree_meta %>%
      select(seq_unique, seq_name) %>%       
      separate_rows(seq_name, sep = ":") %>% 
      left_join(seqs_meta %>% select(seq_name, c_call), by = "seq_name") %>%  
      group_by(seq_unique, c_call) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = c_call, values_from = count, values_fill = 0)
    
    write.csv(isotype_counts, glue("48_GCtree/gctree_meta_{version}/{sample}_isotype_counts.csv"), row.names = FALSE)
    
    # sample_clean_fol_counts.csv
    sample_clean_fol_counts <- gctree_meta %>%
      select(seq_unique, seq_name) %>%       
      separate_rows(seq_name, sep = ":") %>% 
      left_join(seqs_meta %>% select(seq_name, sample_clean_fol), by = "seq_name") %>%  
      group_by(seq_unique, sample_clean_fol) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = sample_clean_fol, values_from = count, values_fill = 0)

    write.csv(sample_clean_fol_counts, glue("48_GCtree/gctree_meta_{version}/{sample}_sample_clean_fol_counts.csv"), row.names = FALSE)

  }
  
}

