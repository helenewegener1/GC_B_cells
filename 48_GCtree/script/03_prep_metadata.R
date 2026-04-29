library(Biostrings) # writing fasta files 
library(msa) # mulitple sequence alignment 
library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

resolve_LC_list_germlined <- readRDS("45_immcantation/out/rds/resolve_LC_list_germlined.rds")

patients <- names(resolve_LC_list_germlined)

fasta_path <- "48_GCtree/fasta"
fasta_files <- list.files(fasta_path)

# ------------------------------------------------------------------------------
# Meta data file writing
# ------------------------------------------------------------------------------



clone_nrs <- 1:5

for (HH in patients){
  
  # HH <- "HH117"
  HH_spec_clones_vj <- resolve_LC_list_germlined[[HH]]
  
  for (clone_nr in clone_nrs){
    
    # Get top "clone_nr" clone from given sample
    # clone_nr <- 5
    filename <- grep(glue("{HH}_clone_nr_{clone_nr}"), fasta_files, value = TRUE)
    fasta <- readDNAStringSet(filepath = glue("{fasta_path}/{filename}"))
    clone <- str_extract(filename, "\\d+_\\d+(?=\\.fasta)")
    
    # Extract metadata
    seqs_meta <- HH_spec_clones_vj %>% filter(clone_subgroup_id == clone & locus == "IGH") %>% select(L1_annotation, c_call)
    
    # Map seq_names on meta data
    seq_names <- names(fasta)[1:length(fasta)-1]
    length(seq_names) == nrow(seqs_meta)
    seqs_meta$seq_name <- seq_names
    
    # Read idmap.txt file
    sample <- str_split_i(filename, "\\.", 1)
    idmap <- read.csv(glue("48_GCtree/out/{sample}/idmap.txt"), header = FALSE, col.names = c("seq_unique", "seq_name")) 
    
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
        c_call = paste(unique(c_call), collapse = ":")
      )
    
    # gctree_meta$L1_annotation %>% table()
    # gctree_meta$c_call %>% table()
    
    write.csv(gctree_meta, glue("48_GCtree/gctree_meta/{sample}_gctree_meta.txt"), row.names = FALSE, quote = FALSE)
    
    
  }
  
}

