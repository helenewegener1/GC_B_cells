getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
library(scRepertoire)
library(gtools)

# Load data
seurat_obj_list <- readRDS("11_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds")

# Investigate non-unique chains across contigs which should be filtered on umi. 

# Subset cells with BCR respectively. 
bcr_mask <- lapply(names(seurat_obj_list), function(x) {"bcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with BCR data 
bcr_seurat_obj_list <- seurat_obj_list[bcr_mask]

# Loading Data into scRepertoire
b_contigs.list <- list()
for (sample_name in names(bcr_seurat_obj_list)){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  # Load contig annotation file for sample 
  b_contigs <- read.csv(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/filtered_contig_annotations.csv"))
  
  # Load Seurat object
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # If sample is muliplexed, split contigs into list of contig files. 
  if ("manual_ADT_ID" %in% colnames(seurat_obj@meta.data)){
    
    b_contigs <- createHTOContigList(b_contigs, 
                                     seurat_obj, 
                                     group.by = "manual_ADT_ID")
    
    # Have name include sample_name
    names(b_contigs) <- paste(sample_name, names(b_contigs), sep = "_")
    
    # Merge list with b_contigs.list
    b_contigs.list <- c(b_contigs.list, b_contigs)
    
  } else {
    
    # Append contiguous annotation of non-multiplexed sample to b_contigs.list
    b_contigs.list[[sample_name]] <- b_contigs
    
  }
  
}

# Check names of contig list 
names(b_contigs.list)

# Here we have one row per contig
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 536 - N contigs
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$contig_id %>% str_split_i("_", 1) %>% unique() %>% length() # 295 - N cells
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3 %>% unique() %>% length() # 486 - amino acid sequence of CDR3
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt %>% unique() %>% length() # 500 - nucleotide sequence of CDR3 (More unique because silent mutations)

# ------------------------------------------------------------------------------
# Understanding start
# ------------------------------------------------------------------------------

# df <- b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`
# 
# # Count chains per barcode
# chain_counts <- df %>%
#   group_by(barcode, chain) %>%
#   summarise(n_chains = n(), .groups = "drop") %>%
#   pivot_wider(names_from = chain, 
#               values_from = n_chains, 
#               values_fill = 0)
# 
# # If you have IGH, IGK, IGL columns, create a summary
# # Adjust column names based on what you actually have
# chain_summary <- chain_counts %>%
#   mutate(
#     total_heavy = if("IGH" %in% names(.)) IGH else 0,
#     total_light = rowSums(select(., matches("IGK|IGL")), na.rm = TRUE),
#     chain_config = paste0("IGH:", total_heavy, " Light:", total_light)
#   ) %>% 
#   count(chain_config) %>%
#   arrange(desc(n))
# 
# ggplot(chain_summary, aes(x = reorder(chain_config, n), y = n)) +
#   geom_col(fill = "steelblue") +
#   geom_text(aes(label = n), hjust = -0.2) +
#   coord_flip() +
#   labs(
#     title = "Distribution of Chain Configurations",
#     x = "Chain Configuration",
#     y = "Number of Cells"
#   ) +
#   theme_minimal(base_size = 10) +
#   theme(plot.title = element_text(face = "bold"))


# What combineBCR does
## Keeps all 295 cells (even those with only heavy or only light) - from plot --> 208 + 62 + 10= 280
## For the 5 cells with 2+ heavy chains: selects one (highest UMI if filterMulti=TRUE)
## For the 13 cells with 2+ light chains: selects one (highest UMI if filterMulti=TRUE)
## Assigns heavy to cdr3_aa1/nt1 and light to cdr3_aa2/nt2
## Cells missing a chain will have NA in the corresponding columns

# ------------------------------------------------------------------------------
# Understanding end
# ------------------------------------------------------------------------------

# combineBCR - one line per cell
# Combine using the default similarity clustering
combined.BCR.NOTfiltered <- combineBCR(b_contigs.list,
                                       samples = names(b_contigs.list), 
                                       # filterNonproductive = FALSE, 
                                       filterMulti = FALSE 
) 

combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 295 - N cells 
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa1 %>% unique() %>% length() # 218 - heavy chain amino acid sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt1 %>% unique() %>% length() # 219 - heavy chain nucleotide sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa2 %>% unique() %>% length() # 256 - ligth chain amino acid sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt2 %>% unique() %>% length() # 267 - ligth chain nucleotide sequence

# Combine using the default similarity clustering
combined.BCR.filtered <- combineBCR(b_contigs.list,
                                    samples = names(b_contigs.list), 
                                    removeNA = TRUE,
                                    # threshold = 0.85, # Default is 0.85. Oliver used default. 
                                    # filterNonproductive = TRUE, # Default. Removes non-productive contigs , keeping only functional receptor chains. 
                                    filterMulti = TRUE # Default. For cells with more than one heavy or light chain detected, this automatically selects the chain with the highest UMI count and discards the others. 
) 

combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 223 - N cells 
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa1 %>% unique() %>% length() # 209 - heavy chain amino acid sequence
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt1 %>% unique() %>% length() # 210 - heavy chain nucleotide sequence
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa2 %>% unique() %>% length() # 200 - ligth chain amino acid sequence
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt2 %>% unique() %>% length() # 208 - ligth chain nucleotide sequence

# How many receptors with NA in any chains 
is.na(combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGH) %>% table() # 62
is.na(combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGLC) %>% table() # 10 

is.na(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGH) %>% table() # 0 
is.na(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGLC) %>% table() # 0 

head(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`, n = 5)

# Adding Variables: sample_high_level
combined.BCR.filtered <- addVariable(combined.BCR.filtered,
                                     variable.name = "sample_high_level",
                                     variables = str_split_i(names(combined.BCR.filtered), "_", 1))

# Adding Variables: Inflamed / non-inflamed
patient_number <- str_split_i(names(b_contigs.list), "-", 1)
patient_number <- ifelse(patient_number == "HH117", "HH117_Crohns", "HH117_Control")

combined.BCR.filtered <- addVariable(combined.BCR.filtered,
                                     variable.name = "patient",
                                     variables = patient_number)

# Adding Variables: Broad cell type


# The CTstrict column contains cluster IDs (e.g., "cluster.1")
head(combined.BCR.filtered[[1]][, c("barcode", "CTstrict", "IGH", "cdr3_aa1")])

saveRDS(combined.BCR.filtered, "20_VDJ/out/combined.BCR.filtered.rds")
# combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.rds")

# ------------------------------------------------------------------------------
# Follicle frequence after filtering 
# ------------------------------------------------------------------------------

combined.BCR.names <- names(combined.BCR.filtered)

fol_freq_list <- list()

sheet_names <- list(
  "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = "HH117",
  "HH119-SI-PP-CD19-Pool1" = "HH119-CD19-Pool1",        
  "HH119-SI-PP-CD19-Pool2" = "HH119-CD19-Pool2",               
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = "HH119-GC-AND-PB-AND-TFH-Pool1",
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = "HH119-GC-AND-PB-AND-TFH-Pool2"
)

for (combined.BCR.name in combined.BCR.names){
  
  # combined.BCR.name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1"
  # combined.BCR.name <- "HH117-SILP-INF-PC"
  combined.BCR.filtered.sample <- combined.BCR.filtered[[combined.BCR.name]]
  
  if (str_detect(combined.BCR.name, "_", negate = TRUE)){
    next
  }
  
  sample_name <- combined.BCR.name %>% str_split_i("_", 1)
  sample_name_sheet_name <- sheet_names[[sample_name]]
  fol_name <- combined.BCR.name %>% str_split_i("_", 2)
  
  fol_freq_line <- c(fol_name, nrow(combined.BCR.filtered.sample))
  
  fol_freq_list[[sample_name_sheet_name]] <- rbind(fol_freq_list[[sample_name_sheet_name]], fol_freq_line)
  
}

# Get in order 

fol_freq_df_list <- lapply(fol_freq_list, function(x) {
  as.data.frame(x) %>% 
    mutate(
      fol_num = as.numeric(str_extract(V1, "(?<=Fol-)\\d+")),
      order_val = case_when(
        V1 == "Doublet" ~ Inf,      # Doublet goes second to last
        V1 == "Negative" ~ Inf + 1, # Negative goes last
        TRUE ~ fol_num                # Fol-X ordered by number
      )
    ) %>%
    arrange(order_val) %>%
    select(-fol_num, -order_val) %>% 
    rename(Fol = V1, Count = V2) %>% 
    remove_rownames()
})

# Export xlsx file 
out_file <- "20_VDJ/table/BCR_filtered_fol_freq.xlsx"

# Use openxlsx::write.xlsx, which takes the named list and writes
# each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
openxlsx::write.xlsx(
  x = fol_freq_df_list,
  file = out_file,
  overwrite = TRUE # Overwrite the file if it already exists
)
