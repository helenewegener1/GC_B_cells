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
seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds") # Broad annotation
seurat_obj_ADT_demulti <- readRDS("10_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds") # Follicol information 

# Investigate non-unique chains across contigs which should be filtered on umi. 

# Subset cells with tcr respectively. 
tcr_mask <- lapply(names(seurat_obj_ADT_demulti), function(x) {"tcr_v_gene_contig_1" %in% colnames(seurat_obj_ADT_demulti[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with tcr data 
tcr_seurat_obj_list <- seurat_obj_ADT_demulti[tcr_mask]

# Loading Data into scRepertoire
t_contigs.list <- list()
for (sample_name in names(tcr_seurat_obj_list)){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  # Load contig annotation file for sample 
  t_contigs <- read.csv(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_t/filtered_contig_annotations.csv"))
  
  # Load Seurat object
  seurat_obj <- seurat_obj_ADT_demulti[[sample_name]]
  
  # If sample is muliplexed, split contigs into list of contig files. 
  if ("manual_ADT_ID" %in% colnames(seurat_obj@meta.data)){
    
    t_contigs <- createHTOContigList(t_contigs, 
                                     seurat_obj, 
                                     group.by = "manual_ADT_ID")
    
    # Have name include sample_name
    names(t_contigs) <- paste(sample_name, names(t_contigs), sep = "_")
    
    # Merge list with t_contigs.list
    t_contigs.list <- c(t_contigs.list, t_contigs)
    
  } else {
    
    # Append contiguous annotation of non-multiplexed sample to t_contigs.list
    t_contigs.list[[sample_name]] <- t_contigs
    
  }
  
}

# Check names of contig list 
names(t_contigs.list)

# Here we have one row per contig
t_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 100 - N contigs
t_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$contig_id %>% str_split_i("_", 1) %>% unique() %>% length() # 49 - N cells
t_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3 %>% unique() %>% length() # 91 - amino acid sequence of CDR3
t_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt %>% unique() %>% length() # 91 - nucleotide sequence of CDR3 (More unique because silent mutations)

# ------------------------------------------------------------------------------
# Understanding start
# ------------------------------------------------------------------------------

# df <- t_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`
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


# What combineTCR does

# ------------------------------------------------------------------------------
# Understanding end
# ------------------------------------------------------------------------------

# combineTCR - one line per cell
# Combine using the default similarity clustering
combined.TCR.NOTfiltered <- combineTCR(t_contigs.list,
                                       samples = names(t_contigs.list), 
                                       # filterNonproductive = FALSE, 
                                       filterMulti = FALSE 
) 

combined.TCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 49 - N cells 
combined.TCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa1 %>% unique() %>% length() # 41 - TCR1 (alpha chain) amino acid sequence
combined.TCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt1 %>% unique() %>% length() # 41 - TCR1 (alpha chain) nucleotide sequence
combined.TCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa2 %>% unique() %>% length() # 44 - TCR2 (beta chain) amino acid sequence
combined.TCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt2 %>% unique() %>% length() # 44 - TCR2 (beta chain) sequence

# Combine using the default similarity clustering
combined.TCR.filtered <- combineTCR(t_contigs.list,
                                    samples = names(t_contigs.list), 
                                    removeNA = TRUE,
                                    # filterNonproductive = TRUE, # Default. Removes non-productive contigs , keeping only functional receptor chains. 
                                    filterMulti = TRUE # Default. For cells with more than one heavy or light chain detected, this automatically selects the chain with the highest UMI count and discards the others. 
) 

combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 43 - N cells 
combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa1 %>% unique() %>% length() # 40 - TCR1 (alpha chain) amino acid sequence
combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt1 %>% unique() %>% length() # 40 - TCR1 (alpha chain) nucleotide sequence
combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa2 %>% unique() %>% length() # 39 - TCR2 (beta chain) amino acid sequence
combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt2 %>% unique() %>% length() # 39 - TCR2 (beta chain) sequence

# How many receptors with NA in any chains 
is.na(combined.TCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$TCR1) %>% table() # 6
is.na(combined.TCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$TCR2) %>% table() # 0

is.na(combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$TCR1) %>% table() # 0 
is.na(combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$TCR2) %>% table() # 0 

head(combined.TCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`, n = 5)

# Adding Variables for Plotting: sample_high_level
combined.TCR.filtered <- addVariable(combined.TCR.filtered,
                                     variable.name = "sample_high_level",
                                     variables = str_split_i(names(combined.TCR.filtered), "_", 1))

# Adding Variables for Plotting: Inflamed / non-inflamed
patient_number <- str_split_i(names(t_contigs.list), "-", 1)
patient_number <- ifelse(patient_number == "HH117", "HH117_Crohns", "HH117_Control")

combined.TCR.filtered <- addVariable(combined.TCR.filtered,
                                     variable.name = "patient",
                                     variables = patient_number)

# The CTstrict column contains cluster IDs (e.g., "cluster.1")
head(combined.TCR.filtered[[1]][, c("barcode", "CTstrict", "TCR1", "cdr3_aa1")])

saveRDS(combined.TCR.filtered, "20_VDJ/out/combined.TCR.filtered.rds")

# ------------------------------------------------------------------------------
# Follicle frequence after filtering 
# ------------------------------------------------------------------------------

combined.TCR.names <- names(combined.TCR.filtered)

fol_freq_list <- list()

sheet_names <- list(
  "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH" = "HH117",
  "HH119-SI-PP-CD19-Pool1" = "HH119-CD19-Pool1",        
  "HH119-SI-PP-CD19-Pool2" = "HH119-CD19-Pool2",               
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = "HH119-GC-AND-PB-AND-TFH-Pool1",
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = "HH119-GC-AND-PB-AND-TFH-Pool2"
)

for (combined.TCR.name in combined.TCR.names){
  
  # combined.TCR.name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1"
  # combined.TCR.name <- "HH117-SILP-INF-PC"
  combined.TCR.filtered.sample <- combined.TCR.filtered[[combined.TCR.name]]
  
  if (str_detect(combined.TCR.name, "_", negate = TRUE)){
    next
  }
  
  sample_name <- combined.TCR.name %>% str_split_i("_", 1)
  sample_name_sheet_name <- sheet_names[[sample_name]]
  fol_name <- combined.TCR.name %>% str_split_i("_", 2)
  
  fol_freq_line <- c(fol_name, nrow(combined.TCR.filtered.sample))
  
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
out_file <- "20_VDJ/table/TCR_filtered_fol_freq.xlsx"

# Use openxlsx::write.xlsx, which takes the named list and writes
# each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
openxlsx::write.xlsx(
  x = fol_freq_df_list,
  file = out_file,
  overwrite = TRUE # Overwrite the file if it already exists
)
