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
# seurat_obj_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_singlets_list.rds")

# ------------------------------------------------------------------------------
# LOAD BCR DATA
# ------------------------------------------------------------------------------

# Subset cells with BCR respectively. 
bcr_mask <- lapply(names(seurat_obj_list), function(x) {"bcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()

# Extract seurat objects with BCR data 
bcr_seurat_obj_list <- seurat_obj_list[bcr_mask]

# Loading Data into scRepertoire
b_contigs.list <- list() # stratified by follicles
for (sample_name in names(bcr_seurat_obj_list)){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  
  # Load contig annotation file for sample 
  b_contigs <- read.csv(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/filtered_contig_annotations.csv"))
  
  # Load Seurat object
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  # If sample is muliplexed, split contigs into list of contig files and added to b_contigs.list
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
# Understanding combineBCR
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
# combineBCR - no filtering 
# ------------------------------------------------------------------------------

# combineBCR - one line per cell
# Combine using the default similarity clustering
combined.BCR.NOTfiltered <- combineBCR(b_contigs.list,
                                       samples = names(b_contigs.list), 
                                       filterNonproductive = FALSE,
                                       filterMulti = FALSE 
) 

# One row per cell
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 295 - N cells 
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa1 %>% unique() %>% length() # 218 - heavy chain amino acid sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt1 %>% unique() %>% length() # 219 - heavy chain nucleotide sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_aa2 %>% unique() %>% length() # 256 - ligth chain amino acid sequence
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt2 %>% unique() %>% length() # 267 - ligth chain nucleotide sequence

# ------------------------------------------------------------------------------
# combineBCR - Filtering 
# ------------------------------------------------------------------------------

# Combine using the default similarity clustering
combined.BCR.filtered <- combineBCR(b_contigs.list,
                                    samples = names(b_contigs.list), 
                                    removeNA = TRUE,
                                    threshold = 0.85, # Default is 0.85. Oliver used default.
                                    filterNonproductive = TRUE, # Default. Removes non-productive contigs , keeping only functional receptor chains.
                                    filterMulti = TRUE # Default. For cells with more than one heavy or light chain detected, this automatically selects the chain with the highest UMI count and discards the others. 
) 

# One row per cell
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

# head(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`, n = 5)

# # Also for the no-follicle-stratification object (13 items)
# combined.BCR.filtered.13 <- combineBCR(b_contigs.list.13,
#                                        samples = names(b_contigs.list.13), 
#                                        removeNA = TRUE,
#                                        threshold = 0.85, # Default is 0.85. Oliver used default.
#                                        # filterNonproductive = TRUE, # Default. Removes non-productive contigs , keeping only functional receptor chains. 
#                                        filterMulti = TRUE # Default. For cells with more than one heavy or light chain detected, this automatically selects the chain with the highest UMI count and discards the others. 
# ) 



# ------------------------------------------------------------------------------
# Filter combined.BCR.filtered based on seurat object 
# ------------------------------------------------------------------------------

combined.BCR.filtered.clean <- list()

for (sample_name_fol in names(combined.BCR.filtered)) {
  
  # sample_name_fol <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13"
  # sample_name_fol <- "HH117-SI-MILF-INF-HLADR-AND-CD19"
  # sample_name_fol <- "HH117-SILP-INF-PC"
  
  sample_name <- sample_name_fol %>% str_split_i("_", 1)
  Fol_name <- sample_name_fol %>% str_split_i("_", 2)
  
  # Get barcodes from Seurat (cells that passed QC)
  seurat_obj <- seurat_obj_list[[sample_name]]
  seurat_obj$barcode <- paste(sample_name_fol, rownames(seurat_obj[[]]), sep = "_")
  
  # Filter for manual_ADT_ID
  if (str_detect(sample_name_fol, "_Fol-")){
    seurat_obj <- subset(seurat_obj, manual_ADT_ID == Fol_name)
  }

  # Get BCR data
  bcr_data <- combined.BCR.filtered[[sample_name_fol]]
  
  # Filter barcodes in combined.BCR.filtered based on seurat obejct since these cells have been QC-filtered
  bcr_filtered <- bcr_data %>%
    filter(barcode %in% seurat_obj$barcode)
  
  # Check numbers 
  # table(bcr_data$barcode %in% seurat_obj$barcode)
  nrow(bcr_data)
  nrow(bcr_filtered)
  
  # -------------------------
  # Add cell type annotation 
  # -------------------------
  
  seurat_meta <- seurat_obj[[]] %>% select(barcode, celltype_broad)
  
  bcr_filtered_annotated <- left_join(bcr_filtered, seurat_meta, by = "barcode")
  
  combined.BCR.filtered.clean[[sample_name_fol]] <- bcr_filtered_annotated
  
  # -------------------------
  # Add sample_high_level 
  # -------------------------
  
  combined.BCR.filtered.clean[[sample_name_fol]]$sample_high_level <- sample_name
  
  # -------------------------
  # Add patient 
  # -------------------------
  
  patient_number <- str_split_i(sample_name, "-", 1)
  patient <- ifelse(patient_number == "HH117", "HH117_Crohns", "HH117_Control")
  
  combined.BCR.filtered.clean[[sample_name_fol]]$sample_high_level <- patient
  
}

# ------------------------------------------------------------------------------
# Combine Pools 
# ------------------------------------------------------------------------------

combined.BCR.filtered.clean.pool_combine <- combined.BCR.filtered.clean

pools_fols <- list(
  "Pool1" = grep("Pool1", names(combined.BCR.filtered.clean.pool_combine), value = TRUE) %>% str_split_i("_", 2) %>% unique(), 
  "Pool2" = grep("Pool2", names(combined.BCR.filtered.clean.pool_combine), value = TRUE) %>% str_split_i("_", 2) %>% unique()
)

for (pool in names(pools_fols)){
  
  # pool <- "Pool1"
  fols <- pools_fols[[pool]]
  
  # Define samples
  for (fol in fols){
    
    # fol <- "Fol-7"
    
    sample_name_new <- glue("HH119-SI-PP_{fol}")
    
    s1 <- glue("HH119-SI-PP-CD19-{pool}_{fol}")
    s2 <- glue("HH119-SI-PP-GC-AND-PB-AND-TFH-{pool}_{fol}")
    
    # Combine samples 
    combined.BCR.filtered.clean.pool_combine[[sample_name_new]] <- rbind(combined.BCR.filtered.clean.pool_combine[[s1]], combined.BCR.filtered.clean.pool_combine[[s2]])
    
    combined.BCR.filtered.clean.pool_combine[[sample_name_new]]$barcode_old <- combined.BCR.filtered.clean.pool_combine[[sample_name_new]]$barcode
    combined.BCR.filtered.clean.pool_combine[[sample_name_new]]$barcode <- paste(
      sample_name_new, 
      combined.BCR.filtered.clean.pool_combine[[sample_name_new]]$barcode_old %>% str_split_i("_", 3), sep = "_"
      )
    
    combined.BCR.filtered.clean.pool_combine[[sample_name_new]]$sample_old <- combined.BCR.filtered.clean.pool_combine[[sample_name_new]]$sample
    combined.BCR.filtered.clean.pool_combine[[sample_name_new]]$sample <- sample_name_new
    
    # Remove old samples 
    combined.BCR.filtered.clean.pool_combine <- combined.BCR.filtered.clean.pool_combine[
      !names(combined.BCR.filtered.clean.pool_combine) %in% c(s1, s2)
    ]
  }
  
}

length(combined.BCR.filtered.clean.pool_combine)

# ------------------------------------------------------------------------------
# Remove negatives and doublets 
# ------------------------------------------------------------------------------

combined.BCR.filtered.clean.pool_combine.rm_neg_doubs <- combined.BCR.filtered.clean.pool_combine

not_keep <- grep("Negative|Doublet", names(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs), value = TRUE)

combined.BCR.filtered.clean.pool_combine.rm_neg_doubs <- combined.BCR.filtered.clean.pool_combine.rm_neg_doubs[
  !names(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs) %in% not_keep
]

length(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs)

# ------------------------------------------------------------------------------
# SAVE BCR DATA (combined.BCR.filtered.clean)
# ------------------------------------------------------------------------------

saveRDS(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs, "20_VDJ/out/combined.BCR.filtered.clean.rds")
# combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.rds")

# ------------------------------------------------------------------------------
# N cell tracking 
# ------------------------------------------------------------------------------

# Base data for all original samples
df_N_cells <- data.frame(
  sample_fol = names(combined.BCR.NOTfiltered),
  sample_name = names(combined.BCR.NOTfiltered) %>% str_split_i("_", 1),
  combineBCR_NOTfiltered = sapply(combined.BCR.NOTfiltered, nrow),
  combineBCR_filtered = sapply(combined.BCR.filtered, nrow),
  seurat_intersection_filtered = sapply(combined.BCR.filtered.clean, nrow),
  row.names = NULL
)

# Add pool_combined stage
pool_combined_counts <- data.frame(
  sample_fol = names(combined.BCR.filtered.clean.pool_combine),
  pool_combine = sapply(combined.BCR.filtered.clean.pool_combine, nrow),
  row.names = NULL
)

# Removing negatives/doublets
final_counts <- data.frame(
  sample_fol = names(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs),
  final_clean = sapply(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs, nrow),
  row.names = NULL
)

# Check that merging of pools went well
pool_summary <- df_N_cells %>% 
  filter(str_detect(sample_name, "Pool")) %>% 
  mutate(Fol = str_split_i(sample_fol, "_", 2)) %>% 
  filter(!(Fol %in% c("Negative", "Doublet"))) %>% 
  group_by(Fol) %>% 
  summarise(n_cells_total = sum(seurat_intersection_filtered)) %>% 
  arrange(Fol)

pool_summary_final <- final_counts %>% filter(str_detect(sample_fol, "HH119-SI-PP")) %>% arrange(sample_fol)

cbind(pool_summary, pool_summary_final) %>% mutate(adds_up = ifelse(n_cells_total == final_clean, TRUE, FALSE))

# Plot
for (SAMPLE_NAME in unique(df_N_cells$sample_name)){

  # SAMPLE_NAME <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"

  df_N_cells %>%
    filter(sample_name == SAMPLE_NAME & str_detect(sample_fol, "Negative|Doublet", negate = TRUE)) %>%
    pivot_longer(
      cols = -c(sample_fol, sample_name),
      names_to = "stage",
      values_to = "n_cells"
    ) %>%
    mutate(stage = factor(stage, levels = c("combineBCR_NOTfiltered", "combineBCR_filtered", "seurat_intersection_filtered"))) %>%
    ggplot(aes(x = stage, y = n_cells, group = sample_fol, color = sample_fol)) +
    geom_line() +
    geom_point(size = 2) +
    theme_bw() +
    theme(legend.position = "right") +
    labs(
      title = "Cell Filtering Through BCR Processing Pipeline",
      subtitle = SAMPLE_NAME, 
      x = "Processing Stage",
      y = "Number of Cells",
      color = "Sample"
    )

  ggsave(glue("20_VDJ/plot/N_cell_stat/BCR_cell_filtering/BCR_cell_filtering_{SAMPLE_NAME}.png"), width = 12, height = 8)

}

# Final N cells barplot
final_counts %>% ggplot(aes(y = sample_fol, x = final_clean)) + 
  geom_col() + 
  geom_text(aes(label = final_clean), hjust = -0.1, size = 3) +
  theme_bw() + 
  labs(
    title = "Final N cells after all filtering",
    y = "", 
    x = "N cells"
  )

ggsave("20_VDJ/plot/N_cell_stat/BCR_final_N_cells.png", width = 14, height = 15)  

# ------------------------------------------------------------------------------
# HH117 - Fol-1 stats
# ------------------------------------------------------------------------------

table(seurat_obj_list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH`$manual_ADT_ID == "Fol-1")
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
combined.BCR.filtered.clean$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
combined.BCR.filtered.clean.pool_combine$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
combined.BCR.filtered.clean.pool_combine.rm_neg_doubs$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()

# HH117-SILP-INF-PC
ncol(seurat_obj_list$`HH117-SILP-INF-PC`)
b_contigs.list$`HH117-SILP-INF-PC` %>% nrow()
combined.BCR.NOTfiltered$`HH117-SILP-INF-PC` %>% nrow()
combined.BCR.filtered$`HH117-SILP-INF-PC` %>% nrow()
combined.BCR.filtered.clean$`HH117-SILP-INF-PC` %>% nrow()
combined.BCR.filtered.clean.pool_combine$`HH117-SILP-INF-PC` %>% nrow()
combined.BCR.filtered.clean.pool_combine.rm_neg_doubs$`HH117-SILP-INF-PC` %>% nrow()


# ------------------------------------------------------------------------------
# Follicle frequence after filtering 
# ------------------------------------------------------------------------------

combined.BCR.names <- names(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs)

fol_freq_list <- list()

for (combined.BCR.name in combined.BCR.names){
  
  # combined.BCR.name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1"
  # combined.BCR.name <- "HH117-SILP-INF-PC"
  combined.BCR.filtered.clean.sample <- combined.BCR.filtered.clean.pool_combine.rm_neg_doubs[[combined.BCR.name]]
  
  if (str_detect(combined.BCR.name, "_", negate = TRUE)){
    next
  }
  
  sample_name <- combined.BCR.name %>% str_split_i("_", 1)
  sample_name_sheet_name <- sample_name %>% str_split_i("-", 1)
  fol_name <- combined.BCR.name %>% str_split_i("_", 2)
  
  fol_freq_line <- c(fol_name, nrow(combined.BCR.filtered.clean.sample))
  
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
