getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
library(gtools)
library(purrr)
library(tibble)
# BiocManager::install("scRepertoire", force = TRUE) # 2.6.2
# devtools::install_github("BorchLab/scRepertoire", force = TRUE) # 2.6.2
# devtools::install_github("BorchLab/scRepertoire@v2.6.0", force = TRUE)
# remotes::install_github(c("BorchLab/immApex", "BorchLab/scRepertoire"), force = TRUE)

# devtools::install_github("BorchLab/scRepertoire@v2.6.0", force = TRUE)
# devtools::install_github("BorchLab/immApex@v1.4.3", force = TRUE)

# devtools::install_github("BorchLab/scRepertoire@v2.5.0", force = TRUE)

library(scRepertoire)
packageVersion("scRepertoire")
library(immApex)
packageVersion("immApex")

# Load data
# seurat_obj_list <- readRDS("11_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed_all.rds")
# seurat_obj_list <- readRDS("13_add_metadata/out/seurat_obj_prepped_list.rds")
seurat_obj_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

# ------------------------------------------------------------------------------
# LOAD BCR DATA - running on computerome.
# ------------------------------------------------------------------------------
# 
# # Subset cells with BCR respectively. 
# bcr_mask <- lapply(names(seurat_obj_list), function(x) {"bcr_v_gene_contig_1" %in% colnames(seurat_obj_list[[x]]@meta.data)}) %>% unlist()
# 
# # Extract seurat objects with BCR data 
# bcr_seurat_obj_list <- seurat_obj_list[bcr_mask]
# 
# # Loading Data into scRepertoire
# b_contigs.list <- list() # stratified by follicles
# full_sequence_IGH.list <- list() 
# for (sample_name in names(bcr_seurat_obj_list)){
#   
#   # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
#   # sample_name <- "HH119-SILP-PC"      
#   
#   # Load contig annotation file for sample
#   b_contigs <- read.csv(glue("05_run_cellranger/out_v9/res_{sample_name}/outs/per_sample_outs/res_{sample_name}/vdj_b/filtered_contig_annotations.csv"))
#   
#   # Load Seurat object
#   seurat_obj <- seurat_obj_list[[sample_name]]
#   
#   # If sample is muliplexed, split contigs into list of contig files and added to b_contigs.list
#   if ("manual_ADT_ID" %in% colnames(seurat_obj@meta.data)){
#     
#     b_contigs <- createHTOContigList(b_contigs, 
#                                      seurat_obj, 
#                                      group.by = "manual_ADT_ID")
#     
#     # Have name include sample_name
#     names(b_contigs) <- paste(sample_name, names(b_contigs), sep = "_")
#     
#     # Merge list with b_contigs.list
#     b_contigs.list <- c(b_contigs.list, b_contigs)
#     
#     # Get full sequence of IgH
#     full_sequence_IGH.list <- c(full_sequence_IGH.list, lapply(b_contigs, function(x){
#       x %>%
#         filter(
#           chain == "IGH",
#           productive == "true",
#           high_confidence == "true",
#           is_cell == "true"
#         ) %>% 
#         group_by(barcode) %>% 
#         slice_max(umis, n = 1, with_ties = FALSE) %>% 
#         ungroup() %>% 
#         mutate(
#           IGH_full_sequence = paste0(
#             fwr1_nt,
#             cdr1_nt,
#             fwr2_nt,
#             cdr2_nt,
#             fwr3_nt,
#             cdr3_nt,
#             fwr4_nt
#           ),
#           sample = sample %>% str_split_i("_", 2), 
#           barcode = paste(sample, barcode, sep = "_")
#         ) %>% 
#         select(barcode, IGH_full_sequence, umis)
#     }))
#     
#   } else {
#     
#     # Append contiguous annotation of non-multiplexed sample to b_contigs.list
#     b_contigs.list[[sample_name]] <- b_contigs
#   
#     # Get full sequence of IgH
#     full_sequence_IGH.list[[sample_name]] <- b_contigs %>%
#       filter(
#         chain == "IGH",
#         productive == "true",
#         high_confidence == "true",
#         is_cell == "true"
#       ) %>% 
#       group_by(barcode) %>% 
#       slice_max(umis, n = 1, with_ties = FALSE) %>% 
#       ungroup() %>% 
#       mutate(
#         IGH_full_sequence = paste0(
#           fwr1_nt,
#           cdr1_nt,
#           fwr2_nt,
#           cdr2_nt,
#           fwr3_nt,
#           cdr3_nt,
#           fwr4_nt
#         ),
#         sample = sample %>% str_split_i("_", 2), 
#         barcode = paste(sample, barcode, sep = "_")
#       ) %>% 
#       select(barcode, IGH_full_sequence, umis)
# 
#   }
#   
# }

b_contigs.list <- readRDS("20_VDJ/out/b_contigs.list.rds")
full_sequence_IGH.list <- readRDS("20_VDJ/out/full_sequence_IGH.list.rds")

# Check names of contig list 
names(b_contigs.list)
names(full_sequence_IGH.list)

# Here we have one row per contig
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% nrow() # 536 - N contigs
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$contig_id %>% str_split_i("_", 1) %>% unique() %>% length() # 295 - N cells
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3 %>% unique() %>% length() # 486 - amino acid sequence of CDR3
b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$cdr3_nt %>% unique() %>% length() # 500 - nucleotide sequence of CDR3 (More unique because silent mutations)

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
# Filter cells that have heavy chain and compare to no filteringgit s
# ------------------------------------------------------------------------------

combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13` %>% filter(!is.na(IGH)) %>% nrow()

n_cells_prefilt <- lapply(names(combined.BCR.NOTfiltered), function(x) combined.BCR.NOTfiltered[[x]] %>% nrow())
n_cells_w_IGH <- lapply(names(combined.BCR.NOTfiltered), function(x) combined.BCR.NOTfiltered[[x]] %>% filter(!is.na(IGH)) %>% nrow())
 
df <- cbind(n_cells_prefilt, n_cells_w_IGH, "sample" = names(combined.BCR.NOTfiltered)) %>% 
  as_tibble() %>% unnest(c(n_cells_prefilt, n_cells_w_IGH, sample)) %>% 
  mutate(
    sample = sample %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC"),
    percentage_loss = (n_cells_prefilt-n_cells_w_IGH)/n_cells_prefilt * 100 
  ) %>% 
  select(sample, n_cells_prefilt, n_cells_w_IGH, percentage_loss)

xlsx::write.xlsx(df, "20_VDJ/table/BCR_both_vs_heavy_chain.xlsx")

filter_conditions <- list(
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1" = function(df) filter(df, str_detect(sample, "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool1")),
  "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2" = function(df) filter(df, str_detect(sample, "HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2")),
  "HH119-SI-PP-CD19-Pool1" = function(df) filter(df, str_detect(sample, "HH119-SI-PP-CD19-Pool1")),
  "HH119-SI-PP-CD19-Pool2" = function(df) filter(df, str_detect(sample, "HH119-SI-PP-CD19-Pool2")),
  "HH119-SI-PP" = function(df) filter(df, str_detect(sample, "HH119-SI-PP")),
  "other"        = function(df) filter(df, !str_detect(sample, paste(c("HH117-SI-PP", "HH119-SI-PP"), collapse = "|")))
)

for (condition_name in names(filter_conditions)) {
  
  # condition_name <- "other"
  
  df_plot <- df %>%
    filter_conditions[[condition_name]](.) %>%
    filter(!str_detect(sample, "Negative|Doublet")) %>%
    pivot_longer(cols = starts_with("n_cells"),
                 names_to = "filter",
                 values_to = "n_cells")
  
  n_cell_max <- max(df_plot$n_cells)
  
  df_plot %>% 
    ggplot(aes(x = filter, y = n_cells, color = sample, group = sample)) + 
    geom_point() + 
    geom_line() + 
    scale_y_continuous(breaks = seq(0, n_cell_max + 25, by = 100)) +
    labs(
      x = "",
      y = "N cells"
    ) + 
    theme_bw()
  
  ggsave(glue("20_VDJ/plot/N_cell_stat/BCR_both_vs_heavy_chain/BCR_both_vs_heavy_chain_sample_{condition_name}.png"), height = 10, width = 8)
  
}


# ------------------------------------------------------------------------------
# combineBCR - Filtering 
# ------------------------------------------------------------------------------

combined.BCR.filtered_HH117 <- combineBCR(
  b_contigs.list[startsWith(names(b_contigs.list), "HH117")],
  samples = names(b_contigs.list[startsWith(names(b_contigs.list), "HH117")]),
  removeNA = TRUE,
  threshold = 0.85,
  filterNonproductive = TRUE,
  filterMulti = TRUE
)

combined.BCR.filtered_HH119_noFol <- combineBCR(
  b_contigs.list[startsWith(names(b_contigs.list), "HH119") & !grepl("Pool", names(b_contigs.list))],
  samples = names(b_contigs.list[startsWith(names(b_contigs.list), "HH119") & !grepl("Pool", names(b_contigs.list))]),
  removeNA = TRUE,
  threshold = 0.85,
  filterNonproductive = TRUE,
  filterMulti = TRUE
)

combined.BCR.filtered_HH119_pool1 <- combineBCR(
  b_contigs.list[startsWith(names(b_contigs.list), "HH119") & grepl("Pool1", names(b_contigs.list))],
  samples = names(b_contigs.list[startsWith(names(b_contigs.list), "HH119") & grepl("Pool1", names(b_contigs.list))]),
  removeNA = TRUE,
  threshold = 0.85,
  filterNonproductive = TRUE,
  filterMulti = TRUE
)

combined.BCR.filtered_HH119_pool2 <- combineBCR(
  b_contigs.list[startsWith(names(b_contigs.list), "HH119") & grepl("Pool2", names(b_contigs.list))],
  samples = names(b_contigs.list[startsWith(names(b_contigs.list), "HH119") & grepl("Pool2", names(b_contigs.list))]),
  removeNA = TRUE,
  threshold = 0.85,
  filterNonproductive = TRUE,
  filterMulti = TRUE
)

# Combine cells across patients 
combined.BCR.filtered_list <- list(
  HH117 = combined.BCR.filtered_HH117 %>% bind_rows(),
  HH119 = list(
    combined.BCR.filtered_HH119_noFol,
    combined.BCR.filtered_HH119_pool1,
    combined.BCR.filtered_HH119_pool2
  ) %>% flatten() %>% bind_rows()
)

rm(combined.BCR.filtered_HH117, combined.BCR.filtered_HH119_noFol, combined.BCR.filtered_HH119_pool1, combined.BCR.filtered_HH119_pool2)

# Define clone types on heavy chain 
clonalCluster_IGH_list <- lapply(combined.BCR.filtered_list, function(x) {
  clonalCluster(
    list(x),            
    chain = "IGH",
    sequence = "nt",
    threshold = 0.85,
    dist_type = "levenshtein",
    normalize = "length",
    cluster.method = "components",
    cluster.prefix = "cluster.",
    use.V = TRUE,
    use.J = TRUE,
    exportAdjMatrix = FALSE,
    exportGraph = FALSE
  )[[1]]                 
})

# Define clone types on light chain 
clonalCluster_Light_list <- lapply(combined.BCR.filtered_list, function(x) {
  clonalCluster(
    list(x),             
    chain = "Light",
    sequence = "nt",
    threshold = 0.85,
    dist_type = "levenshtein",
    normalize = "length",
    cluster.method = "components",
    cluster.prefix = "cluster.",
    use.V = TRUE,
    use.J = TRUE,
    exportAdjMatrix = FALSE,
    exportGraph = FALSE
  )[[1]]                
})

rm(combined.BCR.filtered_list)

# Combine clonal types for heavy and light chain 
combined.BCR.filtered_clonalCluster <- map2(
  clonalCluster_IGH_list,
  clonalCluster_Light_list,
  function(igh, light) {
    inner_join(igh, light) %>%
      mutate(
        IGH_V_gene     = str_split_i(IGH,  "\\.", 1),
        IGH_J_gene     = str_split_i(IGH,  "\\.", 3),
        IGLC_V_gene    = str_split_i(IGLC, "\\.", 1),
        IGLC_J_gene    = str_split_i(IGLC, "\\.", 2),
        gene_cluster   = paste0(IGH_V_gene, "_", IGH_J_gene, "_", IGLC_V_gene, "_", IGLC_J_gene),
        manual_cluster = paste0(IGH.Cluster, "_", Light.Cluster)
      )
  }
)

# Clean up 
rm(clonalCluster_Light_list, clonalCluster_IGH_list)

# Check
combined.BCR.filtered_clonalCluster$HH119 %>% 
  summarize(n = n(), .by = c("gene_cluster", "manual_cluster")) %>%
  filter(gene_cluster == "IGHV4-34_IGHJ4_IGLV1-40_IGLJ7") 
  arrange(desc(n)) %>% head(n = 20)

# combined.BCR.filtered_clonalCluster[["HH117-SILP-INF-PC"]] %>% filter(is.na(IGH.Cluster)) %>% head()

saveRDS(combined.BCR.filtered_clonalCluster, "20_VDJ/out/combined.BCR.filtered_clonalCluster.rds")
# combined.BCR.filtered_clonalCluster <- readRDS("20_VDJ/out/combined.BCR.filtered_clonalCluster.rds")

combined.BCR.filtered <- combined.BCR.filtered_clonalCluster

################################################################################
# One row per cell
combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13") %>% nrow() # 223 - N cells 
combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13") %>% select(cdr3_aa1) %>% unique() %>% nrow() # 209 - heavy chain amino acid sequence
combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13") %>% select(cdr3_nt1) %>% unique() %>% nrow() # 210 - heavy chain nucleotide sequence
combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13") %>% select(cdr3_aa2) %>% unique() %>% nrow() # 200 - ligth chain amino acid sequence
combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13") %>% select(cdr3_nt2) %>% unique() %>% nrow() # 208 - ligth chain nucleotide sequence

# How many receptors with NA in any chains 
# is.na(combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGH) %>% table() # 62
# is.na(combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGLC) %>% table() # 10 
# 
# combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13", !is.na(IGH))
# is.na(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGH) %>% table() # 0 
# is.na(combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-13`$IGLC) %>% table() # 0 

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

# Bind rows
full_sequence_IGH <- list(
  HH117 = dplyr::bind_rows(Filter(function(df) any(grepl("HH117", df$barcode)), full_sequence_IGH.list)),
  HH119 = dplyr::bind_rows(Filter(function(df) any(grepl("HH119", df$barcode)), full_sequence_IGH.list))
)

combined.BCR.filtered.clean <- list()

for (HH in names(combined.BCR.filtered)) {
  
  # HH <- "HH119"
  
  combined.BCR.filtered_HH <- combined.BCR.filtered[[HH]]
  
  # Clean sample names
  seurat_obj <- subset(seurat_obj_integrated, patient == HH)
  seurat_obj$sample_high_level <- seurat_obj$sample %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  seurat_obj[[]] <- seurat_obj[[]] %>% mutate(sample = ifelse(!is.na(manual_ADT_ID), glue("{sample_high_level}_{manual_ADT_ID}"), sample_high_level))
  
  # Clean sample names 
  combined.BCR.filtered_HH$sample_original <- combined.BCR.filtered_HH$sample
  combined.BCR.filtered_HH$sample <- combined.BCR.filtered_HH$sample %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  combined.BCR.filtered_HH$barcode_original <- combined.BCR.filtered_HH$barcode
  combined.BCR.filtered_HH$barcode_original_noFol <- gsub("_Fol-\\d+", "", combined.BCR.filtered_HH$barcode_original)
  combined.BCR.filtered_HH$barcode <- combined.BCR.filtered_HH$barcode %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  combined.BCR.filtered_HH$sample_high_level <- str_split_i(combined.BCR.filtered_HH$sample, "_", 1)
   
  # Get barcodes from Seurat (cells that passed QC)
  # Update barcodes
  barcode_seurat <- rownames(seurat_obj[[]]) %>% str_extract("(?<=_)[^_]+_[^_]+")
  seurat_obj$barcode <- paste(seurat_obj$sample, barcode_seurat, sep = "_")
  
  # Get barcode suffix
  sample_suffix_map <- seurat_obj[[]] %>%
    select(sample_high_level) %>%
    rownames_to_column("barcode") %>%
    mutate(barcode_suffix = str_split_i(barcode, "_", 3)) %>%
    select(sample_high_level, barcode_suffix) %>%
    distinct()
  
  # Then fix BCR barcodes to match integrated object
  combined.BCR.filtered_HH <- combined.BCR.filtered_HH %>%
    left_join(sample_suffix_map, by = "sample_high_level", relationship = "many-to-many") %>% 
    mutate(
      barcode_integrated = paste0(
        barcode,
        "_", barcode_suffix # Add integration suffix
      )
    ) %>% 
    rename(barcode_original_2 = barcode, 
           barcode = barcode_integrated)
  
  # Look at barcodes
  combined.BCR.filtered_HH$barcode %>% tail()
  seurat_obj$barcode %>% tail()
  table(combined.BCR.filtered_HH$barcode %in% seurat_obj$barcode)
  
  # Filter barcodes in combined.BCR.filtered based on seurat obejct since these cells have been QC-filtered
  combined.BCR.filtered_HH_QC <- combined.BCR.filtered_HH %>%
    filter(barcode %in% seurat_obj$barcode)
  
  # Check numbers 
  table(combined.BCR.filtered_HH$barcode %in% seurat_obj$barcode)
  nrow(combined.BCR.filtered_HH)
  nrow(combined.BCR.filtered_HH_QC)
  
  # -------------------------
  # Add cell type annotation 
  # -------------------------
  
  seurat_meta <- seurat_obj[[]] %>% select(barcode, celltype_broad)
  seurat_meta$barcode %>% tail()
  combined.BCR.filtered_HH_QC$barcode %>% tail()
  
  combined.BCR.filtered_HH_QC <- left_join(combined.BCR.filtered_HH_QC, seurat_meta, by = "barcode")
  
  # -------------------------
  # Add patient 
  # -------------------------
  
  patient <- ifelse(HH == "HH117", "HH117_Crohns", "HH119_Control")

  combined.BCR.filtered_HH_QC$patient <- patient
  
  # -------------------------
  # Remove Pool information in sample and barcode
  # -------------------------
  combined.BCR.filtered_HH_QC$barcode <- combined.BCR.filtered_HH_QC$barcode %>% str_remove_all("-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-CD19-Pool1|-CD19-Pool2")
  combined.BCR.filtered_HH_QC$sample <- combined.BCR.filtered_HH_QC$sample %>% str_remove_all("-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-CD19-Pool1|-CD19-Pool2")
  combined.BCR.filtered_HH_QC$sample_high_level <- combined.BCR.filtered_HH_QC$sample %>% str_remove_all("-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-CD19-Pool1|-CD19-Pool2")
  
  # unique(combined.BCR.filtered_HH$sample)
  # unique(combined.BCR.filtered_HH$sample_high_level)
  
  # -------------------------
  # Remove negatives and doublets
  # -------------------------
  
  combined.BCR.filtered_HH_QC <- combined.BCR.filtered_HH_QC %>% filter(!str_detect(sample, "Negative|Doublet"))
  
  # -------------------------
  # Add full sequence 
  # -------------------------
  
  full_sequence_IGH_HH <- full_sequence_IGH[[HH]]
  
  full_sequence_IGH_HH$barcode_original_noFol <- full_sequence_IGH_HH$barcode
  
  # Check
  table(full_sequence_IGH_HH$barcode_original_noFol %in% combined.BCR.filtered_HH_QC$barcode_original_noFol)
  combined.BCR.filtered_HH_QC$barcode_original_noFol %>% length()
  combined.BCR.filtered_HH_QC$barcode_original_noFol %>% unique() %>% length() # same as ^^
  
  combined.BCR.filtered_HH_QC <- combined.BCR.filtered_HH_QC %>% left_join(full_sequence_IGH_HH, by = "barcode_original_noFol")
  
  # -------------------------
  # Add genes and Ig class
  # -------------------------
  
  combined.BCR.filtered_HH_QC <- combined.BCR.filtered_HH_QC %>% 
    mutate(
      IGH_V_gene = IGH %>% str_split_i("\\.", 1),
      IGH_J_gene = IGH %>% str_split_i("\\.", 3),
      IGLC_V_gene = IGLC %>% str_split_i("\\.", 1),
      IGLC_J_gene = IGLC %>% str_split_i("\\.", 2),
      gene_cluster   = paste0(IGH_V_gene, "_", IGH_J_gene, "_", IGLC_V_gene, "_", IGLC_J_gene),
      Ig_class = CTgene %>% str_split_i("_", 1) %>% str_split_i("\\.", 4),
      Ig_class = na_if(Ig_class, "")
  )

  # -------------------------
  # Add sample_clean
  # -------------------------
  
  combined.BCR.filtered_HH_QC <- combined.BCR.filtered_HH_QC %>% 
    mutate(sample_clean = sample %>% str_split_i("_", 1))
  
  # -------------------------
  # Export dataframe
  # -------------------------

  combined.BCR.filtered.clean[[HH]] <- combined.BCR.filtered_HH_QC
  
}

saveRDS(combined.BCR.filtered.clean, "20_VDJ/out/combined.BCR.filtered.clean.rds")
# combined.BCR.joined <- readRDS("20_VDJ/out/combined.BCR.joined.rds")

# ------------------------------------------------------------------------------
# Inspect top clones 
# ------------------------------------------------------------------------------

combined.BCR.filtered$HH119 %>%
  summarise(n = n(), .by = c("manual_cluster", "gene_cluster")) %>% 
  filter(!str_detect(manual_cluster, "NA")) %>% 
  arrange(desc(n)) %>% 
  head(20)

# Are there clones with the same genes?
combined.BCR.filtered$HH119  %>%
  # filter(IGH_V_gene == "IGHV4-34" & IGH_J_gene == "IGHJ4") %>% 
  filter(gene_cluster == "IGHV4-34_IGHJ4_IGLV1-40_IGLJ7") %>% 
  summarise(n = n(), .by = c("manual_cluster", "gene_cluster")) %>% 
  arrange(desc(n)) %>% 
  head(15)

combined.BCR.filtered$HH117 %>%
  summarise(n = n(), .by = c("manual_cluster", "gene_cluster")) %>% 
  filter(!str_detect(manual_cluster, "NA")) %>% 
  arrange(desc(n)) %>% 
  head(20)

# ------------------------------------------------------------------------------
# N cell tracking 
# ------------------------------------------------------------------------------

# Base data for all original samples
# df_N_cells <- data.frame(
#   sample_fol = names(combined.BCR.NOTfiltered),
#   sample_name = names(combined.BCR.NOTfiltered) %>% str_split_i("_", 1),
#   combineBCR_NOTfiltered = sapply(combined.BCR.NOTfiltered, nrow),
#   combineBCR_filtered = sapply(combined.BCR.filtered, nrow),
#   seurat_intersection_filtered = sapply(combined.BCR.filtered.clean, nrow),
#   row.names = NULL
# )

# Add pool_combined stage
# pool_combined_counts <- data.frame(
#   sample_fol = names(combined.BCR.filtered.clean),
#   pool_combine = sapply(combined.BCR.filtered.clean.pool_combine, nrow),
#   row.names = NULL
# )

# Removing negatives/doublets
# final_counts <- data.frame(
#   sample_fol = names(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs),
#   final_clean = sapply(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs, nrow),
#   row.names = NULL
# )

final_counts <- combined.BCR.filtered.clean %>%
  bind_rows() %>% 
  summarize(n = n(), .by = "sample")

# # Check that merging of pools went well
# pool_summary <- df_N_cells %>% 
#   filter(str_detect(sample_name, "Pool")) %>% 
#   mutate(Fol = str_split_i(sample_fol, "_", 2)) %>% 
#   filter(!(Fol %in% c("Negative", "Doublet"))) %>% 
#   group_by(Fol) %>% 
#   summarise(n_cells_total = sum(seurat_intersection_filtered)) %>% 
#   arrange(Fol)
# 
# pool_summary_final <- final_counts %>% filter(str_detect(sample_fol, "HH119-SI-PP")) %>% arrange(sample_fol)
# 
# cbind(pool_summary, pool_summary_final) %>% mutate(adds_up = ifelse(n_cells_total == final_clean, TRUE, FALSE))
# 
# # Plot
# for (SAMPLE_NAME in unique(df_N_cells$sample_name)){
# 
#   # SAMPLE_NAME <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
# 
#   df_N_cells %>%
#     filter(sample_name == SAMPLE_NAME & str_detect(sample_fol, "Negative|Doublet", negate = TRUE)) %>%
#     pivot_longer(
#       cols = -c(sample_fol, sample_name),
#       names_to = "stage",
#       values_to = "n_cells"
#     ) %>%
#     mutate(stage = factor(stage, levels = c("combineBCR_NOTfiltered", "combineBCR_filtered", "seurat_intersection_filtered"))) %>%
#     ggplot(aes(x = stage, y = n_cells, group = sample_fol, color = sample_fol)) +
#     geom_line() +
#     geom_point(size = 2) +
#     theme_bw() +
#     theme(legend.position = "right") +
#     labs(
#       title = "Cell Filtering Through BCR Processing Pipeline",
#       subtitle = SAMPLE_NAME, 
#       x = "Processing Stage",
#       y = "Number of Cells",
#       color = "Sample"
#     )
# 
#   ggsave(glue("20_VDJ/plot/N_cell_stat/BCR_cell_filtering/BCR_cell_filtering_{SAMPLE_NAME}.png"), width = 12, height = 8)
# 
# }

# Final N cells barplot
final_counts %>% ggplot(aes(y = sample, x = n)) + 
  geom_col() + 
  geom_text(aes(label = n), hjust = -0.1, size = 3) +
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

# table(seurat_obj_list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH`$manual_ADT_ID == "Fol-1")
# b_contigs.list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
# combined.BCR.NOTfiltered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
# combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
# combined.BCR.filtered.clean$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
# combined.BCR.filtered.clean.pool_combine$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
# combined.BCR.filtered.clean.pool_combine.rm_neg_doubs$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow()
# 
# # HH117-SILP-INF-PC
# ncol(seurat_obj_list$`HH117-SILP-INF-PC`)
# b_contigs.list$`HH117-SILP-INF-PC` %>% nrow()
# combined.BCR.NOTfiltered$`HH117-SILP-INF-PC` %>% nrow()
# combined.BCR.filtered$`HH117-SILP-INF-PC` %>% nrow()
# combined.BCR.filtered.clean$`HH117-SILP-INF-PC` %>% nrow()
# combined.BCR.filtered.clean.pool_combine$`HH117-SILP-INF-PC` %>% nrow()
# combined.BCR.filtered.clean.pool_combine.rm_neg_doubs$`HH117-SILP-INF-PC` %>% nrow()
# 
# 
# # ------------------------------------------------------------------------------
# # Follicle frequence after filtering 
# # ------------------------------------------------------------------------------
# 
# combined.BCR.names <- names(combined.BCR.filtered.clean.pool_combine.rm_neg_doubs)
# 
# fol_freq_list <- list()
# 
# for (combined.BCR.name in combined.BCR.names){
#   
#   # combined.BCR.name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1"
#   # combined.BCR.name <- "HH117-SILP-INF-PC"
#   combined.BCR.filtered.clean.sample <- combined.BCR.filtered.clean.pool_combine.rm_neg_doubs[[combined.BCR.name]]
#   
#   if (str_detect(combined.BCR.name, "_", negate = TRUE)){
#     next
#   }
#   
#   sample_name <- combined.BCR.name %>% str_split_i("_", 1)
#   sample_name_sheet_name <- sample_name %>% str_split_i("-", 1)
#   fol_name <- combined.BCR.name %>% str_split_i("_", 2)
#   
#   fol_freq_line <- c(fol_name, nrow(combined.BCR.filtered.clean.sample))
#   
#   fol_freq_list[[sample_name_sheet_name]] <- rbind(fol_freq_list[[sample_name_sheet_name]], fol_freq_line)
#   
# }
# 
# # Get in order 
# 
# fol_freq_df_list <- lapply(fol_freq_list, function(x) {
#   as.data.frame(x) %>% 
#     mutate(
#       fol_num = as.numeric(str_extract(V1, "(?<=Fol-)\\d+")),
#       order_val = case_when(
#         V1 == "Doublet" ~ Inf,      # Doublet goes second to last
#         V1 == "Negative" ~ Inf + 1, # Negative goes last
#         TRUE ~ fol_num                # Fol-X ordered by number
#       )
#     ) %>%
#     arrange(order_val) %>%
#     select(-fol_num, -order_val) %>% 
#     rename(Fol = V1, Count = V2) %>% 
#     remove_rownames()
# })
# 
# # Export xlsx file 
# out_file <- "20_VDJ/table/BCR_filtered_fol_freq.xlsx"
# 
# # Use openxlsx::write.xlsx, which takes the named list and writes
# # each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
# openxlsx::write.xlsx(
#   x = fol_freq_df_list,
#   file = out_file,
#   overwrite = TRUE # Overwrite the file if it already exists
# )

