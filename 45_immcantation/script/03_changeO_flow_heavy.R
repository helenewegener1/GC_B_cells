# Load libraries
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
# suppressPackageStartupMessages(library(dowser)) # Needs to be installed 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(scoper))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(shazam))
library(tibble)
library(patchwork)
library(forcats)
library(glue)
library(stringr)

packageVersion("airr")
packageVersion("alakazam")
packageVersion("scoper")
packageVersion("shazam")

# Following this Immcantation flow:
# https://immcantation.readthedocs.io/en/latest/getting_started/10x_tutorial.html

# ------------------------------------------------------------------------------
# Get sample names
# ------------------------------------------------------------------------------

base_path <- "45_immcantation/out"
files <- list.files(base_path, full.names = FALSE)
sample_names <- files[1:length(files)-1]

# Read all samples, adding sample_id and subject_id
bcr_data_tmp <- lapply(sample_names, function(s) {
  # s <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  f <- file.path(base_path, s, paste0(s, "_heavy_germ-pass.tsv"))
  db <- airr::read_rearrangement(f)
  db$sample_id <- s
  db$subject_id <- sub("-(SI|SILP|CO|COLP).*", "", s)  # extracts HH117 or HH119
  # make sequence and cell IDs unique across samples
  db$sequence_id <- paste0(s, "_", db$sequence_id)
  db$cell_id <- paste0(s, "_", db$cell_id)
  return(db)
}) %>% bind_rows()

# split by patient
bcr_data <- list(
  "HH117" = bcr_data_tmp %>% filter(subject_id == "HH117"),
  "HH119" = bcr_data_tmp %>% filter(subject_id == "HH119")
)

cat(paste("HH117:", nrow(bcr_data$HH117), "sequences\n"))
cat(paste("HH119:", nrow(bcr_data$HH119), "sequences\n"))



# ------------------------------------------------------------------------------
# Check V/D/J gene call consistency
# ------------------------------------------------------------------------------
# 
# # Filter out inconsistent sequences:
# # Filter out contigs that differ in heavy/light/kappa/lambda chain across V, J and C gene. 
# 
# bcr_data_filt <- lapply(bcr_data, function(x) {
#   
#   df <- x %>% 
#     filter(grepl("^IGHV", v_call) & grepl("^IGHJ", j_call))
#   
#     # Filter out contigs that differ in heavy/light/kappa/lambda chain across V, J and C gene. 
#     # dplyr::filter(
#     #   (grepl("^IGHV", v_call) & grepl("^IGHJ", j_call) & grepl("^IGH[MGADE]", c_call)) |
#     #     (grepl("^IGKV", v_call) & grepl("^IGKJ", j_call) & grepl("^IGKC", c_call)) |
#     #     (grepl("^IGLV", v_call) & grepl("^IGLJ", j_call) & grepl("^IGLC", c_call))
#     # )
#   
#   return(df)
#   
# }
# )
# 
# # Rows in the data after filtering V/J/C calls inconsistent with the respective locus
# cat(paste("HH117:", nrow(bcr_data_filt$HH117), "sequences\n"))
# cat(paste("HH119:", nrow(bcr_data_filt$HH119), "sequences\n"))
# 
# bcr_data$HH117[!(bcr_data$HH117$sequence_id %in% bcr_data_filt$HH117$sequence_id), ] %>% 
#   select(v_call, j_call, c_call)

# ------------------------------------------------------------------------------
# Remove non-productive sequences
# ------------------------------------------------------------------------------

bcr_data$HH117$productive %>% table()
bcr_data$HH119$productive %>% table()

# bcr_data <- lapply(bcr_data, function(x) {
# 
#   df <- x %>% filter(productive)
# 
#   cat(paste("There are", nrow(x), "rows in the data.\n"))
# 
#   return(df)
# 
# }
# )
# 
# cat(paste("HH117:", nrow(bcr_data_filt$HH117), "sequences\n"))
# cat(paste("HH119:", nrow(bcr_data_filt$HH119), "sequences\n"))

# ------------------------------------------------------------------------------
# Handle multiple heavy chains 
# ------------------------------------------------------------------------------

# Filter cells based on multiple heavy chain
bcr_data_qc <- lapply(bcr_data, function(x) {
  
  df <- x %>%
    group_by(cell_id) %>%
    arrange(desc(umi_count), .by_group = TRUE) %>%
    mutate(
      n_heavy = n(),
      dominant = case_when(
        n_heavy == 1 ~ TRUE,                          # only one heavy chain, keep it
        umi_count[1] >= 2 * umi_count[2] ~ row_number() == 1,  # dominant contig has 2x UMIs
        TRUE ~ FALSE                                  # ambiguous, drop all contigs for this cell
      )
    ) %>%
    filter(dominant) %>%
    select(-n_heavy, -dominant) %>%
    ungroup()
  
  return(df)
  
}
)

# Rows in the data before filtering 
cat(paste("HH117:", nrow(bcr_data$HH117), "sequences\n"))
cat(paste("HH119:", nrow(bcr_data$HH119), "sequences\n"))

# Rows in the data after filtering 
cat(paste("HH117:", nrow(bcr_data_qc$HH117), "sequences\n"))
cat(paste("HH119:", nrow(bcr_data_qc$HH119), "sequences\n"))

# ------------------------------------------------------------------------------
# Add cell type annotation 
# ------------------------------------------------------------------------------

seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")
patients <- names(bcr_data_qc)

# Filter cells based on multiple heavy chain
bcr_data_qc_annot <- lapply(patients, function(HH) {
  
  # Define data
  # HH <- "HH117"
  
  # -------------------
  # Standardize cell IDs
  # -------------------
  
  df <- bcr_data_qc[[HH]]
  seurat_obj <- subset(seurat_integrated, patient == HH)
  
  # Check cell_ids/barcodes
  rownames(seurat_obj) %>% head()
  df$cell_id %>% head()
  
  # Add sample name to seurat_obj
  seurat_obj$cell_id_seurat <- paste0(seurat_obj$sample, "_", colnames(seurat_obj))
  seurat_obj$sample_id <- seurat_obj$sample
  seurat_obj$cell_id_seurat %>% head()
  
  # Add more meta data 
  df$sample_clean <- df$sample_id %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")
  df$cell_id_noFol <- gsub("_Fol-\\d+", "", df$cell_id)
  
  # Get barcode suffix
  sample_suffix_map <- seurat_obj[[]] %>%
    mutate(barcode_suffix = str_split_i(cell_id_seurat, "_", 3)) %>%
    select(sample_id, barcode_suffix) %>%
    distinct()
  
  # Then fix BCR barcodes to match integrated object
  df <- df %>%
    left_join(sample_suffix_map, by = "sample_id", relationship = "many-to-many") %>% 
    mutate(
      cell_id_seurat = paste0(
        cell_id_noFol,
        "_", barcode_suffix # Add integration suffix
      )
    )
  
  # Check 
  df$cell_id_seurat %>% sort() %>% tail()
  seurat_obj$cell_id_seurat %>% sort() %>% tail()
  
  table(df$cell_id_seurat %in% seurat_obj$cell_id_seurat)
  # df$cell_id_seurat[!(df$cell_id_seurat %in% seurat_obj$cell_id_seurat)]
  
  # -------------------
  # Transfer cell type annotations into the BCR data
  # -------------------
  
  seurat_meta <- seurat_obj[[]] %>% 
    select(cell_id_seurat, celltype_broad, L1_annotation, manual_ADT_class, manual_ADT_ID, manual_ADT_full_ID)
  
  df <- df %>% left_join(seurat_meta, by = "cell_id_seurat")
  
  # -------------------
  # Remove cells without GEX data
  # -------------------
  
  df <- df %>% filter(!is.na(celltype_broad))
  
  # -------------------
  # Remove negative follicles and doublets 
  # -------------------
  
  df <- df %>% filter(!(manual_ADT_class %in% c("Negative", "Doublet")))
  # df <- df %>% mutate(ifelse(manual_ADT_class %in% c("Negative", "Doublet"), NA, manual_ADT_class)) # We want this 
  
  # -------------------
  # Combine pools of HH119-SI-PP
  # -------------------
  
  if (HH == "HH119"){
    
    df$sample_clean <- df$sample_clean %>% str_remove_all("-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2")
    
  }
  
  # -------------------
  # Add patient information 
  # -------------------
  
  df$patient_id <- HH
  
  # -------------------
  # Remove TFH cells and GC B cells from the LP and contamination 
  # -------------------
  
  LP_samples <- grep("LP", df$sample_clean, value = TRUE) %>% unique()
  
  # table(df$L1_annotation == "Tfh_cells")
  # table(df$L1_annotation == "GC_B_cells" & df$sample_clean %in% LP_samples)
  
  # df <- df %>% filter(
  #   (celltype_broad != "Tfh_like_cells") & !(df$celltype_broad == "GC_B_cells" & df$sample_clean %in% LP_samples)
  # )
  # nrow(df)
  
  df <- df %>% 
    mutate(L1_annotation = ifelse(L1_annotation == "GC_Bcells", "GC_B_cells", L1_annotation)) %>% 
    filter(
      (L1_annotation != "Tfh_cells"), 
      !(L1_annotation == "GC_B_cells" & sample_clean %in% LP_samples), 
      (str_detect(L1_annotation, "Contamination", negate = TRUE))
    )

  # nrow(df)
  
  # -------------------
  # Return final df
  # -------------------
  
  return(df)
  
}) %>% 
  setNames(patients)


# Rows in the data after subsetting with seurat object  
cat(paste("HH117:", nrow(bcr_data_qc_annot$HH117), "sequences\n"))
cat(paste("HH119:", nrow(bcr_data_qc_annot$HH119), "sequences\n"))

# ------------------------------------------------------------------------------
# Summary of filtering
# ------------------------------------------------------------------------------

# Rows in the data before filtering 
cat(paste("HH117:", nrow(bcr_data$HH117), "sequences\n"))
cat(paste("HH119:", nrow(bcr_data$HH119), "sequences\n"))

# Rows in the data after filtering 
cat(paste("HH117:", nrow(bcr_data_qc$HH117), "sequences\n"))
cat(paste("HH119:", nrow(bcr_data_qc$HH119), "sequences\n"))

# Rows in the data after subsetting with seurat object  
cat(paste("HH117:", nrow(bcr_data_qc_annot$HH117), "sequences\n"))
cat(paste("HH119:", nrow(bcr_data_qc_annot$HH119), "sequences\n"))

# ------------------------------------------------------------------------------
# Export bcr_data_qc_annot
# ------------------------------------------------------------------------------

saveRDS(bcr_data_qc_annot, "45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")
# bcr_data_qc_annot <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")

# ------------------------------------------------------------------------------
# Summary follicles and cell types
# ------------------------------------------------------------------------------

source("10_broad_annotation/script/color_palette.R")
patients <- names(bcr_data_qc_annot)

lapply(patients, function(HH){
  
  # HH <- "HH119"
  
  bcr_data_qc_annot[[HH]] %>% 
    ggplot(aes(x = sample_clean, fill = L1_annotation)) +
    geom_bar() + 
    scale_fill_manual(values = L1_colors) + 
    theme_bw() + 
    labs(
      x = "", 
      y = "Count", 
      title = glue ("{HH}: N cells across samples")
    ) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(glue("45_immcantation/plot/{HH}_N_cells_across_samples.png"), width = 12, height = 7)
  
  bcr_data_qc_annot[[HH]] %>% 
    filter(!is.na(manual_ADT_ID)) %>% 
    ggplot(aes(x = manual_ADT_ID, fill = L1_annotation)) + 
    geom_bar() + 
    scale_fill_manual(values = L1_colors) + 
    theme_bw() + 
    labs(
      x = "", 
      y = "Count", 
      title = glue("{HH}: N cells across follicles of SI-PP")
    ) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  ggsave(glue("45_immcantation/plot/{HH}_N_cells_across_follicles.png"), width = 14, height = 7)
  
})

