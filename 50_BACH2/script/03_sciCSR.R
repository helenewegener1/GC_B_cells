library(reticulate)
use_condaenv("scicsr_new", required = TRUE)
py_config()

library(sciCSR, quietly = TRUE)
source("50_BACH2/script/sciCSR_functions_fixed.R")
library(Seurat)
library(glue)
source("10_broad_annotation/script/color_palette.R")
library(patchwork)
library(tidyverse)

# Run this if this is the first time you set up sciCSR on 
# your machine; this set up the dependencies in python
# prepare_sciCSR() 

# Following tutorial: https://josephng-bio.org/sciCSR/articles/csr.html

# ------------------------------------------------------------------------------
# Prep paths 
# ------------------------------------------------------------------------------

cellranger_out <- "05_run_cellranger/out_v9/"
samples_all <- list.files(cellranger_out) %>% str_remove("res_") 
samples <- grep("SI-PP", samples_all, value = TRUE)
samples <- grep("GC", samples, value = TRUE)

outdir <- "50_BACH2/plot/03_sciCSR"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# # ------------------------------------------------------------------------------
# # Load Seurat object
# # ------------------------------------------------------------------------------
# 
# seurat_obj <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs_annotated.rds")
# 
# # Seurat samples have ocm information in name
# seurat_obj_samples <- seurat_obj$sample %>% unique()
# 
# # Subset seurat object to SI-PP-GC samples
# seurat_obj_samples <- grep("SI-PP", seurat_obj_samples, value = TRUE)
# seurat_obj_samples <- grep("GC", seurat_obj_samples, value = TRUE)
# seurat_obj <- subset(seurat_obj, sample %in% seurat_obj_samples)
# 
# # Remove Tfh cells
# seurat_obj <- subset(seurat_obj, L1_annotation != "Tfh_cells")
# 
# # Reintegrate
# res <- 0.3
# seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
# seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
# seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
# seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
# seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE, cluster.name = "subset_unintegrated_clusters")
# seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:10, verbose = FALSE, reduction.name = "subset_umap.unintegrated")
# 
# DefaultAssay(seurat_obj)
# Reductions(seurat_obj)
# seurat_obj_integrated <- seurat_obj
# 
# DefaultAssay(seurat_obj_integrated) <- "RNA"
# Layers(seurat_obj_integrated[["RNA"]])
# Idents(seurat_obj_integrated) %>% unique()
# 
# Reductions(seurat_obj_integrated)
# 
# seurat_obj_integrated <- IntegrateLayers(
#   object = seurat_obj_integrated, 
#   method = HarmonyIntegration,
#   orig.reduction = "pca", 
#   new.reduction = "subset_RNA_harmony",
#   # group.by.vars = "patient",
#   assay = "RNA", 
#   verbose = FALSE
# )
# 
# seurat_obj_integrated <- FindNeighbors(seurat_obj_integrated, reduction = "subset_RNA_harmony", dims = 1:30)
# seurat_obj_integrated <- FindClusters(seurat_obj_integrated, resolution = 0.3)
# seurat_obj_integrated <- RunUMAP(seurat_obj_integrated, reduction = "pca", dims = 1:10, verbose = FALSE, reduction.name = "subset_umap.harmony")
# 
# Reductions(seurat_obj_integrated)
# 
# DimPlot(seurat_obj_integrated, reduction = "subset_umap.harmony")
# DimPlot(seurat_obj_integrated, reduction = "subset_umap.harmony", group.by = "sample")
# 
# DimPlot(seurat_obj_integrated, group.by = "patient", reduction = "subset_umap.harmony", cols = patient_color_values)
# ggsave(glue("{outdir}/DimPlot_patient.png"))
# 
# DimPlot(seurat_obj_integrated, group.by = "L1_annotation", reduction = "subset_umap.harmony", cols = L1_colors)
# ggsave(glue("{outdir}/DimPlot_L1_annotation.png"))
# 
# # Large CRC cluster included for now 
# 
# # ------------------------------------------------------------------------------
# # Load VDJ data
# # ------------------------------------------------------------------------------
# 
# collapsed_list <- list()
# for (sample in samples){
#   
#   # sample <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
#   # sample <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB" 
#   
#   dirs <- glue("{cellranger_out}/res_{sample}/outs/per_sample_outs")
#   
#   dirs_in_outsdir <- list.dirs(dirs, recursive = FALSE, full.names = FALSE)
#   
#   if (length(dirs_in_outsdir) == 1 & dirs_in_outsdir[[1]] == glue("res_{sample}")){
#     # ocm = FALSE
#     vdj_paths <- glue("{dirs}/res_{sample}/vdj_b/filtered_contig_annotations.csv")
#     vquest_paths <- glue("45_immcantation/out/{sample}/{sample}_db-pass.tsv")
#   } else if (length(dirs_in_outsdir) > 1) {
#     # ocm = TRUE
#     vdj_paths <- glue("{dirs}/{dirs_in_outsdir}/vdj_b/filtered_contig_annotations.csv")
#     vquest_paths <- glue("45_immcantation/out/{sample}_{dirs_in_outsdir}/{sample}_{dirs_in_outsdir}_db-pass.tsv")
#   }
#   
#   for (i in 1:length(dirs_in_outsdir)){
#     
#     # this is the filtered_contig_annotations.csv file from cellranger vdj
#     vdj <- read.csv(vdj_paths[[i]], stringsAsFactors = FALSE)
#     vdj$contig_id <- glue("{sample}_{vdj$contig_id}")
#     
#     # this is the AIRR-format output from IMGT/HighV-Quest
#     vquest <- read.table(vquest_paths[[i]], sep = "\t", stringsAsFactors = FALSE, header = TRUE)
#     vquest$sequence_id <- glue("{sample}_{vquest$sequence_id}")
#     
#     # The % sequence identity to germline V allele is the column 'v_identity'
#     # We want to merge this column into cellranger vdj output 
#     # using the sequence IDs as keys
#     vdj <- merge(vdj, vquest[, c("sequence_id", "v_identity")],
#                  by.x = "contig_id", by.y = "sequence_id", 
#                  all.x = TRUE, all.y = FALSE, sort = FALSE)
#     
#     # collapseBCR annotate for each cell as singlet/doublet etc. based on the number of heavy & light transcripts associated with each cell barcode.
#     collapsed <- collapseBCR(vdj, format = "10X")
#     
#     # Add sample 
#     collapsed$sample <- sample
#     
#     # Add OCM is avaiable 
#     if (length(dirs_in_outsdir) > 1){
#       
#       ocm <- dirs_in_outsdir[[i]]
#       collapsed$ocm <- ocm
#       sample_out <- glue("{sample}_{ocm}")
#       
#     } else{
#       sample_out <- sample
#     }
#     
#     collapsed$sample_name <- sample_out 
#     
#     # Add to list 
#     collapsed_list[[sample_out]] <- collapsed
#     
#   }
#   
# }
# 
# names(collapsed_list)
# 
# collapsed_all <- bind_rows(collapsed_list)
# 
# table(collapsed_all$ocm, useNA = "always")
# 
# # We can add annotations from the Seurat Object into the VDJ table
# # e.g. we add 'Status' (WT/KD) from Seruat object into the repertoire 
# # annotations
# # collapsed <- AddCellMetaToVDJ(
# #   vdj = collapsed,
# #   SeuratObj = seurat_obj_integrated,
# #   metadata_col = c("Status")
# # )
# 
# # We can add data the other way round, i.e. select features from
# # repertoire and add them to the Seurat Object
# 
# # First make sure we have a column in the vdj table which corresponds to
# # sample names - this column must be called 'sample_name'
# # In this case it is the column 'donor' - first rename this to 'sample_name'
# # colnames(collapsed)[which(colnames(collapsed) == "donor")] <- "sample_name"
# 
# # Make sure barcode column and link collapsed_all and seurat_obj_integrated
# collapsed_all$barcode <- glue("{collapsed_all$sample_name}_{collapsed_all$barcode}")
# seurat_obj_integrated$barcode <- glue("{seurat_obj_integrated$sample}_{colnames(seurat_obj_integrated)}") %>% str_remove("_\\d+")
# 
# # Check overlap 
# table(collapsed_all$barcode %in% seurat_obj_integrated$barcode)
# collapsed_all$barcode[collapsed_all$barcode %in% seurat_obj_integrated$barcode] %>% str_remove("_[ATCG]{16}-\\d+$") %>% table()
# 
# # Define 
# keep_columns = c("v_gene", "d_gene", "j_gene", "c_gene", "v_identity", "full_length", 
#                  "productive", "cdr3", "cdr3_nt", "reads", "umis", "barcode", "ocm")
# 
# keep_columns %in% colnames(collapsed_all) %>% table()
# collapsed_all %>% select(keep_columns) %>% head()
# 
# seurat_obj_integrated_comb <- combineBCR(
#   vdj = collapsed_all, SeuratObj = seurat_obj_integrated,
#   # list the columns from vdj you wish to add to the Seurat object down here
#   keep_columns = keep_columns,
#   seurat_cell_name_column = "barcode",
#   seurat_sample_column = "sample"
# )
# 
# # plot DimPlot but colour by isotype from the vdj data
# # to confirm the data merging has happened
# DimPlot(seurat_obj_integrated_comb, group.by = "IGH_c_gene", reduction = "subset_umap.harmony", cols = isotype_colors_custom)
# ggsave(glue("{outdir}/DimPlot_isotype.png"))
# 
# # ------------------------------------------------------------------------------
# # Begin sciCSR workflow 
# # ------------------------------------------------------------------------------
# 
# # pattern for files 
# # "05_run_cellranger/out_v9/res_HH117-SILP-INF-PC/outs/per_sample_outs/res_HH117-SILP-INF-PC/count/sample_alignments.bam"
# 
# # We first need to load the genomic coordinates of the V, D and J genes
# # change to 'human_definitions' if you are working on human data
# # # data("mouse_definitions")
# # data("human_definitions")
# 
# # we need to supply the list of BAM files associated with the scRNA-seq
# # in this case we have five libraries, a subset of the BAM file corresponding
# # to the cells are available in the package
# 
# # # Check for OCMs
# # bam_paths <- list()
# # for (sample in samples){
# # 
# #   # sample <- "HH151-SI-PP-nonINF-MEM-AND-GC-AND-TFH-AND-PB"
# #   # sample <- "HH117-SILP-INF-PC"
# # 
# #   dirs <- glue("{cellranger_out}/res_{sample}/outs/per_sample_outs")
# # 
# #   dirs_in_outsdir <- list.dirs(dirs, recursive = FALSE, full.names = FALSE)
# # 
# #   if (length(dirs_in_outsdir) == 1 & dirs_in_outsdir[[1]] == glue("res_{sample}")){
# #     # ocm = FALSE
# #     bam_path <- glue("{dirs}/res_{sample}/count/sample_alignments.bam")
# #   } else if (length(dirs_in_outsdir) > 1) {
# #     # ocm = TRUE
# #     bam_path <- glue("{dirs}/{dirs_in_outsdir}/count/sample_alignments.bam")
# #   }
# # 
# #   bam_paths <- c(bam_paths, bam_path)
# # 
# # }
# # 
# # bamfiles <- unlist(bam_paths)
# # 
# # # Check that all files exists
# # stopifnot(all(file.exists(bamfiles)))
# # 
# # # loop through these BAM files
# # subset_IGH <- lapply(bamfiles, function(bamfile){
# #   cat(paste0(bamfile, " ...\n"))
# #   # These three functions count productive and sterile
# #   # reads from the BAM file and give a count matrix
# #   out <- getIGHmapping(bamfile, human_definitions)
# #   # out2 <- getIGHreadType(out$read_count)
# #   out2 <- getIGHreadType_fixed(out$read_count)
# #   out3 <- summariseIGHreads(out2, human_definitions)
# #   # 'out3' is a count matrix of productive/sterile transcript 
# #   # of each isotype per cell
# #   out3  
# # })
# # 
# # 
# # # Add names to subset_IGH
# # subset_IGH_samples <- bamfiles %>% str_split_i("/", 4) %>% str_remove("res_")
# # subset_IGH_ocms <- bamfiles %>% str_split_i("/", 7)
# # subset_IGH_names <- glue("{subset_IGH_samples}_{subset_IGH_ocms}") %>% str_remove("_?res_.*$")
# # names(subset_IGH) <- subset_IGH_names
# # 
# # saveRDS(subset_IGH, "50_BACH2/out/subset_IGH.rds")
# subset_IGH <- readRDS("50_BACH2/out/subset_IGH.rds")
# 
# # subset_IGH$
# 
# # Match subset_IGH and seurat_obj_integrated_comb barcodes 
# x <- seurat_obj_integrated_comb$sample
# y <- colnames(seurat_obj_integrated_comb) %>% str_split_i("_", 2)
# 
# table(x, y)
# 
# # Sample names are the same 
# table(unique(x) %in% names(subset_IGH))
# 
# # Repair barcodes 
# subset_IGH2 <- repairBarcode(
#   subset_IGH, 
#   seurat_obj_integrated_comb,
#   sample_names = names(subset_IGH), 
#   seurat_sample_column = "sample"
# )
# 
# # Check that repair of barcode was correct
# table_x_y <- table(x, y) %>% as.data.frame() %>% filter(Freq != 0)
# 
# for (sample in names(subset_IGH)){
# 
#   print(sample)
#   subset_IGH2_y <- subset_IGH2[[sample]] %>% rownames() %>% str_split_i("_", 2) %>% unique()
#   y_check <- table_x_y %>% filter(x == sample) %>% pull(y)
# 
#   print(subset_IGH2_y == y_check)
# 
# }
# 
# # combine these individual count matrices
# subset_IGH2 <- do.call("rbind", subset_IGH2)
# 
# # remove cells which are not in the Seurat object but happen to have
# # observed productive/sterile transcripts
# subset_IGH2 <- subset_IGH2[which(rownames(subset_IGH2) %in% Cells(seurat_obj_integrated_comb)), ]
# 
# # Now that we have this productive/sterile transcript count matrix, we would want to merge this count matrix as a new ‘assay’ in the Seurat object. 
# # Occassionally, there might be a few B cells do not have any productive/sterile reads (in this case 2 cells). 
# # This will cause problems in merging the data. 
# # We can add them back to the matrix and proceed with the creation of this new assay, using the mergeIgHCountsToSeurat function:
# seurat_obj_integrated_comb <- mergeIgHCountsToSeurat(
#     igh_counts = subset_IGH2, 
#     SeuratObj = seurat_obj_integrated_comb,
#     # We will store the productive/sterile count matrix as a new
#     # assay called 'IGHC'
#     assay_name = "IGHC"
#   )
# 
# # Take a look 
# # ‘-S’: these correspond to the sterile transcripts;
# # ‘-P’: these correspond to the productive transcripts;
# # ‘-C’: these corresponds to transcripts which cannot be classified as sterile or productive
# slot(seurat_obj_integrated_comb[["IGHC"]], "counts")[1:9, 1:3]
# 
# # ------------------------------------------------------------------------------
# # Add finer annotations  
# # ------------------------------------------------------------------------------
# 
# # Load annotations 
# GCB_meta_GL_new <- readRDS("00_data/GCB_meta_GL_new.rds")
# 
# GCB_meta_GL_new$barcode <- glue("{GCB_meta_GL_new$Sample_name}_{rownames(GCB_meta_GL_new)}") %>% str_remove("_\\d+")
# 
# seurat_obj_integrated_comb$barcode 
# GCB_meta_GL_new$barcode
# 
# seurat_obj_integrated_comb$barcode  %>% str_remove("_[ATCG]{16}-\\d+$") %>% unique()
# GCB_meta_GL_new$barcode %>% str_remove("_[ATCG]{16}-\\d+$") %>% unique()
# 
# table(GCB_meta_GL_new$barcode %in% seurat_obj_integrated_comb$barcode)
# GCB_meta_GL_new$barcode[!(GCB_meta_GL_new$barcode %in% seurat_obj_integrated_comb$barcode)] %>% str_remove("_[ATCG]{16}-\\d+$") %>% table()
# 
# # Prep to add metadata 
# meta_add <- seurat_obj_integrated_comb[[]] %>%
#   tibble::rownames_to_column("cell") %>%
#   select(cell, barcode, L1_annotation) %>%
#   left_join(
#     GCB_meta_GL_new %>% select(barcode, L2_GCB_annotation, L3_GCB_annotation), 
#     by = "barcode"
#   ) %>%
#   tibble::column_to_rownames("cell") %>%
#   select(-barcode) %>% 
#   mutate(
#     L_annotation_Helene = case_when(
#       !is.na(L3_GCB_annotation) ~ L3_GCB_annotation, 
#       .default = L1_annotation 
#     )
#   ) 
# 
# # Add meta data to seurat object 
# seurat_obj_integrated_comb <- AddMetaData(seurat_obj_integrated_comb, meta_add)
# 
# saveRDS(seurat_obj_integrated_comb, "50_BACH2/out/seurat_obj_integrated_comb.rds")
seurat_obj_integrated_comb <- readRDS("50_BACH2/out/seurat_obj_integrated_comb.rds")

L_annotation_Helene_colors <- c(L1_colors, L3_GC_pretty_colors)

DimPlot(seurat_obj_integrated_comb, group.by = "L_annotation_Helene", reduction = "subset_umap.harmony", cols = L_annotation_Helene_colors)
ggsave(glue("{outdir}/DimPlot_L_annotation_Helene.png"))

# ------------------------------------------------------------------------------
# Look at results 
# ------------------------------------------------------------------------------

# Normalisation; 'IGHC' is the assay that holds the 
# sterile/productive count data
seurat_obj_integrated_comb <- NormalizeData(seurat_obj_integrated_comb, assay = "IGHC")

genes <- rownames(seurat_obj_integrated_comb[["IGHC"]])
# Let's say we want to plot only the sterile and productive counts
genes <- genes[which(grepl("-[SP]$", genes))]

DotPlot(seurat_obj_integrated_comb, assay = "IGHC", features = genes, group.by = "L_annotation_Helene") +
  RotatedAxis() + 
  plot_annotation(
    caption = "S: sterile transcripts, P: productive transcripts"
  )

ggsave(glue("{outdir}/IGHC_dotplot_annotation.png"))

# ‘-S’: these correspond to the sterile transcripts;
# ‘-P’: these correspond to the productive transcripts;
# ‘-C’: these corresponds to transcripts which cannot be classified as sterile or productive

# ------------------------------------------------------------------------------
# Inferring transitions using CSR and SHM information
# 
# sciCSR calculates a “CSR potential” which ranks cells from a naive state to a memory/plasma cell state, 
# based on their expression pattern of the sterile and productive IgH transcripts.
# ------------------------------------------------------------------------------

# Prep before running 
table(seurat_obj_integrated_comb$IGH_c_gene, seurat_obj_integrated_comb$L1_annotation)
table(is.na(seurat_obj_integrated_comb$IGH_c_gene))

# Merge IgD and IgM to IgM
seurat_obj_integrated_comb$IGH_c_gene[seurat_obj_integrated_comb$IGH_c_gene == "IGHD"] <- "IGHM"

# BCR types to keep 
keep_types <- c("singlet", "singlet_H", "multi_LC_diff_class", "multi_LC_same_class")
cells_keep <- colnames(seurat_obj_integrated_comb)[seurat_obj_integrated_comb$bcr_type %in% keep_types]
length(cells_keep)

seurat_obj_integrated_comb <- subset(seurat_obj_integrated_comb, cells = cells_keep)

# calculate CSR potential
seurat_obj_integrated_comb <- getCSRpotential(
  SeuratObj = seurat_obj_integrated_comb,
  # Here we use the isotype information from the merged 
  # repertoire annotations. Column 'IGH_c_gene' stores 
  # this information
  ighc_slot = "data",
  c_gene_anno_name = "IGH_c_gene", 
  reference_based = "human"
)

# get SHM frequency (= 1 - IGH_v_identity)
# we have merged the germline identity information from the
# repertoire to the Seurat object (column 'IGH_v_identity')
seurat_obj_integrated_comb <- getSHM(seurat_obj_integrated_comb, v_identity_anno_name = "IGH_v_identity")

# Prepare a copy for export: joined layers, scaled data, cells with an isotype
obj_h5ad <- seurat_obj_integrated_comb
obj_h5ad[["RNA"]] <- JoinLayers(obj_h5ad[["RNA"]])
obj_h5ad <- ScaleData(obj_h5ad, assay = "RNA")
obj_h5ad <- subset(obj_h5ad, cells = colnames(obj_h5ad)[!is.na(obj_h5ad$isotype)])

# Write the h5ad (write_h5ad_for_sciCSR() from functions script)
# convertSeuratToH5ad(seurat_obj_integrated_comb, assays = c("RNA"), "50_BACH2/out/seurat_obj_integrated_comb.h5ad")
h5ad_file <- "50_BACH2/out/seurat_obj_integrated_comb.h5ad"
write_h5ad_for_sciCSR(obj_h5ad, h5ad_file)

# Levels = Transitions from and to... 
celltype_order <- c(
  "Naive_B_cells",
  "Unconventional_B_cells",
  "GC_B_cells",
  "DZ_cycling",
  "DZ",
  "LZ",
  "LZ_transitioning",
  "BACH2_cycling",
  "BACH2+",
  "Pre_MBC_cycling",
  "Pre_MBC",
  "Memory_B_cells",
  "Pre_PB",
  "PCs"
)

splitAnnData(
  anndata_file = "50_BACH2/out/seurat_obj_integrated_comb.h5ad",
  split.by = "L_annotation_Helene", levels = celltype_order,
  conda_env = "scicsr_new"
)

# fit transitions based on CSR
# First Wild-Type
g_wt_csr <- fitTransitionModel(
  anndata_file = "50_BACH2/out/seurat_obj_integrated_comb.h5ad",
  mode = "pseudotime", pseudotime_key = "csr_pot",
  conda_env = "scicsr_new",
  do_pca = FALSE, do_neighbors = FALSE
)




# This information is stored in the 'Status' column in the 
# Seurat object metadata
# splitAnnData(
#   anndata_file = "hong_sampled_bcells_assay-RNA.h5ad",
#   split.by = "L_annotation_Helene", levels = c("WT", "IL23-/-")
# )
# 
# seurat_obj_integrated_comb$L_annotation_Helene %>% unique()


