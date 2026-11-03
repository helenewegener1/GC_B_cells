# library(reticulate)
# use_condaenv("scicsr_new", required = TRUE)
# py_config()

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

outdir <- "50_BACH2/plot/05_sciCSR"
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
# samples <- seurat_obj_integrated_comb$sample_clean %>% unique() 

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

# Merge IgD and IgM to IgM
seurat_obj_integrated_comb$IGH_c_gene[seurat_obj_integrated_comb$IGH_c_gene == "IGHD"] <- "IGHM"


DotPlot(seurat_obj_integrated_comb, assay = "IGHC", features = genes, group.by = "L_annotation_Helene") +
  RotatedAxis() + 
  plot_annotation(
    caption = "S: sterile transcripts, P: productive transcripts"
  )
ggsave(glue("{outdir}/IGHC_SP_dotplot_annotation.png"))

# ‘-S’: these correspond to the sterile transcripts;
# ‘-P’: these correspond to the productive transcripts;
# ‘-C’: these corresponds to transcripts which cannot be classified as sterile or productive

DotPlot(seurat_obj_integrated_comb, assay = "IGHC",
        features = rownames(seurat_obj_integrated_comb[["IGHC"]]),
        group.by = "L_annotation_Helene") + RotatedAxis()
ggsave(glue("{outdir}/IGHC_SPC_dotplot_annotation.png"))

# ------------------------------------------------------------------------------
# Check fractions of transcript types in total and across isotypes. 
# ------------------------------------------------------------------------------

# Fraction of cells without any IGHC count, per cell type and per sample
ighc_total <- colSums(LayerData(seurat_obj_integrated_comb, assay = "IGHC", layer = "counts"))
seurat_obj_integrated_comb$ighc_total <- ighc_total[colnames(seurat_obj_integrated_comb)]
round(tapply(seurat_obj_integrated_comb$ighc_total == 0,
             seurat_obj_integrated_comb$L_annotation_Helene, mean), 2)
table(seurat_obj_integrated_comb$sample, no_ighc = seurat_obj_integrated_comb$ighc_total == 0)

# Fraction of S, P and C transcripts per cell type.
ighc <- LayerData(seurat_obj_integrated_comb, assay = "IGHC", layer = "counts")
ct <- seurat_obj_integrated_comb$L_annotation_Helene
type_sums <- sapply(c("S", "P", "C"), function(t)
  tapply(colSums(ighc[grepl(paste0("-", t, "$"), rownames(ighc)), ]), ct, sum))

prop.table(type_sums, margin = 1) %>% 
  as.data.frame() %>% 
  na.omit() %>% 
  rownames_to_column("cell_type") %>% 
  pivot_longer(cols = !cell_type) %>% 
  ggplot(aes(x = cell_type, y = value, fill = name)) + 
  geom_col() + 
  theme_bw() + 
  labs(
    title = "Fraction of types of transcripts per cell type", 
    x = "", 
    y = "Fraction", 
    fill = "Transcript type",
    subtitle = "S: sterile transcripts, P: productive transcripts, C: unclassified"
  ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(
    breaks = scales::breaks_width(0.2),
    # minor_breaks = scales::breaks_width(minor_breaks_width),
    limits = c(0, 1)
  )
ggsave(glue("{outdir}/fraction_of_transcript_type.png"))

## IgG1 specifically, per cell type: fraction of cells with any IGHG1-S/P/C count
isotypes <- rownames(ighc) %>% str_remove("-[SPC]") %>% unique()

for (isotype in isotypes){
  # isotype <- "IGHM"
  l <- paste0(isotype, c("-S", "-P", "-C"))
  sapply(l, function(f) tapply(ighc[f, ] > 0, ct, mean)) %>% 
    as.data.frame() %>% 
    na.omit() %>% 
    rownames_to_column("cell_type") %>% 
    pivot_longer(cols = !cell_type) %>% 
    ggplot(aes(x = cell_type, y = value, fill = name)) + 
    geom_col() + 
    theme_bw() + 
    labs(
      title = glue("{isotype} fraction of types of transcripts per cell type"), 
      x = "", 
      y = "Fraction", 
      fill = "Transcript type",
      subtitle = "S: sterile transcripts, P: productive transcripts, C: unclassified"
    ) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(
      breaks = scales::breaks_width(0.2),
      # minor_breaks = scales::breaks_width(minor_breaks_width),
      limits = c(0, 1)
    )
  
  ggsave(glue("{outdir}/{isotype}_fraction_of_transcript_type.png"))
  
}

isotypes <- rownames(ighc) %>% str_remove("-[SPC]$") %>% unique()

df_frac <- map_dfr(isotypes, function(isotype) {
  l <- paste0(isotype, c("-S", "-P", "-C"))
  sapply(l, function(f) tapply(ighc[f, ] > 0, ct, mean)) %>%
    as.data.frame() %>%
    na.omit() %>%
    rownames_to_column("cell_type") %>%
    pivot_longer(cols = !cell_type, names_to = "feature", values_to = "fraction") %>%
    mutate(isotype = isotype,
           transcript_type = str_extract(feature, "[SPC]$"))
}) %>%
  mutate(
    isotype = factor(isotype, levels = isotypes),               # IGH locus order
    cell_type = factor(cell_type),
    transcript_type = factor(transcript_type, levels = c("S", "P", "C"))
  )

ggplot(df_frac, aes(x = cell_type, y = fraction, fill = transcript_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  facet_wrap(~ isotype, ncol = 2) +
  theme_bw() +
  labs(
    title = "Fraction of cells with sterile, productive or unclassified transcripts",
    x = "", y = "Fraction of cells", fill = "Transcript type",
    subtitle = "S: sterile transcripts, P: productive transcripts, C: unclassified"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = scales::breaks_width(0.2), limits = c(0, 1))

ggsave(glue("{outdir}/fraction_of_transcript_type_per_isotype.png"),
       width = 10, height = 7.5)


# ------------------------------------------------------------------------------
# Sterile isotypes - filter in what we can use
# For BCR 
# Isotype switch order: IGHM → IGHD → IGHG3 → IGHG1 → IGHA1 → IGHG2 → IGHG4 → GHE → IGHA2
# ------------------------------------------------------------------------------

# Define switch order
locus <- c("IGHM", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")

# Load clone data 
resolve_LC <- readRDS(glue("45_immcantation/out/rds/06_resolve_LC_germlined.rds"))
patients <- resolve_LC %>% names() 
df_heavy <- resolve_LC %>% bind_rows() %>% filter(locus == "IGH")

df_heavy$cell_id %>% table()
df_heavy$cell_id %>% str_split_i("_", 1) %>% unique()

df_heavy_meta <- df_heavy %>% select(cell_id, c_call, c_call_grouped, clone_subgroup_id_90_similarity) %>% 
  mutate(
    c_call_sciCSR = ifelse(c_call_grouped == "IGHM/D", "IGHM", c_call_grouped)
  )

# table(df_heavy_meta$cell_id %in% L1_annotation_metadata_clean$cell_id) 
# df_heavy_meta <- df_heavy_meta %>% left_join(L1_annotation_metadata_clean, by = "cell_id")

# Filter seurat_obj to only include cells that are also in df_heavy
seurat_obj_integrated_comb$cell_id <- colnames(seurat_obj_integrated_comb) %>% str_remove("_\\d+")
seurat_obj_integrated_comb$cell_id_sample <- seurat_obj_integrated_comb$sample %>% str_split_i("_", 1)
seurat_obj_integrated_comb$cell_id <- glue("{seurat_obj_integrated_comb$cell_id_sample}_{seurat_obj_integrated_comb$cell_id}")

table(seurat_obj_integrated_comb$cell_id %in% df_heavy$cell_id)

# Get sterile counts per cell from the IGHC assay
ighc <- LayerData(seurat_obj_integrated_comb, assay = "IGHC", layer = "counts")
df_sterile <- as.matrix(ighc[grep("-S$", rownames(ighc)), , drop = FALSE]) %>%
  t() %>%
  as.data.frame()

if (table(rownames(df_sterile) == colnames(seurat_obj_integrated_comb))){
  df_sterile$cell_id <- seurat_obj_integrated_comb$cell_id
}

# Get annotations 
seurat_obj_integrated_comb_meta <- seurat_obj_integrated_comb@meta.data %>% 
  select(cell_id, L1_annotation, L3_GCB_annotation, L_annotation_Helene)

# Combine and reshape: one row per cell and sterile transcript
df_sterile_cells <- df_heavy_meta %>%
  inner_join(df_sterile, by = "cell_id") %>%
  inner_join(seurat_obj_integrated_comb_meta, by = "cell_id") %>%
  pivot_longer(cols = ends_with("-S"), names_to = "sterile_feature", values_to = "n_molecules") %>%
  mutate(
    patient = cell_id %>% str_split_i("-", 1),
    sterile_isotype = str_remove(sterile_feature, "-S$"),
    type_switch = case_when(
      match(c_call_sciCSR, locus) == match(sterile_isotype, locus) ~ "same",
      match(c_call_sciCSR, locus) < match(sterile_isotype, locus) ~ "sequential",
      match(c_call_sciCSR, locus) > match(sterile_isotype, locus) ~ "reverse"
    ),
    informative = match(c_call_sciCSR, locus) < match(sterile_isotype, locus)   # upstream cells only
  ) %>% 
  filter(!is.na(c_call_sciCSR))

saveRDS(df_sterile_cells, "50_BACH2/out/df_sterile_cells.rds")

# How many cells does an informative switch (following the isotype switch order)
df_sterile_cells$informative %>% table(useNA = "always")
df_sterile_cells$type_switch %>% table(useNA = "always")

outdir_diagnosis <- glue("{outdir}/diagnosis")
dir.create(outdir_diagnosis, recursive = TRUE, showWarnings = FALSE)

# How many molecules 
df_sterile_cells %>% 
  ggplot(aes(x = n_molecules)) + 
  geom_histogram() + 
  facet_wrap(vars(type_switch)) + 
  scale_y_log10() + 
  theme_bw() + 
  labs(
    title = "N molecule counts"
  )

ggsave(glue("{outdir_diagnosis}/n_molecules_histogram.png"), width = 14, height = 5)

# Investigate non-informative rows 
plot_df <- df_sterile_cells %>%
  mutate(
    c_call_sciCSR = factor(c_call_sciCSR, levels = locus),
    sterile_feature = factor(sterile_feature, levels = paste0(locus, "-S")),
    type_switch = factor(type_switch, levels = c("sequential", "same", "reverse"))
  )

# 1. Summed molecules
plot_df %>%
  group_by(c_call_sciCSR, sterile_feature, type_switch) %>%
  summarise(sum = sum(n_molecules), .groups = "drop") %>%
  ggplot(aes(x = c_call_sciCSR, y = sterile_feature, fill = sum)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sum), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  facet_wrap(vars(type_switch)) +
  labs(
    x = "Productive isotype (BCR)", y = "Sterile transcript",
    fill = "Total molecule count",
    title = "Productive vs sterile transcripts",
    subtitle = "Numbers are summed molecules, not cells"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(glue("{outdir_diagnosis}/sterile_vs_productive_molecules.png"), width = 14, height = 5)

# 2. Cells with at least n_molecules_min sterile molecules
for (n_molecules_min in c(1, 2, 3)){
  
  # n_molecules_min <- 2
  plot_df %>%
    filter(n_molecules >= n_molecules_min) %>%
    count(c_call_sciCSR, sterile_feature, type_switch) %>%
    ggplot(aes(x = c_call_sciCSR, y = sterile_feature, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 4) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    facet_wrap(vars(type_switch), drop = FALSE) +
    labs(
      x = "Productive isotype (BCR)", y = "Sterile transcript",
      fill = "Cell count",
      title = "Productive vs sterile transcripts",
      subtitle = glue("Numbers are cells with at least {n_molecules_min} sterile molecules")
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(glue("{outdir_diagnosis}/sterile_vs_productive_cells_min{n_molecules_min}.png"), width = 14, height = 5)
  
}

# Investigate non-informative rows - per cell type 
plot_df %>%
  group_by(c_call_sciCSR, sterile_feature, type_switch, L_annotation_Helene) %>%
  summarise(sum = sum(n_molecules), .groups = "drop") %>%
  filter(sum > 0) %>%
  ggplot(aes(x = c_call_sciCSR, y = sterile_feature, fill = sum)) +
  geom_tile(aes(colour = type_switch), linewidth = 0.6) +
  geom_text(aes(label = sum), size = 2.5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_colour_manual(values = c(sequential = "forestgreen", same = "grey50", reverse = "purple")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  facet_wrap(vars(L_annotation_Helene), ncol = 4) +
  labs(
    x = "Productive isotype (BCR)", y = "Sterile transcript",
    fill = "Total molecule count", colour = "Switch type",
    title = "Productive vs sterile transcripts per cell type",
    subtitle = "Numbers are summed molecules, not cells"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave(glue("{outdir_diagnosis}/cell_types_sterile_vs_productive_molecules.png"), width = 12, height = 7)


# 2. Cells with at least n_molecules_min sterile molecules
for (n_molecules_min in c(1, 2, 3)){
  
  # n_molecules_min <- 2
  plot_df %>%
    filter(n_molecules >= n_molecules_min) %>%
    count(c_call_sciCSR, sterile_feature, type_switch, L_annotation_Helene) %>%
    ggplot(aes(x = c_call_sciCSR, y = sterile_feature, fill = n)) +
    geom_tile(aes(colour = type_switch), linewidth = 0.6) +
    geom_text(aes(label = n), size = 2.5) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    scale_colour_manual(values = c(sequential = "forestgreen", same = "grey50", reverse = "purple")) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    facet_wrap(vars(L_annotation_Helene), ncol = 4) +
    labs(
      x = "Productive isotype (BCR)", y = "Sterile transcript",
      fill = "Cell count", colour = "Switch type",
      title = "Productive vs sterile transcripts per cell type",
      subtitle = glue("Numbers are cells with at least {n_molecules_min} sterile molecules")
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(glue("{outdir_diagnosis}/cell_types_sterile_vs_productive_cells_min{n_molecules_min}.png"), width = 12, height = 7)
  
}

# Test thresholds for n_molecules 
thresholds <- 0:10

df_thresh <- map_dfr(thresholds, function(t) {
  df_sterile_cells %>%
    filter(n_molecules >= t) %>%
    count(type_switch) %>%
    mutate(min_molecules = t)
}) %>%
  mutate(
    type_switch = factor(type_switch, levels = c("sequential", "same", "reverse")),
    min_molecules = factor(min_molecules, levels = thresholds,
                           labels = c("All", paste0("\u2265", thresholds[-1])))
  )

switch_colors <- c(sequential = "forestgreen", same = "grey50", reverse = "purple")

ggplot(df_thresh, aes(x = min_molecules, y = n, fill = type_switch)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = switch_colors, na.value = "grey85") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 100000),
                     labels = scales::comma,
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Minimum sterile molecules", y = "Cell–transcript pairs (pseudo-log)",
       fill = "Switch type") +
  theme_bw(base_size = 12)

ggsave(glue("{outdir_diagnosis}/sterile_switch_type_by_threshold.png"), width = 11, height = 5.5)

# Clean 
df_sterile_cells_inf <- df_sterile_cells %>% filter(informative)

df_sterile_cells_inf$c_call_sciCSR %>% table()
df_sterile_cells_inf$patient %>% table()
df_sterile_cells_inf$L_annotation_Helene %>% table()
table(df_sterile_cells_inf$L_annotation_Helene, df_sterile_cells_inf$L1_annotation)

# Inspect 
df_sterile_cells_inf %>% filter(c_call_sciCSR == "IGHM" & n_molecules > 1) %>% count(patient, sterile_isotype)

# Heatmap of current and sterile isotype - per cell type 
cell_types <- df_sterile_cells_inf$L_annotation_Helene %>% unique()

lapply(cell_types, function(cell_type){
  
  # cell_type <- "DZ_cycling"
  df_plot <- df_sterile_cells_inf %>% 
    filter(
      n_molecules >= 1,
      L_annotation_Helene == cell_type
    ) %>% 
    select(c_call_sciCSR, sterile_isotype, patient) %>% 
    table() %>% 
    as.data.frame()
  # df_plot
  
  ggplot(df_plot, aes(x = c_call_sciCSR, y = sterile_isotype, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      x = "Current isotype (BCR)",
      y = "Sterile isptype (sciCSR)",
      fill = "Cell count", 
      title = glue("{cell_type}: Current and sterile isotype")
    ) +
    theme_bw() + 
    facet_wrap(vars(patient))
  
  ggsave(glue("{outdir}/current_to_sterile_isotype/current_to_sterile_isotype_{cell_type}.png"), height = 3.5, width = 11)
  
})

# Heatmap of current and sterile isotype - per cell type 
# cell_types <- df_sterile_cells_inf$L_annotation_Helene %>% unique()

HH <- "HH119"
clone_nr <- 1
top_clone <- df_sterile_cells_inf %>%
  filter(
    patient == HH,
    L1_annotation == "GC_B_cells"
  ) %>%
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>%
  slice(clone_nr) %>%
  pull(clone_subgroup_id_90_similarity)

# HH <- "HH117"
# clone_nr <- 3
# top_clone <- "1881_1"

# cell_type <- "DZ_cycling"
df_plot <- df_sterile_cells_inf %>% 
  filter(
    n_molecules >= 1, 
    clone_subgroup_id_90_similarity == top_clone, 
    patient == HH
  ) %>% 
  select(c_call_sciCSR, sterile_isotype, patient) %>% 
  table() %>% 
  as.data.frame()

ggplot(df_plot, aes(x = c_call_sciCSR, y = sterile_isotype, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Current isotype (BCR)",
    y = "Sterile isptype (sciCSR)",
    fill = "Cell count", 
    title = glue("{HH} top {clone_nr} clone ({top_clone}): Current and sterile isotype")
  ) +
  theme_bw() + 
  facet_wrap(vars(patient))

ggsave(glue("{outdir}/current_to_sterile_isotype/current_to_sterile_isotype_{HH}_top{clone_nr}clone.png"), height = 3.5, width = 5.5)

# Get cell ids 
df_plot_cell_ids <- df_sterile_cells_inf %>% 
  filter(
    n_molecules >= 1, 
    clone_subgroup_id_90_similarity == top_clone, 
    patient == HH
  ) %>% 
  pull(cell_id)

# Get cell ids for cells that show sterile transcript of IgE
df_IGES <- df_sterile_cells_inf %>% 
  filter(
    clone_subgroup_id_90_similarity == top_clone, 
    patient == HH, 
    sterile_feature == "IGHE-S"
  ) %>% 
  filter(n_molecules > 0) 

df_IGES %>% 
  ggplot(aes(y = n_molecules, x = 1)) + 
  geom_jitter(alpha = 0.5) + 
  theme_bw() + 
  labs(
    title = "HH119 large clone: N molecules of IGHE-S across cells",
    subtitle = "N molecules > 0 are shown"
  )
ggsave(glue("{outdir}/HH119_large_clone/IGES_Nmolecules_distribution.png"))

df_IGE_cells <- df_IGES$cell_id 
saveRDS(df_IGE_cells, "50_BACH2/out/HH119_large_clone_IGES_cells.rds")
  
# -----
# Map nodes on tree to cell ids. 
# -----
library(ape)
library(ggnewscale)
seq_dir <- readRDS("48_GCtree/out/rds/seq_dir.rds")

# Get clone files 
clone_full_name <- glue("{HH}_clone_nr_{clone_nr}_clone_{top_clone}")

sequence_ids <- seq_dir[[clone_full_name]]
# table(df_clone$sequence_alignment == sequence_ids)

# filter for clone 
df_clone <- df_heavy %>% 
  filter(
    patient_id == HH, 
    clone_subgroup_id_90_similarity == top_clone
  ) %>% 
  mutate(
    sequence_id_gctree = names(sequence_ids)[match(sequence_alignment, sequence_ids)]
  )

# Get clone files
tree_file <- glue("48_GCtree/plot/{clone_full_name}/{clone_full_name}.inference.1.nk")
idmap_file <- glue("48_GCtree/out/{clone_full_name}/idmap.txt")
seq_cell_id_dir <- readRDS("48_GCtree/out/rds/seq_cell_id_dir.rds") 
seq_cell_id_df <- data.frame(
  sequence_id = names(seq_cell_id_dir[[clone_full_name]]), 
  cell_id = seq_cell_id_dir[[clone_full_name]]
)

# Load gctree
tree <- read.tree(tree_file)

# Distance between ALL nodes (tips + internal)
node_labels <- c(tree$tip.label, tree$node.label)
node_dist <- dist.nodes(tree)
rownames(node_dist) <- colnames(node_dist) <- node_labels

# Nodes that count as observed: sequences plus the germline
is_observed <- function(x) str_detect(x, "^(seq[0-9]+|GL)$")
obs_nodes <- node_labels[is_observed(node_labels)]
obs_dist <- node_dist[obs_nodes, obs_nodes]

# Parent of each node: tree$edge holds (parent index, child index) pairs
edge_df <- as_tibble(tree$edge, .name_repair = "minimal") %>%
  set_names(c("parent_idx", "child_idx")) %>%
  mutate(parent = node_labels[parent_idx], child = node_labels[child_idx])
parent_of <- set_names(edge_df$parent, edge_df$child)

# Walk up until the next observed node (a sequence or the germline)
get_nearest_observed_ancestor <- function(node, parent_of) {
  current <- parent_of[node]
  while (!is.na(current) && !is_observed(current)) {
    current <- parent_of[current]
  }
  if (is.na(current)) NA_character_ else unname(current)
}

ancestor_map <- set_names(
  map_chr(obs_nodes, ~ get_nearest_observed_ancestor(.x, parent_of)),
  obs_nodes
)

# One row per observed node: its nearest observed ancestor and the distance between them
ancestor_df <- tibble(node = obs_nodes, ancestor = ancestor_map[obs_nodes]) %>%
  mutate(distance_to_ancestor = if_else(
    is.na(ancestor), NA_real_,
    node_dist[cbind(node, coalesce(ancestor, node))]
  ))

# Expand collapsed nodes back to individual sequence_ids (germline kept)
idmap <- read_csv(idmap_file, col_names = c("node", "sequence_id"), show_col_types = FALSE) %>%
  filter(node %in% obs_nodes) %>%
  separate_longer_delim(sequence_id, delim = ":")

# One row per sequence: tree node -> GCtree sequence ID -> cell metadata
node_seq_df <- idmap %>%
  left_join(
    df_clone %>% select(sequence_id, cell_id, c_call, L1_annotation, sample_clean_fol),
    by = "sequence_id"
  ) %>%
  left_join(ancestor_df, by = "node")

# One row per node: everything it contains, collapsed into readable strings
node_table <- tibble(node = node_labels) %>%
  mutate(node_type = case_when(
    node == "GL" ~ "germline",
    str_detect(node, "^seq[0-9]+$") ~ "observed",
    TRUE ~ "unobserved ancestor"
  )) %>%
  left_join(
    node_seq_df %>%
      group_by(node) %>%
      summarise(
        n_sequences  = n_distinct(sequence_id),
        sequence_ids = paste(sort(unique(sequence_id)), collapse = ", "),
        cell_ids     = paste(sort(unique(na.omit(cell_id))), collapse = ", "),
        isotypes     = paste(sort(unique(na.omit(c_call))), collapse = ", "),
        cell_types   = paste(sort(unique(na.omit(L1_annotation))), collapse = ", "),
        samples      = paste(sort(unique(na.omit(sample_clean_fol))), collapse = ", "),
        .groups = "drop"
      ),
    by = "node"
  ) %>%
  left_join(ancestor_df, by = "node")

write_tsv(node_table, glue("{outdir}/{clone_full_name}_node_table.tsv"))


node_ancestor <- tibble(node = obs_nodes, ancestor = ancestor_map[obs_nodes]) %>%
  left_join(
    idmap,
    by = "node"
  ) %>% 
  left_join(
    seq_cell_id_df, by = "sequence_id"
  ) %>% 
  left_join(
    df_clone %>% select(cell_id, c_call_grouped, L1_annotation, sample_clean_fol),
    by = "cell_id"
  )
node_ancestor

node_ancestor %>% filter(ancestor == "GL")
node_ancestor %>% filter(ancestor == "seq33") # 10 LP PCs


# Nodes with sterile transcripts 
node_ancestor %>% 
  filter(cell_id %in% df_plot_cell_ids)

# node to sequence_id map
node_to_seq_id <- node_ancestor %>% 
  filter(cell_id %in% df_plot_cell_ids) %>% 
  select(node, sequence_id)

# Get sterile transcripts of the nodes
df_sterile_cells_inf %>% filter(cell_id %in% df_plot_cell_ids, n_molecules >= 1) %>% 
  filter(cell_id %in% df_plot_cell_ids) %>% 
  left_join(seq_cell_id_df, by = "cell_id") %>% 
  left_join(node_to_seq_id, by = "sequence_id") %>% 
  select(node, sequence_id, c_call_sciCSR, sterile_isotype) 
  



# # ------------------------------------------------------------------------------
# # Inferring transitions using CSR and SHM information
# # 
# # sciCSR calculates a “CSR potential” which ranks cells from a naive state to a memory/plasma cell state, 
# # based on their expression pattern of the sterile and productive IgH transcripts.
# # ------------------------------------------------------------------------------
# 
# # Prep before running 
# table(seurat_obj_integrated_comb$IGH_c_gene, seurat_obj_integrated_comb$L1_annotation)
# table(seurat_obj_integrated_comb$IGH_c_gene, seurat_obj_integrated_comb$L_annotation_Helene)
# table(is.na(seurat_obj_integrated_comb$IGH_c_gene))
# 
# # BCR types to keep 
# keep_types <- c("singlet", "singlet_H", "multi_LC_diff_class", "multi_LC_same_class")
# cells_keep <- colnames(seurat_obj_integrated_comb)[seurat_obj_integrated_comb$bcr_type %in% keep_types]
# length(cells_keep)
# 
# seurat_obj_integrated_comb <- subset(seurat_obj_integrated_comb, cells = cells_keep)
# 
# # calculate CSR potential
# seurat_obj_integrated_comb <- getCSRpotential(
#   SeuratObj = seurat_obj_integrated_comb,
#   # Here we use the isotype information from the merged 
#   # repertoire annotations. Column 'IGH_c_gene' stores 
#   # this information
#   ighc_slot = "data",
#   c_gene_anno_name = "IGH_c_gene", 
#   reference_based = "human"
# )
# 
# # get SHM frequency (= 1 - IGH_v_identity)
# # we have merged the germline identity information from the
# # repertoire to the Seurat object (column 'IGH_v_identity')
# seurat_obj_integrated_comb <- getSHM(seurat_obj_integrated_comb, v_identity_anno_name = "IGH_v_identity")
# 
# # Prepare a copy for export: joined layers, scaled data, cells with an isotype
# obj_h5ad <- seurat_obj_integrated_comb
# obj_h5ad[["RNA"]] <- JoinLayers(obj_h5ad[["RNA"]])
# obj_h5ad <- ScaleData(obj_h5ad, assay = "RNA")
# obj_h5ad <- subset(obj_h5ad, cells = colnames(obj_h5ad)[!is.na(obj_h5ad$isotype)])
# 
# # Write the h5ad (write_h5ad_for_sciCSR() from functions script)
# # convertSeuratToH5ad(seurat_obj_integrated_comb, assays = c("RNA"), "50_BACH2/out/seurat_obj_integrated_comb.h5ad")
# h5ad_file <- "50_BACH2/out/seurat_obj_integrated_comb.h5ad"
# write_h5ad_for_sciCSR(obj_h5ad, h5ad_file)
# 
# # Levels = Transitions from and to... 
# celltype_order <- c(
#   "Naive_B_cells",
#   "Unconventional_B_cells",
#   "GC_B_cells",
#   "DZ_cycling",
#   "DZ",
#   "LZ",
#   "LZ_transitioning",
#   "BACH2_cycling",
#   "BACH2+",
#   "Pre_MBC_cycling",
#   "Pre_MBC",
#   "Memory_B_cells",
#   "Pre_PB",
#   "PCs"
# )
# 
# splitAnnData(
#   anndata_file = "50_BACH2/out/seurat_obj_integrated_comb.h5ad",
#   split.by = "L_annotation_Helene", levels = celltype_order,
#   conda_env = "scicsr_new"
# )
# 
# # fit transitions based on CSR
# # First Wild-Type
# g_wt_csr <- fitTransitionModel(
#   anndata_file = "50_BACH2/out/seurat_obj_integrated_comb.h5ad",
#   mode = "pseudotime", pseudotime_key = "csr_pot",
#   conda_env = "scicsr_new",
#   do_pca = FALSE, do_neighbors = FALSE
# )
# 
# 
# 
# 
# # This information is stored in the 'Status' column in the 
# # Seurat object metadata
# # splitAnnData(
# #   anndata_file = "hong_sampled_bcells_assay-RNA.h5ad",
# #   split.by = "L_annotation_Helene", levels = c("WT", "IL23-/-")
# # )
# # 
# # seurat_obj_integrated_comb$L_annotation_Helene %>% unique()
# 
# 
