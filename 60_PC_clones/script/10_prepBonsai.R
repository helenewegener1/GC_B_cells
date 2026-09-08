library(Seurat)
library(data.table)
library(glue) 

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# Load seurat object
seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs_May.rds")

# Subset seurat object to LP PCs
seurat_integrated_PC <- subset(seurat_integrated, str_detect(sample, "LP") & L1_annotation == "PCs")

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)
patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

# Load data and filter for LP PCs
resolve_LC_list <- lapply(resolve_LC_files, function(x){
  readRDS(glue("45_immcantation/out/rds/{x}")) %>% 
    filter(locus == "IGH" & L1_annotation == "PCs" & str_detect(sample_clean, "LP"))
}) %>% 
  setNames(patients)

resolve_LC_list <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  resolve_LC_list[[HH]] %>%
    group_by(clone_subgroup_id_90_similarity, locus) %>%
    ungroup()
  
}) %>% setNames(patients)

df_both <- bind_rows(resolve_LC_list)

# ------------------------------------------------------------------------------
# Prep outdir
# ------------------------------------------------------------------------------

outdir <- "60_PC_clones/bonsai_input/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Get metadata
# ------------------------------------------------------------------------------

# Update cell_id
seurat_integrated_PC[[]] <- seurat_integrated_PC[[]] %>% 
  mutate(
    cell_id = glue("{sample}_{rownames(.)}") %>% str_remove("_\\d")
  )

seurat_integrated_PC$cell_id %>% head()
df_both$cell_id %>% head()

table(seurat_integrated_PC$cell_id %in% df_both$cell_id)

# Subset to cells in df_both 
seurat_integrated_PC <- subset(seurat_integrated_PC, cell_id %in% df_both$cell_id)

table(seurat_integrated_PC$cell_id %in% df_both$cell_id)

# Add clone information to seurat object 
seurat_integrated_PC[[]] <- seurat_integrated_PC[[]] %>% left_join(df_both, by = "cell_id")

# Summarise 
PC_meta <- seurat_integrated_PC[[]] %>%
  select(
    cell_id, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo,
    scDblFinder.class, scDblFinder.score, 
    S.Score, G2M.Score, Phase, 
    patient, condition, inflammed, tissue, sample_clean, sample,
    clone_id_90_similarity, clone_subgroup_90_similarity, clone_subgroup_id_90_similarity, c_call, c_call_grouped
  ) %>% 
  rename(CellID = cell_id) 

rownames(PC_meta) <- NULL

nrow(PC_meta)

fwrite(PC_meta, glue("{outdir}/PC_meta_bonsai.tsv.gz"), sep = "\t", quote = FALSE)

# ------------------------------------------------------------------------------
# Get count matrix 
# ------------------------------------------------------------------------------

Layers(seurat_integrated_PC)

seurat_integrated_PC_merge <- JoinLayers(seurat_integrated_PC)

PC_count_matrix <- SeuratObject::LayerData(seurat_integrated_PC_merge, assay = "RNA", layer = "counts")
dim(PC_count_matrix)

colnames(PC_count_matrix) <- PC_meta$CellID

PC_count_matrix_dt <- as.data.table(as.matrix(PC_count_matrix), keep.rownames = "GeneID")

PC_count_matrix_dt[1:5, 1:5]
ncol(PC_count_matrix_dt)
fwrite(PC_count_matrix_dt, glue("{outdir}/PC_counts_bonsai.tsv.gz"), sep = "\t", quote = FALSE)


