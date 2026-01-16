library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)

source("08_seurat_QC_filtering/script/functions.R")

cellRversion <- "v8"

# Load data
seurat_obj_list <- readRDS(glue("07_seurat_QC/out/seurat_obj_QC_{cellRversion}.rds"))

# Initialize filtered list
seurat_obj_QC_filtered_list <- rep(0, length(seurat_obj_list)) %>% as.list()
names(seurat_obj_QC_filtered_list) <- names(seurat_obj_list)

# Check samples
names(seurat_obj_list)

# HH117-SI-MILF-INF-HLADR-AND-CD19 looks different than the others

# Extra plots
# for (sample_name in names(seurat_obj_list)){
#   
#   seurat_obj <- seurat_obj_list[[sample_name]]
#   seurat_obj <- pre_filter_pipeline(seurat_obj) # calculate the percent.mt and percent.hb
#   
#   FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
#   FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#   VlnPlot(seurat_obj, features = "percent.hb", layer = "counts")
#   
# }

################################### Rough QC ################################### 

sample_name <- names(seurat_obj_list)[[1]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name,
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[2]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name,
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[3]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[4]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[5]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[6]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[7]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 7500 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[8]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[9]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 7500 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[10]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[11]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[12]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj)

################################################################################

sample_name <- names(seurat_obj_list)[[13]]
seurat_obj <- seurat_obj_list[[sample_name]]

Idents(seurat_obj) <- "orig.ident"

# Filter on doublets and calculate QC metrics + plot pre_filter plots
seurat_obj <- pre_filter_pipeline(seurat_obj)

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 20)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, sample_name = sample_name, 
        version = "filtered", filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered)

rm(seurat_obj_list)

########################################## Export list of filtered Seurat objects ##########################################

names(seurat_obj_QC_filtered_list)
saveRDS(seurat_obj_QC_filtered_list, glue("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_list_{cellRversion}.rds"))

samples_names <- names(seurat_obj_QC_filtered_list)

memory.limit(size = 70000)

# Subset and save only singlets and only doublets in seperate objects 
seurat_obj_QC_filtered_singlets_list <- rep(0, length(seurat_obj_QC_filtered_list)) %>% as.list()
names(seurat_obj_QC_filtered_singlets_list) <- names(seurat_obj_QC_filtered_list)

seurat_obj_QC_filtered_doublets_list <- rep(0, length(seurat_obj_QC_filtered_list)) %>% as.list()
names(seurat_obj_QC_filtered_doublets_list) <- names(seurat_obj_QC_filtered_list)

for (sample_name in samples_names){
  
  # sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  seurat_obj <- seurat_obj_QC_filtered_list[[sample_name]]
  print(sample_name)
  
  # singlets 
  seurat_obj_singlets <- subset(seurat_obj, subset = scDblFinder.class == "singlet")
  seurat_obj_QC_filtered_singlets_list[[sample_name]] <- seurat_obj_singlets
  print(glue("N singlets: {ncol(seurat_obj_singlets)}"))
  
  # doublets
  seurat_obj_doublets <- subset(seurat_obj, subset = scDblFinder.class == "doublet")
  seurat_obj_QC_filtered_doublets_list[[sample_name]] <- seurat_obj_doublets
  print(glue("N doublets: {ncol(seurat_obj_doublets)}"))
  
  # Doublet percentage
  percentage_doublet <- round((ncol(seurat_obj_doublets)/ncol(seurat_obj) ) * 100, 2)
  print(glue("Doublet percentage: {percentage_doublet}%"))
  
  rm(seurat_obj, seurat_obj_singlets, seurat_obj_doublets)
  print("--------------------------------------")
  
}

saveRDS(seurat_obj_QC_filtered_singlets_list, glue("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_singlets_list_{cellRversion}.rds"))
saveRDS(seurat_obj_QC_filtered_doublets_list, glue("08_seurat_QC_filtering/out/seurat_obj_QC_filtered_doublets_list_{cellRversion}.rds"))

# seurat_obj_QC_filtered_doublets_list$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH`$scDblFinder.class %>% length()

