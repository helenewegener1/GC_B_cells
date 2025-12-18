getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(glue)
library(ggplot2)
library(patchwork)
library(readxl)
# library(ztable)
library(pheatmap)

# Load data
# seurat_obj_list <- readRDS("07_seurat_QC/out/seurat_obj_QC.rds")
seurat_obj_list <- readRDS("08_seurat_QC_filtering/out/seurat_obj_roughQC_list.rds")
seurat_obj <- seurat_obj_list[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

Gina_seurat_obj <- readRDS("00_data/Gina_HH117_PP_broadAnn.rds")

# Explore data
dim(seurat_obj@assays$RNA)
dim(Gina_seurat_obj@assays$RNA)

DimPlot(Gina_seurat_obj, group.by = )

# Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 20, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

DimPlot(seurat_obj)

# Clusters
Gina_seurat_obj@meta.data$Celltype 
Gina_seurat_obj


# Gina seurat object 

Gina_seurat_obj@meta.data$Celltype
Gina_seurat_obj@meta.data$ADT_maxID
Gina_seurat_obj@meta.data$ADT_classification
Gina_seurat_obj@meta.data$ADT_classification.global

x_1 <- table(Gina_seurat_obj@meta.data$Celltype, Gina_seurat_obj@meta.data$ADT_classification.global)
x_2 <- table(Gina_seurat_obj@meta.data$Celltype, Gina_seurat_obj@meta.data$ADT_classification)
x_3 <- table(Gina_seurat_obj@meta.data$Celltype, Gina_seurat_obj@meta.data$ADT_maxID)

heatmap(x_1, Colv = NA, Rowv = NA)
# heatmap(x_2, Colv = NA, Rowv = NA)
heatmap(x_3, Colv = NA, Rowv = NA)

# Plot
x_1 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("global_class") %>% 
  pivot_longer(cols = c("Mem.B cells", "Tfh cells", "GC.B cells", "Naive.B cells", "PBs", "DCs")) %>% 
  ggplot(aes(y = name, x = value, fill = global_class)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  scale_fill_manual(values = c("blue4", "red3", "lightblue")) + 
  labs(x = "Count", y = "")
ggsave("07_seurat_QC/plot/ADT_explore/global_class.png", width = 12, height = 8)


# x_2 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("class") %>% 
#   # pivot_longer(cols = c("Mem.B cells", "Tfh cells", "GC.B cells", "Naive.B cells", "PBs", "DCs")) %>% 
#   ggplot(aes(x = DCs, y = class)) + 
#   geom_col(position = "dodge") + 
#   theme_bw()  
#   # scale_fill_manual(values = c("blue4", "red3", "lightblue"))

x_3 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("maxID") %>% 
  ggplot(aes(x = DCs, y = maxID)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  labs(x = "Count")
ggsave("07_seurat_QC/plot/ADT_explore/maxID_DCs.png", width = 12, height = 8)


x_3 %>% t() %>% as.data.frame.matrix() %>% rownames_to_column("maxID") %>% 
  pivot_longer(cols = c("Mem.B cells", "Tfh cells", "GC.B cells", "Naive.B cells", "PBs", "DCs")) %>%
  ggplot(aes(y = name, x = value, fill = maxID)) + 
  geom_col(position = "dodge", ) + 
  theme_bw() + 
  labs(x = "Count", y = "")
ggsave("07_seurat_QC/plot/ADT_explore/maxID_all_cells.png", width = 12, height = 8)


# Doublets
df_doublets <- Gina_seurat_obj@meta.data %>% 
  filter(ADT_classification.global == "Doublet") %>% 
  select(starts_with("ADT"))

pdf("07_seurat_QC/plot/ADT_explore/doublets_heatmap.pdf", width = 12, height = 8)
table(df_doublets$ADT_maxID, df_doublets$ADT_secondID) %>% 
  heatmap(Colv = NA, Rowv = NA, main = "ADT_maxID vs ADT_secondID for Doublets", xlab = "ADT_maxID", ylab = "ADT_secondID")
dev.off()

Gina_seurat_obj@meta.data$Celltype %>% table()

#
HTOHeatmap(Gina_seurat_obj, assay = "ADT", ncells = 5000)


# Investigate 

df <- Gina_seurat_obj[["ADT"]]$counts %>% 
  as.data.frame() %>% 
  rownames_to_column("follicles") %>% 
  pivot_longer(cols = colnames(Gina_seurat_obj), names_to = "cell", values_to = "value")

df_meta_data <- Gina_seurat_obj@meta.data %>% select(ADT_classification, ADT_maxID, ADT_classification.global) %>% rownames_to_column("cell")

df <- df %>% left_join(df_meta_data, by = "cell", relationship = "many-to-many")

head(df)

cell_doublets <- Gina_seurat_obj@meta.data %>% filter(ADT_classification.global == "Doublet") %>% rownames()

cell_subset_list <- list(
  x1 = cell_doublets[1:5],
  x2 = cell_doublets[5:10], 
  x3 = cell_doublets[10:15]
)

x <- "x3"

df %>% 
  filter(cell %in% cell_subset_list[[x]]) %>%
  filter(ADT_classification.global == "Doublet") %>%
  ggplot(aes(x = follicles, 
             y = value,
             color = ADT_classification)) + 
  # geom_col(position = "dodge") + 
  geom_point() +
  geom_line(aes(group = cell)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  labs(
    x = "", 
    title = "5 cells where ADT_classification.global = 'Doublet'"
  )

ggsave(glue("07_seurat_QC/plot/ADT_explore/doublets_lineplot_{x}.png"), width = 12, height = 8)

# log value 
df %>% 
  filter(cell %in% cell_subset_list[[x]]) %>%
  filter(ADT_classification.global == "Doublet") %>%
  ggplot(aes(x = follicles, 
             y = log(value),
             color = ADT_classification)) + 
  # geom_col(position = "dodge") + 
  geom_point() +
  geom_line(aes(group = cell)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  labs(
    x = "", 
    title = "5 cells where ADT_classification.global = 'Doublet'"
  )

ggsave(glue("07_seurat_QC/plot/ADT_explore/doublets_lineplot_{x}_log.png"), width = 12, height = 8)



# 
DefaultAssay(Gina_seurat_obj)
Idents(Gina_seurat_obj) <- Gina_seurat_obj@meta.data$ADT_maxID
FeatureScatter(Gina_seurat_obj, feature1 = "Fol-1", feature2 = "Fol-2")
FeatureScatter(Gina_seurat_obj, feature1 = "Fol-6", feature2 = "Fol-8")
FeatureScatter(Gina_seurat_obj, feature1 = "Fol-2", feature2 = "Fol-9")

# ADT QC 
VlnPlot(Gina_seurat_obj, features = "nCount_ADT", layer = "counts")
Gina_seurat_obj_subset <- subset(Gina_seurat_obj, subset = nCount_ADT > 100 & nCount_ADT < quantile(Gina_seurat_obj$nCount_ADT, 0.99))

# N cells before and after filtering
ncol(Gina_seurat_obj[["ADT"]])
ncol(Gina_seurat_obj_subset[["ADT"]])

Gina_seurat_obj_subset <- NormalizeData(Gina_seurat_obj_subset, assay = "ADT", normalization.method = "CLR")
Gina_seurat_obj_subset <- HTODemux(Gina_seurat_obj_subset, assay = "ADT", positive.quantile = 0.999)

Gina_seurat_obj_subset@meta.data$ADT_classification.global %>% table()
Gina_seurat_obj_subset@meta.data$ADT_maxID %>% table()

# Pre-filter 
Gina_seurat_obj@meta.data$ADT_classification.global %>% table()
Gina_seurat_obj@meta.data$ADT_maxID %>% table()

# # MULTIseqDemux
# Gina_seurat_obj_multi <- MULTIseqDemux(Gina_seurat_obj, assay = "ADT", quantile = 0.99)
# 
# Gina_seurat_obj_multi@meta.data$MULTI_classification %>% table()
# Gina_seurat_obj_multi@meta.data$MULTI_ID %>% table()
