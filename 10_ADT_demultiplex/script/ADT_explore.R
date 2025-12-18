getwd()

library(SeuratObject)
library(DropletUtils)
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

seurat_obj_list_QC <- readRDS("07_seurat_QC/out/seurat_obj_QC.rds")
seurat_obj_QC <- seurat_obj_list_QC[["HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"]]

# Explore data
dim(seurat_obj@assays$RNA)
dim(Gina_seurat_obj@assays$RNA)

# DimPlot(Gina_seurat_obj, group.by = )

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

df_meta_data <- Gina_seurat_obj@meta.data %>% select(ADT_classification, ADT_maxID, ADT_classification.global, Celltype) %>% rownames_to_column("cell")

df <- df %>% left_join(df_meta_data, by = "cell", relationship = "many-to-many")

head(df)

cell_doublets <- Gina_seurat_obj@meta.data %>% 
  filter(ADT_classification.global == "Doublet" & Celltype != "DCs") %>% 
  rownames()

cell_subset_list <- list(
  x1 = cell_doublets[1:15],
  x2 = cell_doublets[5:10], 
  x3 = cell_doublets[10:15]
)

x <- "x1"

df %>% 
  filter(cell %in% cell_subset_list[[x]]) %>%
  filter(ADT_classification.global == "Doublet") %>%
  ggplot(aes(x = follicles, 
             y = value,
             color = cell)) + 
  # geom_col(position = "dodge") + 
  geom_point() +
  geom_line(aes(group = cell)) +
  facet_wrap(vars(ADT_classification)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
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

# Doublets investigate
df %>% arrange(cell) %>% filter(ADT_classification.global == "Doublet") %>% head(n=20)

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

############## Filter out DCs ############## 

# Add scDblFinder info
scDblFinder_info <- seurat_obj_QC$scDblFinder.class[colnames(Gina_seurat_obj)]
Gina_seurat_obj <- AddMetaData(Gina_seurat_obj, metadata = scDblFinder_info, col.name = "scDblFinder.class")

# Filter out DCs
Gina_seurat_obj_DCs <- subset(Gina_seurat_obj, subset = Celltype == "DCs")
Gina_seurat_obj_DCs$Celltype %>% table()

Gina_seurat_obj_other <- subset(Gina_seurat_obj, subset = Celltype != "DCs")
Gina_seurat_obj_other$Celltype %>% table()

# Demultiplex
Gina_seurat_obj_other <- NormalizeData(Gina_seurat_obj_other, assay = "ADT", normalization.method = "CLR")
Gina_seurat_obj_other <- HTODemux(Gina_seurat_obj_other, assay = "ADT", positive.quantile = 0.9999)

Gina_seurat_obj$ADT_classification.global %>% table()

# When removing the 617 DCs we get 529 fewer negatives. 
Gina_seurat_obj_other$ADT_classification.global %>% table()

# MULTIseqDemux gives many negatives 
Gina_seurat_obj_other <- MULTIseqDemux(Gina_seurat_obj_other, assay = "ADT", quantile = 0.8)
Gina_seurat_obj_other$MULTI_ID %>% table()

# Compare methods 
table(Gina_seurat_obj_other$ADT_classification.global, Gina_seurat_obj_other$MULTI_ID)

# Doublets of demultiplexing VS doublets of scDblFinder
table(Gina_seurat_obj_other$ADT_classification.global, Gina_seurat_obj_other$scDblFinder.class)
table(Gina_seurat_obj_other$MULTI_ID, Gina_seurat_obj_other$scDblFinder.class)

# PLOT DIFFERENCE 

df_other <- Gina_seurat_obj_other[["ADT"]]$counts %>% 
  as.data.frame() %>% 
  rownames_to_column("follicles") %>% 
  pivot_longer(cols = colnames(Gina_seurat_obj_other), names_to = "cell", values_to = "value")

df_meta_data_other <- Gina_seurat_obj_other@meta.data %>% 
  select(MULTI_ID, MULTI_classification, ADT_maxID, ADT_classification.global, ADT_classification, Celltype) %>% rownames_to_column("cell")

df_other <- df_other %>% left_join(df_meta_data_other, by = "cell", relationship = "many-to-many")
head(df_other)

df_other %>% 
  # filter(cell %in% cell_subset_list[[x]]) %>%
  filter(ADT_classification.global == "Doublet") %>%
  arrange(cell) %>% 
  head(n = 10*18) %>% 
  ggplot(aes(x = follicles, 
             y = value,
             color = ADT_classification)) + 
  # geom_col(position = "dodge") + 
  geom_point() +
  geom_line(aes(group = cell)) +
  facet_wrap(vars(ADT_maxID)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  labs(
    x = "", 
    title = "5 cells where ADT_classification.global = 'Doublet'"
  )

ggsave(glue("07_seurat_QC/plot/ADT_explore/MULTI_ID_Negative.png"), width = 12, height = 8)

df_other %>% 
  # filter(cell %in% cell_subset_list[[x]]) %>%
  filter(MULTI_ID == "Negative") %>%
  arrange(cell) %>% 
  head(n = 6*18) %>% 
  ggplot(aes(x = follicles, 
             y = value,
             color = ADT_classification)) + 
  # geom_col(position = "dodge") + 
  geom_point() +
  geom_line(aes(group = cell)) +
  facet_wrap(vars(MULTI_classification)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  labs(
    x = "", 
    title = "5 cells where MULTI_ID = 'Negative'"
  )

ggsave(glue("07_seurat_QC/plot/ADT_explore/MULTI_ID_Negative.png"), width = 12, height = 8)

# Summary of mean CLR-normalized values per hashtag

# Get CLR-normalized data
adt_data <- GetAssayData(Gina_seurat_obj_other, assay = "ADT", slot = "data")
dim(adt_data)

# Compute mean CLR value per hashtag
adt_means <- Matrix::rowMeans(adt_data)

# Summarize in a table
adt_summary <- data.frame(
  ADT = names(adt_means),
  mean_CLR = adt_means,
  median_CLR = apply(adt_data, 1, median),
  max_CLR = apply(adt_data, 1, max)
) %>%
  arrange(desc(mean_CLR))

adt_summary

# Weak tags are the ones with low CLR mean 
weak_tags <- adt_summary %>% arrange(mean_CLR) %>% head(6) %>% pull(ADT)
RidgePlot(Gina_seurat_obj_other, assay = "ADT", features = weak_tags, ncol = 3)

adt_summary_arrange <- adt_summary %>% arrange(ADT)

#### ADT range
Gina_seurat_obj_other[["ADT"]]$counts %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "counts") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  ggplot(aes(x = log(counts), fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(vars(ADT_CLR_mean)) + 
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("07_seurat_QC/plot/ADT_explore/fol_density_counts.png"), width = 10, height = 7)

adt_long <- Gina_seurat_obj_other_singlets[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  rownames_to_column("cell") %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}'))

# Compute 99th quantile per ADT
quantiles <- adt_long %>%
  group_by(ADT) %>%
  summarize(q99 = quantile(CLR, 0.99))

n_above_quantile <- adt_long %>% left_join(quantiles, by = "ADT") %>% ungroup()
n_above_quantile <- n_above_quantile %>% mutate(is_q99 = ifelse(CLR >= q99, TRUE, FALSE))

n_above_quantile %>% group_by(ADT) %>% dplyr::count(is_q99) %>% filter(is_q99 == TRUE)

n_above_quantile %>% filter(is_q99) %>% select(cell)

n_q99 <- (ncol(Gina_seurat_obj_other_singlets) / 100) * (100 - 99) 
n_q99 <- ceiling(n_q99)

n_above_quantile %>% 
  ggplot(aes(x = CLR, fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  geom_vline(aes(xintercept = q99), linetype = "dashed") +
  facet_wrap(vars(ADT_CLR_mean)) +
  theme_bw() + 
  theme(legend.position = "none") + 
  labs(subtitle = glue("q99, n = {n_q99}"))

ggsave(glue("07_seurat_QC/plot/ADT_explore/fol_density_CLR.png"), width = 10, height = 7)


df_other <- Gina_seurat_obj_other[["ADT"]]$data %>% 
  as.data.frame() %>% 
  rownames_to_column("follicles") %>% 
  pivot_longer(cols = colnames(Gina_seurat_obj_other), names_to = "cell", values_to = "CLR")

df_meta_data_other <- Gina_seurat_obj_other@meta.data %>% 
  select(ADT_classification, ADT_maxID, ADT_classification.global, Celltype, scDblFinder.class) %>% 
  rownames_to_column("cell")

df_other <- df_other %>% left_join(df_meta_data_other, by = "cell", relationship = "many-to-many")

df_other %>% colnames()
head(df_other)

cell_doublets <- Gina_seurat_obj_other@meta.data %>% 
  filter(ADT_classification.global == "Doublet" & Celltype != "DCs" & scDblFinder.class == "singlet") %>% 
  rownames()


df_other %>% 
  arrange(cell) %>% 
  filter(ADT_classification.global == "Doublet") %>%
  head(10*18) %>%
  ggplot(aes(x = follicles, 
             y = CLR,
             color = cell)) + 
  geom_point() +
  geom_line(aes(group = cell)) +
  facet_wrap(vars(ADT_classification)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  labs(
    x = "", 
    title = "Cells where ADT_classification.global = 'Doublet'"
  )

ggsave(glue("07_seurat_QC/plot/ADT_explore/doublets_lineplot_clr.png"), width = 12, height = 8)

# log 
# df_other %>% 
#   arrange(cell) %>% 
#   filter(ADT_classification.global == "Doublet") %>%
#   head(10*18) %>%
#   ggplot(aes(x = follicles, 
#              y = log(value),
#              color = cell)) + 
#   geom_point() +
#   geom_line(aes(group = cell)) +
#   facet_wrap(vars(ADT_classification)) + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         legend.position = "none")+ 
#   labs(
#     x = "", 
#     title = "Log Scale: Cells where ADT_classification.global = 'Doublet'"
#   )
# 
# ggsave(glue("07_seurat_QC/plot/ADT_explore/doublets_lineplot.png"), width = 12, height = 8)


# singlet 
df_other %>% 
  arrange(cell) %>% 
  filter(ADT_classification.global == "Singlet") %>%
  head(10*18) %>%
  ggplot(aes(x = follicles, 
             y = CLR,
             color = cell)) + 
  geom_point() +
  geom_line(aes(group = cell)) +
  facet_wrap(vars(ADT_classification)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  labs(
    x = "", 
    title = "Cells where ADT_classification.global = 'Singlet'"
  )

ggsave(glue("07_seurat_QC/plot/ADT_explore/singlet_lineplot_clr.png"), width = 12, height = 8)


df_other %>% 
  arrange(cell) %>% 
  filter(ADT_classification.global == "Negative") %>%
  head(5*18) %>%
  ggplot(aes(x = follicles, 
             y = CLR,
             color = cell)) + 
  geom_point() +
  geom_line(aes(group = cell)) +
  facet_wrap(vars(ADT_classification)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  labs(
    x = "", 
    title = "Cells where ADT_classification.global = 'Negative'"
  )

ggsave(glue("07_seurat_QC/plot/ADT_explore/negatives_lineplot_clr.png"), width = 12, height = 8)

# test HTODemux thresholds 
Gina_seurat_obj_other_singlets <- subset(Gina_seurat_obj_other, scDblFinder.class == 'singlet')

Gina_seurat_obj_other_singlets <- NormalizeData(Gina_seurat_obj_other_singlets, assay = "ADT", normalization.method = "CLR")
Gina_seurat_obj_other_singlets <- HTODemux(Gina_seurat_obj_other_singlets, assay = "ADT", positive.quantile = 0.9999)

Gina_seurat_obj_other_singlets$ADT_classification.global %>% table()

FeatureScatter(Gina_seurat_obj_other_singlets, feature1 = "Fol-1", feature2 = "Fol-6") + labs(title = "Two better ADTs")
ggsave("07_seurat_QC/plot/ADT_explore/FeatureScatter_good.png", width = 10, height = 8)
FeatureScatter(Gina_seurat_obj_other_singlets, feature1 = "Fol-1", feature2 = "Fol-15") + labs(title = "One good, one bad")
ggsave("07_seurat_QC/plot/ADT_explore/FeatureScatter_good_bad.png", width = 10, height = 8)
FeatureScatter(Gina_seurat_obj_other_singlets, feature1 = "Fol-9", feature2 = "Fol-15") + labs(title = "Two bad ADTs")
ggsave("07_seurat_QC/plot/ADT_explore/FeatureScatter_bad.png", width = 10, height = 8)

# demuxmix
dmm2 <- demuxmix(
  Gina_seurat_obj_other_singlets[["ADT"]]$counts %>% as.matrix, 
  rna = Gina_seurat_obj_other_singlets$nFeature_RNA, 
  model = "auto" # default
)

summary(dmm2)

classes <- dmmClassify(dmm1)
table(classes$HTO)

plotDmmHistogram(dmm1)


# SUBSETTING
Gina_seurat_obj_other_filter <- subset(
  Gina_seurat_obj_other, subset = nCount_ADT > 1000 & nCount_ADT < quantile(Gina_seurat_obj_other$nCount_ADT, 0.95)
)

Gina_seurat_obj_other_filter[["ADT"]]$counts %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "counts") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  ggplot(aes(x = log(counts), fill = ADT_CLR_mean)) + 
  geom_histogram(alpha = 0.5) + 
  facet_wrap(vars(ADT_CLR_mean)) + 
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("07_seurat_QC/plot/ADT_explore/fol_density_log_counts.png"), width = 10, height = 7)

Gina_seurat_obj_other_filter[["ADT"]]$counts %>% rowSums() %>% log()
Gina_seurat_obj_other_filter[["ADT"]]$data %>% rowMeans() 

# Demultiplex filtered
Gina_seurat_obj_other_filter <- NormalizeData(Gina_seurat_obj_other_filter, assay = "ADT", normalization.method = "CLR")
Gina_seurat_obj_other_filter <- HTODemux(Gina_seurat_obj_other_filter, assay = "ADT", positive.quantile = 0.99999)

Gina_seurat_obj_other_filter$ADT_classification.global %>% table()

Gina_seurat_obj_other$ADT_classification.global %>% table()

# CLR manually
x_counts <- Gina_seurat_obj_other_filter[["ADT"]]$counts[, "AAACCAGCAGCCCATA-1"] 
geometric_mean <- prod(x_counts + 1)^(1/length(x_counts))

log2(x_counts + 1) - mean(log2(x_counts + 1))

log((x_counts + 1)/geometric_mean, base = 2)

Gina_seurat_obj_other_filter[["ADT"]]$data[, "AAACCAGCAGCCCATA-1"] 

Gina_seurat_obj_other_filter <- NormalizeData(Gina_seurat_obj_other_filter, assay = "ADT", normalization.method = "CLR")

# CLR ADT summary
# Extract CLR data
clr_data <- GetAssayData(Gina_seurat_obj_other, assay = "ADT", slot = "data")

# Mean, median, max CLR per ADT
summary_stats <- data.frame(
  ADT = rownames(clr_data),
  mean_CLR = rowMeans(clr_data),
  median_CLR = apply(clr_data, 1, median),
  max_CLR = apply(clr_data, 1, max)
)

# Optionally, fraction of cells above a threshold
threshold <- 1  # example
summary_stats$frac_positive <- rowMeans(clr_data > threshold)
summary_stats

fols <- Gina_seurat_obj_other[["ADT"]] %>% rownames()

n <- length(fols)

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    fol1 <- fols[i]
    fol2 <- fols[j]
    
    p <- FeatureScatter(Gina_seurat_obj_other, feature1 = fol1, feature2 = fol2, slot = "data") # data == CLR transformed
    
    # Save plot
    ggsave(
      filename = glue("07_seurat_QC/plot/ADT_explore/FeatureScatter_all/FeatureScatter_{fol1}_{fol2}.png"),
      plot = p,
      width = 10,
      height = 8
    )
  }
}


# Different methods # Different methods rowMeans()
library(DropletUtils) # hashedDrops
library(Seurat) # HTODemux
library(demuxmix) # demuxmix
library(cellhashR) # BFF


# demuxmix
# rna not specifiec 
dmm1 <- demuxmix(
  Gina_seurat_obj_other[["ADT"]]$counts %>% as.matrix, 
  model = "naive" # when no rna
)

summary(dmm1)

classes <- dmmClassify(dmm1)
table(classes$HTO)

plotDmmHistogram(dmm1)

# rna specifiec 
dmm2 <- demuxmix(
  Gina_seurat_obj_other[["ADT"]]$counts %>% as.matrix, 
  rna = Gina_seurat_obj_other$nFeature_RNA, 
  model = "auto" # default
)

summary(dmm2)

# # MULTIseqDemux a lot of negatives 
Gina_seurat_obj_other <- MULTIseqDemux(Gina_seurat_obj_other, assay = "ADT", quantile = 0.80)

Gina_seurat_obj_other@meta.data$MULTI_classification %>% table()
Gina_seurat_obj_other@meta.data$MULTI_ID %>% table()

# Exclude 

Gina_seurat_obj_other_singlets[["ADT"]]$data %>% colSums() %>% range()
Gina_seurat_obj_other_singlets[["ADT"]]$data %>% rowSums() %>% range()


HTOHeatmap(Gina_seurat_obj_other_singlets, assay = "ADT")


# # MULTIseqDemux a lot of negatives 
Gina_seurat_obj_other_singlets <- MULTIseqDemux(Gina_seurat_obj_other_singlets, assay = "ADT",  autoThresh = TRUE)

Gina_seurat_obj_other_singlets@meta.data$MULTI_classification %>% table()
Gina_seurat_obj_other_singlets@meta.data$MULTI_ID %>% table()


####### CONSENSUS OF HTODemux AND HTODemux ####### 

# This object contains singlets (doublets removed with scDoubletfinder and only hashtagged cells (all other than DCs))
Gina_seurat_obj_other_singlets

# HTODemux
Gina_seurat_obj_other_singlets <- NormalizeData(Gina_seurat_obj_other_singlets, assay = "ADT", normalization.method = "CLR")
Gina_seurat_obj_other_singlets <- HTODemux(Gina_seurat_obj_other_singlets, assay = "ADT", positive.quantile = 0.99, kfunc = "clara")

Gina_seurat_obj_other_singlets@meta.data$ADT_classification.global %>% table()

# HTODemux
Gina_seurat_obj_other_singlets <- MULTIseqDemux(Gina_seurat_obj_other_singlets, assay = "ADT", autoThresh = TRUE)
# Gina_seurat_obj_other_singlets <- MULTIseqDemux(Gina_seurat_obj_other_singlets, assay = "ADT",  quantile = 0.3)

Gina_seurat_obj_other_singlets@meta.data$MULTI_ID %>% table()

Gina_seurat_obj_other_singlets@meta.data$MULTI_ID %>% str_replace("Fol-\\d+", "Singlet") %>% table()

# Prep for compare 
Gina_seurat_obj_other_singlets@meta.data$ADT_classification.global %>% table()

Gina_seurat_obj_other_singlets@meta.data$ADT_classification_compare <- lapply(Gina_seurat_obj_other_singlets@meta.data$ADT_classification, function(x) {ifelse(str_detect(x, "_"), "Doublet", x)}) %>% unlist()

Gina_seurat_obj_other_singlets@meta.data$ADT_classification_compare %>% table()

# Compare
table(Gina_seurat_obj_other_singlets@meta.data$ADT_classification_compare, Gina_seurat_obj_other_singlets@meta.data$MULTI_ID)

table(Gina_seurat_obj_other_singlets@meta.data$ADT_classification_compare, Gina_seurat_obj_other_singlets@meta.data$MULTI_ID) %>% 
  heatmap(Colv = NA, Rowv = NA, xlab = "HTODemux", ylab = "MULTIseqDemux")

# consensus barplot 
df_consens <- Gina_seurat_obj_other_singlets@meta.data %>% 
  select(ADT_classification, ADT_classification.global, ADT_classification_compare, ADT_maxID, ADT_secondID, MULTI_ID, MULTI_classification) %>% 
  mutate(
    MULTI_maxID = str_split_i(MULTI_classification, "_", 1),
    MULTI_secondID = str_split_i(MULTI_classification, "_", 2)
  ) %>% 
  mutate(
    ADT_consensus_very_hard = ifelse(ADT_classification == MULTI_classification, ADT_classification, NA),
    ADT_consensus_hard = ifelse(ADT_classification_compare == MULTI_ID, ADT_classification_compare, NA),
    # Since HTODemux detects many doublets and MULTIseqDemux detects more negatives. 
    ADT_consensus_medium = case_when(
      # if HTODemux and MULTIseqDemux agree - keep label
      ADT_classification_compare == MULTI_ID ~ ADT_classification_compare,
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
      ADT_classification.global == "Singlet" & !(MULTI_ID %in% c("Negative", "Doublet")) & ADT_classification != MULTI_maxID ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet IF ADT_maxID == MULTI_ID
      ADT_classification.global == "Doublet" & !(MULTI_ID %in% c("Negative", "Doublet")) & ADT_maxID == MULTI_ID ~ MULTI_ID,
      # else...
      ADT_classification.global == "Doublet" & !(MULTI_ID %in% c("Negative", "Doublet")) & ADT_maxID != MULTI_ID ~ "Doublet",
      # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet IF ADT_classification == MULTI_maxID
      ADT_classification.global == "Singlet" & MULTI_ID == "Negative" & ADT_classification == MULTI_maxID ~ ADT_classification,
      # else...
      ADT_classification.global == "Singlet" & MULTI_ID == "Negative" & ADT_classification != MULTI_maxID ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
      ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical but the second ones aren't.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID != MULTI_secondID ~ ADT_maxID,
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet if the two first and the second follicles are identical.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID == MULTI_secondID ~ "Doublet",
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"

    ),
    ADT_consensus_soft = case_when(
      # if HTODemux and MULTIseqDemux agree - keep label
      ADT_classification_compare == MULTI_ID ~ ADT_classification_compare,
      
      # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
      ADT_classification.global == "Singlet" & !(MULTI_ID %in% c("Negative", "Doublet")) & ADT_classification != MULTI_maxID ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet
      ADT_classification.global == "Doublet" & !(MULTI_ID %in% c("Negative", "Doublet"))  ~ MULTI_ID,
      # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet
      ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ ADT_classification,
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
      ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
      
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical.
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID  ~ ADT_maxID,
      # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
      ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"
      
    )
  )

head(df_consens) 

table(df_consens$ADT_consensus_very_hard, useNA = "always") # same doublet 
table(df_consens$ADT_consensus_hard, useNA = "always")
table(df_consens$ADT_consensus_medium, useNA = "always")
table(df_consens$ADT_consensus_soft, useNA = "always")

df_consens %>% filter(!is.na(ADT_consensus_soft)) %>% nrow()
df_consens %>% filter(!is.na(ADT_consensus_medium)) %>% nrow()

df_consens$ADT_consensus_hard %>% str_replace("Fol-\\d+", "Singlet") %>% table(useNA = "ifany")
df_consens$ADT_consensus_medium %>% str_replace("Fol-\\d+", "Singlet") %>% table(useNA = "ifany")
df_consens$ADT_consensus_soft %>% str_replace("Fol-\\d+", "Singlet") %>% table(useNA = "ifany")

df_consens %>% 
  filter(MULTI_ID == "Negative") %>% 
  ggplot(aes(x = ADT_classification_compare)) + 
  geom_bar() +
  theme_bw() + 
  labs(title = "MULTIseqDemux = Negative")

ggsave(glue("07_seurat_QC/plot/ADT_explore/consensus_MULTIseqDemux_Negative.png"), width = 10, height = 7)


df_consens %>% 
  filter(ADT_classification.global == "Doublet") %>% 
  ggplot(aes(x = MULTI_ID)) + 
  geom_bar() +
  theme_bw() + 
  labs(title = "HTODemux = Doublet")

ggsave(glue("07_seurat_QC/plot/ADT_explore/consensus_HTODemux_Doublet.png"), width = 10, height = 7)

df_consens %>% 
  filter(ADT_classification.global == "Doublet") %>% 
  ggplot(aes(x = ADT_maxID, fill = ADT_secondID)) + 
  geom_bar() +
  theme_bw() + 
  labs(title = "HTODemux = Doublet")

ggsave(glue("07_seurat_QC/plot/ADT_explore/consensus_HTODemux_Doublet.png"), width = 10, height = 7)


# pheatmap(dplyr::select(df_consens, ADT_classification, MULTI_classification), cluster_rows = FALSE, cluster_cols = FALSE)

# table(df_consens$ADT_classification, df_consens$MULTI_classification) %>% heatmap(Colv = NA, Rowv = NA)




############################# Sanity check for ADT - pre QC #############################

for (sample_name in names(seurat_obj_list)){
  
  seurat_obj <- seurat_obj_list[[sample_name]]
  
  if ("ADT" %in% names(seurat_obj)){
    
    Idents(seurat_obj)
    Idents(seurat_obj) <- "ADT_maxID"
    RidgePlot(seurat_obj, assay = "ADT", features = rownames(seurat_obj[["ADT"]]))
    ggsave(glue("07_seurat_QC/plot/ADT_RidgePlot/RidgePlot_{sample_name}_preQC.pdf"), width = 16, height = 20)
    
  }
  
}
