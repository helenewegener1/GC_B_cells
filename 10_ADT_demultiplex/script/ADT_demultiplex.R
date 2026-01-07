getwd()

library(SeuratObject)
library(DropletUtils)
library(Seurat)
library(tidyverse)
library(glue)
library(patchwork)
library(readxl)
# library(ztable)
library(pheatmap)
library(ggsci)

# Load data
seurat_obj_nonDC_list <- readRDS("09_seurat_QC_clusters/out/seurat_obj_nonDC_list.rds")

# Get sample
sample_name <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
seurat_obj <- seurat_obj_nonDC_list[[sample_name]]

DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) + NoLegend()

# seurat_obj$scDblFinder.class %>% table()
# seurat_obj$DC_bool %>% table()

# N cells 
n_cells <- seurat_obj %>% ncol()
n_doublets <- seurat_obj$scDblFinder.class

# CLR normalization of ADT
seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")

############################ Run demultplexing tools ########################### 

# Demultiplexing with HTODemux and MULTIseqDemux
seurat_obj <- HTODemux(seurat_obj, assay = "ADT", positive.quantile = 0.999)
seurat_obj$ADT_classification.global %>% table()

seurat_obj <- MULTIseqDemux(seurat_obj, assay = "ADT", autoThresh = FALSE, quantile = 0.7) # Default
seurat_obj$MULTI_ID_0.7 <- seurat_obj$MULTI_ID
seurat_obj$MULTI_classification_0.7 <- seurat_obj$MULTI_classification
seurat_obj$MULTI_ID_0.7 %>% table()

seurat_obj <- MULTIseqDemux(seurat_obj, assay = "ADT", autoThresh = TRUE)
seurat_obj$MULTI_ID_autoThresh <- seurat_obj$MULTI_ID
seurat_obj$MULTI_classification_autoThresh <- seurat_obj$MULTI_classification
seurat_obj$MULTI_ID_autoThresh %>% table()

############################ Add comparable columns ############################

# Mutate columns so result from HTODemux and MULTIseqDemux become more comparable 
seurat_obj@meta.data$ADT_ID <- lapply(seurat_obj@meta.data$ADT_classification, function(x) {ifelse(str_detect(x, "_"), "Doublet", x)}) %>% unlist()
seurat_obj$MULTI_classification.global <- seurat_obj$MULTI_ID %>% str_replace("Fol-\\d+", "Singlet")

# Look at results
seurat_obj$ADT_ID %>% table()
seurat_obj$MULTI_ID %>% table()

seurat_obj$ADT_classification.global %>% table()
seurat_obj$MULTI_classification.global %>% table()

##################################### Plot ##################################### 

# Plot Heatmap
HTOHeatmap(seurat_obj, assay = "ADT", ncells = 5000)

# Plot density of CLR values for each ADT
# Get CLR-normalized data
adt_data <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")
dim(adt_data)

# Compute mean CLR value per hashtag
adt_means <- Matrix::rowMeans(adt_data)

# Summarize in a table
adt_summary <- data.frame(
  ADT = names(adt_means),
  mean_CLR = adt_means,
  median_CLR = apply(adt_data, 1, median),
  max_CLR = apply(adt_data, 1, max),
  q99_CLR = apply(adt_data, 1, function(x) quantile(x, 0.99))
) %>%
  arrange(desc(mean_CLR))

# Pull CLR ADT matrix
clr_counts <- seurat_obj[["ADT"]]$data
expr_counts <- seurat_obj[["ADT"]]$counts

rowSums(clr_counts)
rowSums(expr_counts)

expr <- as.data.frame(rowSums(expr_counts)) %>% rownames_to_column("ADT")
colnames(expr) <- c("ADT", "summed_expr_ADT")
expr

adt_summary <- adt_summary %>% left_join(expr)
rownames(adt_summary) <- adt_summary$ADT

adt_summary %>% arrange(ADT)

# Plot log count value range
seurat_obj[["ADT"]]$counts %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "counts") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}'),
         summed_expr_ADT = glue('{ADT}: {adt_summary[ADT,"summed_expr_ADT"]}')
         ) %>% 
  ggplot(aes(x = log(counts), fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(vars(summed_expr_ADT)) + 
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/ADT_log_counts_density.png"), width = 10, height = 7)

# log CLR range
seurat_obj[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  ggplot(aes(x = log(CLR), fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(vars(ADT_CLR_mean)) +
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/ADT_log_CLR_values_density.png"), width = 10, height = 7)

# CLR range
seurat_obj[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  ggplot(aes(x = CLR, fill = ADT_CLR_mean)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(vars(ADT_CLR_mean)) +
  theme_bw() + 
  theme(legend.position = "none") 

ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/ADT_CLR_values_density.png"), width = 10, height = 7)

# CLR value range
ADT_exprs <- seurat_obj[["ADT"]]$counts %>% rowSums() 
fols <- rownames(seurat_obj[["ADT"]]$data)

for (fol in fols){
  
  # fol <- "Fol-11"
  
  # ADT_exp <- format(ADT_exprs[fol], big.mark = ",")
  
  # log counts
  seurat_obj[["ADT"]]$counts %>% t() %>% as.data.frame() %>%
    pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "counts") %>% 
    filter(ADT == fol) %>%
    ggplot(aes(x = log(counts))) + 
    geom_density(alpha = 0.5) + 
    scale_x_continuous(
      breaks = seq(0, 10, by = 0.5),        # labeled ticks
      minor_breaks = seq(0, 10, by = 0.1)   # unlabeled splits
    ) +
  theme_bw() + 
    theme(legend.position = "none") + 
    labs(title = fol, x = "LogNormalzied counts")
  
  ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/LogNormalized_counts_per_fol/{fol}_ADT_values_density.png"), width = 14, height = 6)
  
  # # log CLR range
  # seurat_obj[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  #   pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  #   mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  #   filter(ADT == fol) %>% 
  #   ggplot(aes(x = CLR, fill = ADT_CLR_mean)) + 
  #   geom_density(alpha = 0.5) + 
  #   scale_x_continuous(trans = "log", breaks = scales::log_breaks(n = 40)) + 
  #   # scale_x_continuous(
  #   #   breaks = seq(-5, 5, by = 0.2),  # adjust 0.2 as needed
  #   #   sec.axis = sec_axis(~exp(.), name = "Original CLR")
  #   # ) +
  #   theme_bw() + 
  #   theme(legend.position = "none") +
  #   labs(title = fol,
  #        subtitle = glue("Summed ADT counts: {ADT_exp}"))
  # 
  # ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/CLR_values_per_fol/{fol}_ADT_log_CLR_values_density.png"), width = 16, height = 6)
  # 
  # # log CLR range
  # seurat_obj[["ADT"]]$data %>% t() %>% as.data.frame() %>% 
  #   pivot_longer(cols = starts_with("Fol"), names_to = "ADT", values_to = "CLR") %>% 
  #   mutate(ADT_CLR_mean = glue('{ADT} CLR_mean: {round(adt_summary[ADT,"mean_CLR"], 2)}')) %>% 
  #   filter(ADT == fol) %>% 
  #   ggplot(aes(x = CLR, fill = ADT_CLR_mean)) + 
  #   geom_density(alpha = 0.5) + 
  #   scale_x_continuous(trans = "log", breaks = scales::log_breaks(n = 40)) + 
  #   # scale_x_continuous(
  #   #   breaks = seq(-5, 5, by = 0.2),  # adjust 0.2 as needed
  #   #   sec.axis = sec_axis(~exp(.), name = "Original CLR")
  #   # ) +
  #   theme_bw() + 
  #   theme(legend.position = "none") +
  #   labs(title = fol,
  #        subtitle = glue("Summed ADT counts: {ADT_exp}"))
  # 
  # ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/CLR_values_per_fol/{fol}_ADT_log_CLR_values_density.png"), width = 16, height = 6)
  # 
}

############################## Manual threshold ################################ 

# three thresholds <- cut_close_to_noise, cut_in_middle, cut_close_to_signal

cut_close_to_noise_thresholds <- c(
  "Fol-1" = 1.75,
  "Fol-2" = 1.8,
  "Fol-3" = 1.5,
  "Fol-4" = 1.4,
  "Fol-5" = 1.3,
  "Fol-6" = 1.75,
  "Fol-7" = 1.7,
  "Fol-8" = 1.4,
  "Fol-9" = 2.25,
  "Fol-10" = 2,
  "Fol-11" = 1.4,
  "Fol-12" = 1.2,
  "Fol-13" = 1.2,
  "Fol-14" = 1.5,
  "Fol-15" = 1.75,
  "Fol-16" = 2.5,
  "Fol-17" = 1.5,  
  "Fol-18" = 1.25
)

cut_in_middle_thresholds <- c(
  "Fol-1" = 2,
  "Fol-2" = 2.1,
  "Fol-3" = 2,
  "Fol-4" = 1.74,
  "Fol-5" = 1.5,
  "Fol-6" = 2,
  "Fol-7" = 2.4,
  "Fol-8" = 1.75,
  "Fol-9" = 2.75,
  "Fol-10" = 2.75,
  "Fol-11" = 1.75,
  "Fol-12" = 1.4,
  "Fol-13" = 1.3,
  "Fol-14" = 2.5,
  "Fol-15" = 2.75,
  "Fol-16" = 3.50,
  "Fol-17" = 2,
  "Fol-18" = 1.5
)

cut_close_to_signal_thresholds <- c(
  "Fol-1" = 2.5,
  "Fol-2" = 2.2,
  "Fol-3" = 2.75,
  "Fol-4" = 2.25,
  "Fol-5" = 1.6,
  "Fol-6" = 2.3,
  "Fol-7" = 3.3,
  "Fol-8" = 2,
  "Fol-9" = 3.75,
  "Fol-10" = 3.25,
  "Fol-11" = 2,
  "Fol-12" = 1.6,
  "Fol-13" = 1.4,
  "Fol-14" = 3,
  "Fol-15" = 3.5,
  "Fol-16" = 4.5,
  "Fol-17" = 2.5,
  "Fol-18" = 1.75
)

# ---------- UNIVERSAL FUNCTION ----------
manual_call <- function(clr_matrix, thresholds) {
  
  # Boolean matrix: ADT > threshold ?
  pos_matrix <- sweep(clr_matrix, 1, thresholds, FUN = ">")
  
  # How many ADTs positive per cell
  pos_counts <- colSums(pos_matrix)
  
  # Initialize calls
  cell_id <- rep("Negative", ncol(clr_matrix))
  
  # For singlets: pick the ADT with highest CLR *among positives*
  is_singlet <- pos_counts == 1
  if (any(is_singlet)) {
    # Pick exactly the row where pos_matrix is TRUE
    singlet_ADT <- apply(pos_matrix[, is_singlet, drop = FALSE], 2,
                         function(x) rownames(pos_matrix)[which(x)])
    cell_id[is_singlet] <- singlet_ADT
  }
  
  # Doublets
  cell_id[pos_counts > 1] <- "Doublet"
  
  # Classification column
  classification <- ifelse(cell_id == "Negative", "Negative",
                           ifelse(cell_id == "Doublet", "Doublet", "Singlet"))
  
  return(list(id = cell_id, class = classification))
  
}


# ---------- APPLY TO ALL 3 THRESHOLD SETS ----------
# 1) close to noise
res_noise <- manual_call(clr_counts, cut_close_to_noise_thresholds)
seurat_obj$manual_ADT_ID_cut_close_to_noise <- res_noise$id
seurat_obj$manual_ADT_class_cut_close_to_noise <- res_noise$class

seurat_obj$manual_ADT_class_cut_close_to_noise %>% table()

# 2) middle
res_middle <- manual_call(clr_counts, cut_in_middle_thresholds)
seurat_obj$manual_ADT_ID_cut_in_middle <- res_middle$id
seurat_obj$manual_ADT_class_cut_in_middle <- res_middle$class

seurat_obj$manual_ADT_class_cut_in_middle %>% table()

# 3) close to signal
res_signal <- manual_call(clr_counts, cut_close_to_signal_thresholds)
seurat_obj$manual_ADT_ID_cut_close_to_signal <- res_signal$id
seurat_obj$manual_ADT_class_cut_close_to_signal <- res_signal$class

seurat_obj$manual_ADT_class_cut_close_to_signal %>% table()

# clr_counts <- seurat_obj[["ADT"]]$data
# 
# # Initialize empty matrix to store positive calls
# pos_matrix <- matrix(0, nrow = nrow(clr_counts), ncol = ncol(clr_counts),
#                      dimnames = dimnames(clr_counts))
# 
# # Loop through hashtags and apply threshold
# for (h in rownames(clr_counts)) {
#   pos_matrix[h, ] <- clr_counts[h, ] > cut_close_to_noise_thresholds[h]
# }
# 
# # Count positives per cell
# pos_counts <- colSums(pos_matrix)
# 
# # Assign singlet/doublet/negative
# cell_class <- rep("Negative", ncol(clr_counts))
# cell_class[pos_counts == 1] <- rownames(pos_matrix)[apply(pos_matrix, 2, which.max)][pos_counts == 1]
# cell_class[pos_counts > 1] <- "Doublet"
# 
# cell_class %>% table()
# 
# # Insert manual ADT demultplexing in seurat object
# seurat_obj$manual_ADT_ID_cut_close_to_noise_thresholds <- cell_class
# seurat_obj$manual_ADT_classification_cut_close_to_noise_thresholds <- seurat_obj$manual_ADT_ID %>% str_replace("Fol-\\d+", "Singlet")

################################# LARS METHOD ################################## 

# Normalize ranges
# seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
# seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "LogNormalize")
# seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "RC")

ADT_counts <- seurat_obj@assays$ADT$counts # Raw counts
ADT_CLR <- seurat_obj@assays$ADT$data # CLR normalized 

# Data wrangle
ADT_counts_t <- ADT_counts %>% as.matrix() %>% t() %>% 
  as.data.frame() %>% rownames_to_column("Cell") %>% 
  pivot_longer(cols = starts_with("Fol"),
               names_to = "Fol", 
               values_to = "Count")

# Show the idea
log_ratios <- ADT_counts_t %>% 
  filter(Cell == "AAACCAAAGACGAACG-1") %>%
  # group_by(Cell) %>% 
  arrange(desc(Count)) %>% 
  mutate(
    Fol = factor(Fol, levels = Fol),
    next_Count = lead(Count),
    log2_ratio = log2(next_Count / Count),
    label_y = pmin(Count, next_Count) * 1.05,  # slightly above lower bar
    label_x = as.numeric(Fol) + 0.5             # between bars
  )

# plot
log_ratios %>% 
  ggplot(aes(x = reorder(Fol, -Count), y = Count)) + 
  geom_col() +
  geom_label(
    data = log_ratios %>% filter(!is.na(log2_ratio)),
    aes(
      x = label_x,
      y = label_y,
      label = round(log2_ratio, 2)
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(x = "Fol", caption = "Log-ratio of 0: same value\nLog-ratio of -1: half the value") + 
  theme_bw()


############################### Compare to tools ###############################

seurat_obj$ADT_ID %>% table()
seurat_obj$MULTI_ID_0.7 %>% table()
seurat_obj$MULTI_ID_autoThresh %>% table()

df_compare_demux <- list(
  # DIY = cell_class %>% table(),
  manual_ADT_ID_cut_close_to_noise = seurat_obj$manual_ADT_ID_cut_close_to_noise %>% table(),
  manual_ADT_ID_cut_in_middle = seurat_obj$manual_ADT_ID_cut_in_middle %>% table(),
  manual_ADT_ID_cut_close_to_signal = seurat_obj$manual_ADT_ID_cut_close_to_signal %>% table(),
  HTO_0.999 = seurat_obj$ADT_ID %>% table(),
  MULTI_ID_0.7 = seurat_obj$MULTI_ID_0.7 %>% table(), 
  MULTI_ID_autoThresh = seurat_obj$MULTI_ID_autoThresh %>% table()
) %>%
  bind_rows(.id = "Method") %>%         # Method column
  as.data.frame()

df_compare_demux_long <- df_compare_demux %>% pivot_longer(cols = !Method,
                                                           names_to = "class", 
                                                           values_to = "count")

df_compare_demux <- df_compare_demux %>% mutate(N_drop_outs = Doublet + Negative)

text_box = glue(
  "
      N drop out
      Manual cut close to noise: {df_compare_demux[df_compare_demux$Method == 'manual_ADT_ID_cut_close_to_noise', 'N_drop_outs']}
      Manual cut close to middle: {df_compare_demux[df_compare_demux$Method == 'manual_ADT_ID_cut_in_middle', 'N_drop_outs']}
      Manual cut close to signal: {df_compare_demux[df_compare_demux$Method == 'manual_ADT_ID_cut_close_to_signal', 'N_drop_outs']}
      HTO_0.999: {df_compare_demux[df_compare_demux$Method == 'HTO_0.999', 'N_drop_outs']}
      MULTI_ID_0.7: {df_compare_demux[df_compare_demux$Method == 'MULTI_ID_0.7', 'N_drop_outs']}
      MULTI_ID_autoThresh: {df_compare_demux[df_compare_demux$Method == 'MULTI_ID_autoThresh', 'N_drop_outs']}
      "
)

# Bar plot comparing tools 
library(wesanderson)
df_compare_demux_long %>% 
  ggplot(aes(y = class, x = count, fill = Method)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  # scale_fill_manual(values = c("hotpink", "blue4", "green4", "red4"))
  scale_fill_manual(values = c(wes_palette("Moonrise3", n = 5), wes_palette("Moonrise2", n = 2)[2])) +
  annotate("text", x=2500, y=10, label= text_box)

ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/ADT_barplot_manual_VS_tool.png"), width = 10, height = 7)

# Bar plot only manual 
text_box = glue(
  "
      N drop out
      Manual cut close to noise: {df_compare_demux[df_compare_demux$Method == 'manual_ADT_ID_cut_close_to_noise', 'N_drop_outs']}
      Manual cut close to middle: {df_compare_demux[df_compare_demux$Method == 'manual_ADT_ID_cut_in_middle', 'N_drop_outs']}
      Manual cut close to signal: {df_compare_demux[df_compare_demux$Method == 'manual_ADT_ID_cut_close_to_signal', 'N_drop_outs']}
  "
)

df_compare_demux_long %>% 
  filter(startsWith(Method, "manual")) %>% 
  ggplot(aes(y = class, x = count, fill = Method)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  # scale_fill_manual(values = c("hotpink", "blue4", "green4", "red4"))
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5)[2:4]) +
  annotate("text", x=1500, y=10, label= text_box)

ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/ADT_barplot_manual.png"), width = 10, height = 7)

################################################################################ 
####################### Doublets: ADT VS DoubletFinder ######################### 
################################################################################

table(seurat_obj$scDblFinder.class) 
table(seurat_obj$manual_ADT_class_cut_close_to_signal) 
table(seurat_obj$scDblFinder.class, seurat_obj$manual_ADT_class_cut_close_to_signal) 

df <- as.data.frame(
  table(
    scDblFinder = seurat_obj$scDblFinder.class,
    ADT_class   = seurat_obj$manual_ADT_class_cut_close_to_signal
  )
)

ggplot(df, aes(x = ADT_class, y = scDblFinder, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "grey90", high = "steelblue") +
  labs(
    x = "Manual ADT classification",
    y = "scDblFinder classification",
    fill = "Cell count"
  ) +
  theme_minimal() + 
  theme(legend.position = "none")

ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/ADT_vs_scDblFinder_doublets.png"), width = 7, height = 4)

# New column to identify ADT singlets and scDblFinder doublets
# Potential communicating cells (!!!)
seurat_obj[[]] <- seurat_obj[[]] %>% 
  mutate(
    ADT_singlet_scDblFinder_doublets = ifelse(scDblFinder.class == "doublet" & manual_ADT_class_cut_close_to_signal == "Singlet", TRUE, FALSE)
  )

seurat_obj$ADT_singlet_scDblFinder_doublets %>% table()

DimPlot(seurat_obj, group.by = "ADT_singlet_scDblFinder_doublets", order = TRUE) + 
  labs(subtitle = sample_name, caption = "Potential communicating cells") 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_ADT_singlet_scDblFinder_doublets.png"), width = 8, height = 8)

################################################################################
################################### DimPlot #################################### 
################################################################################

seurat_obj$manual_ADT_ID_cut_close_to_noise

# Negatives seem to cluster together
DimPlot(seurat_obj, group.by = "manual_ADT_class_cut_close_to_signal") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_ADT_class.png"), width = 8, height = 8)

DimPlot(seurat_obj, group.by = "manual_ADT_ID_cut_close_to_noise") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_ADT_ID.png"), width = 8, height = 8)

FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA") + theme_classic()
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/FeatureScatter.png"), width = 8, height = 6)

seurat_obj@meta.data$percent.ribo %>% hist()

DimPlot(seurat_obj, group.by = "scDblFinder.class") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_scDblFinder.class.png"), width = 8, height = 8)

FeaturePlot(seurat_obj, features = "nFeature_RNA") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_nFeature_RNA.png"), width = 8, height = 8)

FeaturePlot(seurat_obj, features = "nCount_RNA") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_nCount_RNA.png"), width = 8, height = 8)

FeaturePlot(seurat_obj, features = "percent.mt") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_percent.mt.png"), width = 8, height = 8)

FeaturePlot(seurat_obj, features = "percent.ribo") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_percent.ribo.png"), width = 8, height = 8)

# seurat_obj$percent.ribo %>% range()

FeaturePlot(seurat_obj, features = "sce_contamination") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_sce_contamination.png"), width = 8, height = 8)

DimPlot(seurat_obj, group.by = "Phase") + 
  labs(subtitle = sample_name) 
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/UMAP_Phase.png"), width = 8, height = 8)

FeaturePlot(
  seurat_obj,
  features = c("MALAT1", "FOS"),
  reduction = "umap"
)

VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
        group.by = "manual_ADT_class_cut_close_to_noise", ncol = 2)
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/ADT_VlnPlot.png"), width = 12, height = 10)

# DimPlot(seurat_obj, group.by = "ADT_classification.global")
# DimPlot(seurat_obj, group.by = "seurat_clusters")
# seurat_obj$seurat_clusters

######################## Investigate individual cells ########################## 

# cell_id <- colnames(seurat_obj)[6] # Singlet 
# cell_id <- colnames(seurat_obj)[3] # HTO Doublet, MultiSeq singlet
# cell_id <- colnames(seurat_obj)[10] # HTO Singlet, MultiSeq negative
cell_id <- colnames(seurat_obj)[32] # Double doublet (Fol 18)
# cell_id <- colnames(seurat_obj)[1001] # Negative
# cell_id <- colnames(seurat_obj)[2005] # Fol 18 

cell_ADT <- seurat_obj@meta.data[cell_id, ]$ADT_classification
cell_MULTI_0.7 <- seurat_obj@meta.data[cell_id, ]$MULTI_classification_0.7
cell_MULTI_autoThresh <- seurat_obj@meta.data[cell_id, ]$MULTI_classification_autoThresh
cell_manual_ADT <- seurat_obj@meta.data[cell_id, ]$manual_ADT_ID
adt_summary <- adt_summary %>% rownames_to_column("Fol")

one_cell <- seurat_obj[["ADT"]]$data[,cell_id] %>% as.data.frame() %>% rownames_to_column("Fol")
colnames(one_cell) <- c("Fol", cell_id)

one_cell <- left_join(one_cell, adt_summary, by = "Fol")

one_cell %>% 
  ggplot(aes(x = Fol, y = !!sym(cell_id))) + 
  geom_col() + 
  geom_point(aes(x = Fol, y = max_CLR), color = "red2", alpha = 0.5) + 
  geom_point(aes(x = Fol, y = q99_CLR), color = "blue2", alpha = 0.5) + 
  theme_bw() + 
  labs(caption = glue("DIY: {cell_manual_ADT}\nHTODemux: {cell_ADT}\nMULTIseqDemux 0.7: {cell_MULTI_0.7}\nMULTIseqDemux autoThresh: {cell_MULTI_autoThresh}"))

###################### CONSENSUS OF HTODemux AND HTODemux ######################
# 
# Let's not do this
# Lars says that we really need to understand if it makes sense.
# Need to understand excatly how the tools work and a good reason why one might call more doublets and another more negatives. 
# 
# seurat_obj$ADT_ID
# 
# # Compare result of methods 
# table(seurat_obj$ADT_ID, seurat_obj$MULTI_ID)
# 
# table(seurat_obj$ADT_ID, seurat_obj$MULTI_ID) %>% 
#   heatmap(Colv = NA, Rowv = NA, xlab = "HTODemux", ylab = "MULTIseqDemux")
# ggsave(glue("10_ADT_demultiplex/plot/compare_methods_labels_heatmap.png"), width = 8, height = 7)
# 
# # Consensus data used to plot and extract consensus column for final seurat object 
# df_consens <- seurat_obj@meta.data %>% 
#   select(ADT_classification, ADT_classification.global, ADT_ID, ADT_maxID, ADT_secondID, MULTI_ID, MULTI_classification, MULTI_classification.global) %>% 
#   # Extract MULTI_maxID and MULTI_secondID to compare to HTODemux annotations 
#   mutate(
#     MULTI_maxID = str_split_i(MULTI_classification, "_", 1),
#     MULTI_secondID = str_split_i(MULTI_classification, "_", 2)
#   ) %>% 
#   mutate(
#     
#     # Everything same, also same exact doublets, else NA
#     ADT_consensus_very_hard = ifelse(ADT_classification == MULTI_classification, ADT_classification, NA),
#     
#     # Same singlet, negative and doublet detection, else NA
#     ADT_consensus_hard = ifelse(ADT_ID == MULTI_ID, ADT_ID, NA),
#     
#     # Since HTODemux detects many doublets and MULTIseqDemux detects more negatives. 
#     ADT_consensus_medium = case_when(
#       # if HTODemux and MULTIseqDemux agree - keep label
#       ADT_ID == MULTI_ID ~ ADT_ID,
#       
#       # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
#       ADT_classification.global == "Singlet" & MULTI_classification.global == "Singlet" & ADT_classification != MULTI_maxID ~ "Negative",
#       
#       # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet IF ADT_maxID == MULTI_ID
#       ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet" & ADT_maxID == MULTI_ID ~ MULTI_ID,
#       # else if ADT_maxID != MULTI_ID
#       ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet" & ADT_maxID != MULTI_ID ~ "Doublet",
#       
#       # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet IF ADT_classification == MULTI_maxID
#       # MULTI_maxID will be Negative if MULTI_ID is Negative 
#       # ADT_classification.global == "Singlet" & MULTI_ID == "Negative" & ADT_classification == MULTI_maxID ~ ADT_classification,
#       # else if ADT_maxID != MULTI_ID
#       # ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ "Negative",
#       ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ ADT_classification,
#       
#       # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
#       ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
#       
#       # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical but the second ones aren't.
#       ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID != MULTI_secondID ~ ADT_maxID,
#       # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet if the two first and the second follicles are identical.
#       ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID & ADT_secondID == MULTI_secondID ~ "Doublet",
#       # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
#       ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"
# 
#     ),
#     ADT_consensus_soft = case_when(
#       # if HTODemux and MULTIseqDemux agree - keep label
#       ADT_ID == MULTI_ID ~ ADT_ID,
#       
#       # if HTODemux detects a singlet and MULTIseqDemux detects a singlet - let it be a negative if the top follicles are NOT identical
#       ADT_classification.global == "Singlet" & MULTI_classification.global == "Singlet" & ADT_classification != MULTI_maxID ~ "Negative",
#       
#       # if HTODemux detects a doublet and MULTIseqDemux detects a singlet - let it be a singlet
#       ADT_classification.global == "Doublet" & MULTI_classification.global == "Singlet"  ~ MULTI_ID,
#       # if HTODemux detects a singlet and MULTIseqDemux detects a negative - let it be a singlet
#       ADT_classification.global == "Singlet" & MULTI_ID == "Negative" ~ ADT_classification,
#       
#       # if HTODemux detects a doublet and MULTIseqDemux detects a negative - let it be a negative
#       ADT_classification.global == "Doublet" & MULTI_ID == "Negative" ~ "Negative",
#       
#       # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a singlet if the two first follicles are identical.
#       ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID == MULTI_maxID  ~ ADT_maxID,
#       # if HTODemux detects a doublet and MULTIseqDemux detects a doublet - let it be a doublet
#       ADT_classification.global == "Doublet" & MULTI_ID == "Doublet" & ADT_maxID != MULTI_maxID ~ "Doublet"
#       
#     )
#   )
# 
# head(df_consens) 
# 
# # Tables
# table(df_consens$ADT_consensus_very_hard, useNA = "always") # same doublet 
# table(df_consens$ADT_consensus_hard, useNA = "always")
# table(df_consens$ADT_consensus_medium, useNA = "always")
# table(df_consens$ADT_consensus_soft, useNA = "always")
# 
# df_consens$ADT_classification.global %>% table()
# df_consens$MULTI_classification.global %>% table()
# 
# df_consens$ADT_consensus_hard.global <- df_consens$ADT_consensus_hard %>% str_replace("Fol-\\d+", "Singlet") 
# df_consens$ADT_consensus_medium.global <- df_consens$ADT_consensus_medium %>% str_replace("Fol-\\d+", "Singlet") 
# df_consens$ADT_consensus_soft.global <- df_consens$ADT_consensus_soft %>% str_replace("Fol-\\d+", "Singlet")
# 
# table(df_consens$ADT_consensus_hard.global, useNA = "always")
# table(df_consens$ADT_consensus_medium.global, useNA = "always")
# table(df_consens$ADT_consensus_soft.global, useNA = "always")
# 
# # When MULTIseqDemux detects Negative, what does HTODemux detect?
# df_consens %>% 
#   filter(MULTI_ID == "Negative") %>% 
#   ggplot(aes(x = ADT_ID)) + 
#   geom_bar() +
#   theme_bw() + 
#   labs(title = "MULTIseqDemux = Negative", 
#        subtitle = "When MULTIseqDemux detects Negative, what does HTODemux detect?")
# 
# ggsave("10_ADT_demultiplex/plot/MULTIseqDemux_Negative_barplot.png", width = 10, height = 7)
# 
# # When HTODemux detects Doublet, what is the max and second tag?
# df_consens %>% 
#   filter(ADT_classification.global == "Doublet") %>% 
#   ggplot(aes(x = ADT_maxID, 
#              color = ADT_secondID)) + 
#   geom_bar() +
#   theme_bw() + 
#   labs(title = "HTODemux = Doublet", 
#        subtitle = "When HTODemux detects Doublet, what is the max and second tag?")
# 
# ggsave("10_ADT_demultiplex/plot/HTODemux_Doublet.png", width = 10, height = 7)
# 
# #################### Add consensus colums to seurat object #####################
# 
# seurat_obj@meta.data %>% dim()
# df_consens %>% dim()
# 
# seurat_obj@meta.data$ADT_consensus_very_hard <- df_consens$ADT_consensus_very_hard
# seurat_obj@meta.data$ADT_consensus_hard <- df_consens$ADT_consensus_hard
# seurat_obj@meta.data$ADT_consensus_medium <- df_consens$ADT_consensus_medium
# seurat_obj@meta.data$ADT_consensus_soft<- df_consens$ADT_consensus_soft
# 
# seurat_obj@meta.data$ADT_consensus_hard.global<- df_consens$ADT_consensus_hard.global
# seurat_obj@meta.data$ADT_consensus_medium.global<- df_consens$ADT_consensus_medium.global
# seurat_obj@meta.data$ADT_consensus_soft.global <- df_consens$ADT_consensus_soft.global
# 
# ############################### ADT Ridge Plots ################################ 
# 
# idents <- c("ADT_maxID", "ADT_ID", "ADT_consensus_medium", "ADT_consensus_soft")
#   
# for (ident in idents){
#   
#   Idents(seurat_obj) <- ident
#   RidgePlot(seurat_obj, assay = "ADT", features = rownames(seurat_obj[["ADT"]])) 
#   ggsave(glue("10_ADT_demultiplex/plot/RidgePlot_{ident}.png"), width = 16, height = 20)
#   
# }

Idents(seurat_obj) <- "manual_ADT_ID_cut_close_to_signal"
RidgePlot(seurat_obj, assay = "ADT", features = rownames(seurat_obj[["ADT"]]))
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/demultiplexing/RidgePlot.png"), width = 16, height = 20)

# 
# ################################### DimPlot ####################################
# 
# DimPlot(seurat_obj, group.by = "ADT_classification.global")
# DimPlot(seurat_obj, group.by = "MULTI_classification.global")
# 
# # Consensus
# DimPlot(seurat_obj, group.by = "ADT_consensus_hard.global")
# DimPlot(seurat_obj, group.by = "ADT_consensus_medium.global")
# DimPlot(seurat_obj, group.by = "ADT_consensus_soft.global")
# 
# seurat_obj$ADT_consensus_soft

################################### tSNE ####################################

Idents(seurat_obj) <- "manual_ADT_class_cut_close_to_signal"
# seurat_obj$manual_ADT_class_cut_close_to_signal

# First, we will remove negative cells from the object
seurat_obj.subset <- subset(seurat_obj, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(seurat_obj.subset) <- "ADT"
seurat_obj.subset <- ScaleData(seurat_obj.subset, features = rownames(seurat_obj.subset), verbose = FALSE)
seurat_obj.subset <- RunPCA(seurat_obj.subset, features = rownames(seurat_obj.subset), approx = FALSE)
seurat_obj.subset <- RunUMAP(seurat_obj.subset, dims = 1:8, perplexity = 100)
seurat_obj.subset <- RunTSNE(seurat_obj.subset, dims = 1:8, perplexity = 100)

DimPlot(seurat_obj.subset, reduction = "umap") + labs(title = "Manual ADT classification")
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/DimPlot_Doublet_Singlet_umap.png"), width = 9, height = 7)

DimPlot(seurat_obj.subset, reduction = "tsne") + labs(title = "Manual ADT classification")
ggsave(glue("10_ADT_demultiplex/plot/{sample_name}/DimPlot_Doublet_Singlet_tsne.png"), width = 9, height = 7)

####################### Export ADT demultiplexed objects ####################### 

# When doing this for multiple samples, make a list with objects. 

saveRDS(seurat_obj, "10_ADT_demultiplex/out/seurat_obj_ADT_demultiplexed.rds")





