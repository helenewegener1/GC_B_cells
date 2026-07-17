library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC_3_definitions", rds_files, value = TRUE)
patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
patients

resolve_LC_list <- lapply(resolve_LC_files, function(x) readRDS(glue("45_immcantation/out/rds/{x}"))) %>% 
  setNames(patients)

# Look at clone IDs
grep("clone", colnames(resolve_LC_list$HH117), value = TRUE)

# clone_subgroup_id_90_similarity

# resolve_LC <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds"))
# table(resolve_LC$locus)
# 
# df_heavy <- resolve_LC %>% filter(locus == "IGH")
# 
# nrow(df_heavy)

# Load seurat object
seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

outdir <- glue("45_immcantation/plot/18_90_similarity/05_general_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

resolve_LC_list <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  resolve_LC_list[[HH]] %>%
    group_by(clone_subgroup_id_90_similarity, locus) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)

# ==============================================================================
# Export BCR meta data to Gina 
# ==============================================================================

# meta_4_Gina_list <- lapply(patients, function(HH){
# 
#   # HH <- "HH119"
# 
#   seurat_obj <- subset(seurat_integrated, patient == HH)
#   resolve_LC_HH <- resolve_LC_list[[HH]] %>% filter(locus == "IGH")
# 
#   # Check IDs
#   seurat_obj %>% colnames() %>% head()
#   resolve_LC_HH$cell_id %>% head()
# 
#   seurat_obj %>% colnames() %>% length()
#   resolve_LC_HH$cell_id %>% length()
# 
#   (seurat_obj %>% colnames() %>% length()) == (seurat_obj %>% colnames() %>% unique() %>% length())
#   (resolve_LC_HH$cell_id %>% length()) == (resolve_LC_HH$cell_id %>% unique() %>% length())
# 
# 
#   # Wrangle IDs
#   seurat_ids <- seurat_obj %>% colnames()
#   LC_ids <- resolve_LC_HH$cell_id_seurat %>% str_remove(".*?_")
# 
#   # # IDs test
#   # seurat_ids_sub <- seurat_ids %>% str_split_i("_", 2)
#   # table(seurat_ids_sub, seurat_obj$sample_clean)
#   #
#   # LC_ids_sub <- LC_ids %>% str_split_i("_", 2)
#   # table(LC_ids_sub, resolve_LC_HH$sample_clean)
#   # # End
# 
#   (seurat_ids %>% length()) == (seurat_ids %>% unique() %>% length())
#   (LC_ids %>% length()) == (LC_ids %>% unique() %>% length())
# 
#   table(LC_ids %in% seurat_ids)
# 
#   # Prep for merge
#   resolve_LC_HH_meta <- resolve_LC_HH %>%
#     mutate(cell_id_seurat_clean = str_remove(cell_id_seurat, ".*?_")) %>%
#     select(cell_id_seurat_clean, c_call, clone_subgroup_id_90_similarity)
# 
#   # Merge and create final meta data for Gina
#   meta_4_Gina <- seurat_obj[[]] %>%
#     select(manual_ADT_class, manual_ADT_ID, manual_ADT_full_ID, sample, L1_annotation) %>%
#     rownames_to_column("cell_id_seurat_clean") %>%
#     left_join(resolve_LC_HH_meta, by = "cell_id_seurat_clean") %>%
#     column_to_rownames("cell_id_seurat_clean")
# 
# 
#   # table(meta_4_Gina$L1_annotation, meta_4_Gina$c_call, useNA = "always")
# 
#   # Check
#   meta_4_Gina %>% count(clone_subgroup_id_90_similarity, sort = TRUE)  %>% head()
# 
#   return(meta_4_Gina)
# 
# }) %>% setNames(patients)
# 
# 
# saveRDS(meta_4_Gina_list %>% bind_rows(), "45_immcantation/out/rds/meta_4_Gina_list_90_similarity.rds")

# meta <- readRDS("45_immcantation/out/rds/meta_4_Gina_list.rds")
#
# meta %>% filter(patient_id == "HH117")
#
# resolve_LC_list$HH117$L1_annotation %>% table()


# ==============================================================================
# Percentage of data with BCR data of total data 
# ==============================================================================

outdir1 <- glue("{outdir}/bcr_available/")
dir.create(outdir1, recursive = TRUE, showWarnings = FALSE)

# meta <- readRDS("45_immcantation/out/rds/meta_4_Gina_list_90_similarity.rds")

B_cell_subsets <- list(
  "all B cells" = c("GC_Bcells", "Memory_Bcells", "Naive_Bcells", "PCs", "Unconventional_Bcells"), 
  "GC B cells" = "GC_Bcells",
  "Memory B cells" = "Memory_Bcells",
  "PCs" = "PCs"
)

for (subset in names(B_cell_subsets)){
  
  # subset <- "all B cells"
  
  df_plot <- meta %>% 
    filter(L1_annotation %in% B_cell_subsets[[subset]]) %>% 
    mutate(
      has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE),
      patient_id = str_split_i(sample, "-", 1)
    ) 
  
  df_count <- df_plot %>% 
    count(patient_id) 
  
  df_plot %>% 
    ggplot(aes(x = patient_id, fill = has_bcr)) + 
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("grey", "#7FB069")) + 
    geom_text(
      data = df_count,
      aes(x = patient_id, y = 0.95, label = glue("{n} cells")),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    theme_bw() + 
    labs(
      title = glue("Percentage of {subset} that have BCR data available"), 
      x = "Patient ID", 
      y = "Percentage of cells", 
      fill = "Has BCR data"
    ) 
  
  png_string <- str_replace_all(subset, " ", "_")
  ggsave(glue("{outdir1}/bcr_available_{png_string}.png"))
  
}

# ==============================================================================
# N B cells 
# ==============================================================================

outdir2 <- glue("{outdir}/N_B_cells/")
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

# meta <- readRDS("45_immcantation/out/rds/meta_4_Gina_list_90_similarity.rds")

B_cell_subsets <- list(
  "all B cells" = c("GC_Bcells", "Memory_Bcells", "Naive_Bcells", "PCs", "Unconventional_Bcells"), 
  "GC B cells" = "GC_Bcells",
  "Memory B cells" = "Memory_Bcells",
  "PCs" = "PCs"
)

# all cells 
for (subset in names(B_cell_subsets)){
  
  # subset <- "all B cells"
  # subset <- "GC B cells"
    
  df_plot <- meta %>% 
    filter(L1_annotation %in% B_cell_subsets[[subset]]) %>% 
    mutate(
      has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE),
      patient_id = str_split_i(sample, "-", 1)
    )
  
  df_plot %>% 
    ggplot(aes(x = patient_id, fill = has_bcr)) + 
    geom_bar() +
    scale_fill_manual(values = c("grey", "#7FB069")) + 
    theme_bw() + 
    labs(
      title = glue("N {subset}"), 
      x = "Patient ID", 
      y = "N cells", 
      fill = "Has BCR data"
    ) 
  
  png_string <- str_replace_all(subset, " ", "_")
  ggsave(glue("{outdir2}/N_{png_string}.png"))
  
}

# Across tissues
for (subset in names(B_cell_subsets)){
  
  # subset <- "all B cells"
  # subset <- "GC B cells"
    
  df_plot <- meta %>% 
    filter(
      L1_annotation %in% B_cell_subsets[[subset]]
    ) %>% 
    mutate(
      sample_plot = sample %>% str_remove_all("-HLADR-AND-CD19|-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-AND-GC-AND-TFH"),
      has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE)
    )
  
  df_plot %>% 
    ggplot(aes(x = sample_plot, fill = has_bcr)) + 
    geom_bar() +
    scale_fill_manual(values = c("grey", "#7FB069")) + 
    theme_bw() + 
    labs(
      title = glue("N {subset}"), 
      x = "Sample", 
      y = "N cells",
      fill = "Has BCR data"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  png_string <- str_replace_all(subset, " ", "_")
  ggsave(glue("{outdir2}/N_{png_string}_across_tissues.png"))
  
}

# ==============================================================================
# GC B cells in PP 
# ==============================================================================

outdir3 <- glue("{outdir}/PP_GB_cells/")
dir.create(outdir3, recursive = TRUE, showWarnings = FALSE)

# N GC B cells in PPs
df_plot <- meta %>% 
  filter(
    L1_annotation == "GC_Bcells"
  ) %>% 
  mutate(
    sample_plot = sample %>% str_remove_all("-HLADR-AND-CD19|-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-AND-GC-AND-TFH"),
    has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE)
  ) %>% 
  filter(
    sample_plot %in% c("HH117-SI-PP-nonINF", "HH119-SI-PP")
  )

df_plot %>% 
  ggplot(aes(x = sample_plot, fill = has_bcr)) + 
  geom_bar() +
  scale_fill_manual(values = c("grey", "#7FB069")) + 
  theme_bw() + 
  labs(
    title = glue("N GC B cells in Peyer's patches"), 
    subtitle = "ACTUALLY not quite correct since I filter out BCR data for PP cells that is NA in ADT", 
    x = "Sample", 
    y = "N cells",
    fill = "Has BCR data"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(glue("{outdir3}/N_GC_B_cells_in_PP.png"))

# ADT avail
df_plot <- meta %>% 
  filter(
    L1_annotation == "GC_Bcells"
  ) %>% 
  mutate(
    sample_plot = sample %>% str_remove_all("-HLADR-AND-CD19|-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-AND-GC-AND-TFH"),
    has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE),
    has_ADT = ifelse(manual_ADT_class == "Singlet", TRUE, FALSE)
  ) %>% 
  filter(
    sample_plot %in% c("HH117-SI-PP-nonINF", "HH119-SI-PP")
  )

df_plot %>% 
  ggplot(aes(x = sample_plot, fill = has_ADT)) + 
  geom_bar() +
  scale_fill_manual(values = c("grey", "#9A69B0")) + 
  theme_bw() + 
  labs(
    title = glue("N GC B cells in Peyer's patches"), 
    subtitle = "The GC B cells we were able to demultiplex and hence determine their follicle location", 
    x = "Sample", 
    y = "N cells",
    fill = "Has ADT data"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(glue("{outdir3}/N_GC_B_cells_in_PP_ADT.png"))

# N cells per follicle across patients
df_plot <- meta %>% 
  filter(
    L1_annotation == "GC_Bcells"
  ) %>% 
  mutate(
    sample_plot = sample %>% str_remove_all("-HLADR-AND-CD19|-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-AND-GC-AND-TFH"),
    has_ADT = ifelse(manual_ADT_class == "Singlet", TRUE, FALSE),
    has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE)
  ) %>% 
  filter(
    sample_plot %in% c("HH117-SI-PP-nonINF", "HH119-SI-PP"),
    has_ADT
  ) %>% 
  mutate(
    follicle = manual_ADT_ID %>% str_split_i("-", 2) %>% as.integer()
  )

PP_samples <- df_plot$sample_plot %>% unique()

for (PP_sample in PP_samples) {
  
  # PP_sample <- "HH119-SI-PP"

  df_plot %>% 
    filter(sample_plot == PP_sample) %>% 
    count(follicle, has_bcr) %>% 
    ggplot(aes(x = follicle, y = n)) + 
    geom_col(aes(fill = has_bcr)) +
    scale_fill_manual(values = c("grey", "#7FB069")) + 
    geom_text(
      data = df_plot %>% filter(sample_plot == PP_sample) %>% count(follicle),
      aes(label = n), size = 3, vjust = -0.5
    ) + 
    theme_bw() + 
    labs(
      title = glue("{PP_sample}: N GC B cells in Peyer's patch follicles"), 
      x = "Follicle", 
      y = "N cells"
    ) + 
    scale_x_continuous(
      breaks = function(x) seq(1, ceiling(max(x)), by = 1),
      limits = c(0.5, NA),
      expand = c(0, 0.5),
      minor_breaks = scales::breaks_width(1)
    )
  
  ggsave(glue("{outdir3}/N_GC_B_cells_in_PP_follicles_{PP_sample}.png"), width = 10)

}

# ==============================================================================
# N clones
# ==============================================================================

outdir4 <- glue("{outdir}/N_clones/")
dir.create(outdir4, recursive = TRUE, showWarnings = FALSE)

# N GC B cells in PPs
df_plot <- meta %>% 
  mutate(
    sample_plot = sample %>% str_remove_all("-HLADR-AND-CD19|-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-AND-GC-AND-TFH"),
    follicle = manual_ADT_ID %>% str_split_i("-", 2) %>% as.integer()
  ) %>% 
  filter(
    L1_annotation == "GC_Bcells",
    sample_plot %in% c("HH117-SI-PP-nonINF", "HH119-SI-PP"),
    !is.na(manual_ADT_ID), 
    !is.na(clone_subgroup_id_90_similarity)
  ) 

df_plot %>% 
  select(sample_plot, clone_subgroup_id_90_similarity) %>% 
  distinct() %>% 
  count(sample_plot) %>% 
  ggplot(aes(x = sample_plot, y = n)) + 
  geom_col() + 
  theme_bw() + 
  labs(
    title = "N clones for PP samples", 
    y = "N clones"
  )

ggsave(glue("{outdir4}/N_clones_per_sample.png"))
  
# N clones for PP follicles 
PP_samples <- df_plot$sample_plot %>% unique()

for (PP_sample in PP_samples) {
  
  # PP_sample <- "HH119-SI-PP"
  df_plot %>% 
    filter(sample_plot == PP_sample) %>% 
    select(follicle, clone_subgroup_id_90_similarity) %>% 
    distinct() %>% 
    count(follicle) %>% 
    ggplot(aes(x = follicle, y = n)) + 
    geom_col() + 
    geom_text(
      aes(label = n), size = 3, vjust = -0.5
    ) + 
    scale_x_continuous(
      breaks = function(x) seq(1, ceiling(max(x)), by = 1),
      limits = c(0.5, NA),
      expand = c(0, 0.5),
      minor_breaks = scales::breaks_width(1)
    ) +
    theme_bw() + 
    labs(
      title = glue("{PP_sample}: N clones for PP follicles"), 
      x = "Follicle", 
      y = "N clones"
    ) 
  
  ggsave(glue("{outdir4}/N_clones_per_follicle_{PP_sample}.png", width = 10))
    
}


# ==============================================================================
# Clone size graph (Freds plot)
# ==============================================================================

outdir5 <- glue("{outdir}/clone_size/")
dir.create(outdir5, recursive = TRUE, showWarnings = FALSE)

B_cell_subsets <- list(
  "all B cells" = c("GC_Bcells", "Memory_Bcells", "Naive_Bcells", "PCs", "Unconventional_Bcells"), 
  "GC B cells" = "GC_Bcells",
  "Memory B cells" = "Memory_Bcells",
  "PCs" = "PCs"
)

# all cells 
for (subset in names(B_cell_subsets)){

  # subset <- "GC B cells"
  
  # Freqency of clone size
  df_plot <- meta %>% 
    mutate(
      sample_plot = sample %>% str_remove_all("-HLADR-AND-CD19|-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-CD19-Pool1|-CD19-Pool2|-GC-AND-PB-AND-TFH-Pool1|-GC-AND-PB-AND-TFH-Pool2|-AND-GC-AND-TFH"),
      follicle = manual_ADT_ID %>% str_split_i("-", 2) %>% as.integer(),
      has_bcr = ifelse(!is.na(clone_subgroup_id_90_similarity), TRUE, FALSE),
      patient_id = sample %>% str_split_i("-", 1)
    ) %>% 
    filter(
      has_bcr,
      L1_annotation %in% B_cell_subsets[[subset]]
    ) %>% 
    count(patient_id, clone_subgroup_id_90_similarity) %>% 
    dplyr::rename(clone_size = n) %>% 
    mutate(
      clone_size_group = case_when(
        clone_size == 1 ~ "Singleton",
        clone_size > 1 & clone_size <= 5 ~ "2-5", 
        clone_size > 5 & clone_size <= 10 ~ "6-10",
        clone_size > 10 & clone_size <= 20 ~ "11-20",
        clone_size > 20 & clone_size <= 50 ~ "21-50",
        clone_size > 50 & clone_size <= 100 ~ "51-100",
        clone_size > 100 ~ "100+"
      ),
      clone_size_group = factor(clone_size_group, levels = c("Singleton", "2-5", "6-10", "11-20", "21-50", "51-100", "100+"))
    ) %>% 
    count(patient_id, clone_size_group)
  
  df_plot %>% 
    # filter(clone_size_group != "Singleton") %>% 
    ggplot(aes(x = clone_size_group, y = n, color = patient_id)) + 
    geom_point(alpha = 0.5, size = 2) + 
    geom_line(aes(group = patient_id)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(
      title = glue("Frequency of clone sizes of {subset}"),
      x = "Clone size (N cells)", 
      y = "N clones", 
      color = "Patient ID"
    )
  
  png_string <- str_replace_all(subset, " ", "_")
  ggsave(glue("{outdir5}/freq_of_clones_size_{png_string}.png"))
  
}

# ==============================================================================
# All cells: Summary of cell types in follicles.
# ==============================================================================

outdir_1 <- glue("{outdir}/Follicle_cell_types")
dir.create(outdir_1, recursive = TRUE)

lapply(patients, function(HH){
  
  # HH <- "HH117"
  p <- patient_names[[HH]]
  
  seurat_obj <- subset(seurat_integrated, patient == HH)
  
  # Define LP samples
  LP_samples <- grep("LP", seurat_obj[[]]$sample_clean, value = TRUE) %>% unique()
  
  # How many Tfh cells with BCR?
  seurat_obj[[]] %>% filter((!is.na(bcr_productive_contig_1) & !is.na(bcr_productive_contig_2) & L1_annotation == "Tfh_cells")) 
  
  # Clean meta data and prep for plotting 
  seurat_meta_clean <- seurat_obj[[]] %>%  
    mutate(L1_annotation = ifelse(L1_annotation == "GC_Bcells", "GC_B_cells", L1_annotation)) %>% 
    filter(
      (str_detect(L1_annotation, "Contamination", negate = TRUE)), # Remove contamination
      !(!is.na(bcr_productive_contig_1) & !is.na(bcr_productive_contig_2) & L1_annotation == "Tfh_cells"), # Remove Tfh cells with BCR
      !(L1_annotation == "GC_B_cells" & sample_clean %in% LP_samples), # Remove GC B cells in LP samples
    ) %>% mutate(
      sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-")),
      sample_clean_plot = fct_infreq(sample_clean_plot) #%>% fct_rev()
    ) %>% 
    add_count(sample_clean_plot, name = "Count") 
  
  # Across follicles 
  # HH_fol_sample_clean <- seurat_meta_clean %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))
  
  # Count
  if (HH == "HH117"){
    width <- 12 
  } else if (HH == "HH119"){
    width <- 15.5
  }
  png(glue("{outdir_1}/{HH}_N_cells_across_follicles.png"), width = width, height = 7, res = 1000, units = "in")
  
  print(
    seurat_meta_clean %>% 
      filter(!is.na(manual_ADT_ID)) %>% 
      mutate(
        manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
      ) %>% 
      ggplot(aes(x = manual_ADT_ID_plot, fill = L1_annotation)) +
      geom_bar() + 
      scale_fill_manual(
        values = L1_colors, 
        labels = cell_type_names
      ) + 
      scale_x_continuous(
        breaks = function(x) seq(1, ceiling(max(x)), by = 1),
        limits = c(0.5, NA),
        expand = c(0, 0.5)
      ) + 
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Count", 
        # title = glue ("{p}: {HH_fol_sample_clean} follicles"),
        title = glue ("{p}\nPeyer's patch follicles"),
        fill = "Cell type"
      ) + 
      theme(
        plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)
      )
  )
  
  dev.off()
  
})

#   
# 
# # Facet wrap
# 
# # Define LP samples
# LP_samples <- grep("LP", seurat_integrated[[]]$sample_clean, value = TRUE) %>% unique()
# 
# # How many Tfh cells with BCR?
# seurat_integrated[[]] %>% filter((!is.na(bcr_productive_contig_1) & !is.na(bcr_productive_contig_2) & L1_annotation == "Tfh_cells")) 
# 
# # Clean meta data and prep for plotting 
# seurat_meta_clean <- seurat_integrated[[]] %>%  
#   mutate(L1_annotation = ifelse(L1_annotation == "GC_Bcells", "GC_B_cells", L1_annotation)) %>% 
#   filter(
#     (str_detect(L1_annotation, "Contamination", negate = TRUE)), # Remove contamination
#     !(!is.na(bcr_productive_contig_1) & !is.na(bcr_productive_contig_2) & L1_annotation == "Tfh_cells"), # Remove Tfh cells with BCR
#     !(L1_annotation == "GC_B_cells" & sample_clean %in% LP_samples), # Remove GC B cells in LP samples
#   ) %>% mutate(
#     sample_clean_plot = sample_clean %>% str_remove_all(glue("{HH}-")),
#     sample_clean_plot = fct_infreq(sample_clean_plot) #%>% fct_rev()
#   ) %>% 
#   add_count(sample_clean_plot, name = "Count") 
# 
# # Across follicles 
# HH_fol_sample_clean <- seurat_meta_clean %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))
# 
# # Count
# if (HH == "HH117"){
#   width <- 12 
# } else if (HH == "HH119"){
#   width <- 15
# }
# png(glue("{outdir_1}/{HH}_N_cells_across_follicles.png"), width = width, height = 7, res = 1000, units = "in")
# 
# print(
#   seurat_meta_clean %>% 
#     filter(!is.na(manual_ADT_ID)) %>% 
#     mutate(
#       manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
#     ) %>% 
#     ggplot(aes(x = manual_ADT_ID_plot, fill = L1_annotation)) +
#     geom_bar() + 
#     scale_fill_manual(
#       values = L1_colors, 
#       labels = cell_type_names
#     ) + 
#     scale_x_continuous(
#       breaks = function(x) seq(1, ceiling(max(x)), by = 1),
#       limits = c(0.5, NA),
#       expand = c(0, 0.5)
#     ) + 
#     facet_wrap(vars(patient), drop = TRUE) +
#     theme_classic() +
#     labs(
#       x = "Follicle number", 
#       y = "Count", 
#       title = glue ("Cell types across {HH_fol_sample_clean} follicles"),
#       fill = "Cell type"
#     ) + 
#     theme(
#       plot.title = element_text(face = "bold", size = 26),
#       axis.title = element_text(size = 20),
#       axis.text = element_text(size = 16),
#       legend.title = element_text(size = 20),
#       legend.text = element_text(size = 16)
#     )
# )
# 
# dev.off()

# ==============================================================================
# G B cells: Summary of follicles and isotypes
# ==============================================================================

outdir_2 <- glue("{outdir}/Follicle_GC_B_cells_isotypes")
dir.create(outdir_2, recursive = TRUE)

lapply(patients, function(HH){
  
  # HH <- "HH119"
  p <- patient_names[[HH]]
  
  plot_df <- resolve_LC_list[[HH]] %>% 
    filter(
      locus == "IGH",
      !is.na(manual_ADT_ID), 
      L1_annotation == "GC_B_cells",
      !is.na(c_call)
    ) %>% 
    mutate(
      manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
    ) %>%
    add_count(manual_ADT_ID_plot, name = "Count") 
  
  # Across follicles 
  # HH_fol_sample_clean <- plot_df %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))
  
  # Isotype
  
  ## Freq
  if (HH == "HH117"){
    width <- 12 
  } else if (HH == "HH119"){
    width <- 15.5
  }
  
  png(glue("{outdir_2}/{HH}_Isotype_freq_across_follicles.png"), width = width, height = 7, res = 1000, units = "in")
  
  print(
    plot_df %>% 
      filter(
        !is.na(manual_ADT_ID), 
        L1_annotation == "GC_B_cells",
        !is.na(c_call)
      ) %>% 
      ggplot(aes(x = manual_ADT_ID_plot, fill = c_call)) + 
      geom_bar(position = "fill") + 
      geom_text(
        aes(x = manual_ADT_ID_plot, y = 1.02, label = Count)
      ) + 
      scale_fill_manual(values = isotype_colors_custom) +
      scale_y_continuous(labels = scales::percent) +
      scale_x_continuous(
        breaks = function(x) seq(1, ceiling(max(x)), by = 1),
        limits = c(0.5, NA),
        expand = c(0, 0.5)
      ) + 
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Frequency", 
        title = glue("{p}\nGC B cells from Peyer's patch follicles"),
        fill = "Isotype"
      ) + 
      theme(
        plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)
      )
  )
  
  dev.off()
  
})


# ------------------------------------------------------------------------------
# top GC clones 
# ------------------------------------------------------------------------------

# Define top GC clone
top_GC_clones <- lapply(patients, function(HH) {
  
  # find clones that have GC cells in at least 2 different sample_ids
  GC_clones <- resolve_LC_list[[HH]] %>%
    filter(locus == "IGH") %>% 
    filter(L1_annotation == "GC_B_cells") %>%
    count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
    head(10) %>% 
    pull(clone_subgroup_id_90_similarity)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# N clones barplot
# ------------------------------------------------------------------------------

plot_df <- resolve_LC_list %>%
  bind_rows() %>% 
  filter(
    locus == "IGH",
    L1_annotation == "GC_B_cells", 
    !is.na(clone_subgroup_id_90_similarity)
    # !is.na(manual_ADT_ID)
  ) %>% 
  mutate(
    manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
  ) %>%
  count(patient_id, clone_subgroup_id_90_similarity) %>% 
  dplyr::rename(clone_size = n) %>% 
  mutate(
    clone_size_group = case_when(
      clone_size == 1 ~ "Singletons", 
      clone_size > 1 & clone_size <= 5 ~ "2-5", 
      clone_size > 5 & clone_size <= 10 ~ "6-10",
      clone_size > 10 & clone_size <= 20 ~ "11-20",
      clone_size > 20 & clone_size <= 50 ~ "21-50",
      clone_size > 50 & clone_size <= 100 ~ "51-100",
      clone_size > 100 ~ "100+"
    ),
    clone_size_group = factor(clone_size_group, levels = c("Singletons", "2-5", "6-10", "11-20", "21-50", "51-100", "100+"))
  ) %>% 
  count(patient_id, clone_size_group)
  # select(patient_id, clone_subgroup_id_90_similarity, clone_size_group) %>%
  # distinct() %>%
  # count(patient_id)

plot_df %>% 
  ggplot(aes(x = patient_id, y = n, fill = clone_size_group)) +
  geom_col() + 
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  theme_bw() + 
  labs(
    title = "N GC B cell clones per patient",
    subtitle = "Cells where follicle could not be determined are included",
    x = "Patient",
    y = "N clones", 
    fill = "Clone size group"
  )

ggsave(glue("{outdir}/N_clones_per_patient.png"))


# ==============================================================================
# G B cells: Summary of follicles and isotypes - CRC without the big clone
# ==============================================================================

HH <- "HH119"
p <- patient_names[[HH]]
large_clone <- top_GC_clones[[HH]][[1]]
# large_clone <- top_GC_clones[[HH]][c(1,2)]

plot_df <- resolve_LC_list[[HH]] %>% 
  filter(
    locus == "IGH",
    !is.na(manual_ADT_ID), 
    L1_annotation == "GC_B_cells",
    !is.na(c_call),
    !(clone_subgroup_id_90_similarity %in% large_clone)
  ) %>% 
  mutate(
    manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
  ) %>%
  add_count(manual_ADT_ID_plot, name = "Count") 


# Isotype
## Freq
png(glue("{outdir_2}/{HH}_Isotype_freq_across_follicles_rm_large_clone.png"), width = 15.5, height = 7, res = 1000, units = "in")
# png(glue("{outdir_2}/{HH}_Isotype_freq_across_follicles_rm_large_clone_2.png"), width = 15.5, height = 7, res = 1000, units = "in")

print(
  plot_df %>% 
    filter(
      !is.na(manual_ADT_ID), 
      L1_annotation == "GC_B_cells",
      !is.na(c_call)
    ) %>% 
    ggplot(aes(x = manual_ADT_ID_plot, fill = c_call)) + 
    geom_bar(position = "fill") + 
    geom_text(
      aes(x = manual_ADT_ID_plot, y = 1.02, label = Count)
    ) + 
    scale_fill_manual(values = isotype_colors_custom) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(
      breaks = function(x) seq(1, ceiling(max(x)), by = 1),
      limits = c(0.5, NA),
      expand = c(0, 0.5)
    ) + 
    theme_classic() +
    labs(
      x = "Follicle number", 
      y = "Frequency", 
      title = glue("{p}\nGC B cells from Peyer's patch follicles - Largest clones removed"),
      # title = glue("{p}\nGC B cells from Peyer's patch follicles - Two largest clones removed"),
      fill = "Isotype"
    ) + 
    theme(
      plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16)
    )
)

dev.off()

# ==============================================================================
# Frequency of top clone per follicle 
# ==============================================================================

outdir_6 <- glue("{outdir}/Follicle_GC_B_cells_freq_barplot")
dir.create(outdir_6, recursive = TRUE)

n_clones <- 10

lapply(patients, function(HH){
  
  # HH <- "HH119"
  p <- patient_names[[HH]]
  
  # Subset clones
  top_GC_clones_subset <- top_GC_clones[[HH]][c(1:n_clones)]
  
  plot_df <- resolve_LC_list[[HH]] %>% 
    filter(
      locus == "IGH", 
      L1_annotation == "GC_B_cells",
      !is.na(manual_ADT_ID)
    ) %>% 
    mutate(
      manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer(),
      clone_subgroup_id_90_similarity_plot = ifelse(clone_subgroup_id_90_similarity %in% top_GC_clones_subset, clone_subgroup_id_90_similarity, "other"),
      clone_subgroup_id_90_similarity_plot = factor(clone_subgroup_id_90_similarity_plot, levels = c(top_GC_clones_subset, "other"))
    ) %>%
    add_count(manual_ADT_ID_plot, name = "Count") 
  
  # Across follicles 
  # HH_fol_sample_clean <- plot_df %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))
  
  # Define clone colors 
  clone_colors <- list(
    "#E05C8A", "#66CC55", "#5588DD", "#EE9944", "#AA3377",
    "#44BBAA", "#CC6644", "#4499CC", "#AACC33", "#9955BB",
    # "#FF0000", "#0000FF", "#00CC00", "#FF6600", "#9900CC",
    # "#00CCCC", "#FF0099", "#996600", "#0099FF", "#669900",
    "grey85"
  ) %>% setNames(c(top_GC_clones_subset, "other"))
  
  # Define clone names
  clone_names <- c(paste("Clone", 1:n_clones), "Other") %>% as.list() %>% setNames(c(top_GC_clones_subset, "other"))
  
  # N clones 
  N_clones_per_fol <- plot_df %>%
    filter(
      !is.na(manual_ADT_ID)
    ) %>%
    mutate(
      manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
    ) %>%
    group_by(manual_ADT_ID_plot) %>%
    count(clone_subgroup_id_90_similarity) %>%
    count(manual_ADT_ID_plot) %>%
    ungroup() %>%
    complete(
      manual_ADT_ID_plot = seq(min(manual_ADT_ID_plot), max(manual_ADT_ID_plot)),
      fill = list(n = 0)
    ) 
  
  # colnames(N_clones_per_fol) <- c("Follicle", "N clones")
  # 
  # ggtexttable(N_clones_per_fol, rows = NULL, theme = ttheme("classic"))
  # # grid.text(
  # #   glue("{p}: N clones per follicle"),
  # #   x = 0.50, y = 0.97,          # adjust position as needed
  # #   gp = gpar(fontsize = 20, fontface = "bold")
  # # )
  # ggsave(glue("{outdir_6}/{HH}_N_clones_table.png"), dpi = 1000, height = 10)
  # 
  # N clones 
  
  if (HH == "HH117"){
    width <- 12 
  } else if (HH == "HH119"){
    width <- 15.5
  }
  
  png(glue("{outdir_6}/{HH}_N_{n_clones}.png"), width = width, height = 7, units = "in", res = 1000)
  
  print(
    plot_df %>%
      filter(!is.na(manual_ADT_ID)) %>%
      ggplot(aes(x = manual_ADT_ID_plot)) + 
      geom_bar(aes(fill = clone_subgroup_id_90_similarity_plot), position = "fill") + 
      # geom_text(
      #   # data = N_clones_per_fol, 
      #   aes(x = manual_ADT_ID_plot, y = 1.02, label = Count)
      # ) +
      geom_text(
        data = N_clones_per_fol,
        aes(x = manual_ADT_ID_plot, y = 1.02, label = n)
      ) +
      scale_fill_manual(
        values = clone_colors, 
        labels = clone_names
      ) + 
      scale_x_continuous(
        breaks = function(x) seq(1, ceiling(max(x)), by = 1),
        limits = c(0.5, NA),
        expand = c(0, 0.5)
      ) + 
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Frequency", 
        # title = glue("{p}: Top 10 clones across GC B cells in {HH_fol_sample_clean} follicles"),
        # title = glue("{p}\nTop 10 clones across GC B cells from Peyer's patch follicles"),
        title = glue("{p}\nTop {n_clones} GC B cell clones in Peyer's patch follicles"),
        # subtitle = glue("Top {n_clones} clones highlighted and number of clones with in each follicle is stated on top of the bars"),
        fill = "Clone"
      ) + 
      theme(
        plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)
      )
  )
  
  dev.off()
  
  
  
})

# ==============================================================================
# Frequency of top clone per follicle - junction sequence 
# ==============================================================================

# outdir_6 <- glue("{outdir}/15_poster_figures/Follicle_GC_B_cells_freq_barplot")
# dir.create(outdir_6, recursive = TRUE)
# 
n_clones <- 10

clone_colors_all <- list(
  "HH117" = c(
    "#E05C8A", "#66CC55", "#5588DD", "#EE9944", "#AA3377",
    "#44BBAA", "#CC6644", "#4499CC", "#AACC33", "#9955BB",
    "grey85"
  ), 
  "HH119" = c(
    "#00CCCC", "#FF0099", "#996600", "#0099FF", "#669900",
    "#FF0000", "#0000FF", "#00CC00", "#FF6600", "#9900CC",
    "grey85"
  )
) 

lapply(patients, function(HH){
  
  # HH <- "HH119"
  p <- patient_names[[HH]]
  
  # Subset clones
  top_GC_clones_subset <- top_GC_clones[[HH]][c(1:n_clones)]
  
  plot_df <- resolve_LC_list[[HH]] %>% 
    filter(
      locus == "IGH", 
      L1_annotation == "GC_B_cells",
      !is.na(manual_ADT_ID)
    ) %>% 
    mutate(
      manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer(),
      clone_subgroup_id_90_similarity_plot = ifelse(clone_subgroup_id_90_similarity %in% top_GC_clones_subset, clone_subgroup_id_90_similarity, "other"),
      clone_subgroup_id_90_similarity_plot = factor(clone_subgroup_id_90_similarity_plot, levels = c(top_GC_clones_subset, "other"))
    ) %>%
    add_count(manual_ADT_ID_plot, name = "Count") 
  
  # Across follicles 
  # HH_fol_sample_clean <- plot_df %>% filter(!is.na(manual_ADT_ID)) %>% pull(sample_clean) %>% unique() %>% str_remove(glue("{HH}-"))
  
  # Define clone colors 
  clone_colors <- clone_colors_all[[HH]] %>% setNames(c(top_GC_clones_subset, "other"))
  
  # Define clone names
  # clone_names <- c(paste("Clone", 1:n_clones), "Other") %>% as.list() %>% setNames(c(top_GC_clones_subset, "other"))
  
  # Define majority junction sequence as clone name 
  clone_names <- resolve_LC_list[[HH]] %>% 
    filter(
      locus == "IGH", clone_subgroup_id_90_similarity %in% top_GC_clones_subset
    ) %>% 
    count(clone_subgroup_id_90_similarity, junction, sort = TRUE) %>% 
    group_by(clone_subgroup_id_90_similarity) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(-n) %>% 
    deframe() %>% 
    as.list()
  
  clone_names <- c(clone_names, "other"= "Other")
  
  
  # N clones 
  N_clones_per_fol <- plot_df %>%
    filter(
      !is.na(manual_ADT_ID)
    ) %>%
    mutate(
      manual_ADT_ID_plot = str_split_i(manual_ADT_ID, "-", 2) %>% as.integer()
    ) %>%
    group_by(manual_ADT_ID_plot) %>%
    count(clone_subgroup_id_90_similarity) %>%
    count(manual_ADT_ID_plot) %>%
    ungroup() %>%
    complete(
      manual_ADT_ID_plot = seq(min(manual_ADT_ID_plot), max(manual_ADT_ID_plot)),
      fill = list(n = 0)
    ) 
  
  # colnames(N_clones_per_fol) <- c("Follicle", "N clones")
  # 
  # ggtexttable(N_clones_per_fol, rows = NULL, theme = ttheme("classic"))
  # # grid.text(
  # #   glue("{p}: N clones per follicle"),
  # #   x = 0.50, y = 0.97,          # adjust position as needed
  # #   gp = gpar(fontsize = 20, fontface = "bold")
  # # )
  # ggsave(glue("{outdir_6}/{HH}_N_clones_table.png"), dpi = 1000, height = 10)
  # 
  # N clones 
  
  if (HH == "HH117"){
    width <- 15
  } else if (HH == "HH119"){
    width <- 20
  }
  
  png(glue("{outdir_6}/{HH}_N_{n_clones}_sequences.png"), width = width, height = 7, units = "in", res = 1000)
  
  print(
    plot_df %>%
      filter(!is.na(manual_ADT_ID)) %>%
      ggplot(aes(x = manual_ADT_ID_plot)) + 
      geom_bar(aes(fill = clone_subgroup_id_90_similarity_plot), position = "fill") + 
      # geom_text(
      #   # data = N_clones_per_fol, 
      #   aes(x = manual_ADT_ID_plot, y = 1.02, label = Count)
      # ) +
      geom_text(
        data = N_clones_per_fol,
        aes(x = manual_ADT_ID_plot, y = 1.02, label = n)
      ) +
      scale_fill_manual(
        values = clone_colors, 
        labels = clone_names
      ) + 
      scale_x_continuous(
        breaks = function(x) seq(1, ceiling(max(x)), by = 1),
        limits = c(0.5, NA),
        expand = c(0, 0.5)
      ) + 
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      labs(
        x = "Follicle number", 
        y = "Frequency", 
        # title = glue("{p}: Top 10 clones across GC B cells in {HH_fol_sample_clean} follicles"),
        # title = glue("{p}\nTop 10 clones across GC B cells from Peyer's patch follicles"),
        title = glue("{p}\nTop {n_clones} GC B cell clones in Peyer's patch follicles"),
        # subtitle = glue("Top {n_clones} clones highlighted and number of clones with in each follicle is stated on top of the bars"),
        fill = "Clone"
      ) + 
      theme(
        plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16)
        # legend.title = element_text(size = 20),
        # legend.text = element_text(size = 16)
      )
  )
  
  dev.off()
  
  
  
})
