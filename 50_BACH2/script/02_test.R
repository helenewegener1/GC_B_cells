library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

L3_GCB_annotation_2 <- readRDS("00_data/GCB_meta_GL.rds")
L3_GCB_annotation_4 <- readRDS("00_data/GCB_meta_GL_new.rds")

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC\\.", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 2)) %>% unlist()
patients

# Load both patients
df_all <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/05_{HH}_resolve_LC.rds")) %>%
    filter(
      locus == "IGH"
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows()

# Define largest CRC clone
large_crc_clone <- df_all %>%
  filter(patient_id == "HH119", L1_annotation == "GC_B_cells") %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  head(1) %>% 
  pull(clone_subgroup_id_90_similarity)

# Create outdit
outdir <- glue("50_BACH2/plot/02_test")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Wrangle data
# ------------------------------------------------------------------------------

L3_GCB_annotation_4$cell_id <- glue("{L3_GCB_annotation_4$Sample_name}_{rownames(L3_GCB_annotation_4)}") %>% str_remove("_\\d+")

df_all$cell_id 
L3_GCB_annotation_4$cell_id

df_all$cell_id %>% str_split_i("_", 1) %>% unique()
L3_GCB_annotation_4$cell_id %>% str_split_i("_", 1) %>% unique()

table(L3_GCB_annotation_4$cell_id %in% df_all$cell_id)
L3_GCB_annotation_4$cell_id[!(L3_GCB_annotation_4$cell_id %in% df_all$cell_id)]

L3_GCB_annotation_4_clean <- L3_GCB_annotation_4 %>% 
  select(cell_id, L2_GCB_annotation, L3_GCB_annotation) %>% 
  rename(L2_GCB_annotation_4 = L2_GCB_annotation, L3_GCB_annotation_4 = L3_GCB_annotation) 
rownames(L3_GCB_annotation_4_clean) <- NULL

df_all <- df_all %>% left_join(L3_GCB_annotation_4_clean, by = "cell_id")

L3_GCB_annotation_2_clean <- L3_GCB_annotation_2 %>% 
  select(cell_id, L2_GCB_annotation, L3_GCB_annotation) %>% 
  rename(L2_GCB_annotation_2 = L2_GCB_annotation, L3_GCB_annotation_2 = L3_GCB_annotation) 
rownames(L3_GCB_annotation_2_clean) <- NULL

df_all <- df_all %>% left_join(L3_GCB_annotation_2_clean, by = "cell_id")

saveRDS(df_all, "50_BACH2/out/df_all.rds")

# ------------------------------------------------------------------------------
# Overview
# ------------------------------------------------------------------------------

df_all_GC <- df_all %>% filter(L1_annotation == "GC_B_cells")
table(df_all_GC$patient, df_all_GC$L3_GCB_annotation_4)

table(df_all_GC$patient_id, df_all_GC$L3_GCB_annotation_4, useNA = "always")

table(df_all_GC$c_call_grouped, df_all_GC$L3_GCB_annotation_4, useNA = "always")

# Largest clone
df_all %>% 
  filter(patient_id == "HH119", clone_subgroup_id_90_similarity == large_crc_clone) %>% 
  pull(L3_GCB_annotation_4) %>% 
  table()


# # 03_heavy_bcr_data_qc_annot.rds
# data <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds") 
# 
# HH <- "HH119"
# GC_B_cells <- L3_GCB_annotation_clean %>% filter(!is.na(L3_GCB_annotation), L3_GCB_annotation == "BACH2+")
# df_heavy <- data[[HH]] %>% filter(locus == "IGH") 
# 
# table(GC_B_cells$cell_id %in% df_heavy$cell_id)
# 
# cell_id_heavy <- GC_B_cells$cell_id[GC_B_cells$cell_id %in% df_heavy$cell_id]
# length(cell_id_heavy)
# 
# # 04_bcr_heavy_light.rds
# all_combined_list <- readRDS("45_immcantation/out/rds/04_bcr_heavy_light.rds")
# 
# HH <- "HH119"
# GC_B_cells <- L3_GCB_annotation_clean %>% filter(!is.na(L3_GCB_annotation), L3_GCB_annotation == "BACH2+")
# df_heavy <- all_combined_list[[HH]] %>% filter(locus == "IGH") 
# df_light <- all_combined_list[[HH]] %>% filter(locus != "IGH") 
# 
# table(GC_B_cells$cell_id %in% df_heavy$cell_id)
# table(GC_B_cells$cell_id %in% df_light$cell_id)
# 
# cell_id_heavy <- GC_B_cells$cell_id[GC_B_cells$cell_id %in% df_heavy$cell_id]
# cell_id_light <- GC_B_cells$cell_id[GC_B_cells$cell_id %in% df_light$cell_id]
# 
# length(cell_id_heavy)
# intersect(cell_id_heavy, cell_id_light) %>% length()

# df <- all_combined_list$HH119 %>% left_join(L3_GCB_annotation_clean, by = "cell_id")


# ------------------------------------------------------------------------------
# Clones
# ------------------------------------------------------------------------------

df_all %>% count(clone_subgroup_id_90_similarity, sort = TRUE) %>% head(5)

df_all %>% 
  filter(L3_GCB_annotation_4 == "BACH2+") %>% 
  count(patient_id, clone_subgroup_id_90_similarity, sort = TRUE) 

df_all %>% 
  filter(L3_GCB_annotation_4 == "BACH2+") %>% 
  count(patient_id, clone_subgroup_id_90_similarity, c_call_grouped, sort = TRUE) %>% 
  head(20)


# ------------------------------------------------------------------------------
# Compare old and new annotations 
# ------------------------------------------------------------------------------

old <- L3_GCB_annotation_2 %>% select(cell_id, Patient, L1_annotation, L3_GCB_annotation)
new <- L3_GCB_annotation_4 %>% select(cell_id, Patient, L1_annotation, L3_GCB_annotation)

old$L1_annotation %>% table()
new$L1_annotation %>% table()

old$L3_GCB_annotation %>% table()
new$L3_GCB_annotation %>% table()

table(old$Patient, old$L1_annotation)
table(new$Patient, new$L1_annotation)

table(old$Patient, old$L3_GCB_annotation)
table(new$L3_GCB_annotation, new$Patient)

old_bach2_id <- old %>% filter(L3_GCB_annotation == "BACH2+") %>% pull(cell_id)
new_bach2_id <- new %>% filter(L3_GCB_annotation == "BACH2+") %>% pull(cell_id)

old_bach2_cyc_id <- old %>% filter(L3_GCB_annotation == "BACH2_cycling") %>% pull(cell_id)
new_bach2_cyc_id <- new %>% filter(L3_GCB_annotation == "BACH2_cycling") %>% pull(cell_id)

length(old_bach2_id)
length(new_bach2_id)

table(new_bach2_id %in% old_bach2_id)
table(old_bach2_id %in% new_bach2_id)

table(new_bach2_id %in% old_bach2_cyc_id)
table(old_bach2_cyc_id %in% new_bach2_id)

table(new_bach2_cyc_id %in% old_bach2_id)
table(old_bach2_id %in% new_bach2_cyc_id)

old_bach2_id[!(old_bach2_id %in% new_bach2_id)] %>% str_split_i("_", 1) %>% table()
old_bach2_id[!(old_bach2_id %in% new_bach2_id)] %>% str_split_i("_", 1) %>% str_split_i("-", 1) %>% table()

# Combined 
old_id <- c(old_bach2_id, old_bach2_cyc_id)
new_id <- c(new_bach2_id, new_bach2_cyc_id)

table(old_id %in% new_id)
table(new_id %in% old_id)

old_id[!(old_id %in% new_id)] %>% str_split_i("_", 1) %>% str_split_i("-", 1) %>% table()

new_id[!(new_id %in% old_id)] %>% str_split_i("_", 1) %>% str_split_i("-", 1) %>% table()

table(new$patient_id, new$L3_GCB_annotation)

# ------------------------------------------------------------------------------
# Compare old and new annotations in df_all
# ------------------------------------------------------------------------------

# HH <- "HH119"

order <- df_all$L3_GCB_annotation_2 %>% unique() %>% as.character() %>% na.exclude() %>% sort()

for (HH in c("HH117", "HH119")){
  
  df_HH <- df_all %>% filter(patient_id == HH)
  
  # table(df_HH$L3_GCB_annotation_2, df_HH$L3_GCB_annotation_4, useNA = "always")
  
  df <- as.data.frame(
    table(
      annotation_2_patients = df_HH$L3_GCB_annotation_2,
      annotation_4_patients = df_HH$L3_GCB_annotation_4, 
      useNA = "always"
    )
  ) %>% 
    mutate(
      # Aviod super high number of NA/NA to mess with the colors 
      Freq = ifelse(is.na(annotation_2_patients) & is.na(annotation_4_patients), -1, Freq) %>% as.integer(),
      annotation_2_patients = factor(annotation_2_patients, levels = order),
      annotation_4_patients = factor(annotation_4_patients, levels = order) 
    )
  
  ggplot(df, aes(x = annotation_2_patients, y = annotation_4_patients, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 5, color = "grey80") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      x = "Annotation (2 patients)",
      y = "Annotation (4 patients)",
      fill = "Cell count", 
      title = glue("{HH}: Old VS new GC B cell subset annotation"), 
      caption = "NA in both annotation is marked -1 even though the number is high."
    ) +
    theme_minimal(base_size = 14)
  
  ggsave(glue("{outdir}/{HH}_annotation_differences_2_4.png"), width = 15, height = 7.5)
  
}


# Investigate the interesting cases by looking at all cells 

# HH119, annotation_2 == BACH2+, annotation_4 == NA
subset <- df_all %>% filter(patient_id == "HH119", L3_GCB_annotation_2 == "BACH2+", is.na(L3_GCB_annotation_4)) 

nrow(subset)

subset %>% select(L1_annotation) %>% table()

subset %>% count(clone_subgroup_id_90_similarity, sort = TRUE)
large_crc_clone

# Conclusion: 150 of the BACH2 cells from annotation_2 are NA in annotation_4 AND these cells are of the large clone and are removed during annotation_4. 

# HH119, annotation_2 == BACH2+, annotation_4 == NA
subset <- df_all %>% filter(patient_id == "HH119", L3_GCB_annotation_2 == "BACH2+", L3_GCB_annotation_4 == "LZ_transitioning") 

nrow(subset)

subset %>% select(L1_annotation) %>% table()

subset %>% count(clone_subgroup_id_90_similarity, sort = TRUE)
large_crc_clone

# Conclusion: 25 of the BACH2 cells from annotation_2 are LZ_transitioning in annotation_4 AND these cells are of the large clone. 


# Look at the annotation of the large clone
subset <- df_all %>% filter(clone_subgroup_id_90_similarity == large_crc_clone)

# table(subset$L3_GCB_annotation_2, subset$L3_GCB_annotation_4, tab)

df <- as.data.frame(
  table(
    annotation_2_patients = subset$L3_GCB_annotation_2,
    annotation_4_patients = subset$L3_GCB_annotation_4, 
    useNA = "always"
  )
) %>% 
  mutate(
    # Aviod super high number of NA/NA to mess with the colors 
    Freq = ifelse(is.na(annotation_2_patients) & is.na(annotation_4_patients), -1, Freq) %>% as.integer(),
    annotation_2_patients = factor(annotation_2_patients, levels = order),
    annotation_4_patients = factor(annotation_4_patients, levels = order) 
  )
  

ggplot(df, aes(x = annotation_2_patients, y = annotation_4_patients, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, color = "grey80") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Annotation (2 patients)",
    y = "Annotation (4 patients)",
    fill = "Cell count", 
    title = glue("{HH}, large clone: Old VS new GC B cell subset annotation"), 
    caption = "NA in both annotation is marked -1 even though the number is high."
  ) +
  theme_minimal(base_size = 14)

ggsave(glue("{outdir}/{HH}_rm_large_clone_annotation_differences_2_4.png"), width = 15, height = 7.5)

  