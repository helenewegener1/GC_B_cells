library(glue)
library(tidyverse)
library(alakazam)
library(scatterpie)
library(patchwork)

source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC_3_definitions", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
patients

# Prep output
outdir = glue("45_immcantation/plot/18_90_similarity/08_clonal_sharing/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# For each clone, count presents in follicles (GC B cells)
# ------------------------------------------------------------------------------

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>%
    filter(
      locus == "IGH"
      # !is.na(manual_ADT_full_ID)
    ) %>%
    mutate(patient = HH)
}) %>% bind_rows()

n_min_cells <- 2

df_n_fol <- df_both %>% 
  filter(
    L1_annotation == "GC_B_cells",
    !is.na(manual_ADT_full_ID)
  ) %>% # should it only be across GC B cells?
  group_by(clone_subgroup_id_90_similarity, patient_id) %>% 
  count(manual_ADT_full_ID) %>% 
  filter(n >= n_min_cells) %>%
  count(manual_ADT_full_ID) %>% 
  count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
  dplyr::rename(n_follicles = n) %>% 
  ungroup()

# N 
df_N <- df_n_fol %>%
  select(-clone_subgroup_id_90_similarity) %>% 
  count(patient_id, n_follicles, sort = TRUE)

df_n_fol %>% 
  ggplot(aes(x = patient_id, y = n_follicles)) + 
  geom_jitter(width = 0.20, height = 0) + 
  geom_text(
    data = df_N,
    aes(x = patient_id, y = n_follicles, label = glue("{n} clones")),
    size = 3,
    position = position_nudge(x = 0.35)
  ) +
  theme_bw() + 
  scale_y_continuous(
    breaks = scales::breaks_width(5),
    minor_breaks = scales::breaks_width(1)
  ) + 
  labs(
    x = "Patient ID",
    y = "N follicles",
    title = "Clonal sharing across follicles",
    subtitle = glue("For each clone (GC B cells), how many follicles is it present in (at least {n_min_cells} cells)")
  )

ggsave(glue("{outdir}/clones_in_n_follicles.png"), width = 8, height = 6)

# # ------------------------------------------------------------------------------
# # Different isotypes in different follicles? (GC B cells)
# # ------------------------------------------------------------------------------
# 
# df_shared_clones <- df_n_fol %>% filter(n_follicles > 1)
# 
# # ---- 1. prep counts per clone x follicle x isotype ----
# 
# plot_data <- df_both %>% 
#   filter(
#     !is.na(c_call),
#     L1_annotation == "GC_B_cells",
#     !is.na(manual_ADT_full_ID)
#   ) %>% 
#   mutate(manual_ADT_ID_plot = str_remove(manual_ADT_full_ID, "Fol-") %>% as.integer() %>% as.factor()) %>% 
#   inner_join(df_shared_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
#   count(patient_id, clone_subgroup_id_90_similarity, manual_ADT_ID_plot, c_call)
# 
# clone_id_lookup <- plot_data %>% 
#   distinct(patient_id, clone_subgroup_id_90_similarity) %>% 
#   arrange(patient_id, clone_subgroup_id_90_similarity) %>% 
#   group_by(patient_id) %>% 
#   mutate(y = row_number()) %>% 
#   ungroup()
# 
# pie_data <- plot_data %>% 
#   mutate(x = as.numeric(as.character(manual_ADT_ID_plot))) %>% 
#   left_join(clone_id_lookup, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
#   pivot_wider(
#     id_cols = c(patient_id, clone_subgroup_id_90_similarity, x, y),
#     names_from = c_call, 
#     values_from = n, 
#     values_fill = 0
#   )
# 
# isotype_cols <- setdiff(names(pie_data), c("patient_id", "clone_subgroup_id_90_similarity", "x", "y"))
# 
# # ---- 2. one plot per patient, each with its own axis labels and coord_equal ----
# 
# make_pie_plot <- function(HH) {
#   
#   # HH <- "HH117"
#   data_sub <- pie_data %>% filter(patient_id == HH)
#   lookup_sub <- clone_id_lookup %>% filter(patient_id == HH)
#   
#   ggplot() + 
#     geom_scatterpie(
#       data = data_sub, 
#       aes(x = x, y = y, r = 0.35),
#       cols = isotype_cols, 
#       color = NA
#     ) + 
#     scale_fill_manual(
#       values = isotype_colors_custom,
#       breaks = c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2",  "IGHG3", "IGHG4", "IGHE")
#     ) + 
#     scale_x_continuous(
#       # breaks = sort(unique(data_sub$x)),
#       breaks = 1:max(data_sub$x),
#       minor_breaks = scales::breaks_width(1)
#     ) +
#     scale_y_continuous(
#       breaks = lookup_sub$y,
#       labels = lookup_sub$clone_subgroup_id_90_similarity,
#       minor_breaks = scales::breaks_width(1)
#     ) +
#     coord_equal() + 
#     labs(
#       x = "Follicle", y = "Clone", fill = "Isotype", 
#       title = glue("{HH}: Isotype usage in follicle-shared clones")
#     ) + 
#     theme_bw()
# }
# 
# # patient_names
# 
# for (HH in names(patient_names)){
#   # HH <- "HH119"
#   make_pie_plot(HH)
#   ggsave(glue("{outdir}/{HH}_isotype_follicle_shared_clones.png"), width = 10, height = 8)
#   
# }

# ------------------------------------------------------------------------------
# Different isotypes in different follicles and sites? (All cells )
# ------------------------------------------------------------------------------

# ---- 1. build an ordered site lookup per patient ----
site_lookup <- df_both %>% 
  # filter(
  #   !is.na(c_call)
  # ) %>% 
  distinct(patient_id, sample_clean_fol) %>% 
  mutate(
    follicle_num = str_extract(sample_clean_fol, "(?<=Fol-)\\d+") %>% as.integer(),
    is_follicle = !is.na(follicle_num)
  ) %>% 
  arrange(patient_id, is_follicle, follicle_num, sample_clean_fol) %>% 
  group_by(patient_id) %>% 
  mutate(x = row_number()) %>% 
  ungroup() %>% 
  mutate(
    sample_clean_fol_plot = sample_clean_fol %>% str_remove("HH11\\d+-") %>% str_remove("SI-PP_|SI-PP-nonINF_")
  ) %>% 
  select(patient_id, sample_clean_fol_plot, x)
  

# ---- 2. prep counts, joining in x from the lookup instead of parsing numerically ----

plot_data <- df_both %>% 
  inner_join(df_shared_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  mutate(
    sample_clean_fol_plot = sample_clean_fol %>% str_remove("HH11\\d+-") %>% str_remove("SI-PP_|SI-PP-nonINF_")
  ) %>% 
  count(patient_id, clone_subgroup_id_90_similarity, sample_clean_fol_plot, c_call) %>% 
  left_join(site_lookup, by = c("patient_id", "sample_clean_fol_plot"))

clone_id_lookup <- plot_data %>% 
  distinct(patient_id, clone_subgroup_id_90_similarity) %>% 
  arrange(patient_id, clone_subgroup_id_90_similarity) %>% 
  group_by(patient_id) %>% 
  mutate(y = row_number()) %>% 
  ungroup()

isotype_cols <- setdiff(names(pie_data), c("patient_id", "clone_subgroup_id_90_similarity", "x", "y", "total_n"))

pie_data <- plot_data %>% 
  left_join(clone_id_lookup, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  pivot_wider(
    id_cols = c(patient_id, clone_subgroup_id_90_similarity, x, y),
    names_from = c_call, 
    values_from = n, 
    values_fill = 0
  ) %>% 
  # add a total_n column to pie_data first
  mutate(total_n = rowSums(across(all_of(isotype_cols))))

# ---- 3. plot function, now labeling x with site names instead of numbers ----

make_pie_plot <- function(patient_name) {
  
  data_sub <- pie_data %>% filter(patient_id == patient_name)
  lookup_sub <- clone_id_lookup %>% filter(patient_id == patient_name)
  x_lookup_sub <- site_lookup %>% filter(patient_id == patient_name)
  
  ggplot() + 
    geom_scatterpie(
      data = data_sub, 
      aes(x = x, y = y, r = 0.35),
      cols = isotype_cols, 
      color = NA
    ) + 
    geom_text(
      data = data_sub,
      aes(x = x, y = y, label = total_n),
      size = 2.5,
      nudge_y = -0.5  # push label below the pie so it doesn't overlap the wedges
    ) +
    scale_fill_manual(
      values = isotype_colors_custom,
      breaks = c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2",  "IGHG3", "IGHG4", "IGHE")
    ) + 
    scale_x_continuous(
      breaks = x_lookup_sub$x, 
      labels = x_lookup_sub$sample_clean_fol_plot,
      minor_breaks = scales::breaks_width(1)
    ) + 
    scale_y_continuous(
      breaks = lookup_sub$y, 
      labels = lookup_sub$clone_subgroup_id_90_similarity,
      minor_breaks = scales::breaks_width(1)
    ) + 
    coord_equal() + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(
      x = "Site", y = "Clone", fill = "Isotype", 
      title = glue("{HH}: Isotype usage in follicle-shared clones - all cells")
    )
}

# patient_names
width <- 12
for (HH in names(patient_names)){
  # HH <- "HH117"
  make_pie_plot(HH)
  
  if (HH == "HH119"){
    width <- width + 2
  }
  
  ggsave(glue("{outdir}/{HH}_isotype_follicle_shared_clones_all_cells.png"), width = width, height = 8)
  
}

# ------------------------------------------------------------------------------
# Different isotypes in different follicles and sites? (GC B cells )
# ------------------------------------------------------------------------------

# ---- 1. build an ordered site lookup per patient ----
site_lookup <- df_both %>% 
  filter(
    # !is.na(c_call),
    L1_annotation == "GC_B_cells"
  ) %>% 
  distinct(patient_id, sample_clean_fol) %>% 
  mutate(
    follicle_num = str_extract(sample_clean_fol, "(?<=Fol-)\\d+") %>% as.integer(),
    is_follicle = !is.na(follicle_num)
  ) %>% 
  arrange(patient_id, is_follicle, follicle_num, sample_clean_fol) %>% 
  group_by(patient_id) %>% 
  mutate(x = row_number()) %>% 
  ungroup() %>% 
  mutate(
    sample_clean_fol_plot = sample_clean_fol %>% str_remove("HH11\\d+-") %>% str_remove("SI-PP_|SI-PP-nonINF_")
  ) %>% 
  select(patient_id, sample_clean_fol_plot, x)


# ---- 2. prep counts, joining in x from the lookup instead of parsing numerically ----

plot_data <- df_both %>% 
  filter(L1_annotation == "GC_B_cells") %>% 
  inner_join(df_shared_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  mutate(
    sample_clean_fol_plot = sample_clean_fol %>% str_remove("HH11\\d+-") %>% str_remove("SI-PP_|SI-PP-nonINF_")
  ) %>% 
  count(patient_id, clone_subgroup_id_90_similarity, sample_clean_fol_plot, c_call) %>% 
  left_join(site_lookup, by = c("patient_id", "sample_clean_fol_plot"))

clone_id_lookup <- plot_data %>% 
  distinct(patient_id, clone_subgroup_id_90_similarity) %>% 
  arrange(patient_id, clone_subgroup_id_90_similarity) %>% 
  group_by(patient_id) %>% 
  mutate(y = row_number()) %>% 
  ungroup()

isotype_cols <- setdiff(names(pie_data), c("patient_id", "clone_subgroup_id_90_similarity", "x", "y", "total_n"))

pie_data <- plot_data %>% 
  left_join(clone_id_lookup, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
  pivot_wider(
    id_cols = c(patient_id, clone_subgroup_id_90_similarity, x, y),
    names_from = c_call, 
    values_from = n, 
    values_fill = 0
  ) %>% 
  # add a total_n column to pie_data first
  mutate(total_n = rowSums(across(all_of(isotype_cols))))

# ---- 3. plot function, now labeling x with site names instead of numbers ----

make_pie_plot <- function(patient_name) {
  
  data_sub <- pie_data %>% filter(patient_id == patient_name)
  lookup_sub <- clone_id_lookup %>% filter(patient_id == patient_name)
  x_lookup_sub <- site_lookup %>% filter(patient_id == patient_name)
  
  ggplot() + 
    geom_scatterpie(
      data = data_sub, 
      aes(x = x, y = y, r = 0.35),
      cols = isotype_cols, 
      color = NA
    ) + 
    geom_text(
      data = data_sub,
      aes(x = x, y = y, label = total_n),
      size = 2.5,
      nudge_y = -0.5  # push label below the pie so it doesn't overlap the wedges
    ) +
    scale_fill_manual(
      values = isotype_colors_custom,
      breaks = c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2",  "IGHG3", "IGHG4", "IGHE")
    ) + 
    scale_x_continuous(
      breaks = x_lookup_sub$x, 
      labels = x_lookup_sub$sample_clean_fol_plot,
      minor_breaks = scales::breaks_width(1)
    ) + 
    scale_y_continuous(
      breaks = lookup_sub$y, 
      labels = lookup_sub$clone_subgroup_id_90_similarity,
      minor_breaks = scales::breaks_width(1)
    ) + 
    coord_equal() + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(
      x = "Site", y = "Clone", fill = "Isotype", 
      title = glue("{HH}: Isotype usage in follicle-shared clones - GC B cells")
    )
}

# patient_names
width <- 12
for (HH in names(patient_names)){
  # HH <- "HH117"
  make_pie_plot(HH)
  
  if (HH == "HH119"){
    width <- width + 2
  }
  ggsave(glue("{outdir}/{HH}_isotype_follicle_shared_clones_GC_B_cells.png"), width = width, height = 8)
  
}

