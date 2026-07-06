library(glue)
library(tidyverse)
library(alakazam)

source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

rds_files <- list.files("45_immcantation/out/rds") 
resolve_LC_files <- grep("resolve_LC_3_definitions", rds_files, value = TRUE)

patients <- lapply(resolve_LC_files, function(x) str_split_i(x, "_", 1)) %>% unlist()
patients

# HH <- "HH119"

# Read rds
# df_heavy <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>% 
#   filter(
#     locus == "IGH",
#     !is.na(manual_ADT_full_ID)
#   )

# df_heavy$clone_subgroup_id_90_similarity
# df_heavy$manual_ADT_full_ID

# Remove largest clone as it "takes all the signal"
# df_heavy <- df_heavy %>% filter(clone_subgroup_id_90_similarity != "20693_1")
# extra <- "_largest_removed"

# Prep output
outdir = glue("45_immcantation/plot/18_90_similarity/08_clonal_sharing/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# Prep colors 
HH_samples <- names(sample_clean_plot_colors) %>% str_subset(glue("^{HH}")) %>% str_subset("Fol")
HH_samples_colors <- sample_clean_plot_colors[HH_samples]
names(HH_samples_colors) <- names(HH_samples_colors) %>% str_split_i("_", 2)
HH_samples_colors

# ------------------------------------------------------------------------------
# For each clone, count presents in follicles 
# ------------------------------------------------------------------------------

# Load both patients
df_both <- lapply(patients, function(HH) {
  readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds")) %>%
    filter(
      locus == "IGH",
      !is.na(manual_ADT_full_ID)
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
  rename(n_follicles = n)

# N 
df_N <- df_n_fol %>%
  ungroup() %>% 
  select(-clone_subgroup_id_90_similarity) %>% 
  count(patient_id, n_follicles, sort = TRUE)

df_n_fol %>% 
  left_join(df_N_clones, by = c("patient_id", "clone_subgroup_id_90_similarity")) %>% 
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


