# Load libraries
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
# suppressPackageStartupMessages(library(dowser)) # Needs to be installed 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(scoper))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(shazam))
library(tibble)
library(patchwork)
library(forcats)
library(glue)
library(stringr)

packageVersion("airr")
packageVersion("alakazam")
packageVersion("scoper")
packageVersion("shazam")

# Following this Immcantation flow:
# https://immcantation.readthedocs.io/en/latest/getting_started/10x_tutorial.html

# ------------------------------------------------------------------------------
# Load annotations 
# ------------------------------------------------------------------------------

L3_GCB_annotation <- readRDS("00_data/GCB_meta_GL.rds")
nrow(L3_GCB_annotation)

# Create outdir
outdir <- "50_BACH2/plot/01_bcr_avail/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Get sample names - Heavy chain 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

base_path <- "45_immcantation/out"
files <- list.files(base_path, full.names = FALSE)
sample_names <- files[1:length(files)-1]

# Read all samples, adding sample_id and subject_id
bcr_data_tmp <- lapply(sample_names, function(s) {
  # s <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  f <- file.path(base_path, s, paste0(s, "_heavy_germ-pass.tsv"))
  db <- airr::read_rearrangement(f)
  db$sample_id <- s
  db$subject_id <- sub("-(SI|SILP|CO|COLP).*", "", s)  # extracts HH117 or HH119
  # make sequence and cell IDs unique across samples
  db$sequence_id <- paste0(s, "_", db$sequence_id)
  db$cell_id <- paste0(s, "_", db$cell_id)
  return(db)
}) %>% bind_rows()

# Add annotations 
table(L3_GCB_annotation$cell_id %in% bcr_data_tmp$cell_id)

table(L3_GCB_annotation$L3_GCB_annotation)
L3_GCB_annotation[L3_GCB_annotation$cell_id %in% bcr_data_tmp$cell_id, ] %>% pull(L3_GCB_annotation) %>% table()

bcr_data_tmp <- bcr_data_tmp %>% left_join(L3_GCB_annotation %>% select(cell_id, L3_GCB_annotation, manual_ADT_ID), by = "cell_id")
bcr_data_tmp$L3_GCB_annotation %>% table()

bcr_data_tmp$cell_id %>% length()
bcr_data_tmp$cell_id %>% unique() %>% length()

# split by patient
bcr_data <- list(
  "HH117" = bcr_data_tmp %>% filter(subject_id == "HH117"),
  "HH119" = bcr_data_tmp %>% filter(subject_id == "HH119"),
  "HH151" = bcr_data_tmp %>% filter(subject_id == "HH151"),
  "HH153" = bcr_data_tmp %>% filter(subject_id == "HH153")
)

patients <- names(bcr_data)

for (HH in patients){
  cat(paste(HH, "-", nrow(bcr_data[[HH]]), "sequences\n"))
}


# ------------------------------------------------------------------------------
# Remove non-productive sequences
# ------------------------------------------------------------------------------

for (HH in patients){
  print(paste(HH, "-", bcr_data[[HH]]$productive %>% table()))
}

# ------------------------------------------------------------------------------
# Handle multiple heavy chains 
# ------------------------------------------------------------------------------

# Visualize UMIs of multiple contigs (heavy chains)
# lapply(patients, function(HH){
#   
#   # HH <- "HH153"
#   
#   umi_pairs <- bcr_data[[HH]] %>% 
#     group_by(cell_id) %>% 
#     arrange(desc(umi_count), .by_group = TRUE) %>% 
#     filter(n() == 2) %>% 
#     summarise(
#       top_umi = umi_count[1], 
#       second_umi = umi_count[2],
#       ratio = umi_count[1] / umi_count[2],
#       .groups = "drop"
#     )
#   
#   ggplot(umi_pairs, aes(x = top_umi, y = second_umi, color = ratio >= 2)) + 
#     geom_jitter(alpha = 0.4, width = 0.1, height = 0.1) + 
#     scale_x_log10() + 
#     scale_y_log10() + 
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") + 
#     theme_bw() + 
#     labs(
#       title = glue("{HH}: Top vs second heavy chain contig UMI count"), 
#       x = "Top contig UMI (log scale)", y = "Second contig UMI (log scale)",
#       color = "Ratio ≥ 2"
#     )
#   
#   ggsave(glue("{outdir}/{HH}_two_contigs.png"))
#   
#   
# })

min_second_umi_noise <- 3  # from the flat bottom band in the plot; adjust if you want it tighter/looser

bcr_data_qc <- lapply(bcr_data, function(x) {
  
  # x <- bcr_data[["HH117"]]
  
  df <- x %>%
    group_by(cell_id) %>%
    arrange(desc(umi_count), .by_group = TRUE) %>%
    mutate(
      n_heavy = n(),
      same_rearrangement = n_distinct(v_call, j_call, junction) == 1,
      # only IGHM+IGHD specifically counts as benign co-expression
      is_md_pair = n_heavy == 2 & all(c_call %in% c("IGHM", "IGHD")) & n_distinct(c_call) == 2,
      dominant = case_when(
        n_heavy == 1 ~ TRUE,
        same_rearrangement & is_md_pair ~ c_call == "IGHM",
        # everything else with 2 heavy contigs (different rearrangements, or same rearrangement but not IGHM/IGHD) -> same ratio/noise logic
        n_heavy == 2 & umi_count[2] <= min_second_umi_noise ~ row_number() == 1,
        n_heavy == 2 & umi_count[2] > min_second_umi_noise & umi_count[1] >= 2 * umi_count[2] ~ row_number() == 1,
        TRUE ~ FALSE
      )
    ) %>%
    filter(dominant) %>%
    select(-n_heavy, -same_rearrangement, -is_md_pair, -dominant) %>%
    ungroup() %>% 
    mutate(
      c_call_grouped = if_else(c_call %in% c("IGHM", "IGHD"), "IGHM/D", c_call)
    )
  
  return(df)
  
})

# Rows in the data before filtering 
for (HH in patients){
  cat(paste(HH, "-", nrow(bcr_data[[HH]]), "sequences\n"))
}

# Rows in the data after filtering 
for (HH in patients){
  cat(paste(HH, "-", nrow(bcr_data_qc[[HH]]), "sequences\n"))
}

# Rows in the data after filtering 
for (HH in patients){
  cat(paste(HH, "-", length(bcr_data_qc$HH$cell_id) == unique(length(bcr_data_qc$HH$cell_id)), "that sequences are unique\n"))
}

# ------------------------------------------------------------------------------
# BCR availability stats - heavy chain 
# ------------------------------------------------------------------------------

# GEX cells
table(L3_GCB_annotation$L3_GCB_annotation, L3_GCB_annotation$Patient)
# L3_GCB_annotation[L3_GCB_annotation$cell_id %in% bcr_data_tmp$cell_id, ] %>% select(L3_GCB_annotation, Patient) %>% table()

# BCR QC
bcr_data_qc$HH117$L3_GCB_annotation %>% table()
bcr_data_qc$HH119$L3_GCB_annotation %>% table()

# bcr_data_qc$HH119 %>%
#   filter(
#     # !is.na(c_call),
#     # !is.na(manual_ADT_ID), 
#     !(manual_ADT_ID %in% c("Negative", "Doublet"))
#   ) %>% 
#   pull(L3_GCB_annotation) %>% 
#   table()

# bcr_data_qc$HH119 %>% filter(!is.na(c_call)) %>% pull(L3_GCB_annotation) %>% table()

# Bar plot: x-axis = L3_GCB_annotation, facet_wrap on patient, dodge barplot 
# GEX counts
gex_counts <- L3_GCB_annotation %>%
  filter(
    Patient %in% names(bcr_data_qc), 
    !(manual_ADT_ID %in% c("Negative", "Doublet"))
  ) %>%
  count(Patient, L3_GCB_annotation, name = "n") %>%
  mutate(Source = "GEX")

# BCR (post-QC) counts
bcr_counts <- bcr_data_qc %>%
  imap_dfr(~ .x %>%
             filter(
               !(manual_ADT_ID %in% c("Negative", "Doublet"))
             )
           count(L3_GCB_annotation, name = "n") %>%
             mutate(Patient = .y)) %>%
  mutate(Source = "BCR (QC)")

# BCR (post-QC) counts with non-NA c_call
bcr_ccall_counts <- bcr_data_qc %>%
  imap_dfr(~ .x %>%
             filter(
               !is.na(c_call), 
               !(manual_ADT_ID %in% c("Negative", "Doublet")),
             ) %>%
             count(L3_GCB_annotation, name = "n") %>%
             mutate(Patient = .y)) %>%
  mutate(Source = "BCR (QC, c_call)")

# BCR (post-QC) counts with non-NA ADT_manual_ID
bcr_ADT_counts <- bcr_data_qc %>%
  imap_dfr(~ .x %>%
             filter(
               !is.na(c_call),
               !(manual_ADT_ID %in% c("Negative", "Doublet")),
               !is.na(manual_ADT_ID)
             ) %>%
             count(L3_GCB_annotation, name = "n") %>%
             mutate(Patient = .y)) %>%
  mutate(Source = "BCR (QC, ADT)")

source_levels <- c("GEX", "BCR (QC)", "BCR (QC, c_call)", "BCR (QC, ADT)")

plot_df <- bind_rows(gex_counts, bcr_counts, bcr_ccall_counts, bcr_ADT_counts) %>%
  mutate(Source = factor(Source, levels = source_levels)) %>%
  complete(Patient, L3_GCB_annotation, Source, fill = list(n = 0)) %>% 
  filter(
    !is.na(L3_GCB_annotation), 
    Patient %in% c("HH117", "HH119")
  )

ggplot(plot_df, aes(x = L3_GCB_annotation, y = n, fill = Source)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.85),
            vjust = -0.3, size = 2.5) +
  facet_wrap(~ Patient, ncol = 1) +
  scale_fill_manual(values = c("GEX" = "grey60",
                               "BCR (QC)" = "#2C7FB8",
                               "BCR (QC, c_call)" = "#F28E2B",
                               "BCR (QC, ADT)" = "forestgreen")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey90")) + 
  labs(
    title = "Heavy chain GC B cell filtering", 
    subtitle = "ADT Negative and Doublet removed.",
    x = "L3 GCB annotation", y = "Number of cells", fill = NULL
  )

ggsave(glue("{outdir}/bcr_avail_heavy_chain.png"), width = 9, height = 8)

# Save data 
stats_data <- list(
  "heavy" = list(
    "bcr_data_tmp" = bcr_data_tmp,
    "bcr_data" = bcr_data, 
    "bcr_data_qc" = bcr_data_qc
  )
)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Get sample names - Light chain 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Clean up 
rm(bcr_data_tmp, bcr_data, bcr_data_qc)

base_path <- "45_immcantation/out"
files <- list.files(base_path, full.names = FALSE)
sample_names <- files[1:length(files)-1]

# Read all samples, adding sample_id and subject_id
bcr_data_tmp <- lapply(sample_names, function(s) {
  # s <- "HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
  f <- file.path(base_path, s, paste0(s, "_light_germ-pass.tsv"))
  db <- airr::read_rearrangement(f)
  db$sample_id <- s
  db$subject_id <- sub("-(SI|SILP|CO|COLP).*", "", s)  # extracts HH117 or HH119
  # make sequence and cell IDs unique across samples
  db$sequence_id <- paste0(s, "_", db$sequence_id)
  db$cell_id <- paste0(s, "_", db$cell_id)
  return(db)
}) %>% bind_rows()

# Add annotations 
table(L3_GCB_annotation$cell_id %in% bcr_data_tmp$cell_id)

table(L3_GCB_annotation$L3_GCB_annotation)
L3_GCB_annotation[L3_GCB_annotation$cell_id %in% bcr_data_tmp$cell_id, ] %>% pull(L3_GCB_annotation) %>% table()

bcr_data_tmp <- bcr_data_tmp %>% left_join(L3_GCB_annotation %>% select(cell_id, L3_GCB_annotation, manual_ADT_ID), by = "cell_id")
bcr_data_tmp$L3_GCB_annotation %>% table()

bcr_data_tmp$cell_id %>% length()
bcr_data_tmp$cell_id %>% unique() %>% length()

# split by patient
bcr_data <- list(
  "HH117" = bcr_data_tmp %>% filter(subject_id == "HH117"),
  "HH119" = bcr_data_tmp %>% filter(subject_id == "HH119"),
  "HH151" = bcr_data_tmp %>% filter(subject_id == "HH151"),
  "HH153" = bcr_data_tmp %>% filter(subject_id == "HH153")
)

patients <- names(bcr_data)

for (HH in patients){
  cat(paste(HH, "-", nrow(bcr_data[[HH]]), "sequences\n"))
}


# ------------------------------------------------------------------------------
# Remove non-productive sequences
# ------------------------------------------------------------------------------

for (HH in patients){
  print(paste(HH, "-", bcr_data[[HH]]$productive %>% table()))
}

# ------------------------------------------------------------------------------
# Handle multiple light chains 
# ------------------------------------------------------------------------------

# Filter cells based on multiple light chain
bcr_data_qc <- lapply(bcr_data, function(x) {
  
  df <- x %>%
    group_by(cell_id) %>%
    arrange(desc(umi_count), .by_group = TRUE) %>%
    mutate(
      n_light = n(),
      dominant = case_when(
        n_light == 1 ~ TRUE,                          # only one light chain, keep it
        umi_count[1] >= 2 * umi_count[2] ~ row_number() == 1,  # dominant contig has 2x UMIs
        TRUE ~ FALSE                                  # ambiguous, drop all contigs for this cell
      )
    ) %>%
    filter(dominant) %>%
    select(-n_light, -dominant) %>%
    ungroup()
  
  return(df)
  
}
)

# Rows in the data before filtering 
for (HH in patients){
  cat(paste(HH, "-", nrow(bcr_data[[HH]]), "sequences\n"))
}

# Rows in the data after filtering 
for (HH in patients){
  cat(paste(HH, "-", nrow(bcr_data_qc[[HH]]), "sequences\n"))
}

# Rows in the data after filtering 
for (HH in patients){
  cat(paste(HH, "-", length(bcr_data_qc$HH$cell_id) == unique(length(bcr_data_qc$HH$cell_id)), "that sequences are unique\n"))
}

# ------------------------------------------------------------------------------
# BCR availability stats - light chain 
# ------------------------------------------------------------------------------

# GEX cells
table(L3_GCB_annotation$L3_GCB_annotation, L3_GCB_annotation$Patient)
# L3_GCB_annotation[L3_GCB_annotation$cell_id %in% bcr_data_tmp$cell_id, ] %>% select(L3_GCB_annotation, Patient) %>% table()

# BCR QC
bcr_data_qc$HH117$L3_GCB_annotation %>% table()
bcr_data_qc$HH119$L3_GCB_annotation %>% table()

# Add a column with cells that are not NA in c_call of bcr_data_qc

# Bar plot: x-axis = L3_GCB_annotation, y = counts, facet_wrap on patient
# GEX counts
gex_counts <- L3_GCB_annotation %>%
  filter(
    Patient %in% names(bcr_data_qc), 
    !(manual_ADT_ID %in% c("Negative", "Doublet"))
  ) %>%
  count(Patient, L3_GCB_annotation, name = "n") %>%
  mutate(Source = "GEX")

# BCR (post-QC) counts
bcr_counts <- bcr_data_qc %>%
  imap_dfr(~ .x %>%
             filter(
               !(manual_ADT_ID %in% c("Negative", "Doublet"))
             ) %>% 
             count(L3_GCB_annotation, name = "n") %>%
             mutate(Patient = .y)) %>%
  mutate(Source = "BCR (QC)")

# BCR (post-QC) counts with non-NA c_call
bcr_ccall_counts <- bcr_data_qc %>%
  imap_dfr(~ .x %>%
             filter(
               !is.na(c_call), 
               !(manual_ADT_ID %in% c("Negative", "Doublet")),
             ) %>%
             count(L3_GCB_annotation, name = "n") %>%
             mutate(Patient = .y)) %>%
  mutate(Source = "BCR (QC, c_call)")

# BCR (post-QC) counts with non-NA ADT_manual_ID
bcr_ADT_counts <- bcr_data_qc %>%
  imap_dfr(~ .x %>%
             filter(
               !is.na(c_call),
               !(manual_ADT_ID %in% c("Negative", "Doublet")),
               !is.na(manual_ADT_ID)
             ) %>%
             count(L3_GCB_annotation, name = "n") %>%
             mutate(Patient = .y)) %>%
  mutate(Source = "BCR (QC, ADT)")

source_levels <- c("GEX", "BCR (QC)", "BCR (QC, c_call)", "BCR (QC, ADT)")

plot_df <- bind_rows(gex_counts, bcr_counts, bcr_ccall_counts, bcr_ADT_counts) %>%
  mutate(Source = factor(Source, levels = source_levels)) %>%
  complete(Patient, L3_GCB_annotation, Source, fill = list(n = 0)) %>% 
  filter(
    !is.na(L3_GCB_annotation), 
    Patient %in% c("HH117", "HH119")
  )

ggplot(plot_df, aes(x = L3_GCB_annotation, y = n, fill = Source)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.85),
            vjust = -0.3, size = 2.5) +
  facet_wrap(~ Patient, ncol = 1) +
  scale_fill_manual(values = c("GEX" = "grey60",
                               "BCR (QC)" = "#2C7FB8",
                               "BCR (QC, c_call)" = "#F28E2B",
                               "BCR (QC, ADT)" = "forestgreen")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey90")) + 
  labs(
    title = "Light chain GC B cell filtering",
    subtitle = "ADT Negative and Doublet removed.",
    x = "L3 GCB annotation", y = "Number of cells", fill = NULL
  )

ggsave(glue("{outdir}/bcr_avail_light_chain.png"), width = 9, height = 8)

# Save data 
stats_data_light <- list(
  "light" = list(
    "bcr_data_tmp" = bcr_data_tmp,
    "bcr_data" = bcr_data, 
    "bcr_data_qc" = bcr_data_qc
  )
)

stats_data <- c(stats_data, stats_data_light)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Stats for combined heavy and light chain 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ADT filter applied to all counts (NA in manual_ADT_ID is kept)
filter_adt <- function(df) {
  df %>% filter(!(manual_ADT_ID %in% c("Negative", "Doublet")))
}

# Count unique cells per annotation and patient
count_cells <- function(df_list, source) {
  df_list %>%
    imap_dfr(~ .x %>%
               distinct(cell_id, L3_GCB_annotation) %>%
               count(L3_GCB_annotation, name = "n") %>%
               mutate(Patient = .y)) %>%
    mutate(Source = source)
}

# ADT filtered heavy and light QC data
heavy_qc <- map(stats_data$heavy$bcr_data_qc, filter_adt)
light_qc <- map(stats_data$light$bcr_data_qc, filter_adt)

# GEX counts (ADT filtered)
gex_counts <- L3_GCB_annotation %>%
  filter_adt() %>%
  count(Patient, L3_GCB_annotation, name = "n") %>%
  mutate(Source = "GEX")

# Cells passing QC for both heavy and light chain (per patient)
paired_cells <- imap(heavy_qc, ~ intersect(.x$cell_id, light_qc[[.y]]$cell_id))

# Paired data: heavy + light rows for paired cells only (Immcantation format)
bcr_paired <- imap(heavy_qc, ~ bind_rows(
  .x %>% filter(cell_id %in% paired_cells[[.y]]),
  light_qc[[.y]] %>% filter(cell_id %in% paired_cells[[.y]])
))

# Paired cells, using heavy chain rows (one row per cell, heavy c_call = isotype)
paired_heavy_qc <- imap(heavy_qc, ~ .x %>% filter(cell_id %in% paired_cells[[.y]]))

# Cumulative filtering steps
paired_ccall  <- map(paired_heavy_qc, ~ filter(.x, !is.na(c_call)))
paired_adt_id <- map(paired_ccall, ~ filter(.x, !is.na(manual_ADT_ID)))

source_levels_paired <- c("GEX", "Heavy (QC)", "Light (QC)", "Paired (QC)",
                          "Paired (QC, heavy c_call)", "Paired (QC, heavy c_call, ADT ID)")

plot_df <- bind_rows(
  gex_counts,
  count_cells(heavy_qc, "Heavy (QC)"),
  count_cells(light_qc, "Light (QC)"),
  count_cells(paired_heavy_qc, "Paired (QC)"),
  count_cells(paired_ccall, "Paired (QC, heavy c_call)"),
  count_cells(paired_adt_id, "Paired (QC, heavy c_call, ADT ID)")
) %>%
  mutate(Source = factor(Source, levels = source_levels_paired)) %>%
  complete(Patient, L3_GCB_annotation, Source, fill = list(n = 0)) %>%
  filter(
    !is.na(L3_GCB_annotation),
    Patient %in% c("HH117", "HH119")
  )

ggplot(plot_df, aes(x = L3_GCB_annotation, y = n, fill = Source)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 2) +
  facet_wrap(~ Patient, ncol = 1) +
  scale_fill_manual(values = c("GEX" = "grey60",
                               "Heavy (QC)" = "#2C7FB8",
                               "Light (QC)" = "#41AB5D",
                               "Paired (QC)" = "#F28E2B",
                               "Paired (QC, heavy c_call)" = "#B15928",
                               "Paired (QC, heavy c_call, ADT ID)" = "#6A3D9A")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey90")) +
  labs(
    title = "Paired heavy and light chain GC B cell filtering",
    subtitle = "ADT Negative and Doublet removed from all counts.",
    x = "L3 GCB annotation", y = "Number of cells", fill = NULL
  )

ggsave(glue("{outdir}/bcr_avail_paired.png"), width = 12, height = 8)

# Save
stats_data$paired <- list(
  "paired_cells" = paired_cells,
  "bcr_paired" = bcr_paired,
  "bcr_paired_adt_id" = imap(bcr_paired, ~ .x %>% filter(cell_id %in% paired_adt_id[[.y]]$cell_id))
)