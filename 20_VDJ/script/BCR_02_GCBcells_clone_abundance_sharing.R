getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(readxl)
library(tidyr)
library(tibble)
library(purrr)

source("10_broad_annotation/script/color_palette.R")
# celltype_colors

combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")
# combined.BCR.joined <- readRDS("20_VDJ/out/combined.BCR.joined.rds")

# names(combined.BCR.filtered)
# names(combined.BCR.joined)

names(combined.BCR.filtered) <- names(combined.BCR.filtered) %>% 
  str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")

# Filter for germinal center B cells 
combined.BCR.filtered_GCB <- lapply(names(combined.BCR.filtered), function(x){
  combined.BCR.filtered[[x]] %>% filter(celltype_broad == "GC_B_cells")
}) %>% setNames(names(combined.BCR.filtered))

# Check number of GC B cells in each sample. 
lapply(names(combined.BCR.filtered_GCB), function(x) combined.BCR.filtered_GCB[[x]] %>% nrow())

# ------------------------------------------------------------------------------
# Clonal abundance per sample
# ------------------------------------------------------------------------------
# 
# # Count of clone sizes
# clone_abundance <- lapply(combined.BCR.filtered, function(x) x %>% 
#   summarise(abundance = n(), .by = c(CTstrict)) %>% arrange(desc(abundance)))
# 
# # Check numbers
# combined.BCR.filtered$`HH117-SI-PP-nonINF_Fol-1` %>% nrow() # N cells
# combined.BCR.filtered$`HH117-SI-PP-nonINF_Fol-1`$CTstrict %>% unique() %>% length() # N unique clones
# clone_abundance$`HH117-SI-PP-nonINF_Fol-1` %>% nrow() # N unique clones
# 
# # Export as excel 
# out_file <- "20_VDJ/table/BCR_clonal_abundance.xlsx"
# 
# # Use openxlsx::write.xlsx, which takes the named list and writes
# # each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
# openxlsx::write.xlsx(
#   x = clone_abundance,
#   file = out_file, 
#   overwrite = TRUE # Overwrite the file if it already exists
# )

# ------------------------------------------------------------------------------
# Number of unique clones per sample
# ------------------------------------------------------------------------------

# Count unique clones per sample
n_unique_clones <- data.frame(
  sample = names(combined.BCR.filtered_GCB),
  n_clones = sapply(combined.BCR.filtered_GCB, function(x) length(unique(x$CTstrict))),
  n_cells = sapply(combined.BCR.filtered_GCB, nrow)
) %>%
  mutate(
    patient = str_extract(sample, "HH\\d+"),
    tissue = case_when(
      str_detect(sample, "SI-PP") ~ "PP",
      str_detect(sample, "SI-MILF") ~ "MILF",
      str_detect(sample, "SILP") ~ "SILP",
      TRUE ~ "Other"
    ),
    condition = case_when(
      str_detect(sample, "nonINF") ~ "non-inflamed",
      str_detect(sample, "INF") ~ "inflamed",
      TRUE ~ "NA"
    ),
    sample_type = str_extract(sample, "Fol-\\d+|MILF|SILP|PP")
  )

# Plot 1: Bar plot of unique clones per sample
ggplot(n_unique_clones, aes(x = sample, y = n_clones)) +
  geom_col() +
  geom_text(aes(label = n_clones), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Number of Unique Clones per Sample",
    subtitle = "Only GC B cells",
    x = "",
    y = "Number of unique clones"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave("20_VDJ/plot/BCR_clonal_sharing_GCB/BCR_final_N_unique_clones_GCB.png", 
       width = 14, height = 15, dpi = 300)

# ------------------------------------------------------------------------------
# Number of cells colored by CTstrict 
# ------------------------------------------------------------------------------

for (sample_name in c("HH117-SI-PP", "HH119-SI-PP")){

  # sample_name <- "HH117-SI-PP"
  # sample_name <- "HH119-SI-PP"
  
  # Subset sample
  combined.BCR.filtered_subset <- combined.BCR.filtered_GCB[str_detect(names(combined.BCR.filtered_GCB), sample_name)]
  combined.BCR.filtered_subset <- bind_rows(combined.BCR.filtered_subset)
  
  combined.BCR.filtered_subset <- combined.BCR.filtered_subset %>% 
    mutate(CTstrict_abundance = n(), .by = c(sample, CTstrict)) %>% 
    mutate(
      Fol = str_split_i(sample, "_", 2),
      CTstrict_above_2 = ifelse(CTstrict_abundance >= 2, CTstrict, NA),
      CTstrict_above_5 = ifelse(CTstrict_abundance >= 5, CTstrict, NA)
    )
  
  combined.BCR.filtered_subset %>% 
    ggplot(aes(
      x = Fol,
      fill = CTstrict_above_2
      # fill = CTstrict_above_5
    )) +
    geom_bar() + 
    labs(
      title = "N cells colored by CTstrict",
      subtitle = sample_name,
      x = "",
      y = "N cells", 
      caption = "CTstrict with abundance >= 2 are colored"
    ) + 
    theme_bw() + 
    theme(legend.position = "none")
  
  ggsave(glue("20_VDJ/plot/BCR_clonal_sharing_GCB/BCR_N_cells_CTstrict_{sample_name}_GCB.png"), 
         width = 14, height = 7, dpi = 300)
  
}

# ------------------------------------------------------------------------------
# Clonal abundance per sample per broad cell type 
# ------------------------------------------------------------------------------

# Count of clone sizes
clone_abundance_celltype <- lapply(combined.BCR.filtered_GCB, function(x) x %>% 
                                     summarise(abundance = n(), .by = c(CTstrict, celltype_broad)) %>% 
                                     mutate(abundance_total = sum(abundance), .by = CTstrict) %>% 
                                     arrange(desc(abundance_total), desc(abundance))
                                   )

# Check numbers
combined.BCR.filtered_GCB$`HH117-SI-PP-nonINF_Fol-1` %>% nrow() # N GCB cells
combined.BCR.filtered_GCB$`HH117-SI-PP-nonINF_Fol-1` %>% select(CTstrict, celltype_broad) %>% distinct() %>% nrow() # N unique clones
clone_abundance_celltype$`HH117-SI-PP-nonINF_Fol-1` %>% nrow() # N unique clones

# Plot
for (sample_name in names(clone_abundance_celltype)){
  
  # sample_name <- "HH117-SI-PP-nonINF_Fol-13"
  
  min_abundance <- ifelse(str_detect(sample_name, "Fol"), 1, 7)
  
  # Filter data
  plot_data <- clone_abundance_celltype[[sample_name]] %>% 
    filter(abundance_total > min_abundance)
  
  # Create plot
  plot_data %>% 
    ggplot(aes(x = abundance, y = reorder(CTstrict, abundance_total), fill = celltype_broad)) + 
    geom_col() + 
    geom_text(
      data = plot_data %>% select(CTstrict, abundance_total) %>% distinct(),
      aes(x = abundance_total, y = CTstrict, label = abundance_total),  # Added x and y here!
      hjust = -0.1, 
      size = 3, 
      inherit.aes = FALSE
    ) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = celltype_colors) + 
    labs(
      title = "Clonal abundance by cell type",
      subtitle = sample_name, 
      x = "Clone size (total cells)",
      y = "",
      fill = "Cell type",
      caption = glue("Clones with > {min_abundance} cells shown")
    ) + 
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  ggsave(glue("20_VDJ/plot/BCR_clonal_sharing_GCB/{sample_name}_celltypes_GCB.png"), width = 20, height = 10)
  
}

# ------------------------------------------------------------------------------
# TFH with BCR 
# ------------------------------------------------------------------------------

TFH_w_BCR <- sapply(clone_abundance_celltype, function(x) x %>% filter(celltype_broad == "Tfh_like_cells") %>% nrow())

TFH_w_BCR %>% data.frame()

# ------------------------------------------------------------------------------
# Track clone manually  
# ------------------------------------------------------------------------------

# clone <- "IGH:Cluster.4529.IGHV4-34_IGLC:Cluster.172.IGLV1-40"
# 
# clone_abundance_celltype$`HH119-SILP` %>% filter(CTstrict == clone)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# CLONAL SHARING
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Heatmap of shared clones  
# ------------------------------------------------------------------------------

# Create a matrix of shared clones between samples
# Get all unique clones per sample
clone_by_sample <- lapply(names(combined.BCR.filtered_GCB), function(sample_name) {
  combined.BCR.filtered_GCB[[sample_name]] %>% 
    pull(CTstrict) %>% 
    unique()
})
names(clone_by_sample) <- names(combined.BCR.filtered_GCB)

# Create pairwise sharing matrix
sample_names <- names(clone_by_sample)
n_samples <- length(sample_names)

sharing_matrix <- matrix(0, nrow = n_samples, ncol = n_samples,
                         dimnames = list(sample_names, sample_names))

for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    shared_clones <- intersect(clone_by_sample[[i]], clone_by_sample[[j]])
    sharing_matrix[i, j] <- length(shared_clones)
  }
}

# Convert to long format for ggplot
sharing_df <- sharing_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("sample1") %>% 
  pivot_longer(-sample1, names_to = "sample2", values_to = "n_shared")

# Plot
for (patient in c("HH117", "HH119")){
  
  # patient <- "HH119"
  
  sharing_df_filtered <- sharing_df %>% 
    filter(
      str_detect(sample1, patient) & str_detect(sample2, patient)
    ) #%>% 
    # filter(sample1 < sample2) # Delete top part of diagonal (including diagonal)
  
  midpoint <- max(sharing_df_filtered$n_shared)/2
  
  sharing_df_filtered %>% 
    ggplot(aes(x = sample2, y = sample1, fill = n_shared)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_shared), size = 3) +
    scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue", 
                         midpoint = midpoint) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(fill = "Shared\nclones",
         title = patient
         ) +
    coord_fixed() 
  
  ggsave(glue("20_VDJ/plot/BCR_clonal_sharing_GCB/BCR_shared_clones_heatmap_{patient}.png"), width = 14, height = 13, dpi = 300)

}

# ------------------------------------------------------------------------------
# Circos plot
# ------------------------------------------------------------------------------

library(circlize)
library(RColorBrewer)

# Get clones per sample
clone_by_sample <- lapply(names(combined.BCR.filtered), function(sample_name) {
  combined.BCR.filtered[[sample_name]] %>% 
    select(CTstrict) %>% 
    distinct() %>%
    mutate(sample = sample_name)
})

# Combine all
all_clones <- bind_rows(clone_by_sample)

# Find clones present in multiple samples
clone_counts <- all_clones %>%
  count(CTstrict, name = "n_samples") %>%
  filter(n_samples > 1)  # Only shared clones

# Get pairwise links
links <- all_clones %>%
  filter(CTstrict %in% clone_counts$CTstrict) %>%
  inner_join(., ., by = "CTstrict", relationship = "many-to-many") %>%
  filter(sample.x < sample.y) %>%  # Avoid duplicates
  count(sample.x, sample.y, name = "n_shared")

for (patient in c("HH117", "HH119")) {
  
  patient <- "HH117"
  
  # Filter links for this patient
  links_patient <- links %>% 
    filter(str_detect(sample.x, patient) & str_detect(sample.y, patient))
  
  # Skip if no links
  if (nrow(links_patient) == 0) {
    message(glue("No shared clones for {patient}, skipping"))
    next
  }
  
  # Set up colors by tissue/condition
  all_samples <- unique(c(links_patient$sample.x, links_patient$sample.y))
  sample_colors <- setNames(
    colorRampPalette(brewer.pal(9, "Set1"))(length(all_samples)),
    all_samples
  )
  
  # Plot
  png(glue("20_VDJ/plot/BCR_clonal_sharing/BCR_shared_clones_circos_{patient}.png"), 
      width = 12, height = 12, units = "in", res = 500)
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.01, 0.01))
  
  chordDiagram(
    links_patient %>% select(sample.x, sample.y, n_shared),
    grid.col = sample_colors,
    transparency = 0.4,
    directional = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(
      track.height = 0.2
    )
  )
  
  # Add sample names with smaller text
  circos.trackPlotRegion(
    track.index = 1, 
    panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(mean(xlim), ylim[1], sector.name, 
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0.5), cex = 0.6)
    }, 
    bg.border = NA
  )
  
  title(glue("B Cell Clonal Sharing Network: {patient}"), cex.main = 1.5)
  
  circos.clear()
  dev.off()
  
}


# ------------------------------------------------------------------------------
# Heatmap of shared clones - Follicles 
# ------------------------------------------------------------------------------

# Needs preprocessing code from other heatmap 

# Plot
for (patient in c("HH117", "HH119")){
  
  # patient <- "HH119"
  
  sharing_df_filtered <- sharing_df %>% 
    filter(
      str_detect(sample1, patient) & str_detect(sample2, patient) &
        str_detect(sample1, "Fol") & str_detect(sample2, "Fol")
    ) %>% 
    filter(sample1 < sample2) # Delete top part of diagonal (including diagonal)
  
  midpoint <- max(sharing_df_filtered$n_shared)/2
  
  sharing_df_filtered %>% 
    ggplot(aes(x = sample2, y = sample1, fill = n_shared)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_shared), size = 3) +
    scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue", 
                         midpoint = midpoint) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(fill = "Shared\nclones",
         title = patient
    ) +
    coord_fixed() 
  
  ggsave(glue("20_VDJ/plot/BCR_clonal_sharing_GCB/BCR_shared_clones_heatmap_{patient}_Fol.png"), width = 10, height = 8, dpi = 300)
  
}

# ------------------------------------------------------------------------------
# Circos plot - muted 
# ------------------------------------------------------------------------------

# Needs preprocessing code from other Circos 

for (patient in c("HH117", "HH119")) {
  
  # patient <- "HH117"
  
  # Filter links for this patient
  links_patient <- links %>% 
    filter(str_detect(sample.x, patient) & str_detect(sample.y, patient))
  
  # Skip if no links
  if (nrow(links_patient) == 0) {
    message(glue("No shared clones for {patient}, skipping"))
    next
  }
  
  # Set up colors by tissue/condition
  all_samples <- unique(c(links_patient$sample.x, links_patient$sample.y))
  sample_colors <- setNames(
    ifelse(str_detect(all_samples, "Fol"), "grey", "darkblue"),
    all_samples
  )
  
  # if (patient == "HH117"){
  #   sample_colors[["HH117-SI-MILF-INF"]] <- "darkblue"
  #   sample_colors[["HH117-SI-MILF-nonINF"]] <- "darkgreen"
  #   sample_colors[["HH117-SILP-INF"]] <- "lightblue"
  #   sample_colors[["HH117-SILP-nonINF"]] <- "lightgreen"
  # } else if (patient == "HH119") {
  #   
  # }

  
  # Plot
  png(glue("20_VDJ/plot/BCR_clonal_sharing/BCR_shared_clones_circos_{patient}_muted.png"), 
      width = 12, height = 12, units = "in", res = 500)
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.01, 0.01))
  
  chordDiagram(
    links_patient %>% select(sample.x, sample.y, n_shared),
    grid.col = sample_colors,
    transparency = 0.4,
    directional = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(
      track.height = 0.2
    )
  )
  
  # Add sample names with smaller text
  circos.trackPlotRegion(
    track.index = 1, 
    panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(mean(xlim), ylim[1], sector.name, 
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0.5), cex = 0.6)
    }, 
    bg.border = NA
  )
  
  title(glue("B Cell Clonal Sharing Network: {patient}"), cex.main = 1.5)
  
  circos.clear()
  dev.off()
  
}

