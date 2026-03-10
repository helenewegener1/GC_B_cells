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
library(Biostrings) # writing fasta files 
library(msa) # mulitple sequence alignment 

# celltype_colors
source("10_broad_annotation/script/color_palette.R")

combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")

# ------------------------------------------------------------------------------
# Clonal abundance per sample
# ------------------------------------------------------------------------------

# Count of clone sizes
clone_abundance <- lapply(combined.BCR.filtered, function(x){
  
  x %>%
    filter(!str_detect(manual_cluster, "NA")) %>% 
    summarise(abundance = n(), .by = c(manual_cluster, sample)) %>% 
    arrange(desc(abundance))
  
})

# Check numbers
## N cells
combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF_Fol-1") %>% nrow() 
## N unique clones
combined.BCR.filtered$HH117 %>% filter(sample == "HH117-SI-PP-nonINF_Fol-1" & !str_detect(manual_cluster, "NA") & !is.na(manual_cluster)) %>% nrow()
## N unique clones
clone_abundance$HH117 %>% filter(sample == "HH117-SI-PP-nonINF_Fol-1") %>% pull(abundance) %>% sum()

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
# Number of unique CT per sample
# ------------------------------------------------------------------------------

combined.BCR.filtered_all <- combined.BCR.filtered %>%
  bind_rows() %>% 
  filter(!str_detect(manual_cluster, "NA") & !is.na(manual_cluster))

n_unique_CT <- combined.BCR.filtered_all %>% 
  summarize(n = n(), .by = c("sample", "manual_cluster")) %>% 
  mutate(total_cells = sum(n), .by = "sample")

# Plot 1: Bar plot of unique clones per sample
ggplot() +
  geom_col(
    data = n_unique_CT, 
    aes(x = sample, y = n, fill = manual_cluster)
  ) +
  geom_text(
    data = n_unique_CT %>% select(sample, total_cells) %>% distinct(), 
    aes(x = sample, y = total_cells, label = total_cells), 
    hjust = -0.1, 
    size = 3 
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Number of Unique CT per Sample",
    caption = "Colored by CT", 
    x = "",
    y = "Number of unique CT"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )

ggsave("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_final_N_unique_CTs.png", 
       width = 14, height = 15, dpi = 300)

# ------------------------------------------------------------------------------
# Number of unique GCB CT per sample
# ------------------------------------------------------------------------------

GCB_CT_clones <- combined.BCR.filtered_all %>% 
  filter(celltype_broad == "GC_B_cells") %>% 
  select(manual_cluster) %>% 
  distinct() %>% 
  pull()
  
n_unique_GCB_CT <- combined.BCR.filtered_all %>% 
  filter(manual_cluster %in% GCB_CT_clones) %>% 
  summarize(n = n(), .by = c("sample", "manual_cluster")) %>% 
  mutate(total_cells = sum(n), .by = "sample")

# Plot 1: Bar plot of unique clones per sample
ggplot() +
  geom_col(
    data = n_unique_GCB_CT, 
    aes(x = sample, y = n, fill = manual_cluster)
  ) +
  geom_text(
    data = n_unique_GCB_CT %>% select(sample, total_cells) %>% distinct(), 
    aes(x = sample, y = total_cells, label = total_cells), 
    hjust = -0.1, 
    size = 3 
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Number of Unique GCB CT per Sample",
    caption = "Colored by CT", 
    x = "",
    y = "Number of unique GCB CT"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )

ggsave("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_final_N_unique_GCB_CTs.png", 
       width = 14, height = 15, dpi = 300)

#BCR_final_N_cells.png", width = 14, height = 15) 

# ------------------------------------------------------------------------------
# Number of cells colored by CTstrict 
# ------------------------------------------------------------------------------

for (sample_name in c("HH117-SI-PP-nonINF", "HH119-SI-PP")){

  # sample_name <- "HH117-SI-PP-nonINF"
  # sample_name <- "HH119-SI-PP"
  
  # Subset sample
  combined.BCR.filtered_subset <- combined.BCR.filtered_all %>% filter(sample_clean == sample_name)
  
  combined.BCR.filtered_subset <- combined.BCR.filtered_subset %>% 
    mutate(CT_abundance = n(), .by = c(sample, manual_cluster)) %>% 
    mutate(
      Fol = str_split_i(sample, "_", 2),
      CT_above_2 = ifelse(CT_abundance >= 2, manual_cluster, NA),
      CT_above_5 = ifelse(CT_abundance >= 5, manual_cluster, NA)
    )
  
  combined.BCR.filtered_subset %>% 
    ggplot(aes(
      x = Fol,
      # fill = CT_above_2
      fill = CT_above_5
    )) +
    geom_bar() + 
    labs(
      title = "N cells colored by CT",
      subtitle = sample_name,
      x = "",
      y = "N cells", 
      caption = "CT with abundance >= 5 are colored"
    ) + 
    theme_bw() + 
    theme(legend.position = "none")
  
  ggsave(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_N_cells_CT_{sample_name}.png"), 
         width = 14, height = 7, dpi = 300)
  
}

# ------------------------------------------------------------------------------
# Clonal abundance per sample per broad cell type 
# ------------------------------------------------------------------------------

# Count of clone sizes
clone_abundance_celltype <- combined.BCR.filtered_all %>% 
  summarise(abundance = n(), .by = c(sample, manual_cluster, celltype_broad)) %>% 
  mutate(abundance_total = sum(abundance), .by = manual_cluster) %>% 
  arrange(desc(abundance_total), desc(abundance))

# Get top clones
patients <- names(combined.BCR.filtered)
top_GC_clones <- lapply(patients, function(x){
  
  combined.BCR.filtered[[x]] %>%
    filter(celltype_broad == "GC_B_cells" & !is.na(manual_cluster) & !str_detect(manual_cluster, "NA")) %>% 
    summarise(n = n(), .by = "manual_cluster") %>% 
    arrange(desc(n)) %>% 
    head(n = 10) %>% 
    pull(manual_cluster)
  
}) %>% 
  setNames(patients)

# Plot
for (HH in patients){
  
  for (CT_nr in 1:length(top_GC_clones[[HH]])){
    
    # HH <- "HH119"
    # CT_nr <- 1
    CT <- top_GC_clones[[HH]][[CT_nr]]
    
    # ------------------------------------------------------------------------------
    # Bar plot of clonal abundance per sample per broad cell type 
    # ------------------------------------------------------------------------------
 
    # Filter data
    plot_data <- clone_abundance_celltype %>% 
      filter(manual_cluster == CT) %>% 
      mutate(n_cells_CT_per_sample = sum(abundance), .by = "sample")
    
    df_n_cells_CT_per_sample <- plot_data %>% 
      select(sample, n_cells_CT_per_sample) %>% 
      distinct()
    
    n_cells_total <- plot_data$abundance_total %>% unique()
    genes <- combined.BCR.filtered[[HH]] %>% filter(manual_cluster == CT) %>% pull(gene_cluster) %>% unique()
    
    # Create plot
    plot_data %>% 
      ggplot(aes(x = abundance, y = reorder(sample, n_cells_CT_per_sample), fill = celltype_broad)) + 
      geom_col() + 
      geom_text(
        data = df_n_cells_CT_per_sample,
        aes(x = n_cells_CT_per_sample, y = sample, label = n_cells_CT_per_sample),  
        hjust = -0.1, 
        size = 3, 
        inherit.aes = FALSE
      ) + 
      scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
      scale_fill_manual(values = celltype_colors) + 
      labs(
        title = "Clonal abundance by sample and cell type",
        subtitle = glue("{HH}: {CT}\nGenes: {genes}"), 
        caption = glue("N cells total: {n_cells_total}\nCT number: {CT_nr}"),
        x = "Clone size (total cells)",
        y = "",
        fill = "Cell type"
      ) + 
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 8),
        legend.position = "right"
      )
    
    ggsave(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/CT_individual/{HH}_CT{CT_nr}_abundance.png"), width = 20, height = 10)
    
    # ------------------------------------------------------------------------------
    # Heatmap of sample and isotype of one CT
    # ------------------------------------------------------------------------------
    
    combined.BCR.filtered[[HH]] %>% 
      filter(manual_cluster == CT) %>% 
      group_by(sample, Ig_class) %>%
      summarise(Count = n(), .groups = "drop") %>%
      complete(sample, Ig_class, fill = list(Count = 0)) %>%
      ggplot(aes(x = sample, y = Ig_class, fill = Count)) + 
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(Count > 0, Count, "0")), 
                color = "white", size = 2.5) +
      scale_fill_gradient(low = "#c8d8e8", high = "#0d2a4e",
                          limits = c(0, NA)) +
      labs(
        x = NULL, 
        y = "Ig Class", 
        fill = "Count",  
        subtitle = glue("{HH}: {CT}\nGenes: {genes}"),
        caption = glue("N cells total: {n_cells_total}\nCT number: {CT_nr}")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    ggsave(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/CT_individual/{HH}_CT{CT_nr}_Igclass_vs_sample.png"), width = 12, height = 6)
    
  }
  
  # ------------------------------------------------------------------------------
  # N clones per CT in fasta file 
  # ------------------------------------------------------------------------------

  # HH <- "HH119"
  # CT_nr <- 1
  # CT <- top_GC_clones[[HH]][[CT_nr]]
  
  clone_df <- combined.BCR.filtered[[HH]] %>% 
    filter(manual_cluster == CT) %>% 
    summarise(abundance = n(), .by = "IGH_full_sequence") %>% 
    arrange(desc(abundance)) 
  
  # N unique clones
  # nrow(clone_df)
  
  # Get sequences and name with abundance
  seqs <- DNAStringSet(clone_df$IGH_full_sequence)
  names(seqs) <- clone_df$abundance
  
  # Add germline to rest 
  # seqs_final <- c(seqs, germline_dna)
  
  # Align sequence (msa package works on computerome)
  ## Run alignment (ClustalW, ClustalOmega, or MUSCLE available)
  # alignment <- msa(seqs_final, method = "ClustalOmega")
  
  ## Convert to DNAStringSet for downstream use
  # alignment_dna <- as(alignment, "DNAStringSet")
  
  # Export as FASTA file
  # writeXStringSet(alignment_dna, filepath = glue("20_VDJ/fasta/{HH}_clone_{clone_nr}_aligned_dna.fasta"))
  writeXStringSet(seqs, filepath = glue("20_VDJ/fasta/{HH}_clone_{CT_nr}_{CT}_rna.fasta"))
  
}

# ------------------------------------------------------------------------------
# Distance between sequences in CT
# ------------------------------------------------------------------------------

library(stringdist)

HH <- "HH117"

# Calculate max pairwise Levenshtein distance within each manual cluster
cluster_distances <- combined.BCR.filtered[[HH]] %>%
  filter(manual_cluster %in% top_GC_clones[[HH]]) %>%
  group_by(manual_cluster) %>%
  filter(n() > 1) %>%
  summarize(
    n_cells = n(),
    n_unique_seqs = n_distinct(cdr3_nt1, na.rm = TRUE),
    max_dist = {
      seqs <- unique(cdr3_nt1[!is.na(cdr3_nt1)])
      if (length(seqs) > 1) {
        dm <- stringdistmatrix(seqs, seqs, method = "lv")
        lens <- nchar(seqs)
        max_len <- outer(lens, lens, pmax)
        max(dm / max_len)
      } else NA_real_
    },
    similarity = 1-max_dist,
    .groups = "drop"
  ) %>%
  arrange(desc(max_dist))

cluster_distances

# Inspect
combined.BCR.filtered[[HH]] %>%
  filter(manual_cluster == "cluster.2231_cluster.474")

# ------------------------------------------------------------------------------
# TFH with BCR 
# ------------------------------------------------------------------------------

clone_abundance_celltype %>% filter(celltype_broad == "Tfh_like_cells") %>% data.frame()

# ------------------------------------------------------------------------------
# Clonal abundance per sample
# ------------------------------------------------------------------------------

HHs <- combined.BCR.filtered_all$patient %>% unique()

lapply(HHs, function(HH){
  
  # HH <- "HH119_Control"
  
  # Count of clone sizes
  HH_clone_abundance_celltype <- combined.BCR.filtered_all %>% 
    filter(patient == HH) %>% 
    summarise(abundance = n(), .by = c(manual_cluster)) %>% 
    arrange(desc(abundance))
  
  HH_clone_abundance_celltype %>% 
    filter(abundance > 31) %>% 
    ggplot(aes(x = abundance, y = reorder(manual_cluster, abundance))) + 
    geom_col() + 
    geom_text(
      aes(x = abundance, y = manual_cluster, label = abundance),  
      hjust = -0.1, 
      size = 3, 
      inherit.aes = FALSE
    ) + 
    labs(
      title = HH,
      y = "", 
      x = "N cells", 
      caption = "Only CTs with N cells > 20 is shown"
    ) + 
    theme_bw() 
  
  ggsave(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_N_cells_in_CT_{HH}.png"), width = 10, height = 10)
  
  
})


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# CLONAL SHARING
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# # ------------------------------------------------------------------------------
# # Heatmap of shared clones  
# # ------------------------------------------------------------------------------
# 
# # Create a matrix of shared clones between samples
# # Get all unique clones per sample
# clone_by_sample <- lapply(names(combined.BCR.filtered), function(HH) {
#   combined.BCR.filtered[[HH]] %>% 
#     pull(manual_cluster) %>% 
#     unique()
# }) %>% 
#   setNames(names(combined.BCR.filtered))
# 
# HH <- "HH117"
# 
# # Create pairwise sharing matrix
# sample_names <- names(clone_by_sample)
# n_samples <- length(sample_names)
# 
# sharing_matrix <- matrix(0, nrow = n_samples, ncol = n_samples,
#                          dimnames = list(sample_names, sample_names))
# 
# for (i in 1:n_samples) {
#   for (j in 1:n_samples) {
#     shared_clones <- intersect(clone_by_sample[[i]], clone_by_sample[[j]])
#     sharing_matrix[i, j] <- length(shared_clones)
#   }
# }
# 
# # Convert to long format for ggplot
# sharing_df <- sharing_matrix %>% 
#   as.data.frame() %>% 
#   rownames_to_column("sample1") %>% 
#   pivot_longer(-sample1, names_to = "sample2", values_to = "n_shared")
# 
# # Plot
# for (patient in c("HH117", "HH119")){
#   
#   # patient <- "HH117"
#   
#   sharing_df_filtered <- sharing_df %>% 
#     filter(
#       str_detect(sample1, patient) & str_detect(sample2, patient)
#     ) #%>% 
#     # filter(sample1 < sample2) # Delete top part of diagonal (including diagonal)
#   
#   midpoint <- max(sharing_df_filtered$n_shared)/2
#   
#   sharing_df_filtered %>% 
#     ggplot(aes(x = sample2, y = sample1, fill = n_shared)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = n_shared), size = 3) +
#     scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue", 
#                          midpoint = midpoint) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#       axis.text.y = element_text(size = 10),
#       axis.title = element_blank(),
#       panel.grid = element_blank()
#     ) +
#     labs(fill = "Shared\nclones",
#          title = patient
#          ) +
#     coord_fixed() 
#   
#   ggsave(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_shared_clones_heatmap_{patient}.png"), width = 14, height = 13, dpi = 300)
# 
# }
# 
# # ------------------------------------------------------------------------------
# # Circos plot
# # ------------------------------------------------------------------------------
# 
# library(circlize)
# library(RColorBrewer)
# 
# # Get clones per sample
# clone_by_sample <- lapply(names(combined.BCR.filtered), function(sample_name) {
#   combined.BCR.filtered[[sample_name]] %>% 
#     select(CTstrict) %>% 
#     distinct() %>%
#     mutate(sample = sample_name)
# })
# 
# # Combine all
# all_clones <- bind_rows(clone_by_sample)
# 
# # Find clones present in multiple samples
# clone_counts <- all_clones %>%
#   count(CTstrict, name = "n_samples") %>%
#   filter(n_samples > 1)  # Only shared clones
# 
# # Get pairwise links
# links <- all_clones %>%
#   filter(CTstrict %in% clone_counts$CTstrict) %>%
#   inner_join(., ., by = "CTstrict", relationship = "many-to-many") %>%
#   filter(sample.x < sample.y) %>%  # Avoid duplicates
#   count(sample.x, sample.y, name = "n_shared")
# 
# for (patient in c("HH117", "HH119")) {
#   
#   patient <- "HH117"
#   
#   # Filter links for this patient
#   links_patient <- links %>% 
#     filter(str_detect(sample.x, patient) & str_detect(sample.y, patient))
#   
#   # Skip if no links
#   if (nrow(links_patient) == 0) {
#     message(glue("No shared clones for {patient}, skipping"))
#     next
#   }
#   
#   # Set up colors by tissue/condition
#   all_samples <- unique(c(links_patient$sample.x, links_patient$sample.y))
#   sample_colors <- setNames(
#     colorRampPalette(brewer.pal(9, "Set1"))(length(all_samples)),
#     all_samples
#   )
#   
#   # Plot
#   png(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_shared_clones_circos_{patient}.png"), 
#       width = 12, height = 12, units = "in", res = 500)
#   
#   circos.clear()
#   circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.01, 0.01))
#   
#   chordDiagram(
#     links_patient %>% select(sample.x, sample.y, n_shared),
#     grid.col = sample_colors,
#     transparency = 0.4,
#     directional = 0,
#     annotationTrack = "grid",
#     preAllocateTracks = list(
#       track.height = 0.2
#     )
#   )
#   
#   # Add sample names with smaller text
#   circos.trackPlotRegion(
#     track.index = 1, 
#     panel.fun = function(x, y) {
#       xlim = get.cell.meta.data("xlim")
#       ylim = get.cell.meta.data("ylim")
#       sector.name = get.cell.meta.data("sector.index")
#       circos.text(mean(xlim), ylim[1], sector.name, 
#                   facing = "clockwise", niceFacing = TRUE,
#                   adj = c(0, 0.5), cex = 0.6)
#     }, 
#     bg.border = NA
#   )
#   
#   title(glue("B Cell Clonal Sharing Network: {patient}"), cex.main = 1.5)
#   
#   circos.clear()
#   dev.off()
#   
# }
# 
# 
# # ------------------------------------------------------------------------------
# # Heatmap of shared clones - Follicles 
# # ------------------------------------------------------------------------------
# 
# # Needs preprocessing code from other heatmap 
# 
# # Plot
# for (patient in c("HH117", "HH119")){
#   
#   # patient <- "HH119"
#   
#   sharing_df_filtered <- sharing_df %>% 
#     filter(
#       str_detect(sample1, patient) & str_detect(sample2, patient) &
#         str_detect(sample1, "Fol") & str_detect(sample2, "Fol")
#     ) %>% 
#     filter(sample1 < sample2) # Delete top part of diagonal (including diagonal)
#   
#   midpoint <- max(sharing_df_filtered$n_shared)/2
#   
#   sharing_df_filtered %>% 
#     ggplot(aes(x = sample2, y = sample1, fill = n_shared)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = n_shared), size = 3) +
#     scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue", 
#                          midpoint = midpoint) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#       axis.text.y = element_text(size = 10),
#       axis.title = element_blank(),
#       panel.grid = element_blank()
#     ) +
#     labs(fill = "Shared\nclones",
#          title = patient
#     ) +
#     coord_fixed() 
#   
#   ggsave(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_shared_clones_heatmap_{patient}_Fol.png"), width = 10, height = 8, dpi = 300)
#   
# }
# 
# # ------------------------------------------------------------------------------
# # Circos plot - muted 
# # ------------------------------------------------------------------------------
# 
# # Needs preprocessing code from other Circos 
# 
# for (patient in c("HH117", "HH119")) {
#   
#   # patient <- "HH117"
#   
#   # Filter links for this patient
#   links_patient <- links %>% 
#     filter(str_detect(sample.x, patient) & str_detect(sample.y, patient))
#   
#   # Skip if no links
#   if (nrow(links_patient) == 0) {
#     message(glue("No shared clones for {patient}, skipping"))
#     next
#   }
#   
#   # Set up colors by tissue/condition
#   all_samples <- unique(c(links_patient$sample.x, links_patient$sample.y))
#   sample_colors <- setNames(
#     ifelse(str_detect(all_samples, "Fol"), "grey", "darkblue"),
#     all_samples
#   )
#   
#   # if (patient == "HH117"){
#   #   sample_colors[["HH117-SI-MILF-INF"]] <- "darkblue"
#   #   sample_colors[["HH117-SI-MILF-nonINF"]] <- "darkgreen"
#   #   sample_colors[["HH117-SILP-INF"]] <- "lightblue"
#   #   sample_colors[["HH117-SILP-nonINF"]] <- "lightgreen"
#   # } else if (patient == "HH119") {
#   #   
#   # }
# 
#   
#   # Plot
#   png(glue("20_VDJ/plot/BCR_02_clonal_abundance_sharing/BCR_shared_clones_circos_{patient}_muted.png"), 
#       width = 12, height = 12, units = "in", res = 500)
#   
#   circos.clear()
#   circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.01, 0.01))
#   
#   chordDiagram(
#     links_patient %>% select(sample.x, sample.y, n_shared),
#     grid.col = sample_colors,
#     transparency = 0.4,
#     directional = 0,
#     annotationTrack = "grid",
#     preAllocateTracks = list(
#       track.height = 0.2
#     )
#   )
#   
#   # Add sample names with smaller text
#   circos.trackPlotRegion(
#     track.index = 1, 
#     panel.fun = function(x, y) {
#       xlim = get.cell.meta.data("xlim")
#       ylim = get.cell.meta.data("ylim")
#       sector.name = get.cell.meta.data("sector.index")
#       circos.text(mean(xlim), ylim[1], sector.name, 
#                   facing = "clockwise", niceFacing = TRUE,
#                   adj = c(0, 0.5), cex = 0.6)
#     }, 
#     bg.border = NA
#   )
#   
#   title(glue("B Cell Clonal Sharing Network: {patient}"), cex.main = 1.5)
#   
#   circos.clear()
#   dev.off()
#   
# }

