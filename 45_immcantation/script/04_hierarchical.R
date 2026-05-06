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

bcr_data_qc_annot <- readRDS("45_immcantation/out/rds/03_heavy_bcr_data_qc_annot.rds")

# ------------------------------------------------------------------------------
# Clonal analysis by Hierarchical clustering
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Identify clonal threshold
# Determine clonal clustering threshold: sequences which are under (closer/more similar) this cut-off are clonally related.
# ------------------------------------------------------------------------------

# -------------------
# Automatic: Cross patients 
# -------------------

# Combine BCR data across patients 
bcr_data <- bind_rows(bcr_data_qc_annot)
bcr_data$locus %>% unique()
bcr_data$patient_id %>% table()

# Calculate nearest-neighbor Hamming distance distribution across patients 
dist_crossSubj <- distToNearest(
  bcr_data,
  sequenceColumn = "junction",
  nproc = 1, 
  cross = "patient_id",
  cellIdColumn="cell_id"
)

dist_crossSubj

# -------------------
# Automatic: Per patient
# -------------------

patients <- names(bcr_data_qc_annot)

list_thresholds <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  HH_thresholds <- list()
  
  # -------------------
  # Manual
  # -------------------
  
  bcr_data_HH <- bcr_data_qc_annot[[HH]]
  # bcr_data_HH$locus %>% unique() # Only heavy chain
  
  # Calculate nearest-neighbor Hamming distance distribution
  # Per sequence, what is the Hamming distance to the nearest-neighbor (1 distance per sequence)
  dist_nearest <- distToNearest(
    bcr_data_HH,
    sequenceColumn = "junction",
    cellIdColumn="cell_id" # Important for invoking in run in single-cell mode
  )
  
  # Check how many NAs
  dist_nearest$dist_nearest %>% is.na() %>% table()
  
  # Nearest-neighbor Hamming distance distribution/histogram
  dist_nearest %>% 
    filter(!is.na(dist_nearest)) %>% 
    ggplot(aes(x = dist_nearest)) + 
    geom_histogram(color = "white", binwidth = 0.02) +
    labs(
      title = glue("{HH}: Nearest-neighbor Hamming distance distribution"), 
      x = "Hamming distance", 
      y = "Count"
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    theme_bw()
  
  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram.png"), width = 11, height = 6.5)
  
  # -------------------
  # Automatic - GMM 
  # -------------------
  
  # find threshold for cloning automatically
  threshold_output <- shazam::findThreshold(
    dist_nearest$dist_nearest,
    method = "gmm",
    model = "gamma-norm",
    cutoff = "user",
    spc = 0.995  # specificity
    # spc = 0.99 # slight improve # HH119 - 0.04334775
  )

  threshold <- threshold_output@threshold
  threshold # HH119 - 0.02838713
  HH_thresholds["gmm"] <- threshold

  plot(threshold_output, binwidth = 0.02, silent = TRUE) +
    theme(axis.title = element_text(size = 12)) +
    plot_annotation(
      title = glue("{HH}: Nearest-neighbor Hamming distance distribution"),
      subtitle = glue("GMM, gamma-norm model, specificity = 0.995 threshold: {threshold}")
    )

  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_automatic_gmm_threshold.png"), width = 12, height = 6.5)
  
  # -------------------
  # Automatic - Density
  # -------------------
  
  # find threshold for cloning automatically using density 
  # density finds the valley between the two modes directly
  threshold_output <- shazam::findThreshold(
    dist_nearest$dist_nearest,
    method = "density" 
  )
  
  threshold <- threshold_output@threshold
  threshold # HH119 - 0.1601825
  HH_thresholds["density"] <- threshold
  
  plot(threshold_output, binwidth = 0.02, silent = TRUE) +
    theme(axis.title = element_text(size = 12)) + 
    plot_annotation(
      title = glue("{HH}: Nearest-neighbor Hamming distance distribution"), 
      subtitle = glue("Automatic density threshold: {threshold}")
    )
  
  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_automatic_density_threshold.png"), width = 12, height = 6.5)
  
  # -------------------
  # GMM - cross patient 
  # -------------------
  
  # find threshold for cloning automatically and initialize the Gaussian fit
  # parameters of the nearest-neighbor

  # distance of inter (between) clones using cross subjects distribution of distance to nearest
  threshold_output <- shazam::findThreshold(
    dist_nearest$dist_nearest,
    method = "gmm",
    model = "gamma-norm",
    cross = dist_crossSubj$cross_dist_nearest,
    cutoff = "user",
    spc = 0.995
  )

  threshold_withcross <- threshold_output@threshold
  threshold_withcross
  HH_thresholds["gmm_cross"] <- threshold

  plot(threshold_output, binwidth = 0.02, silent = TRUE) +
    theme(axis.title = element_text(size = 12)) +
    plot_annotation(
      title = glue("{HH} cross patient: Nearest-neighbor Hamming distance distribution"),
      subtitle = glue("GMM, gamma-norm model, specificity = 0.995 threshold: {threshold}")
    )

  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_automatic_gmm_threshold_cross_patient.png"), width = 12, height = 6.5)

  # -------------------
  # Density - cross patient 
  # -------------------
  
  # find threshold for cloning automatically using density 
  # density finds the valley between the two modes directly
  threshold_output <- shazam::findThreshold(
    dist_nearest$dist_nearest,
    method = "density",
    cross = dist_crossSubj$cross_dist_nearest
  )
  
  threshold_withcross <- threshold_output@threshold
  threshold_withcross # 0.1601825
  HH_thresholds["density_cross"] <- threshold
  
  # plot the threshold along the density plot
  plot(threshold_output, binwidth = 0.02, silent = TRUE) +
    theme(axis.title = element_text(size = 12)) + 
    plot_annotation(
      title = glue("{HH} cross patient: Nearest-neighbor Hamming distance distribution"), 
      subtitle = glue("Automatic density threshold: {threshold}")
    )
  
  ggsave(glue("45_immcantation/plot/hier_distance_distributions/{HH}_nearest_neighbor_Hamming_distance_histogram_automatic_density_threshold_cross_patient.png"), width = 12, height = 6.5)
  
  return(HH_thresholds)
  
}) %>% setNames(patients)

list_thresholds
saveRDS(list_thresholds, "45_immcantation/out/rds/04_list_thresholds.rds")

# ------------------------------------------------------------------------------
# Define clonal groups
# ------------------------------------------------------------------------------

# Mean threshold for all patients
# threshold <- list_thresholds %>% unlist() %>% mean()

# Patient-specific threshold
hier <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  bcr_data_HH <- bcr_data_qc_annot[[HH]]
  
  hier_HH <- hierarchicalClones(
    bcr_data_HH,
    cell_id = "cell_id",
    threshold = list_thresholds[[HH]]$density,
    only_heavy = TRUE,
    split_light = FALSE,
    summarize_clones = FALSE
  )
  
  return(hier_HH)

}) %>% setNames(patients)

# Add metadata
hier <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  hier_HH <- hier[[HH]]
  
  hier_HH <- hier_HH %>%
    mutate(
      sample_clean_fol = ifelse(!is.na(manual_ADT_ID), paste(sample_clean, manual_ADT_ID, sep = "_"), sample_clean)
    )
  
  return(hier_HH)
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Clonal overview
# ------------------------------------------------------------------------------

# Inspect V and J call across samples and clones 
# IGHJ4*01 and IGHJ4*02 are the same gene from different alleles 
# Some genes have 1 alleles varies and some have many. They can be looked up in a database. 

# Allelic variants are of the same gene at the same genomic locus.
# They differ only by a few nucleotides (SNPs, rarely small indels)

# If a sequence aligns equally well to multiple allelic variations, they are all stated seperated by commas. 

# -------------------
# HH119
# -------------------

HH <- "HH119"
hier[[HH]] %>% 
  count(clone_id, sort = TRUE)
# count(clone_id, v_call, j_call, sort = TRUE)

top_clone <- hier[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

hier[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean_fol, v_call, j_call, sort = TRUE)

# Several J genes in one clone??? Potential reasons:  
# misassignments
# IGHJ4 --> extensive SHM --> sequence looks more like IGHJ5 now

# How large is the problem? - not that big :)
hier[[HH]] %>%
  filter(clone_id == top_clone) %>%
  mutate(j_gene = sub("\\*.*", "", j_call)) %>%  # strip allele, keep gene
  count(j_gene) %>%
  mutate(pct = n/sum(n)*100)

# Same genes in other clones
# Good, means that CDR3/junction region is being taken into account
hier[[HH]] %>% 
  filter(str_detect(v_call, "IGHV4-34") & str_detect(j_call, "IGHJ4")) %>% 
  count(clone_id, sort = TRUE)

hier[[HH]] %>% 
  filter(clone_id == 34978) %>% 
  count(v_call, j_call, sort = TRUE)

# -------------------
# HH117
# -------------------

HH <- "HH117"
hier[[HH]] %>% 
  count(clone_id, v_call, j_call, sort = TRUE)

top_clone <- hier[[HH]] %>% count(clone_id, v_call, j_call, sort = TRUE) %>% 
  head(1) %>% pull(clone_id)

hier[[HH]] %>% 
  filter(clone_id == top_clone) %>% 
  count(sample_clean_fol, v_call, j_call, sort = TRUE)

# ------------------------------------------------------------------------------
# Define top clones
# ------------------------------------------------------------------------------

top_GC_clones_hier <- lapply(patients, function(HH) {
  
  # find clones that contain at least 1 GC cell
  GC_clones <- hier[[HH]] %>%
    filter(celltype_broad == "GC_B_cells") %>%
    pull(clone_id) %>%
    unique()
  
  # rank those clones by total size (all cell types) and take top 10
  hier[[HH]] %>%
    filter(clone_id %in% GC_clones) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = 10) %>%
    pull(clone_id)
  
}) %>% setNames(patients)

# -------------------
# Look at top clones novj
# -------------------

lapply(patients, function(HH){
  
  # HH <- "HH119"
  # HH <- "HH117"
  
  hier[[HH]] %>%
    filter(clone_id %in% top_GC_clones_hier[[HH]]) %>% 
    mutate(
      v_gene = v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = v_call
    ) %>% 
    mutate(
      j_gene = j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", "), .by = j_call
    ) %>% 
    count(clone_id, v_gene, j_gene, sort = TRUE) 
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Visualize top clones
# ------------------------------------------------------------------------------

source("10_broad_annotation/script/color_palette.R")

for (HH in patients){
  
  # HH <- "HH119"
  HH_top_clones <- top_GC_clones_hier[[HH]]
  
  for (clone_nr in 1:length(HH_top_clones)){
    
    # clone_nr <- 1
    clone <- HH_top_clones[clone_nr]
    
    plot_df <- hier[[HH]] %>% 
      filter(clone_id == clone) %>% 
      mutate(sample_clean_fol = fct_infreq(sample_clean_fol) %>% fct_rev()) %>%
      add_count(sample_clean_fol, name = "Count")
    
    # Meta data
    n_cells <- plot_df %>% nrow()
    v_gene <- plot_df$v_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    j_gene <- plot_df$j_call %>% str_split_i("\\*", 1) %>% unique() %>% paste0(collapse = ", ")
    
    plot_df %>%   
      ggplot(aes(y = sample_clean_fol, fill = celltype_broad)) +
      geom_bar() + 
      scale_fill_manual(values = celltype_colors) + 
      geom_text(aes(x = Count, label = Count), 
                hjust = -0.2, color = "black",
                stat = "unique") +
      labs(
        title = glue("{HH}: Top {clone_nr} GCB clone"),
        subtitle = glue("Clone ID: {clone}"),
        caption = glue("N cells: {n_cells}\nV gene: {v_gene}\n J gene: {j_gene}"),
        y = ""
      ) +
      theme_bw()
    
    ggsave(glue("45_immcantation/plot/GC_clones_hier/{HH}_clone_nr_{clone_nr}_across_samples_and_cell_types.png"), width = 15, height = 8.5)
    
  }
}

# ------------------------------------------------------------------------------
# Export hier
# ------------------------------------------------------------------------------

saveRDS(hier, "45_immcantation/out/rds/hier_clones.rds")

# ------------------------------------------------------------------------------
# Visualize clonal abundance - I think we have too many clones for these plots to make sense
# ------------------------------------------------------------------------------
# 
# # calculate and plot the rank-abundance curve
# abund <- estimateAbundance(
#   hier,
#   group = "sample_id", 
#   nboot = 100
# )
# 
# abund_plot <- plot(abund, silent=T)
# abund_plot
# 
# abund_plot + facet_wrap("sample_id", ncol = 3)
# 
# # get clone sizes using dplyr functions
# clone_sizes <- countClones(
#   hier,
#   groups = "sample_clean"
# )
# 
# clone_sizes$sample_clean %>% unique()
# 
# # plot cells per clone
# clone_sizes %>% 
#   filter(sample_clean == "HH117-SI-PP-nonINF") %>%
#   # filter(sample_clean == "HH117-SILP-INF") %>%
#   # filter(sample_clean == "HH119-SI-PP") %>%
#   ggplot(aes(x = seq_count)) +
#   geom_bar() +
#   # facet_wrap("sample_clean", ncol = 3) +
#   labs(x = "Sequences per clone") +
#   theme_bw()
# 
# # calculate and plot the rank-abundance curve
# div <- alphaDiversity(
#   hier,
#   group = "sample_clean", 
#   nboot = 100
# )
# plot(div, silent = TRUE) + facet_wrap("sample_clean", ncol = 3)

# ------------------------------------------------------------------------------
# Create germlines
# ------------------------------------------------------------------------------

# In terminal 
# cd /home/people/helweg/ciir/people/helweg/projects/packages
# git clone https://bitbucket.org/kleinstein/immcantation.git
# ./immcantation/scripts/fetch_imgtdb.sh 

# read in IMGT data if downloaded on your own (above)
# references <- readIMGT(dir = "/home/people/helweg/ciir/people/helweg/projects/packages/human/vdj")
# readIMGT function does not exist 

# Claude suggested this instead:

# vdj_files <- list.files("/home/people/helweg/ciir/people/helweg/projects/packages/human/vdj", 
#                         full.names = TRUE, pattern = "\\.fasta$")
# 
# library(Biostrings)
# references <- lapply(vdj_files, readDNAStringSet)
# names(references) <- gsub("\\.fasta$", "", basename(vdj_files))
# 
# ## [1] "Read in 1243 from 17 fasta files"
# 
# # reconstruct germlines
# # createGermlines function also does not exist...
# results <- createGermlines(
#   hier, 
#   references, 
#   fields = "patient",
#   nproc = 1
# )

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------






