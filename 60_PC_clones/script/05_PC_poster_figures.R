library(glue)
library(tidyverse)
library(UpSetR)
library(grid)
library(scatterpie)
source("10_broad_annotation/script/color_palette.R")

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

df_both <- readRDS("45_immcantation/out/rds/resolve_LC_90_similarity_germlined.rds")
patients <- names(df_both)

# Look at clone IDs
grep("clone", colnames(df_both$HH117), value = TRUE)
# clone_subgroup_id_90_similarity

# Load seurat object
seurat_integrated <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")

outdir <- "60_PC_clones/plot/PC_poster_figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

get_majority <- function(calls) {
  genes <- unlist(strsplit(calls, ","))
  tab <- table(genes)
  max_count <- max(tab)
  paste(names(tab[tab == max_count]), collapse = ",")
}

df_both <- lapply(patients, function(HH){
  
  # HH <- "HH119"
  df_both[[HH]] %>%
    group_by(clone_subgroup_id_90_similarity, locus) %>%
    mutate(
      v_call_majority = get_majority(v_call),
      j_call_majority = get_majority(j_call)
    ) %>%
    ungroup()
  
}) %>% setNames(patients)

# ------------------------------------------------------------------------------
# Define top PC clones 
# ------------------------------------------------------------------------------

PC_clones <- list()

for (HH in patients){
  
  # HH <- "HH117"
  
  df_HH <- df_both[[HH]]
  
  LP_sites <- grep("LP", unique(df_HH$sample_clean_fol), value = TRUE)
  
  for (site in LP_sites){
    
    site_clones <- df_HH %>% 
      filter(
        locus == "IGH", 
        L1_annotation == "PCs", 
        sample_clean_fol == sym(site)
      ) %>% 
      dplyr::count(clone_subgroup_id_90_similarity, sort = TRUE) %>% 
      head(20) %>% 
      pull(clone_subgroup_id_90_similarity)
    
    PC_clones[[site]] <- site_clones
    
  }
  
}

PC_clones

# ------------------------------------------------------------------------------
# Different isotypes in different follicles and sites? (All cells )
# ------------------------------------------------------------------------------

outdir_1 <- glue("{outdir}/clonal_sharing")
dir.create(outdir_1, recursive = TRUE)

for (site in names(PC_clones)){
  
  # site <- "HH117-SILP-nonINF"
  HH <- str_split_i(site, "-", 1)
  
  df_plot <- df_both[[HH]] %>% 
    filter(
      L1_annotation == "PCs", 
      clone_subgroup_id_90_similarity %in% PC_clones[[site]]
    )
  
  # ---- build x/y lookups so pies have numeric coordinates ----
  
  x_lookup <- df_plot %>% 
    distinct(sample_clean) %>% 
    arrange(sample_clean) %>% 
    mutate(x = row_number())
  
  y_lookup <- df_plot %>% 
    distinct(clone_subgroup_id_90_similarity) %>% 
    arrange(clone_subgroup_id_90_similarity) %>% 
    mutate(y = row_number())
  
  # ---- counts per clone x sample x isotype, pivoted wide for scatterpie ----
  
  pie_data <- df_plot %>% 
    count(sample_clean, clone_subgroup_id_90_similarity, c_call) %>% 
    left_join(x_lookup, by = "sample_clean") %>% 
    left_join(y_lookup, by = "clone_subgroup_id_90_similarity") %>% 
    pivot_wider(
      id_cols = c(sample_clean, clone_subgroup_id_90_similarity, x, y),
      names_from = c_call, 
      values_from = n, 
      values_fill = 0
    )
  
  isotype_cols <- setdiff(names(pie_data), c("sample_clean", "clone_subgroup_id_90_similarity", "x", "y"))
  
  pie_data <- pie_data %>% 
    mutate(total_n = rowSums(across(all_of(isotype_cols))))
  
  # ---- plot ----
  ggplot() + 
    geom_scatterpie(
      data = pie_data, 
      aes(x = x, y = y, r = 0.35),
      cols = isotype_cols, 
      color = NA
    ) + 
    geom_text(
      data = pie_data,
      aes(x = x, y = y, label = total_n),
      size = 2.5,
      nudge_y = -0.5
    ) +
    coord_equal() + 
    scale_fill_manual(values = isotype_colors_custom, na.value = "grey80") + 
    scale_x_continuous(breaks = x_lookup$x, labels = x_lookup$sample_clean, minor_breaks = scales::breaks_width(1)) +
    scale_y_continuous(breaks = y_lookup$y, labels = y_lookup$clone_subgroup_id_90_similarity, minor_breaks = scales::breaks_width(1)) +
    coord_equal() + 
    labs(x = "Sample", y = "Clone", fill = "Isotype", title = glue("{site} largest clones: PCs only")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(glue("{outdir_1}/{site}_isotype_follicle_shared_clones.png"), width = 5.5, height = 8)
  
}


# ==============================================================================
# Frequency of top clone per follicle 
# ==============================================================================

outdir_6 <- glue("{outdir}/PC_clones_freq_barplot")
dir.create(outdir_6, recursive = TRUE)

n_clones <- 10

# Visualize PC top clones
for (site in names(PC_clones)){
  
  # site <- "HH117-SILP-INF"
  HH <- site %>% str_split_i("-", 1)
  df_HH <- df_both[[HH]]
  
  # Subset clones
  top_PC_clones <- PC_clones[[site]][c(1:n_clones)]

  plot_df <- df_both[[HH]] %>% 
    filter(
      locus == "IGH", 
      L1_annotation == "PCs",
      str_detect(sample_clean_fol, "LP")
    ) %>% 
    mutate(
      clone_subgroup_id_90_similarity_plot = ifelse(clone_subgroup_id_90_similarity %in% top_PC_clones, clone_subgroup_id_90_similarity, "other"),
      clone_subgroup_id_90_similarity_plot = factor(clone_subgroup_id_90_similarity_plot, levels = c(top_PC_clones, "other"))
    ) %>%
    add_count(sample_clean_fol, name = "Count") 

  # Define clone colors 
  clone_colors <- list(
    # "#E05C8A", "#66CC55", "#5588DD", "#EE9944", "#AA3377",
    # "#44BBAA", , "#4499CC", "#AACC33", "#9955BB",
    "#FF0000", "#0000FF", "#CC6644", "#FF6600", "#9900CC",
    "#00CCCC", "#FF0099", "#996600", "#0099FF", "#669900",
    "grey85"
  ) %>% setNames(c(top_PC_clones, "other"))
  
  # Define clone names
  clone_names <- c(paste("Clone", 1:n_clones), "Other") %>% as.list() %>% setNames(c(top_PC_clones, "other"))
  
  # N clones 
  N_clones_per_fol <- plot_df %>% 
    group_by(sample_clean_fol) %>%
    count(clone_subgroup_id_90_similarity) %>%
    count(sample_clean_fol) %>%
    ungroup()
  
  if (HH == "HH117"){
    width <- 12 
  } else if (HH == "HH119"){
    width <- 15.5
  }
  
  png(glue("{outdir_6}/{site}_N_{n_clones}.png"), width = width, height = 7, units = "in", res = 1000)
  
  print(
    plot_df %>%
      ggplot(aes(x = sample_clean_fol)) + 
      geom_bar(aes(fill = clone_subgroup_id_90_similarity_plot)) + 
      # geom_text(
      #   # data = N_clones_per_fol, 
      #   aes(x = manual_ADT_ID_plot, y = 1.02, label = Count)
      # ) +
      geom_text(
        data = N_clones_per_fol,
        aes(x = sample_clean_fol, y = 5000, label = glue("N clones total: {n}"))
      ) +
      scale_fill_manual(
        values = clone_colors, 
        labels = clone_names
      ) + 
      # scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      labs(
        x = "Site", 
        y = "Cell count", 
        # title = glue("{p}: Top 10 clones across GC B cells in {HH_fol_sample_clean} follicles"),
        # title = glue("{p}\nTop 10 clones across GC B cells from Peyer's patch follicles"),
        title = glue("{p}\nTop {n_clones} PC clones in LP sites"),
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
  
}


