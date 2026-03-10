# getwd()

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

# celltype_colors
source("10_broad_annotation/script/color_palette.R")

combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")

outdir <- "20_VDJ/plot/BCR_06_GermlineSequence"
# dir.create(outdir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Get clones within GC B cells 
# ------------------------------------------------------------------------------

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

CT <- top_GC_clones$HH119[[1]]
combined.BCR.filtered$HH119 %>% 
  filter(manual_cluster == CT) %>% 
  nrow()

combined.BCR.filtered$HH119 %>% 
  filter(manual_cluster == CT) %>% 
  select(IGH_V_gene, IGH_J_gene, IGLC_V_gene, IGLC_J_gene) %>% distinct()

# ------------------------------------------------------------------------------
# Look at top clones across samples and cell types. 
# ------------------------------------------------------------------------------

for (patient in patients){
  
  # patient <- "HH119"
  
  for (clone_nr in 1:length(top_GC_clones[[patient]])){
    
    # clone_nr <- 1
    
    # CT <- top_GC_clones$HH119_Control[1]
    CT <- top_GC_clones[[patient]][clone_nr]
    
    # Barplot with CT across follicles  
    df_plot <- combined.BCR.filtered[[patient]] %>% 
      filter(manual_cluster == CT) %>%
      # mutate(sample = sample %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")) %>% 
      summarize(n = n(), .by = c(sample, celltype_broad))
    
    n_cells_total <- df_plot$n %>% sum()
    
    # Access genes of CT
    CT_genes <- combined.BCR.filtered[[patient]] %>% 
      filter(manual_cluster == CT) %>% 
      mutate(
        IGH_V_gene = IGH %>% str_split_i("\\.", 1),
        IGH_J_gene = IGH %>% str_split_i("\\.", 3),
        IGLC_V_gene = IGLC %>% str_split_i("\\.", 1),
        IGLC_J_gene = IGLC %>% str_split_i("\\.", 2),
      ) %>% 
      select(IGH_V_gene, IGH_J_gene, IGLC_V_gene, IGLC_J_gene) %>% 
      distinct() %>% 
      paste(collapse = "_")
    
    df_plot %>% 
      ggplot(aes(x = sample, y = n, fill = celltype_broad)) + 
      geom_col() + 
      scale_fill_manual(values = celltype_colors) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = "", 
        y = "N cells",
        title = glue("{patient}: {CT}"),
        subtitle = CT_genes, 
        caption = glue("Clone nr {clone_nr}\nN cells total: {n_cells_total}")
      )
    
    
    ggsave(glue("{outdir}/{patient}_clone{clone_nr}_barplot.png"), width = 14, height = 8)
    
  }
  
}

# ------------------------------------------------------------------------------
# Access abundance of one CT
# ------------------------------------------------------------------------------

patient <- "HH119"
CT <- top_GC_clones[[patient]][[1]]

combined.BCR.filtered$HH119 %>% 
  filter(manual_cluster == CT) %>% 
  summarize(n = n(), .by = c(sample, celltype_broad)) %>% 
  arrange(desc(n)) %>% 
  head(10)

# ------------------------------------------------------------------------------
# Access genes of one CT
# ------------------------------------------------------------------------------

# Get the V and J genes from this CT
CT_genes <- combined.BCR.filtered[[patient]] %>% 
  filter(manual_cluster == CT) %>% 
  mutate(
    IGH_V_gene = IGH %>% str_split_i("\\.", 1),
    IGH_J_gene = IGH %>% str_split_i("\\.", 3),
    IGLC_V_gene = IGLC %>% str_split_i("\\.", 1),
    IGLC_J_gene = IGLC %>% str_split_i("\\.", 2),
  ) %>% 
  select(IGH_V_gene, IGH_J_gene, IGLC_V_gene, IGLC_J_gene) %>% 
  distinct()

CT_genes

IGH_v_gene <- CT_genes$IGH_V_gene
IGH_j_gene <- CT_genes$IGH_J_gene

IGL_v_gene <- CT_genes$IGLC_V_gene
IGL_j_gene <- CT_genes$IGLC_J_gene


# ------------------------------------------------------------------------------
# Find germline sequence in Naive 
# ------------------------------------------------------------------------------

# Get the V and J genes from this CT
CT_genes <- combined.BCR.filtered[[patient]] %>% 
  filter(manual_cluster == CT) %>% 
  mutate(
    IGH_V_gene = IGH %>% str_split_i("\\.", 1),
    IGH_J_gene = IGH %>% str_split_i("\\.", 3),
    IGLC_V_gene = IGLC %>% str_split_i("\\.", 1),
    IGLC_J_gene = IGLC %>% str_split_i("\\.", 2),
  ) %>% 
  select(IGH_V_gene, IGH_J_gene, IGLC_V_gene, IGLC_J_gene) %>% 
  distinct()
  
CT_genes

IGH_v_gene <- CT_genes$IGH_V_gene
IGH_j_gene <- CT_genes$IGH_J_gene

IGL_v_gene <- CT_genes$IGLC_V_gene
IGL_j_gene <- CT_genes$IGLC_J_gene


# In 01
# V gene sequence: cdr1+fwr2+cdr2+fwr3
# J gene sequence: fwr4
# What about light chain? Do we want to define clones from this?

combined.BCR.filtered[[patient]] %>% 
  filter(celltype_broad == "Naïve_memory_B_cells") %>% 
  mutate(
    IGH_V_gene = IGH %>% str_split_i("\\.", 1),
    IGH_J_gene = IGH %>% str_split_i("\\.", 3),
    IGLC_V_gene = IGLC %>% str_split_i("\\.", 1),
    IGLC_J_gene = IGLC %>% str_split_i("\\.", 2),
  ) %>% 
  filter(IGH_V_gene == IGH_v_gene) %>% 
  select(IGH_full_sequence) %>% 
  distinct()




