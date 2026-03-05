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

source("10_broad_annotation/script/color_palette.R")
# celltype_colors

combined.BCR.joined <- readRDS("20_VDJ/out/combined.BCR.joined.rds")
combined.BCR.joined_all <- combined.BCR.joined %>% bind_rows()

# outdir <- "20_VDJ/plot/BCR_06_GermlineSequence"
# dir.create(outdir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Get clones within GC B cells 
# ------------------------------------------------------------------------------

patients <- combined.BCR.joined_all$patient %>% unique()

top_GC_clones <- lapply(patients, function(x){
  
  combined.BCR.joined_all %>% 
    filter(patient == x & celltype_broad == "GC_B_cells" & !is.na(CTstrict)) %>% 
    summarise(n = n(), .by = "CTstrict") %>% 
    arrange(desc(n)) %>% 
    head(n = 10) %>% 
    pull(CTstrict)
  
}) %>% 
  setNames(patients)

# ------------------------------------------------------------------------------
# Look at top clones across samples and cell types. 
# ------------------------------------------------------------------------------

for (patient in patients){
  
  # patient <- "HH117_Crohns"
  # patient <- "HH119_Control"
  
  for (clone_nr in 1:length(top_GC_clones[[patient]])){
    
    # clone_nr <- 1
    
    # CT <- top_GC_clones$HH119_Control[1]
    CT <- top_GC_clones[[patient]][clone_nr]
    
    # Barplot with CT across follicles  
    df_plot <- combined.BCR.joined_all %>% 
      # filter(CTstrict == CT & str_detect(sample, "Fol")) %>%
      filter(CTstrict == CT) %>%
      mutate(sample = sample %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")) %>% 
      summarize(n = n(), .by = c(sample, celltype_broad))
    
    n_cells_total <- df_plot$n %>% sum()
      
    df_plot %>% 
      ggplot(aes(x = sample, y = n, fill = celltype_broad)) + 
      geom_col() + 
      scale_fill_manual(values = celltype_colors) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = "", 
        y = "N cells",
        title = patient,
        subtitle = CT, 
        caption = glue("Clone nr {clone_nr}\nN cells total: {n_cells_total}")
      )
    
    
    ggsave(glue("{outdir}/{patient}_clone{clone_nr}_barplot.png"), width = 14, height = 8)
    
  }
  
}

# ------------------------------------------------------------------------------
# Find germline sequence in Naive 
# ------------------------------------------------------------------------------

CT <- "cluster.7833_IGLV1-40.TGCCAGTCTTATGACACCAGACTGAGGGCCACTGTGTTC"
# CT <- "cluster.1"
patient <- "HH119_Control"

# Get the V and J genes from this CT
CT_genes <- combined.BCR.joined_all %>% 
  filter(CTstrict == CT, patient == patient) %>% 
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

combined.BCR.joined_all %>% 
  filter(patient == patient, celltype_broad == "Naïve_memory_B_cells") %>% 
  mutate(
    IGH_V_gene = IGH %>% str_split_i("\\.", 1),
    IGH_J_gene = IGH %>% str_split_i("\\.", 3),
    IGLC_V_gene = IGLC %>% str_split_i("\\.", 1),
    IGLC_J_gene = IGLC %>% str_split_i("\\.", 2),
  ) %>% 
  filter(IGH_V_gene == IGH_v_gene) %>% 
  select(IGH_full_sequence) %>% 
  distinct()




