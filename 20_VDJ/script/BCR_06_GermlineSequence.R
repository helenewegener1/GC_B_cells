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

outdir <- "20_VDJ/plot/BCR_06_GermlineSequence"
dir.create(outdir)

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
# Look at top clones across samples
# ------------------------------------------------------------------------------

for (patient in patients){
  
  # patient <- "HH117_Crohns"
  # patient <- "HH119_Control"
  
  for (clone_nr in 1:length(top_GC_clones[[patient]])){
    
    # clone_nr <- 1
    
    # CT <- top_GC_clones$HH119_Control[1]
    CT <- top_GC_clones[[patient]][clone_nr]
    
    # Barplot with CT across follicles  
    combined.BCR.joined_all %>% 
      # filter(CTstrict == CT & str_detect(sample, "Fol")) %>%
      filter(CTstrict == CT) %>%
      mutate(sample = sample %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")) %>% 
      summarize(n = n(), .by = c(sample, celltype_broad)) %>% 
      ggplot(aes(x = sample, y = n, fill = celltype_broad)) + 
      geom_col() + 
      scale_fill_manual(values = celltype_colors) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = "", 
        y = "N cells",
        title = patient,
        subtitle = CT
      )
    
    ggsave(glue("{outdir}/{patient}_{CT}_barplot.png"), width = 14, height = 7)
    
  }
  
}











