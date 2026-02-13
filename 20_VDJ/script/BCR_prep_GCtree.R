getwd()

library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(readxl)

combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.clean.rds")

names(combined.BCR.filtered)

# ------------------------------------------------------------------------------
# Clonal abundance per sample
# ------------------------------------------------------------------------------

# Count of clone sizes
clone_abundance <- lapply(combined.BCR.filtered, function(x) x %>% 
  summarise(abundance = n(), .by = c(CTstrict)) %>% arrange(desc(abundance)))

# Check numbers
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow() # N cells
combined.BCR.filtered$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1`$CTstrict %>% unique() %>% length() # N unique clones
clone_abundance$`HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH_Fol-1` %>% nrow() # N unique clones

# Plot
for (sample_name in names(clone_abundance)){

  min_abundance <- ifelse(str_detect(sample_name, "Fol"), 1, 7)
  
  clone_abundance[[sample_name]] %>% 
    filter(abundance > min_abundance) %>% 
    ggplot(aes(x = abundance, y = reorder(CTstrict, abundance))) + 
    geom_col() + 
    geom_text(aes(label = abundance), hjust = -0.1, size = 3) + 
    labs(
      title = "Abundance of clones (CTstrict)",
      subtitle = sample_name, 
      caption = "Only abundance > 1 is included."
    ) + 
    theme_bw() 
  
  ggsave(glue("20_VDJ/plot/BCR_clonal_abundance/{sample_name}.png"), width = 20, height = 10)
  
}

# Export as excel 
out_file <- "20_VDJ/table/BCR_clonal_abundance.xlsx"

# Prep names
names(clone_abundance) <- names(clone_abundance) %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")

# Use openxlsx::write.xlsx, which takes the named list and writes
# each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
openxlsx::write.xlsx(
  x = clone_abundance,
  file = out_file, 
  overwrite = TRUE # Overwrite the file if it already exists
)
