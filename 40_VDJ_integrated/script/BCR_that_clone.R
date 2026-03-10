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
library(scRepertoire)

# Load data
seurat_obj <- readRDS("30_seurat_integration/out/seurat_integrated_10PCs.rds")
combined.BCR.joined <- readRDS("20_VDJ/out/combined.BCR.joined.rds") # with full length seq information 
# combined.BCR.filtered_all <- readRDS("20_VDJ/out/combined.BCR.filtered.clean_all.rds") # with full length seq information 

# nrow(combined.BCR.filtered_all)

clones <- readRDS("40_VDJ_integrated/out/top_5_clones_per_patient.rds")

# Define clone of interest 
CT <- "cluster.7833_IGLV1-40.TGCCAGTCTTATGACACCAGACTGAGGGCCACTGTGTTC"


combined.BCR.joined <- bind_rows(combined.BCR.joined)

combined.BCR.joined %>% nrow()
combined.BCR.joined %>% distinct(barcode, .keep_all = TRUE) %>% nrow()
combined.BCR.joined[duplicated(combined.BCR.joined$barcode), ]

# Adjust barcode 
combined.BCR.joined$barcode %>% tail()
colnames(seurat_obj) %>% tail()

# ------------------------------------------------------------------------------
# Fix barcodes 
# ------------------------------------------------------------------------------

# Check what number corresponds to what sample
# Create a lookup: sample name -> suffix number
sample_suffix_map <- seurat_obj[[]] %>%
  select(sample) %>%
  rownames_to_column("barcode") %>%
  mutate(suffix = str_split_i(barcode, "_", 3)) %>%
  select(sample, suffix) %>%
  distinct() %>% 
  rename(sample_high_level = sample)

sample_suffix_map

colnames(seurat_obj) %>% head()
combined.BCR.joined$barcode %>% tail()

# Then fix BCR barcodes to match integrated object
combined.BCR.joined <- combined.BCR.joined %>%
  left_join(sample_suffix_map, by = "sample_high_level") %>%
  mutate(
    barcode_integrated = paste0(
      barcode,
      "_", suffix                            # Add integration suffix
    ) %>% str_remove("_Fol-\\d+")
  ) %>% 
  rename(barcode_rm = barcode, 
         barcode = barcode_integrated)

# Check matches
table(colnames(seurat_obj) %in% combined.BCR.joined$barcode) # TRUE = 76175

combined.BCR.joined %>% nrow()
combined.BCR.joined %>% distinct(barcode, .keep_all = TRUE) %>% nrow()
combined.BCR.joined[duplicated(combined.BCR.joined$barcode), ]

seurat_obj[[]] %>% filter(!is.na(CTstrict)) %>% rownames() %>% length() # 76175
combined.BCR.joined$barcode %>% length() # 76181

# Exclude cells that are not in seurat object 
combined.BCR.joined <- combined.BCR.joined %>%
  filter(barcode %in% rownames(filter(seurat_obj[[]], !is.na(CTstrict))))

nrow(combined.BCR.joined)

# combined.BCR.joined$barcode[!(combined.BCR.joined$barcode %in% rownames(filter(seurat_obj[[]], !is.na(CTstrict))))] 

# Remove all instances of duplicated barcodes (not just the second occurrence)
combined.BCR.joined <- combined.BCR.joined %>% distinct(barcode, .keep_all = TRUE) # Remove Doublets #TODO: Look more into in 01 script 

# ------------------------------------------------------------------------------

# N cells of each clone

# Combine combined.BCR.filtered_sample and seurat_obj
seurat_obj_BCR <- combineExpression(
  combined.BCR.joined,
  seurat_obj,
  cloneCall = "strict",
  proportion = TRUE
)

for (C in unlist(clones)){
  
  n_cells <- seurat_obj_BCR[[]] %>% filter(CTstrict == C) %>% nrow()
  print(glue("{C}: {n_cells}"))
  
}

# Check 
table(seurat_obj_BCR[[]]$patient)
table(seurat_obj_BCR[[]]$patient, seurat_obj_BCR[[]]$celltype_broad) # HH119 GC B cells 20242
table(seurat_obj_BCR[[]]$CTstrict, seurat_obj_BCR[[]]$celltype_broad) %>% data.frame() %>% filter(Var1 == CT) %>% arrange(desc(Freq)) %>% head()
table(seurat_obj_BCR[[]]$CTstrict == CT)

# DimPlot of clone 
seurat_obj_BCR[[]] <- seurat_obj_BCR[[]] %>% mutate(is_clone = ifelse(CTstrict == CT, TRUE, FALSE))

DimPlot(seurat_obj_BCR, group.by = "is_clone") + labs(caption = CT)
ggsave("40_VDJ_integrated/plot/DimPlot_is_clone.png", width = 8, height = 7)

DimPlot(seurat_obj_BCR, group.by = "is_clone", split.by = "patient") + labs(caption = CT)
ggsave("40_VDJ_integrated/plot/DimPlot_is_clone_split.by_patient.png", width = 13, height = 7)

# Table of that clone across cell types
# seurat_obj_BCR_CT <- seurat_obj_BCR[[]] %>% 
#   filter(CTstrict == CT)
# 
# seurat_obj_BCR_CT

# ------------------------------------------------------------------------
# Bill plot 
# Would be nice to get a graph from each showing size of clones over the two data set with size on Y axis.
# ------------------------------------------------------------------------

# N cells patient 
n_cells_117 <- combined.BCR.joined %>% filter(patient == "HH117_Crohns") %>% nrow() # 19414
n_cells_119 <- combined.BCR.joined %>% filter(patient == "HH117_Control") %>% nrow() # 56761

# N clone types 
n_CTs_117 <- combined.BCR.joined %>% filter(patient == "HH117_Crohns") %>% select(CTstrict) %>% unique() %>% nrow() # 9263
n_CTs_119 <- combined.BCR.joined %>% filter(patient == "HH117_Control") %>% select(CTstrict) %>% unique() %>% na.omit() %>% nrow() # 28481

# N clones in clone types 
combined.BCR.joined %>% 
  filter(patient == "HH117_Crohns") %>% 
  group_by(CTstrict) %>% 
  summarize(n = n()) %>% select(n) %>% pull() %>% table() 

combined.BCR.joined %>% 
  filter(patient == "HH117_Control") %>% 
  group_by(CTstrict) %>% 
  summarize(n = n()) %>% select(n) %>% pull() %>% table() 

# plot
max_117 <- 20
combined.BCR.joined %>% 
  filter(patient == "HH117_Crohns") %>% 
  group_by(CTstrict) %>% 
  summarize(n = n()) %>% 
  filter(n >= max_117) %>% 
  ggplot(aes(y = reorder(CTstrict, n), x = n)) + 
  geom_col() +
  geom_text(aes(label = n), hjust = -0.1, size = 3) + 
  labs(
    subtitle = "HH117_Crohns",
    caption = glue("Only showing clonal types >= {max_117} clones\nN cells: {n_cells_117}\nN CT: {n_CTs_117}"),
    y = "CTstrict"
  ) + 
  theme_bw()

ggsave("20_VDJ/plot/BCR_clonal_abundance/N_CTstrict_HH117.png", height = 10, width = 10)

max_119 <- 40
combined.BCR.joined %>% 
  filter(patient == "HH117_Control") %>% 
  group_by(CTstrict) %>% 
  summarize(n = n()) %>% 
  filter(n >= max_119) %>% 
  ggplot(aes(y = reorder(CTstrict, n), x = n)) + 
  geom_col() +
  geom_text(aes(label = n), hjust = -0.1, size = 3) + 
  labs(
    subtitle = "HH119_Control",
    caption = glue("Only showing clonal types >= {max_119} clones\nN cells: {n_cells_119}\nN CT: {n_CTs_119}"),
    y = "CTstrict"
  ) + 
  theme_bw()
  
ggsave("20_VDJ/plot/BCR_clonal_abundance/N_CTstrict_HH119.png", height = 10, width = 12)






# ------------------------------------------------------------------------
# Define clones in clonetype 
# ------------------------------------------------------------------------

clones
# CT <- "IGH:Cluster.4529.IGHV4-34_IGLC:Cluster.172.IGLV1-40"
CT <- "IGH:Cluster.370.IGHV1-8_IGLC:Cluster.1.IGKV1-5"


df_clone_nt <- combined.BCR.joined %>% 
  filter(CTstrict == CT) %>% 
  mutate( 
    IGH_Vgene = IGH %>% str_split_i("\\.", 1),
    IGH_Jgene = IGH %>% str_split_i("\\.", 3),
    IGH_CDR3 = cdr3_nt1, # Heavy
    # IGLC_CDR3_nt = cdr3_nt2, # Light 
    Ig_class = CTgene %>% str_split_i("_", 1) %>% str_split_i("\\.", 4)
  )

# Heavy chain: CTscrict, full seqeunce  abundance
clone_df <- df_clone_nt %>% 
  select(-barcode) %>% 
  summarise(abundance = n(), .by = c(CTstrict, IGH_Vgene, IGH_Jgene, IGH_full_sequence, IGH_CDR3)) %>% 
  arrange(desc(abundance))

clone_df

# unique clones 
nrow(clone_df)

sum(clone_df$abundance)

####
#### PLOT
####

clone_df %>% 
  ggplot(aes(y = IGH_CDR3, x = abundance)) + 
  geom_col()

####
# ISOTYPES 
####

table(df_clone_nt$sample) %>% data.frame() %>% arrange(desc(Freq))

table(df_clone_nt$Ig_class) %>% data.frame() %>% arrange(desc(Freq))

table(df_clone_nt$Ig_class, df_clone_nt$sample)

table(df_clone_nt$Ig_class, df_clone_nt$sample) %>%
  as.data.frame() %>%
  mutate(Var2 = Var2 %>% str_remove_all("-HLADR-AND-CD19-AND-GC-AND-TFH|-CD19-AND-GC-AND-PB-AND-TFH|-HLADR-AND-CD19|-PC")) %>% 
  ggplot(aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white", size = 3) +
  scale_fill_gradient(low = "#a8c5da", high = "#08306b") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  labs(x = "", y = "Ig Class", fill = "Count")

ggsave("20_VDJ/plot/BCR_thatclone/Ig_class_VS_sample.png", width = 13, height = 6)

# Look into 
# combined.BCR.filtered.clean_all <- combined.BCR.filtered.clean_all %>% 
#   mutate(Ig_class = CTgene %>% str_split_i("_", 1) %>% 
#            str_split_i("\\.", 4))

library(Biostrings) # writing fasta files 

seqs <- DNAStringSet(clone_df$IGH_full_sequence)
names(seqs) <- glue("Clone:{clone_df$CTstrict}|CDR3nt:{clone_df$IGH_CDR3}|Abundance:{clone_df$abundance}")

sprintf(">Clone:%s|CDR3nt:%s|Abundance:%s", 
        clone_df$CTstrict,
        clone_df$IGH_CDR3,
        clone_df$abundance)    

writeXStringSet(seqs, filepath = glue("20_VDJ/fasta/clone_{CT}.fasta"))

