setwd("~/Documents/projects/project_Bcells/Gina_project/")

library(tidyverse)
library(xlsx)

#################################### HH117 ##################################### 

# Load data
HH117_meta <- xlsx::read.xlsx("00_data/HH117 Calculation Hashtag Information_latestVer.xlsx", 
                              sheetName = "Hashtag Information HH117", 
                              startRow = 3)

# Clean up 
HH117_meta_clean <- HH117_meta %>% 
  select(Follicle, Cite.Seq..Hashing.AB, starts_with("Sequence")) %>% # Select relevant columns 
  filter(str_starts(Cite.Seq..Hashing.AB, "Hashtag")) %>% # Select relevant rows 
  mutate(Follicle = str_trim(Follicle), # Remove white space
         Sequence = coalesce(!!!select(., starts_with("Sequence")))) %>% # collapse hashtag sequence columns into one
  select(Follicle, Cite.Seq..Hashing.AB, Sequence) %>% 
  rename("id" = Cite.Seq..Hashing.AB, 
         "name" = Follicle, 
         "sequence" = Sequence)

# Rename id (hashtag IDs)
id_tmp <- str_replace(HH117_meta_clean$id , " ", "_")
id_new <- paste0(id_tmp, "_TotalSeqC")
HH117_meta_clean$id <- id_new

# Rename name (demultiplexed sub sample name)
name_new <- HH117_meta_clean$name %>% 
  str_replace("Gina: ", "") %>% 
  str_replace_all(" ", "") %>% 
  str_replace_all("\\.", "_")

HH117_meta_clean$name <- name_new

# Add new columns 
HH117_meta_clean$read <- "R2"
HH117_meta_clean$pattern <- "5PNNNNNNNNNN(BC)"
HH117_meta_clean$feature_type <- "Antibody Capture"

# Order columns 
# output header: id,name,read,pattern,sequence,feature_type
HH117_meta_clean <- HH117_meta_clean %>% select(id, name, read, pattern, sequence, feature_type)

# Export as csv
# output: feature_ref.csv
write.csv(HH117_meta_clean, "01_prep_feature_ref/out/feature_ref_HH117.csv", quote = FALSE, row.names = FALSE)

rm(HH117_meta, HH117_meta_clean)

################################ HH119 - Pool 1 ################################ 

# Load data
HH119_meta <- xlsx::read.xlsx("00_data/HH119Calculation Hashtag Information_latestVersion.xlsx", 
                              sheetName = "Hashtag Information HH119", 
                              startRow = 3)

# Row index (start) of Pool 2
pool_2_start_index <- min(which(HH119_meta$Follicle == "Gina: Fol.pool 2"))

# Clean up 
HH119_pool_1_meta_clean <- HH119_meta %>%
  select(Follicle, Cite.Seq..Hashing.AB, starts_with("Sequence")) %>% # Select relevant columns 
  slice_head(n = pool_2_start_index - 1) %>% # Only include Pool 1 rows 
  filter(str_starts(Cite.Seq..Hashing.AB, "Hashtag")) %>% # Select relevant rows 
  mutate(Follicle = str_trim(Follicle), # Remove white space
         Sequence = coalesce(!!!select(., starts_with("Sequence")))) %>% # collapse hashtag sequence columns into one
  select(Follicle, Cite.Seq..Hashing.AB, Sequence) %>% 
  rename("id" = Cite.Seq..Hashing.AB, 
         "name" = Follicle, 
         "sequence" = Sequence)

# Rename id (hashtag IDs)
id_tmp <- str_replace(HH119_pool_1_meta_clean$id , " ", "_")
id_new <- paste0(id_tmp, "_TotalSeqC")
HH119_pool_1_meta_clean$id <- id_new

# Rename name (demultiplexed sub sample name)
name_new <- HH119_pool_1_meta_clean$name %>% 
  str_replace_all(" ", "") %>% 
  str_replace_all("\\.", "_")

HH119_pool_1_meta_clean$name <- name_new

# Add new columns 
HH119_pool_1_meta_clean$read <- "R2"
HH119_pool_1_meta_clean$pattern <- "5PNNNNNNNNNN(BC)"
HH119_pool_1_meta_clean$feature_type <- "Antibody Capture"

# Order columns 
# output header: id,name,read,pattern,sequence,feature_type
HH119_pool_1_meta_clean <- HH119_pool_1_meta_clean %>% select(id, name, read, pattern, sequence, feature_type)

# Export as csv
# output: feature_ref.csv
write.csv(HH119_pool_1_meta_clean, "01_prep_feature_ref/out/feature_ref_HH119_pool_1.csv", quote = FALSE, row.names = FALSE)

rm(HH119_pool_1_meta_clean)


################################ HH119 - Pool 2 ################################ 

# Clean up 
HH119_pool_2_meta_clean <- HH119_meta %>%
  select(Follicle, Cite.Seq..Hashing.AB, starts_with("Sequence")) %>% # Select relevant columns 
  slice(pool_2_start_index:n()) %>% # Only include Pool 2 rows 
  filter(str_starts(Cite.Seq..Hashing.AB, "Hashtag")) %>% # Select relevant rows 
  mutate(Follicle = str_trim(Follicle), # Remove white space
         Sequence = coalesce(!!!select(., starts_with("Sequence")))) %>% # collapse hashtag sequence columns into one
  select(Follicle, Cite.Seq..Hashing.AB, Sequence) %>% 
  rename("id" = Cite.Seq..Hashing.AB, 
         "name" = Follicle, 
         "sequence" = Sequence)

# Rename id (hashtag IDs)
id_tmp <- str_replace(HH119_pool_2_meta_clean$id , " ", "_")
id_new <- paste0(id_tmp, "_TotalSeqC")
HH119_pool_2_meta_clean$id <- id_new

# Rename name (demultiplexed sub sample name)
name_new <- HH119_pool_2_meta_clean$name %>% 
  str_replace_all(" ", "") %>% 
  str_replace_all("\\.", "_")

HH119_pool_2_meta_clean$name <- name_new

# Add new columns 
HH119_pool_2_meta_clean$read <- "R2"
HH119_pool_2_meta_clean$pattern <- "5PNNNNNNNNNN(BC)"
HH119_pool_2_meta_clean$feature_type <- "Antibody Capture"

# Order columns 
# output header: id,name,read,pattern,sequence,feature_type
HH119_pool_2_meta_clean <- HH119_pool_2_meta_clean %>% select(id, name, read, pattern, sequence, feature_type)

# Export as csv
# output: feature_ref.csv
write.csv(HH119_pool_2_meta_clean, "01_prep_feature_ref/out/feature_ref_HH119_pool_2.csv", quote = FALSE, row.names = FALSE)

rm(HH119_pool_2_meta_clean)
