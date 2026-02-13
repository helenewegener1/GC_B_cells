library(tidyverse)
library(readxl)
library(wesanderson)

 # ------------------------------------------------------------------------------
# Cell frequency with BCR after filtering (No fol information)
# ------------------------------------------------------------------------------
# combined.BCR.filtered <- readRDS("20_VDJ/out/combined.BCR.filtered.rds")
# source("20_VDJ/script/functions.R")
# 
# # Apply the function
# combined.BCR.combined_follicles <- combine_follicles(combined.BCR.filtered)
# 
# cell_counts_combined <- data.frame(
#   sample = names(combined.BCR.combined_follicles),
#   n_cells = sapply(combined.BCR.combined_follicles, nrow),
#   type = "Combined"
# )
# 
# # Plot comparison
# ggplot(cell_counts_combined, aes(x = sample, y = n_cells)) +
#   geom_col(fill = "steelblue") +
#   geom_text(aes(label = n_cells), vjust = -0.5, size = 3) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 8)) +
#   labs(
#     title = "Cell Counts After combineBCR",
#     x = "",
#     y = "Number of cells"
#   )
# 
# ggsave("20_VDJ/plot/N_cell_stat/N_cells_combined_BCR_filtered.png", width = 12, height = 7)

# ------------------------------------------------------------------------------
# Cell frequency after filtering - per follicle
# ------------------------------------------------------------------------------

# Load data 
BCR_sheet_names <- readxl::excel_sheets("10_ADT_demultiplex/table/fol_freq.xlsx")
TCR_sheet_names <- readxl::excel_sheets("20_VDJ/table/TCR_filtered_fol_freq.xlsx")

initial <- list()
BCR_filt <- list()
TCR_filt <- list()

for (sheet in BCR_sheet_names){
  
  initial[[sheet]] <- read_xlsx("10_ADT_demultiplex/table/fol_freq.xlsx", sheet = sheet)
  BCR_filt[[sheet]] <- read_xlsx("20_VDJ/table/BCR_filtered_fol_freq.xlsx", sheet = sheet)
  
  if (sheet %in% TCR_sheet_names){
    TCR_filt[[sheet]] <- read_xlsx("20_VDJ/table/TCR_filtered_fol_freq.xlsx", sheet = sheet)
  }
  
}

# Update
# Add identifiers if needed
initial_clean <- lapply(names(initial), function(x) initial[[x]] %>% mutate(version = "initial", sample = x))
BCR_filt_clean <- lapply(names(BCR_filt), 
                         function(x) BCR_filt[[x]] %>% mutate(version = "BCR_filt", sample = x, Count = as.integer(Count)))
TCR_filt_clean <- lapply(names(TCR_filt), 
                         function(x) TCR_filt[[x]] %>% mutate(version = "TCR_filt", sample = x, Count = as.integer(Count)))

# Combine both lists into one data frame
combined_df <- bind_rows(c(initial_clean, BCR_filt_clean, TCR_filt_clean))
combined_df <- combined_df %>% mutate(version = factor(version, levels = c("initial", "BCR_filt", "TCR_filt")))

## Test
combined_df %>% filter(Fol == "Fol-1", sample == "HH117")

combined_df %>% filter(sample == "HH119-CD19-Pool1") %>% view()

# Plot
for (sample_name in BCR_sheet_names){
  
  # sample <- "HH117"
  # sample_name <- "HH119-CD19-Pool2"
  
  combined_df %>% 
    filter(sample == sample_name) %>% 
    ggplot(aes(x = Fol, 
               y = Count,
               fill = version)) +
    geom_col(position = "dodge") + 
    scale_fill_manual(values = c("#CCC591FF", "#798E87FF", "#C27D38FF")) + 
    # scale_fill_manual(values = c("initial" = "#CCC591FF", "BCR_filt" = "#798E87FF", "TCR_filt" = "#C27D38FF")) + 
    labs(
      title = sample_name,
      y = "N cells"
    ) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(glue("20_VDJ/plot/N_cell_stat/{sample_name}.png"), width = 12, height = 7)
  
}

