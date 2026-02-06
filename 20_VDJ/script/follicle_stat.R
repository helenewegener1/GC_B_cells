library(tidyverse)
library(readxl)

# Load data 
sheet_names <- readxl::excel_sheets("10_ADT_demultiplex/table/fol_freq.xlsx")

initial <- list()
BCR_filt <- list()
# TCR_filt <- list()

for (sheet in sheet_names){
  
  initial[[sheet]] <- read_xlsx("10_ADT_demultiplex/table/fol_freq.xlsx", sheet = sheet)
  BCR_filt[[sheet]] <- read_xlsx("20_VDJ/table/BCR_filtered_fol_freq.xlsx", sheet = sheet)
  # TCR_filt[[sheet]] <- read_xlsx("20_VDJ/table/TCR_filtered_fol_freq.xlsx", sheet = sheet)
  
}

# Update
# Add identifiers if needed
initial_clean <- lapply(names(initial), function(x) initial[[x]] %>% mutate(version = "initial", sample = x))
BCR_filt_clean <- lapply(names(BCR_filt), 
                         function(x) BCR_filt[[x]] %>% mutate(version = "BCR_filt", sample = x, Count = as.integer(Count)))

# Combine both lists into one data frame
combined_df <- bind_rows(c(initial_clean, BCR_filt_clean))

## Test
combined_df %>% filter(Fol == "Fol-1", sample == "HH117")

combined_df %>% filter(sample == "HH119-CD19-Pool1") %>% view()

# Plot
for (sample_name in sheet_names){
  
  # sample <- "HH117"
  # sample_name <- "HH119-CD19-Pool2"
  
  combined_df %>% 
    filter(sample == sample_name) %>% 
    ggplot(aes(x = Fol, 
               y = Count,
               fill = version)) +
    geom_col(position = "dodge") + 
    scale_fill_manual(values = c("#798E87FF", "#C27D38FF")) + 
    labs(
      title = sample_name
    ) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(glue("20_VDJ/plot/fol_stat/{sample_name}.png"), width = 12, height = 7)
  
}

