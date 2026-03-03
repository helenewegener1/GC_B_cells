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
library(readr)


combined.BCR.filtered.clean_all <- readRDS("20_VDJ/out/combined.BCR.filtered.clean_all.rds")

combined.BCR.filtered.clean_all <- combined.BCR.filtered.clean_all %>% 
  mutate(Ig_class = CTgene %>% str_split_i("_", 1) %>% 
           str_split_i("\\.", 4))


# combined.BCR.filtered.clean_all %>% filter(Ig_class == "NA") %>% select(CTgene)
combined.BCR.filtered.clean_all$Ig_class %>% table(useNA = "always")

# # ------------------------------------------------------------------------------
# # Classes per sample
# # ------------------------------------------------------------------------------
# 
# sample_names <- unique(combined.BCR.filtered.clean_all$sample_clean)
# 
# for (sample_name in sample_names){
#   
#   # sample_name <- "HH117-SILP-INF"
#   
#   df <- combined.BCR.filtered.clean_all %>% filter(sample_clean == sample_name)
#   
#   n_cells <- nrow(df)
#   
#   df %>% 
#     ggplot(aes(x = Ig_class)) + 
#     geom_bar() + 
#     theme_bw() + 
#     labs(
#       title = "Abundance of Ig classes",
#       subtitle = sample_name,
#       caption = glue("N cells: {n_cells}")
#     )
#   
#   ggsave(glue("20_VDJ/plot/BCR_IgClassesAbundance/count/BCR_IgClassesAbundance_{sample_name}_count.png"), width = 7.5, height = 5.5)
#   
#   df %>% 
#     ggplot(aes(x = Ig_class, y = after_stat(count/sum(count)*100))) + 
#     geom_bar() + 
#     theme_bw() + 
#     labs(
#       title = "Abundance of Ig classes",
#       subtitle = sample_name,
#       y = "Percentage (%)",
#       caption = glue("N cells: {n_cells}")
#     )
#   
#   ggsave(glue("20_VDJ/plot/BCR_IgClassesAbundance/percentage/BCR_IgClassesAbundance_{sample_name}_percentage.png"), width = 7.5, height = 5.5)
#   
# }

# ------------------------------------------------------------------------------
# Colors 
# ------------------------------------------------------------------------------

combined.BCR.filtered.clean_all$Ig_class %>% unique()
ig_class_colors <- c("IGHA1" = "#2166AC", "IGHA2" = "#6BAED6", 
                     "IGHG1" = "#1A7A3C", "IGHG2" = "#4DAF4A", "IGHG3" = "#A1D76A", "IGHG4" = "#D9F0A3", 
                     "IGHM" = "#C1392B", "IGHD" = "#E8735A", "IGHE" = "#9E3A8C", "NA" = "#AAAAAA")


# ------------------------------------------------------------------------------
# Classes per sample - combined plot
# ------------------------------------------------------------------------------

# Count plot
combined.BCR.filtered.clean_all %>%
  ggplot(aes(x = sample_clean, fill = Ig_class)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  # scale_fill_manual(values = ig_class_colors, na.value = "#AAAAAA") +
  scale_fill_manual(values = ig_class_colors) +
  labs(
    title = "Abundance of Ig classes",
    x = "",
    y = "Count",
    fill = "Ig Class"
  )

ggsave("20_VDJ/plot/BCR_IgClassesAbundance/BCR_IgClassesAbundance_combined_count.png", width = 10, height = 5.5)

# Percentage plot
combined.BCR.filtered.clean_all %>%
  group_by(sample_clean, Ig_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample_clean) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = sample_clean, y = pct, fill = Ig_class)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_fill_manual(values = ig_class_colors, na.value = "#AAAAAA") +
  scale_fill_manual(values = ig_class_colors) +
  labs(
    title = "Abundance of Ig classes",
    x = "",
    y = "Percentage (%)",
    fill = "Ig Class"
  )

ggsave("20_VDJ/plot/BCR_IgClassesAbundance/BCR_IgClassesAbundance_combined_percentage.png", width = 10, height = 5.5)

# ------------------------------------------------------------------------------
# Classes per fol for PP samples
# ------------------------------------------------------------------------------

sample_names <- combined.BCR.filtered.clean_all$sample_clean %>% unique()
fol_samples <- grep("PP", sample_names, value = TRUE)

lapply(fol_samples, function(x) {
  
  # Count plot
  combined.BCR.filtered.clean_all %>%
    filter(sample_clean == x) %>% 
    mutate(fol = sample %>% str_split_i("_", 2)) %>% 
    ggplot(aes(x = fol, fill = Ig_class)) +
    geom_bar() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
    # scale_fill_manual(values = ig_class_colors, na.value = "#AAAAAA") +
    scale_fill_manual(values = ig_class_colors) +
    labs(
      title = "Abundance of Ig classes",
      subtitle = x,
      x = "",
      y = "Count",
      fill = "Ig Class"
    )
  
  ggsave(glue("20_VDJ/plot/BCR_IgClassesAbundance/BCR_IgClassesAbundance_{x}_count.png"), width = 10, height = 5.5)
  
  # Percentage plot
  combined.BCR.filtered.clean_all %>%
    filter(sample_clean == x) %>% 
    mutate(fol = sample %>% str_split_i("_", 2)) %>% 
    group_by(fol, Ig_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(fol) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ggplot(aes(x = fol, y = pct, fill = Ig_class)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # scale_fill_manual(values = ig_class_colors, na.value = "#AAAAAA") +
    scale_fill_manual(values = ig_class_colors) +
    labs(
      title = "Abundance of Ig classes",
      subtitle = x,
      x = "",
      y = "Percentage (%)",
      fill = "Ig Class"
    )
  
  ggsave(glue("20_VDJ/plot/BCR_IgClassesAbundance/BCR_IgClassesAbundance_{x}_percentage.png"), width = 10, height = 5.5)
  
  
  
})


