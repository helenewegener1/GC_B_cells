library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

LC_list <- readRDS("45_immcantation/out/rds/resolve_LC_list.rds")

HH <- "HH119"
align_free <- read.delim(glue("47_alignment_free_clone_identification/out/{HH}_clone_output.tsv"))

LC <- LC_list[[HH]]

LC %>% filter(locus == "IGH") %>% nrow()
align_free %>% nrow()

# ------------------------------------------------------------------------------
# Inspect clones
# ------------------------------------------------------------------------------

LC %>% filter(locus == "IGH") %>% count(clone_id, sort = TRUE) %>% head()
align_free %>% count(CLONE, sort = TRUE) %>% head()

# Take SCOPer's largest clone and check if Lindenbaum agrees
# LC_clone_cell_id <- LC %>% filter(clone_id == "578") %>% pull(cell_id)
LC_clone_cell_id <- LC %>% filter(clone_id == "4516") %>% pull(cell_id)

align_free %>% filter(ID %in% LC_clone_cell_id) %>% 
  count(CLONE, sort = TRUE) %>% head(n = 10)

# ------------------------------------------------------------------------------
# Visualize 
# ------------------------------------------------------------------------------

outdir <- glue("47_alignment_free_clone_identification/plot/LC_vs_align_free")
dir.create(outdir, showWarnings = FALSE)

# Clean dfs
LC_clean <- LC %>% filter(locus == "IGH") %>% select(cell_id, clone_id)
align_free_clean <- align_free %>%  select(ID, CLONE)

# Define top clones
top_clones_names <- list(
  "LC" = LC_clean %>% count(clone_id, sort = TRUE) %>% head(10) %>% pull(clone_id),
  "align_free" = align_free_clean %>% count(CLONE, sort = TRUE) %>% head(10) %>% pull(CLONE) %>% as.character()
)

top_clones_count <- list(
  "LC" = LC_clean %>% count(clone_id, sort = TRUE) %>% head(10) %>% pull(n),
  "align_free" = align_free_clean %>% count(CLONE, sort = TRUE) %>% head(10) %>% pull(n) 
)

names(top_clones_names$LC) <- paste("Clone", 1:10, "\nN cells: ", top_clones_count$LC)
names(top_clones_names$align_free) <- paste("Clone", 1:10, "\nN cells: ", top_clones_count$align_free)


# Make table
compare <- table(LC_clean$clone_id, align_free_clean$CLONE) %>% 
  as.data.frame() %>%
  setNames(c("LC_clones", "align_free_clones", "Count")) 
 
# LC clones in align_free 
min_count <- 4
compare %>% 
  filter(LC_clones %in% top_clones_names$LC & Count > min_count) %>% 
  ggplot(aes(x = LC_clones, y = align_free_clones, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Count), size = 3, color = "black") +
  scale_fill_gradient(low = "#EEF3FB", high = "#185FA5") +
  scale_x_discrete(
    limits = top_clones_names$LC,
    labels = setNames(names(top_clones_names$LC), top_clones_names$LC)
  ) + 
  labs(
    x     = "LC clones",
    y     = "Alignment free clones",
    title = glue("{HH}: Top 10 LC clusters VS Alignment free clones"),
    caption = glue("Minimum count: {min_count}")
  ) +
  theme_minimal(base_size = 10) + 
  theme(
    # axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "none"
  )

ggsave(glue("{outdir}/{HH}_LC_VS_align_free_clones.png"), width = 10, height = 9)

# align_free clones in LC
min_count <- 4
compare %>% 
  filter(align_free_clones %in% top_clones_names$align_free & Count > min_count) %>% 
  ggplot(aes(x = align_free_clones, y = LC_clones, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Count), size = 3, color = "black") +
  scale_fill_gradient(low = "#EEF3FB", high = "#185FA5") +
  scale_x_discrete(
    limits = top_clones_names$align_free,
    labels = setNames(names(top_clones_names$align_free), top_clones_names$align_free)
  ) + 
  labs(
    x     = "Alignment free clones",
    y     = "LC clones",
    title = glue("{HH}: Top 10 Alignment free VS LC clusters clones"),
    caption = glue("Minimum count: {min_count}")
  ) +
  theme_minimal(base_size = 10) + 
  theme(
    # axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "none"
  )

ggsave(glue("{outdir}/{HH}_align_free_VS_LC_clones.png"), width = 10, height = 11)



