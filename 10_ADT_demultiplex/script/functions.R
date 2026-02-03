
rules_text <- 'top_signal < 1  --> "Negative"
n_above_1 == 1 --> "Singlet",
n_above_1 > 1 & ratio_12 < 0.5 --> "Singlet"
TRUE --> "Doublet"
'

# plot
plot_ADT_cells <- function(df, cell_nr){
  
  # cell_nr <- 18
  cell_name <- colnames(seurat_obj)[cell_nr]
  
  label_y_text <- df %>% 
    filter(Cell == cell_name) %>% select(top_signal) %>% pull() %>% unique() / 2
  
  p <- df %>%
    filter(Cell == cell_name) %>%
    ggplot(aes(x = reorder(Fol, -Count), y = Count)) +
    geom_col() +
    geom_label(aes(label = round(Count, 2)), size = 2) +
    geom_label(
      data = df %>% filter(Cell == cell_name) %>% filter(!is.na(Count)),
      aes(
        x = label_x,
        y = label_y,
        label = round(ratio, 2)
      ),
      inherit.aes = FALSE,
      size = 3
    ) +
    labs(
      title = glue("Cell {cell_nr}: ADT counts corrected by thershold"),
      y = "Normalized counts (counts divided by zero-point-count)",
      x = "Fol",
      subtitle = "Normalized count of 1: same value as zero-point",
      # subtitle = "Log-ratio of 0: same value\nLog-ratio of -1: half the value",
      caption = glue("Dominant ADT: {df$dominant_ADT_full_ID[df$Cell == cell_name][1]}
                     HG ADT: {df$HG_ADT_full_ID[df$Cell == cell_name][1]}
                     ")
    ) +
    annotate("text", x=13, y=label_y_text, label= rules_text) +
    theme_bw()
  
  ggsave(filename = glue("{outdir_plot}/{cell_nr}.png"),
         p,
         width = 13,
         height = 7)
  
  return(p)
  
}
