

seurat_obj

seurat_obj$scDblFinder.score
seurat_obj$orig.ident

Idents(seurat_obj) <- "orig.ident"

seurat_obj[[]] %>% 
  ggplot(aes(y = scDblFinder.score, x = Phase)) + 
  geom_violin() + 
  geom_point(aes(color = scDblFinder.class), alpha = 0.5) +
  # scale_color_manual(values = c("grey", "red", "blue")) + 
  theme_bw() 


# VlnPlot(seurat_obj, features = "scDblFinder.score", split.by = )
