library(scoper)
library(dplyr)

# Following this SCOPer tutorial: 
# https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

spec_clones_novj <- readRDS("45_immcantation/out/rds/spec_clones_novj.rds")
spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")
hier <- readRDS("45_immcantation/out/rds/hier_clones.rds")

patients <- names(hier)

# ------------------------------------------------------------------------------
# Compare spetral to hierical clustering
# ------------------------------------------------------------------------------

HH <- "HH119"

HH_spec_clones_vj <- spec_clones_vj[[HH]]@db
HH_spec_clones_novj <- spec_clones_novj[[HH]]@db
HH_hier <- hier %>% filter(patient_id == HH)

# Not in the same order...
table(HH_spec_clones_vj$cell_id %in% HH_spec_clones_novj$cell_id) 
table(HH_spec_clones_vj$cell_id == HH_spec_clones_novj$cell_id) 
table(HH_spec_clones_vj$cell_id %in% HH_hier$cell_id) 
table(HH_spec_clones_vj$cell_id == HH_hier$cell_id) 

length(HH_spec_clones_vj$cell_id)
length(HH_spec_clones_novj$cell_id)
length(HH_hier$cell_id)

# reorder novj to match the order of vj
HH_spec_clones_novj <- HH_spec_clones_novj[match(HH_spec_clones_vj$cell_id, HH_spec_clones_novj$cell_id), ]
HH_hier <- HH_hier[match(HH_spec_clones_vj$cell_id, HH_hier$cell_id), ]

# Compare largest clones
HH_hier %>% count(clone_id, v_call, j_call, sort = TRUE) # HH119 IGHV4-34 clone is top but rest are different from the spetral clones

HH_spec_clones_novj %>% count(clone_id, v_call, j_call, sort = TRUE)

HH_spec_clones_vj %>% count(clone_id, v_call, j_call, sort = TRUE)

# Compare N unique clones 
HH_hier$clone_id %>% unique() %>% length() # More unique clones 
HH_spec_clones_novj$clone_id %>% unique() %>% length()
HH_spec_clones_vj$clone_id %>% unique() %>% length()

# Compare N clones of 1 cell
HH_hier$clone_id %>% table() %>% as.data.frame() %>% filter(Freq == 1) %>% nrow() # More clones of the size of 1 cell
HH_spec_clones_novj$clone_id %>% table() %>% as.data.frame() %>% filter(Freq == 1) %>% nrow()
HH_spec_clones_vj$clone_id %>% table() %>% as.data.frame() %>% filter(Freq == 1) %>% nrow()

# Compare N clones of 2 cells
HH_hier$clone_id %>% table() %>% as.data.frame() %>% filter(Freq == 2) %>% nrow()
HH_spec_clones_novj$clone_id %>% table() %>% as.data.frame() %>% filter(Freq == 2) %>% nrow()
HH_spec_clones_vj$clone_id %>% table() %>% as.data.frame() %>% filter(Freq == 2) %>% nrow()

# Claude compare
comparison <- HH_spec_clones_vj %>%
  select(cell_id, clone_vj = clone_id, v_call, j_call) %>%
  left_join(
    HH_spec_clones_novj %>% select(cell_id, clone_novj = clone_id),
    by = "cell_id"
  ) %>%
  left_join(
    HH_hier %>% select(cell_id, clone_hier = clone_id),
    by = "cell_id"
  )

# Do cells in the same vj clone also end up in the same novj and hier clone?
comparison %>%
  group_by(clone_vj) %>%
  summarise(
    n_cells = n(),
    n_novj_clones = n_distinct(clone_novj), # N unique clones that clone_vj splits into hence 1 = perfect agreement
    n_hier_clones = n_distinct(clone_hier), # N unique clones that clone_hier splits into hence 1 = perfect agreement
    dominant_v = names(sort(table(v_call), decreasing = TRUE))[1],
    dominant_j = names(sort(table(j_call), decreasing = TRUE))[1]
  ) %>%
  arrange(desc(n_cells))

# Do cells in the same novj clone also end up in the same vj and hier clone?
comparison %>%
  group_by(clone_novj) %>%
  summarise(
    n_cells = n(),
    n_vj_clones = n_distinct(clone_vj), 
    n_hier_clones = n_distinct(clone_hier), 
    dominant_v = names(sort(table(v_call), decreasing = TRUE))[1],
    dominant_j = names(sort(table(j_call), decreasing = TRUE))[1]
  ) %>%
  arrange(desc(n_cells))

# Do cells in the same hier clone also end up in the same vj and novj clone?
# Hier makes more smaller clones. 
comparison %>%
  group_by(clone_hier) %>%
  summarise(
    n_cells = n(),
    n_vj_clones = n_distinct(clone_vj),
    n_novj_clones = n_distinct(clone_novj), 
    dominant_v = names(sort(table(v_call), decreasing = TRUE))[1],
    dominant_j = names(sort(table(j_call), decreasing = TRUE))[1]
  ) %>%
  arrange(desc(n_cells))

# Junctions

# Figure out which on to go with. 
# hierical clustering done with mean threshold for patients - should this had been done seperatly? Test
# Continue with Immcantation Tutorials


