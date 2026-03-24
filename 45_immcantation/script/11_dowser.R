library(dowser)
library(alakazam)
library(dplyr)
library(ggplot2)
library(shazam)
library(glue)

# Following: https://dowser.readthedocs.io/en/stable/vignettes/Germlines-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")

HH <- "HH119"
HH_spec_clones_vj <- spec_clones_vj[[HH]]

nrow(HH_spec_clones_vj)

# Check that neccessary columns are present
HH_spec_clones_vj$sequence_alignment %>% head()
HH_spec_clones_vj$germline_alignment_d_mask %>% head()

# ------------------------------------------------------------------------------
# Construct clonal germlines
# ------------------------------------------------------------------------------

# Download (IMGT) germline reference database
# In terminal
# git clone https://github.com/immcantation
# cd immcantation/scripts
# fetch_imgtdb.sh -o germlines # Run script to obtain IMGT gapped sequences

# Read in IMGT-gapped sequences
# references <- readIMGT(dir = "../packages/immcantation/scripts/germlines/human/vdj")

# remove germline alignment columns for this example
db <- select(HH_spec_clones_vj, -"germline_alignment", -"germline_alignment_d_mask")

# Reconstruct germline sequences
HH_spec_clones_vj <- createGermlines(db, references, nproc=1)

# Check germline of first row
HH_spec_clones_vj$germline_alignment_d_mask[1]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Build Lineage Trees
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Format clones
# ------------------------------------------------------------------------------

# Top clone
clone_nr <- 5
clone <- HH_spec_clones_vj %>% count(clone_id, sort = TRUE) %>% head(1) %>% pull(clone_id)

# Subset data for this example
HH_spec_clones_vj <- HH_spec_clones_vj[HH_spec_clones_vj$clone_id == clone,]

# Add count for identical clones - used for tipsize of tree
HH_spec_clones_vj <- HH_spec_clones_vj %>%
  group_by(clone_id, sequence_alignment) %>%
  mutate(n_identical = n()) %>%
  ungroup()

# Process example data using default settings
clones <- formatClones(HH_spec_clones_vj, text_fields = c("c_call", "celltype_broad"), num_fields=c("n_identical"))

print(clones)

# ------------------------------------------------------------------------------
# Build trees
# ------------------------------------------------------------------------------

# Maximum parsimony trees using phangorn.
# clones <- getTrees(clones, nproc=1)

# Maximum likelihood trees using phangorn
clones <- getTrees(clones, build="pml")

# ------------------------------------------------------------------------------
# Plot tree
# ------------------------------------------------------------------------------

# plotTrees(clones)

plotTrees(
  clones, 
  tips="c_call", 
  tipsize="n_identical",
  title = FALSE
)[[1]] + 
  plot_annotation(
    title = glue("{HH}: Clone number {clone_nr}")
  )

ggsave(glue("45_immcantation/plot/11_dowser/{HH}_dowser_tree_clone_{clone_nr}_c_call.png"), width = 15, height = 25)

plotTrees(
  clones, 
  tips="celltype_broad", 
  tipsize="n_identical",
  title = FALSE
)[[1]] + 
  plot_annotation(
    title = glue("{HH}: Clone number {clone_nr}")
  )

ggsave(glue("45_immcantation/plot/11_dowser/{HH}_dowser_tree_clone_{clone_nr}_celltype.png"), width = 15, height = 25)

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
