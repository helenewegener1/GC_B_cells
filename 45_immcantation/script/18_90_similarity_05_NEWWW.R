library(alakazam)
library(dplyr)
library(scales)
library(glue)

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

HH <- "HH119"

resolve_LC <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds"))
table(resolve_LC$locus)

df_heavy <- resolve_LC %>% filter(locus == "IGH")

nrow(df_heavy)

outdir <- "45_immcantation/plot/18_90_similarity/05_V4_34_zoom/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Expression of V4-34 gene 
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
