library(alakazam)
library(dplyr)
library(ggplot2)
library(shazam)
library(glue)

# Following: https://shazam.readthedocs.io/en/stable/vignettes/Mutation-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

HH <- "HH119"

resolve_LC <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds"))
table(resolve_LC$locus)

df_heavy <- resolve_LC %>% filter(locus == "IGH")

nrow(df_heavy)

outdir <- "45_immcantation/plot/18_90_similarity/03_mutational_load/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Calculate the counts and frequencies of mutations over the entire sequence
# ------------------------------------------------------------------------------

# Mutations can be calculated as either counts or frequencies, 
# divided into replacement (R) and silent (S) mutations, 
# and subset into FWR and CDR specific mutations

# ----------------------
# Calculate R and S mutation counts
# ----------------------
db_obs <- observedMutations(df_heavy, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
# Show new mutation count columns
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>% head(n=4)

# ----------------------
# Calculate R and S mutation frequencies
# ----------------------
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>% head(n=4)

# ----------------------
# Specifying the combine=TRUE argument will aggregate all mutation columns into a single value.
# Calculate combined R and S mutation frequencies
# ----------------------
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            combine=TRUE,
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq")) %>% head(n=4)

# ----------------------
# Plot
# ----------------------

db_obs %>% 
  filter(!is.na(c_call)) %>% 
  ggplot(aes(x=c_call, y=mu_freq, fill=c_call)) +
  geom_boxplot() + 
  labs(title="Total mutations", 
       x="Isotype", y="Mutation frequency") +
  theme_bw() + 
  facet_wrap(vars(sample_clean)) + 
  theme(legend.position = "none")

ggsave(glue("45_immcantation/plot/08_mutational_load/{HH}_total_mutations_across_isotypes_and_sample_clean.png"), width = 15, height = 8)

# Continues...