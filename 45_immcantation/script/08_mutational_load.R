library(alakazam)
library(dplyr)
library(ggplot2)
library(shazam)
library(glue)

# Following: https://shazam.readthedocs.io/en/stable/vignettes/Mutation-Vignette/

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

# Subset data
db <- subset(HH_spec_clones_vj)# & sample_clean == "")

# ------------------------------------------------------------------------------
# Calculate the counts and frequencies of mutations over the entire sequence
# ------------------------------------------------------------------------------

# Mutations can be calculated as either counts or frequencies, 
# divided into replacement (R) and silent (S) mutations, 
# and subset into FWR and CDR specific mutations

# ----------------------
# Calculate R and S mutation counts
# ----------------------
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
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