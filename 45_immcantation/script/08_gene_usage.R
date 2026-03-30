library(alakazam)
library(dplyr)
library(scales)
library(glue)

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

spec_clones_vj <- readRDS("45_immcantation/out/rds/spec_clones_vj.rds")

HH <- "HH119"
HH_spec_clones_vj <- spec_clones_vj[[HH]]

nrow(HH_spec_clones_vj)

# ------------------------------------------------------------------------------
# Study gene usage
# ------------------------------------------------------------------------------

# Quantify usage at the gene level
gene <- countGenes(
  HH_spec_clones_vj, 
  gene="v_call", 
  groups="sample_clean", 
  mode="gene"
)

gene

# Plot V gene usage in the IGHVX family by sample
that_gene <- "IGHV4"

gene %>%
  mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
  filter(getFamily(gene) == that_gene) %>% 
  ggplot(aes(x=gene, y=seq_freq)) +
  theme_bw() +
  ggtitle(glue("{HH}: {that_gene} Usage")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  geom_point(aes(color=sample_clean), size=5, alpha=0.8, position=position_dodge(width = 0.1))


# Quantify V family usage by sample
that_gene_call <- "v_call"

# Plot V family usage by sample
countGenes(HH_spec_clones_vj, gene=that_gene_call, groups="sample_clean", mode="family") %>% 
  ggplot(aes(x=gene, y=seq_freq)) +
  theme_bw() +
  ggtitle(glue("{HH}: Family Usage of {that_gene_call}")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  geom_point(aes(color=sample_clean), size=5, alpha=0.8, position=position_dodge(width = 0.1))


# ------------------------------------------------------------------------------
# Study gene usage - more groups - count
# ------------------------------------------------------------------------------

# Quantify V family clonal usage by sample and isotype
that_gene_call <- "v_call"

family <- countGenes(HH_spec_clones_vj, gene=that_gene_call, groups=c("sample_clean", "c_call"), 
                     clone="clone_id", mode="family")
head(family, n=4)

# Subset to IGHM and IGHG for plotting
# Plot V family clonal usage by sample and isotype
family %>% 
  filter(c_call %in% c("IGHM", "IGHG1")) %>% 
  ggplot(aes(x=gene, y=clone_freq)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  geom_point(aes(color=sample_clean), size=5, alpha=0.8, position=position_dodge(width = 0.1)) +
  facet_grid(. ~ c_call)

# ------------------------------------------------------------------------------
# Study gene usage - more groups - copy number
# ------------------------------------------------------------------------------

# # Quantify V family clonal usage by sample and isotype
# that_gene_call <- "v_call"
# 
# family <- countGenes(HH_spec_clones_vj, gene=that_gene_call, groups=c("sample_clean", "c_call"), 
#                      copy="duplicate_count")
# head(family, n=4)
# 
# # Subset to IGHM and IGHG for plotting
# family <- filter(family, c_call %in% c("IGHM", "IGHG"))
# # Plot V family copy abundance by sample and isotype
# g4 <- ggplot(family, aes(x=gene, y=copy_freq)) +
#   theme_bw() +
#   ggtitle("Copy Number") +
#   theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
#   ylab("Percent of repertoire") +
#   xlab("") +
#   scale_y_continuous(labels=percent) +
#   scale_color_brewer(palette="Set1") +
#   geom_point(aes(color=sample_id), size=5, alpha=0.8, position=position_dodge(width = 0.1)) +
#   facet_grid(. ~ c_call)
# plot(g4)