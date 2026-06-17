library(alakazam)
library(dplyr)
library(scales)
library(glue)

# Following: https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette/

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

HH <- "HH117"

resolve_LC <- readRDS(glue("45_immcantation/out/rds/{HH}_resolve_LC_3_definitions.rds"))
table(resolve_LC$locus)

df_heavy <- resolve_LC %>% filter(locus == "IGH")

nrow(df_heavy)

outdir <- "45_immcantation/plot/18_90_similarity/02_gene_usage/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Study gene usage
# ------------------------------------------------------------------------------

# Quantify usage at the gene level
gene <- countGenes(
  df_heavy, 
  gene="v_call", 
  groups="sample_clean", 
  mode="gene"
)

gene

# Quantify V family usage by sample
for (that_gene_call in c("v_call", "j_call")){
  # that_gene_call <- "v_call"
  
  # Plot V family usage by sample
  countGenes(df_heavy, gene=that_gene_call, groups="sample_clean", mode="family") %>% 
    ggplot(aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle(glue("{HH}: Family Usage of {that_gene_call}")) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    geom_point(aes(color=sample_clean), size=5, alpha=0.8, position=position_dodge(width = 0.1))
  
  ggsave(glue("{outdir}/{HH}_{that_gene_call}.png"))
  
}


# Plot V gene usage in the IGHVX family by sample
that_gene <- "IGHV3"

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

ggsave(glue("{outdir}/{HH}_{that_gene}.png"))


# ------------------------------------------------------------------------------
# Study gene usage - more groups - count
# ------------------------------------------------------------------------------

# Quantify V family clonal usage by sample and isotype
that_gene_call <- "v_call"

family <- countGenes(df_heavy, gene=that_gene_call, groups=c("sample_clean", "c_call"), 
                     clone="clone_subgroup_id_90_similarity", mode="family")
head(family, n=4)

# Subset to IGHM and IGHG for plotting
# Plot V family clonal usage by sample and isotype
family %>% 
  # filter(c_call %in% c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4")) %>% 
  ggplot(aes(x=gene, y=clone_freq)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(
    axis.text.x=element_text(angle=45, hjust=1, vjust=1), 
    legend.position = "bottom"
  ) +
  ylab("Percent of repertoire") +
  xlab("") +
  geom_point(aes(color=sample_clean), size=5, alpha=0.8, position=position_dodge(width = 0.1)) +
  facet_grid(. ~ c_call)

ggsave(glue("{outdir}/{HH}_{that_gene_call}_vs_isotype.png"), width = 15)

# ------------------------------------------------------------------------------
# IGHV4-34 
# ------------------------------------------------------------------------------

this_gene <- "IGHV4-34"

df_heavy %>% filter(v_call_no_allele == this_gene) %>% count(sample_clean, c_call, sort = T) %>% 
  ggplot(aes(x = c_call, y = n, color = sample_clean)) + 
  geom_point(size = 4, alpha = 0.7) + 
  theme_bw() + 
  labs(
    title = glue("N cells with {this_gene} across isotypes")
  )

ggsave(glue("{outdir}/{HH}_{this_gene}_vs_isotype.png"))

df_heavy %>% filter(v_call_no_allele == this_gene) %>% count(sample_clean, L1_annotation, sort = T) %>% 
  ggplot(aes(x = L1_annotation, y = n, color = sample_clean)) + 
  geom_point(size = 4, alpha = 0.7) + 
  theme_bw() + 
  labs(
    title = glue("N cells with {this_gene} across isotypes")
  )

ggsave(glue("{outdir}/{HH}_{this_gene}_vs_L1_annotation.png"))

# ------------------------------------------------------------------------------
# Study gene usage - more groups - copy number
# ------------------------------------------------------------------------------

# # Quantify V family clonal usage by sample and isotype
# that_gene_call <- "v_call"
# 
# family <- countGenes(df_heavy, gene=that_gene_call, groups=c("sample_clean", "c_call"),
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