clones <- read.table("50_MiXCR/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2.clones.tsv", 
                     header=TRUE, sep="\t")

cell_clones <- read.table("50_MiXCR/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2.clone.groups_IG.tsv", 
                          header=TRUE, sep="\t")

# How many cells have both heavy and light chain?
table(!is.na(cell_clones$IGH.primary.aaSeqCDR3) & 
        (!is.na(cell_clones$IGL.primary.aaSeqCDR3) | !is.na(cell_clones$IGK.primary.aaSeqCDR3)))

# Isotype distribution
table(cell_clones$IGH.primary.bestCHit)

# Chain usage (IGL vs IGK)
table(!is.na(cell_clones$IGL.primary.aaSeqCDR3), 
      !is.na(cell_clones$IGK.primary.aaSeqCDR3))


cell_clones$IGH.primary.bestVHit
cell_clones$IGL.primary.bestVHit
cell_clones$IGL.primary.bestVHit


clones$allVHitsWithScore


cell_clones[which.max(cell_clones$groupReadCount), ]

cell_clones %>%
  arrange(desc(groupReadCount)) %>%
  select(cellGroup, groupReadCount, groupUniqueCellCount, 
         IGH.primary.bestVHit, IGH.primary.bestJHit, IGH.primary.bestCHit,
         IGH.primary.aaSeqCDR3, 
         IGL.primary.aaSeqCDR3, IGK.primary.aaSeqCDR3) %>%
  head(10)

# Ranked by read count
clones %>%
  arrange(desc(readCount)) %>%
  select(cloneId, topChains, readCount, uniqueMoleculeCount,
         allVHitsWithScore, allJHitsWithScore, aaSeqCDR3) %>%
  head(10)
