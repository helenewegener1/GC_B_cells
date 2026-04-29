#!/bin/bash

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/48_GCtree
DATA_DIR=${WD}/fasta
OUT_DIR=${WD}/out
PLOT_DIR=${WD}/plot

# Get samples
sample_list=$(ls $DATA_DIR | cut -d "." -f1)

# # Run gctree
# for sample in $sample_list; do

  # sample=HH117_clone_nr_10_clone_2587_1
  sample=HH117_clone_nr_6_clone_1278_1

  # Make sample sepecific outdir
  OUT_DIR_SAMPLE=${OUT_DIR}/$sample
  mkdir -p $OUT_DIR_SAMPLE
  cd $OUT_DIR_SAMPLE
  
  # Clean up before new run 
  rm -f outfile outtree

  # Deduplication and sequence abundances
  deduplicate $DATA_DIR/${sample}.fasta \
  --root GL \
  --abundance_file abundances.csv \
  --idmapfile idmap.txt > deduplicated.phylip

  # Parsimony trees
  mkconfig deduplicated.phylip dnapars > dnapars.cfg
  dnapars < dnapars.cfg > dnapars.log
  
  # Make plotting sepecific outdir
  PLOT_DIR_SAMPLE=${PLOT_DIR}/$sample
  mkdir -p $PLOT_DIR_SAMPLE

  # Gctree Ranking
  xvfb-run -a gctree infer outfile abundances.csv --root GL --frame 1 --verbose --outbase $PLOT_DIR_SAMPLE/$sample

# done

