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

  sample=HH117_clone_nr_10_clone_2587_1

  # Make sample sepecific outdir
  OUT_DIR_SAMPLE=${OUT_DIR}/$sample
  mkdir -p $OUT_DIR_SAMPLE
  cd $OUT_DIR_SAMPLE
  
  # Clean up before new run 
  rm outfile outtree

  # Deduplication and sequence abundances
  deduplicate $DATA_DIR/${sample}.fasta \
  --root GL \
  --abundance_file abundances.csv \
  --idmapfile idmap.txt > deduplicated.phylip

  # Parsimony trees
  mkconfig deduplicated.phylip dnapars > dnapars.cfg
  dnapars < dnapars.cfg > dnapars.log

  # Gctree Ranking
  # gctree infer outfile abundances.csv --root GL --frame 1 --verbose
  xvfb-run -a gctree infer outfile abundances.csv --root GL --frame 1 --verbose --outbase $PLOT_DIR/$sample

# done

# 
# # Extra?
# export QT_QPA_PLATFORM=offscreen
# export XDG_RUNTIME_DIR=/tmp/runtime-runner
# export MPLBACKEND=agg
# export MPLCONFIGDIR=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/.matplotlib
# 
# gctree infer outfile abundances.csv --root GL --frame 1 --verbose

# export CONDA_PKGS_DIRS=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/.conda_pkgs
# 
# /services/tools/anaconda3/2024.06-1/bin/conda create --prefix /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/48_GCtree/script/gctree_env_310 python=3.10
