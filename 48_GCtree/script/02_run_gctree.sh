#!/bin/bash
set -e

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/48_GCtree
# DATA_DIR=${WD}/fasta/GC_clones
# OUT_DIR=${WD}/out
# PLOT_DIR=${WD}/plot

# # Get samples
# # sample_list=$(ls $DATA_DIR | cut -d "." -f1)
# sample_list=(
#   HH119_clone_nr_1_clone_4516_1
#   # HH117_clone_nr_4_clone_735_1
#   # HH119_clone_nr_8_clone_2062_1
#   # HH117_clone_nr_5_clone_617_1
#   # HH119_clone_nr_9_clone_2388_1
#   # HH117_clone_nr_6_clone_1012_1
#   # HH119_clone_nr_3_clone_2466_1
#   # HH119_clone_nr_2_clone_5791_1
# )

DATA_DIR=${WD}/fasta/gmm_threshold_GC_clones
OUT_DIR=${WD}/out_gmm_threshold
PLOT_DIR=${WD}/plot_gmm_threshold

# sample_list=$(ls $DATA_DIR | cut -d "." -f1)
sample_list=(
  
  HH117_clone_nr_10_clone_2169_1
  HH117_clone_nr_9_clone_1910_1
  HH117_clone_nr_8_clone_6115_1
  HH117_clone_nr_7_clone_5941_1
  HH117_clone_nr_6_clone_1320_1
  HH117_clone_nr_5_clone_2301_1
  HH117_clone_nr_4_clone_3709_1
  HH117_clone_nr_3_clone_1849_1
  HH117_clone_nr_2_clone_2628_1
  HH117_clone_nr_1_clone_4221_1
  
  HH119_clone_nr_10_clone_7913_1
  HH119_clone_nr_9_clone_25158_1
  HH119_clone_nr_8_clone_3869_1
  HH119_clone_nr_7_clone_8372_1
  HH119_clone_nr_6_clone_23124_1
  HH119_clone_nr_5_clone_15287_1
  HH119_clone_nr_4_clone_12120_1
  HH119_clone_nr_3_clone_28075_1
  HH119_clone_nr_2_clone_22387_1
  # HH119_clone_nr_1_clone_22355_1
  
)

# # Run gctree
# for sample in $sample_list; do
for sample in "${sample_list[@]}"; do

  echo "Processing $sample..."

  # sample=HH117_clone_nr_10_clone_2587_1
  # sample=HH117_clone_nr_6_clone_1278_1

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
  xvfb-run -a gctree infer outfile abundances.csv --root GL --frame 1 --verbose --outbase $PLOT_DIR_SAMPLE/${sample}

  echo "Processing of $sample is complete!"

done

