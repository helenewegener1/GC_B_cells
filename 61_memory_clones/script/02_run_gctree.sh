#!/bin/bash
set -e

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/61_PC_memory
# DATA_DIR=${WD}/fasta/GC_clones
# OUT_DIR=${WD}/out
# PLOT_DIR=${WD}/plot

# DATA_DIR=${WD}/fasta/gmm_threshold_GC_clones
# OUT_DIR=${WD}/out_gmm_threshold
# PLOT_DIR=${WD}/plot_gmm_threshold

DATA_DIR=${WD}/fasta
OUT_DIR=${WD}/out
PLOT_DIR=${WD}/plot

# sample_list=$(ls $DATA_DIR | cut -d "." -f1)
sample_list=(

  HH117_clone_nr_1_clone_6680_1
  HH117_clone_nr_2_clone_1788_1
  HH117_clone_nr_3_clone_4448_1
  HH117_clone_nr_4_clone_2392_1
  HH117_clone_nr_5_clone_2703_1
  HH117_clone_nr_6_clone_9105_1
  HH117_clone_nr_7_clone_3932_1
  HH117_clone_nr_8_clone_4601_1
  HH117_clone_nr_9_clone_5150_1
  HH117_clone_nr_10_clone_5208_1
  HH117_clone_nr_11_clone_522_1
  HH117_clone_nr_12_clone_8608_1
  HH117_clone_nr_13_clone_10202_1
  HH117_clone_nr_14_clone_10223_1
  HH117_clone_nr_15_clone_1502_1
  HH117_clone_nr_16_clone_5243_1
  HH117_clone_nr_17_clone_8213_1
  HH117_clone_nr_18_clone_9573_1
  HH117_clone_nr_19_clone_9628_1
  HH117_clone_nr_20_clone_10056_1
  HH119_clone_nr_1_clone_19030_1
  HH119_clone_nr_2_clone_11649_1
  HH119_clone_nr_3_clone_17493_1
  HH119_clone_nr_4_clone_21605_1
  HH119_clone_nr_5_clone_28521_1
  HH119_clone_nr_6_clone_5519_1
  HH119_clone_nr_7_clone_16678_1
  HH119_clone_nr_8_clone_16073_1
  HH119_clone_nr_9_clone_17526_1
  HH119_clone_nr_10_clone_18616_1
  HH119_clone_nr_11_clone_7441_1
  HH119_clone_nr_12_clone_13948_1
  HH119_clone_nr_13_clone_16668_1
  HH119_clone_nr_14_clone_9837_1
  HH119_clone_nr_15_clone_10865_1
  HH119_clone_nr_16_clone_24292_1
  HH119_clone_nr_17_clone_13034_1
  HH119_clone_nr_18_clone_10721_1
  HH119_clone_nr_19_clone_16040_1
  HH119_clone_nr_20_clone_17902_1

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

