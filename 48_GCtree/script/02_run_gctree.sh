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

# DATA_DIR=${WD}/fasta/gmm_threshold_GC_clones
# OUT_DIR=${WD}/out_gmm_threshold
# PLOT_DIR=${WD}/plot_gmm_threshold

DATA_DIR=${WD}/fasta/
OUT_DIR=${WD}/out/
PLOT_DIR=${WD}/plot/

# This works to run all clones
# sample_list=$(ls $DATA_DIR | cut -d "." -f1)
# echo $sample_list

sample_list=(

  HH117_clone_nr_1_clone_4252_1
  HH117_clone_nr_2_clone_2654_1
  HH117_clone_nr_3_clone_1854_1
  HH117_clone_nr_4_clone_3755_1
  HH117_clone_nr_5_clone_2335_1
  HH117_clone_nr_6_clone_1307_1
  HH117_clone_nr_7_clone_5945_1
  HH117_clone_nr_8_clone_6111_1
  HH117_clone_nr_9_clone_1922_1
  HH117_clone_nr_10_clone_2197_1
  HH117_clone_nr_11_clone_2835_1
  HH117_clone_nr_12_clone_5340_1
  HH117_clone_nr_13_clone_6254_1
  HH117_clone_nr_14_clone_6680_1
  HH117_clone_nr_15_clone_688_1
  HH117_clone_nr_16_clone_10074_1
  HH117_clone_nr_17_clone_9559_1
  HH117_clone_nr_18_clone_5832_1
  HH117_clone_nr_19_clone_6058_1
  HH117_clone_nr_20_clone_7375_1
  
  # HH119_clone_nr_1_clone_21605_1
  HH119_clone_nr_2_clone_27055_1
  HH119_clone_nr_3_clone_11748_1
  HH119_clone_nr_4_clone_7647_1
  HH119_clone_nr_5_clone_14825_1
  HH119_clone_nr_6_clone_22333_1
  HH119_clone_nr_7_clone_8071_1
  HH119_clone_nr_8_clone_8404_1
  HH119_clone_nr_9_clone_3761_1
  HH119_clone_nr_10_clone_24273_1
  HH119_clone_nr_11_clone_8356_1
  HH119_clone_nr_12_clone_180_1
  HH119_clone_nr_13_clone_28685_1
  HH119_clone_nr_14_clone_15980_1
  HH119_clone_nr_15_clone_10216_1
  HH119_clone_nr_16_clone_10133_1
  HH119_clone_nr_17_clone_8759_1
  HH119_clone_nr_18_clone_25753_1
  HH119_clone_nr_19_clone_26371_1
  HH119_clone_nr_20_clone_10141_1

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

