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

DATA_DIR=${WD}/fasta/90_similarity
OUT_DIR=${WD}/out_90_similarity
PLOT_DIR=${WD}/plot_90_similarity

# sample_list=$(ls $DATA_DIR | cut -d "." -f1)
sample_list=(

  HH117_clone_nr_10_clone_2457_1, 
  HH117_clone_nr_1_clone_500_1, 
  HH117_clone_nr_2_clone_2885_1, 
  HH117_clone_nr_3_clone_2018_1, 
  HH117_clone_nr_4_clone_3953_1, 
  HH117_clone_nr_5_clone_2588_1, 
  HH117_clone_nr_6_clone_1571_1, 
  HH117_clone_nr_7_clone_5554_1, 
  HH117_clone_nr_8_clone_5725_1, 
  HH117_clone_nr_9_clone_2211_1,
  
  HH119_clone_nr_10_clone_23719_1, 
  # HH119_clone_nr_1_clone_20693_1, 
  # HH119_clone_nr_2_clone_27035_1, 
  HH119_clone_nr_3_clone_11978_1, 
  HH119_clone_nr_4_clone_7879_1, 
  HH119_clone_nr_5_clone_14912_1, 
  HH119_clone_nr_6_clone_21402_1, 
  HH119_clone_nr_7_clone_8286_1, 
  HH119_clone_nr_8_clone_8625_1, 
  HH119_clone_nr_9_clone_3913_1

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

