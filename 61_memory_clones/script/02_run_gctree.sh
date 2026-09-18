#!/bin/bash
set -e

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/61_memory_clones
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

  HH117_clone_nr_10_clone_546_1
  HH117_clone_nr_11_clone_8630_1
  HH117_clone_nr_12_clone_10216_1
  HH117_clone_nr_13_clone_1529_1
  HH117_clone_nr_14_clone_5261_1
  HH117_clone_nr_15_clone_8234_1
  HH117_clone_nr_16_clone_9610_1
  HH117_clone_nr_17_clone_9665_1
  HH117_clone_nr_18_clone_10097_1
  HH117_clone_nr_19_clone_1087_1
  HH117_clone_nr_1_clone_6696_1
  HH117_clone_nr_20_clone_1166_1
  HH117_clone_nr_2_clone_4458_1
  HH117_clone_nr_3_clone_2423_1
  HH117_clone_nr_4_clone_2714_1
  HH117_clone_nr_5_clone_9144_1
  HH117_clone_nr_6_clone_3949_1
  HH117_clone_nr_7_clone_4609_1
  HH117_clone_nr_8_clone_5169_1
  HH117_clone_nr_9_clone_5226_1
  HH119_clone_nr_10_clone_7523_1
  HH119_clone_nr_11_clone_14174_1
  HH119_clone_nr_12_clone_16124_1
  HH119_clone_nr_13_clone_9959_1
  HH119_clone_nr_14_clone_13246_1
  HH119_clone_nr_15_clone_10940_1
  HH119_clone_nr_16_clone_24453_1
  HH119_clone_nr_17_clone_10843_1
  HH119_clone_nr_18_clone_16262_1
  HH119_clone_nr_19_clone_18286_1
  HH119_clone_nr_1_clone_11791_1
  HH119_clone_nr_20_clone_25347_1
  HH119_clone_nr_2_clone_19239_1
  HH119_clone_nr_3_clone_17557_1
  HH119_clone_nr_4_clone_28841_1
  HH119_clone_nr_5_clone_5410_1
  HH119_clone_nr_6_clone_16134_1
  HH119_clone_nr_7_clone_18860_1
  HH119_clone_nr_8_clone_16295_1
  HH119_clone_nr_9_clone_17590_1
  HH151_clone_nr_1_clone_428_1
  HH151_clone_nr_2_clone_5497_1
  HH151_clone_nr_3_clone_5744_1
  HH151_clone_nr_4_clone_7292_1
  HH151_clone_nr_5_clone_8245_1
  HH153_clone_nr_1_clone_1128_1
  HH153_clone_nr_2_clone_219_1

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

