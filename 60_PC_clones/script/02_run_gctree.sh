#!/bin/bash
set -e

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/60_PC_clones

DATA_DIR=${WD}/fasta
OUT_DIR=${WD}/out
PLOT_DIR=${WD}/plot

# sample_list=$(ls $DATA_DIR | cut -d "." -f1)
sample_list=(
  
  HH117-SILP-INF_clone_nr_1_clone_8843_1
  HH117-SILP-INF_clone_nr_2_clone_9213_1
  HH117-SILP-INF_clone_nr_3_clone_6758_1
  HH117-SILP-INF_clone_nr_4_clone_5794_1
  HH117-SILP-INF_clone_nr_5_clone_4207_1
  HH117-SILP-INF_clone_nr_6_clone_9790_1
  HH117-SILP-INF_clone_nr_7_clone_6895_1
  HH117-SILP-INF_clone_nr_8_clone_7500_1
  HH117-SILP-INF_clone_nr_9_clone_2015_1
  HH117-SILP-INF_clone_nr_10_clone_142_1
  HH117-SILP-nonINF_clone_nr_1_clone_9213_1
  HH117-SILP-nonINF_clone_nr_2_clone_6758_1
  HH117-SILP-nonINF_clone_nr_3_clone_8750_1
  HH117-SILP-nonINF_clone_nr_4_clone_8565_1
  HH117-SILP-nonINF_clone_nr_5_clone_9492_1
  HH117-SILP-nonINF_clone_nr_6_clone_9232_1
  HH117-SILP-nonINF_clone_nr_7_clone_4207_1
  HH117-SILP-nonINF_clone_nr_8_clone_5794_1
  HH117-SILP-nonINF_clone_nr_9_clone_6938_1
  HH117-SILP-nonINF_clone_nr_10_clone_8843_1
  HH119-COLP_clone_nr_1_clone_24419_1
  HH119-COLP_clone_nr_2_clone_4664_1
  HH119-COLP_clone_nr_3_clone_4528_1
  HH119-COLP_clone_nr_4_clone_9_1
  HH119-COLP_clone_nr_5_clone_17441_1
  HH119-COLP_clone_nr_6_clone_2415_1
  HH119-COLP_clone_nr_7_clone_14852_1
  HH119-COLP_clone_nr_8_clone_24807_1
  HH119-COLP_clone_nr_9_clone_24782_1
  HH119-COLP_clone_nr_10_clone_28532_1
  HH119-SILP_clone_nr_1_clone_27931_1
  HH119-SILP_clone_nr_2_clone_15801_1
  HH119-SILP_clone_nr_3_clone_4112_1
  HH119-SILP_clone_nr_4_clone_14426_1
  HH119-SILP_clone_nr_5_clone_24269_1
  HH119-SILP_clone_nr_6_clone_28233_1
  HH119-SILP_clone_nr_7_clone_344_1
  HH119-SILP_clone_nr_8_clone_24246_2
  HH119-SILP_clone_nr_9_clone_19014_1
  HH119-SILP_clone_nr_10_clone_11471_1

)

# # Run gctree
# for sample in $sample_list; do
for sample in "${sample_list[@]}"; do

  echo "Processing $sample..."

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

