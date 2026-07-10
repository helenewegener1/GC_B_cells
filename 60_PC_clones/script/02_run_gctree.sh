#!/bin/bash
set -e

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/60_PC_clones
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
  
  # HH117-SILP-INF_clone_nr_10_clone_718_1      
  # HH117-SILP-INF_clone_nr_1_clone_8985_1      
  # HH117-SILP-INF_clone_nr_2_clone_455_1     
  # HH117-SILP-INF_clone_nr_3_clone_6469_1      
  # HH117-SILP-INF_clone_nr_4_clone_9677_1      
  # HH117-SILP-INF_clone_nr_5_clone_6896_1      
  # HH117-SILP-INF_clone_nr_6_clone_7587_1      
  # HH117-SILP-INF_clone_nr_7_clone_2275_1      
  # HH117-SILP-INF_clone_nr_8_clone_24_1    
  # HH117-SILP-INF_clone_nr_9_clone_9150_1 
  # HH117-SILP-nonINF_clone_nr_10_clone_8985_1
  # HH117-SILP-nonINF_clone_nr_1_clone_6469_1 
  # HH117-SILP-nonINF_clone_nr_2_clone_8879_1 
  # HH117-SILP-nonINF_clone_nr_3_clone_455_1 
  # HH117-SILP-nonINF_clone_nr_4_clone_8700_1 
  # HH117-SILP-nonINF_clone_nr_5_clone_9624_1 
  # # HH117-SILP-nonINF_clone_nr_6_clone_9347_1 
  # HH117-SILP-nonINF_clone_nr_7_clone_9677_1 
  # HH117-SILP-nonINF_clone_nr_8_clone_6961_1 
  # # HH117-SILP-nonINF_clone_nr_9_clone_8888_1
  HH119-COLP_clone_nr_10_clone_24340_1
  HH119-COLP_clone_nr_1_clone_23994_1
  HH119-COLP_clone_nr_2_clone_4885_1
  HH119-COLP_clone_nr_3_clone_2482_1
  HH119-COLP_clone_nr_4_clone_4747_1
  HH119-COLP_clone_nr_5_clone_28868_1
  HH119-COLP_clone_nr_6_clone_23862_1
  HH119-COLP_clone_nr_7_clone_14958_1
  HH119-COLP_clone_nr_8_clone_24368_1
  HH119-COLP_clone_nr_9_clone_20840_1
  HH119-SILP_clone_nr_10_clone_1800_1
  HH119-SILP_clone_nr_1_clone_27895_1
  HH119-SILP_clone_nr_2_clone_15879_1
  HH119-SILP_clone_nr_3_clone_4315_1
  HH119-SILP_clone_nr_4_clone_14406_1
  HH119-SILP_clone_nr_5_clone_23715_1
  HH119-SILP_clone_nr_6_clone_28184_1
  HH119-SILP_clone_nr_7_clone_20693_1
  HH119-SILP_clone_nr_8_clone_362_1
  HH119-SILP_clone_nr_9_clone_23831_2

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

