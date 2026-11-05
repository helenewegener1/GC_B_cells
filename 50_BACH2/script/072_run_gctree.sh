#!/bin/bash
set -e

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_BACH2
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
  
  # "HH117_clone_nr_1_clone_1881_1"
  # "HH117_clone_nr_2_clone_1951_1"
  # "HH117_clone_nr_3_clone_2667_1"
  # "HH117_clone_nr_4_clone_4261_1"
  # "HH117_clone_nr_5_clone_5139_1"
  # "HH117_clone_nr_6_clone_5959_1"
  # "HH117_clone_nr_7_clone_4066_1"
  # "HH117_clone_nr_8_clone_5608_1"
  # "HH117_clone_nr_9_clone_5834_1"
  # "HH117_clone_nr_10_clone_6270_1"
  # "HH117_clone_nr_11_clone_715_1"
  # "HH117_clone_nr_12_clone_7509_1"
  # "HH117_clone_nr_13_clone_9817_1"
  # 
  # "HH119_clone_nr_1_clone_27383_1"
  # "HH119_clone_nr_2_clone_4494_1"
  # "HH119_clone_nr_3_clone_29181_1"
  # "HH119_clone_nr_4_clone_8598_1"
  # "HH119_clone_nr_5_clone_15039_1"
  # "HH119_clone_nr_6_clone_26043_1"
  # "HH119_clone_nr_7_clone_8506_1"
  # "HH119_clone_nr_8_clone_22295_1"
  # "HH119_clone_nr_9_clone_17585_1"
  # "HH119_clone_nr_10_clone_8457_1"
  # "HH119_clone_nr_11_clone_22776_1"
  # "HH119_clone_nr_12_clone_26698_1"
  # "HH119_clone_nr_13_clone_28739_1"
  # "HH119_clone_nr_14_clone_3804_1"
  # "HH119_clone_nr_15_clone_8177_1"
  # "HH119_clone_nr_16_clone_24431_1"
  # "HH119_clone_nr_17_clone_28674_1"
  # "HH119_clone_nr_18_clone_1830_1"
  # "HH119_clone_nr_19_clone_24483_1"
  # "HH119_clone_nr_20_clone_8591_1"
  # "HH119_clone_nr_21_clone_903_1"
  # "HH119_clone_nr_22_clone_19562_1"
  # "HH119_clone_nr_23_clone_22252_1"
  # "HH119_clone_nr_24_clone_23494_1"
  # "HH119_clone_nr_25_clone_25333_1"
  # "HH119_clone_nr_26_clone_10201_1"
  # "HH119_clone_nr_27_clone_17585_2"
  # "HH119_clone_nr_28_clone_17807_1"
  # "HH119_clone_nr_29_clone_26121_1"
  # "HH119_clone_nr_30_clone_2812_1"
  # "HH119_clone_nr_31_clone_5037_1"
  # "HH119_clone_nr_32_clone_5652_1"
  # "HH119_clone_nr_33_clone_5836_1"
  # "HH119_clone_nr_34_clone_8685_1"
  # 
  # "HH151_clone_nr_1_clone_3084_1"
  # "HH151_clone_nr_2_clone_4747_1"
  # "HH151_clone_nr_3_clone_6909_1"

  "HH153_clone_nr_1_clone_188_1"
  "HH153_clone_nr_2_clone_3105_1"
  "HH153_clone_nr_3_clone_2204_1"
  "HH153_clone_nr_4_clone_125_1"
  "HH153_clone_nr_5_clone_3058_1"
  "HH153_clone_nr_6_clone_2209_1"
  "HH153_clone_nr_7_clone_3107_1"
  "HH153_clone_nr_8_clone_1399_1"
  "HH153_clone_nr_9_clone_3100_1"
  "HH153_clone_nr_10_clone_3482_1"
  "HH153_clone_nr_11_clone_234_1"
  "HH153_clone_nr_12_clone_324_1"
  "HH153_clone_nr_13_clone_149_1"
  "HH153_clone_nr_14_clone_3526_1"
  "HH153_clone_nr_15_clone_2512_1"
  "HH153_clone_nr_16_clone_1283_1"
  "HH153_clone_nr_17_clone_215_1"
  "HH153_clone_nr_18_clone_1841_1"
  "HH153_clone_nr_19_clone_2551_1"
  "HH153_clone_nr_20_clone_1221_1"
  "HH153_clone_nr_21_clone_1559_1"
  "HH153_clone_nr_22_clone_2045_1"
  "HH153_clone_nr_23_clone_1227_1"
  "HH153_clone_nr_24_clone_1268_1"
  "HH153_clone_nr_25_clone_1402_1"
  "HH153_clone_nr_26_clone_2601_1"
  "HH153_clone_nr_27_clone_2934_1"
  
)

# Track samples that fail at any step
failed_samples=()

fail_sample() {
  echo "WARNING: $2 failed for $1, skipping." >&2
  failed_samples+=("$1 ($2)")
}

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
  if ! deduplicate $DATA_DIR/${sample}.fasta \
       --root GL \
       --abundance_file abundances.csv \
       --idmapfile idmap.txt > deduplicated.phylip; then
    fail_sample "$sample" "deduplicate"; continue
  fi

  # Parsimony trees
  if ! { mkconfig deduplicated.phylip dnapars > dnapars.cfg && dnapars < dnapars.cfg > dnapars.log; }; then
    fail_sample "$sample" "dnapars"; continue
  fi

  # Make plotting sepecific outdir
  PLOT_DIR_SAMPLE=${PLOT_DIR}/$sample
  mkdir -p $PLOT_DIR_SAMPLE

  # Gctree Ranking
  if ! xvfb-run -a gctree infer outfile abundances.csv --root GL --frame 1 --verbose --outbase $PLOT_DIR_SAMPLE/${sample}; then
    fail_sample "$sample" "gctree infer"; continue
  fi

  echo "Processing of $sample is complete!"

done

# Summary of failed samples
if [ ${#failed_samples[@]} -gt 0 ]; then
  echo ""
  echo "${#failed_samples[@]} sample(s) failed:"
  printf '  %s\n' "${failed_samples[@]}"
  printf '%s\n' "${failed_samples[@]}" > ${OUT_DIR}/failed_samples.txt
else
  echo "All samples processed successfully."
fi
