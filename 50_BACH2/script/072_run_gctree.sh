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
  
  "HH117_clone_nr_1_clone_1881_1"
  "HH117_clone_nr_2_clone_1951_1"
  "HH117_clone_nr_3_clone_2008_1"
  "HH117_clone_nr_4_clone_2667_1"
  "HH117_clone_nr_5_clone_3265_1"
  "HH117_clone_nr_6_clone_4066_1"
  "HH117_clone_nr_7_clone_4261_1"
  "HH117_clone_nr_8_clone_5139_1"
  "HH117_clone_nr_9_clone_5608_1"
  "HH117_clone_nr_10_clone_5834_1"
  "HH117_clone_nr_11_clone_5959_1"
  "HH117_clone_nr_12_clone_6270_1"
  "HH117_clone_nr_13_clone_6752_2"
  "HH117_clone_nr_14_clone_715_1"
  "HH117_clone_nr_15_clone_7509_1"
  "HH117_clone_nr_16_clone_7527_2"
  "HH117_clone_nr_17_clone_9464_4"
  "HH117_clone_nr_18_clone_9817_1"
  
  "HH119_clone_nr_1_clone_10201_1"
  "HH119_clone_nr_2_clone_10275_1"
  "HH119_clone_nr_3_clone_10613_1"
  "HH119_clone_nr_4_clone_11895_2"
  "HH119_clone_nr_5_clone_15039_1"
  "HH119_clone_nr_6_clone_1605_1"
  "HH119_clone_nr_7_clone_17441_1"
  "HH119_clone_nr_8_clone_17441_2"
  "HH119_clone_nr_9_clone_17585_1"
  "HH119_clone_nr_10_clone_17585_2"
  "HH119_clone_nr_11_clone_17731_1"
  "HH119_clone_nr_12_clone_17807_1"
  "HH119_clone_nr_13_clone_1830_1"
  "HH119_clone_nr_14_clone_19551_1"
  "HH119_clone_nr_15_clone_19562_1"
  "HH119_clone_nr_16_clone_21854_1"
  "HH119_clone_nr_17_clone_22252_1"
  "HH119_clone_nr_18_clone_22295_1"
  "HH119_clone_nr_19_clone_22776_1"
  "HH119_clone_nr_20_clone_23207_1"
  "HH119_clone_nr_21_clone_23494_1"
  "HH119_clone_nr_22_clone_24431_1"
  "HH119_clone_nr_23_clone_24483_1"
  "HH119_clone_nr_24_clone_25333_1"
  "HH119_clone_nr_25_clone_26043_1"
  "HH119_clone_nr_26_clone_26121_1"
  "HH119_clone_nr_27_clone_26262_1"
  "HH119_clone_nr_28_clone_26365_1"
  "HH119_clone_nr_29_clone_26418_1"
  "HH119_clone_nr_30_clone_26698_1"
  "HH119_clone_nr_31_clone_27383_1"
  "HH119_clone_nr_32_clone_27583_1"
  "HH119_clone_nr_33_clone_2812_1"
  "HH119_clone_nr_34_clone_28674_1"
  "HH119_clone_nr_35_clone_28704_2"
  "HH119_clone_nr_36_clone_28739_1"
  "HH119_clone_nr_37_clone_29181_1"
  "HH119_clone_nr_38_clone_3804_1"
  "HH119_clone_nr_39_clone_4494_1"
  "HH119_clone_nr_40_clone_5037_1"
  "HH119_clone_nr_41_clone_5652_1"
  "HH119_clone_nr_42_clone_5836_1"
  "HH119_clone_nr_43_clone_6398_2"
  "HH119_clone_nr_44_clone_6457_1"
  "HH119_clone_nr_45_clone_7418_1"
  "HH119_clone_nr_46_clone_8177_1"
  "HH119_clone_nr_47_clone_8457_1"
  "HH119_clone_nr_48_clone_8506_1"
  "HH119_clone_nr_49_clone_8591_1"
  "HH119_clone_nr_50_clone_8598_1"
  "HH119_clone_nr_51_clone_8685_1"
  "HH119_clone_nr_52_clone_8991_1"
  "HH119_clone_nr_53_clone_9038_1"
  "HH119_clone_nr_54_clone_903_1"
  
  "HH151_clone_nr_1_clone_3084_1"
  "HH151_clone_nr_2_clone_4747_1"
  "HH151_clone_nr_3_clone_5221_2"
  "HH151_clone_nr_4_clone_6909_1"
  
  "HH153_clone_nr_1_clone_102_1"
  "HH153_clone_nr_2_clone_1221_1"
  "HH153_clone_nr_3_clone_1224_1"
  "HH153_clone_nr_4_clone_1227_1"
  "HH153_clone_nr_5_clone_1234_1"
  "HH153_clone_nr_6_clone_1237_1"
  "HH153_clone_nr_7_clone_125_1"
  "HH153_clone_nr_8_clone_1268_1"
  "HH153_clone_nr_9_clone_1283_1"
  "HH153_clone_nr_10_clone_1338_1"
  "HH153_clone_nr_11_clone_1399_1"
  "HH153_clone_nr_12_clone_1401_1"
  "HH153_clone_nr_13_clone_1402_1"
  "HH153_clone_nr_14_clone_1414_1"
  "HH153_clone_nr_15_clone_1457_3"
  "HH153_clone_nr_16_clone_149_1"
  "HH153_clone_nr_17_clone_1559_1"
  "HH153_clone_nr_18_clone_1609_1"
  "HH153_clone_nr_19_clone_1680_1"
  "HH153_clone_nr_20_clone_1731_1"
  "HH153_clone_nr_21_clone_1830_1"
  "HH153_clone_nr_22_clone_1841_1"
  "HH153_clone_nr_23_clone_1865_1"
  "HH153_clone_nr_24_clone_188_1"
  "HH153_clone_nr_25_clone_188_2"
  "HH153_clone_nr_26_clone_2045_1"
  "HH153_clone_nr_27_clone_2055_1"
  "HH153_clone_nr_28_clone_215_1"
  "HH153_clone_nr_29_clone_2171_1"
  "HH153_clone_nr_30_clone_2204_1"
  "HH153_clone_nr_31_clone_2209_1"
  "HH153_clone_nr_32_clone_2285_1"
  "HH153_clone_nr_33_clone_2326_1"
  "HH153_clone_nr_34_clone_234_1"
  "HH153_clone_nr_35_clone_2469_1"
  "HH153_clone_nr_36_clone_2512_1"
  "HH153_clone_nr_37_clone_2551_1"
  "HH153_clone_nr_38_clone_2601_1"
  "HH153_clone_nr_39_clone_2708_1"
  "HH153_clone_nr_40_clone_270_1"
  "HH153_clone_nr_41_clone_2717_1"
  "HH153_clone_nr_42_clone_2745_1"
  "HH153_clone_nr_43_clone_2793_1"
  "HH153_clone_nr_44_clone_2907_1"
  "HH153_clone_nr_45_clone_2934_1"
  "HH153_clone_nr_46_clone_2967_1"
  "HH153_clone_nr_47_clone_3058_1"
  "HH153_clone_nr_48_clone_3100_1"
  "HH153_clone_nr_49_clone_3105_1"
  "HH153_clone_nr_50_clone_3107_1"
  "HH153_clone_nr_51_clone_324_1"
  "HH153_clone_nr_52_clone_3344_1"
  "HH153_clone_nr_53_clone_3452_1"
  "HH153_clone_nr_54_clone_3482_1"
  "HH153_clone_nr_55_clone_3526_1"
  "HH153_clone_nr_56_clone_397_1"
  "HH153_clone_nr_57_clone_400_1"
  "HH153_clone_nr_58_clone_483_2"
  "HH153_clone_nr_59_clone_4_1"
  "HH153_clone_nr_60_clone_879_1"
  
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

