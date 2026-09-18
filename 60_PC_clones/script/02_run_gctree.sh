#!/bin/bash
set -e

# Define directories 
WD=/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/60_PC_clones

DATA_DIR=${WD}/fasta
OUT_DIR=${WD}/out
PLOT_DIR=${WD}/plot

# sample_list=$(ls $DATA_DIR | cut -d "." -f1)
sample_list=(
  
  HH117-SILP-INF_clone_nr_10_clone_169_1
  HH117-SILP-INF_clone_nr_11_clone_4805_1
  HH117-SILP-INF_clone_nr_12_clone_488_1
  HH117-SILP-INF_clone_nr_13_clone_2477_1
  HH117-SILP-INF_clone_nr_14_clone_9145_1
  HH117-SILP-INF_clone_nr_15_clone_6575_1
  HH117-SILP-INF_clone_nr_16_clone_9593_1
  HH117-SILP-INF_clone_nr_17_clone_4511_1
  HH117-SILP-INF_clone_nr_18_clone_506_1
  HH117-SILP-INF_clone_nr_19_clone_3791_1
  HH117-SILP-INF_clone_nr_1_clone_8865_1
  HH117-SILP-INF_clone_nr_20_clone_8248_1
  HH117-SILP-INF_clone_nr_2_clone_9227_1
  HH117-SILP-INF_clone_nr_3_clone_6771_1
  HH117-SILP-INF_clone_nr_4_clone_2044_1
  HH117-SILP-INF_clone_nr_5_clone_5857_1
  HH117-SILP-INF_clone_nr_6_clone_4217_1
  HH117-SILP-INF_clone_nr_7_clone_6913_1
  HH117-SILP-INF_clone_nr_8_clone_7511_1
  HH117-SILP-INF_clone_nr_9_clone_9813_1
  HH117-SILP-nonINF_clone_nr_10_clone_8865_1
  HH117-SILP-nonINF_clone_nr_11_clone_2044_1
  HH117-SILP-nonINF_clone_nr_12_clone_9012_1
  HH117-SILP-nonINF_clone_nr_13_clone_10134_1
  HH117-SILP-nonINF_clone_nr_14_clone_8183_1
  HH117-SILP-nonINF_clone_nr_15_clone_8760_1
  HH117-SILP-nonINF_clone_nr_16_clone_3852_1
  HH117-SILP-nonINF_clone_nr_17_clone_8248_1
  HH117-SILP-nonINF_clone_nr_18_clone_5261_1
  HH117-SILP-nonINF_clone_nr_19_clone_6334_1
  HH117-SILP-nonINF_clone_nr_1_clone_9227_1
  HH117-SILP-nonINF_clone_nr_20_clone_7656_1
  HH117-SILP-nonINF_clone_nr_2_clone_6771_1
  HH117-SILP-nonINF_clone_nr_3_clone_8771_1
  HH117-SILP-nonINF_clone_nr_4_clone_9502_1
  HH117-SILP-nonINF_clone_nr_5_clone_8588_1
  HH117-SILP-nonINF_clone_nr_6_clone_4217_1
  HH117-SILP-nonINF_clone_nr_7_clone_5857_1
  HH117-SILP-nonINF_clone_nr_8_clone_9247_1
  HH117-SILP-nonINF_clone_nr_9_clone_6956_1
  HH119-COLP_clone_nr_10_clone_23464_1
  HH119-COLP_clone_nr_11_clone_28568_1
  HH119-COLP_clone_nr_12_clone_28852_1
  HH119-COLP_clone_nr_13_clone_11639_1
  HH119-COLP_clone_nr_14_clone_25750_1
  HH119-COLP_clone_nr_15_clone_26143_1
  HH119-COLP_clone_nr_16_clone_28017_1
  HH119-COLP_clone_nr_17_clone_3650_1
  HH119-COLP_clone_nr_18_clone_1335_1
  HH119-COLP_clone_nr_19_clone_2149_1
  HH119-COLP_clone_nr_1_clone_24582_1
  HH119-COLP_clone_nr_20_clone_25699_1
  HH119-COLP_clone_nr_2_clone_4718_1
  HH119-COLP_clone_nr_3_clone_4576_1
  HH119-COLP_clone_nr_4_clone_14_1
  HH119-COLP_clone_nr_5_clone_2472_1
  HH119-COLP_clone_nr_6_clone_17505_1
  HH119-COLP_clone_nr_7_clone_15065_1
  HH119-COLP_clone_nr_8_clone_24963_1
  HH119-COLP_clone_nr_9_clone_24974_1
  HH119-SILP_clone_nr_10_clone_11664_1
  HH119-SILP_clone_nr_11_clone_24338_1
  HH119-SILP_clone_nr_12_clone_24405_1
  HH119-SILP_clone_nr_13_clone_12918_1
  HH119-SILP_clone_nr_14_clone_17590_1
  HH119-SILP_clone_nr_15_clone_21586_1
  HH119-SILP_clone_nr_16_clone_22404_1
  HH119-SILP_clone_nr_17_clone_5624_1
  HH119-SILP_clone_nr_18_clone_12704_1
  HH119-SILP_clone_nr_19_clone_17998_1
  HH119-SILP_clone_nr_1_clone_28238_1
  HH119-SILP_clone_nr_20_clone_11580_1
  HH119-SILP_clone_nr_2_clone_15994_1
  HH119-SILP_clone_nr_3_clone_4157_1
  HH119-SILP_clone_nr_4_clone_14647_1
  HH119-SILP_clone_nr_5_clone_24427_1
  HH119-SILP_clone_nr_6_clone_28541_1
  HH119-SILP_clone_nr_7_clone_346_1
  HH119-SILP_clone_nr_8_clone_24405_2
  HH119-SILP_clone_nr_9_clone_19224_1
  HH151-SILP-INF_clone_nr_10_clone_532_1
  HH151-SILP-INF_clone_nr_11_clone_8511_1
  HH151-SILP-INF_clone_nr_12_clone_366_1
  HH151-SILP-INF_clone_nr_13_clone_4131_1
  HH151-SILP-INF_clone_nr_14_clone_2935_1
  HH151-SILP-INF_clone_nr_15_clone_2599_1
  HH151-SILP-INF_clone_nr_16_clone_4747_1
  HH151-SILP-INF_clone_nr_17_clone_8000_1
  HH151-SILP-INF_clone_nr_18_clone_920_1
  HH151-SILP-INF_clone_nr_19_clone_1292_1
  HH151-SILP-INF_clone_nr_1_clone_8016_1
  HH151-SILP-INF_clone_nr_20_clone_656_1
  HH151-SILP-INF_clone_nr_2_clone_2697_1
  HH151-SILP-INF_clone_nr_3_clone_6425_1
  HH151-SILP-INF_clone_nr_4_clone_6197_1
  HH151-SILP-INF_clone_nr_5_clone_6261_1
  HH151-SILP-INF_clone_nr_6_clone_1826_1
  HH151-SILP-INF_clone_nr_7_clone_8753_1
  HH151-SILP-INF_clone_nr_8_clone_2691_1
  HH151-SILP-INF_clone_nr_9_clone_440_1
  HH151-SILP-nonINF_clone_nr_10_clone_2114_1
  HH151-SILP-nonINF_clone_nr_11_clone_3650_1
  HH151-SILP-nonINF_clone_nr_12_clone_1292_1
  HH151-SILP-nonINF_clone_nr_13_clone_3868_1
  HH151-SILP-nonINF_clone_nr_14_clone_4040_1
  HH151-SILP-nonINF_clone_nr_15_clone_4960_1
  HH151-SILP-nonINF_clone_nr_16_clone_6197_1
  HH151-SILP-nonINF_clone_nr_17_clone_8277_1
  HH151-SILP-nonINF_clone_nr_18_clone_4131_1
  HH151-SILP-nonINF_clone_nr_19_clone_5377_1
  HH151-SILP-nonINF_clone_nr_1_clone_8124_1
  HH151-SILP-nonINF_clone_nr_20_clone_5675_1
  HH151-SILP-nonINF_clone_nr_2_clone_2977_1
  HH151-SILP-nonINF_clone_nr_3_clone_6091_1
  HH151-SILP-nonINF_clone_nr_4_clone_8016_1
  HH151-SILP-nonINF_clone_nr_5_clone_366_1
  HH151-SILP-nonINF_clone_nr_6_clone_6817_1
  HH151-SILP-nonINF_clone_nr_7_clone_6059_1
  HH151-SILP-nonINF_clone_nr_8_clone_2017_1
  HH151-SILP-nonINF_clone_nr_9_clone_8044_1
  HH153-SILP-INF_clone_nr_10_clone_2908_1
  HH153-SILP-INF_clone_nr_11_clone_2995_1
  HH153-SILP-INF_clone_nr_12_clone_3303_1
  HH153-SILP-INF_clone_nr_13_clone_1526_1
  HH153-SILP-INF_clone_nr_14_clone_1772_1
  HH153-SILP-INF_clone_nr_15_clone_1894_1
  HH153-SILP-INF_clone_nr_16_clone_1958_1
  HH153-SILP-INF_clone_nr_17_clone_1985_1
  HH153-SILP-INF_clone_nr_18_clone_2010_1
  HH153-SILP-INF_clone_nr_19_clone_3037_1
  HH153-SILP-INF_clone_nr_1_clone_1515_1
  HH153-SILP-INF_clone_nr_20_clone_3105_1
  HH153-SILP-INF_clone_nr_2_clone_1655_1
  HH153-SILP-INF_clone_nr_3_clone_2204_1
  HH153-SILP-INF_clone_nr_4_clone_2080_1
  HH153-SILP-INF_clone_nr_5_clone_717_1
  HH153-SILP-INF_clone_nr_6_clone_568_1
  HH153-SILP-INF_clone_nr_7_clone_1387_1
  HH153-SILP-INF_clone_nr_8_clone_2564_1
  HH153-SILP-INF_clone_nr_9_clone_2693_1
  HH153-SILP-nonINF_clone_nr_10_clone_2519_1
  HH153-SILP-nonINF_clone_nr_11_clone_1113_1
  HH153-SILP-nonINF_clone_nr_12_clone_1452_1
  HH153-SILP-nonINF_clone_nr_13_clone_2910_1
  HH153-SILP-nonINF_clone_nr_14_clone_3105_1
  HH153-SILP-nonINF_clone_nr_15_clone_423_1
  HH153-SILP-nonINF_clone_nr_16_clone_1020_1
  HH153-SILP-nonINF_clone_nr_17_clone_1396_1
  HH153-SILP-nonINF_clone_nr_18_clone_1598_1
  HH153-SILP-nonINF_clone_nr_19_clone_1666_1
  HH153-SILP-nonINF_clone_nr_1_clone_1655_1
  HH153-SILP-nonINF_clone_nr_20_clone_135_1
  HH153-SILP-nonINF_clone_nr_2_clone_882_1
  HH153-SILP-nonINF_clone_nr_3_clone_1958_1
  HH153-SILP-nonINF_clone_nr_4_clone_1772_1
  HH153-SILP-nonINF_clone_nr_5_clone_818_1
  HH153-SILP-nonINF_clone_nr_6_clone_1410_1
  HH153-SILP-nonINF_clone_nr_7_clone_1822_1
  HH153-SILP-nonINF_clone_nr_8_clone_1526_1
  HH153-SILP-nonINF_clone_nr_9_clone_2209_1

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

