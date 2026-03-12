#!/bin/bash

OUT_DIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/mixcr_output"
FIGS_DIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/mixcr_figures"

sample="BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"

mkdir -p ${FIGS_DIR}

# mixcr exportQc align \
#       ${OUT_DIR}/*.clns \
#       ${FIGS_DIR}/alignQc.pdf
 
# Allele inference     
mixcr findAlleles \
      -f \
      --report ${OUT_DIR}/${sample}.findAlleles.report.txt \
      --export-alleles-mutations ${OUT_DIR}/${sample}_alleles.tsv \
      --export-library ${OUT_DIR}/${sample}_alleles.json \
      --output-template {file_dir_path}/{file_name}.reassigned.clns \
      ${OUT_DIR}/*.clns
      
# Export clones
# mixcr exportClones \ 
#       ${OUT_DIR}/DCh_T1_pbmc_H_2.reassigned.clns \
#       ${OUT_DIR}/DCh_T1_pbmc_H_2.txt
