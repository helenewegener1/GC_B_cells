#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=24:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/run_mixcr.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/run_mixcr.log
#PBS -N run_mixcr

module load tools jdk/24.0.2 openjdk/23 java/17-openjdk jre/17-openjdk mixcr/4.7.0

# mixcr activate-license 
# E-ZBKUJLUYAGTUQQLPMXPDVOHYKMVBWUFNATIMAYZBPLPQSLIU

FASTQ_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"
OUT_DIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/mixcr_output"
SAMPLESHEET="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/samplesheet.csv"
# SAMPLESHEET="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/samplesheet_test.csv"
TMP_DIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/tmp"

mkdir -p ${OUT_DIR}
mkdir -p ${TMP_DIR}

while IFS=',' read -r sampleid prefix; do
echo "Submitting: ${sampleid}"

  # Concatenate lanes - MiXCR needs this 
  cat ${FASTQ_DIR}/${prefix}_L00{1,2,3,4}_R1_001.fastq.gz > ${TMP_DIR}/${sampleid}_R1.fastq.gz
  cat ${FASTQ_DIR}/${prefix}_L00{1,2,3,4}_R2_001.fastq.gz > ${TMP_DIR}/${sampleid}_R2.fastq.gz

  mixcr analyze 10x-sc-xcr-vdj \
    --species hsa \
    ${TMP_DIR}/${sampleid}_R1.fastq.gz \
    ${TMP_DIR}/${sampleid}_R2.fastq.gz \
    ${OUT_DIR}/${sampleid} \
    --threads 16
                
done < <(tail -n +2 ${SAMPLESHEET})

echo "All jobs submitted!"
