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

FASTQ_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"
OUT_DIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/mixcr_output"
# SAMPLESHEET="samplesheet.csv"
SAMPLESHEET="samplesheet_test.csv"

mkdir -p ${OUT_DIR}

while IFS=',' read -r sampleid prefix; do
echo "Submitting: ${sampleid}"

sbatch \
--job-name=mixcr_${sampleid} \
--cpus-per-task=16 \
--mem=64G \
--time=12:00:00 \
--output=${OUT_DIR}/${sampleid}.log \
--wrap="
            module load tools jdk/24.0.2 mixcr/4.7.0 &&
            mixcr analyze 10x-vdj-bcr \
                ${FASTQ_DIR}/${prefix}_L001_R1_001.fastq.gz \
                ${FASTQ_DIR}/${prefix}_L001_R2_001.fastq.gz \
                ${FASTQ_DIR}/${prefix}_L002_R1_001.fastq.gz \
                ${FASTQ_DIR}/${prefix}_L002_R2_001.fastq.gz \
                ${FASTQ_DIR}/${prefix}_L003_R1_001.fastq.gz \
                ${FASTQ_DIR}/${prefix}_L003_R2_001.fastq.gz \
                ${FASTQ_DIR}/${prefix}_L004_R1_001.fastq.gz \
                ${FASTQ_DIR}/${prefix}_L004_R2_001.fastq.gz \
                ${OUT_DIR}/${sampleid} \
                --threads 16
        "
done < <(tail -n +2 ${SAMPLESHEET})

echo "All jobs submitted!"
