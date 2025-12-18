#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=02:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/Gina_project/02_prep_fastq/script/prep_fastq.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/Gina_project/02_prep_fastq/script/prep_fastq.log
#PBS -N prep_fastq

# Copy fastq files into structure expected by cellranger multi

# Define working directory 
WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project"
SAMPLE_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"
OUTPUT_BASE_DIR="02_prep_fastq/out"

# Navigate to your desired output location
cd "$WD"

# Use read -a (read into an array) to populate samples_array
read -r -d '' -a samples_array < <(ls "${SAMPLE_DIR}" | awk -F '_' '{print $1}' | sort | uniq)

for ID in "${samples_array[@]}"; do

    echo "--- Processing sample: ${ID} ---"
    
    # Create the directory (-p ensures parent directories are created if needed)
    mkdir -p "${OUTPUT_BASE_DIR}/${ID}"
    
    # Copy fastq files into structure expected by cellranger multi
    cp "${SAMPLE_DIR}"/*"${ID}"*.fastq.gz "${OUTPUT_BASE_DIR}/${ID}"
    
done

echo "--- COMPLETE ---"

