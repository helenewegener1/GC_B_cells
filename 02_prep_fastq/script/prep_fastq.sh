#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=02:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/02_prep_fastq/script/prep_fastq.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/02_prep_fastq/script/prep_fastq.log
#PBS -N prep_fastq

# Copy fastq files into structure expected by cellranger multi
 
# Define working directory
WD="/home/people/helweg/ciir/people/helweg/projects/GC_B_cells"
OUTPUT_BASE_DIR="02_prep_fastq/out"
 
# Each batch is "SAMPLE_DIR|SAMPLE_PATTERN" - pattern is an extended-regex
# passed to grep to select which sample prefixes to pull from that directory.
BATCHES=(
    "/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF|HH117|HH119"
    "/home/projects/dtu_00062/data/KU10/FASTQ_KU10|HH151|HH153"
)
 
# Navigate to your desired output location
cd "$WD"
 
for BATCH in "${BATCHES[@]}"; do
    # Split "SAMPLE_DIR|PAT1|PAT2" on '|'
    IFS='|' read -r SAMPLE_DIR PAT1 PAT2 <<< "${BATCH}"
    SAMPLE_PATTERN="${PAT1}|${PAT2}"
 
    echo "=== Batch: ${SAMPLE_DIR} (${SAMPLE_PATTERN}) ==="
 
    # Use read -a (read into an array) to populate samples_array
    read -r -d '' -a samples_array < <(ls "${SAMPLE_DIR}" | awk -F '_' '{print $1}' | sort -u | grep -E "${SAMPLE_PATTERN}")
 
    for ID in "${samples_array[@]}"; do
        echo "--- Processing sample: ${ID} ---"
 
        # Create the directory (-p ensures parent directories are created if needed)
        mkdir -p "${OUTPUT_BASE_DIR}/${ID}"
 
        # Copy fastq files into structure expected by cellranger multi
        cp "${SAMPLE_DIR}"/*"${ID}"*.fastq.gz "${OUTPUT_BASE_DIR}/${ID}"
    done
done
 
echo "--- COMPLETE ---"


# # Define working directory 
# WD="/home/people/helweg/ciir/people/helweg/projects/GC_B_cells"
# # SAMPLE_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF" # HH117 & HH119 
# SAMPLE_DIR="/home/projects/dtu_00062/data/KU10/FASTQ_KU10" # HH151 & HH153
# OUTPUT_BASE_DIR="02_prep_fastq/out"
# 
# # Navigate to your desired output location
# cd "$WD"
# 
# # Use read -a (read into an array) to populate samples_array
# # read -r -d '' -a samples_array < <(ls "${SAMPLE_DIR}" | awk -F '_' '{print $1}' | sort | uniq) # HH117 & HH119 
# read -r -d '' -a samples_array < <(ls "${SAMPLE_DIR}" | awk -F '_' '{print $1}' | sort -u | grep -E 'HH151|HH153') # HH151 & HH153
# 
# for ID in "${samples_array[@]}"; do
# 
#     echo "--- Processing sample: ${ID} ---"
#     
#     # Create the directory (-p ensures parent directories are created if needed)
#     mkdir -p "${OUTPUT_BASE_DIR}/${ID}"
#     
#     # Copy fastq files into structure expected by cellranger multi
#     cp "${SAMPLE_DIR}"/*"${ID}"*.fastq.gz "${OUTPUT_BASE_DIR}/${ID}"
#     
# done
# 
# echo "--- COMPLETE ---"

