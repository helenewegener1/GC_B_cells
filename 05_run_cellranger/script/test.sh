#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=24:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/Gina_project/05_run_cellranger/script/test.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/Gina_project/05_run_cellranger/script/test.log
#PBS -N run_cellranger

# Define working directory (which is the out dir here)
WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project/05_run_cellranger/out"
SAMPLE_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"
#CONFIG_DIR="/home/people/helweg/ciir/people/helweg/projects/Gina_project/04_prep_config/out"

module load tools
module load cellranger/9.0.1

# Navigate to your desired output location
cd "$WD"

# List of all your unique sample identifiers (e.g., HH117-SI-PP-...)
# Use read -a (read into an array) to populate samples_array
# read -r -d '' -a samples_array < <(ls "$SAMPLE_DIR" | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort | uniq)

# Run cellranger multi for all IDs
# for ID in "${samples_array[@]}"; do

ID="HH117-SI-PP-nonINF-HLADR-AND-CD19-AND-GC-AND-TFH"
    
echo "--- Starting analysis for sample: ${ID} ---"
# echo "${OUT_DIR}/${ID}_Multi_Analysis"
# echo "${CONFIG_DIR}/multi_config_${ID}.csv"

echo "ID variable content (should be clean):"
echo "$ID" | cat -v

cellranger multi \
    --id="res_${ID}" \
    --csv="../../04_prep_config/out/multi_config_${ID}.csv" \
    --localcores=32 \
    --localmem=128
    
# done
