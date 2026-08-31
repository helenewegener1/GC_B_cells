#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=70:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/05_run_cellranger/script/run_cellranger_v9.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/05_run_cellranger/script/run_cellranger_v9.log
#PBS -N run_cellranger_v9

# Define working directory (which is the out dir here)
WD="/home/people/helweg/ciir/people/helweg/projects/GC_B_cells/05_run_cellranger/out_v9"
CONFIG_DIR="/home/people/helweg/ciir/people/helweg/projects/GC_B_cells/04_prep_config/out"

# Each batch is "SAMPLE_DIR|PATIENT1|PATIENT2" - the two patient codes that
# live in that FASTQ directory. All 4 patients across both directories.
BATCHES=(
    "/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF|HH117|HH119"
    "/home/projects/dtu_00062/data/KU10/FASTQ_KU10|HH151|HH153"
)

# --- Run a subset of patients instead of all 4 ---
# Submit with:  qsub -v PATIENT=HH117 run_cellranger_v9.sh          (one patient)
#           or  qsub -v PATIENT=HH117,HH151 run_cellranger_v9.sh    (several, comma-separated, no spaces)
# Leave unset (or "ALL") to run all 4 patients.
PATIENT="${PATIENT:-ALL}"

module load tools
module load cellranger/9.0.1

# Navigate to your desired output location
cd "$WD"

for BATCH in "${BATCHES[@]}"; do

  IFS='|' read -r SAMPLE_DIR PAT1 PAT2 <<< "${BATCH}"

  # Of this batch's two patients, keep only the ones that were requested
  MATCHING=""
  for P in "${PAT1}" "${PAT2}"; do
    if [ "${PATIENT}" == "ALL" ] || [[ ",${PATIENT}," == *",${P},"* ]]; then
      MATCHING="${MATCHING}${MATCHING:+|}${P}"
    fi
  done

  # Neither of this batch's patients was requested - skip the directory entirely
  if [ -z "${MATCHING}" ]; then
    continue
  fi
  SAMPLE_PATTERN="${MATCHING}"

  read -r -d '' -a samples_array < <(ls "$SAMPLE_DIR" | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort -u | grep -E "${SAMPLE_PATTERN}")

  # Run cellranger multi for each matching sample ID
  for ID in "${samples_array[@]}"; do

    echo "--- Starting analysis for sample: ${ID} ---"

    cellranger multi \
        --id="res_${ID}" \
        --csv="${CONFIG_DIR}/multi_config_${ID}.csv" \
        --localcores=32 \
        --localmem=128

  done
done

