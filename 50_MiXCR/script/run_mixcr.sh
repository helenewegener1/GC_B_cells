# Define paths 
FASTQ_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"
OUT_DIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/mixcr_output"
SAMPLESHEET="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/samplesheet.csv"
WHITELIST="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/00_data/3M-february-2018.txt"
# SAMPLESHEET="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/samplesheet_test.csv"
TMP_DIR="/scratch/helweg/mixcr_tmp"

mkdir -p ${OUT_DIR}
mkdir -p ${TMP_DIR}

while IFS=',' read -r sampleid prefix; do
  echo "Submitting: ${sampleid}"

  cat ${FASTQ_DIR}/${prefix}_L00{1,2,3,4}_R1_001.fastq.gz > ${TMP_DIR}/${sampleid}_R1.fastq.gz
  cat ${FASTQ_DIR}/${prefix}_L00{1,2,3,4}_R2_001.fastq.gz > ${TMP_DIR}/${sampleid}_R2.fastq.gz

  # mixcr analyze 10x-sc-xcr-vdj \
  mixcr analyze 10x-sc-xcr-vdj-gemx-v3 \
    --species hsa \
    # --tag-pattern "^(CELL:N{16})(UMI:N{12})\^(R2:*)" \
    # --set-whitelist CELL=file:${WHITELIST} \
    ${TMP_DIR}/${sampleid}_R1.fastq.gz \
    ${TMP_DIR}/${sampleid}_R2.fastq.gz \
    ${OUT_DIR}/${sampleid} \
    --threads 16 &   # <-- run in background

  # Wait every 2 samples before starting next pair
  if (( $(jobs -r | wc -l) >= 2 )); then
    wait
  fi

done < <(tail -n +2 ${SAMPLESHEET})
wait  # wait for last batch

echo "All jobs submitted!"
