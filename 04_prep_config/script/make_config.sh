#!/bin/bash
# ==============================================================================
# Sample Structure Analysis Script (Concise Version)
# Generates customized Cell Ranger multi config files based on detected libraries.
# ==============================================================================

# Define working directory
WD="/home/people/helweg/ciir/people/helweg/projects/GC_B_cells/04_prep_config"
CONFIG_DIR="out"

# Each batch is "SAMPLE_DIR|PATTERN1|PATTERN2" - the patterns are an
# extended-regex passed to grep to select which sample IDs to process from
# that directory (in case other, unrelated samples also live there).
BATCHES=(
    "/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF|HH117|HH119"
    "/home/projects/dtu_00062/data/KU10/FASTQ_KU10|HH151|HH153"
)

# --- Setup ---

cd "$WD"
mkdir -p "${CONFIG_DIR}"

# Print header
echo "Sample_ID,GEX_Present,ADT_Present,TCR_Present,BCR_Present,Required_Template"

# --- Analysis Loop ---

for BATCH in "${BATCHES[@]}"; do

  # Split "SAMPLE_DIR|PAT1|PAT2" on '|'
  IFS='|' read -r SAMPLE_DIR PAT1 PAT2 <<< "${BATCH}"
  SAMPLE_PATTERN="${PAT1}|${PAT2}"

  # List of all your unique sample identifiers (e.g., HH117-SI-PP-...) in
  # this batch's directory. Use read -a (read into an array) to populate
  # samples_array.
  read -r -d '' -a samples_array < <(ls "$SAMPLE_DIR" | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort -u | grep -E "${SAMPLE_PATTERN}")

  for ID in "${samples_array[@]}"; do

    # --- 1. Check for library files using single-line conditional assignment ---
    # (Sets status to 'Yes' if files exist, 'No' otherwise)
    GEX_STATUS=$(ls "${SAMPLE_DIR}/GEX-${ID}"*.fastq.gz > /dev/null 2>&1 && echo "Yes" || echo "No")
    ADT_STATUS=$(ls "${SAMPLE_DIR}/ADT-${ID}"*.fastq.gz > /dev/null 2>&1 && echo "Yes" || echo "No")
    TCR_STATUS=$(ls "${SAMPLE_DIR}/TCR-${ID}"*.fastq.gz > /dev/null 2>&1 && echo "Yes" || echo "No")
    BCR_STATUS=$(ls "${SAMPLE_DIR}/BCR-${ID}"*.fastq.gz > /dev/null 2>&1 && echo "Yes" || echo "No")

    TEMPLATE_NAME="NO_MATCHING_TEMPLATE"

    # --- 2. Template Selection Logic (Dynamic Construction) ---

    if [ "${GEX_STATUS}" == "Yes" ] && [ "${BCR_STATUS}" == "Yes" ]; then

      # Start template name with mandatory libraries
      TEMPLATE_BASE="multi_config_GEX_BCR"

      # Append TCR if present
      if [ "${TCR_STATUS}" == "Yes" ]; then TEMPLATE_BASE="${TEMPLATE_BASE}_TCR"; fi

      # Append ADT if present
      if [ "${ADT_STATUS}" == "Yes" ]; then TEMPLATE_BASE="${TEMPLATE_BASE}_ADT"; fi

      # Add suffix
      TEMPLATE_NAME="${TEMPLATE_BASE}_template.csv"

    else
      TEMPLATE_NAME="ERROR: GEX or BCR files missing"

    fi

    # --- 3. CONFIG FILE GENERATION ---

    FEATURE_SUFFIX=""

    if [ "${TEMPLATE_NAME}" != "NO_MATCHING_TEMPLATE" ] && [[ "${TEMPLATE_NAME}" != ERROR* ]]; then

      OUTPUT_CSV="${CONFIG_DIR}/multi_config_${ID}.csv"
      TEMPLATE_PATH="script/${TEMPLATE_NAME}"

      # 3a. Initial substitution: copy template, replace SAMPLE_PREFIX and
      # SAMPLE_DIR (must run first). SAMPLE_DIR contains '/' so use '#' as
      # the sed delimiter instead of the default '/'.
      sed -e "s/SAMPLE_PREFIX/${ID}/g" -e "s#SAMPLE_DIR#${SAMPLE_DIR}#g" "${TEMPLATE_PATH}" > "${OUTPUT_CSV}"

      # 3b. Conditional substitution: Replace FEATURE_REF using 'sed -i' on the new file
      if [ "${ADT_STATUS}" == "Yes" ]; then

        # Check for the most specific case first: HH119 and Pool2
        if [[ "${ID}" == *"HH119"* ]] && [[ "${ID}" == *"Pool2"* ]]; then
          FEATURE_SUFFIX="HH119_pool_2"

        # Check for generic HH119 (implies Pool 1 if Pool 2 was not matched)
        elif [[ "${ID}" == *"HH119"* ]]; then
          FEATURE_SUFFIX="HH119_pool_1"

        # Check for HH117
        elif [[ "${ID}" == *"HH117"* ]]; then
          FEATURE_SUFFIX="HH117"

        # Check for HH151 (OCM + hashtag - one shared feature ref for the run)
        elif [[ "${ID}" == *"HH151"* ]]; then
          FEATURE_SUFFIX="HH151"

        # Check for HH153 (OCM + hashtag - one shared feature ref for the run)
        elif [[ "${ID}" == *"HH153"* ]]; then
          FEATURE_SUFFIX="HH153"
        fi

        # Execute in-place replacement if a reference suffix was determined
        if [ -n "${FEATURE_SUFFIX}" ]; then
          # Use sed -i (in-place) to modify the file created in step 3a
          sed -i "s/FEATURE_SUFFIX/${FEATURE_SUFFIX}/g" "${OUTPUT_CSV}"
        fi
      fi

    # Note: Removed "Generated config file" echo for brevity
    fi

    # Print the CSV line result
    echo "${ID},${GEX_STATUS},${ADT_STATUS},${TCR_STATUS},${BCR_STATUS},${TEMPLATE_NAME},${FEATURE_SUFFIX}"

  done
done

echo "--------------------------------------------------------"
echo "Analysis complete. Customized config files are in the '${CONFIG_DIR}' directory."
