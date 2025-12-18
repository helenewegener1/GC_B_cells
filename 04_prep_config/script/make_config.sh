#!/bin/bash
# ==============================================================================
# Sample Structure Analysis Script (Concise Version)
# Generates customized Cell Ranger multi config files based on detected libraries.
# ==============================================================================

# Define working directory 
WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project/04_prep_config"
SAMPLE_DIR="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"
CONFIG_DIR="out" 

# --- Setup ---

cd "$WD"
mkdir -p "${CONFIG_DIR}"

# Print header
echo "Sample_ID,GEX_Present,ADT_Present,TCR_Present,BCR_Present,Required_Template"

# --- Analysis Loop ---

# List of all your unique sample identifiers (e.g., HH117-SI-PP-...)
# Use read -a (read into an array) to populate samples_array
read -r -d '' -a samples_array < <(ls "$SAMPLE_DIR" | awk -F '_' '{print $1}' | awk -F '-' '{sub($1 FS, "", $0); print $0}' | sort | uniq)

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
  
  if [ "${TEMPLATE_NAME}" != "NO_MATCHING_TEMPLATE" ] && [[ "${TEMPLATE_NAME}" != ERROR* ]]; then
  
    OUTPUT_CSV="${CONFIG_DIR}/multi_config_${ID}.csv"
    TEMPLATE_PATH="script/${TEMPLATE_NAME}"
    
    # 3a. Initial substitution: Copy template and replace SAMPLE_PREFIX (must run first)
    sed "s/SAMPLE_PREFIX/${ID}/g" "${TEMPLATE_PATH}" > "${OUTPUT_CSV}"

    # 3b. Conditional substitution: Replace FEATURE_REF using 'sed -i' on the new file
    FEATURE_SUFFIX=""
    
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
        fi
        
        # Execute in-place replacement if a reference path was determined
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

echo "--------------------------------------------------------"
echo "Analysis complete. Customized config files are in the '${CONFIG_DIR}' directory."
