#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=64GB
#PBS -l walltime=12:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_immcantation/script/run_add_light_chain.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_immcantation/script/run_add_light_chain.log
#PBS -N run_add_light_chain

module load tools singularity/4.3.0

cd /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_immcantation/script

# Get change-O from singularity
./01_get_changeo_from_singularity.sh

# Mounting the data I need 
export SINGULARITY_BIND="/home/projects/dtu_00062:/home/projects/dtu_00062,/scratch/helweg:/scratch/helweg"

IMAGE="/scratch/helweg/singularity/immcantation_suite-4.7.0.sif"

while IFS= read -r SAMPLE_NAME; do

  # Output directory
  OUTDIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_immcantation/out/${SAMPLE_NAME}"
  
  # Light chain correction (run AFTER SCOPer in R has added clone_id)
  # This step expects the heavy chain file to already have clone_id column
  # Run this after exporting SCOPer results back to TSV from R
  singularity exec $IMAGE light_cluster.py \
      -d ${OUTDIR}/${SAMPLE_NAME}_heavy_germ-pass_clone-pass.tsv \
      -e ${OUTDIR}/${SAMPLE_NAME}_light_parse-select.tsv \
      -o ${OUTDIR}/${SAMPLE_NAME}_10X_clone-pass.tsv \
      --doublets count
    
done < <(ls /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_immcantation/out/ | sed 's/^res_//' | grep '^HH')

