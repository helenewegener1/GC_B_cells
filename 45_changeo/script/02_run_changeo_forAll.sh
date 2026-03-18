#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=64GB
#PBS -l walltime=24:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_changeo/script/run_changeo.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_changeo/script/run_changeo.log
#PBS -N run_changeo

module load tools singularity/4.3.0

# Change path 
cd /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_changeo/script

# Get change-O from singularity
./01_get_changeo_from_singularity.sh

# Mounting the data I need 
export SINGULARITY_BIND="/home/projects/dtu_00062:/home/projects/dtu_00062,/scratch/helweg:/scratch/helweg"

SAMPLE_NAMES=$(ls /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/05_run_cellranger/out_v9/ | sed 's/^res_//')

for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
    ./02_run_changeo.sh $SAMPLE_NAME
done

