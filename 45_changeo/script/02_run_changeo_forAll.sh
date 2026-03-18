#!/bin/bash
mkdir -p /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_changeo/logs

SAMPLE_NAMES=$(ls /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/05_run_cellranger/out_v9/ | sed 's/^res_//')

for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
    ./02_run_changeo.sh $SAMPLE_NAME
done

