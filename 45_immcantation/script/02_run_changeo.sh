#!/bin/bash

# Mounting the data I need 
# export SINGULARITY_BIND="/home/projects/dtu_00062:/home/projects/dtu_00062,/scratch/helweg:/scratch/helweg"

# Files
IMAGE="/scratch/helweg/singularity/immcantation_suite-4.7.0.sif"

# SAMPLE_NAME="HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"
SAMPLE_NAME=$1   # pass sample name as argument

CELLRANGER_OUT_PATH="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/05_run_cellranger/out_v9/res_${SAMPLE_NAME}/outs/per_sample_outs/res_${SAMPLE_NAME}/vdj_b"
FILTERED_CONTIGS_PATH="${CELLRANGER_OUT_PATH}/filtered_contig.fasta"
FILTERED_CONTIGS_ANNOTATIONS_PATH="${CELLRANGER_OUT_PATH}/filtered_contig_annotations.csv"

# Output directory
OUTDIR="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_immcantation/out/${SAMPLE_NAME}"
mkdir -p $OUTDIR

# Step 1: V(D)J assignment with IgBLAST
singularity exec $IMAGE AssignGenes.py igblast \
    -s $FILTERED_CONTIGS_PATH \
    -b /usr/local/share/igblast \
    --organism human --loci ig --format blast \
    --outdir $OUTDIR --outname $SAMPLE_NAME --nproc 4

# Step 2: Parse into AIRR format, merging Cell Ranger barcode annotations
singularity exec $IMAGE MakeDb.py igblast \
    -i ${OUTDIR}/${SAMPLE_NAME}_igblast.fmt7 \
    -s $FILTERED_CONTIGS_PATH \
    -r /usr/local/share/germlines/imgt/human/vdj/ \
    --10x $FILTERED_CONTIGS_ANNOTATIONS_PATH \
    --extended --outdir $OUTDIR --outname $SAMPLE_NAME

# Step 3: Filter to productive sequences
singularity exec $IMAGE ParseDb.py select \
    -d ${OUTDIR}/${SAMPLE_NAME}_db-pass.tsv \
    -f productive -u T \
    --outdir $OUTDIR --outname $SAMPLE_NAME

# Step 4: Split heavy and light chains
singularity exec $IMAGE ParseDb.py select \
    -d ${OUTDIR}/${SAMPLE_NAME}_parse-select.tsv \
    -f locus -u IGH --regex \
    --outname ${SAMPLE_NAME}_heavy --outdir $OUTDIR

singularity exec $IMAGE ParseDb.py select \
    -d ${OUTDIR}/${SAMPLE_NAME}_parse-select.tsv \
    -f locus -u "IG[LK]" --regex \
    --outname ${SAMPLE_NAME}_light --outdir $OUTDIR

# Step 5: Reconstruct germlines
# singularity exec $IMAGE CreateGermlines.py \
#     -d ${OUTDIR}/${SAMPLE_NAME}_heavy_parse-select.tsv \
#     -r /usr/local/share/germlines/imgt/human/vdj/ \
#     --cloned --outdir $OUTDIR --outname ${SAMPLE_NAME}_heavy
    
singularity exec $IMAGE CreateGermlines.py \
    -d ${OUTDIR}/${SAMPLE_NAME}_heavy_parse-select.tsv \
    -r /usr/local/share/germlines/imgt/human/vdj/ \
    --outdir $OUTDIR --outname ${SAMPLE_NAME}_heavy
    
  