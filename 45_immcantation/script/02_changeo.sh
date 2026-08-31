#!/bin/bash

# Mounting the data I need 
# export SINGULARITY_BIND="/home/projects/dtu_00062:/home/projects/dtu_00062,/scratch/helweg:/scratch/helweg"

# Files
IMAGE="/scratch/helweg/singularity/immcantation_suite-4.7.0.sif"

# SAMPLE_NAME="HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2"
SAMPLE_NAME=$1   # pass sample name as argument

CELLRANGER_BASE="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/05_run_cellranger/out_v9/res_${SAMPLE_NAME}/outs/per_sample_outs"
IMMCANTATION_OUT_BASE="/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/45_immcantation/out"

# --- Detect OCM pools -------------------------------------------------------
# Non-OCM samples have a single per_sample_outs subdir named "res_${SAMPLE_NAME}".
# OCM samples have one subdir per OCM pool (whatever sample_id values were used
# in the [samples] section of the multi config, e.g. Blue/Red/Green/Yellow), and
# each pool needs its own Change-O run - same detection logic as
# 06_seurat_load/script/load_seurat_v9.R.
mapfile -t POOL_DIRS < <(ls -1 "${CELLRANGER_BASE}")

OCM_LABELS=()
if [ "${#POOL_DIRS[@]}" -eq 1 ] && [ "${POOL_DIRS[0]}" == "res_${SAMPLE_NAME}" ]; then
    # Not OCM - single pool, keep the plain sample name
    OCM_LABELS=("${SAMPLE_NAME}")
else
    # OCM - one pool per detected subdir, labeled SAMPLE_pool
    for OCM in "${POOL_DIRS[@]}"; do
        OCM_LABELS+=("${SAMPLE_NAME}_${OCM}")
    done
fi

# --- Run Change-O for each OCM ---------------------------------------------
for i in "${!POOL_DIRS[@]}"; do

    POOL_DIR="${POOL_DIRS[$i]}"
    OUTNAME="${OCM_LABELS[$i]}"

    CELLRANGER_OUT_PATH="${CELLRANGER_BASE}/${POOL_DIR}/vdj_b"
    FILTERED_CONTIGS_PATH="${CELLRANGER_OUT_PATH}/filtered_contig.fasta"
    FILTERED_CONTIGS_ANNOTATIONS_PATH="${CELLRANGER_OUT_PATH}/filtered_contig_annotations.csv"

    OUTDIR="${IMMCANTATION_OUT_BASE}/${OUTNAME}"
    mkdir -p "$OUTDIR"

    echo "=== Processing ${OUTNAME} (${CELLRANGER_OUT_PATH}) ==="

    # Step 1: V(D)J assignment with IgBLAST
    singularity exec $IMAGE AssignGenes.py igblast \
        -s $FILTERED_CONTIGS_PATH \
        -b /usr/local/share/igblast \
        --organism human --loci ig --format blast \
        --outdir $OUTDIR --outname $OUTNAME --nproc 4

    # Step 2: Parse into AIRR format, merging Cell Ranger barcode annotations
    singularity exec $IMAGE MakeDb.py igblast \
        -i ${OUTDIR}/${OUTNAME}_igblast.fmt7 \
        -s $FILTERED_CONTIGS_PATH \
        -r /usr/local/share/germlines/imgt/human/vdj/ \
        --10x $FILTERED_CONTIGS_ANNOTATIONS_PATH \
        --extended --outdir $OUTDIR --outname $OUTNAME

    # Step 3: Filter to productive sequences
    singularity exec $IMAGE ParseDb.py select \
        -d ${OUTDIR}/${OUTNAME}_db-pass.tsv \
        -f productive -u T \
        --outdir $OUTDIR --outname $OUTNAME

    # Step 4: Split heavy and light chains
    singularity exec $IMAGE ParseDb.py select \
        -d ${OUTDIR}/${OUTNAME}_parse-select.tsv \
        -f locus -u IGH --regex \
        --outname ${OUTNAME}_heavy --outdir $OUTDIR

    singularity exec $IMAGE ParseDb.py select \
        -d ${OUTDIR}/${OUTNAME}_parse-select.tsv \
        -f locus -u "IG[LK]" --regex \
        --outname ${OUTNAME}_light --outdir $OUTDIR

    # Step 5: Reconstruct germlines
    singularity exec $IMAGE CreateGermlines.py \
        -d ${OUTDIR}/${OUTNAME}_heavy_parse-select.tsv \
        -r /usr/local/share/germlines/imgt/human/vdj/ \
        --outdir $OUTDIR --outname ${OUTNAME}_heavy

    singularity exec $IMAGE CreateGermlines.py \
        -d ${OUTDIR}/${OUTNAME}_light_parse-select.tsv \
        -r /usr/local/share/germlines/imgt/human/vdj/ \
        --outdir $OUTDIR --outname ${OUTNAME}_light

done
