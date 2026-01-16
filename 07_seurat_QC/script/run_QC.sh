#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=24:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/07_seurat_QC/script/run_QC.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/07_seurat_QC/script/run_QC.log
#PBS -N run_QC

module load tools
module load gcc/14.2.0
module load intel/basekit/INITIALIZE/2023.0.0
module load intel/basekit/mkl/2023.0.0
module load R/4.5.0

Rscript /home/people/helweg/ciir/people/helweg/projects/GC_B_cells/07_seurat_QC/script/QC_v8.R
# Rscript /home/people/helweg/ciir/people/helweg/projects/GC_B_cells/07_seurat_QC/script/QC_v9.R
