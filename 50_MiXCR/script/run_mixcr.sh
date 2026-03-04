#!/bin/bash
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=128GB
#PBS -l walltime=24:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/run_mixcr.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/GC_B_cells/50_MiXCR/script/run_mixcr.log
#PBS -N run_mixcr

module load tools jdk/24.0.2 openjdk/23 java/17-openjdk jre/17-openjdk mixcr/4.7.0

proj_dir="/home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF"

mixcr analyze 10x-vdj-bcr \
    $proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L001_R1_001.fastq.gz,$proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L002_R1_001.fastq.gz,$proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L003_R1_001.fastq.gz,$proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L004_R1_001.fastq.gz \
    $proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L001_R2_001.fastq.gz,$proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L002_R2_001.fastq.gz,$proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L003_R2_001.fastq.gz,$proj_dir/BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2_S28_L004_R2_001.fastq.gz \
    BCR-HH119-SI-PP-GC-AND-PB-AND-TFH-Pool2
    
    