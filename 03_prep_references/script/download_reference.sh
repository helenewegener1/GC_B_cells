#!/bin/sh
#PBS -W group_list=dtu_00062 -A dtu_00062 # CIIR
#PBS -l nodes=1:ppn=32:thinnode
#PBS -l gpus=0
#PBS -l mem=32GB
#PBS -l walltime=01:00:00
#PBS -e /home/projects/dtu_00062/people/helweg/projects/Gina_project/03_prep_references/script/download_cellranger_reference.err
#PBS -o /home/projects/dtu_00062/people/helweg/projects/Gina_project/03_prep_references/script/download_cellranger_reference.log
#PBS -N download_cellranger_reference

WD="/home/people/helweg/ciir/people/helweg/projects/Gina_project/03_prep_references/out"

# Navigate to your desired output location
cd "$WD"

# References downloaded from: https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads

# Human reference (GRCh38) - 2024-A
# md5sum: a7b5b7ceefe10e435719edc1a8b8b2fa

wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

# Human V(D)J reference (GRCh38)
# 65b5b033723b07fc1bb5375e5645761c

wget "https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz"
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz

