## Data 
Human samples. 

HH117 is a patient with Crohns Disease where tissue have been taken from a inflamed (INF) and non-inflamed (nonINF) site.

HH119 is a patient with Colorectal Cancer where the sample have been taken from a non-tumor cite and is therefore used as the best healthy control we can get.

The samples from SI-PP are separated into follicles that each are labeled with a hashtag.

### Patient 
**HH117**: patient with Crohns Disease (CD)

**HH119**: patient with Colorectal Cancer (CRC) ~ Control 

### Inflammation Status
**nonINF**: non-inflamed

**INF**: inflamed 

### Library Types (Sequencing Data)
**GEX**: Gene expression 

**ADT**: Antibody Capture

**BCR**: B cell receptor 

**TCR**: T cell receptor 

### Sites in the Intestine
**SI-PP**: Small Intestine - Peyer's Patch

**SI-MILF**: Small Intestine - Muscularis layer of the Intestinal Lamina Fibrinogen-rich area

**SILP**: Small Intestine - Lamina Propria

**COLP**: Colon - Lamina Propria

**CO-SMILF**: Colon - Submucosal layer or Muscularis layer of the Intestinal Lamina Fibrinogen-rich area

### Cell Types
**CD19**: Cluster of Differentiation 19 (Marker for B-cells)

**GC**: Germinal Center (B-cells)

**TFH**: T Follicular Helper (T-cells)

**PB**: Plasmablast

**HLADR**: Human Leukocyte Antigen – DR isotype (MHC Class II molecule)

**PC**: Plasma Cell

## Structure of analysis 
fastq files on computerome: /home/projects/dtu_00062/data/KU09/FASTQ_ku09_mkfastq/outs/fastq_path/HKL3YDSXF

**00_data**: Excel files with information about hashtag antibodies from each 

**01_prep_feature_ref**: Using the excel files in 00_data to make feature_ref files for each hashtag-combination. 

**02_prep_fastq**: fastq files are copied from the original directory into a structure where all fastq files from the same sample are in one folder. This might not be needed anyways...

**03_prep_references**: Download human genome reference and human vdj reference for B and T cell receptors.

**04_prep_config**: From self-made template files, there is created a config file per sample depending on which files exist for the given sample. 

**05_run_cellranger**: Running cellranger multi using the config files. 




