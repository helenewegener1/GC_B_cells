#!/usr/bin/env python3

# Create environment: python3.11 -m venv gctree_env_local
# source gctree_env_local/bin/activate
# gctree_env_local/bin/pip install gctree

# Following: https://matsen.group/gctree/rendering-demo.html

print("Importing packages...")
import gctree
import pickle
import numpy as np
import os
import pandas as pd
import re

# For computerome run
os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["XDG_RUNTIME_DIR"] = "/tmp/runtime-runner"
os.environ["MPLBACKEND"] = "agg"
os.environ["MPLCONFIGDIR"] = "/home/projects/dtu_00062/people/helweg/projects/GC_B_cells/.matplotlib"

samples_dict = {
  
  "HH117-SILP-INF_clone_nr_10_clone_169_1": "HH117 SI-LP-INF: 10. largest PC clone",
  "HH117-SILP-INF_clone_nr_11_clone_4805_1": "HH117 SI-LP-INF: 11. largest PC clone",
  "HH117-SILP-INF_clone_nr_12_clone_488_1": "HH117 SI-LP-INF: 12. largest PC clone",
  "HH117-SILP-INF_clone_nr_13_clone_2477_1": "HH117 SI-LP-INF: 13. largest PC clone",
  "HH117-SILP-INF_clone_nr_14_clone_9145_1": "HH117 SI-LP-INF: 14. largest PC clone",
  "HH117-SILP-INF_clone_nr_15_clone_6575_1": "HH117 SI-LP-INF: 15. largest PC clone",
  "HH117-SILP-INF_clone_nr_16_clone_9593_1": "HH117 SI-LP-INF: 16. largest PC clone",
  "HH117-SILP-INF_clone_nr_17_clone_4511_1": "HH117 SI-LP-INF: 17. largest PC clone",
  "HH117-SILP-INF_clone_nr_18_clone_506_1": "HH117 SI-LP-INF: 18. largest PC clone",
  "HH117-SILP-INF_clone_nr_19_clone_3791_1": "HH117 SI-LP-INF: 19. largest PC clone",
  "HH117-SILP-INF_clone_nr_1_clone_8865_1": "HH117 SI-LP-INF: 1. largest PC clone",
  "HH117-SILP-INF_clone_nr_20_clone_8248_1": "HH117 SI-LP-INF: 20. largest PC clone",
  "HH117-SILP-INF_clone_nr_2_clone_9227_1": "HH117 SI-LP-INF: 2. largest PC clone",
  "HH117-SILP-INF_clone_nr_3_clone_6771_1": "HH117 SI-LP-INF: 3. largest PC clone",
  "HH117-SILP-INF_clone_nr_4_clone_2044_1": "HH117 SI-LP-INF: 4. largest PC clone",
  "HH117-SILP-INF_clone_nr_5_clone_5857_1": "HH117 SI-LP-INF: 5. largest PC clone",
  "HH117-SILP-INF_clone_nr_6_clone_4217_1": "HH117 SI-LP-INF: 6. largest PC clone",
  "HH117-SILP-INF_clone_nr_7_clone_6913_1": "HH117 SI-LP-INF: 7. largest PC clone",
  "HH117-SILP-INF_clone_nr_8_clone_7511_1": "HH117 SI-LP-INF: 8. largest PC clone",
  "HH117-SILP-INF_clone_nr_9_clone_9813_1": "HH117 SI-LP-INF: 9. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_10_clone_8865_1": "HH117 SI-LP-nonINF: 10. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_11_clone_2044_1": "HH117 SI-LP-nonINF: 11. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_12_clone_9012_1": "HH117 SI-LP-nonINF: 12. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_13_clone_10134_1": "HH117 SI-LP-nonINF: 13. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_14_clone_8183_1": "HH117 SI-LP-nonINF: 14. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_15_clone_8760_1": "HH117 SI-LP-nonINF: 15. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_16_clone_3852_1": "HH117 SI-LP-nonINF: 16. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_17_clone_8248_1": "HH117 SI-LP-nonINF: 17. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_18_clone_5261_1": "HH117 SI-LP-nonINF: 18. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_19_clone_6334_1": "HH117 SI-LP-nonINF: 19. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_1_clone_9227_1": "HH117 SI-LP-nonINF: 1. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_20_clone_7656_1": "HH117 SI-LP-nonINF: 20. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_2_clone_6771_1": "HH117 SI-LP-nonINF: 2. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_3_clone_8771_1": "HH117 SI-LP-nonINF: 3. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_4_clone_9502_1": "HH117 SI-LP-nonINF: 4. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_5_clone_8588_1": "HH117 SI-LP-nonINF: 5. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_6_clone_4217_1": "HH117 SI-LP-nonINF: 6. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_7_clone_5857_1": "HH117 SI-LP-nonINF: 7. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_8_clone_9247_1": "HH117 SI-LP-nonINF: 8. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_9_clone_6956_1": "HH117 SI-LP-nonINF: 9. largest PC clone",
  "HH119-COLP_clone_nr_10_clone_23464_1": "HH119 COL-LP: 10. largest PC clone",
  "HH119-COLP_clone_nr_11_clone_28568_1": "HH119 COL-LP: 11. largest PC clone",
  "HH119-COLP_clone_nr_12_clone_28852_1": "HH119 COL-LP: 12. largest PC clone",
  "HH119-COLP_clone_nr_13_clone_11639_1": "HH119 COL-LP: 13. largest PC clone",
  "HH119-COLP_clone_nr_14_clone_25750_1": "HH119 COL-LP: 14. largest PC clone",
  "HH119-COLP_clone_nr_15_clone_26143_1": "HH119 COL-LP: 15. largest PC clone",
  "HH119-COLP_clone_nr_16_clone_28017_1": "HH119 COL-LP: 16. largest PC clone",
  "HH119-COLP_clone_nr_17_clone_3650_1": "HH119 COL-LP: 17. largest PC clone",
  "HH119-COLP_clone_nr_18_clone_1335_1": "HH119 COL-LP: 18. largest PC clone",
  "HH119-COLP_clone_nr_19_clone_2149_1": "HH119 COL-LP: 19. largest PC clone",
  "HH119-COLP_clone_nr_1_clone_24582_1": "HH119 COL-LP: 1. largest PC clone",
  "HH119-COLP_clone_nr_20_clone_25699_1": "HH119 COL-LP: 20. largest PC clone",
  "HH119-COLP_clone_nr_2_clone_4718_1": "HH119 COL-LP: 2. largest PC clone",
  "HH119-COLP_clone_nr_3_clone_4576_1": "HH119 COL-LP: 3. largest PC clone",
  "HH119-COLP_clone_nr_4_clone_14_1": "HH119 COL-LP: 4. largest PC clone",
  "HH119-COLP_clone_nr_5_clone_2472_1": "HH119 COL-LP: 5. largest PC clone",
  "HH119-COLP_clone_nr_6_clone_17505_1": "HH119 COL-LP: 6. largest PC clone",
  "HH119-COLP_clone_nr_7_clone_15065_1": "HH119 COL-LP: 7. largest PC clone",
  "HH119-COLP_clone_nr_8_clone_24963_1": "HH119 COL-LP: 8. largest PC clone",
  "HH119-COLP_clone_nr_9_clone_24974_1": "HH119 COL-LP: 9. largest PC clone",
  "HH119-SILP_clone_nr_10_clone_11664_1": "HH119 SI-LP: 10. largest PC clone",
  "HH119-SILP_clone_nr_11_clone_24338_1": "HH119 SI-LP: 11. largest PC clone",
  "HH119-SILP_clone_nr_12_clone_24405_1": "HH119 SI-LP: 12. largest PC clone",
  "HH119-SILP_clone_nr_13_clone_12918_1": "HH119 SI-LP: 13. largest PC clone",
  "HH119-SILP_clone_nr_14_clone_17590_1": "HH119 SI-LP: 14. largest PC clone",
  "HH119-SILP_clone_nr_15_clone_21586_1": "HH119 SI-LP: 15. largest PC clone",
  "HH119-SILP_clone_nr_16_clone_22404_1": "HH119 SI-LP: 16. largest PC clone",
  "HH119-SILP_clone_nr_17_clone_5624_1": "HH119 SI-LP: 17. largest PC clone",
  "HH119-SILP_clone_nr_18_clone_12704_1": "HH119 SI-LP: 18. largest PC clone",
  "HH119-SILP_clone_nr_19_clone_17998_1": "HH119 SI-LP: 19. largest PC clone",
  "HH119-SILP_clone_nr_1_clone_28238_1": "HH119 SI-LP: 1. largest PC clone",
  "HH119-SILP_clone_nr_20_clone_11580_1": "HH119 SI-LP: 20. largest PC clone",
  "HH119-SILP_clone_nr_2_clone_15994_1": "HH119 SI-LP: 2. largest PC clone",
  "HH119-SILP_clone_nr_3_clone_4157_1": "HH119 SI-LP: 3. largest PC clone",
  "HH119-SILP_clone_nr_4_clone_14647_1": "HH119 SI-LP: 4. largest PC clone",
  "HH119-SILP_clone_nr_5_clone_24427_1": "HH119 SI-LP: 5. largest PC clone",
  "HH119-SILP_clone_nr_6_clone_28541_1": "HH119 SI-LP: 6. largest PC clone",
  "HH119-SILP_clone_nr_7_clone_346_1": "HH119 SI-LP: 7. largest PC clone",
  "HH119-SILP_clone_nr_8_clone_24405_2": "HH119 SI-LP: 8. largest PC clone",
  "HH119-SILP_clone_nr_9_clone_19224_1": "HH119 SI-LP: 9. largest PC clone",
  "HH151-SILP-INF_clone_nr_10_clone_532_1": "HH151 SI-LP-INF: 10. largest PC clone",
  "HH151-SILP-INF_clone_nr_11_clone_8511_1": "HH151 SI-LP-INF: 11. largest PC clone",
  "HH151-SILP-INF_clone_nr_12_clone_366_1": "HH151 SI-LP-INF: 12. largest PC clone",
  "HH151-SILP-INF_clone_nr_13_clone_4131_1": "HH151 SI-LP-INF: 13. largest PC clone",
  "HH151-SILP-INF_clone_nr_14_clone_2935_1": "HH151 SI-LP-INF: 14. largest PC clone",
  "HH151-SILP-INF_clone_nr_15_clone_2599_1": "HH151 SI-LP-INF: 15. largest PC clone",
  "HH151-SILP-INF_clone_nr_16_clone_4747_1": "HH151 SI-LP-INF: 16. largest PC clone",
  "HH151-SILP-INF_clone_nr_17_clone_8000_1": "HH151 SI-LP-INF: 17. largest PC clone",
  "HH151-SILP-INF_clone_nr_18_clone_920_1": "HH151 SI-LP-INF: 18. largest PC clone",
  "HH151-SILP-INF_clone_nr_19_clone_1292_1": "HH151 SI-LP-INF: 19. largest PC clone",
  "HH151-SILP-INF_clone_nr_1_clone_8016_1": "HH151 SI-LP-INF: 1. largest PC clone",
  "HH151-SILP-INF_clone_nr_20_clone_656_1": "HH151 SI-LP-INF: 20. largest PC clone",
  "HH151-SILP-INF_clone_nr_2_clone_2697_1": "HH151 SI-LP-INF: 2. largest PC clone",
  "HH151-SILP-INF_clone_nr_3_clone_6425_1": "HH151 SI-LP-INF: 3. largest PC clone",
  "HH151-SILP-INF_clone_nr_4_clone_6197_1": "HH151 SI-LP-INF: 4. largest PC clone",
  "HH151-SILP-INF_clone_nr_5_clone_6261_1": "HH151 SI-LP-INF: 5. largest PC clone",
  "HH151-SILP-INF_clone_nr_6_clone_1826_1": "HH151 SI-LP-INF: 6. largest PC clone",
  "HH151-SILP-INF_clone_nr_7_clone_8753_1": "HH151 SI-LP-INF: 7. largest PC clone",
  "HH151-SILP-INF_clone_nr_8_clone_2691_1": "HH151 SI-LP-INF: 8. largest PC clone",
  "HH151-SILP-INF_clone_nr_9_clone_440_1": "HH151 SI-LP-INF: 9. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_10_clone_2114_1": "HH151 SI-LP-nonINF: 10. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_11_clone_3650_1": "HH151 SI-LP-nonINF: 11. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_12_clone_1292_1": "HH151 SI-LP-nonINF: 12. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_13_clone_3868_1": "HH151 SI-LP-nonINF: 13. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_14_clone_4040_1": "HH151 SI-LP-nonINF: 14. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_15_clone_4960_1": "HH151 SI-LP-nonINF: 15. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_16_clone_6197_1": "HH151 SI-LP-nonINF: 16. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_17_clone_8277_1": "HH151 SI-LP-nonINF: 17. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_18_clone_4131_1": "HH151 SI-LP-nonINF: 18. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_19_clone_5377_1": "HH151 SI-LP-nonINF: 19. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_1_clone_8124_1": "HH151 SI-LP-nonINF: 1. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_20_clone_5675_1": "HH151 SI-LP-nonINF: 20. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_2_clone_2977_1": "HH151 SI-LP-nonINF: 2. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_3_clone_6091_1": "HH151 SI-LP-nonINF: 3. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_4_clone_8016_1": "HH151 SI-LP-nonINF: 4. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_5_clone_366_1": "HH151 SI-LP-nonINF: 5. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_6_clone_6817_1": "HH151 SI-LP-nonINF: 6. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_7_clone_6059_1": "HH151 SI-LP-nonINF: 7. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_8_clone_2017_1": "HH151 SI-LP-nonINF: 8. largest PC clone",
  "HH151-SILP-nonINF_clone_nr_9_clone_8044_1": "HH151 SI-LP-nonINF: 9. largest PC clone",
  "HH153-SILP-INF_clone_nr_10_clone_2908_1": "HH153 SI-LP-INF: 10. largest PC clone",
  "HH153-SILP-INF_clone_nr_11_clone_2995_1": "HH153 SI-LP-INF: 11. largest PC clone",
  "HH153-SILP-INF_clone_nr_12_clone_3303_1": "HH153 SI-LP-INF: 12. largest PC clone",
  "HH153-SILP-INF_clone_nr_13_clone_1526_1": "HH153 SI-LP-INF: 13. largest PC clone",
  "HH153-SILP-INF_clone_nr_14_clone_1772_1": "HH153 SI-LP-INF: 14. largest PC clone",
  "HH153-SILP-INF_clone_nr_15_clone_1894_1": "HH153 SI-LP-INF: 15. largest PC clone",
  "HH153-SILP-INF_clone_nr_16_clone_1958_1": "HH153 SI-LP-INF: 16. largest PC clone",
  "HH153-SILP-INF_clone_nr_17_clone_1985_1": "HH153 SI-LP-INF: 17. largest PC clone",
  "HH153-SILP-INF_clone_nr_18_clone_2010_1": "HH153 SI-LP-INF: 18. largest PC clone",
  "HH153-SILP-INF_clone_nr_19_clone_3037_1": "HH153 SI-LP-INF: 19. largest PC clone",
  "HH153-SILP-INF_clone_nr_1_clone_1515_1": "HH153 SI-LP-INF: 1. largest PC clone",
  "HH153-SILP-INF_clone_nr_20_clone_3105_1": "HH153 SI-LP-INF: 20. largest PC clone",
  "HH153-SILP-INF_clone_nr_2_clone_1655_1": "HH153 SI-LP-INF: 2. largest PC clone",
  "HH153-SILP-INF_clone_nr_3_clone_2204_1": "HH153 SI-LP-INF: 3. largest PC clone",
  "HH153-SILP-INF_clone_nr_4_clone_2080_1": "HH153 SI-LP-INF: 4. largest PC clone",
  "HH153-SILP-INF_clone_nr_5_clone_717_1": "HH153 SI-LP-INF: 5. largest PC clone",
  "HH153-SILP-INF_clone_nr_6_clone_568_1": "HH153 SI-LP-INF: 6. largest PC clone",
  "HH153-SILP-INF_clone_nr_7_clone_1387_1": "HH153 SI-LP-INF: 7. largest PC clone",
  "HH153-SILP-INF_clone_nr_8_clone_2564_1": "HH153 SI-LP-INF: 8. largest PC clone",
  "HH153-SILP-INF_clone_nr_9_clone_2693_1": "HH153 SI-LP-INF: 9. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_10_clone_2519_1": "HH153 SI-LP-nonINF: 10. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_11_clone_1113_1": "HH153 SI-LP-nonINF: 11. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_12_clone_1452_1": "HH153 SI-LP-nonINF: 12. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_13_clone_2910_1": "HH153 SI-LP-nonINF: 13. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_14_clone_3105_1": "HH153 SI-LP-nonINF: 14. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_15_clone_423_1": "HH153 SI-LP-nonINF: 15. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_16_clone_1020_1": "HH153 SI-LP-nonINF: 16. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_17_clone_1396_1": "HH153 SI-LP-nonINF: 17. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_18_clone_1598_1": "HH153 SI-LP-nonINF: 18. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_19_clone_1666_1": "HH153 SI-LP-nonINF: 19. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_1_clone_1655_1": "HH153 SI-LP-nonINF: 1. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_20_clone_135_1": "HH153 SI-LP-nonINF: 20. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_2_clone_882_1": "HH153 SI-LP-nonINF: 2. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_3_clone_1958_1": "HH153 SI-LP-nonINF: 3. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_4_clone_1772_1": "HH153 SI-LP-nonINF: 4. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_5_clone_818_1": "HH153 SI-LP-nonINF: 5. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_6_clone_1410_1": "HH153 SI-LP-nonINF: 6. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_7_clone_1822_1": "HH153 SI-LP-nonINF: 7. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_8_clone_1526_1": "HH153 SI-LP-nonINF: 8. largest PC clone",
  "HH153-SILP-nonINF_clone_nr_9_clone_2209_1": "HH153 SI-LP-nonINF: 9. largest PC clone",
  
}

########################################################################################
# COLOR BY ANYTHING
########################################################################################

# Define colors 
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cairosvg
from PIL import Image
import io

# Prep colors
color_list = {
    "L1_annotation": {
        "Tfh_cells":             "#E8608A",
        "Naive_Bcells":          "#D4C420",
        "Memory_Bcells":         "#2AAAC8",
        "GC_B_cells":            "#E08C20",
        "PCs":                   "#C42030",
        "Unconventional_Bcells": "#8855CC",
    },
    "c_call_grouped": {
        "IGHA1": "#FF7F00",
        "IGHA2": "#E31A1C",
        "IGHM/D":  "#1F78B4",
        # "IGHD":  "#00E5FF",
        "IGHG1": "#3F007D",
        "IGHG2": "#54278F",
        "IGHG3": "#756BB1",
        "IGHG4": "#9E9AC8"
    },
    "sample_clean_fol": {
    # HH117 samples
    "HH117-SI-MILF-INF":          "#C42030",
    "HH117-SILP-INF":             "#E08C20",

    "HH117-SI-MILF-nonINF":       "#2AAAC8",
    "HH117-SILP-nonINF":          "#1A6090",
    
    "HH117-SI-PP-nonINF_Fol-1":   "#8855CC",
    "HH117-SI-PP-nonINF_Fol-2":   "#A066DD",
    "HH117-SI-PP-nonINF_Fol-3":   "#55CC77",
    "HH117-SI-PP-nonINF_Fol-4":   "#44BB66",
    "HH117-SI-PP-nonINF_Fol-5":   "#33AA55",
    "HH117-SI-PP-nonINF_Fol-6":   "#D4C420",
    "HH117-SI-PP-nonINF_Fol-7":   "#C4B410",
    "HH117-SI-PP-nonINF_Fol-8":   "#E8608A",
    "HH117-SI-PP-nonINF_Fol-9":   "#D85070",
    "HH117-SI-PP-nonINF_Fol-10":  "#EF19EC",
    "HH117-SI-PP-nonINF_Fol-11":  "#CC10CC",
    "HH117-SI-PP-nonINF_Fol-12":  "#FF8C00",
    "HH117-SI-PP-nonINF_Fol-13":  "#20B2AA",
    "HH117-SI-PP-nonINF_Fol-14":  "#10A090",
    "HH117-SI-PP-nonINF_Fol-15":  "#6A5ACD",
    "HH117-SI-PP-nonINF_Fol-16":  "#7B68EE",
    "HH117-SI-PP-nonINF_Fol-17":  "#228B22",
    "HH117-SI-PP-nonINF_Fol-18":  "#32CD32",

    # HH119 samples
    "HH119-CO-SMILF":      "#C42030",
    "HH119-COLP":          "#E08C20",
    "HH119-SILP":          "#1A6090",
    "HH119-SI-MILF":       "#2AAAC8",

    "HH119-SI-PP_Fol-1":   "#8855CC",
    "HH119-SI-PP_Fol-2":   "#A066DD",
    "HH119-SI-PP_Fol-3":   "#55CC77",
    "HH119-SI-PP_Fol-4":   "#44BB66",
    "HH119-SI-PP_Fol-5":   "#33AA55",
    "HH119-SI-PP_Fol-6":   "#D4C420",
    "HH119-SI-PP_Fol-7":   "#C4B410",
    "HH119-SI-PP_Fol-8":   "#E8608A",
    "HH119-SI-PP_Fol-9":   "#D85070",
    "HH119-SI-PP_Fol-10":  "#EF19EC",
    "HH119-SI-PP_Fol-11":  "#CC10CC",
    "HH119-SI-PP_Fol-12":  "#FF8C00",
    "HH119-SI-PP_Fol-13":  "#20B2AA",
    "HH119-SI-PP_Fol-14":  "#10A090",
    "HH119-SI-PP_Fol-15":  "#6A5ACD",
    "HH119-SI-PP_Fol-16":  "#7B68EE",
    "HH119-SI-PP_Fol-17":  "#228B22",
    "HH119-SI-PP_Fol-18":  "#32CD32",
    "HH119-SI-PP_Fol-19":  "#C42030",
    "HH119-SI-PP_Fol-20":  "#E08C20",
    "HH119-SI-PP_Fol-21":  "#2AAAC8",
    "HH119-SI-PP_Fol-22":  "#1A6090",
    "HH119-SI-PP_Fol-23":  "#8855CC",
    "HH119-SI-PP_Fol-24":  "#55CC77",
    "HH119-SI-PP_Fol-25":  "#D4C420",
    "HH119-SI-PP_Fol-26":  "#E8608A",
    "HH119-SI-PP_Fol-27":  "#EF19EC",
    "HH119-SI-PP_Fol-28":  "#FF8C00",
    "HH119-SI-PP_Fol-29":  "#20B2AA",
    "HH119-SI-PP_Fol-30":  "#6A5ACD",
    "HH119-SI-PP_Fol-31":  "#228B22",
    "HH119-SI-PP_Fol-32":  "#A066DD",
    "HH119-SI-PP_Fol-33":  "#44BB66",
    "HH119-SI-PP_Fol-34":  "#C4B410",


    }
    
}

var_translate = {
  "L1_annotation": "Cell type", 
  "c_call_grouped": "Isotype",
  "sample_clean_fol": "Sample"
}

label_translate = {
    "L1_annotation": {
        "Tfh_cells":             "Tfh cells",
        "Naive_Bcells":          "Naive B cells",
        "Memory_Bcells":         "Memory B cells",
        "GC_B_cells":            "GC B cells",
        "PCs":                   "Plasma cells",
        "Unconventional_Bcells": "Unconventional B cells",
    },
    "c_call_grouped": {},        # no translation needed
    "sample_clean_fol": {
      "SI-PP": "PP",
      "SILP": "Ileum LP",
      "COLP": "Colon LP",
      "SI-MILF": "M-ILF",
      "CO-SMILF": "SM-ILF",
      "SILP-nonINF": "Ileum LP (non-INF)",
      "SILP-INF": "Ileum LP (INF)",
      "SI-MILF-nonINF": "M-ILF (non-INF)",
      "SI-MILF-INF": "M-ILF (INF)",
      "SI-PP-nonINF": "PP (non-INF)",
    }  
}

# sample_clean_fol translator
#def format_label(label, HH):
#    # Remove patient prefix e.g. "HH117-"
#    label = label.replace(f"{HH}-", "")
#    label = label.replace(f"SILP", "SI-LP")
#    # Replace "Fol-1" with "Follicle 1"
#    label = re.sub(r".*_Fol-(\d+)", r"Follicle \1", label)
#    return label

# sample_clean_fol translator
def format_label(label, HH):
    # Remove patient prefix e.g. "HH117-"
    label = label.replace(f"{HH}-", "")
    # Replace "Fol-1" with "Follicle 1"
    label = re.sub(r".*_Fol-(\d+)", r"Follicle \1", label)
    # Translate using label_translate if available
    label = label_translate.get("sample_clean_fol", {}).get(label, label)
    return label

for sample, sample_name in samples_dict.items():
  
  print(sample, sample_name)

  HH = sample.split("_")[0]
  plot_path = f"./../plot/{sample}"
  p_file = f"{plot_path}/{sample}.inference.1.p"
  # print(p_file)
  
  # Load an inferred CollapsedTree object from a pickle file (such as are output by the gctree CLI)
  print("Loading data...")
  with open(p_file, "rb") as f:
      tree = pickle.load(f)
  
  # Prep dir for new plots 
  custom_plot_path = f"{plot_path}/costum_trees"
  os.makedirs(custom_plot_path, exist_ok=True)
  
  # The default tree rendering
  tree.render(f"{custom_plot_path}/{sample}_default.png")
  
  # ########################################################################################
  # # Tree size
  # ########################################################################################
  # # Rendering arguments
  # tree.render(f"{custom_plot_path}/{sample}_scale50.png", scale = 50)
  
  # # Change width of tree
  # tree.render(f"{custom_plot_path}/{sample}_scale30_margin20.png", scale = 30, branch_margin=20)
  
  # # Set defult node size
  # tree.render(f"{custom_plot_path}/{sample}_nodesize5.png", node_size=5)

  ########################################################################################
  # Add meta data
  ########################################################################################
  
  # Read meta data
  df_meta = pd.read_csv(f"../gctree_meta/{sample}_gctree_meta.txt",     
                   sep=",",          
                   header=0)
  
  print(df_meta.head())
  
  def plot_tree(tree, df_meta, var, sample, custom_plot_path, color_list, counts_dir=None):
      
      # Get color dict for this variable
      color_dict = color_list[var]
      color_map_with_na = {**color_dict, "NA": "#C8C8C8"}
  
      # Build feature dicts
      dict_str = df_meta.set_index("seq_unique")[var].to_dict()
      dict_int = df_meta.set_index("seq_unique")[f"{var}_int"].to_dict()
  
      # Add features to nodes
      for node in tree.tree.traverse():
          if node.name in dict_str:
              node.add_feature(var, dict_str[node.name])
              node.add_feature(f"{var}_int", dict_int[node.name])
          else:
              node.add_feature(var, "NA")
              node.add_feature(f"{var}_int", np.nan)
  
      # Load counts file if provided (for pie charts with correct abundances)
      counts = None
      if counts_dir is not None:
          counts_path = f"{counts_dir}/{sample}_{var}_counts.csv"
          if os.path.exists(counts_path):
              counts = pd.read_csv(counts_path, index_col=0)
  
      # Build colormap
      def get_color(node, default="#C8C8C8"):
          annotation = getattr(node, var, "NA")
          if annotation is None or annotation == "NA":
              return default
          if counts is not None and node.name in counts.index:
              row = counts.loc[node.name]
              result = {
                  color_map_with_na[col]: count
                  for col, count in row.items()
                  if count > 0 and col in color_map_with_na
              }
              if result:
                  return result
          if ":" in str(annotation):
              parts = [p for p in str(annotation).split(":") if p in color_dict]
              if not parts:
                  return default
              abundance_per_part = node.abundance / len(parts)
              return {color_dict[p]: abundance_per_part for p in parts}
          return color_dict.get(annotation, default)
  
      colormap = {node.name: get_color(node) for node in tree.tree.traverse()}
  
      # Render tree
      svg_path = f"{custom_plot_path}/{sample}_{var}.svg"
      png_path = f"{custom_plot_path}/{sample}_{var}.png"
      tree.render(svg_path, colormap=colormap, scale=15, branch_margin=20)
      cairosvg.svg2png(url=svg_path, write_to=png_path, dpi=150)
  
      # Build legend
      # present_labels = set(
      #     part
      #     for node in tree.tree.traverse()
      #     for part in getattr(node, var, "NA").split(":")
      #     if part != "NA" and part in color_dict
      # )
      present_labels = set(
          part
          for node in tree.tree.traverse()
          for part in str(getattr(node, var, "NA")).split(":")
          if part != "NA" and part != "nan" and part in color_dict
      )
  
      translate = label_translate.get(var, {})
      
      patches = [
          mpatches.Patch(
            color=color_dict[label], 
            label=translate.get(label, format_label(label, HH))  # translate if available, else format_label
          )
          for label in present_labels
          if label in color_dict
      ]
  
      # Add NA to legend if present
      has_na = (
          counts["NA"].sum() > 0 if counts is not None and "NA" in counts.columns
          else any(
              "NA" in getattr(node, var, "").split(":")
              for node in tree.tree.traverse()
          )
      )
      if has_na:
          patches.append(mpatches.Patch(color="#C8C8C8", label="NA"))
  
      # Add legend and title to plot
      img = Image.open(png_path)
      fig, ax = plt.subplots(figsize=(img.width/200, img.height/200))
      ax.imshow(img)
      ax.axis("off")
      ax.legend(handles=patches, loc="upper left", title=var_translate[var], fontsize=22, title_fontsize=24)
      ax.set_title(sample_name, fontsize=26, fontweight="bold")
      plt.savefig(png_path, bbox_inches="tight", dpi=150)
      plt.close()
  
      print(f"Saved: {png_path}")
  
  # Color by L1_annotation
  plot_tree(tree, df_meta, "L1_annotation", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")
  
  # Color by c_call_grouped (with counts for correct pie charts)
  plot_tree(tree, df_meta, "c_call_grouped", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")
  
  # Color by sample_clean_fol (with counts for correct pie charts)
  plot_tree(tree, df_meta, "sample_clean_fol", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")
