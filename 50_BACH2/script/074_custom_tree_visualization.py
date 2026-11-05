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
  # "HH117_clone_nr_1_clone_1881_1": "HH117: Clone 1",
  # "HH117_clone_nr_2_clone_1951_1": "HH117: Clone 2",
  # "HH117_clone_nr_3_clone_2667_1": "HH117: Clone 3",
  # "HH117_clone_nr_4_clone_4261_1": "HH117: Clone 4",
  # "HH117_clone_nr_5_clone_5139_1": "HH117: Clone 5",
  # "HH117_clone_nr_6_clone_5959_1": "HH117: Clone 6",
  # "HH117_clone_nr_7_clone_4066_1": "HH117: Clone 7",
  # # "HH117_clone_nr_8_clone_5608_1": "HH117: Clone 8", # One sequence --> no tree
  "HH117_clone_nr_9_clone_5834_1": "HH117: Clone 9",
  "HH117_clone_nr_10_clone_6270_1": "HH117: Clone 10",
  "HH117_clone_nr_11_clone_715_1": "HH117: Clone 11",
  "HH117_clone_nr_12_clone_7509_1": "HH117: Clone 12",
  "HH117_clone_nr_13_clone_9817_1": "HH117: Clone 13",

  # "HH119_clone_nr_1_clone_27383_1": "HH119: Clone 1",
  # "HH119_clone_nr_2_clone_4494_1": "HH119: Clone 2",
  # "HH119_clone_nr_3_clone_29181_1": "HH119: Clone 3",
  # "HH119_clone_nr_4_clone_8598_1": "HH119: Clone 4",
  # "HH119_clone_nr_5_clone_15039_1": "HH119: Clone 5",
  # "HH119_clone_nr_6_clone_26043_1": "HH119: Clone 6",
  # "HH119_clone_nr_7_clone_8506_1": "HH119: Clone 7",
  # "HH119_clone_nr_8_clone_22295_1": "HH119: Clone 8",
  # # "HH119_clone_nr_9_clone_17585_1": "HH119: Clone 9", # One sequence --> no tree
  "HH119_clone_nr_10_clone_8457_1": "HH119: Clone 10",
  "HH119_clone_nr_11_clone_22776_1": "HH119: Clone 11",
  "HH119_clone_nr_12_clone_26698_1": "HH119: Clone 12",
  "HH119_clone_nr_13_clone_28739_1": "HH119: Clone 13",
  "HH119_clone_nr_14_clone_3804_1": "HH119: Clone 14",
  "HH119_clone_nr_15_clone_8177_1": "HH119: Clone 15",
  "HH119_clone_nr_16_clone_24431_1": "HH119: Clone 16",
  "HH119_clone_nr_17_clone_28674_1": "HH119: Clone 17",
  "HH119_clone_nr_18_clone_1830_1": "HH119: Clone 18",
  "HH119_clone_nr_19_clone_24483_1": "HH119: Clone 19",
  "HH119_clone_nr_20_clone_8591_1": "HH119: Clone 20",
  "HH119_clone_nr_21_clone_903_1": "HH119: Clone 21",
  "HH119_clone_nr_22_clone_19562_1": "HH119: Clone 22",
  "HH119_clone_nr_23_clone_22252_1": "HH119: Clone 23",
  "HH119_clone_nr_24_clone_23494_1": "HH119: Clone 24",
  "HH119_clone_nr_25_clone_25333_1": "HH119: Clone 25",
  "HH119_clone_nr_26_clone_10201_1": "HH119: Clone 26",
  "HH119_clone_nr_27_clone_17585_2": "HH119: Clone 27",
  "HH119_clone_nr_28_clone_17807_1": "HH119: Clone 28",
  "HH119_clone_nr_29_clone_26121_1": "HH119: Clone 29",
  "HH119_clone_nr_30_clone_2812_1": "HH119: Clone 30",
  "HH119_clone_nr_31_clone_5037_1": "HH119: Clone 31",
  "HH119_clone_nr_32_clone_5652_1": "HH119: Clone 32",
  "HH119_clone_nr_33_clone_5836_1": "HH119: Clone 33",
  "HH119_clone_nr_34_clone_8685_1": "HH119: Clone 34",

  "HH151_clone_nr_1_clone_3084_1": "HH151: Clone 1",
  "HH151_clone_nr_2_clone_4747_1": "HH151: Clone 2",
  "HH151_clone_nr_3_clone_6909_1": "HH151: Clone 3",

  "HH153_clone_nr_1_clone_188_1": "HH153: Clone 1",
  "HH153_clone_nr_2_clone_3105_1": "HH153: Clone 2",
  "HH153_clone_nr_3_clone_2204_1": "HH153: Clone 3",
  "HH153_clone_nr_4_clone_125_1": "HH153: Clone 4",
  "HH153_clone_nr_5_clone_3058_1": "HH153: Clone 5",
  "HH153_clone_nr_6_clone_2209_1": "HH153: Clone 6",
  "HH153_clone_nr_7_clone_3107_1": "HH153: Clone 7",
  "HH153_clone_nr_8_clone_1399_1": "HH153: Clone 8",
  "HH153_clone_nr_9_clone_3100_1": "HH153: Clone 9",
  "HH153_clone_nr_10_clone_3482_1": "HH153: Clone 10",
  "HH153_clone_nr_11_clone_234_1": "HH153: Clone 11",
  "HH153_clone_nr_12_clone_324_1": "HH153: Clone 12",
  "HH153_clone_nr_13_clone_149_1": "HH153: Clone 13",
  "HH153_clone_nr_14_clone_3526_1": "HH153: Clone 14",
  "HH153_clone_nr_15_clone_2512_1": "HH153: Clone 15",
  "HH153_clone_nr_16_clone_1283_1": "HH153: Clone 16",
  "HH153_clone_nr_17_clone_215_1": "HH153: Clone 17",
  "HH153_clone_nr_18_clone_1841_1": "HH153: Clone 18",
  "HH153_clone_nr_19_clone_2551_1": "HH153: Clone 19",
  "HH153_clone_nr_20_clone_1221_1": "HH153: Clone 20",
  "HH153_clone_nr_21_clone_1559_1": "HH153: Clone 21",
  "HH153_clone_nr_22_clone_2045_1": "HH153: Clone 22",
  "HH153_clone_nr_23_clone_1227_1": "HH153: Clone 23",
  "HH153_clone_nr_24_clone_1268_1": "HH153: Clone 24",
  "HH153_clone_nr_25_clone_1402_1": "HH153: Clone 25",
  "HH153_clone_nr_26_clone_2601_1": "HH153: Clone 26",
  "HH153_clone_nr_27_clone_2934_1": "HH153: Clone 27",
}



# for sample in samples:


# Define samples 
# sample = "HH117_clone_nr_1_clone_4221_1"
# sample_name = "Crohn's Diease: Largest GC B cell clone"

# sample = "HH119_clone_nr_2_clone_5791_1"
# sample_name = "Colorectal Cancer: Second largest clone with GC B cells across samples"


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
    
   "L1_L3_annotation": {
      # L1 annotation
      "Tfh_cells":              "#E8608A",
      "Naive_B_cells":          "#D4C420",
      "Memory_B_cells":         "#2AAAC8",
      "GC_B_cells":             "#E08C20",
      "PCs":                    "#C42030",
      "Unconventional_B_cells": "#8855CC",
  
      # GC B cell subsets
      "DZ_cycling":             "#6A98C4",
      "DZ":                     "#A9CDE0",
      "Pre_MBC_cycling":        "#A485C6",
      "LZ":                     "#7BA270",
      "LZ_transitioning":       "#A3C293",
      "BACH2+":                 "#E77F6D",
      "BACH2_cycling":          "#BA6052",
      "Pre_MBC":                "#5BA4A4",
      "Pre_PB":                 "#745973",
  
      # Clusters
      "0":                      "#E41A1C",
      "1":                      "#FF8C00",
      "2":                      "#4DAF4A",
      "3":                      "#00A896",
      "4":                      "#377EB8",
      "5":                      "#984EA3",
      "6":                      "#E7298A",
    },
    
    
    # "c_call": {
    #     "IGHA1": "#FF7F00",
    #     "IGHA2": "#E31A1C",
    #     "IGHM":  "#1F78B4",
    #     "IGHD":  "#00E5FF",
    #     "IGHG1": "#3F007D",
    #     "IGHG2": "#54278F",
    #     "IGHG3": "#756BB1",
    #     "IGHG4": "#9E9AC8"
    # },
    "c_call_grouped": {
        "IGHA1": "#FF7F00",
        "IGHA2": "#E31A1C",
        "IGHM/D":  "#1F78B4",
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
    
        # HH153
    "HH153-SILP-INF":           "#1A6090",
    "HH153-SILP-nonINF":        "#2AAAC8",
    "HH153-SI-PP-nonINF_Fol-1":  "#228B22",
    "HH153-SI-PP-nonINF_Fol-2":  "#A066DD",
    "HH153-SI-PP-nonINF_Fol-3":  "#44BB66",
    "HH153-SI-PP-nonINF_Fol-4":  "#C4B410",
    "HH153-SI-PP-nonINF_Fol-5":  "#1F77B4",
    "HH153-SI-PP-nonINF_Fol-6":  "#E377C2",
    "HH153-SI-PP-nonINF_Fol-7":  "#17BECF",
    "HH153-SI-PP-nonINF_Fol-8":  "#D62728",
    "HH153-SI-PP-nonINF_Fol-9":  "#8C564B",
    "HH153-SI-PP-nonINF_Fol-10": "#FF7F0E",
    "HH153-SI-PP-nonINF_Fol-11": "#7B68EE",
    "HH153-SI-PP-nonINF_Fol-12": "#9ACD32",
    "HH153-SI-PP-nonINF_Fol-13": "#FF69B4",
    "HH153-SI-PP-nonINF_Fol-14": "#00CED1",
    "HH153-SI-PP-nonINF_Fol-15": "#B8860B",
    "HH153-SI-PP-nonINF_Fol-16": "#4B0082",
    "HH153-SI-PP-nonINF_Fol-17": "#2E8B57",
    "HH153-SI-PP-nonINF_Fol-18": "#FF4500",
    "HH153-SI-PP-nonINF_Fol-19": "#4682B4",
    "HH153-SI-PP-nonINF_Fol-20": "#DA70D6",

    # HH151
    "HH151-SILP-INF":           "#1A6090",
    "HH151-SILP-nonINF":        "#2AAAC8",
    "HH151-SI-PP-nonINF_Fol-1":  "#228B22",
    "HH151-SI-PP-nonINF_Fol-2":  "#A066DD",
    "HH151-SI-PP-nonINF_Fol-3":  "#44BB66",
    "HH151-SI-PP-nonINF_Fol-4":  "#C4B410",
    "HH151-SI-PP-nonINF_Fol-5":  "#1F77B4",
    "HH151-SI-PP-nonINF_Fol-6":  "#E377C2",
    "HH151-SI-PP-nonINF_Fol-7":  "#17BECF",
    "HH151-SI-PP-nonINF_Fol-8":  "#D62728",
    "HH151-SI-PP-nonINF_Fol-9":  "#8C564B",
    "HH151-SI-PP-nonINF_Fol-10": "#FF7F0E",
    "HH151-SI-PP-nonINF_Fol-11": "#7B68EE",
    "HH151-SI-PP-nonINF_Fol-12": "#9ACD32",
    "HH151-SI-PP-nonINF_Fol-13": "#FF69B4",
    "HH151-SI-PP-nonINF_Fol-14": "#00CED1",
    "HH151-SI-PP-nonINF_Fol-15": "#B8860B",
    "HH151-SI-PP-nonINF_Fol-16": "#4B0082",
    "HH151-SI-PP-nonINF_Fol-17": "#2E8B57",
    "HH151-SI-PP-nonINF_Fol-18": "#FF4500",
    "HH151-SI-PP-nonINF_Fol-19": "#4682B4",
    "HH151-SI-PP-nonINF_Fol-20": "#DA70D6",


    }
    
}

var_translate = {
  "L1_annotation": "Cell type", 
  "L1_L3_annotation": "Cell type (GC B subsets + PC clusters)", 
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
    "L1_L3_annotation": {
      
     "Tfh_cells":             "Tfh cells",
        "Naive_Bcells":          "Naive B cells",
        "Memory_Bcells":         "Memory B cells",
        "GC_B_cells":            "GC B cells",
        "PCs":                   "Plasma cells",
        "Unconventional_Bcells": "Unconventional B cells",
  
      # GC B cell subsets
      "DZ_cycling":             "DZ cycling",      
      "DZ":                     "DZ",              
      "Pre_MBC_cycling":        "Pre MBC cycling", 
      "LZ":                     "LZ",             
      "LZ_transitioning":       "LZ transitioning",
      "BACH2+":                 "BACH2+",          
      "BACH2_cycling":          "BACH2 cycling",   
      "Pre_MBC":                "Pre MBC",         
      "Pre_PB":                 "Pre PB",         
  
      # Clusters
      "0":                      "PC cluster 1",
      "1":                      "PC cluster 2",
      "2":                      "PC cluster 3",
      "3":                      "PC cluster 4",
      "4":                      "PC cluster 5",
      "5":                      "PC cluster 6",
      "6":                      "PC cluster 7",
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
  plot_tree(tree, df_meta, "L1_annotation", sample, custom_plot_path, color_list, counts_dir=f"../gctree_meta")
  
  # Color by L3_annotation
  plot_tree(tree, df_meta, "L1_L3_annotation", sample, custom_plot_path, color_list, counts_dir=f"../gctree_meta")
  
  # Color by c_call (with counts for correct pie charts)
  plot_tree(tree, df_meta, "c_call_grouped", sample, custom_plot_path, color_list, counts_dir=f"../gctree_meta")
  
  # Color by sample_clean_fol (with counts for correct pie charts)
  plot_tree(tree, df_meta, "sample_clean_fol", sample, custom_plot_path, color_list, counts_dir=f"../gctree_meta")
