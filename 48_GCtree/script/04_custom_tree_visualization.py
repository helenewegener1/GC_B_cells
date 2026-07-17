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
  
  #"HH117_clone_nr_1_clone_500_1":  "CD: Largest GC B cell clone",
  #"HH117_clone_nr_2_clone_2885_1":  "CD: 2. largest GC B cell clone",
  #"HH117_clone_nr_3_clone_2018_1":  "Crohn's disease: 3. largest GC B cell clone",
  #"HH117_clone_nr_4_clone_3953_1":  "CD: 4. largest GC B cell clone",
  #"HH117_clone_nr_5_clone_2588_1":  "CD: 5. largest GC B cell clone",
  # "HH117_clone_nr_6_clone_1571_1":  "CD: 6. largest GC B cell clone",
  # "HH117_clone_nr_7_clone_5554_1":  "CD: 7. largest GC B cell clone",
  # "HH117_clone_nr_8_clone_5725_1":  "CD: 8. largest GC B cell clone",
  # "HH117_clone_nr_9_clone_2211_1":  "CD: 9. largest GC B cell clone",
  # "HH117_clone_nr_10_clone_2457_1": "CD: 10. largest GC B cell clone",
  "HH117_clone_nr_11_clone_5904_1": "CD: 11. largest GC B cell clone",
  "HH117_clone_nr_12_clone_7177_1": "CD: 12. largest GC B cell clone",
  "HH117_clone_nr_13_clone_6361_1": "CD: 13. largest GC B cell clone",
  "HH117_clone_nr_14_clone_972_1":  "CD: 14. largest GC B cell clone",
  "HH117_clone_nr_15_clone_9669_1": "CD: 15. largest GC B cell clone",
  "HH117_clone_nr_16_clone_9994_1": "CD: 16. largest GC B cell clone",
  "HH117_clone_nr_17_clone_5459_1": "CD: 17. largest GC B cell clone",
  "HH117_clone_nr_18_clone_5680_1": "CD: 18. largest GC B cell clone",
  "HH117_clone_nr_19_clone_7466_1": "CD: 19. largest GC B cell clone",
  "HH117_clone_nr_20_clone_2140_1": "CD: 20. largest GC B cell clone",
  
  # # "HH119_clone_nr_1_clone_22355_1": "CRC: Largest GC B cell clone",
  # "HH119_clone_nr_2_clone_27035_1": "CRC: 2. largest GC B cell clone",
  # "HH119_clone_nr_3_clone_11978_1": "Colorectal cancer: 3. largest GC B cell clone",
  # "HH119_clone_nr_4_clone_7879_1": "CRC: 4. largest GC B cell clone",
  # "HH119_clone_nr_5_clone_14912_1": "CRC: 5. largest GC B cell clone",
  # "HH119_clone_nr_6_clone_21402_1": "CRC: 6. largest GC B cell clone",
  # "HH119_clone_nr_7_clone_8286_1":  "CRC: 7. largest GC B cell clone",
  # "HH119_clone_nr_8_clone_8625_1":  "CRC: 8. largest GC B cell clone",
  # "HH119_clone_nr_9_clone_3913_1": "CRC: 9. largest GC B cell clone",
  # "HH119_clone_nr_10_clone_23719_1": "CRC: 10. largest GC B cell clone",
  "HH119_clone_nr_11_clone_8576_1": "CRC: 11. largest GC B cell clone",
  "HH119_clone_nr_12_clone_145_1": "CRC: 12. largest GC B cell clone",
  "HH119_clone_nr_13_clone_8374_1": "CRC: 13. largest GC B cell clone",
  "HH119_clone_nr_14_clone_16064_1": "CRC: 14. largest GC B cell clone",
  "HH119_clone_nr_15_clone_10590_1": "CRC: 15. largest GC B cell clone",
  "HH119_clone_nr_16_clone_8980_1": "CRC: 16. largest GC B cell clone",
  "HH119_clone_nr_17_clone_25709_1": "CRC: 17. largest GC B cell clone",
  "HH119_clone_nr_18_clone_26412_1": "CRC: 18. largest GC B cell clone",
  "HH119_clone_nr_19_clone_10599_1": "CRC: 19. largest GC B cell clone",
  "HH119_clone_nr_20_clone_27222_1": "CRC: 20. largest GC B cell clone",
}


# for sample in samples:


# Define samples 
# sample = "HH117_clone_nr_1_clone_4221_1"
# sample_name = "Crohn's Diease: Largest GC B cell clone"

# sample = "HH119_clone_nr_2_clone_5791_1"
# sample_name = "Colorectal Cancer: Second largest clone with GC B cells across samples"

# version = ""
# version = "gmm_threshold"
version = "90_similarity"

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
    "c_call": {
        "IGHA1": "#FF7F00",
        "IGHA2": "#E31A1C",
        "IGHM":  "#1F78B4",
        "IGHD":  "#00E5FF",
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
  "c_call": "Isotype",
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
    "c_call": {},        # no translation needed
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
  plot_path = f"./../plot_{version}/{sample}"
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
  df_meta = pd.read_csv(f"../gctree_meta_{version}/{sample}_gctree_meta.txt",     
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
  plot_tree(tree, df_meta, "L1_annotation", sample, custom_plot_path, color_list, counts_dir=f"../gctree_meta_{version}")
  
  # Color by c_call (with counts for correct pie charts)
  plot_tree(tree, df_meta, "c_call", sample, custom_plot_path, color_list, counts_dir=f"../gctree_meta_{version}")
  
  # Color by sample_clean_fol (with counts for correct pie charts)
  plot_tree(tree, df_meta, "sample_clean_fol", sample, custom_plot_path, color_list, counts_dir=f"../gctree_meta_{version}")
