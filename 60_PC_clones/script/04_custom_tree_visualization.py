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

  # "HH117-SILP-INF_clone_nr_10_clone_718_1": "CD SI-LP-INF: 10. largest PC clone",
  # "HH117-SILP-INF_clone_nr_1_clone_8985_1": "CD SI-LP-INF: 1. largest PC clone",
  # "HH117-SILP-INF_clone_nr_2_clone_455_1": "CD SI-LP-INF: 2. largest PC clone",
  # "HH117-SILP-INF_clone_nr_3_clone_6469_1": "CD SI-LP-INF: 3. largest PC clone",
  # "HH117-SILP-INF_clone_nr_4_clone_9677_1": "CD SI-LP-INF: 4. largest PC clone",
  # "HH117-SILP-INF_clone_nr_5_clone_6896_1": "CD SI-LP-INF: 5. largest PC clone",
  # # "HH117-SILP-INF_clone_nr_6_clone_7587_1": "CD SI-LP-INF: 6. largest PC clone",
  # "HH117-SILP-INF_clone_nr_7_clone_2275_1": "CD SI-LP-INF: 7. largest PC clone",
  # "HH117-SILP-INF_clone_nr_8_clone_24_1": "CD SI-LP-INF: 8. largest PC clone",
  # "HH117-SILP-INF_clone_nr_9_clone_9150_1": "CD SI-LP-INF: 9. largest PC clone",
  # "HH117-SILP-nonINF_clone_nr_10_clone_8985_1": "CD SI-LP-nonINF: 10. largest PC clone",
  # "HH117-SILP-nonINF_clone_nr_1_clone_6469_1": "CD SI-LP-nonINF: 1. largest PC clone",
  # "HH117-SILP-nonINF_clone_nr_2_clone_8879_1": "CD SI-LP-nonINF: 2. largest PC clone",
  # "HH117-SILP-nonINF_clone_nr_3_clone_455_1": "CD SI-LP-nonINF: 3. largest PC clone",
  # "HH117-SILP-nonINF_clone_nr_4_clone_8700_1": "CD SI-LP-nonINF: 4. largest PC clone",
  # "HH117-SILP-nonINF_clone_nr_5_clone_9624_1": "CD SI-LP-nonINF: 5. largest PC clone",
  # # "HH117-SILP-nonINF_clone_nr_6_clone_9347_1": "CD SI-LP-nonINF: 6. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_7_clone_9677_1": "CD SI-LP-nonINF: 7. largest PC clone",
  "HH117-SILP-nonINF_clone_nr_8_clone_6961_1": "CD SI-LP-nonINF: 8. largest PC clone",
  # "HH117-SILP-nonINF_clone_nr_9_clone_8888_1": "CD SI-LP-nonINF: 9. largest PC clone",

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
  
  # Color by c_call (with counts for correct pie charts)
  plot_tree(tree, df_meta, "c_call", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")
  
  # Color by sample_clean_fol (with counts for correct pie charts)
  plot_tree(tree, df_meta, "sample_clean_fol", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")
