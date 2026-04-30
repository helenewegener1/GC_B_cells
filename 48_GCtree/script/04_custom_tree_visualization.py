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

# Define paths 
sample = "HH117_clone_nr_5_clone_1856_1"
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

# ########################################################################################
# # Color nodes
# ########################################################################################

# colormap = tree.feature_colormap("abundance")
# tree.render(f"{custom_plot_path}/{sample}_color_abundance.png" , colormap=colormap)


########################################################################################
# Add meta data
########################################################################################

# Read meta data
df_meta = pd.read_csv(f"../gctree_meta/{sample}_gctree_meta.txt",     
                 sep=",",          
                 header=0)

print(df_meta.head())

# Define colors 
from matplotlib.colors import ListedColormap
import matplotlib as mpl

########################################################################################
# Color by L1 annotation 
########################################################################################

# # # Make dict
# # L1_dict = df_meta.set_index("seq_unique")["L1_annotation_int"].to_dict()

# # Add features
# # for node in tree.tree.traverse():
# #     if node.name in L1_dict:
# #         # print(node.name)
# #         # print(L1_dict[node.name])
# #         node.add_feature("L1_annotation_int", L1_dict[node.name])
# #     else:
# #         node.add_feature("L1_annotation_int", np.nan)

# L1_dict_str = df_meta.set_index("seq_unique")["L1_annotation"].to_dict()
# L1_dict_int = df_meta.set_index("seq_unique")["L1_annotation_int"].to_dict()

# for node in tree.tree.traverse():
#     if node.name in L1_dict_str:
#         print(f"Adding feature to {node.name}: {L1_dict_str[node.name]}")  # check it's hitting
#         node.add_feature("L1_annotation", L1_dict_str[node.name])
#         node.add_feature("L1_annotation_int", L1_dict_int[node.name])
#     else:
#         node.add_feature("L1_annotation", "NA")
#         node.add_feature("L1_annotation_int", np.nan)


# # Define your colors
# color_dict = {
#     "Tfh_cells":             "#E8608A",
#     "Naive_Bcells":          "#D4C420",
#     "Memory_Bcells":         "#2AAAC8",
#     "GC_B_cells":            "#E08C20",
#     "PCs":                   "#C42030",
#     "Unconventional_Bcells": "#8855CC",
# }

# # Build mapping (cell type --> number) from df_meta
# category_mapping = dict(zip(
#     df_meta["L1_annotation"],
#     df_meta["L1_annotation_int"].astype(int)
# ))
# # e.g. {"PCs": 2, "Memory_Bcells": 1}

# category_mapping_inv = {v: k for k, v in category_mapping.items()}
# # e.g. {2: "PCs", 1: "Memory_Bcells"}

# # Build colormap
# # def get_color(node, color_dict, category_mapping_inv, default="#808080"):
# #     val = getattr(node, "L1_annotation_int", np.nan)
# #     if val is None or (isinstance(val, float) and np.isnan(val)):
# #         return default
# #     label = category_mapping_inv.get(int(val), None)
# #     if label is None:
# #         return default
# #     return color_dict.get(label, default)

# def get_color(node, color_dict, default="#808080"):
#     annotation = getattr(node, "L1_annotation", "NA")
    
#     if annotation is None or annotation == "NA":
#         return default
    
#     if ":" in str(annotation):
#         parts = [p for p in str(annotation).split(":") if p in color_dict]
#         if not parts:
#             return default
#         # Use abundance per part instead of frequency
#         abundance_per_part = node.abundance / len(parts)
#         return {color_dict[p]: abundance_per_part for p in parts}
    
#     return color_dict.get(annotation, default)

# colormap = {
#     node.name: get_color(node, color_dict)
#     for node in tree.tree.traverse()
# }

# print(colormap)


# # RENDER TREE
# # colormap = tree.feature_colormap("L1_annotation_int", cmap="plasma")
# # tree.render(f"{custom_plot_path}/{sample}_L1_annotation_int.svg" , colormap=colormap, scale = 10, branch_margin=2)
# # print(help(tree.feature_colormap))

# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# import cairosvg
# from PIL import Image
# import io

# # Render tree to SVG
# svg_path = f"{custom_plot_path}/{sample}_L1_annotation_int.svg"
# png_path = f"{custom_plot_path}/{sample}_L1_annotation_int.png"
# tree.render(svg_path, colormap=colormap, scale=20, branch_margin=10)

# # Convert SVG to PNG
# cairosvg.svg2png(url=svg_path, write_to=png_path, dpi=150)

# # # Build legend patches only for categories present in this tree
# # present_labels = set(
# #     category_mapping_inv.get(int(getattr(node, "L1_annotation_int", np.nan)), None)
# #     for node in tree.tree.traverse()
# #     if not (isinstance(getattr(node, "L1_annotation_int", np.nan), float) and 
# #             np.isnan(getattr(node, "L1_annotation_int", np.nan)))
# # )
# # present_labels.discard(None)

# # patches = [
# #     mpatches.Patch(color=color_dict[label], label=label)
# #     for label in present_labels
# #     if label in color_dict
# # ]

# # Get present labels from the tree
# present_labels = set(
#     part
#     for node in tree.tree.traverse()
#     for part in getattr(node, "L1_annotation", "NA").split(":")
#     if part != "NA" and part in color_dict
# )

# print("Present labels:", present_labels)

# # Build legend patches
# patches = [
#     mpatches.Patch(color=color_dict[label], label=label)
#     for label in present_labels
#     if label in color_dict
# ]

# # Open PNG and add legend
# img = Image.open(png_path)
# fig, ax = plt.subplots(figsize=(img.width/200, img.height/200))
# ax.imshow(img)
# ax.axis("off")
# ax.legend(handles=patches, loc="upper left", title="Cell type", fontsize=14, title_fontsize=16)
# ax.set_title(f"GC tree of {sample}", fontsize=20, fontweight="bold")
# plt.savefig(png_path, bbox_inches="tight", dpi=150)
# plt.close()

# ########################################################################################
# # COLOR BY ISOTYPE
# ########################################################################################

# c_call_dict_str = df_meta.set_index("seq_unique")["c_call"].to_dict()
# c_call_dict_int = df_meta.set_index("seq_unique")["c_call_int"].to_dict()

# for node in tree.tree.traverse():
#     if node.name in c_call_dict_str:
#         print(f"Adding feature to {node.name}: {c_call_dict_str[node.name]}")  # check it's hitting
#         node.add_feature("c_call", c_call_dict_str[node.name])
#         node.add_feature("c_call_int", c_call_dict_int[node.name])
#     else:
#         node.add_feature("c_call", "NA")
#         node.add_feature("c_call_int", np.nan)


# # Define your colors
# color_dict = {
#     "IGHA1": "#E8608A",
#     "IGHA2": "#D4C420",
#     "IGHM":  "#2AAAC8",
#     "IGHD":  "#E08C20",
#     "IGHG1": "#C42030",
#     "IGHG2": "#8855CC",
#     "IGHG3": "#55CC77",
#     "IGHG4": "#EF19EC"
# }

# # Build mapping (cell type --> number) from df_meta
# category_mapping = dict(zip(
#     df_meta["c_call"],
#     df_meta["c_call_int"].astype(int)
# ))
# # e.g. {"PCs": 2, "Memory_Bcells": 1}

# category_mapping_inv = {v: k for k, v in category_mapping.items()}
# # e.g. {2: "PCs", 1: "Memory_Bcells"}

# # Build colormap
# isotype_counts = pd.read_csv(f"../gctree_meta/{sample}_isotype_counts.csv",
#                              index_col=0,
#                              sep=",",          
#                              header=0)

# print(isotype_counts.head())

# color_map_with_na = {**color_dict, "NA": "#808080"}

# def get_color(node, color_dict, isotype_counts, default="#808080"):
#     annotation = getattr(node, "c_call", "NA")
    
#     if annotation is None or annotation == "NA":
#         return default
    
#     if node.name in isotype_counts.index:
#         row = isotype_counts.loc[node.name]
#         result = {
#             color_map_with_na[col]: count
#             for col, count in row.items()
#             if count > 0 and col in color_map_with_na
#         }
#         if result:
#             return result
    
#     return color_dict.get(annotation, default)

# colormap = {
#     node.name: get_color(node, color_dict, isotype_counts)
#     for node in tree.tree.traverse()
# }

# # print(colormap)

# # RENDER TREE
# # Render tree to SVG
# svg_path = f"{custom_plot_path}/{sample}_c_call_int.svg"
# png_path = f"{custom_plot_path}/{sample}_c_call_int.png"
# tree.render(svg_path, colormap=colormap, scale=20, branch_margin=20)

# # Convert SVG to PNG
# cairosvg.svg2png(url=svg_path, write_to=png_path, dpi=150)

# # Get present labels from the tree
# present_labels = set(
#     part
#     for node in tree.tree.traverse()
#     for part in getattr(node, "c_call", "NA").split(":")
#     if part != "NA" and part in color_dict
# )

# print("Present labels:", present_labels)

# # Check if any node has NA counts
# has_na = isotype_counts["NA"].sum() > 0 if "NA" in isotype_counts.columns else False

# # Build legend patches
# patches = [
#     mpatches.Patch(color=color_dict[label], label=label)
#     for label in present_labels
#     if label in color_dict
# ]

# # Add NA patch if present
# if has_na:
#     patches.append(mpatches.Patch(color="#808080", label="NA"))

# # Open PNG and add legend
# img = Image.open(png_path)
# fig, ax = plt.subplots(figsize=(img.width/200, img.height/200))
# ax.imshow(img)
# ax.axis("off")
# ax.legend(handles=patches, loc="upper left", title="Cell type", fontsize=14, title_fontsize=16)
# ax.set_title(f"GC tree of {sample}", fontsize=20, fontweight="bold")
# plt.savefig(png_path, bbox_inches="tight", dpi=150)
# plt.close()



########################################################################################
# COLOR BY ANYTHING
########################################################################################

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
        "IGHA1": "#E8608A",
        "IGHA2": "#D4C420",
        "IGHM":  "#2AAAC8",
        "IGHD":  "#E08C20",
        "IGHG1": "#C42030",
        "IGHG2": "#8855CC",
        "IGHG3": "#55CC77",
        "IGHG4": "#EF19EC"
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

def plot_tree(tree, df_meta, var, sample, custom_plot_path, color_list, counts_dir=None):
    
    # Get color dict for this variable
    color_dict = color_list[var]
    color_map_with_na = {**color_dict, "NA": "#808080"}

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
    def get_color(node, default="#808080"):
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
    tree.render(svg_path, colormap=colormap, scale=20, branch_margin=20)
    cairosvg.svg2png(url=svg_path, write_to=png_path, dpi=150)

    # Build legend
    present_labels = set(
        part
        for node in tree.tree.traverse()
        for part in getattr(node, var, "NA").split(":")
        if part != "NA" and part in color_dict
    )

    patches = [
        mpatches.Patch(color=color_dict[label], label=label)
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
        patches.append(mpatches.Patch(color="#808080", label="NA"))

    # Add legend and title to plot
    img = Image.open(png_path)
    fig, ax = plt.subplots(figsize=(img.width/200, img.height/200))
    ax.imshow(img)
    ax.axis("off")
    ax.legend(handles=patches, loc="upper left", title=var, fontsize=14, title_fontsize=16)
    ax.set_title(f"GC tree of {sample}", fontsize=20, fontweight="bold")
    plt.savefig(png_path, bbox_inches="tight", dpi=150)
    plt.close()

    print(f"Saved: {png_path}")

# Color by L1_annotation
plot_tree(tree, df_meta, "L1_annotation", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")

# Color by c_call (with counts for correct pie charts)
plot_tree(tree, df_meta, "c_call", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")

# Color by sample_clean_fol (with counts for correct pie charts)
plot_tree(tree, df_meta, "sample_clean_fol", sample, custom_plot_path, color_list, counts_dir="../gctree_meta")