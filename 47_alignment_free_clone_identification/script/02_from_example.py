#!/usr/bin/env python3

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

import matplotlib
matplotlib.use('Agg')

import sys
sys.path.append("/path/to/Lindenbaum2020/")

# from core import *
from core_updated import *
import pandas as pd
import numpy as np

HH = "HH119"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
repertoire = pd.read_csv(f"../out/{HH}_repertoire.tsv", sep="\t")
table_neg  = pd.read_csv(f"../out/{HH}_negation.tsv",   sep="\t")

print("Data loaded!")

# ------------------------------------------------------------------------------
# Checking data
# ------------------------------------------------------------------------------

print("Checking data...")
print("")
print("Head of repertoire:")
print(repertoire.head())
print("")
print("Head of table_neg:")
print(table_neg.head())

print("")
print("")

df = pd.DataFrame(repertoire.SEQUENCE[:].values)
print("Head of df:")
print(df.head())

print("")

df2=pd.DataFrame(table_neg["SEQUENCE"].values )
print("Head of df2:")
print(df2.head())

# ------------------------------------------------------------------------------
# Run alignment-free clonal clustering
# ------------------------------------------------------------------------------
# W_l: window length - last W_l nucleotides used from each sequence
#      set based on your trimmed sequence length distribution
# per: percentile of negation distance distribution used as threshold
#      0.1 = 10th percentile (from paper default)

print("Running alignement_free_clone...")

result = alignement_free_clone(
    repertoire = repertoire,
    table_neg  = table_neg["SEQUENCE"],  # core.py uses .values on this
    W_l        = 300,                    # adjust based on your length distribution
    per        = 0.1
)

print("alignement_free_clone finished!")

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

print("Saving restults...")

result.to_csv(f"../out/{HH}_clone_output.tsv", sep="\t", index=False)
print(f"Done. {result.CLONE.nunique()} clones identified in {len(result)} sequences.")

