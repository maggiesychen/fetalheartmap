#!/usr/bin/env python
import argparse
import gzip
import os

import muon as mu
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

# %%
# --- Toggle between single sublib testing and batch mode ---
TEST_MODE = False  # set to False for batch submission

# for test mode, just run as python script
if TEST_MODE:
    sublib = "1"  # change this to test different sublibs
else:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sublib", type=str, required=True)
    args = parser.parse_args()
    sublib = args.sublib

input_path = f"/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/2.qcfiltering/float_output/qc_filtered_data/300k_sublib{sublib}_qc_filtered_gex_and_guide.h5mu"
out_dir = f"/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3b.cellranger_format/sublib{sublib}/filtered_feature_bc_matrix/"

os.makedirs(out_dir, exist_ok=True)
print(f"Processing sublib {sublib}")
print(f"Input: {input_path}")
print(f"Output: {out_dir}")

# %%
# load
mdata = mu.read(input_path)
gex = mdata["GEX"].copy()
guide = mdata["guide"].copy()

# converting everything to int here! rounding, not truncating
gex.X = np.round(gex.X).astype(np.int32)
guide.X = np.round(guide.X).astype(np.int32)

print(f"GEX: {gex.n_obs} cells x {gex.n_vars} genes")
print(f"Guide: {guide.n_obs} cells x {guide.n_vars} guides")
print(f"GEX dtype: {gex.X.dtype}")
print(f"Guide dtype: {guide.X.dtype}")


# %%
# sort barcodes consistently, append sublib label
sorted_barcodes = sorted(gex.obs_names)
gex = gex[sorted_barcodes, :]
guide = guide[sorted_barcodes, :]
assert (gex.obs_names == guide.obs_names).all(), "Barcode mismatch after sorting!"
print(f"Barcodes sorted and matched: {len(sorted_barcodes)}")

# verify matrix column order matches barcodes
assert list(gex.obs_names) == sorted_barcodes, "GEX row order doesn't match sorted_barcodes!"
assert list(guide.obs_names) == sorted_barcodes, "Guide row order doesn't match sorted_barcodes!"

# %%
# barcodes.tsv.gz — append sublib to each barcode
print("Barcodes before appending:")
print("Sample barcodes:", sorted_barcodes[:3])

with gzip.open(out_dir + "barcodes.tsv.gz", "wt") as f:
    for bc in sorted_barcodes:
        f.write(f"{bc}_sublib{sublib}\n")
print("Written: barcodes.tsv.gz")

# %%
# features.tsv.gz — genes then guides
# check var columns
print("GEX var columns:", gex.var.columns.tolist())
print("GEX var_names sample:", gex.var_names[:5].tolist())
print("GEX var shape:", gex.var.shape)
print(gex.var.head())

print("Guide var_names sample:", guide.var_names[:5].tolist())
print("Guide var shape:", guide.var.shape)
print(guide.var.head())

# %%
# strip version number from Ensembl IDs e.g. ENSG00000227232.6 -> ENSG00000227232
gene_ids = [g.split(".")[0] for g in gex.var_names.tolist()]

# use symbol where available, fall back to stripped Ensembl ID
# also replace symbols that are just unversioned Ensembl IDs (failed mapping)
if "symbol" in gex.var.columns:
    gene_names = [
        symbol if (isinstance(symbol, str) and not symbol.startswith("ENSG")) else eid
        for symbol, eid in zip(gex.var["symbol"].tolist(), gene_ids)
    ]
else:
    gene_names = gene_ids

gene_features = pd.DataFrame({
    "id": gene_ids,
    "name": gene_names,
    "type": "Gene Expression"
})

# guide features: id, name, type
guide_features = pd.DataFrame({
    "id": guide.var_names.tolist(),
    "name": guide.var_names.tolist(),
    "type": "CRISPR Guide Capture"
})

features = pd.concat([gene_features, guide_features], ignore_index=True)
with gzip.open(out_dir + "features.tsv.gz", "wt") as f:
    features.to_csv(f, sep="\t", header=False, index=False)
print(f"Written: features.tsv.gz ({len(gene_features)} genes + {len(guide_features)} guides)")

# checking features --> lots of printing
print("\n--- Features check ---")
print(f"Total features: {len(features)}")
print(f"Gene features: {len(gene_features)}")
print(f"Guide features: {len(guide_features)}")
print("\nFirst 5 gene features:")
print(features.head())
print("\nLast 5 gene features:")
print(features.iloc[len(gene_features)-5:len(gene_features)])
print("\nFirst 5 guide features:")
print(features.iloc[len(gene_features):len(gene_features)+5])
print("\nLast 5 guide features:")
print(features.tail())
print("\nAny NaN in features:")
print(features.isnull().sum())
print("\nFeature types:")
print(features["type"].value_counts())



# %%
# matrix.mtx.gz — stack genes on top of guides, features x barcodes
gene_matrix = sp.csc_matrix(gex.X).T    # genes x cells
guide_matrix = sp.csc_matrix(guide.X).T  # guides x cells
combined = sp.vstack([gene_matrix, guide_matrix])
print(f"Combined matrix shape: {combined.shape} (features x cells)")

# eliminating 0s to maintain sparse matrix format (added 4/10/26)
combined.eliminate_zeros()  # <-- ADD THIS LINE
print("Eliminated 0 entries!")

with gzip.open(out_dir + "matrix.mtx.gz", "wb") as f:
    sio.mmwrite(f, combined)
print("Written: matrix.mtx.gz")

# %%
# sanity checks
print(f"\n--- Sanity checks ---")
print(f"n_barcodes: {len(sorted_barcodes)}")
print(f"n_features: {len(features)}")
print(f"matrix shape: {combined.shape}")
assert combined.shape[0] == len(features), "Feature count mismatch!"
assert combined.shape[1] == len(sorted_barcodes), "Barcode count mismatch!"

# --- Additional sanity checks ---

# 1. feature count matches matrix rows
assert combined.shape[0] == len(features), f"Feature count mismatch: matrix has {combined.shape[0]} rows but features has {len(features)}"

# 2. barcode count matches matrix columns
assert combined.shape[1] == len(sorted_barcodes), f"Barcode count mismatch: matrix has {combined.shape[1]} cols but {len(sorted_barcodes)} barcodes"

# 3. spot check gene block (first 100 genes, first 100 cells)
gene_block_sample = combined[:100, :100].toarray()
gex_sample = sp.csc_matrix(gex.X).T[:100, :100].toarray()
assert np.allclose(gene_block_sample, gex_sample), "Gene block spot check failed!"

# 4. spot check guide block (first 100 guides, first 100 cells)
guide_block_sample = combined[gex.n_vars:gex.n_vars+100, :100].toarray()
guide_sample = sp.csc_matrix(guide.X).T[:100, :100].toarray()
assert np.allclose(guide_block_sample, guide_sample), "Guide block spot check failed!"

# 5. no negative values
assert combined.data.min() >= 0, "Negative values found in matrix!"

# 6. check features file has right number of each type
n_gene_features = (features["type"] == "Gene Expression").sum()
n_guide_features = (features["type"] == "CRISPR Guide Capture").sum()
assert n_gene_features == gex.n_vars, f"Gene feature count mismatch: {n_gene_features} vs {gex.n_vars}"
assert n_guide_features == guide.n_vars, f"Guide feature count mismatch: {n_guide_features} vs {guide.n_vars}"

print(f"Gene features: {n_gene_features}")
print(f"Guide features: {n_guide_features}")
print(f"Matrix shape: {combined.shape}")
print(f"Matrix dtype: {combined.dtype}")
print(f"Matrix nnz: {combined.nnz}")
print("All sanity checks passed.")
print("All checks passed.")