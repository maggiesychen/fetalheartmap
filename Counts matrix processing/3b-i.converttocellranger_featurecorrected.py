#!/usr/bin/env python
import gzip
import os

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

# %%
# --- Toggle between single sublib testing and batch mode ---
TEST_MODE = False  # set to True to test a single sublib

if TEST_MODE:
    sublib = "1"
else:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sublib", type=str, required=True)
    args = parser.parse_args()
    sublib = args.sublib

base_input = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/2.qcfiltering/float_output/qc_filtered_data"
base_output = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3b-i.cellranger_format_featurecorrected"
master_gene_list_path = f"{base_output}/master_gene_ids.txt"

input_path = f"{base_input}/300k_sublib{sublib}_qc_filtered_gex_and_guide.h5mu"
out_dir = f"{base_output}/sublib{sublib}/filtered_feature_bc_matrix/"

os.makedirs(out_dir, exist_ok=True)
print(f"Processing sublib {sublib}")
print(f"Input: {input_path}")
print(f"Output: {out_dir}")

# %%
# Pass 1: build master gene list across all sublibs (only if not already saved)
# This is saved to disk so it persists across separate Slurm jobs
if not os.path.exists(master_gene_list_path):
    print("Building master gene list from all sublibs...")
    all_gene_ids = set()
    for s in ["1", "2", "3", "4", "5", "6", "7"]:
        path_s = f"{base_input}/300k_sublib{s}_qc_filtered_gex_and_guide.h5mu"
        mdata_s = mu.read(path_s)
        gene_ids_s = [g.split(".")[0] for g in mdata_s["GEX"].var_names.tolist()]
        all_gene_ids.update(gene_ids_s)
        print(f"  sublib{s}: {len(gene_ids_s)} genes")
    master_gene_ids = sorted(list(all_gene_ids))
    # save to disk
    with open(master_gene_list_path, "w") as f:
        for g in master_gene_ids:
            f.write(g + "\n")
    print(f"Master gene list saved: {len(master_gene_ids)} genes → {master_gene_list_path}")
else:
    with open(master_gene_list_path, "r") as f:
        master_gene_ids = [line.strip() for line in f.readlines()]
    print(f"Master gene list loaded: {len(master_gene_ids)} genes from {master_gene_list_path}")

master_gene_ids_set = set(master_gene_ids)

# %%
# load sublib
mdata = mu.read(input_path)
gex = mdata["GEX"].copy()
guide = mdata["guide"].copy()

# strip version numbers from Ensembl IDs BEFORE reindexing
gex.var_names = pd.Index([g.split(".")[0] for g in gex.var_names.tolist()])

# converting everything to int here! rounding, not truncating
gex.X = np.round(gex.X).astype(np.int32)
guide.X = np.round(guide.X).astype(np.int32)

print(f"GEX before reindex: {gex.n_obs} cells x {gex.n_vars} genes")
print(f"Guide: {guide.n_obs} cells x {guide.n_vars} guides")

# %%
# reindex GEX to master gene list — adds missing genes as zeros
missing_genes = [g for g in master_gene_ids if g not in gex.var_names]
print(f"Missing genes in this sublib: {len(missing_genes)}")

if missing_genes:
    empty = ad.AnnData(
        X=sp.csc_matrix((gex.n_obs, len(missing_genes))),
        obs=gex.obs,
        var=pd.DataFrame(index=missing_genes)
    )
    gex = ad.concat([gex, empty], axis=1)

# reorder to match master list exactly
gex = gex[:, master_gene_ids]
print(f"GEX after reindex: {gex.n_obs} cells x {gex.n_vars} genes")
assert gex.n_vars == len(master_gene_ids), "Reindex failed — gene count mismatch!"
assert list(gex.var_names) == master_gene_ids, "Reindex failed — gene order mismatch!"

# %%
# sort barcodes consistently, append sublib label
sorted_barcodes = sorted(gex.obs_names)
gex = gex[sorted_barcodes, :]
guide = guide[sorted_barcodes, :]
assert (gex.obs_names == guide.obs_names).all(), "Barcode mismatch after sorting!"
print(f"Barcodes sorted and matched: {len(sorted_barcodes)}")

# %%
# barcodes.tsv.gz
with gzip.open(out_dir + "barcodes.tsv.gz", "wt") as f:
    for bc in sorted_barcodes:
        f.write(f"{bc}_sublib{sublib}\n")
print("Written: barcodes.tsv.gz")

# %%
# features.tsv.gz — genes (master list) then guides
# gene names: use symbol where available, fall back to Ensembl ID
if "symbol" in gex.var.columns:
    gene_names = [
        symbol if (isinstance(symbol, str) and not symbol.startswith("ENSG")) else eid
        for symbol, eid in zip(gex.var["symbol"].fillna("").tolist(), master_gene_ids)
    ]
else:
    gene_names = master_gene_ids

gene_features = pd.DataFrame({
    "id": master_gene_ids,
    "name": gene_names,
    "type": "Gene Expression"
})

guide_features = pd.DataFrame({
    "id": guide.var_names.tolist(),
    "name": guide.var_names.tolist(),
    "type": "CRISPR Guide Capture"
})

features = pd.concat([gene_features, guide_features], ignore_index=True)

# check for duplicate feature IDs
assert features["id"].is_unique, "Duplicate feature IDs found!"

with gzip.open(out_dir + "features.tsv.gz", "wt") as f:
    features.to_csv(f, sep="\t", header=False, index=False)
print(f"Written: features.tsv.gz ({len(gene_features)} genes + {len(guide_features)} guides)")

# save reference features for cross-sublib check
reference_features_path = f"{base_output}/reference_features.txt"
if not os.path.exists(reference_features_path):
    with open(reference_features_path, "w") as f:
        for fid in features["id"].tolist():
            f.write(fid + "\n")
    print(f"Reference features saved to {reference_features_path}")
else:
    with open(reference_features_path, "r") as f:
        reference_features = [line.strip() for line in f.readlines()]
    assert features["id"].tolist() == reference_features, f"Feature mismatch in sublib {sublib}!"
    print(f"sublib{sublib}: features match reference ✓")

# %%
# matrix.mtx.gz — stack genes on top of guides, features x barcodes
gene_matrix = sp.csc_matrix(gex.X).T    # genes x cells
guide_matrix = sp.csc_matrix(guide.X).T  # guides x cells
combined = sp.vstack([gene_matrix, guide_matrix])
print(f"Combined matrix shape: {combined.shape} (features x cells)")

# eliminate explicit zeros (added 4/10/26)
combined.eliminate_zeros()
print(f"After eliminate_zeros — nnz: {combined.nnz}")

with gzip.open(out_dir + "matrix.mtx.gz", "wb") as f:
    sio.mmwrite(f, combined)
print("Written: matrix.mtx.gz")

# %%
# matrix stats check — counts should differ per sublib
print(f"\nsublib{sublib} matrix stats:")
print(f"  shape: {combined.shape}")
print(f"  nnz: {combined.nnz}")
print(f"  total UMI count: {combined.sum()}")
print(f"  mean UMIs per cell: {combined.sum() / combined.shape[1]:.1f}")

# spot check ACTB
if "ENSG00000075624" in master_gene_ids:
    gene_idx = master_gene_ids.index("ENSG00000075624")
    gene_counts = combined[gene_idx, :].toarray().flatten()
    print(f"  ACTB (ENSG00000075624): mean={gene_counts.mean():.2f}, max={gene_counts.max()}, nonzero={np.count_nonzero(gene_counts)}")

# %%
# sanity checks
print(f"\n--- Sanity checks ---")
assert combined.shape[0] == len(features), f"Feature count mismatch: matrix {combined.shape[0]} vs features {len(features)}"
assert combined.shape[1] == len(sorted_barcodes), f"Barcode count mismatch: matrix {combined.shape[1]} vs {len(sorted_barcodes)}"

# spot check a gene known to exist in this sublib (not a missing/zero gene)
existing_genes = [g for g in master_gene_ids if g not in missing_genes]
if len(existing_genes) >= 1:
    test_idx = master_gene_ids.index(existing_genes[0])
    assert np.allclose(
        combined[test_idx, :100].toarray(),
        sp.csc_matrix(gex.X).T[test_idx, :100].toarray()
    ), "Existing gene spot check failed!"
    print(f"Existing gene spot check passed for {existing_genes[0]} ✓")

# spot check guide block
guide_block_sample = combined[len(master_gene_ids):len(master_gene_ids)+100, :100].toarray()
guide_sample = sp.csc_matrix(guide.X).T[:100, :100].toarray()
assert np.allclose(guide_block_sample, guide_sample), "Guide block spot check failed!"

# no negative values
if combined.data.size > 0:
    assert combined.data.min() >= 0, "Negative values found in matrix!"

# feature type counts
n_gene_features = (features["type"] == "Gene Expression").sum()
n_guide_features = (features["type"] == "CRISPR Guide Capture").sum()
assert n_gene_features == len(master_gene_ids), f"Gene feature count mismatch: {n_gene_features} vs {len(master_gene_ids)}"
assert n_guide_features == guide.n_vars, f"Guide feature count mismatch: {n_guide_features} vs {guide.n_vars}"

print(f"Gene features: {n_gene_features}")
print(f"Guide features: {n_guide_features}")
print(f"Matrix shape: {combined.shape}")
print(f"Matrix dtype: {combined.dtype}")
print(f"Matrix nnz: {combined.nnz}")
print("All sanity checks passed.")
print("All checks passed.")