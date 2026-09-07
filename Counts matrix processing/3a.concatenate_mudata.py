#!/usr/bin/env python
# %%
import os

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
from muon import MuData

# %%
# --- Parameters ---
input_dir = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/2.qcfiltering/float_output/qc_filtered_data"
output_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/3a.concat-old-forclustering"
output_name = "300k_combined_all_sublibs_newguidelist"

# If True, make barcodes unique by appending the sublib label (e.g. ACGT-1_sublib1).
# Set to False only if you are certain barcodes are already globally unique.
make_barcodes_unique = True

# %%
# --- Discover all .h5mu files in the input directory ---
h5mu_files = sorted([
    f for f in os.listdir(input_dir) if f.endswith(".h5mu")
])
assert len(h5mu_files) > 0, f"No .h5mu files found in {input_dir}"
print(f"Found {len(h5mu_files)} MuData files:")
for f in h5mu_files:
    print(f"  {f}")

# %%
# --- Load each MuData, tag barcodes with the sublib label, collect ---
gex_adatas = []
guide_adatas = []
sublib_labels = []

for fname in h5mu_files:
    sublib = fname.replace("_qc_filtered_gex_and_guide.h5mu", "")
    fpath = os.path.join(input_dir, fname)
    print(f"\nLoading {fname}  (sublib='{sublib}')")

    mdata = mu.read(fpath)
    gex = mdata["GEX"].copy()
    guide = mdata["guide"].copy()

    print(f"  GEX:   {gex.n_obs} cells x {gex.n_vars} genes")
    print(f"  guide: {guide.n_obs} cells x {guide.n_vars} guides")

    # Tag cells with their sublib of origin before any collision check
    gex.obs["sublib"] = sublib
    guide.obs["sublib"] = sublib

    if make_barcodes_unique:
        gex.obs_names = [f"{bc}_{sublib}" for bc in gex.obs_names]
        guide.obs_names = [f"{bc}_{sublib}" for bc in guide.obs_names]

    gex_adatas.append(gex)
    guide_adatas.append(guide)
    sublib_labels.append(sublib)

# %%
# --- Barcode collision check ---
print("\n--- Barcode collision check ---")

all_gex_barcodes = [bc for a in gex_adatas for bc in a.obs_names]
all_guide_barcodes = [bc for a in guide_adatas for bc in a.obs_names]

gex_duplicates = pd.Series(all_gex_barcodes)[pd.Series(all_gex_barcodes).duplicated()].unique()
guide_duplicates = pd.Series(all_guide_barcodes)[pd.Series(all_guide_barcodes).duplicated()].unique()

if len(gex_duplicates) == 0:
    print("GEX:   no duplicate barcodes")
else:
    print(f"GEX:   {len(gex_duplicates)} duplicate barcodes found!")
    print(f"  Examples: {list(gex_duplicates[:10])}")

if len(guide_duplicates) == 0:
    print("guide: no duplicate barcodes")
else:
    print(f"guide: {len(guide_duplicates)} duplicate barcodes found!")
    print(f"  Examples: {list(guide_duplicates[:10])}")

assert len(gex_duplicates) == 0 and len(guide_duplicates) == 0, (
    "Duplicate barcodes detected. Set make_barcodes_unique=True or investigate the files above."
)

# %%
# --- Concatenate GEX ---
print("\n--- Concatenating GEX ---")
gex_combined = ad.concat(
    gex_adatas,
    join="outer",       # keep the union of genes; fills missing with 0
    label="sublib",
    keys=sublib_labels,
    index_unique=None,  # barcodes are already unique
)
print(f"Combined GEX: {gex_combined.n_obs} cells x {gex_combined.n_vars} genes")

# %%
# --- Concatenate guide ---
print("\n--- Concatenating guide ---")
guide_combined = ad.concat(
    guide_adatas,
    join="outer",
    label="sublib",
    keys=sublib_labels,
    index_unique=None,
)
print(f"Combined guide: {guide_combined.n_obs} cells x {guide_combined.n_vars} guides")

# %%
# --- Sanity checks on the combined objects ---
print("\n--- Post-concatenation sanity checks ---")

# Cell counts should match
assert gex_combined.n_obs == guide_combined.n_obs, (
    f"Cell count mismatch after concat: GEX={gex_combined.n_obs}, guide={guide_combined.n_obs}"
)

# Barcodes should be in the same order
assert (gex_combined.obs_names == guide_combined.obs_names).all(), (
    "Barcode order differs between combined GEX and guide objects."
)

# Confirm sublib breakdown
print("Cells per sublib:")
print(gex_combined.obs["sublib"].value_counts().to_string())

print("\nAll checks passed.")

# %%
# --- Build combined MuData and write ---
mdata_combined = MuData({
    "GEX": gex_combined,
    "guide": guide_combined,
})

os.makedirs(output_path, exist_ok=True)
out_file = os.path.join(output_path, f"{output_name}.h5mu")
mdata_combined.write(out_file)
print(f"\nSaved combined MuData to:\n  {out_file}")
print(f"Final: {gex_combined.n_obs} cells x {gex_combined.n_vars} genes")
# %%
