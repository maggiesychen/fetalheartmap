#!/usr/bin/env python
# %%
import argparse
import os

import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
from muon import MuData

# %%
# --- Parameters ---

##doing this per sublib separately so just change the number
#gex_input_path = "/scratch/users/msychen/20260407_mergedfastq_barcoderanksinflection_300kround3feb_ipscvicd9_newguidelist_analysis/merged:iPSC-VICs_sublibrary_1_merged/kb_all_main_raw/counts_unfiltered_modified/adata.h5ad"       # path to GEX AnnData (.h5ad or .h5mu modality)
#guide_input_path = "/scratch/users/msychen/20260407_mergedfastq_barcoderanksinflection_300kround3feb_ipscvicd9_newguidelist_analysis/merged:splitseq_guide_pool_sublibrary_1_merged/kb_guide_main_raw/counts_unfiltered_modified/adata.h5ad"     # path to guide/gRNA AnnData (.h5ad or .h5mu modality)
#output_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing"          # directory where output will be written
#output_name = "300k_sublib1_umi_filtered_gex_and_guide"  # base name for output .h5mu file
#sublib = "sublib1"  # label for this run; plots saved under plots/QC/<sublib>/

## run this part for batch submission
# --- Parameters ---
parser = argparse.ArgumentParser()
parser.add_argument("--sublib", type=str, required=True)
parser.add_argument("--umi_threshold", type=int, default=1921)
parser.add_argument("--umi_max", type=int, default=100000)
args = parser.parse_args()

sublib = args.sublib
umi_threshold = args.umi_threshold
umi_max = args.umi_max

gex_input_path = f"/scratch/users/msychen/20260407_mergedfastq_barcoderanksinflection_300kround3feb_ipscvicd9_newguidelist_analysis/merged:iPSC-VICs_sublibrary_{sublib}_merged/kb_all_main_raw/counts_unfiltered_modified/adata.h5ad"
guide_input_path = f"/scratch/users/msychen/20260407_mergedfastq_barcoderanksinflection_300kround3feb_ipscvicd9_newguidelist_analysis/merged:splitseq_guide_pool_sublibrary_{sublib}_merged/kb_guide_main_raw/counts_unfiltered_modified/adata.h5ad"
output_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/1.umifiltering/output/20260409-intcorrected-output"
output_name = f"300k_sublib{sublib}_umi_filtered_gex_and_guide"

print(f"Running sublib {sublib} | umi_threshold={umi_threshold} | umi_max={umi_max}")

# %%
os.makedirs(output_path, exist_ok=True)

# %%
# Load GEX and guide AnnData objects
gex_adata = sc.read(gex_input_path)
guide_adata = sc.read(guide_input_path)

# print to check values before any processing
print(gex_adata.X[:5, :5].toarray())
print(guide_adata.X[:5, :5].toarray())

# %%
# Compute per-cell UMI counts and sort cells descending by UMI count
umi_counts = np.array(gex_adata.X.sum(axis=1, dtype=np.float32)).flatten()
gex_adata.obs["total_counts"] = umi_counts

knee_df = pd.DataFrame(
    {
        "total_counts": umi_counts,
        "barcode": gex_adata.obs_names.values,
    }
)
knee_df = knee_df.sort_values("total_counts", ascending=False).reset_index(drop=True)
knee_df["rank"] = knee_df.index + 1

print(f"Total cells before UMI filtering: {len(knee_df)}")
print(
    f"Cells with total_counts >= {umi_threshold}: "
    f"{(knee_df['total_counts'] >= umi_threshold).sum()}"
)

# %%
# Plot UMI rank curve so you can visually confirm the threshold
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(knee_df["rank"], knee_df["total_counts"], linewidth=1, color="steelblue")
ax.axhline(y=umi_threshold, color="red", linestyle="--", label=f"min threshold = {umi_threshold}")
ax.axhline(y=umi_max, color="darkorange", linestyle="--", label=f"max threshold = {umi_max}")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Cell rank (by UMI count)", fontsize=12)
ax.set_ylabel("UMI count", fontsize=12)
ax.set_title("UMI rank curve", fontsize=13)
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_path, f"{output_name}_umi_rank_curve.png"), dpi=150)
plt.close()

# %%
# Subset GEX to cells within the UMI range [umi_threshold, umi_max]
gex_adata = gex_adata[
    (gex_adata.obs["total_counts"] >= umi_threshold)
    & (gex_adata.obs["total_counts"] <= umi_max),
    :,
].copy()
print(f"Cells after UMI filtering (>={umi_threshold} and <={umi_max}): {gex_adata.n_obs}")

# %%
# Get shared barcodes (sorted to ensure consistent order)
shared_barcodes = sorted(
    set(gex_adata.obs_names).intersection(guide_adata.obs_names)
)
print(f"Barcodes present in both GEX and guide after filtering: {len(shared_barcodes)}")

# %%
# subset BOTH to shared barcodes
gex_adata = gex_adata[shared_barcodes, :].copy()
guide_adata = guide_adata[shared_barcodes, :].copy()
print("Both GEX and guide subsetted to shared barcodes")

#%%
# Check layers BEFORE rounding
if all(k in gex_adata.layers for k in ["nascent", "ambiguous", "mature"]):
    layer_sum = gex_adata.layers["nascent"] + gex_adata.layers["ambiguous"] + gex_adata.layers["mature"]
    diff = gex_adata.X - layer_sum  # both still floats here
    print("Max difference between X and layer sum:", np.abs(diff).max())
    print("X equals layer sum:", np.allclose(gex_adata.X.toarray(), layer_sum.toarray()))

# %%
# after subsetting GEX and guide, THEN round everything to int
gex_adata.X = np.round(gex_adata.X).astype(np.int32)
guide_adata.X = np.round(guide_adata.X).astype(np.int32)
print("Counts converted to int32")

# print to check values and make sure they're int
print(gex_adata.X[:5, :5].toarray())
print(guide_adata.X[:5, :5].toarray())

# %%
# Check if layers exist
print("Layers available:", list(gex_adata.layers.keys()))

# now look at layer sum again
if all(k in gex_adata.layers for k in ["nascent", "ambiguous", "mature"]):
    layer_sum = gex_adata.layers["nascent"] + gex_adata.layers["ambiguous"] + gex_adata.layers["mature"]
    diff = gex_adata.X - layer_sum
    print("Max difference between X and layer sum:", np.abs(diff).max())
    print("X equals layer sum:", np.allclose(gex_adata.X.toarray(), layer_sum.toarray()))
else:
    print("No nac layers found — X is likely the direct count matrix")

# %%
# Final sanity checks before writing
print(f"Final GEX: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")
print(f"Final guide: {guide_adata.n_obs} cells x {guide_adata.n_vars} guides")
print(f"GEX dtype: {gex_adata.X.dtype}")
print(f"Guide dtype: {guide_adata.X.dtype}")
print(f"GEX n_obs: {gex_adata.n_obs}")
print(f"Guide n_obs: {guide_adata.n_obs}")
if gex_adata.n_obs == guide_adata.n_obs:
    print(f"GEX barcodes match guide barcodes: {(gex_adata.obs_names == guide_adata.obs_names).all()}")
else:
    print(f"WARNING: GEX and guide have different cell counts — subsetting GEX to shared_barcodes")
print("GEX sample values:")
print(gex_adata.X[:5, :5].toarray())
print("Guide sample values:")
print(guide_adata.X[:5, :5].toarray())

# %%
# Build MuData and write
mudata_dict = {
    "GEX": gex_adata,    # already subsetted and rounded
    "guide": guide_adata, # already subsetted and rounded
}
mdata = MuData(mudata_dict)
out_file = os.path.join(output_path, f"{output_name}.h5mu")
mdata.write(out_file)
print(f"Saved MuData to: {out_file}")