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
output_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/1.umifiltering/output"
output_name = f"300k_sublib{sublib}_umi_filtered_gex_and_guide"

print(f"Running sublib {sublib} | umi_threshold={umi_threshold} | umi_max={umi_max}")

# %%
os.makedirs(output_path, exist_ok=True)

# %%
# Load GEX and guide AnnData objects
gex_adata = sc.read(gex_input_path)
guide_adata = sc.read(guide_input_path)

gex_adata.X = gex_adata.X.astype(np.float32)

# %%
# Compute per-cell UMI counts and sort cells descending by UMI count
umi_counts = np.array(gex_adata.X.sum(axis=1)).flatten()
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
# Plot UMI rank curve
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
# Subset guide AnnData to the same set of barcodes that remain in GEX
shared_barcodes = list(
    set(gex_adata.obs_names).intersection(guide_adata.obs_names)
)
print(
    f"Barcodes present in both GEX and guide after filtering: {len(shared_barcodes)}"
)

# subset guide and GEX to correct barcodes 
mudata_dict = {
    "GEX": gex_adata[shared_barcodes, :].copy(),
    "guide": guide_adata[shared_barcodes, :].copy(),
}
mdata = MuData(mudata_dict)

out_file = os.path.join(output_path, f"{output_name}.h5mu")
mdata.write(out_file)
print(f"Saved MuData to: {out_file}")