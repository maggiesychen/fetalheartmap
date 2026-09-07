#!/usr/bin/env python
# %%
import os
import argparse
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
from muon import MuData
from scipy.stats import median_abs_deviation

# %%
##### uncomment if running one sublib at a time
# --- Parameters ---
#input_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260223-WINTERVICHIGHMOISCREEN-round1-novaseq/countsmatrices/postprocessing/20260317-postprocessing-newguidelist/sublib4_umi_filtered_gex_and_guide.h5mu"
#output_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260223-WINTERVICHIGHMOISCREEN-round1-novaseq/countsmatrices/postprocessing/20260317-postprocessing-newguidelist"
#output_name = "sublib4_qc_filtered_gex_and_guide"
#sublib = "sublib4"  # label for this run; plots saved under plots/QC/<sublib>/

parser = argparse.ArgumentParser()
parser.add_argument("--sublib", type=str, required=True)
args = parser.parse_args()

sublib = args.sublib

# make sure input is from umi filtered h5mu objs
input_path = f"/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/1.umifiltering/output/300k_sublib{sublib}_umi_filtered_gex_and_guide.h5mu"
output_path = "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/2.qcfiltering/float_output"
output_name = f"300k_sublib{sublib}_qc_filtered_gex_and_guide"
sublib_label = f"sublib{sublib}"

print(f"Running QC filtering for sublib {sublib}")



reference = "human"  # "human" or "mouse"

# Path to ensembl→symbol dict (used if var_names are ENSG IDs)
ens2symbol_path = "/oak/stanford/groups/engreitz/Users/opushkar/genome/ensembl_to_symbol_dict_v43.npy"

# QC filtering thresholds
mt_threshold = 15       # max % mitochondrial UMIs per cell
ribo_threshold = 4    # max % ribosomal UMIs per cell --> Olga's script has 5, we used 12 for the first iPSC-VIC screen (winter 2025)
n_mads_total_counts = 4 # MADs from median for log1p_total_counts outlier cutoff
n_mads_n_genes = 4      # MADs from median for log1p_n_genes_by_counts outlier cutoff

filter_outliers = True
filter_cells_by_min_genes = False
min_genes_per_cell = 10
min_counts_per_cell = 500

# per gene filtering thresholds
min_cells_per_gene = 10
min_counts_per_gene = 100 ## changed this to 100 to make less conservative (20260304)

# %%
# --- Output directories ---
plots_dir = os.path.join(output_path, "plots", "QC", sublib_label)
filtered_data_dir = os.path.join(output_path, "qc_filtered_data")
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(filtered_data_dir, exist_ok=True)

# %%
# --- Helper: MAD-based outlier detection ---
def is_outlier(adata, metric, n_mads):
    """Return boolean array: True if cell is an outlier on `metric`."""
    vals = adata.obs[metric].values.astype(float)
    med = np.median(vals)
    mad = median_abs_deviation(vals)
    return pd.Series(
        (vals < med - n_mads * mad) | (vals > med + n_mads * mad),
        index=adata.obs_names,
    )


# %%
# --- Helper: QC violin plots ---
def plot_qc_violins(adata, metrics, thresholds, title_suffix, save_path):
    n = len(metrics)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 5))
    if n == 1:
        axes = [axes]
    labels = {
        "n_genes_by_counts": "# genes per cell",
        "total_counts": "# UMIs per cell",
        "pct_counts_mt": "% mitochondrial counts",
        "pct_counts_ribo": "% ribosomal counts",
    }
    for ax, metric in zip(axes, metrics):
        vals = adata.obs[metric].dropna()
        ax.violinplot(vals, positions=[0], showmedians=True)
        ax.set_xticks([])
        ax.set_ylabel(labels.get(metric, metric), fontsize=11)
        if metric in thresholds:
            ax.axhline(thresholds[metric], color="red", linestyle="--",
                       linewidth=1, label=f"threshold={thresholds[metric]}")
            ax.legend(fontsize=8)
    fig.suptitle(f"QC metrics {title_suffix}", fontsize=13)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# %%
# --- Helper: UMI vs genes scatter ---
def plot_umi_gene_scatter(adata, title_suffix, save_path, color_col=None):
    fig, ax = plt.subplots(figsize=(6, 5))
    x = adata.obs["total_counts"]
    y = adata.obs["n_genes_by_counts"]
    if color_col and color_col in adata.obs.columns:
        sc_obj = ax.scatter(x, y, c=adata.obs[color_col], s=3, alpha=0.4,
                            linewidth=0, rasterized=True, cmap="tab20")
        plt.colorbar(sc_obj, ax=ax, label=color_col)
    else:
        ax.scatter(x, y, s=3, alpha=0.3, linewidth=0,
                   rasterized=True, color="steelblue")
    ax.set_xlabel("# UMIs per cell", fontsize=12)
    ax.set_ylabel("# genes per cell", fontsize=12)
    ax.set_title(f"UMI vs genes {title_suffix}", fontsize=13)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


# %%
# --- Load MuData ---
print(f"Loading MuData from:\n  {input_path}")
mdata = mu.read(input_path)
gex_adata = mdata["GEX"].copy()
guide_adata = mdata["guide"].copy()

# changing everything to int
#print("GEX dtype before:", gex_adata.X.dtype)
#print(gex_adata.X[:5, :5].toarray())
#gex_adata.X = np.round(gex_adata.X).astype(np.int32)
#print("GEX dtype after:", gex_adata.X.dtype)
#print(gex_adata.X[:5, :5].toarray())

# changing everything to int
#print("guide dtype before:", guide_adata.X.dtype)
#print(guide_adata.X[:5, :5].toarray())
#guide_adata.X = np.round(guide_adata.X).astype(np.int32)
#print("guide dtype after:", guide_adata.X.dtype)
#print(guide_adata.X[:5, :5].toarray())


print(f"Loaded: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")

# %%
# --- Annotate MT and ribo genes ---
if reference == "human":
    mt_prefix = "MT-"
    ribo_prefix = ("RPS", "RPL")
else:
    mt_prefix = "Mt-"
    ribo_prefix = ("Rps", "Rpl")

# Resolve gene symbols from var_names (handle ensembl IDs or plain symbols)
var_names_are_ensembl = gex_adata.var_names[0].startswith("ENSG") or \
                        gex_adata.var_names[0].startswith("ENSMUSG")

if var_names_are_ensembl:
    # Try to load ensembl→symbol mapping
    if os.path.exists(ens2symbol_path):
        ens2symbol = np.load(ens2symbol_path, allow_pickle=True).item()
        symbols = (
            pd.Series(gex_adata.var_names)
            .str.split(".").str[0]
            .map(ens2symbol)
        )
    else:
        print(f"Warning: ens2symbol dict not found at {ens2symbol_path}. "
              "MT/ribo annotation may be incomplete.")
        symbols = pd.Series(gex_adata.var_names)
else:
    # var_names are already gene symbols
    if "gene_id" in gex_adata.var.columns:
        gex_adata.var.reset_index(inplace=True)
        gex_adata.var.index = gex_adata.var["gene_id"]
    symbols = pd.Series(gex_adata.var_names, index=gex_adata.var_names)

gex_adata.var["symbol"] = symbols.values
gex_adata.var["mt"] = gex_adata.var["symbol"].str.startswith(mt_prefix).fillna(False)
gex_adata.var["ribo"] = gex_adata.var["symbol"].str.startswith(ribo_prefix).fillna(False)
print(f"MT genes: {gex_adata.var['mt'].sum()}  |  Ribo genes: {gex_adata.var['ribo'].sum()}")

# %%
# --- Calculate QC metrics ---
sc.pp.calculate_qc_metrics(
    gex_adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
)

# %%
# --- Plots before filtering ---
qc_thresholds = {
    "pct_counts_mt": mt_threshold,
    "pct_counts_ribo": ribo_threshold,
}
plot_qc_violins(
    gex_adata,
    metrics=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    thresholds=qc_thresholds,
    title_suffix="(before filtering)",
    save_path=os.path.join(plots_dir, "qc_violin_before_filtering.png"),
)
plot_umi_gene_scatter(
    gex_adata,
    title_suffix="(before filtering)",
    save_path=os.path.join(plots_dir, "umi_gene_scatter_before_filtering.png"),
)

# %%
# --- Outlier detection (MAD-based) ---
print(f"\nBefore filtering: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")
gex_adata.obs["outlier"] = (
    is_outlier(gex_adata, "log1p_total_counts", n_mads_total_counts)
    | is_outlier(gex_adata, "log1p_n_genes_by_counts", n_mads_n_genes)
)
print(f"Outlier cells flagged: {gex_adata.obs['outlier'].sum()}")

# %%
# --- Filter by MT and ribo thresholds ---
gex_adata = gex_adata[gex_adata.obs["pct_counts_mt"] <= mt_threshold, :]
print(f"After MT filter  (≤{mt_threshold}%):   {gex_adata.n_obs} cells")
gex_adata = gex_adata[gex_adata.obs["pct_counts_ribo"] <= ribo_threshold, :]
print(f"After ribo filter (≤{ribo_threshold}%):  {gex_adata.n_obs} cells")

# %%
# --- Filter MAD outliers ---
if filter_outliers:
    gex_adata = gex_adata[~gex_adata.obs["outlier"]].copy()
    print(f"After outlier removal:            {gex_adata.n_obs} cells")
else:
    print("Skipping outlier removal (filter_outliers=False)")

# %%
# --- Optional per-cell min-gene/count filter ---
if filter_cells_by_min_genes:
    sc.pp.filter_cells(gex_adata, min_genes=min_genes_per_cell)
    sc.pp.filter_cells(gex_adata, min_counts=min_counts_per_cell)
    print(f"After min-gene/count filter:      {gex_adata.n_obs} cells")

# %%
# --- Gene filters ---
print(f"\nApplying gene filters: min_cells={min_cells_per_gene}, min_counts={min_counts_per_gene}")
sc.pp.filter_genes(gex_adata, min_cells=min_cells_per_gene)
sc.pp.filter_genes(gex_adata, min_counts=min_counts_per_gene)
print(f"After gene filtering: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")

# %%
# --- Plots after filtering ---
sc.pp.calculate_qc_metrics(
    gex_adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
)
plot_qc_violins(
    gex_adata,
    metrics=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    thresholds=qc_thresholds,
    title_suffix="(after filtering)",
    save_path=os.path.join(plots_dir, "qc_violin_after_filtering.png"),
)
plot_umi_gene_scatter(
    gex_adata,
    title_suffix="(after filtering)",
    save_path=os.path.join(plots_dir, "umi_gene_scatter_after_filtering.png"),
)
print(f"\nFinal dataset: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")

# %%
# --- Subset guide to matching barcodes and write MuData ---
shared_barcodes = list(
    set(gex_adata.obs_names).intersection(guide_adata.obs_names)
)
print(f"Barcodes in both GEX and guide: {len(shared_barcodes)}")

mdata_filtered = MuData({
    "GEX": gex_adata[shared_barcodes, :].copy(),
    "guide": guide_adata[shared_barcodes, :].copy(),
})

out_file = os.path.join(filtered_data_dir, f"{output_name}.h5mu")
mdata_filtered.write(out_file)
print(f"\nSaved filtered MuData to:\n  {out_file}")

#checking counts type 
print(mdata_filtered["GEX"].X.dtype)
print(mdata_filtered["guide"].X.dtype)
# %%
