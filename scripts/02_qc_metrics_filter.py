#!/usr/bin/env python
"""Step 2 -- per-cell and per-gene QC filtering for one sublibrary.

Applies, in order:

1. mitochondrial and ribosomal percentage ceilings;
2. MAD-based outlier removal on ``log1p_total_counts`` and
   ``log1p_n_genes_by_counts``;
3. optional per-cell minimum gene / count filters;
4. per-gene minimum cell / count filters.

MADs are computed on the unfiltered population and applied after the MT/ribo
ceilings, so the cutoffs reflect the pre-filter distribution. Gene filters run
per sublibrary, which is why sublibraries end up with different gene sets --
see ``convert_to_cellranger.py --master-gene-list``.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import muon as mu  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from muon import MuData  # noqa: E402
from scipy.stats import median_abs_deviation  # noqa: E402

from scripts.pipeline_utils import annotate_mito_ribo, load_config

METRIC_LABELS = {
    "n_genes_by_counts": "# genes per cell",
    "total_counts": "# UMIs per cell",
    "pct_counts_mt": "% mitochondrial counts",
    "pct_counts_ribo": "% ribosomal counts",
}
QC_METRICS = list(METRIC_LABELS)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="UMI-filtered .h5mu")
    parser.add_argument("--output", required=True, help="Output .h5mu path")
    parser.add_argument("--config", required=True, help="Workflow config YAML")
    parser.add_argument(
        "--plots-dir", default=None, help="Directory for the QC plots (optional)"
    )
    return parser.parse_args()


def is_outlier(adata, metric, n_mads):
    """Boolean series: True where ``metric`` is > ``n_mads`` MADs from median."""
    vals = adata.obs[metric].values.astype(float)
    med = np.median(vals)
    mad = median_abs_deviation(vals)
    return pd.Series(
        (vals < med - n_mads * mad) | (vals > med + n_mads * mad),
        index=adata.obs_names,
    )


def plot_qc_violins(adata, metrics, thresholds, title_suffix, save_path):
    n = len(metrics)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, metric in zip(axes, metrics):
        vals = adata.obs[metric].dropna()
        ax.violinplot(vals, positions=[0], showmedians=True)
        ax.set_xticks([])
        ax.set_ylabel(METRIC_LABELS.get(metric, metric), fontsize=11)
        if metric in thresholds:
            ax.axhline(
                thresholds[metric],
                color="red",
                linestyle="--",
                linewidth=1,
                label=f"threshold={thresholds[metric]}",
            )
            ax.legend(fontsize=8)
    fig.suptitle(f"QC metrics {title_suffix}", fontsize=13)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


def plot_umi_gene_scatter(adata, title_suffix, save_path):
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        s=3,
        alpha=0.3,
        linewidth=0,
        rasterized=True,
        color="steelblue",
    )
    ax.set_xlabel("# UMIs per cell", fontsize=12)
    ax.set_ylabel("# genes per cell", fontsize=12)
    ax.set_title(f"UMI vs genes {title_suffix}", fontsize=13)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {save_path}")


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["qc_filtering"]

    plots_dir = args.plots_dir
    if plots_dir:
        os.makedirs(plots_dir, exist_ok=True)

    print(f"Loading MuData from:\n  {args.input}")
    mdata = mu.read(args.input)
    gex_adata = mdata["GEX"].copy()
    guide_adata = mdata["guide"].copy()
    print(f"Loaded: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")

    annotate_mito_ribo(
        gex_adata,
        reference=params["reference"],
        ensembl_to_symbol=config["input_paths"].get("ensembl_to_symbol"),
    )

    sc.pp.calculate_qc_metrics(
        gex_adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
    )

    qc_thresholds = {
        "pct_counts_mt": params["mt_threshold"],
        "pct_counts_ribo": params["ribo_threshold"],
    }

    if plots_dir:
        plot_qc_violins(
            gex_adata,
            metrics=QC_METRICS,
            thresholds=qc_thresholds,
            title_suffix="(before filtering)",
            save_path=os.path.join(plots_dir, "qc_violin_before_filtering.png"),
        )
        plot_umi_gene_scatter(
            gex_adata,
            title_suffix="(before filtering)",
            save_path=os.path.join(plots_dir, "umi_gene_scatter_before_filtering.png"),
        )

    # MADs are computed here, on the unfiltered population, and applied below.
    print(f"\nBefore filtering: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")
    gex_adata.obs["outlier"] = is_outlier(
        gex_adata, "log1p_total_counts", params["n_mads_total_counts"]
    ) | is_outlier(
        gex_adata, "log1p_n_genes_by_counts", params["n_mads_n_genes"]
    )
    print(f"Outlier cells flagged: {gex_adata.obs['outlier'].sum()}")

    gex_adata = gex_adata[gex_adata.obs["pct_counts_mt"] <= params["mt_threshold"], :]
    print(f"After MT filter   (<={params['mt_threshold']}%): {gex_adata.n_obs} cells")
    gex_adata = gex_adata[
        gex_adata.obs["pct_counts_ribo"] <= params["ribo_threshold"], :
    ]
    print(f"After ribo filter (<={params['ribo_threshold']}%): {gex_adata.n_obs} cells")

    if params["filter_outliers"]:
        gex_adata = gex_adata[~gex_adata.obs["outlier"]].copy()
        print(f"After outlier removal:            {gex_adata.n_obs} cells")
    else:
        print("Skipping outlier removal (filter_outliers=false)")

    if params["filter_cells_by_min_genes"]:
        sc.pp.filter_cells(gex_adata, min_genes=params["min_genes_per_cell"])
        sc.pp.filter_cells(gex_adata, min_counts=params["min_counts_per_cell"])
        print(f"After min-gene/count filter:      {gex_adata.n_obs} cells")

    print(
        f"\nApplying gene filters: min_cells={params['min_cells_per_gene']}, "
        f"min_counts={params['min_counts_per_gene']}"
    )
    sc.pp.filter_genes(gex_adata, min_cells=params["min_cells_per_gene"])
    sc.pp.filter_genes(gex_adata, min_counts=params["min_counts_per_gene"])
    print(f"After gene filtering: {gex_adata.n_obs} cells x {gex_adata.n_vars} genes")

    sc.pp.calculate_qc_metrics(
        gex_adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
    )
    if plots_dir:
        plot_qc_violins(
            gex_adata,
            metrics=QC_METRICS,
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

    shared_barcodes = sorted(
        set(gex_adata.obs_names).intersection(guide_adata.obs_names)
    )
    print(f"Barcodes in both GEX and guide: {len(shared_barcodes)}")

    mdata_filtered = MuData(
        {
            "GEX": gex_adata[shared_barcodes, :].copy(),
            "guide": guide_adata[shared_barcodes, :].copy(),
        }
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    mdata_filtered.write(args.output)
    print(f"\nSaved filtered MuData to:\n  {args.output}")
    print(f"GEX dtype:   {mdata_filtered['GEX'].X.dtype}")
    print(f"guide dtype: {mdata_filtered['guide'].X.dtype}")


if __name__ == "__main__":
    main()
