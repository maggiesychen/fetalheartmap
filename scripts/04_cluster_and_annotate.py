#!/usr/bin/env python
"""Step 4 -- normalise, cluster and annotate the combined GEX matrix.

Follows the scanpy clustering tutorial
(https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html):
CP10k + log1p, cell-cycle scoring, HVG selection, PCA, a k-NN graph, UMAP, and
Leiden at several resolutions, then marker-panel scoring and dotplots.

The raw counts are preserved in ``layers['counts']`` and HVG selection reads
that layer, because ``flavor='seurat_v3'`` expects counts rather than
log-normalised values.

Cell-cycle regression is **off** by default: the published run did not regress,
and its output is what every downstream step consumes. See the README.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import muon as mu  # noqa: E402
import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

from scripts.pipeline_utils import load_config, load_symbol_map, read_gene_list

matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False
matplotlib.rcParams["svg.fonttype"] = "none"
sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=False)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Combined .h5mu")
    parser.add_argument("--output", required=True, help="Output clustered .h5ad")
    parser.add_argument("--config", required=True, help="Workflow config YAML")
    parser.add_argument(
        "--plots-dir", default=None, help="Directory for the clustering plots"
    )
    return parser.parse_args()


def leiden_key(resolution):
    """Stable obs key for one Leiden resolution, e.g. 0.25 -> leiden_res0_25."""
    return f"leiden_res{resolution:.2f}".replace(".", "_")


def make_savefig(plots_dir):
    """Return a ``savefig(name)`` that no-ops when no plots dir is configured."""

    def savefig(name):
        if not plots_dir:
            plt.close()
            return
        path = os.path.join(plots_dir, name)
        plt.savefig(path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {path}")

    return savefig


def map_var_names_to_symbols(adata, ensembl_to_symbol, symbol_to_ensembl):
    """Switch ``var_names`` from Ensembl IDs to gene symbols, if possible.

    Marker panels are specified by symbol, so lookups need symbols as the
    index. Genes whose symbol cannot be resolved are dropped.
    """
    if not str(adata.var_names[0]).startswith(("ENSG", "ENSMUSG")):
        print("var_names appear to already be gene symbols -- skipping mapping.")
        return adata

    symbol_map = load_symbol_map(ensembl_to_symbol)
    if symbol_map is None:
        print(
            f"Warning: ensembl->symbol dict not found at {ensembl_to_symbol}. "
            "Leaving var_names as Ensembl IDs; marker panels will not resolve."
        )
        return adata

    # Loaded for symmetry with the reverse lookup used by downstream steps.
    load_symbol_map(symbol_to_ensembl)

    adata.var["gene_id"] = adata.var_names.str.split(".").str[0]
    adata.var["symbol"] = adata.var["gene_id"].map(symbol_map)
    adata = adata[:, adata.var["symbol"].notna()].copy()
    adata.var_names = adata.var["symbol"].astype(str)
    adata.var_names_make_unique()
    print(f"After symbol mapping: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def score_cell_cycle(adata, s_genes_file, g2m_genes_file, random_state):
    """Score S and G2M phase from the shipped Tirosh et al. gene lists."""
    s_genes = read_gene_list(s_genes_file)
    g2m_genes = read_gene_list(g2m_genes_file)

    s_present = [g for g in s_genes if g in adata.var_names]
    g2m_present = [g for g in g2m_genes if g in adata.var_names]
    missing = len(s_genes) - len(s_present) + len(g2m_genes) - len(g2m_present)

    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_present,
        g2m_genes=g2m_present,
        random_state=random_state,
    )
    print(
        f"Cell cycle scoring done. S genes found: {len(s_present)}/{len(s_genes)}, "
        f"G2M genes found: {len(g2m_present)}/{len(g2m_genes)} "
        f"({missing} list entries absent from the data)"
    )


def score_marker_panels(adata, marker_panels, random_state):
    """Add one aggregate ``{panel}_score`` per marker panel."""
    for panel_name, genes in marker_panels.items():
        present = [g for g in genes if g in adata.var_names]
        missing = [g for g in genes if g not in adata.var_names]
        if missing:
            print(f"[{panel_name}] genes not found in data: {missing}")
        if not present:
            print(f"[{panel_name}] no genes found -- skipping score.")
            continue
        sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=f"{panel_name}_score",
            random_state=random_state,
        )
        print(f"[{panel_name}] scored {len(present)} genes.")


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["clustering"]
    input_paths = config["input_paths"]
    seed = params["random_seed"]

    np.random.seed(seed)

    plots_dir = args.plots_dir
    if plots_dir:
        os.makedirs(plots_dir, exist_ok=True)
    savefig = make_savefig(plots_dir)

    # --- 1. Load ------------------------------------------------------------
    print("Loading MuData...")
    mdata = mu.read(args.input)
    adata = mdata["GEX"].copy()
    adata.X = adata.X.astype("float32")
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- 2. Gene symbols ----------------------------------------------------
    adata = map_var_names_to_symbols(
        adata,
        input_paths.get("ensembl_to_symbol"),
        input_paths.get("symbol_to_ensembl"),
    )

    # --- 3. QC metrics, recomputed on the combined object -------------------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
    )

    # --- 4. Normalisation ---------------------------------------------------
    adata.layers["counts"] = adata.X.copy()  # preserve counts for seurat_v3 HVGs
    sc.pp.normalize_total(adata, target_sum=params["target_sum"])
    sc.pp.log1p(adata)
    print("Normalization complete.")

    # --- 5. Cell cycle ------------------------------------------------------
    score_cell_cycle(
        adata,
        input_paths["cell_cycle_s_genes"],
        input_paths["cell_cycle_g2m_genes"],
        random_state=seed,
    )

    # --- 6. Highly variable genes -------------------------------------------
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=params["n_top_genes"],
        batch_key=params["hvg_batch_key"] or None,
        flavor=params["hvg_flavor"],
        layer="counts",
    )
    print(f"HVGs selected: {adata.var['highly_variable'].sum()}")
    sc.pl.highly_variable_genes(adata, show=False)
    savefig("hvg_plot.png")

    # --- 6b. Optional cell-cycle regression ---------------------------------
    # Off for the published run. Turning it on also scales the data, which
    # changes PCA and every downstream clustering result.
    if params["regress_cell_cycle"]:
        print("Regressing out S_score and G2M_score, then scaling.")
        sc.pp.regress_out(adata, ["S_score", "G2M_score"])
        sc.pp.scale(adata, max_value=10)

    # --- 7. PCA -------------------------------------------------------------
    sc.tl.pca(
        adata,
        n_comps=params["n_pcs"],
        mask_var="highly_variable",
        random_state=seed,
    )
    sc.pl.pca_variance_ratio(adata, n_pcs=params["n_pcs"], log=True, show=False)
    savefig("pca_variance_ratio.png")
    sc.pl.pca(
        adata,
        color=["sublib", "pct_counts_mt", "total_counts", "phase"],
        ncols=2,
        size=3,
        show=False,
    )
    savefig("pca_qc_overview.png")

    # --- 8. Neighbour graph + UMAP ------------------------------------------
    sc.pp.neighbors(
        adata,
        n_pcs=params["n_pcs_use"],
        n_neighbors=params["n_neighbors"],
        random_state=seed,
    )
    sc.tl.umap(adata, random_state=seed)
    print("UMAP computed.")

    # --- 9. Leiden clustering -----------------------------------------------
    resolutions = params["leiden_resolutions"]
    for resolution in resolutions:
        key = leiden_key(resolution)
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key,
            flavor="igraph",
            n_iterations=2,
            random_state=seed,
        )
        print(f"Leiden res={resolution}: {adata.obs[key].nunique()} clusters")

    sc.pl.umap(
        adata,
        color=[leiden_key(r) for r in resolutions],
        legend_loc="on data",
        ncols=2,
        show=False,
    )
    savefig("umap_leiden_all_resolutions.png")

    # --- 10-11. UMAPs of QC metrics and cell cycle --------------------------
    sc.pl.umap(
        adata,
        color=[
            "sublib",
            "total_counts",
            "log1p_total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
        ],
        ncols=3,
        size=3,
        show=False,
    )
    savefig("umap_qc_metrics.png")

    sc.pl.umap(
        adata, color=["phase", "S_score", "G2M_score"], ncols=3, size=3, show=False
    )
    savefig("umap_cell_cycle.png")

    # --- 12. Marker panels --------------------------------------------------
    marker_panels = params["marker_panels"]
    score_marker_panels(adata, marker_panels, random_state=seed)

    score_cols = [
        f"{panel}_score"
        for panel in marker_panels
        if f"{panel}_score" in adata.obs.columns
    ]
    if score_cols:
        sc.pl.umap(
            adata,
            color=score_cols,
            cmap="bwr",
            vcenter=0,
            ncols=3,
            size=3,
            show=False,
        )
        savefig("umap_marker_scores.png")

    for panel_name, genes in marker_panels.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.pl.umap(adata, color=present, ncols=4, size=3, show=False, vmin=0)
        savefig(f"umap_{panel_name}_individual_genes.png")

    # --- 13. Marker dotplots by cluster -------------------------------------
    dotplot_key = leiden_key(params["dotplot_resolution"])
    if dotplot_key not in adata.obs.columns:
        raise ValueError(
            f"clustering.dotplot_resolution "
            f"({params['dotplot_resolution']}) is not one of "
            f"clustering.leiden_resolutions ({resolutions})"
        )
    for panel_name, genes in marker_panels.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.pl.dotplot(
            adata,
            var_names=present,
            groupby=dotplot_key,
            standard_scale="var",
            show=False,
        )
        savefig(f"dotplot_{panel_name}_by_{dotplot_key}.png")

    # --- 14. Save -----------------------------------------------------------
    adata.var.index.name = "gene_symbol"
    adata.uns["clustering_params"] = {
        key: value for key, value in params.items() if key != "marker_panels"
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"\nSaved clustered AnnData to:\n  {args.output}")
    print(f"Final: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
