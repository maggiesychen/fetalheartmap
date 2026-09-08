#!/usr/bin/env python
"""Steps 3b / 3b-i -- write one sublibrary as a CellRanger feature-barcode matrix.

Emits ``barcodes.tsv.gz``, ``features.tsv.gz`` and ``matrix.mtx.gz`` with genes
stacked above guides, which is the layout SCEPTRE's
``import_data_from_cellranger`` expects.

Counts are float32 up to this point; they are **rounded** (not truncated) to
int32 here. Truncating zeroes every sub-1.0 count and silently drops cells and
genes -- see the README.

Two modes:

* without ``--master-gene-list`` (step 3b): each sublibrary keeps its own gene
  set, so ``features.tsv.gz`` differs between sublibraries;
* with ``--master-gene-list`` (step 3b-i): the sublibrary is reindexed onto the
  shared gene list, zero-filling missing genes, so every sublibrary emits an
  identical ``features.tsv.gz``. This is the mode downstream steps consume.
"""

import argparse
import gzip
import os

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

from scripts.pipeline_utils import load_config


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="QC-filtered .h5mu")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output filtered_feature_bc_matrix/ directory",
    )
    parser.add_argument("--config", required=True, help="Workflow config YAML")
    parser.add_argument(
        "--sublibrary",
        required=True,
        help="Sublibrary label, appended to every barcode",
    )
    parser.add_argument(
        "--master-gene-list",
        default=None,
        help="Optional shared gene ID list from build_master_gene_list.py",
    )
    return parser.parse_args()


def read_gene_ids(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def reindex_to_master(gex, master_gene_ids):
    """Add the genes this sublibrary lost as explicit zeros, then reorder."""
    present = set(gex.var_names)
    missing_genes = [g for g in master_gene_ids if g not in present]
    print(f"Missing genes in this sublibrary: {len(missing_genes)}")

    if missing_genes:
        empty = ad.AnnData(
            X=sp.csc_matrix((gex.n_obs, len(missing_genes))),
            obs=gex.obs,
            var=pd.DataFrame(index=missing_genes),
        )
        gex = ad.concat([gex, empty], axis=1)

    gex = gex[:, master_gene_ids]
    print(f"GEX after reindex: {gex.n_obs} cells x {gex.n_vars} genes")
    assert gex.n_vars == len(master_gene_ids), "Reindex failed -- gene count mismatch"
    assert list(gex.var_names) == master_gene_ids, (
        "Reindex failed -- gene order mismatch"
    )
    return gex, missing_genes


def build_features(gex, guide, gene_ids, config):
    """Gene features stacked above guide features, CellRanger 3-column format."""
    params = config["cellranger_export"]

    # Prefer the mapped symbol; fall back to the Ensembl ID when the mapping
    # failed (symbol is NaN, or is itself an unversioned Ensembl ID).
    if "symbol" in gex.var.columns:
        symbols = gex.var["symbol"].fillna("").tolist()
        gene_names = [
            symbol
            if (isinstance(symbol, str) and symbol and not symbol.startswith("ENSG"))
            else gene_id
            for symbol, gene_id in zip(symbols, gene_ids)
        ]
    else:
        gene_names = list(gene_ids)

    gene_features = pd.DataFrame(
        {"id": gene_ids, "name": gene_names, "type": params["gene_feature_type"]}
    )
    guide_features = pd.DataFrame(
        {
            "id": guide.var_names.tolist(),
            "name": guide.var_names.tolist(),
            "type": params["guide_feature_type"],
        }
    )
    features = pd.concat([gene_features, guide_features], ignore_index=True)

    # Stripping Ensembl version suffixes could in principle collide two IDs.
    assert features["id"].is_unique, "Duplicate feature IDs found!"
    return features, len(gene_features), len(guide_features)


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["cellranger_export"]

    out_dir = args.output_dir.rstrip("/") + "/"
    os.makedirs(out_dir, exist_ok=True)
    print(f"Processing {args.sublibrary}")
    print(f"Input:  {args.input}")
    print(f"Output: {out_dir}")

    mdata = mu.read(args.input)
    gex = mdata["GEX"].copy()
    guide = mdata["guide"].copy()

    if params["strip_ensembl_versions"]:
        gex.var_names = pd.Index([g.split(".")[0] for g in gex.var_names.tolist()])

    # Round, do not truncate: truncation zeroes every sub-1.0 count.
    gex.X = np.round(gex.X).astype(np.int32)
    guide.X = np.round(guide.X).astype(np.int32)

    print(f"GEX:   {gex.n_obs} cells x {gex.n_vars} genes ({gex.X.dtype})")
    print(f"Guide: {guide.n_obs} cells x {guide.n_vars} guides ({guide.X.dtype})")

    missing_genes = []
    if args.master_gene_list:
        master_gene_ids = read_gene_ids(args.master_gene_list)
        print(f"Master gene list: {len(master_gene_ids)} genes")
        gex, missing_genes = reindex_to_master(gex, master_gene_ids)
        gene_ids = master_gene_ids
    else:
        gene_ids = gex.var_names.tolist()

    # Sort barcodes so the matrix column order is deterministic and identical
    # between the GEX and guide blocks.
    sorted_barcodes = sorted(gex.obs_names)
    gex = gex[sorted_barcodes, :]
    guide = guide[sorted_barcodes, :]
    assert (gex.obs_names == guide.obs_names).all(), "Barcode mismatch after sorting!"
    print(f"Barcodes sorted and matched: {len(sorted_barcodes)}")

    with gzip.open(out_dir + "barcodes.tsv.gz", "wt") as handle:
        for barcode in sorted_barcodes:
            handle.write(f"{barcode}_{args.sublibrary}\n")
    print("Written: barcodes.tsv.gz")

    features, n_gene_features, n_guide_features = build_features(
        gex, guide, gene_ids, config
    )
    with gzip.open(out_dir + "features.tsv.gz", "wt") as handle:
        features.to_csv(handle, sep="\t", header=False, index=False)
    print(
        f"Written: features.tsv.gz ({n_gene_features} genes + "
        f"{n_guide_features} guides)"
    )

    gene_matrix = sp.csc_matrix(gex.X).T  # genes x cells
    guide_matrix = sp.csc_matrix(guide.X).T  # guides x cells
    combined = sp.vstack([gene_matrix, guide_matrix])
    print(f"Combined matrix shape: {combined.shape} (features x cells)")

    # Drop the explicit zeros introduced by rounding and by reindexing.
    combined.eliminate_zeros()
    print(f"After eliminate_zeros -- nnz: {combined.nnz}")

    with gzip.open(out_dir + "matrix.mtx.gz", "wb") as handle:
        sio.mmwrite(handle, combined)
    print("Written: matrix.mtx.gz")

    print(f"\n{args.sublibrary} matrix stats:")
    print(f"  shape: {combined.shape}")
    print(f"  nnz: {combined.nnz}")
    print(f"  total UMI count: {combined.sum()}")
    print(f"  mean UMIs per cell: {combined.sum() / combined.shape[1]:.1f}")

    spot_check_gene = params.get("spot_check_gene_id")
    if spot_check_gene and spot_check_gene in gene_ids:
        gene_counts = combined[gene_ids.index(spot_check_gene), :].toarray().flatten()
        print(
            f"  {spot_check_gene}: mean={gene_counts.mean():.2f}, "
            f"max={gene_counts.max()}, nonzero={np.count_nonzero(gene_counts)}"
        )

    print("\n--- Sanity checks ---")
    assert combined.shape[0] == len(features), (
        f"Feature count mismatch: matrix {combined.shape[0]} vs "
        f"features {len(features)}"
    )
    assert combined.shape[1] == len(sorted_barcodes), (
        f"Barcode count mismatch: matrix {combined.shape[1]} vs "
        f"{len(sorted_barcodes)} barcodes"
    )

    # Spot-check a gene that really exists in this sublibrary, not a
    # zero-filled one, against the source matrix.
    gex_transposed = sp.csc_matrix(gex.X).T
    present_genes = [g for g in gene_ids if g not in set(missing_genes)]
    if present_genes:
        test_idx = gene_ids.index(present_genes[0])
        assert np.allclose(
            combined[test_idx, :100].toarray(), gex_transposed[test_idx, :100].toarray()
        ), "Gene block spot check failed!"
        print(f"Gene block spot check passed for {present_genes[0]}")

    guide_block = combined[len(gene_ids) : len(gene_ids) + 100, :100].toarray()
    guide_source = sp.csc_matrix(guide.X).T[:100, :100].toarray()
    assert np.allclose(guide_block, guide_source), "Guide block spot check failed!"
    print("Guide block spot check passed")

    if combined.data.size > 0:
        assert combined.data.min() >= 0, "Negative values found in matrix!"

    assert n_gene_features == len(gene_ids), (
        f"Gene feature count mismatch: {n_gene_features} vs {len(gene_ids)}"
    )
    assert n_guide_features == guide.n_vars, (
        f"Guide feature count mismatch: {n_guide_features} vs {guide.n_vars}"
    )

    print(f"Gene features:  {n_gene_features}")
    print(f"Guide features: {n_guide_features}")
    print(f"Matrix shape:   {combined.shape}")
    print(f"Matrix dtype:   {combined.dtype}")
    print(f"Matrix nnz:     {combined.nnz}")
    print("All sanity checks passed.")


if __name__ == "__main__":
    main()
