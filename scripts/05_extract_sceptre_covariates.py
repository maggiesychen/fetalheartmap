#!/usr/bin/env python
"""Step 5 -- build the extra per-cell covariate table for SCEPTRE.

SCEPTRE's own covariates cover sequencing depth and batch. Cell-state
covariates (cell-cycle scores, ribosomal and mitochondrial fraction, Leiden
cluster) have to be supplied separately, and they must be in **exactly** the
cell order SCEPTRE used at import: the concatenation of
``barcodes.tsv.gz`` across sublibraries, in the order
``import_data_from_cellranger()`` received the directories.

Row *i* of the output therefore corresponds to SCEPTRE integer cell index *i*.
That positional contract is what ``06c_validate_cell_ordering.py`` checks, and
it is safe here because cell order is stable across SCEPTRE runs (only gRNA
order is permuted between runs).

Barcode formats differ between the two sources and are bridged by an infix:

    CellRanger / SCEPTRE : {barcode}_{sublibrary}
    clustered AnnData    : {barcode}{infix}_{sublibrary}

Note on which clustering to read: ``S_score``, ``G2M_score``, ``pct_ribo`` and
``pct_mito`` are identical between the cell-cycle-regressed and
no-regression clusterings (verified 2026-09-08, max abs diff 0), so for the
trans_v2 variant either file reproduces the published covariates. Only
``leiden_res0_25`` differs -- 6 clusters vs 5, with 83.5% of cells relabelled --
so the trans_v3 variant is specific to the clustering named in the config.
"""

import argparse
import gzip
import os

import anndata as ad
import pandas as pd

from scripts.pipeline_utils import load_config


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--barcode-files",
        nargs="+",
        required=True,
        help="barcodes.tsv.gz per sublibrary, in SCEPTRE import order",
    )
    parser.add_argument(
        "--sublibraries",
        nargs="+",
        required=True,
        help="Sublibrary labels, same order as --barcode-files",
    )
    parser.add_argument("--output", required=True, help="Output covariate CSV")
    parser.add_argument("--config", required=True, help="Workflow config YAML")
    return parser.parse_args()


def read_barcodes(path):
    with gzip.open(path, "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["sceptre"]["covariates"]

    if len(args.barcode_files) != len(args.sublibraries):
        raise ValueError("--barcode-files and --sublibraries must be the same length")

    # --- 1. Barcodes in SCEPTRE import order --------------------------------
    print("Reading barcodes in SCEPTRE import order...")
    sceptre_barcodes = []
    for label, path in zip(args.sublibraries, args.barcode_files):
        barcodes = read_barcodes(path)
        print(f"  {label}: {len(barcodes)} barcodes  (first: {barcodes[0]})")
        sceptre_barcodes.extend(barcodes)
    print(f"Total SCEPTRE cells: {len(sceptre_barcodes)}")

    if len(set(sceptre_barcodes)) != len(sceptre_barcodes):
        raise ValueError(
            "Duplicate barcodes across sublibraries -- the sublibrary suffix "
            "is not making them unique, so the positional contract is unsafe."
        )

    # --- 2. Translate to the clustered AnnData's barcode format -------------
    infix = params["barcode_h5ad_infix"]
    h5ad_barcodes = [
        bc.replace(f"_{label}", f"{infix}_{label}", 1)
        for bc, label in (
            (bc, next(s for s in args.sublibraries if bc.endswith(f"_{s}")))
            for bc in sceptre_barcodes
        )
    ]

    # --- 3. Load the clustered AnnData and check coverage -------------------
    h5ad_path = params["clustered_h5ad"]
    print(f"\nLoading clustered AnnData (backed):\n  {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path, backed="r")

    present = set(adata.obs_names)
    n_found = sum(bc in present for bc in h5ad_barcodes)
    print(f"Barcodes matched in the AnnData: {n_found} / {len(h5ad_barcodes)}")
    if n_found != len(h5ad_barcodes):
        missing = [bc for bc in h5ad_barcodes if bc not in present][:5]
        raise RuntimeError(
            f"Only {n_found}/{len(h5ad_barcodes)} barcodes found. Check "
            f"sceptre.covariates.barcode_h5ad_infix. Examples missing: {missing}"
        )

    # --- 4. Extract, in SCEPTRE cell order ----------------------------------
    columns = params["columns"]  # output name -> AnnData obs column
    missing_cols = [src for src in columns.values() if src not in adata.obs.columns]
    if missing_cols:
        raise ValueError(
            f"AnnData obs is missing configured columns: {missing_cols}. "
            f"Available: {sorted(adata.obs.columns)}"
        )

    print("Extracting covariates in SCEPTRE cell order...")
    df = adata.obs.loc[h5ad_barcodes, list(columns.values())].copy()
    df.columns = list(columns.keys())

    # Categorical columns are written as strings; the R side coerces the ones
    # listed in the variant's `categorical_covariates` back to factors.
    for name in df.columns:
        if str(df[name].dtype) == "category" or df[name].dtype == object:
            df[name] = df[name].astype(str)

    # Index by the SCEPTRE barcode for traceability. Position is what the R
    # side relies on; the index is there so mistakes are auditable.
    df.index = sceptre_barcodes
    df.index.name = "sceptre_barcode"

    for name in df.columns:
        n_na = int(df[name].isna().sum())
        print(f"  NAs in {name}: {n_na}")
        if n_na:
            raise ValueError(
                f"{name} has {n_na} missing values; SCEPTRE cannot fit a model "
                "with NA covariates."
            )

    for name in df.columns:
        if df[name].dtype == object:
            counts = df[name].value_counts().sort_index()
            print(f"\n{name} value counts:\n{counts.to_string()}")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output)
    print(
        f"\nSaved {df.shape[0]} cells x {df.shape[1]} covariates to:\n  {args.output}"
    )


if __name__ == "__main__":
    main()
