#!/usr/bin/env python
"""Step 3a -- concatenate the per-sublibrary QC-filtered MuData objects.

Joins all sublibraries on the union of genes (zero-filling genes a sublibrary
lost to its own per-gene filters), tags every cell with its sublibrary of
origin, and asserts that GEX and guide stay row-aligned.

Barcodes are suffixed with the sublibrary label because the split-seq barcode
space is reused across sublibraries -- without the suffix, barcodes collide.
"""

import argparse
import os

import anndata as ad
import muon as mu
import pandas as pd
from muon import MuData

from scripts.pipeline_utils import load_config


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="QC-filtered .h5mu files, one per sublibrary",
    )
    parser.add_argument(
        "--sublibraries",
        nargs="+",
        required=True,
        help="Sublibrary labels, in the same order as --inputs",
    )
    parser.add_argument("--output", required=True, help="Output .h5mu path")
    parser.add_argument("--config", required=True, help="Workflow config YAML")
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["concatenate"]

    if len(args.inputs) != len(args.sublibraries):
        raise ValueError(
            f"--inputs ({len(args.inputs)}) and --sublibraries "
            f"({len(args.sublibraries)}) must have the same length"
        )

    print(f"Concatenating {len(args.inputs)} MuData files:")
    for label, path in zip(args.sublibraries, args.inputs):
        print(f"  {label}: {path}")

    gex_adatas = []
    guide_adatas = []

    for label, path in zip(args.sublibraries, args.inputs):
        print(f"\nLoading {path}  (sublibrary='{label}')")
        mdata = mu.read(path)
        gex = mdata["GEX"].copy()
        guide = mdata["guide"].copy()

        print(f"  GEX:   {gex.n_obs} cells x {gex.n_vars} genes")
        print(f"  guide: {guide.n_obs} cells x {guide.n_vars} guides")

        # Tag cells with their sublibrary of origin before any collision check.
        gex.obs["sublib"] = label
        guide.obs["sublib"] = label

        if params["make_barcodes_unique"]:
            gex.obs_names = [f"{bc}_{label}" for bc in gex.obs_names]
            guide.obs_names = [f"{bc}_{label}" for bc in guide.obs_names]

        gex_adatas.append(gex)
        guide_adatas.append(guide)

    print("\n--- Barcode collision check ---")
    for name, adatas in (("GEX", gex_adatas), ("guide", guide_adatas)):
        barcodes = pd.Series([bc for a in adatas for bc in a.obs_names])
        duplicates = barcodes[barcodes.duplicated()].unique()
        if len(duplicates) == 0:
            print(f"{name}:   no duplicate barcodes")
        else:
            print(f"{name}:   {len(duplicates)} duplicate barcodes found!")
            print(f"  Examples: {list(duplicates[:10])}")
        assert len(duplicates) == 0, (
            f"Duplicate {name} barcodes detected. Set "
            "concatenate.make_barcodes_unique: true, or investigate the inputs."
        )

    print(f"\n--- Concatenating GEX (join={params['join']}) ---")
    gex_combined = ad.concat(
        gex_adatas,
        join=params["join"],
        label="sublib",
        keys=list(args.sublibraries),
        index_unique=None,  # barcodes are already unique
    )
    print(f"Combined GEX: {gex_combined.n_obs} cells x {gex_combined.n_vars} genes")

    print(f"\n--- Concatenating guide (join={params['join']}) ---")
    guide_combined = ad.concat(
        guide_adatas,
        join=params["join"],
        label="sublib",
        keys=list(args.sublibraries),
        index_unique=None,
    )
    print(
        f"Combined guide: {guide_combined.n_obs} cells x "
        f"{guide_combined.n_vars} guides"
    )

    print("\n--- Post-concatenation sanity checks ---")
    assert gex_combined.n_obs == guide_combined.n_obs, (
        f"Cell count mismatch after concat: GEX={gex_combined.n_obs}, "
        f"guide={guide_combined.n_obs}"
    )
    assert (gex_combined.obs_names == guide_combined.obs_names).all(), (
        "Barcode order differs between combined GEX and guide objects."
    )

    print("Cells per sublibrary:")
    print(gex_combined.obs["sublib"].value_counts().to_string())
    print("\nAll checks passed.")

    mdata_combined = MuData({"GEX": gex_combined, "guide": guide_combined})

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    mdata_combined.write(args.output)
    print(f"\nSaved combined MuData to:\n  {args.output}")
    print(f"Final: {gex_combined.n_obs} cells x {gex_combined.n_vars} genes")


if __name__ == "__main__":
    main()
