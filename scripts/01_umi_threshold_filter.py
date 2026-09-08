#!/usr/bin/env python
"""Step 1 -- filter one sublibrary to a UMI window and pair GEX with guides.

Reads the raw kb-python GEX and guide count matrices for a single sublibrary,
keeps cells whose total GEX UMI count falls inside ``[--umi-min, --umi-max]``,
intersects those barcodes with the guide matrix, and writes a two-modality
MuData (``GEX`` and ``guide``).

The lower bound is the barcode-rank inflection point reported by the upstream
perturb-seq pipeline's cell calling; the upper bound removes probable multiplets.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from muon import MuData  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gex-input", required=True, help="Raw GEX adata.h5ad")
    parser.add_argument("--guide-input", required=True, help="Raw guide adata.h5ad")
    parser.add_argument("--output", required=True, help="Output .h5mu path")
    parser.add_argument(
        "--umi-min", type=int, required=True, help="Minimum total GEX UMIs per cell"
    )
    parser.add_argument(
        "--umi-max",
        type=int,
        default=None,
        help="Maximum total GEX UMIs per cell (default: no upper bound)",
    )
    parser.add_argument(
        "--rank-curve-plot",
        default=None,
        help="Optional path for the UMI rank-curve PNG",
    )
    parser.add_argument(
        "--round-counts",
        action="store_true",
        help=(
            "Round X to int32 before writing. The published run did NOT do this "
            "-- rounding happens at CellRanger export instead."
        ),
    )
    parser.add_argument(
        "--check-layer-sums",
        action="store_true",
        help="Assert nascent + ambiguous + mature == X (diagnostic)",
    )
    return parser.parse_args()


def plot_rank_curve(knee_df, umi_min, umi_max, path):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(knee_df["rank"], knee_df["total_counts"], linewidth=1, color="steelblue")
    ax.axhline(
        y=umi_min, color="red", linestyle="--", label=f"min threshold = {umi_min}"
    )
    if umi_max is not None:
        ax.axhline(
            y=umi_max,
            color="darkorange",
            linestyle="--",
            label=f"max threshold = {umi_max}",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cell rank (by UMI count)", fontsize=12)
    ax.set_ylabel("UMI count", fontsize=12)
    ax.set_title("UMI rank curve", fontsize=13)
    ax.legend()
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Saved: {path}")


def check_layer_sums(adata):
    """Confirm the kb-python layers decompose X exactly."""
    required = ("nascent", "ambiguous", "mature")
    missing = [name for name in required if name not in adata.layers]
    if missing:
        print(f"Layer check skipped -- missing layers: {missing}")
        return
    layer_sum = sum(adata.layers[name] for name in required)
    diff = abs(layer_sum - adata.X)
    max_diff = diff.max() if hasattr(diff, "max") else np.max(diff)
    print(f"Layer check: max |nascent + ambiguous + mature - X| = {max_diff}")
    assert max_diff < 1e-3, "kb-python layers do not sum to X"


def main():
    args = parse_args()

    print(
        f"UMI filtering | umi_min={args.umi_min} | "
        f"umi_max={args.umi_max if args.umi_max is not None else 'none'}"
    )

    gex_adata = sc.read(args.gex_input)
    guide_adata = sc.read(args.guide_input)
    gex_adata.X = gex_adata.X.astype(np.float32)

    if args.check_layer_sums:
        check_layer_sums(gex_adata)

    umi_counts = np.array(gex_adata.X.sum(axis=1)).flatten()
    gex_adata.obs["total_counts"] = umi_counts

    knee_df = pd.DataFrame(
        {"total_counts": umi_counts, "barcode": gex_adata.obs_names.values}
    )
    knee_df = knee_df.sort_values("total_counts", ascending=False).reset_index(drop=True)
    knee_df["rank"] = knee_df.index + 1

    print(f"Total cells before UMI filtering: {len(knee_df)}")
    print(
        f"Cells with total_counts >= {args.umi_min}: "
        f"{(knee_df['total_counts'] >= args.umi_min).sum()}"
    )

    if args.rank_curve_plot:
        plot_rank_curve(knee_df, args.umi_min, args.umi_max, args.rank_curve_plot)

    keep = gex_adata.obs["total_counts"] >= args.umi_min
    if args.umi_max is not None:
        keep &= gex_adata.obs["total_counts"] <= args.umi_max
    gex_adata = gex_adata[keep, :].copy()
    print(f"Cells after UMI filtering: {gex_adata.n_obs}")

    # Guide barcodes are a superset of GEX barcodes; keep the intersection so
    # the two modalities stay row-aligned.
    shared_barcodes = sorted(
        set(gex_adata.obs_names).intersection(guide_adata.obs_names)
    )
    print(f"Barcodes present in both GEX and guide: {len(shared_barcodes)}")

    gex_out = gex_adata[shared_barcodes, :].copy()
    guide_out = guide_adata[shared_barcodes, :].copy()

    if args.round_counts:
        gex_out.X = np.round(gex_out.X).astype(np.int32)
        guide_out.X = np.round(guide_out.X).astype(np.int32)
        print("Counts rounded to int32.")

    mdata = MuData({"GEX": gex_out, "guide": guide_out})

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    mdata.write(args.output)
    print(f"Saved MuData to: {args.output}")
    print(f"GEX dtype: {gex_out.X.dtype} | guide dtype: {guide_out.X.dtype}")


if __name__ == "__main__":
    main()
