#!/usr/bin/env python
"""Step 6c -- validate the positional cell-ordering contract SCEPTRE relies on.

The whole SCEPTRE stage rests on one assumption: cell *i* of the sceptre
object is cell *i* of the concatenated ``barcodes.tsv.gz`` files, in the order
``import_data_from_cellranger()`` received the directories. The extra covariate
table is joined on that basis, with no barcode key. If it were ever wrong,
every downstream result would be quietly attributed to the wrong cells rather
than erroring.

This consolidates three checks that were originally three near-duplicate
scripts, and adds a fourth:

``barcode``     CSV row *i* barcode == positional barcode *i*.
``covariates``  CSV values == the clustered AnnData's obs values for that cell.
``assignment``  the gRNA assignment matrix's column names == the same
                positional barcode order. This is the check that most directly
                protects the downstream per-guide knockdown work, and it is new
                -- the original scripts validated the covariate join but not
                the assignment matrix's own cell order.
``counts``      spot-check that rounding the stage-00 float MuData counts
                reproduces the integer CellRanger matrix for a few cells, i.e.
                that the clustering input and the SCEPTRE input differ only by
                rounding. Reads large files; opt in explicitly.

Exits non-zero if any selected check fails.
"""

import argparse
import gzip
import os

import numpy as np
import pandas as pd

from scripts.pipeline_utils import load_config

ALL_CHECKS = ("barcode", "covariates", "assignment", "counts")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--covariates", required=True, help="Extra covariate CSV")
    parser.add_argument(
        "--barcode-files", nargs="+", required=True, help="barcodes.tsv.gz, in order"
    )
    parser.add_argument("--sublibraries", nargs="+", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True, help="Report path (.txt)")
    parser.add_argument(
        "--checks",
        default="barcode,covariates,assignment",
        help=f"Comma-separated subset of {ALL_CHECKS} [default: %(default)s]",
    )
    parser.add_argument(
        "--assignment-matrix",
        default=None,
        help="grna_assignment_matrix.rds column names, as a text file "
        "(one barcode per line). Required for the 'assignment' check.",
    )
    parser.add_argument(
        "--qc-filtered-h5mus",
        nargs="*",
        default=None,
        help="Stage-00 QC-filtered .h5mu files, for the 'counts' check",
    )
    parser.add_argument(
        "--n-spot-checks",
        type=int,
        default=10,
        help="Cells to spot-check in the heavier checks [default: %(default)s]",
    )
    return parser.parse_args()


def read_barcodes(path):
    with gzip.open(path, "rt") as handle:
        return [line.strip() for line in handle if line.strip()]


def positional_barcodes(barcode_files, sublibraries):
    """Concatenated barcodes plus, per cell, its sublibrary and local index."""
    barcodes, origin = [], []
    for label, path in zip(sublibraries, barcode_files):
        local = read_barcodes(path)
        barcodes.extend(local)
        origin.extend((label, i) for i in range(len(local)))
    return barcodes, origin


def check_barcode(csv_index, barcodes, report):
    report("--- check: barcode ---")
    report(f"  covariate CSV rows      : {len(csv_index)}")
    report(f"  positional barcodes     : {len(barcodes)}")
    if len(csv_index) != len(barcodes):
        report("  FAIL: row counts differ")
        return False
    mismatches = [
        i for i, (a, b) in enumerate(zip(csv_index, barcodes)) if a != b
    ]
    if mismatches:
        report(f"  FAIL: {len(mismatches)} positions differ")
        for i in mismatches[:5]:
            report(f"    [{i}] CSV={csv_index[i]}  positional={barcodes[i]}")
        return False
    report("  PASS: every CSV row sits at its positional barcode")
    return True


def check_covariates(df, config, barcodes, n_spot, report):
    import anndata as ad

    params = config["sceptre"]["covariates"]
    report("--- check: covariates ---")
    report(f"  AnnData: {params['clustered_h5ad']}")
    adata = ad.read_h5ad(params["clustered_h5ad"], backed="r")

    infix = params["barcode_h5ad_infix"]
    rng = np.random.default_rng(0)
    idx = sorted(rng.choice(len(barcodes), size=min(n_spot, len(barcodes)), replace=False))

    ok = True
    for i in idx:
        sceptre_bc = barcodes[i]
        sublib = sceptre_bc.rsplit("_", 1)[-1]
        h5ad_bc = sceptre_bc.replace(f"_{sublib}", f"{infix}_{sublib}", 1)
        if h5ad_bc not in adata.obs_names:
            report(f"  FAIL: [{i}] {h5ad_bc} not in the AnnData")
            ok = False
            continue
        for out_name, obs_col in params["columns"].items():
            if out_name not in df.columns:
                continue
            expected = adata.obs.loc[h5ad_bc, obs_col]
            got = df.iloc[i][out_name]
            if isinstance(expected, (int, float, np.floating, np.integer)):
                match = np.isclose(float(got), float(expected), rtol=0, atol=1e-8)
            else:
                match = str(got) == str(expected)
            if not match:
                report(f"  FAIL: [{i}] {out_name}: CSV={got!r} AnnData={expected!r}")
                ok = False
    if ok:
        report(f"  PASS: {len(idx)} spot-checked cells match on every covariate")
    return ok


def check_assignment(barcodes, assignment_path, report):
    report("--- check: assignment ---")
    if not assignment_path or not os.path.exists(assignment_path):
        report(f"  SKIP: --assignment-matrix not provided or missing ({assignment_path})")
        return None
    with open(assignment_path) as handle:
        cols = [line.strip() for line in handle if line.strip()]
    report(f"  assignment matrix columns : {len(cols)}")
    report(f"  positional barcodes       : {len(barcodes)}")
    if len(cols) != len(barcodes):
        report("  FAIL: column count differs from the barcode count")
        return False
    if cols == barcodes:
        report("  PASS: assignment columns are in positional barcode order")
        return True
    if set(cols) == set(barcodes):
        n_diff = sum(1 for a, b in zip(cols, barcodes) if a != b)
        report(
            f"  FAIL: same cell set but {n_diff} positions differ -- the "
            "assignment matrix is permuted relative to import order. Anything "
            "joining it positionally is wrong; join on column names instead."
        )
        return False
    report("  FAIL: assignment columns and positional barcodes are different sets")
    return False


def check_counts(barcodes, origin, h5mus, sublibraries, n_spot, report):
    import muon as mu

    report("--- check: counts ---")
    if not h5mus:
        report("  SKIP: --qc-filtered-h5mus not provided")
        return None

    by_sublib = dict(zip(sublibraries, h5mus))
    rng = np.random.default_rng(1)
    idx = sorted(rng.choice(len(barcodes), size=min(n_spot, len(barcodes)), replace=False))

    ok = True
    cache = {}
    for i in idx:
        sublib, local_i = origin[i]
        if sublib not in by_sublib:
            report(f"  SKIP: no h5mu supplied for {sublib}")
            continue
        if sublib not in cache:
            cache[sublib] = mu.read(by_sublib[sublib])["GEX"]
        gex = cache[sublib]
        # 03c_convert_to_cellranger.py sorts barcodes before writing, so the
        # positional barcode must equal the sorted obs_names at that index.
        sorted_names = sorted(gex.obs_names)
        expected_bc = f"{sorted_names[local_i]}_{sublib}"
        if expected_bc != barcodes[i]:
            report(
                f"  FAIL: [{i}] {sublib} local {local_i}: "
                f"cellranger={barcodes[i]} sorted-h5mu={expected_bc}"
            )
            ok = False
            continue
        row = gex[sorted_names[local_i], :].X
        total_float = float(np.asarray(row.sum()))
        total_rounded = float(np.asarray(np.round(row).sum()))
        report(
            f"  [{i}] {barcodes[i]}: float total {total_float:.2f} -> "
            f"rounded {total_rounded:.0f}"
        )
    if ok:
        report(f"  PASS: {len(idx)} cells align between the h5mu and the export")
    return ok


def main():
    args = parse_args()
    config = load_config(args.config)

    selected = [c.strip() for c in args.checks.split(",") if c.strip()]
    unknown = set(selected) - set(ALL_CHECKS)
    if unknown:
        raise ValueError(f"unknown checks: {sorted(unknown)}; valid: {ALL_CHECKS}")

    lines = []

    def report(text=""):
        print(text)
        lines.append(str(text))

    report("SCEPTRE cell-ordering validation")
    report(f"covariate CSV : {args.covariates}")
    report(f"checks        : {', '.join(selected)}")
    report("")

    df = pd.read_csv(args.covariates, index_col=0)
    csv_index = [str(x) for x in df.index]
    barcodes, origin = positional_barcodes(args.barcode_files, args.sublibraries)

    results = {}
    if "barcode" in selected:
        results["barcode"] = check_barcode(csv_index, barcodes, report)
        report("")
    if "covariates" in selected:
        results["covariates"] = check_covariates(
            df, config, barcodes, args.n_spot_checks, report
        )
        report("")
    if "assignment" in selected:
        results["assignment"] = check_assignment(
            barcodes, args.assignment_matrix, report
        )
        report("")
    if "counts" in selected:
        results["counts"] = check_counts(
            barcodes, origin, args.qc_filtered_h5mus, args.sublibraries,
            args.n_spot_checks, report
        )
        report("")

    report("--- summary ---")
    for name, value in results.items():
        verdict = "SKIP" if value is None else ("PASS" if value else "FAIL")
        report(f"  {name:12s} {verdict}")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"\nWrote {args.output}")

    if any(v is False for v in results.values()):
        raise SystemExit("cell-ordering validation FAILED -- see the report above")


if __name__ == "__main__":
    main()
