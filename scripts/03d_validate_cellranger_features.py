#!/usr/bin/env python
"""Step 3b-i, pass 3 -- confirm every sublibrary emitted identical features.

SCEPTRE's ``import_data_from_cellranger`` reads a list of
``filtered_feature_bc_matrix`` directories and indexes features positionally,
so it silently produces wrong results if the directories disagree on
``features.tsv.gz``. This rule fails loudly instead.

Writes a short report listing the feature counts per sublibrary and, on
success, the shared feature order that downstream steps can rely on.
"""

import argparse
import gzip
import os


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--features-files",
        nargs="+",
        required=True,
        help="features.tsv.gz files, one per sublibrary",
    )
    parser.add_argument(
        "--sublibraries",
        nargs="+",
        required=True,
        help="Sublibrary labels, in the same order as --features-files",
    )
    parser.add_argument(
        "--output", required=True, help="Output report path (.txt)"
    )
    parser.add_argument(
        "--reference-features",
        default=None,
        help="Optional path to write the shared feature ID list to",
    )
    return parser.parse_args()


def read_features(path):
    """Return the list of (id, name, type) tuples from a features.tsv.gz."""
    with gzip.open(path, "rt") as handle:
        return [tuple(line.rstrip("\n").split("\t")) for line in handle if line.strip()]


def main():
    args = parse_args()

    if len(args.features_files) != len(args.sublibraries):
        raise ValueError(
            "--features-files and --sublibraries must have the same length"
        )

    per_sublibrary = {}
    for label, path in zip(args.sublibraries, args.features_files):
        features = read_features(path)
        per_sublibrary[label] = features
        n_gene = sum(1 for f in features if f[2] == "Gene Expression")
        n_guide = sum(1 for f in features if f[2] == "CRISPR Guide Capture")
        print(
            f"{label}: {len(features)} features "
            f"({n_gene} genes + {n_guide} guides) -- {path}"
        )

    reference_label = args.sublibraries[0]
    reference = per_sublibrary[reference_label]

    mismatches = []
    for label, features in per_sublibrary.items():
        if features == reference:
            continue
        if len(features) != len(reference):
            mismatches.append(
                f"{label}: {len(features)} features vs {len(reference)} "
                f"in {reference_label}"
            )
            continue
        first_diff = next(
            i for i, (a, b) in enumerate(zip(features, reference)) if a != b
        )
        mismatches.append(
            f"{label}: first difference at row {first_diff}: "
            f"{features[first_diff]} vs {reference[first_diff]} "
            f"in {reference_label}"
        )

    lines = [
        "CellRanger feature consistency check",
        f"reference sublibrary: {reference_label}",
        f"features per sublibrary: {len(reference)}",
        "",
    ]
    for label, features in per_sublibrary.items():
        n_gene = sum(1 for f in features if f[2] == "Gene Expression")
        n_guide = sum(1 for f in features if f[2] == "CRISPR Guide Capture")
        lines.append(f"{label}\t{len(features)}\t{n_gene}\t{n_guide}")
    lines.append("")

    if mismatches:
        lines.append("FAILED")
        lines.extend(mismatches)
    else:
        lines.append("PASSED: all sublibraries emit identical features.tsv.gz")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    if mismatches:
        raise SystemExit(
            "features.tsv.gz differs between sublibraries -- rerun the export "
            "with --master-gene-list so SCEPTRE can index features positionally."
        )

    if args.reference_features:
        os.makedirs(
            os.path.dirname(os.path.abspath(args.reference_features)), exist_ok=True
        )
        with open(args.reference_features, "w") as handle:
            for feature in reference:
                handle.write(feature[0] + "\n")
        print(f"Shared feature IDs written to {args.reference_features}")


if __name__ == "__main__":
    main()
