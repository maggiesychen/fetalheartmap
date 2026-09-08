#!/usr/bin/env python
"""Step 3b-i, pass 1 -- build the gene list shared by every sublibrary.

``qc_metrics_filter.py`` applies its per-gene filters one sublibrary at a time,
so the sublibraries end up with different gene sets. Anything that reads the
CellRanger-format matrices as a set (SCEPTRE's
``import_data_from_cellranger``, which requires identical ``features.tsv.gz``
across directories) needs them reconciled first.

This writes the sorted union of unversioned Ensembl gene IDs across all
sublibraries. ``convert_to_cellranger.py --master-gene-list`` then reindexes
each sublibrary onto it, zero-filling the genes that sublibrary lost.
"""

import argparse
import os

import muon as mu


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="QC-filtered .h5mu files, one per sublibrary",
    )
    parser.add_argument(
        "--output", required=True, help="Output master gene ID list (one per line)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Building master gene list from {len(args.inputs)} sublibraries...")
    all_gene_ids = set()
    for path in args.inputs:
        mdata = mu.read(path)
        # Strip Ensembl version suffixes so IDs are comparable across
        # sublibraries and match what convert_to_cellranger.py writes.
        gene_ids = [g.split(".")[0] for g in mdata["GEX"].var_names.tolist()]
        all_gene_ids.update(gene_ids)
        print(f"  {os.path.basename(path)}: {len(gene_ids)} genes")

    master_gene_ids = sorted(all_gene_ids)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as handle:
        for gene_id in master_gene_ids:
            handle.write(gene_id + "\n")

    print(f"Master gene list: {len(master_gene_ids)} genes -> {args.output}")


if __name__ == "__main__":
    main()
