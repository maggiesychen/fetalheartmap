#!/usr/bin/env python
"""Step 7e -- combined percent-knockdown bar figure.

Merges the per-gene summary from 07d (column order, asterisks, headline
knockdown) with the per-gRNA table from 07c (the dots), converts each gRNA's
log2 fold change to percent knockdown, and draws bars = mean of that gene's own
dots with SEM whiskers.

Reads only those two CSVs, so it runs in seconds and can be re-styled without
recomputing anything.

Percent knockdown is ``100 * (1 - 2^log2fc)``, a **concave** transform. That
has two consequences this script enforces rather than assumes:

  * the bar must equal the mean of the dots actually plotted, not the percent
    knockdown of the mean log2 fold change -- both are reported so the
    difference is visible;
  * by Jensen's inequality ``mean(%KD) <= %KD(mean log2FC)``, checked below.

Four runtime guards are preserved from the published script: excluded genes
must exist, every summary gene must have per-gRNA rows, the bar must match the
mean of its dots, and the Jensen inequality must hold.
"""

import argparse
import os

import numpy as np
import pandas as pd

import figstyle as fs  # noqa: F401  -- sets the Agg backend on import
import matplotlib.pyplot as plt  # noqa: E402

from scripts.pipeline_utils import load_config  # noqa: E402

TOL = 1e-9


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, help="07d per-gene summary CSV")
    parser.add_argument("--per-guide", required=True, help="07c per-gRNA CSV")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-figure", required=True, help="Figure stem")
    parser.add_argument("--output-source-data", required=True)
    return parser.parse_args()


def pct_knockdown(log2fc):
    return 100.0 * (1.0 - 2.0**np.asarray(log2fc, dtype=float))


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["per_guide_knockdown"]
    plot_params = params["plots"]

    summary = pd.read_csv(args.summary)
    per_guide = pd.read_csv(args.per_guide)

    # --- Guard 1: excluded genes must exist -------------------------------
    exclude = list(plot_params["exclude_genes"])
    missing = [g for g in exclude if g not in set(summary["gene_symbol"])]
    if missing:
        raise ValueError(
            f"plots.exclude_genes names genes absent from the summary: {missing}. "
            "Either the gene set changed or the names are stale."
        )
    summary = summary[~summary["gene_symbol"].isin(exclude)].copy()
    print(f"Excluded {len(exclude)} genes: {', '.join(exclude)}")

    summary = summary.sort_values("log2fc_tss_vs_ntc_guides").reset_index(drop=True)
    genes = summary["gene_symbol"].tolist()
    print(f"Plotting {len(genes)} genes")

    tss = per_guide[per_guide["series"] == "TSS"].copy()
    ntc = per_guide[per_guide["series"] == "NTC"].copy()
    tss["pct_kd"] = pct_knockdown(tss["log2_fold_change"])
    ntc["pct_kd"] = pct_knockdown(ntc["log2_fold_change"])

    # --- Guard 2: every plotted gene needs dots ---------------------------
    have_dots = set(tss["gene_symbol"])
    no_dots = [g for g in genes if g not in have_dots]
    if no_dots:
        raise ValueError(f"no per-gRNA rows for: {no_dots}")

    rows = []
    for gene in genes:
        vals = tss.loc[tss["gene_symbol"] == gene, "pct_kd"].dropna().to_numpy()
        lfcs = (
            tss.loc[tss["gene_symbol"] == gene, "log2_fold_change"]
            .dropna()
            .to_numpy()
        )
        summary_row = summary[summary["gene_symbol"] == gene].iloc[0]
        mean_pct = float(vals.mean())
        pct_of_mean = float(pct_knockdown(lfcs.mean()))

        # --- Guard 4: Jensen's inequality --------------------------------
        if mean_pct > pct_of_mean + 1e-6:
            raise AssertionError(
                f"{gene}: mean(%KD)={mean_pct:.6f} exceeds "
                f"%KD(mean log2FC)={pct_of_mean:.6f}, which a concave "
                "transform cannot do -- check the inputs."
            )

        rows.append(
            {
                "gene_symbol": gene,
                "n_tss_guides_plotted": int(vals.size),
                "n_ntc_guides_plotted": int(
                    ntc.loc[ntc["gene_symbol"] == gene, "pct_kd"].notna().sum()
                ),
                "perguide_mean_pct_kd": mean_pct,
                "perguide_sem_pct_kd": (
                    float(vals.std(ddof=1) / np.sqrt(vals.size))
                    if vals.size > 1 else np.nan
                ),
                "pct_kd_of_mean_log2fc": pct_of_mean,
                "perguide_mean_log2fc": float(lfcs.mean()),
                "pct_knockdown": summary_row.get("pct_knockdown", np.nan),
                "log2fc_tss_vs_ntc_guides": summary_row["log2fc_tss_vs_ntc_guides"],
                "p_mannwhitney": summary_row["p_mannwhitney"],
                "q_bh": summary_row["q_bh"],
            }
        )
    plot_df = pd.DataFrame(rows)

    ntc_all = ntc["pct_kd"].dropna().to_numpy()
    lo_pct, hi_pct = np.percentile(ntc_all, plot_params["ntc_percentiles"])
    plot_df["ntc_p_lo_pct_kd"] = lo_pct
    plot_df["ntc_p_hi_pct_kd"] = hi_pct

    # --- Figure -----------------------------------------------------------
    fs.set_style(base_font_size=12)
    rng = np.random.RandomState(0)  # per-call reset keeps jitter deterministic

    n = len(genes)
    fig, ax = plt.subplots(figsize=(0.50 * n + 2.1, 4.9))

    ax.axhspan(
        lo_pct, hi_pct, color=fs.NTC, alpha=0.16, lw=0,
        label=f"non-targeting {plot_params['ntc_percentiles'][0]}"
              f"-{plot_params['ntc_percentiles'][1]}th pct",
    )
    ax.axhline(0, color=fs.INK_MUTED, lw=0.8, ls="--")

    heights = plot_df["perguide_mean_pct_kd"].to_numpy()
    errors = plot_df["perguide_sem_pct_kd"].to_numpy()
    ax.bar(
        range(n), heights, width=0.66, color=fs.BAR_GREY, alpha=0.55, lw=0, zorder=2,
    )
    ax.errorbar(
        range(n), heights, yerr=errors, fmt="none",
        ecolor=fs.TSS_DARK, elinewidth=1.0, capsize=2.0, zorder=4,
    )

    for i, gene in enumerate(genes):
        vals = tss.loc[tss["gene_symbol"] == gene, "pct_kd"].dropna().to_numpy()
        ax.scatter(
            np.full(vals.size, i) + rng.uniform(-0.20, 0.20, size=vals.size),
            vals, s=7, color=fs.DOT_INK, alpha=0.75, linewidth=0, zorder=5,
        )
        # --- Guard 3: bar height == mean of the plotted dots -------------
        if abs(vals.mean() - heights[i]) > TOL:
            raise AssertionError(
                f"{gene}: bar {heights[i]:.12f} != mean of dots "
                f"{vals.mean():.12f}"
            )

    span = float(np.nanmax(heights) - min(0.0, lo_pct)) or 1.0
    top = float(np.nanmax([np.nanmax(heights), hi_pct]))
    bottom = min(0.0, lo_pct, float(np.nanmin(tss["pct_kd"])))
    ax.set_ylim(bottom - 0.06 * span, top + 0.14 * span)

    for i, q in enumerate(plot_df["q_bh"]):
        marker = fs.stars(q) or "n.s."
        ax.annotate(
            marker, xy=(i, top + 0.025 * span), ha="center",
            fontsize=8, color=fs.INK,
        )

    ax.set_xticks(range(n))
    labels = ax.set_xticklabels(
        genes, rotation=90, style="italic", fontsize=plot_params["gene_label_pt"]
    )
    # --- Guard 5: label geometry is what the figure spec asks for --------
    for label in labels:
        assert label.get_rotation() == 90
        assert abs(label.get_fontsize() - plot_params["gene_label_pt"]) < 1e-9

    ax.set_ylabel("Percent knockdown", color=fs.INK_SOFT)
    ax.set_title(
        f"Per-gRNA knockdown of {n} target genes", color=fs.INK, fontsize=12
    )
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.margins(x=0.01)
    ax.annotate(
        fs.STAR_LEGEND, xy=(0.0, 1.0), xycoords="axes fraction",
        xytext=(0, 12), textcoords="offset points",
        fontsize=7, color=fs.INK_SOFT,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output_figure)), exist_ok=True)
    fs.save(fig, args.output_figure)
    print(f"Wrote {args.output_figure}.svg / .png")

    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_source_data)), exist_ok=True
    )
    plot_df.to_csv(args.output_source_data, index=False)
    print(f"Wrote {args.output_source_data}")

    sig = plot_df[plot_df["q_bh"] < plot_params["fdr_threshold"]]
    print(
        f"\n{len(sig)} / {n} genes significant at q < "
        f"{plot_params['fdr_threshold']}"
    )
    print(
        "Mean per-gRNA knockdown across plotted genes: "
        f"{plot_df['perguide_mean_pct_kd'].mean():.1f}%"
    )
    print("All guards passed.")


if __name__ == "__main__":
    main()
