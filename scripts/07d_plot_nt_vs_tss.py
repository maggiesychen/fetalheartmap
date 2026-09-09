#!/usr/bin/env python
"""Step 7d -- non-targeting vs TSS pseudobulk expression, per target gene.

The model-free counterpart to 07c. Tests each gene's TSS gRNAs against the
non-targeting gRNAs on **CP10K** (not log2 fold change), so the result depends
on the count data and the assignment matrix but on no model at all.

Also cross-checks the two estimators against each other and reports the
correlation between the CP10K-derived and Poisson-derived per-gene fold
changes. Disagreement here would mean one of them is picking up a modelling
artefact.

Exports the per-gene summary that ``07e_plot_pct_knockdown_bar.py`` uses for
its asterisks, its column ordering and its headline knockdown value. Note that
BH correction is applied over whichever genes are in *this* table -- changing
the gene set changes the multiple-testing universe and therefore the stars.
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

import figstyle as fs  # noqa: F401  -- sets the Agg backend on import
import matplotlib.pyplot as plt  # noqa: E402

from scripts.pipeline_utils import benjamini_hochberg, load_config  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True, help="07b per-guide TSV")
    parser.add_argument("--targeting", required=True, help="07a targeting TSV")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-summary", required=True, help="Per-gene summary CSV")
    parser.add_argument("--output-figure", required=True, help="Figure stem")
    parser.add_argument(
        "--value-column",
        default="cp10k",
        help="Expression column to plot and test [default: %(default)s]",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["per_guide_knockdown"]
    plot_params = params["plots"]
    value = args.value_column

    expr = pd.read_csv(args.expression, sep="\t")
    if value not in expr.columns:
        raise ValueError(f"{value!r} not in {args.expression}")
    print(f"Expression rows: {len(expr)}  (testing on {value})")

    # Poisson per-gene mean, for the cross-check column only.
    targeting = pd.read_csv(args.targeting, sep="\t")
    on_target = targeting[
        targeting["is_pos_ctrl"].astype(bool)
        & (targeting["grna_target"] == targeting["response_id"])
    ]
    poisson_mean = (
        on_target.groupby("gene_symbol")["log2_fold_change"].mean().to_dict()
    )

    genes = sorted(expr["gene_symbol"].dropna().unique())
    rows = []
    for gene in genes:
        block = expr[expr["gene_symbol"] == gene]
        tss = block.loc[block["series"] == "TSS", value].to_numpy()
        ntc = block.loc[block["series"] == "NTC", value].to_numpy()
        pool = block["nt_pool_cp10k"].dropna()
        pool_value = float(pool.iloc[0]) if len(pool) else np.nan

        if tss.size >= 2 and ntc.size >= 2:
            p = stats.mannwhitneyu(tss, ntc, alternative="two-sided").pvalue
        else:
            p = np.nan

        tss_mean = float(np.mean(tss)) if tss.size else np.nan
        ntc_mean = float(np.mean(ntc)) if ntc.size else np.nan
        log2fc = (
            float(np.log2(tss_mean / ntc_mean))
            if np.isfinite(tss_mean) and np.isfinite(ntc_mean) and ntc_mean > 0
            else np.nan
        )
        rows.append(
            {
                "gene_symbol": gene,
                "n_tss_guides": int(tss.size),
                "n_ntc_guides": int(ntc.size),
                "nt_pool_cp10k": pool_value,
                "ntc_guide_mean_cp10k": ntc_mean,
                "tss_guide_mean_cp10k": tss_mean,
                "log2fc_tss_vs_ntc_guides": log2fc,
                "pct_knockdown": (
                    100 * (1 - tss_mean / ntc_mean)
                    if np.isfinite(tss_mean) and np.isfinite(ntc_mean) and ntc_mean > 0
                    else np.nan
                ),
                "p_mannwhitney": p,
                "poisson_mean_log2fc": poisson_mean.get(gene, np.nan),
            }
        )

    summary = pd.DataFrame(rows)
    summary["q_bh"] = benjamini_hochberg(summary["p_mannwhitney"])

    # --- Cross-check the two independent estimators -----------------------
    both = summary.dropna(subset=["log2fc_tss_vs_ntc_guides", "poisson_mean_log2fc"])
    if len(both) >= 3:
        r = stats.pearsonr(
            both["log2fc_tss_vs_ntc_guides"], both["poisson_mean_log2fc"]
        )
        print(
            f"\nCP10K vs Poisson log2FC across {len(both)} genes: "
            f"Pearson r = {r.statistic:.3f} (p = {r.pvalue:.3g})"
        )
        if r.statistic < 0.8:
            print(
                "  WARNING: the two estimators disagree (r < 0.8). One of them "
                "is likely picking up an artefact; investigate before using "
                "either for a figure."
            )
    n_sig = int((summary["q_bh"] < plot_params["fdr_threshold"]).sum())
    print(
        f"Significant at q < {plot_params['fdr_threshold']}: {n_sig} / "
        f"{int(summary['p_mannwhitney'].notna().sum())} tested"
    )

    # --- Small multiples --------------------------------------------------
    fs.set_style(plot_params["base_font_size"])
    rng = np.random.default_rng(plot_params["random_seed"])

    ordered = summary.sort_values("log2fc_tss_vs_ntc_guides")["gene_symbol"].tolist()
    ncol = 5
    nrow = int(np.ceil(len(ordered) / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(7.2, 1.42 * nrow + 0.70), squeeze=False
    )

    for ax_i, gene in enumerate(ordered):
        ax = axes[ax_i // ncol][ax_i % ncol]
        block = expr[expr["gene_symbol"] == gene]
        series = [
            ("NTC", block.loc[block["series"] == "NTC", value].to_numpy(),
             fs.NTC, fs.NTC_DARK, 3.0, 0.30),
            ("TSS", block.loc[block["series"] == "TSS", value].to_numpy(),
             fs.TSS, fs.TSS_DARK, 9.0, 0.80),
        ]
        for x, (label, vals, colour, dark, size, alpha) in enumerate(series):
            if not vals.size:
                continue
            ax.bar(x, vals.mean(), width=0.64, color=colour, alpha=0.35, lw=0)
            ax.scatter(
                np.full(vals.size, x) + rng.uniform(-0.19, 0.19, size=vals.size),
                vals, s=size, color=dark, alpha=alpha, linewidth=0, zorder=3,
            )
        row = summary[summary["gene_symbol"] == gene].iloc[0]
        top = max(
            (s[1].max() for s in series if s[1].size), default=1.0
        )
        ax.set_ylim(0, top * 1.40)
        marker = fs.stars(row["q_bh"]) or "n.s."
        ax.annotate(
            marker, xy=(0.5, top * 1.13), ha="center", fontsize=7, color=fs.INK,
        )
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["NT", "TSS"], fontsize=6)
        ax.set_title(gene, fontsize=7, style="italic", color=fs.INK)
        ax.tick_params(labelsize=6)

    for ax_i in range(len(ordered), nrow * ncol):
        axes[ax_i // ncol][ax_i % ncol].axis("off")

    fig.supylabel(f"{value} (pseudobulk)", fontsize=8, color=fs.INK_SOFT)
    fig.suptitle(
        f"Target-gene expression, non-targeting vs TSS gRNAs ({len(ordered)} genes)",
        fontsize=10, color=fs.INK,
    )
    fig.tight_layout(rect=(0.02, 0, 1, 0.97))

    os.makedirs(os.path.dirname(os.path.abspath(args.output_figure)), exist_ok=True)
    fs.save(fig, args.output_figure)
    print(f"\nWrote {args.output_figure}.svg / .png")

    os.makedirs(os.path.dirname(os.path.abspath(args.output_summary)), exist_ok=True)
    summary.to_csv(args.output_summary, index=False)
    print(f"Wrote {args.output_summary} ({len(summary)} genes)")


if __name__ == "__main__":
    main()
