#!/usr/bin/env python
"""Step 7c -- per-gRNA TSS knockdown against the non-targeting null.

Takes the two Poisson fold-change tables from 07a and, per target gene, tests
whether that gene's TSS gRNAs are shifted relative to the non-targeting gRNAs
(two-sided Mann-Whitney on log2 fold change, BH-corrected across genes).

Also exports the per-gRNA table that ``07e_plot_pct_knockdown_bar.py``
consumes as its dots -- that export, not the figure, is what the rest of the
stage depends on.

The test is on log2 fold change, deliberately. Percent knockdown is a concave
transform of it, so rank-based tests give the same answer but means do not;
keeping the test on the log scale means the p-values and the plotted bars are
computed on scales that agree.
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
    parser.add_argument("--targeting", required=True, help="07a targeting TSV")
    parser.add_argument("--non-targeting", required=True, help="07a non-targeting TSV")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-per-guide", required=True, help="Per-gRNA dot CSV")
    parser.add_argument("--output-summary", required=True, help="Per-gene summary CSV")
    parser.add_argument("--output-figure", required=True, help="Figure stem (no suffix)")
    return parser.parse_args()


def load_panel(config):
    """Panel gene symbols, and the subset flagged as expressed in the cell type."""
    params = config["per_guide_knockdown"]
    table = pd.read_csv(params["chd_gene_table"], sep="\t", dtype=str)
    name_col = params["chd_gene_column"]
    flag_col = params["vic_expressed_column"]
    if name_col not in table.columns:
        raise ValueError(f"{name_col!r} not in {params['chd_gene_table']}")
    genes = table[name_col].str.strip()
    if flag_col in table.columns:
        flagged = table[flag_col].str.strip() == params["vic_expressed_value"]
        focus = genes[flagged].tolist()
    else:
        focus = genes.tolist()
    return genes.tolist(), focus


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["per_guide_knockdown"]
    plot_params = params["plots"]

    all_genes, focus_genes = load_panel(config)
    print(f"Panel: {len(all_genes)} genes, {len(focus_genes)} flagged as expressed")

    # --- TSS gRNAs: on-target rows only -----------------------------------
    targeting = pd.read_csv(args.targeting, sep="\t")
    on_target = targeting[
        targeting["is_pos_ctrl"].astype(bool)
        & (targeting["grna_target"] == targeting["response_id"])
    ].copy()
    on_target["series"] = "TSS"
    print(f"On-target TSS rows: {len(on_target)}")

    nt = pd.read_csv(args.non_targeting, sep="\t").copy()
    nt["series"] = "NTC"
    print(f"Non-targeting rows: {len(nt)}")

    combined = pd.concat([on_target, nt], ignore_index=True)
    combined = combined[combined["gene_symbol"].isin(all_genes)]
    combined = combined.dropna(subset=["log2_fold_change"])
    combined["in_main_figure"] = combined["gene_symbol"].isin(focus_genes)

    # --- Per-gene Mann-Whitney on log2 fold change ------------------------
    rows = []
    for gene in all_genes:
        block = combined[combined["gene_symbol"] == gene]
        tss = block.loc[block["series"] == "TSS", "log2_fold_change"].to_numpy()
        ntc = block.loc[block["series"] == "NTC", "log2_fold_change"].to_numpy()
        if tss.size >= 2 and ntc.size >= 2:
            p = stats.mannwhitneyu(tss, ntc, alternative="two-sided").pvalue
        else:
            p = np.nan
        mean_lfc = float(np.mean(tss)) if tss.size else np.nan
        rows.append(
            {
                "gene_symbol": gene,
                "n_tss_guides": int(tss.size),
                "n_ntc_guides": int(ntc.size),
                "tss_mean_log2fc": mean_lfc,
                "ntc_mean_log2fc": float(np.mean(ntc)) if ntc.size else np.nan,
                # Convert once, from the mean log2FC -- never the mean of
                # per-guide percentages (concave transform).
                "pct_knockdown_of_mean_log2fc": (
                    100 * (1 - 2**mean_lfc) if np.isfinite(mean_lfc) else np.nan
                ),
                "p_mannwhitney": p,
                "in_main_figure": gene in focus_genes,
            }
        )
    summary = pd.DataFrame(rows)
    summary["q_bh"] = benjamini_hochberg(summary["p_mannwhitney"])

    n_sig = int((summary["q_bh"] < plot_params["fdr_threshold"]).sum())
    print(
        f"Significant at q < {plot_params['fdr_threshold']}: "
        f"{n_sig} / {summary['p_mannwhitney'].notna().sum()} tested"
    )

    # --- Figure: per-gene dot strips against the NT band ------------------
    fs.set_style(plot_params["base_font_size"])
    rng = np.random.default_rng(plot_params["random_seed"])

    plotted = summary[summary["in_main_figure"]].copy()
    plotted = plotted.sort_values("tss_mean_log2fc")
    genes = plotted["gene_symbol"].tolist()

    fig, ax = plt.subplots(figsize=(max(6.0, 0.30 * len(genes) + 2.0), 4.2))

    ntc_all = combined.loc[combined["series"] == "NTC", "log2_fold_change"].to_numpy()
    lo, hi = np.percentile(ntc_all, plot_params["ntc_percentiles"])
    ax.axhspan(
        lo, hi, color=fs.NTC, alpha=0.18, lw=0,
        label=f"non-targeting {plot_params['ntc_percentiles'][0]}"
              f"-{plot_params['ntc_percentiles'][1]}th pct",
    )
    ax.axhline(0, color=fs.INK_MUTED, lw=0.8, ls="--")
    for fraction in plot_params["kd_reference_fractions"]:
        y = np.log2(fraction)
        ax.axhline(y, color=fs.INK_MUTED, lw=0.6, ls=":")
        ax.annotate(
            f"{100 * (1 - fraction):.0f}% KD",
            xy=(1.0, y), xycoords=("axes fraction", "data"),
            xytext=(-2, 2), textcoords="offset points",
            ha="right", fontsize=6.5, color=fs.INK_SOFT,
        )

    for i, gene in enumerate(genes):
        vals = combined[
            (combined["gene_symbol"] == gene) & (combined["series"] == "TSS")
        ]["log2_fold_change"].to_numpy()
        if not vals.size:
            continue
        jitter = rng.uniform(-0.20, 0.20, size=vals.size)
        ax.scatter(
            np.full(vals.size, i) + jitter, vals,
            s=9, color=fs.TSS, edgecolor=fs.SURFACE, linewidth=0.4, zorder=3,
        )
        ax.plot([i - 0.30, i + 0.30], [vals.mean()] * 2,
                color=fs.TSS_DARK, lw=1.4, zorder=4)

    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, style="italic", fontsize=7)
    ax.set_ylabel("log$_2$ fold change vs null model", color=fs.INK_SOFT)
    ax.set_xlabel("")
    ax.set_title(
        f"Per-gRNA on-target knockdown, {len(genes)} target genes",
        color=fs.INK, fontsize=10,
    )
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    ax.margins(x=0.01)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_figure)), exist_ok=True)
    fs.save(fig, args.output_figure)
    print(f"Wrote {args.output_figure}.svg / .png")

    # --- Exports ----------------------------------------------------------
    for path in (args.output_per_guide, args.output_summary):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    per_guide_cols = [
        "grna_id", "grna_target", "response_id", "gene_symbol",
        "n_trt", "log2_fold_change", "series", "in_main_figure",
    ]
    combined[per_guide_cols].to_csv(args.output_per_guide, index=False)
    summary.to_csv(args.output_summary, index=False)
    print(f"Wrote {args.output_per_guide} ({len(combined)} rows)")
    print(f"Wrote {args.output_summary} ({len(summary)} rows)")


if __name__ == "__main__":
    main()
