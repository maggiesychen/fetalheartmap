"""Shared figure style for the CHD-gene TSS knockdown figure.

Two series only, and each colour does exactly one job:

  TSS      -> blue  #2a78d6  : cells carrying a TSS-targeting gRNA (the perturbation)
  NTC      -> grey  #6b6a66  : cells carrying a non-targeting gRNA (the reference)

Colour choice was validated with the data-viz palette validator
(`validate_palette.js "#2a78d6,#6b6a66" --mode light`): lightness band PASS,
CVD separation PASS (worst adjacent dE 17.6 protan / 12.5 tritan, target >=8),
normal-vision floor PASS (17.6, floor 15), contrast vs surface PASS (both >=3:1).
The one reported FAIL is the *chroma floor* on the grey - i.e. "this reads as
grey". That is deliberate and not a defect here: NTC is a neutral reference
series, not a peer category, so it should recede. Identity is never carried by
colour alone - the NTC bar is directly labelled on the axis and in the legend,
and it is separated from the gene bars by a gap in the category axis.

No value-ramp is applied across genes: every gene bar is the same blue, because
gene identity is nominal and bar length already encodes the effect size.

The combined all-genes bar figure (05_...) is a later, stripped-down panel that
plots only the TSS series and so uses the single-series greys below instead of
the blue/grey pair; 02_, 04_ and 06_ still use the two-series colours.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- series colours ---------------------------------------------------------
TSS = "#2a78d6"          # TSS-targeting gRNAs
TSS_DARK = "#1c5296"     # dot / errorbar ink for the TSS series
NTC = "#6b6a66"          # non-targeting gRNAs (neutral reference)
NTC_DARK = "#46453f"

# --- single-series palette (05_plot_log2fc_combined_bar.py) ------------------
# That figure drops the NTC series entirely, so there is no two-series encoding
# left to preserve and no colour carries identity: one grey bar per gene, black
# dots for the individual gRNAs. The grey is lighter than the NTC grey above
# because the dots sit INSIDE the bars, so bar-to-dot contrast is what has to
# hold: black on BAR_GREY is 5.7:1, on the NTC grey it would be only 3.5:1.
# BAR_GREY still clears 3:1 against the white surface (3.3:1).
BAR_GREY = "#8E8E8A"
DOT_INK = "#111111"

# --- ink / chrome ----------------------------------------------------------
INK = "#111111"
INK_SOFT = "#5A5F66"
INK_MUTED = "#8A9099"
SURFACE = "#FFFFFF"      # the 2px "surface gap/ring" colour around marks
GRID = "#E3E4E1"


def set_style(base_font_size=8):
    """Paper-ready matplotlib defaults.

    svg.fonttype='none' keeps text as real <text> nodes so labels stay editable
    in Illustrator during figure assembly - which is why this folder ships SVG
    rather than PDF.
    """
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": base_font_size,
        "axes.titlesize": base_font_size + 1,
        "axes.labelsize": base_font_size,
        "xtick.labelsize": base_font_size - 1,
        "ytick.labelsize": base_font_size - 1,
        "legend.fontsize": base_font_size - 1,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK_SOFT,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK_SOFT,
        "ytick.color": INK_SOFT,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "text.color": INK,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.transparent": False,
    })


def save(fig, out_stem):
    """PNG for on-screen checking, SVG (editable text) for figure assembly."""
    fig.savefig(f"{out_stem}.png", dpi=400)
    fig.savefig(f"{out_stem}.svg")
    plt.close(fig)
    print(f"  saved {out_stem}.png / .svg", flush=True)


# p-value -> asterisks (same ladder used elsewhere in this project).
STAR_LADDER = ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*"))
STAR_LEGEND = "* q<0.05   ** q<0.01   *** q<0.001   **** q<1e-4"


def stars(q):
    """BH-adjusted p -> asterisk string ('' if not significant / missing)."""
    try:
        q = float(q)
    except (TypeError, ValueError):
        return ""
    if q != q:  # NaN
        return ""
    for cut, s in STAR_LADDER:
        if q < cut:
            return s
    return "ns"
