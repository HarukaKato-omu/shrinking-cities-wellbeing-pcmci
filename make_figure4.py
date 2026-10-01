#!/usr/bin/env python3
"""Figure 4 of the manuscript: dose-response of the within-person
Caregiving Support (t-1) -> Work-Life Balance (t) dependency among older
adults by the municipal annual population change rate.

Reads, from the output directory of validation_analysis.py (steps `dose`
or `dose_extra`):
  V9_figure4_lines.json            fitted coefficient lines (pooled 65+ and 15-64)
  V9_figure4_quintile_slopes.csv   within OLS slopes by rate quintile, 65+ by scale
  V9_dose_response.csv             interaction estimates and cluster p-values (annotation)

Usage:  python make_figure4.py --indir validation --outdir figures
Style matches the manuscript figures (grey plot surface, serif font).
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

COL_RED, COL_BLUE, COL_GREY, SURF = "#c62828", "#1565c0", "#9e9e9e", "#e6e6e6"
THRESH = -0.15                                     # SCiRN criterion, %/yr


def set_font():
    """TeX Gyre Termes if installed, else Times New Roman, else the default serif."""
    for f in font_manager.findSystemFonts():
        if "termes" in Path(f).name.lower():
            font_manager.fontManager.addfont(f)
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("TeX Gyre Termes", "Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"):
        if name in available:
            plt.rcParams.update({"font.family": name, "mathtext.fontset": "custom", "mathtext.rm": name,
                                 "mathtext.it": f"{name}:italic", "mathtext.bf": f"{name}:bold"})
            return name
    plt.rcParams["font.family"] = "serif"
    return "serif"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True, help="output directory of validation_analysis.py")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()
    indir, outdir = Path(args.indir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    set_font()

    lines = json.loads((indir / "V9_figure4_lines.json").read_text(encoding="utf-8"))
    bins = pd.read_csv(indir / "V9_figure4_quintile_slopes.csv")
    dose = pd.read_csv(indir / "V9_dose_response.csv")
    ols = dose[dose.estimator.str.startswith("within OLS")].set_index(["scale", "cohort"])
    p65, p1564 = ols.loc[("pooled", "65+"), "p_cluster_t"], ols.loc[("pooled", "15-64"), "p_cluster_t"]

    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    fig.subplots_adjust(left=0.11, right=0.985, top=0.96, bottom=0.17)
    ax.set_facecolor(SURF)
    for sp in ax.spines.values():
        sp.set_color("black"); sp.set_linewidth(0.8)
    ax.axvline(THRESH, color="#555555", lw=0.8, ls=(0, (4, 3)))
    ax.axhline(0, color="#888888", lw=0.8, ls=":")
    ax.text(THRESH - 0.04, 0.315, "shrinking (SCiRN)  ", ha="right", va="top", fontsize=8, color="#444444")
    ax.text(THRESH + 0.04, 0.315, "  non-shrinking", ha="left", va="top", fontsize=8, color="#444444")

    L = lines["65+"]; g = np.asarray(L["grid"]); s = np.asarray(L["slope"]); se = np.asarray(L["se"])
    lo, hi = L["r_range"]; m = (g >= lo) & (g <= hi)
    ax.fill_between(g[m], (s - 1.96 * se)[m], (s + 1.96 * se)[m], color="#4a4a4a", alpha=0.13, lw=0)
    ax.plot(g[m], s[m], color="#4a4a4a", lw=2.0)
    L2 = lines["15-64"]; g2 = np.asarray(L2["grid"]); s2 = np.asarray(L2["slope"])
    m2 = (g2 >= L2["r_range"][0]) & (g2 <= L2["r_range"][1])
    ax.plot(g2[m2], s2[m2], color=COL_GREY, lw=1.6, ls=(0, (5, 2)))

    for scale, col, mk, dx in (("SMC", COL_RED, "o", 0.0), ("Large", COL_BLUE, "s", 0.03)):
        d = bins[bins.scale == scale]
        for r in d.itertuples():
            ax.plot([r.rate_median + dx] * 2, [r.slope_CS - 1.96 * r.se_cluster, r.slope_CS + 1.96 * r.se_cluster],
                    color=col, lw=1.2, zorder=3)
        ax.plot(d.rate_median + dx, d.slope_CS, mk, color=col, markersize=6.5, markeredgecolor=col,
                markerfacecolor=col, zorder=4)

    ax.set_xlim(-2.45, 1.05); ax.set_ylim(-0.33, 0.33)
    ax.set_xticks([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]); ax.set_yticks([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3])
    ax.tick_params(labelsize=9)
    ax.set_xlabel(r"$\bf{Municipal\ annual\ population\ change\ rate}$ (%/yr, 2017–2023)", fontsize=10)
    ax.set_ylabel(r"$\bf{Within}$-$\bf{person\ coefficient}$:" + "\nCaregiving Support (t-1) → Work-Life Balance (t)", fontsize=9.5)
    ax.text(0.985, 0.04,
            f"65+ (both scales): Caregiving Support (t-1) × rate = {L['b_int']:.3f} (SE {L['se_int']:.3f}), cluster $p$ = {p65:.3f}\n"
            f"15-64 (age-pattern check): {L2['b_int']:.3f} (SE {L2['se_int']:.3f}), $p$ = {p1564:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.2, color="#222222",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=2.5))
    handles = [Line2D([], [], marker="o", color=COL_RED, lw=0, markerfacecolor=COL_RED, markersize=6, label="SMCs, 65+ (rate quintiles)"),
               Line2D([], [], marker="s", color=COL_BLUE, lw=0, markerfacecolor=COL_BLUE, markersize=6, label="Large cities, 65+ (rate quintiles)"),
               Line2D([], [], color="#4a4a4a", lw=2.0, label="Fitted coefficient, 65+ (both scales) with 95% CI"),
               Line2D([], [], color=COL_GREY, lw=1.6, ls=(0, (5, 2)), label="Fitted coefficient, 15-64 (age-pattern check)")]
    ax.legend(handles=handles, loc="upper left", fontsize=8.2, frameon=False, handlelength=2.2, borderaxespad=0.4)
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"figure4.{ext}", dpi=300)
    print(f"written {outdir / 'figure4.pdf'} and figure4.png")


if __name__ == "__main__":
    main()
