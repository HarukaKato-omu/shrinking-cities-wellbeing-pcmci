#!/usr/bin/env python3
"""Figures 2-4 of the article, drawn from aggregate outputs only (repository release v2.1.1).

  Figure 2   pooled lagged dependency networks, eight strata             <- network run
  Figure 3   person-demeaned (within-person) networks in SMCs, by age    <- network run
  Figure 4   Caregiving Support (t-1) -> 13 target domains (t), older adults in SMCs, by the pace of
             population decline: (a) coefficient at the 10th and 90th percentile of the rate,
             (b) interaction with the rate, (c) coefficient along the rate for Work-Life Balance and
             Jobs & Wages                               <- care_links.py

Usage
  python make_figures.py --network-dir NET [--strain-dir STR --care-dir CARE] --outdir figures

The network run is validation_analysis.py --window strict3 --design mci (file N_all_edges.csv).  Figure 4
needs CARE/V10_family_SMC_lag1.csv and CARE/V11_quintile_slopes.csv; STR/V9_figure4_lines.json, if present,
supplies the range of the rate over which the fitted lines are drawn (its central 98%).
Node positions, the out-strength colour scale (0-1.2) and the edge-width scale are the same in every network
panel and are those of the earlier releases (ring order fixed in MANUSCRIPT_RING; the 2-opt routine that
produced it is kept in ring_order_2opt and can be re-run with --layout refit).  Reads no person-level data.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors, font_manager, gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch

MAKE_FIGURES_VERSION = "2.1.3"
HUB, CS, WLB, JW = "Life Satisfaction", "Caregiving Support", "Work-Life Balance", "Jobs & Wages"
DOMAINS = ["Life Satisfaction", "Finances & Assets", "Jobs & Wages", "Housing", "Work-Life Balance", "Health", "Education",
           "Friends & Community", "Trust in Government & Courts", "Air & Water (environment)", "Personal Safety",
           "Parenting Support", "Caregiving Support", "Enjoyment of Life"]
MANUSCRIPT_RING = ["Jobs & Wages", "Caregiving Support", "Friends & Community", "Health", "Enjoyment of Life", "Education",
                   "Housing", "Air & Water (environment)", "Personal Safety", "Parenting Support",
                   "Trust in Government & Courts", "Finances & Assets", "Work-Life Balance"]
# label text, (dx, dy) from the node centre, horizontal and vertical alignment
LABELS = {
    "Life Satisfaction": ("Life\nSatisfaction", (0.0, -0.15), "center", "top"),
    "Jobs & Wages": ("Jobs & Wages", (0.0, 0.13), "center", "bottom"),
    "Caregiving Support": ("Caregiving Support", (0.11, 0.10), "left", "bottom"),
    "Friends & Community": ("Friends &\nCommunity", (0.14, 0.02), "left", "center"),
    "Health": ("Health", (0.16, 0.0), "left", "center"),
    "Enjoyment of Life": ("Enjoyment of Life", (0.30, -0.12), "center", "top"),
    "Education": ("Education", (0.08, -0.11), "left", "top"),
    "Housing": ("Housing", (0.08, -0.12), "center", "top"),
    "Air & Water (environment)": ("Air & Water (environment)", (0.10, -0.12), "right", "top"),
    "Personal Safety": ("Personal Safety", (-0.08, -0.11), "right", "top"),
    "Parenting Support": ("Parenting\nSupport", (-0.14, -0.07), "right", "center"),
    "Trust in Government & Courts": ("Trust in\nGovernment\n& Courts", (-0.125, 0.0), "right", "center"),
    "Finances & Assets": ("Finances &\nAssets", (-0.12, 0.03), "right", "center"),
    "Work-Life Balance": ("Work–Life Balance", (-0.08, 0.11), "right", "bottom"),
}
SCALE_MAX, WIDTH_PER_MCI = 1.2, 17.0
NODE_R0, NODE_R1 = 0.085, 0.075
CURV_RECIP, CURV_AVOID, AVOID_DIST = 0.18, 0.16, 0.17
CMAP, NORM = cm.Oranges, colors.Normalize(vmin=0.0, vmax=SCALE_MAX)
COL_POS, COL_NEG, COL_EMPTY, BAND = "#808080", "#c0392b", "#f5f5f5", "#231815"
PANEL_W, PANEL_H, LEGEND_H = 3.8, 3.14, 0.78
RED, GREY, INK, SURF = "#c62828", "#555555", "#222222", "#f2f2f2"
THRESH = -0.15                                     # SCiRN criterion, %/yr


def label(prefix, scale, cohort):
    return f"{prefix} {scale}, aged {cohort}"


FIG2 = [("(a1)", label("All", "Large cities", "15-64")), ("(a2)", label("All", "Large cities", "65+")),
        ("(b1)", label("Shrinking", "Large cities", "15-64")), ("(b2)", label("Shrinking", "Large cities", "65+")),
        ("(c1)", label("All", "SMCs", "15-64")), ("(c2)", label("All", "SMCs", "65+")),
        ("(d1)", label("Shrinking", "SMCs", "15-64")), ("(d2)", label("Shrinking", "SMCs", "65+"))]
FIG3 = [("(a)", label("All", "SMCs", "15-64")), ("(b)", label("All", "SMCs", "65+")),
        ("(c)", label("Shrinking", "SMCs", "15-64")), ("(d)", label("Shrinking", "SMCs", "65+"))]


def set_font():
    """TeX Gyre Termes if installed, else Times New Roman, else the nearest serif."""
    for f in font_manager.findSystemFonts():
        if "termes" in Path(f).name.lower():
            font_manager.fontManager.addfont(f)
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("TeX Gyre Termes", "Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"):
        if name in available:
            return {"font.family": name, "mathtext.fontset": "custom", "mathtext.rm": name,
                    "mathtext.it": f"{name}:italic", "mathtext.bf": f"{name}:bold"}
    return {"font.family": "serif"}


# ---------------------------------------------------------------------------------------------------------
# networks (Figures 2 and 3)
# ---------------------------------------------------------------------------------------------------------
def edge_dict(path, mode):
    """Reported edges (q < .05 and |MCI| >= .10, as flagged by validation_analysis.py) per stratum."""
    e = pd.read_csv(path, encoding="utf-8-sig")
    e = e[(e["mode"] == mode) & e["reported"].astype(bool)]
    return {m: [(r.source, r.target, float(r.MCI)) for r in g.itertuples()] for m, g in e.groupby("model")}


def ring_order_2opt(pair_w, initial):
    """Deterministic 2-opt descent minimising |MCI|-weighted angular distance (the ring layout routine)."""
    order, n = list(initial), len(initial)

    def cost(o):
        pos = {node: i for i, node in enumerate(o)}
        return sum(w * min(abs(pos[a] - pos[b]), n - abs(pos[a] - pos[b])) for (a, b), w in pair_w.items())

    improved = True
    while improved:
        improved, base = False, cost(order)
        for i in range(n - 1):
            for j in range(i + 1, n):
                cand = list(order)
                cand[i], cand[j] = cand[j], cand[i]
                c = cost(cand)
                if c < base - 1e-12:
                    order, base, improved = cand, c, True
    return order


def refit_ring(edge_sets):
    pair_w = {}
    for edges in edge_sets:
        for u, v, mci in edges:
            if HUB not in (u, v):
                k = tuple(sorted((u, v)))
                pair_w[k] = max(pair_w.get(k, 0.0), abs(mci))
    ring = ring_order_2opt(pair_w, [d for d in DOMAINS if d != HUB])
    print("[layout] ring order:", " | ".join(ring))
    return ring


def ring_positions(ring):
    pos = {HUB: np.array([0.0, 0.0])}
    for i, node in enumerate(ring):
        ang = np.deg2rad(90.0 - i * 360.0 / len(ring))
        pos[node] = np.array([np.cos(ang), np.sin(ang)])
    return pos


def out_strength(edges):
    s = {}
    for u, _v, mci in edges:
        s[u] = s.get(u, 0.0) + abs(mci)
    return s


def node_radius(s):
    return NODE_R0 + NODE_R1 * (min(s, SCALE_MAX) / SCALE_MAX)


def edge_rad(u, v, pairset, pos):
    if (v, u) in pairset:
        return CURV_RECIP
    p, q = pos[u], pos[v]
    d = q - p
    L = float(np.linalg.norm(d))
    for w, xy in pos.items():
        if w in (u, v):
            continue
        t = float(np.dot(xy - p, d)) / (L * L)
        if 0.05 < t < 0.95 and np.linalg.norm(p + t * d - xy) < AVOID_DIST:
            return CURV_AVOID
    return 0.0


def draw_panel(ax, edges, title, pos):
    ax.set_xlim(-1.56, 1.64); ax.set_ylim(-1.30, 1.34)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("black"); sp.set_linewidth(0.8)
    ax.set_title(title, fontsize=10.5, fontweight="bold", loc="left", pad=3)
    strengths, patches = out_strength(edges), {}
    for node, xy in pos.items():
        s = strengths.get(node, 0.0)
        c = Circle(xy, node_radius(s), facecolor=CMAP(NORM(s)) if s > 0 else COL_EMPTY, edgecolor="#222222", linewidth=0.8, zorder=3)
        ax.add_patch(c); patches[node] = c
        text, (dx, dy), ha, va = LABELS[node]
        ax.text(xy[0] + dx, xy[1] + dy, text, ha=ha, va=va, fontsize=8.6, linespacing=0.95, zorder=5)
    pairset = {(u, v) for u, v, _ in edges}
    for u, v, mci in sorted(edges, key=lambda e: abs(e[2])):
        neg = mci < 0
        ax.add_patch(FancyArrowPatch(pos[u], pos[v], patchA=patches[u], patchB=patches[v],
                                     connectionstyle=f"arc3,rad={edge_rad(u, v, pairset, pos)}", arrowstyle="-|>", mutation_scale=9,
                                     linewidth=abs(mci) * WIDTH_PER_MCI, color=COL_NEG if neg else COL_POS,
                                     linestyle=(0, (4, 2)) if neg else "solid", shrinkA=0, shrinkB=0, zorder=2))


def legend_band(fig, fig_h):
    h = LEGEND_H / fig_h
    band = fig.add_axes([0.0, 0.0, 1.0, h]); band.set_facecolor(BAND)
    band.set_xlim(0, 1); band.set_ylim(0, 1); band.set_xticks([]); band.set_yticks([])
    for sp in band.spines.values():
        sp.set_visible(False)
    band.text(0.012, 0.86, "Legend", color="white", fontsize=11, va="top", ha="left")
    band.text(0.26, 0.86, "Out-Strength centrality", color="white", fontsize=11, va="top", ha="center")
    band.text(0.985, 0.86, "Lagged dependency strength", color="white", fontsize=11, va="top", ha="right")
    for y, lw, lab, col, ls in ((0.48, 0.10 * WIDTH_PER_MCI, "|MCI|: 0.10", "white", "solid"),
                                (0.28, 0.15 * WIDTH_PER_MCI, "|MCI|: 0.15", "white", "solid"),
                                (0.08, 0.13 * WIDTH_PER_MCI, "Negative", COL_NEG, (0, (4, 2)))):
        band.plot([0.805, 0.875], [y, y], color=col, lw=lw, ls=ls, solid_capstyle="butt")
        band.text(0.885, y, lab, color="white", fontsize=10.5, va="center", ha="left")
    cax = fig.add_axes([0.18, 0.30 * h, 0.535, 0.24 * h])
    cb = fig.colorbar(cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    cb.ax.tick_params(labelsize=10, colors="white", length=0, pad=2)
    cb.outline.set_edgecolor("white"); cb.outline.set_linewidth(1.2)


def network_figure(edges_by_model, panels, pos, fname, outdir):
    nrow = (len(panels) + 1) // 2
    fig_h = nrow * PANEL_H + LEGEND_H
    fig = plt.figure(figsize=(2 * PANEL_W, fig_h))
    gs = gridspec.GridSpec(nrow, 2, left=0.002, right=0.998, top=1.0 - 0.20 / fig_h, bottom=(LEGEND_H + 0.03) / fig_h, wspace=0.0, hspace=0.075)
    for k, (tag, model) in enumerate(panels):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        strat, age = model.rsplit(", aged ", 1)
        draw_panel(ax, edges_by_model.get(model, []), f"{tag} {strat} (Respondents aged {age})", pos)
    legend_band(fig, fig_h)
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{fname}.{ext}", dpi=300)
    plt.close(fig)
    print(f"[figure] {fname} written ({nrow}x2 panels)")


# ---------------------------------------------------------------------------------------------------------
# Figure 4
# ---------------------------------------------------------------------------------------------------------
def fitted_line(row, grid):
    """Implied coefficient b_CS + b_int * r and its standard error along the rate, recovered from the
    estimates at the 10th and 90th percentile and the interaction (the variance is quadratic in r)."""
    b3, v33 = row.interaction, row.se_interaction ** 2
    r10, r90, v10, v90 = row.rate_p10, row.rate_p90, row.se_at_p10 ** 2, row.se_at_p90 ** 2
    v13 = ((v10 - v90) - (r10 ** 2 - r90 ** 2) * v33) / (2.0 * (r10 - r90))
    v11 = v10 - 2.0 * r10 * v13 - r10 ** 2 * v33
    return row.coef_at_p10 + b3 * (grid - r10), np.sqrt(np.maximum(v11 + 2.0 * grid * v13 + grid ** 2 * v33, 0.0))


def minus(s):
    return s.replace("-", "−")


def figure4(fam_csv, quint_csv, lines_json, outdir, fname="figure4"):
    tab = pd.read_csv(fam_csv, encoding="utf-8-sig").sort_values("interaction", ascending=False).reset_index(drop=True)
    quint = pd.read_csv(quint_csv, encoding="utf-8-sig"); quint = quint[quint.scale == "SMC"]
    k = len(tab); ypos = np.arange(k)[::-1]; r0 = tab.iloc[0]
    fig = plt.figure(figsize=(7.4, 8.15))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.74, 1.0], left=0.262, right=0.945, top=0.972, bottom=0.052, hspace=0.30, wspace=0.06)
    ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1], sharey=ax1); ax3 = fig.add_subplot(gs[1, :])
    for ax in (ax1, ax2, ax3):
        ax.set_facecolor(SURF)
        for sp in ax.spines.values():
            sp.set_color("black"); sp.set_linewidth(1.0)
        ax.tick_params(labelsize=10)
    eb = dict(capsize=3.5, capthick=1.3, elinewidth=1.3)
    off = 0.17
    for ax in (ax1, ax2):
        ax.axvline(0, color="#888888", lw=1.0, ls=":")
        for y in ypos:
            ax.axhline(y - 0.5, color="white", lw=0.9)
    for y, r in zip(ypos, tab.itertuples()):
        ax1.errorbar(r.coef_at_p10, y + off, xerr=1.96 * r.se_at_p10, fmt="o", color=RED, markersize=6.5, **eb)
        ax1.errorbar(r.coef_at_p90, y - off, xerr=1.96 * r.se_at_p90, fmt="o", color=GREY, markerfacecolor="white",
                     markeredgecolor=GREY, markeredgewidth=1.3, markersize=6.5, **eb)
        ax2.errorbar(r.interaction, y, xerr=1.96 * r.se_interaction, fmt="s", color=INK, markersize=6, **eb)
        ax2.text(1.03, y, f"{r.p_cluster:.3f}", transform=ax2.get_yaxis_transform(), va="center", ha="left", fontsize=10)
    ax2.text(1.03, k - 0.15, "p", transform=ax2.get_yaxis_transform(), va="center", ha="left", fontsize=10, fontstyle="italic")
    ax1.text(-0.035, k - 0.05, "Target domain (t)", transform=ax1.get_yaxis_transform(), va="center", ha="right", fontsize=10, fontstyle="italic")
    ax1.set_yticks(ypos); ax1.set_yticklabels(tab.target, fontsize=10)
    for lab in ax1.get_yticklabels():
        if lab.get_text() in (WLB, JW):
            lab.set_fontweight("bold")
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax1.set_ylim(-1.4, k + 0.95)
    ax1.set_xlim(-0.287, 0.164); ax1.set_xticks([-0.2, -0.1, 0.0, 0.1])
    ax2.set_xlim(-0.178, 0.218); ax2.set_xticks([-0.1, 0.0, 0.1, 0.2])
    ax1.set_title("(a) Coefficient under fast and no decline", fontsize=10.4, fontweight="bold", loc="left", pad=3)
    ax2.set_title("(b) Change with the rate (slope in c)", fontsize=10.4, fontweight="bold", loc="left", pad=3)
    ax1.set_xlabel("Within-person coefficient of\nCaregiving Support (t−1)", fontsize=10)
    ax2.set_xlabel("Change in the coefficient per 1 %/yr of\npopulation change (interaction, 95% CI)", fontsize=10)
    ax1.text(-0.005, -0.93, "←  negative dependency", ha="right", va="center", fontsize=10, fontstyle="italic", color="#444444")
    ax2.text(-0.022, -0.93, "←  with growth", ha="right", va="center", fontsize=10, fontstyle="italic", color="#444444")
    ax2.text(0.006, -0.93, "more negative with decline →", ha="left", va="center", fontsize=9.4, fontstyle="italic", color="#444444")
    h1 = [ax1.errorbar([], [], xerr=[], fmt="o", color=RED, markersize=6.5, **eb),
          ax1.errorbar([], [], xerr=[], fmt="o", color=GREY, markerfacecolor="white", markeredgecolor=GREY, markeredgewidth=1.3, markersize=6.5, **eb)]
    ax1.legend(h1, [minus(f"fast decline (rate at P10, {r0.rate_p10:+.2f} %/yr)"), minus(f"no decline (rate at P90, {r0.rate_p90:+.2f} %/yr)")],
               loc="upper right", fontsize=9.6, framealpha=0.95, edgecolor="#bbbbbb", borderpad=0.25, labelspacing=0.2, handletextpad=0.3)
    ax2.legend([ax2.errorbar([], [], xerr=[], fmt="s", color=INK, markersize=6, **eb)], ["interaction with the rate (95% CI)"],
               loc="upper right", fontsize=9.6, framealpha=0.95, edgecolor="#bbbbbb", borderpad=0.25, handletextpad=0.3)
    # (c) coefficient along the rate
    lo, hi = max(float(quint.rate_lo.min()), -2.45), min(float(quint.rate_hi.max()), 0.92)
    if lines_json and Path(lines_json).exists():
        L = json.loads(Path(lines_json).read_text(encoding="utf-8"))
        if "SMC 65+" in L:
            lo, hi = L["SMC 65+"]["r_range"]
    grid = np.linspace(lo, hi, 200)
    ax3.axhline(0, color="#888888", lw=1.0, ls=":")
    ax3.axvline(THRESH, color="#444444", lw=1.3, ls=(0, (5, 3)))
    for rv, name in ((r0.rate_p10, "P10"), (r0.rate_p90, "P90")):
        ax3.axvline(rv, color="#888888", lw=1.0, ls=":")
        ax3.text(rv + 0.03, -0.405, name, ha="left", va="bottom", fontsize=10, color="#555555")
    ax3.text(THRESH - 0.02, 0.322, "shrinking (SCiRN)", ha="right", va="top", fontsize=10, color="#555555")
    ax3.text(THRESH + 0.04, 0.322, "non-shrinking", ha="left", va="top", fontsize=10, color="#555555")
    handles, labels = [], []
    for target, col, mk, face, dx in ((WLB, RED, "o", RED, 0.0), (JW, INK, "s", "white", 0.035)):
        row = tab[tab.target == target].iloc[0]
        line, se = fitted_line(row, grid)
        if target == WLB:
            ax3.fill_between(grid, line - 1.96 * se, line + 1.96 * se, color=RED, alpha=0.14, lw=0)
        ax3.plot(grid, line, color=col, lw=1.8, ls=(0, (5, 2)))
        q = quint[quint.target == target].sort_values("quintile")
        ax3.errorbar(q.rate_median + dx, q.slope_CS, yerr=1.96 * q.se_cluster, fmt=mk, color=col, markerfacecolor=face,
                     markeredgecolor=col, markeredgewidth=1.4, markersize=7, zorder=4, capsize=4, capthick=1.5, elinewidth=1.5)
        handles.append(ax3.errorbar([], [], xerr=[], fmt=mk, color=col, markerfacecolor=face, markeredgecolor=col, markeredgewidth=1.4,
                                    markersize=7, capsize=4, capthick=1.5, elinewidth=1.5))
        labels.append(f"{target} (slope " + minus(f"{row.interaction:+.3f}") + ")")
    ax3.legend(handles, labels, loc="upper left", fontsize=10, framealpha=0.95, edgecolor="#bbbbbb", borderpad=0.3, labelspacing=0.25)
    ax3.set_xlim(-2.55, 1.02); ax3.set_ylim(-0.42, 0.33)
    ax3.set_xticks(np.arange(-2.5, 1.01, 0.5)); ax3.set_yticks(np.arange(-0.4, 0.31, 0.1))
    ax3.xaxis.set_major_formatter(lambda v, _p: minus(f"{v:.1f}")); ax3.yaxis.set_major_formatter(lambda v, _p: minus(f"{v:.1f}"))
    for ax in (ax1, ax2):
        ax.xaxis.set_major_formatter(lambda v, _p: minus(f"{v:.1f}"))
    ax3.set_title("(c) Coefficient along the rate: Work-Life Balance and Jobs & Wages", fontsize=11.5, fontweight="bold", loc="left", pad=3)
    ax3.set_xlabel("Municipal annual population change rate, 2017–2023 (%/yr)", fontsize=10, labelpad=1)
    ax3.set_ylabel("Within-person coefficient of\nCaregiving Support (t−1)", fontsize=10)
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{fname}.{ext}", dpi=300)
    plt.close(fig)
    print(f"[figure] {fname} written (13 targets; n = {int(r0.n):,} lag-1 transitions, {int(r0.n_munis)} municipalities)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network-dir", required=True, help="output directory of the network run (strict3 / mci)")
    ap.add_argument("--strain-dir", default=None, help="output directory of the targeted run (lag1 / lag1only)")
    ap.add_argument("--care-dir", default=None, help="output directory of care_links.py")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--layout", choices=["manuscript", "refit"], default="manuscript")
    a = ap.parse_args(argv)
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    net = Path(a.network_dir) / "N_all_edges.csv"
    P, W = edge_dict(net, "pooled"), edge_dict(net, "within")
    with plt.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42, **set_font()}):
        ring = MANUSCRIPT_RING if a.layout == "manuscript" else refit_ring(
            [P.get(m, []) for _, m in FIG2] + [W.get(m, []) for _, m in FIG3])
        pos = ring_positions(ring)
        network_figure(P, FIG2, pos, "figure2", out)
        network_figure(W, FIG3, pos, "figure3", out)
        fam = Path(a.care_dir) / "V10_family_SMC_lag1.csv" if a.care_dir else None
        qs = Path(a.care_dir) / "V11_quintile_slopes.csv" if a.care_dir else None
        if fam and qs and fam.exists() and qs.exists():
            figure4(fam, qs, Path(a.strain_dir) / "V9_figure4_lines.json" if a.strain_dir else None, out)
        else:
            print("[figure] figure4 skipped: needs --care-dir with V10_family_SMC_lag1.csv and V11_quintile_slopes.csv")
    print(f"[done] {out.resolve()}")


if __name__ == "__main__":
    main()
