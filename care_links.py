#!/usr/bin/env python3
"""Care-work links among older adults aged 65+ (repository release v2.1.1).

Within-person regressions on person-demeaned scores, with the target's own lag, survey-year fixed effects and
municipality-clustered standard errors.  No model here contains a lag-2 term.

1. Screen of the 13 links with Caregiving Support as source.  For every target domain Y the model of Table 8,

       Y(t) = b0 + b1 CS(t-1) + b2 Y(t-1) + b3 CS(t-1) x rate + b4 Y(t-1) x rate + b5 rate + year fixed effects,

   is fitted on all valid lag-1 transitions (sample "lag1", as Table 8) and on three-consecutive-wave windows
   (sample "windows", the sample of the network models), for SMCs, large cities and both scales pooled.
   Reported per target: the CS(t-1) x rate interaction with its cluster p-value, a randomisation p-value (the
   rate is reassigned across municipalities within city scale), a family-wise p-value (max |t| over the 13
   targets under the same reassignments), and the coefficient of CS(t-1) implied at the 10th, 50th and 90th
   percentile of the municipal annual population change rate.
       V10_family_<scale>_<sample>.csv   Figure 4a-b; Supplementary Tables S4(c, d), S7(b), S8(a)
   The Work-Life Balance row on the lag-1 sample must reproduce Table 8; the script prints the check.

2. Coefficient of CS(t-1) in the model without the rate terms, within quintiles of the transitions ranked by
   the rate, for Work-Life Balance and Jobs & Wages.
       V11_quintile_slopes.csv           Supplementary Table S8(b); points of Figure 4c
   Same rule as validation_analysis.step_dose_extra (b), whose Work-Life Balance values it must reproduce
   (checked when --strain-dir is given).

3. The same model without the rate terms, separately in the shrinking and the non-shrinking subset.
       V11_subset_regressions.csv        Supplementary Table S7(b)

Usage
  python care_links.py --xlsx WORKBOOK --outdir CARE [--strain-dir TARGETED]      (or --data-csv and --pop-csv)
Run in the folder that contains validation_analysis.py (>= 2.1.1), or pass --va-dir; its loader builds the
panel, so the samples are those of every other script.  Reads the individual data, writes aggregates only.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CARE_LINKS_VERSION = "2.1.1"
MIN_WAVES = 3
CS, WLB, JW = "Caregiving Support", "Work-Life Balance", "Jobs & Wages"
SLAB = {"SMC": "SMCs", "Large": "Large cities", "pooled": "Both scales"}
# Published values, used only for the printed checks (Table 8; Supplementary Tables S7b and S8b).
TABLE8 = {("SMC", "lag1"): (0.103, 0.049, 1899), ("pooled", "lag1"): (0.115, 0.040, 3564),
          ("Large", "lag1"): (0.119, 0.072, 1630)}
PUBLISHED_SUBSET = {("SMC", JW, "shrinking"): (-0.123, 0.037, 1505), ("SMC", JW, "non-shrinking"): (-0.101, 0.062, 386),
                    ("Large", JW, "shrinking"): (-0.041, 0.039, 1127), ("Large", JW, "non-shrinking"): (0.120, 0.064, 497)}
PUBLISHED_QUINTILE = {("SMC", WLB): [-0.115, -0.159, -0.027, -0.087, 0.085], ("SMC", JW): [-0.096, -0.167, -0.046, -0.169, -0.089]}


def rd(path):
    return pd.read_csv(path, encoding="utf-8-sig")


# ----------------------------------------------------------------------------
# 1. screen of the 13 caregiving-source links
# ----------------------------------------------------------------------------
def stratum(d, scale):
    m = d.age65plus == 1
    if scale != "pooled":
        m &= d.city_scale == scale
    return d[m].sort_values(["respondent_id", "year"]).reset_index(drop=True)


def demean(sub, S):
    n = sub.groupby("respondent_id").year.transform("size")
    sub = sub[n >= MIN_WAVES].copy().reset_index(drop=True)
    sub[S] = sub[S] - sub.groupby("respondent_id")[S].transform("mean")
    return sub.sort_values(["respondent_id", "year"]).reset_index(drop=True)


def valid_index(d, history):
    """Rows t whose previous `history` rows are the same respondent in consecutive years."""
    ids, yrs = d.respondent_id.to_numpy(), d.year.to_numpy()
    idx = np.arange(history, len(d))
    good = np.ones(len(idx), bool)
    for lag in range(1, history + 1):
        good &= (ids[idx] == ids[idx - lag]) & (yrs[idx] - yrs[idx - lag] == lag)
    return idx[good]


def ols_cluster(y, X, starts, g):
    """OLS with CR1 cluster-robust covariance; rows must be sorted by cluster."""
    n, k = X.shape
    bread = np.linalg.inv(X.T @ X)
    b = bread @ (X.T @ y)
    u = y - X @ b
    scores = np.add.reduceat(X * u[:, None], starts, axis=0)
    V = bread @ (scores.T @ scores) @ bread.T * (g / (g - 1)) * ((n - 1) / (n - k))
    return b, V


def design(x, y1, r, year_fe, extra):
    cols = [np.ones(len(x)), x, y1, x * r, y1 * r, r, year_fe]
    if extra is not None:
        cols.append(extra)
    return np.column_stack(cols)


def run(analysis, S, CS, scale, sample, nperm, seed):
    sub = stratum(analysis, scale)
    d = demean(sub, S)
    idx = valid_index(d, 1 if sample == "lag1" else 2)
    muni = d.muni_code.to_numpy()[idx]
    order = np.argsort(muni, kind="stable")                  # sort rows by cluster once
    idx, muni = idx[order], muni[order]
    codes, inv = np.unique(muni, return_inverse=True)
    starts = np.flatnonzero(np.r_[True, inv[1:] != inv[:-1]])
    g = len(codes)
    vals = d[S].to_numpy(float)
    r = 100.0 * d["annual_change_rate"].to_numpy(float)[idx]
    yrs = d.year.to_numpy()[idx]
    lv = sorted(np.unique(yrs))
    year_fe = np.column_stack([(yrs == v).astype(float) for v in lv[1:]]) if len(lv) > 1 else np.empty((len(idx), 0))
    sc = d.city_scale.to_numpy()[idx]
    i = S.index(CS)
    x = vals[idx - 1, i]
    extra = None
    if scale == "pooled":
        large = (sc == "Large").astype(float)
        extra = np.column_stack([large, large * x])
    targets = [j for j in range(len(S)) if j != i]
    q10, q50, q90 = np.quantile(r, [.10, .50, .90])

    def tvec(rr):
        t = np.empty(len(targets))
        fits = []
        for a, j in enumerate(targets):
            X = design(x, vals[idx - 1, j], rr, year_fe, extra)
            b, V = ols_cluster(vals[idx, j], X, starts, g)
            t[a] = b[3] / np.sqrt(V[3, 3])
            fits.append((b, V))
        return t, fits

    t_obs, fits = tvec(r)

    # randomisation: reassign the rate across municipalities within city scale
    rng = np.random.default_rng([int(seed), 10, {"SMC": 1, "Large": 2, "pooled": 3}[scale],
                                 1 if sample == "lag1" else 2])
    muni_rate = np.zeros(g); muni_rate[inv] = r
    muni_scale = np.empty(g, dtype=object); muni_scale[inv] = sc
    blocks = [np.flatnonzero(muni_scale == s) for s in np.unique(muni_scale)]
    null = np.empty((nperm, len(targets)))
    for p in range(nperm):
        perm = muni_rate.copy()
        for m in blocks:
            perm[m] = muni_rate[m][rng.permutation(len(m))]
        null[p], _ = tvec(perm[inv])
        if (p + 1) % 500 == 0:
            print(f"  [{scale} {sample}] randomisation {p + 1}/{nperm}", flush=True)
    absnull = np.abs(null)
    maxnull = absnull.max(axis=1) if nperm else np.array([])

    rows = []
    for a, j in enumerate(targets):
        b, V = fits[a]
        row = {"scale": scale, "sample": sample, "target": S[j], "n": len(idx), "n_munis": g,
               "interaction": b[3], "se_interaction": np.sqrt(V[3, 3]), "t_cluster": t_obs[a],
               "p_cluster": float(2 * stats.t.sf(abs(t_obs[a]), g - 1))}
        if nperm:
            row["p_randomisation"] = float((1 + (absnull[:, a] >= abs(t_obs[a])).sum()) / (nperm + 1))
            row["p_family13"] = float((1 + (maxnull >= abs(t_obs[a])).sum()) / (nperm + 1))
        for name, rv in (("p10", q10), ("p50", q50), ("p90", q90)):
            gvec = np.zeros(len(b)); gvec[1] = 1.0; gvec[3] = rv
            row[f"rate_{name}"] = rv
            row[f"coef_at_{name}"] = float(gvec @ b)
            row[f"se_at_{name}"] = float(np.sqrt(gvec @ V @ gvec))
        rows.append(row)
    out = pd.DataFrame(rows)
    out["rank_abs_t"] = out.t_cluster.abs().rank(ascending=False).astype(int)
    return out.sort_values("interaction", ascending=False).reset_index(drop=True)


# ----------------------------------------------------------------------------
# 2-3. quintile coefficients and subset regressions
# ----------------------------------------------------------------------------
def target_design(va, sub, target):
    """Within-person design on valid lag-1 transitions (as validation_analysis._dose_design, any target)."""
    d = va.demean(sub)
    idx = va.valid_index(d, 1)
    yrs = d.year.to_numpy()[idx]
    lv = sorted(np.unique(yrs))
    year_fe = np.column_stack([(yrs == v).astype(float) for v in lv[1:]]) if len(lv) > 1 else np.empty((len(idx), 0))
    return {"y": d[target].to_numpy(float)[idx], "x": d[CS].to_numpy(float)[idx - 1], "y1": d[target].to_numpy(float)[idx - 1],
            "r": 100.0 * d["annual_change_rate"].to_numpy(float)[idx], "year_fe": year_fe,
            "muni": d.muni_code.to_numpy()[idx], "n": len(idx)}


def slope(va, D, m=None):
    """Coefficient of Caregiving Support (t-1): y on x, own lag and year FE; CR1 cluster SE by municipality."""
    m = np.ones(D["n"], bool) if m is None else m
    fe = D["year_fe"][m]
    fe = fe[:, fe.any(axis=0)]
    X = np.column_stack([np.ones(int(m.sum())), D["x"][m], D["y1"][m], fe])
    b, se, _ = va._ols_cluster(D["y"][m], X, D["muni"][m])
    g = len(np.unique(D["muni"][m]))
    return float(b[1]), float(se[1]), int(m.sum()), g


def estimate(va, analysis, out, scales):
    qrows, srows = [], []
    for scale in scales:
        for target in (WLB, JW):
            D = target_design(va, va.stratum(analysis, scale, True, None), target)
            rr = D["r"]
            e = np.quantile(rr, [0, .2, .4, .6, .8, 1.0])
            for q in range(5):
                m = (rr >= e[q]) & (rr <= e[q + 1]) if q == 4 else (rr >= e[q]) & (rr < e[q + 1])
                b, se, n, g = slope(va, D, m)
                qrows.append({"scale": scale, "target": target, "quintile": q + 1, "rate_lo": float(e[q]), "rate_hi": float(e[q + 1]),
                              "rate_median": float(np.median(rr[m])), "n": n, "n_munis": g, "slope_CS": b, "se_cluster": se})
            for shrink, lab in ((True, "shrinking"), (False, "non-shrinking")):
                D = target_design(va, va.stratum(analysis, scale, True, shrink), target)
                b, se, n, g = slope(va, D)
                srows.append({"scale": scale, "target": target, "subset": lab, "n": n, "n_munis": g, "coef_CS": b, "se_cluster": se,
                              "p_cluster_t": float(2 * stats.t.sf(abs(b / se), g - 1))})
    Q, R = pd.DataFrame(qrows), pd.DataFrame(srows)
    Q.to_csv(out / "V11_quintile_slopes.csv", index=False, encoding="utf-8-sig")
    R.to_csv(out / "V11_subset_regressions.csv", index=False, encoding="utf-8-sig")
    return Q, R


def check_published(Q, R, strain_dir):
    """Printed comparison with the published values and with the targeted run (never stops the run)."""
    bad = 0
    print("  [check] comparison with the published values (meaningful for the survey workbook only)", flush=True)

    def line(name, got, exp, tol=0.0015):
        nonlocal bad
        ok = all(abs(g - x) <= tol for g, x in zip(np.atleast_1d(got), np.atleast_1d(exp)))
        bad += not ok
        print(f"  [check] {name}: {'OK' if ok else 'MISMATCH'}", flush=True)

    for (scale, target, subset), (b, se, n) in PUBLISHED_SUBSET.items():
        r = R[(R.scale == scale) & (R.target == target) & (R.subset == subset)]
        if len(r):
            line(f"S7b {SLAB[scale]}, {subset}: {b:+.3f} ({se:.3f}), n = {n:,}", [r.coef_CS.iloc[0], r.se_cluster.iloc[0], r.n.iloc[0]], [b, se, n])
    for (scale, target), exp in PUBLISHED_QUINTILE.items():
        q = Q[(Q.scale == scale) & (Q.target == target)].sort_values("quintile")
        if len(q):
            line(f"S8b {SLAB[scale]}, {target}, five quintiles", q.slope_CS.to_numpy(), exp)
    f = Path(strain_dir) / "V9_figure4_quintile_slopes.csv" if strain_dir else None
    if f and f.exists():
        v9 = rd(f).sort_values(["scale", "quintile"])
        mine = Q[Q.target == WLB].merge(v9, on=["scale", "quintile"], suffixes=("", "_v9"))
        if len(mine):
            line("Work-Life Balance quintiles equal the targeted run (V9_figure4_quintile_slopes.csv)",
                 np.r_[mine.slope_CS, mine.se_cluster, mine.n], np.r_[mine.slope_CS_v9, mine.se_cluster_v9, mine.n_v9], tol=1e-9)
    print(f"  [check] {'all checks passed' if not bad else str(bad) + ' check(s) failed: do not use these outputs before resolving'}", flush=True)
    return bad == 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx"); ap.add_argument("--data-csv"); ap.add_argument("--pop-csv")
    ap.add_argument("--va-dir", default=str(Path(__file__).resolve().parent), help="folder containing validation_analysis.py")
    ap.add_argument("--strain-dir", default=None, help="output directory of the targeted run (for a cross-check)")
    ap.add_argument("--outdir", default="care")
    ap.add_argument("--scales", default="SMC,pooled,Large", help="comma-separated: SMC, pooled, Large")
    ap.add_argument("--samples", default="lag1,windows", help="comma-separated: lag1, windows")
    ap.add_argument("--nperm", type=int, default=2000, help="randomisation replicates of the 13-target screen")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args(argv)
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(Path(a.va_dir).resolve()))
    import validation_analysis as va
    _, analysis, _ = va.load_frame(xlsx=a.xlsx, data_csv=a.data_csv, pop_csv=a.pop_csv)
    S = list(va.S)
    scales = a.scales.split(",")
    for scale in scales:
        for sample in a.samples.split(","):
            tab = run(analysis, S, CS, scale, sample, a.nperm, a.seed)
            tab.to_csv(out / f"V10_family_{scale}_{sample}.csv", index=False, encoding="utf-8-sig")
            print(f"[care-links] {SLAB[scale]}, 65+, {sample}: n = {tab.n.iloc[0]:,} ({tab.n_munis.iloc[0]} municipalities); "
                  f"rate P10/P50/P90 = {tab.rate_p10.iloc[0]:+.2f}/{tab.rate_p50.iloc[0]:+.2f}/{tab.rate_p90.iloc[0]:+.2f} %/yr", flush=True)
            ref = TABLE8.get((scale, sample))
            if ref:
                w = tab[tab.target == WLB].iloc[0]
                ok = abs(w.interaction - ref[0]) < 0.0015 and abs(w.se_interaction - ref[1]) < 0.0015 and int(w.n) == ref[2]
                print(f"  [check] Table 8, {SLAB[scale]} (Work-Life Balance: {ref[0]:.3f}, SE {ref[1]:.3f}, n = {ref[2]:,}): "
                      f"{'OK' if ok else 'MISMATCH'}", flush=True)
    Q, R = estimate(va, analysis, out, [s for s in scales if s != "pooled"])
    check_published(Q, R, a.strain_dir)
    print(f"[care-links] written to {out.resolve()}", flush=True)


if __name__ == "__main__":
    main()
