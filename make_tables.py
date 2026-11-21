#!/usr/bin/env python3
"""All tables of the article and its supplement (repository release v2.1.1).

Assembles Tables 1-8 and Supplementary Tables S1-S8 from the outputs of the other scripts:

  network run   validation_analysis.py --window strict3 --design mci    -> Tables 3-6, Block A of Table 7
  targeted run  validation_analysis.py --window lag1 --design lag1only  -> Block B of Table 7, Table 8, S2-S4
  care links    care_links.py                                           -> S4(c, d), S7(b), S8

Tables 1-2 (sample structure, descriptive statistics) and the conditioning depths (Table 5, columns (1)-(3);
Supplementary Table S5) are computed from the individual data and are written only when the workbook is given.

Usage
  python make_tables.py --network-dir NET --strain-dir STR --care-dir CARE --outdir tables --xlsx WORKBOOK
  (or --data-csv and --pop-csv; without data inputs the tables that need them are skipped)
Writes one CSV per table and tables.xlsx (one sheet per table, plus a README sheet).  With data inputs it needs
validation_analysis.py (>= 2.1.1) in the same folder; it reads the individual data and writes aggregates only.
"""
import argparse, json, math, sys
from pathlib import Path
import numpy as np, pandas as pd

MAKE_TABLES_VERSION = "2.1.3"
Q, FLOOR, LOOSE = 0.05, 0.10, 0.05
CS, WLB, JW, HUB = "Caregiving Support", "Work-Life Balance", "Jobs & Wages", "Life Satisfaction"
DOMAINS = ["Life Satisfaction", "Finances & Assets", "Jobs & Wages", "Housing", "Work-Life Balance", "Health", "Education",
           "Friends & Community", "Trust in Government & Courts", "Air & Water (environment)", "Personal Safety",
           "Parenting Support", "Caregiving Support", "Enjoyment of Life"]
MAIN8 = [f"{a} {b}, aged {c}" for b in ("Large cities", "SMCs") for a in ("All", "Shrinking") for c in ("15-64", "65+")]
STRATA8 = [(sc, older, sh) for sc in ("Large", "SMC") for sh in (None, True) for older in (False, True)]
SCALES = ("SMCs", "Large cities")
SLAB = {"SMC": "SMCs", "Large": "Large cities", "pooled": "Both scales"}
SPEC_NET = "Full MCI, three-consecutive-wave windows"
SPEC_STR = "Lagged partial correlation (target parents only), all lag-1 transitions"
GROUPS = [("All Large Cities", "Large", None), ("All SMCs", "SMC", None), ("Shrinking Large Cities", "Large", True), ("Shrinking SMCs", "SMC", True)]
MALE_CODE = 1                                   # sex code of men in the survey file
# Published Table 1 rows that Gate 1 of validation_analysis.py does not cover (printed check only).
PUBLISHED_TABLE1 = {"Respondents with >=1 consecutive pair": [4068, 4512, 2578, 3505],
                    "Valid three-consecutive-wave windows (aged 15-64)": [2904, 3202, 1793, 2507],
                    "Valid three-consecutive-wave windows (aged 65+)": [882, 1048, 651, 834],
                    "Sex: male (%)": [49.6, 50.7, 49.8, 50.4], "Age: 65+ (%)": [16.0, 16.6, 15.9, 16.6],
                    "COVID-period waves (%)": [19.0, 18.0, 17.7, 17.7]}
CONTENT = {
    "T1_sample_structure": "Table 1: sample structure",
    "T2_descriptives": "Table 2: mean and total, between-person and within-person SD of the 14 measures (long format)",
    "T2_descriptives_wide": "Table 2 in the layout of the article",
    "T3_edge_counts": "Table 3: number of reported edges per stratum (network run), main and sensitivity criteria",
    "T4_pooled_edges": "Table 4: reported pooled edges (network run)",
    "T5_out_strength": "Table 5 / S6: out-strength centrality of each source domain, pooled models (network run)",
    "T6_within_edges": "Table 6: reported person-demeaned edges with time-permutation null bands (network run)",
    "T6b_null_summary_and_care_links": "Table 6 note: edges under the time-permutation null; caregiving-source links against their null bands",
    "T7_decline_contrast": "Table 7: decline contrast of Caregiving Support -> Work-Life Balance (65+) under the two specifications",
    "T7b_anderson_hsiao": "Table S4(a): Anderson-Hsiao IV in first differences (three-wave histories)",
    "T7c_conditioning_conventions": "Table 7 note / S2: contrast under alternative conditioning conventions, both specifications",
    "T7d_other_contrasts": "Table 7: age-pattern check (15-64); Table S7(a): Caregiving Support -> Jobs & Wages, both specifications",
    "T7e_position_and_period": "Table S3: contrast by position in the run of waves and by year (targeted run)",
    "T7f_window_by_design": "Table S3: contrast by window and design (targeted run)",
    "T8_dose_response": "Table 8 / S4(b): continuous moderation by the municipal population change rate (SMCs first)",
    "T8b_family_test": "Table S4(c): family-wise test over the 182 lagged links (classical t)",
    "T8c_region_robustness": "Table 8 note / S4(d): region fixed effects and region-blocked randomisation",
    "S_two_layers_summary": "Not in the article: the network run and the lag-1-only network layer side by side",
    "T5b_out_strength_both_layers": "Not in the article: out-strength under the full MCI and under the lag-1-only design",
    "S_lag1only_edge_counts": "Not in the article: edge counts under the lag-1-only design",
    "S_lag1only_out_strength": "Not in the article: out-strength under the lag-1-only design",
    "S_depth_negative_pooled_edges": "Table S5: conditioning depth of the negative pooled edges of Table 4",
    "S_depth_life_satisfaction_hub": "Table 5, columns (1)-(5): mean coefficient of the 13 Life Satisfaction -> domain links by conditioning depth",
    "S4c_interaction_placement": "Table S4(c): placement of the interaction among the 182 links, incl. cluster-robust family-wise p",
    "S4d_rate_region_windows": "Table S4(d): rate percentiles, region fixed effects, re-estimation on windows",
    "S7b_jobs_wages_regressions": "Table S7(b): Caregiving Support -> Jobs & Wages, subset regressions and rate interaction",
    "S8a_family_13_targets": "Table S8(a): the interaction model for each of the 13 targets of Caregiving Support (SMCs, 65+)",
    "S8b_quintile_coefficients": "Table S8(b): coefficient within quintiles of the rate, Work-Life Balance and Jobs & Wages",
}


def norm_p(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def fisher(r1, n1, r2, n2):
    z = (math.atanh(r1) - math.atanh(r2)) / math.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    return z, norm_p(z)


def rd(path):
    return pd.read_csv(path, encoding="utf-8-sig")


def jl(path):
    return json.loads(Path(path).read_text(encoding="utf-8")) if Path(path).exists() else None


def edge_counts(e):
    rows = []
    for m in MAIN8:
        r = {"stratum": m}
        for mode in ("pooled", "within"):
            d = e[(e.model == m) & (e["mode"] == mode)]
            main = (d.q < Q) & (d.MCI.abs() >= FLOOR)
            r.update({f"{mode}: n (transitions)": int(d.n_valid.iloc[0]), f"{mode}: edges (q<.05, |r|>=.10)": int(main.sum()),
                      f"{mode}: of which negative": int((main & (d.MCI < 0)).sum()),
                      f"{mode}: q<.05, |r|>=.05": int(((d.q < Q) & (d.MCI.abs() >= LOOSE)).sum()),
                      f"{mode}: q<.05 only": int((d.q < Q).sum()), f"{mode}: p<.05, |r|>=.10": int(((d.p < .05) & (d.MCI.abs() >= FLOOR)).sum())})
        rows.append(r)
    return pd.DataFrame(rows)


def out_strength(e):
    d = e[(e["mode"] == "pooled") & e.reported & e.model.isin(MAIN8)]
    t = d.assign(a=d.MCI.abs()).pivot_table(index="source", columns="model", values="a", aggfunc="sum").reindex(index=DOMAINS, columns=MAIN8).fillna(0.0)
    t.loc["(number of reported edges)"] = [int(((d.model == m)).sum()) for m in MAIN8]
    return t.reset_index().rename(columns={"source": "source domain"})


def null_bands(d):
    out = {}
    for sc in SCALES:
        f = d / f"V3_table6_edges_vs_null_{'SMC' if sc == 'SMCs' else 'Large'}.csv"
        if f.exists():
            out[sc] = rd(f)
    return out


def contrast_rows(d, spec, nb):
    """Decline contrast of Caregiving Support (t-1) -> Work-Life Balance (t), 65+, for one specification."""
    v2, v4 = rd(d / "V2_common_set.csv"), rd(d / "V4_cluster_boot.csv")
    rows = []
    for sc in SCALES:
        key = "SMC" if sc == "SMCs" else "Large"
        a = v2[(v2.scale == sc) & (v2.conditioning == "own (Table 7)")].iloc[0]; b = v4[(v4.scale == sc) & (v4.convention == "own")].iloc[0]
        r = {"specification": spec, "city scale": sc, "shrinking: coefficient": a.shrinking_r, "shrinking: n": int(a.n_s),
             "non-shrinking: coefficient": a.non_shrinking_r, "non-shrinking: n": int(a.n_n), "Fisher z": a.z, "p (Fisher z)": a.p_iid,
             "cluster-bootstrap z": b.z_cluster_boot, "p (cluster bootstrap)": b.p_cluster_boot_wald,
             "municipalities (shrinking / non-shrinking)": f"{int(b.n_munis_s)} / {int(b.n_munis_n)}"}
        if sc in nb:
            t = nb[sc]; t = t[(t.model == "shrinking 65+") & (t.source == CS) & (t.target == WLB)]
            if len(t):
                r.update({"time-permutation null 95% band (shrinking)": f"[{t.null_q025.iloc[0]:.3f}, {t.null_q975.iloc[0]:.3f}]",
                          "time-permutation p (shrinking coefficient)": t.perm_p_two_tail.iloc[0]})
        for mode, lab in (("global", ""), ("region", ", region-blocked")):
            j = jl(d / f"V6_perm_label_{key}_{mode}.json")
            if j:
                r.update({f"label-permutation p: this link{lab}": j["paper__p_single_edge"],
                          f"label-permutation p: 13 caregiving-source links{lab}": j["paper__p_family_CS13"],
                          f"label-permutation p: all 182 links{lab}": j["paper__p_family_all182"]})
        rows.append(r)
    return rows


def other_contrasts(e, spec):
    """Fisher z contrasts from the network edge lists (conditioning sets as selected in each subset model)."""
    rows = []
    for tgt in (WLB, JW):
        for sc in SCALES:
            for age in ("65+", "15-64"):
                g = lambda lab: e[(e["mode"] == "within") & (e.model == f"{lab} {sc}, aged {age}") & (e.source == CS) & (e.target == tgt)].iloc[0]
                s, n = g("Shrinking"), g("Non-shrinking"); z, p = fisher(s.MCI, s.n_valid, n.MCI, n.n_valid)
                rows.append({"specification": spec, "dependency": f"{CS} -> {tgt}", "city scale": sc, "cohort": age,
                             "shrinking: coefficient": s.MCI, "shrinking: q": s.q, "shrinking: n": int(s.n_valid), "shrinking: reported edge": bool(s.reported),
                             "non-shrinking: coefficient": n.MCI, "non-shrinking: n": int(n.n_valid), "Fisher z": z, "p (Fisher z)": p})
    return rows


def article_tables(N, S):
    """Tables 3-8 and the supplementary tables that need only the aggregate outputs of the two runs."""
    en, es = rd(N / "N_all_edges.csv"), rd(S / "N_all_edges.csv")
    for e in (en, es):
        e["reported"] = e["reported"].astype(bool)
    nbn, nbs = null_bands(N), null_bands(S)
    T = {}

    T["T3_edge_counts"] = edge_counts(en)
    d = en[(en["mode"] == "pooled") & en.reported & en.model.isin(MAIN8)].copy()
    d["stratum"] = pd.Categorical(d.model, MAIN8, ordered=True)
    T["T4_pooled_edges"] = d.sort_values(["stratum", "MCI"], ascending=[True, False])[["stratum", "source", "target", "MCI", "q", "n_valid"]].rename(columns={"n_valid": "n"})
    T["T5_out_strength"] = out_strength(en)

    rows = []                                                       # Table 6: reported within-person edges + null bands
    w = en[(en["mode"] == "within") & en.reported & en.model.isin(MAIN8)]
    for r in w.itertuples():
        row = {"stratum": r.model, "source": r.source, "target": r.target, "coefficient": r.MCI, "q": r.q, "n": int(r.n_valid)}
        sc = "SMCs" if "SMCs" in r.model else "Large cities"; mdl = ("shrinking 65+" if r.model.startswith("Shrinking") else "all 65+") if r.model.endswith("65+") else None
        if mdl and sc in nbn:
            t = nbn[sc]; t = t[(t.model == mdl) & (t.source == r.source) & (t.target == r.target)]
            if len(t):
                row.update({"null mean": t.null_mean.iloc[0], "null 95% band": f"[{t.null_q025.iloc[0]:.3f}, {t.null_q975.iloc[0]:.3f}]", "time-permutation p": t.perm_p_two_tail.iloc[0]})
        rows.append(row)
    T["T6_within_edges"] = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["stratum", "source", "target", "coefficient", "q", "n"])
    rows = []
    for sc in SCALES:
        j = jl(N / f"V3_perm_time_{'SMC' if sc == 'SMCs' else 'Large'}.json")
        if j:
            rows.append({"city scale": sc, "replicates": j["nrep"], "share of null replicates with >= 1 reported edge (shrinking 65+)": j["null_share_reps_with_reported_edge_s"],
                         "mean number of reported edges per null replicate": j["null_mean_reported_edges_s"]})
        if sc in nbn:
            t = nbn[sc]; t = t[(t.source == CS) & t.target.isin([WLB, JW])]
            for r in t.itertuples():
                rows.append({"city scale": sc, "model": r.model, "dependency": f"{r.source} -> {r.target}", "coefficient": r.mci_obs, "reported edge": bool(r.reported),
                             "null mean": r.null_mean, "null 95% band": f"[{r.null_q025:.3f}, {r.null_q975:.3f}]", "time-permutation p": r.perm_p_two_tail})
    T["T6b_null_summary_and_care_links"] = pd.DataFrame(rows)

    T["T7_decline_contrast"] = pd.DataFrame(contrast_rows(N, SPEC_NET, nbn) + contrast_rows(S, SPEC_STR, nbs))      # Table 7, two specifications
    v5 = rd(S / "V5_anderson_hsiao.csv"); v5 = v5[v5.specification == "year_fe"]
    T["T7b_anderson_hsiao"] = v5[["scale", "subset", "n_equations", "n_respondents", "n_munis", "gamma_CS_to_WLB", "se_cluster_muni", "p_cluster_t", "first_stage_F_dWLB", "first_stage_F_dCS"]]
    v2 = pd.concat([rd(N / "V2_common_set.csv").assign(specification=SPEC_NET), rd(S / "V2_common_set.csv").assign(specification=SPEC_STR)])
    T["T7c_conditioning_conventions"] = v2[["specification", "scale", "conditioning", "shrinking_r", "n_s", "non_shrinking_r", "n_n", "z", "p_iid", "Z_s", "Z_n"]]
    T["T7d_other_contrasts"] = pd.DataFrame(other_contrasts(en, SPEC_NET) + other_contrasts(es, SPEC_STR))
    if (S / "V8_period_position_contrasts.csv").exists():
        T["T7e_position_and_period"] = rd(S / "V8_period_position_contrasts.csv").sort_values(["scale", "rows"])
    if (S / "V7_window_sensitivity.csv").exists():
        T["T7f_window_by_design"] = rd(S / "V7_window_sensitivity.csv")[["scale", "window", "design", "n_s", "n_n", "mci_s", "mci_n", "z_paper", "p_iid", "Z_s", "Z_n"]]

    if (S / "V9_dose_response.csv").exists():                       # Table 8
        d9 = rd(S / "V9_dose_response.csv")
        order = {("SMC", "65+"): 0, ("SMC", "15-64"): 1, ("Large", "65+"): 2, ("Large", "15-64"): 3, ("pooled", "65+"): 4, ("pooled", "15-64"): 5}
        d9["o"] = [order[(s, c)] for s, c in zip(d9.scale, d9.cohort)]
        T["T8_dose_response"] = d9.sort_values(["o", "estimator"], ascending=[True, False]).drop(columns="o")
        for f, k in (("V9_dose_response_family.csv", "T8b_family_test"), ("V9_dose_response_region.csv", "T8c_region_robustness")):
            if (S / f).exists():
                T[k] = rd(S / f)

    sn, ss = rd(N / "N_network_summary.csv"), rd(S / "N_network_summary.csv")                       # supplement
    rows = []
    for m in MAIN8:
        r = {"stratum": m}
        for lab, sm in (("full MCI, 3-wave", sn), ("lag-1-only, all lag-1", ss)):
            for mode in ("pooled", "within"):
                x = sm[(sm.model == m) & (sm["mode"] == mode)].iloc[0]
                r.update({f"{lab} | {mode}: n": int(x.n_valid), f"{lab} | {mode}: edges": int(x.n_edges), f"{lab} | {mode}: negative": int(x.n_negative)})
                if mode == "pooled":
                    r[f"{lab} | pooled: Life Satisfaction out-strength"] = x.LS_out_strength
        rows.append(r)
    T["S_two_layers_summary"] = pd.DataFrame(rows)
    a, b = out_strength(en).set_index("source domain"), out_strength(es).set_index("source domain")       # Table 5, both layers side by side
    both = pd.concat({"full MCI, 3-wave": a, "lag-1-only, all lag-1": b}, axis=1).swaplevel(axis=1)
    both = both[[(m, l) for m in MAIN8 for l in ("full MCI, 3-wave", "lag-1-only, all lag-1")]]
    both.columns = [f"{m} | {l}" for m, l in both.columns]
    T["T5b_out_strength_both_layers"] = both.reset_index()
    T["S_lag1only_edge_counts"] = edge_counts(es)
    T["S_lag1only_out_strength"] = out_strength(es)
    return T, en, es



# Tables 1-2 (individual data)
def table1(va, frame):
    cols = {}
    for name, scale, shrink in GROUPS:
        sub = va.stratum(frame, scale, None, shrink)
        idx = va.valid_index(sub, 1)
        cols[name] = {
            "Number of municipalities": sub.muni_code.nunique(),
            "Number of shrinking municipalities (SCiRN)": sub.loc[sub.shrinking == 1, "muni_code"].nunique(),
            "Median annual population change (%/yr)": round(100 * sub.groupby("muni_code")["annual_change_rate"].first().median(), 2),
            "Person-years": len(sub), "Unique respondents": sub.respondent_id.nunique(),
            "Respondents with >=1 consecutive pair": int(sub.respondent_id.to_numpy()[idx].size and np.unique(sub.respondent_id.to_numpy()[idx]).size),
            "Valid lag-1 transitions (total)": len(idx),
            "Valid lag-1 transitions (aged 15-64)": len(va.valid_index(va.stratum(frame, scale, False, shrink), 1)),
            "Valid lag-1 transitions (aged 65+)": len(va.valid_index(va.stratum(frame, scale, True, shrink), 1)),
            "Valid three-consecutive-wave windows (aged 15-64)": len(va.valid_index(va.stratum(frame, scale, False, shrink), 2)),
            "Valid three-consecutive-wave windows (aged 65+)": len(va.valid_index(va.stratum(frame, scale, True, shrink), 2)),
            "Sex: male (%)": round(100 * float((sub.sex.round() == MALE_CODE).mean()), 1),
            "Age: 15-64 (%)": round(100 * float((sub.age65plus == 0).mean()), 1),
            "Age: 65+ (%)": round(100 * float((sub.age65plus == 1).mean()), 1),
            "COVID-period waves (%)": round(100 * float(sub.year.isin(va.COVID_YEARS).mean()), 1)}
    t = pd.DataFrame(cols)
    t.index.name = "item"
    return t.reset_index()


def table2(va, frame):
    rows = []
    for name, scale, shrink in GROUPS:
        for older in (False, True):
            sub = va.stratum(frame, scale, older, shrink)
            for c in va.S:
                pm = sub.groupby("respondent_id")[c].mean()
                within = sub[c] - sub.groupby("respondent_id")[c].transform("mean")
                rows.append({"city_group": name, "age_group": "65+" if older else "15-64", "domain": c, "mean": sub[c].mean(),
                             "sd_total": sub[c].std(ddof=1), "sd_between": pm.std(ddof=1), "sd_within": within.std(ddof=1)})
    long = pd.DataFrame(rows)
    long["cell"] = [f"{r.mean:.2f} ({r.sd_total:.2f}; {r.sd_between:.2f}; {r.sd_within:.2f})" for r in long.itertuples()]
    long["block"] = np.where(long.city_group.str.startswith("All"), "All Cities", "Shrinking Cities")
    long["column"] = np.where(long.city_group.str.contains("Large"), "Large cities, ", "SMCs, ") + long.age_group
    wide = long.pivot_table(index=["block", "domain"], columns="column", values="cell", aggfunc="first")
    wide = wide[["Large cities, 15-64", "Large cities, 65+", "SMCs, 15-64", "SMCs, 65+"]].reindex(
        [(b, d) for b in ("All Cities", "Shrinking Cities") for d in va.S]).reset_index()
    return long.drop(columns=["cell", "block", "column"]), wide


def check_table1(t1):
    bad = 0
    print("  [check] comparison with the published values (meaningful for the survey workbook only)", flush=True)
    for item, exp in PUBLISHED_TABLE1.items():
        got = t1.loc[t1.item == item, [g[0] for g in GROUPS]].to_numpy(float).ravel()
        ok = bool(len(got)) and all(abs(g - e) <= 0.051 for g, e in zip(got, exp))
        bad += not ok
        print(f"  [check] Table 1, {item}: {'OK' if ok else 'MISMATCH'}", flush=True)
    print(f"  [check] {'all checks passed' if not bad else str(bad) + ' check(s) failed (Table 1 rows not covered by Gate 1)'}", flush=True)


def depth_tables(va, frame, en, es):
    """Conditioning depths on three-consecutive-wave windows of the pooled strata: (1) zero-order lagged
    correlation, (2) given Y(t-1), (3) given Y(t-1) and X(t-2), (4) full MCI of the network run, (5) the
    person-demeaned full MCI.  Part A: every reported negative pooled edge; Part B: mean over the 13
    Life Satisfaction -> domain links in each of the eight strata."""
    va.DESIGN, va.HISTORY = "mci", 2                       # three-consecutive-wave rows, lag-2 terms allowed
    get = lambda e, m, mode, s, t: float(e[(e.model == m) & (e["mode"] == mode) & (e.source == s) & (e.target == t)].MCI.iloc[0])

    def depth(sub, s, t):
        j, i = va.S.index(t), va.S.index(s)
        r0, n = va.parcorr_masked(sub, s, t, [])
        r1, _ = va.parcorr_masked(sub, s, t, [(j, -1)])
        r2, _ = va.parcorr_masked(sub, s, t, [(j, -1), (i, -2)])
        return n, r0, r1, r2

    rows_a, rows_b = [], []
    for sc, older, sh in STRATA8:
        m = va.model_label(sc, older, sh); sub = va.stratum(frame, sc, older, sh)
        neg = en[(en.model == m) & (en["mode"] == "pooled") & en.reported & (en.MCI < 0)]
        for r in neg.itertuples():
            n, r0, r1, r2 = depth(sub, r.source, r.target)
            rows_a.append({"stratum": m, "source": r.source, "target": r.target, "n (three-wave)": n, "(1) zero-order lagged r": r0,
                           "(2) given Y(t-1)": r1, "(3) given Y(t-1), X(t-2)": r2, "(4) full MCI, pooled": r.MCI, "q": r.q,
                           "number of conditions in (4)": len(str(r.Z).split(";")), "Life Satisfaction in the conditioning set": "Life Satisfaction" in str(r.Z),
                           "(5) full MCI, person-demeaned": get(en, m, "within", r.source, r.target),
                           "lag-1-only, pooled (all lag-1 transitions)": get(es, m, "pooled", r.source, r.target)})
        d = np.array([depth(sub, va.HUB, t)[1:] for t in va.S if t != va.HUB])
        f = np.array([get(en, m, "pooled", va.HUB, t) for t in va.S if t != va.HUB]); w = np.array([get(en, m, "within", va.HUB, t) for t in va.S if t != va.HUB])
        l = np.array([get(es, m, "pooled", va.HUB, t) for t in va.S if t != va.HUB])
        rows_b.append({"stratum": m, "n (three-wave)": depth(sub, va.HUB, va.S[1])[0], "(1) zero-order lagged r, mean": d[:, 0].mean(), "(2) given Y(t-1), mean": d[:, 1].mean(),
                       "(3) given Y(t-1), X(t-2), mean": d[:, 2].mean(), "(4) full MCI pooled, mean": f.mean(), "(4) links with |r| >= .10": int((np.abs(f) >= .10).sum()),
                       "(5) full MCI person-demeaned, mean": w.mean(), "lag-1-only pooled, mean": l.mean(), "lag-1-only: links with |r| >= .10": int((np.abs(l) >= .10).sum())})
    A, B = pd.DataFrame(rows_a), pd.DataFrame(rows_b)
    return A, B



# Supplementary tables on the care-work links (aggregate outputs of care_links.py)
def fam(family_dir, scale, sample):
    f = Path(family_dir) / f"V10_family_{scale}_{sample}.csv" if family_dir else None
    return rd(f) if f and f.exists() else None


def care_tables(strain_dir, family_dir):
    S = Path(strain_dir) if strain_dir else None
    T = {}
    has = lambda name: S is not None and (S / name).exists()

    if has("V9_dose_response_family.csv"):                                   # S4(c)
        d = rd(S / "V9_dose_response_family.csv")
        rows = []
        for scale in ("pooled", "SMC", "Large"):
            r = d[d.scale == scale]
            if not len(r):
                continue
            r = r.iloc[0]
            row = {"sample": f"{SLAB[scale]}, 65+", "OLS t of interaction": r.t_CS_x_rate, "rank among 182": int(r.rank_among_182),
                   "family-wise p, 13 links": r.p_family_CS13, "family-wise p, 182 links": r.p_family_all182,
                   "null max |t| 13 links, median and 97.5th pct": r.null_max_abs_t_CS13_q50_q975,
                   "null max |t| 182 links, median and 97.5th pct": r.null_max_abs_t_all182_q50_q975}
            v = fam(family_dir, scale, "lag1")
            if v is not None and "p_family13" in v.columns:
                row["family-wise p, 13 links, cluster-robust t"] = float(v[v.target == WLB].p_family13.iloc[0])
            rows.append(row)
        T["S4c_interaction_placement"] = pd.DataFrame(rows)

    if has("V9_dose_response.csv"):                                          # S4(d)
        d = rd(S / "V9_dose_response.csv"); d = d[d.estimator.str.startswith("within OLS")]
        reg = rd(S / "V9_dose_response_region.csv") if has("V9_dose_response_region.csv") else None
        rows = []
        for scale, cohort in (("SMC", "65+"), ("SMC", "15-64"), ("Large", "65+"), ("pooled", "65+"), ("pooled", "15-64")):
            r = d[(d.scale == scale) & (d.cohort == cohort)]
            if not len(r):
                continue
            r = r.iloc[0]
            row = {"sample": f"{SLAB[scale]}, {cohort}", "rate P10": r.r_p10, "rate P50": r.r_p50, "rate P90": r.r_p90}
            if reg is not None and cohort == "65+" and len(reg[reg.scale == scale]):
                g = reg[reg.scale == scale].iloc[0]
                row.update({"region FE: interaction": g.beta_interaction_region_FE, "region FE: SE": g.se_region_FE,
                            "region FE: p": g.p_cluster_t_region_FE, "region-blocked randomisation p": g.p_randomisation_region_blocks})
            v = fam(family_dir, scale, "windows") if cohort == "65+" else None
            if v is not None:
                w = v[v.target == WLB].iloc[0]
                row.update({"windows: interaction": w.interaction, "windows: SE": w.se_interaction, "windows: p": w.p_cluster, "windows: n": int(w.n)})
            rows.append(row)
        T["S4d_rate_region_windows"] = pd.DataFrame(rows)

    E = Path(family_dir) if family_dir else Path(".")
    if (E / "V11_subset_regressions.csv").exists():                          # S7(b)
        R = rd(E / "V11_subset_regressions.csv"); R = R[R.target == JW]
        rows = []
        for scale in ("SMC", "Large"):
            row = {"sample": f"{SLAB[scale]}, 65+"}
            for lab in ("shrinking", "non-shrinking"):
                r = R[(R.scale == scale) & (R.subset == lab)]
                if len(r):
                    row.update({f"{lab}: coefficient": r.coef_CS.iloc[0], f"{lab}: SE": r.se_cluster.iloc[0], f"{lab}: n": int(r.n.iloc[0])})
            v = fam(family_dir, scale, "lag1")
            if v is not None:
                j = v[v.target == JW].iloc[0]
                row.update({"Caregiving Support (t-1) x rate": j.interaction, "interaction SE": j.se_interaction,
                            "interaction p": j.p_cluster, "interaction n": int(j.n)})
            if len(row) > 1:
                rows.append(row)
        T["S7b_jobs_wages_regressions"] = pd.DataFrame(rows)

    v = fam(family_dir, "SMC", "lag1")
    if v is not None:                                                        # S8(a)
        t = pd.DataFrame({"target (t)": v.target, "Caregiving Support (t-1) x rate": v.interaction, "SE": v.se_interaction,
                          "cluster p": v.p_cluster, "rate P10": v.rate_p10, "coefficient at P10": v.coef_at_p10,
                          "P10 CI low": v.coef_at_p10 - 1.96 * v.se_at_p10, "P10 CI high": v.coef_at_p10 + 1.96 * v.se_at_p10,
                          "rate P90": v.rate_p90, "coefficient at P90": v.coef_at_p90,
                          "P90 CI low": v.coef_at_p90 - 1.96 * v.se_at_p90, "P90 CI high": v.coef_at_p90 + 1.96 * v.se_at_p90,
                          "n": v.n, "municipalities": v.n_munis})
        T["S8a_family_13_targets"] = t.sort_values("Caregiving Support (t-1) x rate", ascending=False)

    if (E / "V11_quintile_slopes.csv").exists():                             # S8(b)
        Q = rd(E / "V11_quintile_slopes.csv"); Q = Q[Q.scale == "SMC"]
        w, j = (Q[Q.target == t].set_index("quintile") for t in (WLB, JW))
        T["S8b_quintile_coefficients"] = pd.DataFrame({
            "quintile": w.index, "median rate (%/yr)": w.rate_median.to_numpy(), "n": w.n.to_numpy(), "municipalities": w.n_munis.to_numpy(),
            "Work-Life Balance": w.slope_CS.to_numpy(), "Work-Life Balance SE": w.se_cluster.to_numpy(),
            "Jobs & Wages": j.slope_CS.reindex(w.index).to_numpy(), "Jobs & Wages SE": j.se_cluster.reindex(w.index).to_numpy()})

    return T


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network-dir", required=True); ap.add_argument("--strain-dir", required=True)
    ap.add_argument("--care-dir", default=None, help="output directory of care_links.py")
    ap.add_argument("--outdir", default="tables")
    ap.add_argument("--xlsx"); ap.add_argument("--data-csv"); ap.add_argument("--pop-csv")
    a = ap.parse_args(argv); N, S, O = Path(a.network_dir), Path(a.strain_dir), Path(a.outdir); O.mkdir(parents=True, exist_ok=True)
    T = {}
    art, en, es = article_tables(N, S)
    t1 = None
    if a.xlsx or (a.data_csv and a.pop_csv):
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import validation_analysis as va
        frame, _, _ = va.load_frame(xlsx=a.xlsx, data_csv=a.data_csv, pop_csv=a.pop_csv)
        t1 = table1(va, frame)
        T["T1_sample_structure"] = t1
        T["T2_descriptives"], T["T2_descriptives_wide"] = table2(va, frame)
        T.update(art)
        T["S_depth_negative_pooled_edges"], T["S_depth_life_satisfaction_hub"] = depth_tables(va, frame, en, es)
    else:
        T.update(art)
        print("[tables] no data inputs: Tables 1-2 and the conditioning depths (Table 5, Table S5) are skipped", flush=True)
    T.update(care_tables(S, a.care_dir))
    readme = pd.DataFrame({"sheet": list(T), "content": [CONTENT.get(k, "") for k in T]})
    with pd.ExcelWriter(O / "tables.xlsx", engine="openpyxl") as xw:
        readme.to_excel(xw, sheet_name="README", index=False)
        for k, df in T.items():
            num = df.select_dtypes("number").columns
            out = df.copy(); out[num] = out[num].round(4)
            out.to_excel(xw, sheet_name=k[:31], index=False); out.to_csv(O / f"{k}.csv", index=False, encoding="utf-8-sig")
    print(f"[tables] {len(T)} tables written to {O.resolve()} (tables.xlsx and CSV files)", flush=True)
    if t1 is not None:
        check_table1(t1)


if __name__ == "__main__":
    main()
