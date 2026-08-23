#!/usr/bin/env python3
"""
complete_pipeline.py
====================
Complete, self-verifying analysis and reporting pipeline for

  Kato, H. "Population decline and the temporal structure of
  well-being: Contextual strain amid between-person stability in
  Japan's shrinking cities." (Cities, in revision)

Workbook -> panel -> calibration -> validation gates -> PCMCI
(pooled + person-demeaned + composition check) -> decline contrasts
-> full-precision exports -> Tables 1-7 -> Figures 2-4.

USAGE
-----
    python complete_pipeline.py --xlsx /path/to/workbook.xlsx
        [--outdir OUT] [--force] [--jitter-seed N] [--skip-figures]

or paste into one notebook cell after setting, in a prior cell,

    XLSX_PATH_OVERRIDE = "/path/to/workbook.xlsx"

DATA ACCESS
-----------
The individual-level workbook is provided by the Cabinet Office of
Japan for approved academic use and is NOT distributed with this
repository. Researchers can apply for access through the Cabinet
Office's designated data-use procedures (see README). All pipeline
outputs are aggregate-level (edge lists, centralities, contrasts,
counts) and safe to share; nothing person-level is written to disk.

VALUE CODINGS
-------------
  - satisfaction items are raw-coded 1-11 in the workbook; the
    pipeline subtracts 1 to obtain the 0-10 scale used in the survey
    documentation and the manuscript
  - city scale: the set of CityType values classified as "Large"
    (designated + core cities) is identified automatically and must
    reproduce the manuscript count of 122 large-city municipalities
    exactly, otherwise the pipeline stops (see calibrate_and_clean)
  - age: 5-year bin codes; the 65+ cohort corresponds to codes >= 11
  - shrinkage: SCiRN criterion, compound annual change rate of TOTAL
    registered population (Japanese + foreign residents, Basic
    Resident Register) over 2017-2023; municipalities with
    r <= -0.15 %/yr are classified as shrinking (manuscript Eq. 1)
  - sample: rows missing the respondent ID or any analysis variable
    are dropped -> 56,206 person-years from 30,474 respondents
  - n_valid: transition-mask-valid lag-1 count (consecutive years
    within the same respondent)

METHODOLOGICAL NOTES FOR RE-USERS
---------------------------------
  1. Jitter. A negligible jitter (SD = 1e-5, far below the coding
     resolution) is added to the seven covariate codes to avoid
     degenerate ties in ParCorr. To reproduce the manuscript's random
     stream EXACTLY, the fixed protocol is: np.random.seed (default
     42), sequential per-column draws in the order JITTER_ORDER,
     applied BEFORE the final (respondent, year) sort. No reported
     estimate depends on this jitter: re-running with a different
     seed (e.g. --jitter-seed 1) reproduces every reported edge,
     out-strength value, and contrast to the printed precision.
  2. PC_ALPHA = 0.05 is the significance level of the PCMCI condition
     -selection (PC) phase. It affects conditioning-set cardinality
     and is reported, with all other parameters, in
     data/run_metadata.json.
  3. Person-demeaning is applied within each estimation subset, so
     for respondents who moved between subsets during the panel the
     person mean is computed from their waves within that subset.
     Transitions that cross subset boundaries are invalidated by the
     transition mask in any case.
  4. The transition mask is applied to the dependent variable
     (mask_type="y"), so MCI estimates and their analytic
     significance are computed exclusively from valid within-person
     consecutive-year transitions (manuscript Section 3.3).
  5. The contemporaneous (lag-0) layer is exported for verification
     but is unoriented and reported only descriptively in the
     manuscript.

OUTPUTS (all aggregate-level; safe to share)
--------------------------------------------
data/
  edges_pooled.csv, edges_within.csv, contrasts.csv    full precision
  all_edges_pooled.csv, all_edges_within.csv           all 182 lagged
                                     hypotheses per model, incl.
                                     non-significant (within file
                                     includes non-shrinking subsets)
  sensitivity_edge_counts.csv        edge counts under 4 criteria
  composition_check.csv              3-wave-plus pooled models
  contemporaneous_pairs.csv          lag-0 verification
  run_metadata.json                  parameters + software versions
tables/
  table1_sample_structure.csv    table2_descriptives.csv
  table3_edge_counts.csv         table4_pooled_edges.csv
  table5_out_strength.csv        table6_within_edges.csv
  table7_contrasts.csv
figures/
  figure2.(png|pdf)  pooled networks, all eight strata (4 x 2 grid)
  figure3.(png|pdf)  person-demeaned within-person networks (65+)
  figure4.(png|pdf)  decline-contrast forest plot

(Figure 1, the municipal map, is produced in GIS and is not part of
this pipeline.)

SELF-VERIFICATION GATES
-----------------------
  Gate 1  sample structure vs manuscript Table 1 (stops on failure)
  Gate 2  out-strength sentinels vs manuscript Table 5
  Gate 3  within-person edge counts vs manuscript Section 4.3 /
          Table 6 and its note (0/0/0/0 for the 15-64 strata;
          0/1/1/4 for the 65+ strata)
  Gate 4  sensitivity criterion: Life-Satisfaction out-strength zero
          in seven of eight demeaned strata (q < .05, |MCI| >= .05)
  Gate 5  composition check: pooled models on the 3-wave-plus
          subsample retain the Life-Satisfaction hub (out-strength
          0.62-0.90) in all eight strata
  Gate 6  contemporaneous layer: all 91 unoriented lag-0 pairs
          significant in every pooled model

Requirements: Python >= 3.10, numpy, pandas, openpyxl, matplotlib,
tigramite (>= 5.2).  Runtime: roughly 30-90 minutes (28 PCMCI runs).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib import cm, colors, gridspec

import tigramite
from tigramite.data_processing import DataFrame as TigDataFrame
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr

__version__ = "1.0.0"

# ======================================================================
# CONFIG (defaults; overridable from the command line)
# ======================================================================
XLSX_PATH = "workbook.xlsx"              # set via --xlsx or override
OUT_ROOT = Path(".")
DATA_DIR, TABLE_DIR, FIG_DIR = Path("data"), Path("tables"), Path("figures")

PC_ALPHA = 0.05                          # PC-phase selection level
JITTER_SD, JITTER_SEED = 1e-5, 42        # fixed protocol (note 1)
Q_THRESHOLD, MCI_FLOOR = 0.05, 0.10      # main reporting criteria
SENS_FLOOR = 0.05                        # loose sensitivity floor
TAU_MIN, TAU_MAX = 0, 1
MIN_WAVES_DEMEANED = 3
FORCE_PROCEED = False
SKIP_FIGURES = False
EXPORT_ALL_EDGES = True                  # full per-model edge lists

# ======================================================================
# Constants
# ======================================================================
SHEET_ANALYSIS, SHEET_POP = "24_分析データ", "17_Population"
POP_BASE_YEAR, POP_END_YEAR, SHRINK_YEARS = 2017, 2023, 6
SHRINK_THRESHOLD = -0.0015
COVID_YEARS = {2020, 2021}
MISSING_FLAG = 99999.0                   # defensive; frame is complete

ANALYSIS_COLS = [
    "PID", "Life Satisfaction", "Finances & Assets", "Jobs & Wages",
    "Housing", "Work-Life Balance", "Health", "Education",
    "Friends & Community", "Trust in Government & Courts",
    "Air & Water (environment)", "Personal Safety", "Parenting Support",
    "Caregiving Support", "Enjoyment of Life",
    "CityName", "CityCode6", "muni_code", "city_type", "Areakm2",
    "merged_flag", "metro_flag", "pop_japanese", "pop_foreign",
    "pop_total", "inflow", "outflow",
    "sex", "age_category", "year",
    "education", "household_income", "household_size",
]
SATISFACTION_VARS = ANALYSIS_COLS[1:15]
N_SAT = len(SATISFACTION_VARS)
N_PAIRS_LAG0 = N_SAT * (N_SAT - 1) // 2  # 91 unoriented pairs
COVARIATES = ["sex", "age_category", "survey_year", "covid_period",
              "education", "household_income", "household_size"]
# jitter order follows the demographic column order of the original
# estimation protocol (note 1); changing it changes nothing reported:
JITTER_ORDER = ["sex", "age_category", "survey_year",
                "education", "household_income", "household_size",
                "covid_period"]
ALL_VARS = SATISFACTION_VARS + COVARIATES

POP_COLS = ["sy_code", "sy_jp", "pop_year", "region_code", "muni_code",
            "region_name", "japanese", "foreign", "population",
            "inflow", "outflow"]


def model_label(scale, cohort65, shrink):
    """Canonical stratum label used in every table, figure, and gate.

    Labels follow the manuscript convention: fixed group names in
    title case ("All Large Cities", "Shrinking SMCs", ...).
    """
    s = {None: "All", True: "Shrinking", False: "Non-shrinking"}[shrink]
    return (f"{s} {'Large Cities' if scale == 'Large' else 'SMCs'}, "
            f"aged {'65+' if cohort65 else '15-64'}")


# Manuscript Table 1 targets (Gate 1).
TABLE1 = {
    "person_years_total": 56206, "respondents_total": 30474,
    "share65": 16.3,
    "Large": dict(munis=122, shrinking_munis=69, person_years=26714,
                  respondents=14837, median_r=-0.26, trans_total=8218,
                  trans_1564=6205, trans_65=1811),
    "SMC": dict(munis=1214, shrinking_munis=996, person_years=29492,
                respondents=16425, median_r=-1.08, trans_total=9181,
                trans_1564=6818, trans_65=2114),
    "ShrinkLarge": dict(person_years=17506, respondents=10097,
                        median_r=-0.46, trans_total=5228,
                        trans_1564=3859, trans_65=1237),
    "ShrinkSMC": dict(person_years=23116, respondents=12967,
                      median_r=-1.32, trans_total=7178,
                      trans_1564=5301, trans_65=1669),
}
# Manuscript Table 5 sentinels (Gate 2).
TABLE5_SENTINELS = {
    (model_label("Large", True, None), "Life Satisfaction"): 0.727,
    (model_label("SMC", False, None), "Life Satisfaction"): 0.673,
    (model_label("Large", False, True), "Life Satisfaction"): 1.117,
    (model_label("Large", True, True), "Finances & Assets"): 0.402,
    (model_label("SMC", True, True), "Life Satisfaction"): 0.851,
}
# Manuscript Section 4.3 / Table 6 and its note (Gate 3): reported
# within-person edge counts across all eight primary strata ("zero
# edges in five of eight strata, one edge in two strata, and four
# edges in one stratum"; no edge in any 15-64 model).
WITHIN_PROFILE = {
    model_label("Large", False, None): 0,
    model_label("SMC", False, None): 0,
    model_label("Large", False, True): 0,
    model_label("SMC", False, True): 0,
    model_label("Large", True, None): 0,
    model_label("SMC", True, None): 1,
    model_label("Large", True, True): 1,
    model_label("SMC", True, True): 4,
}
# Manuscript Section 4.3 composition check (Gate 5): pooled models on
# the 3-wave-plus subsample retain the Life-Satisfaction hub.
COMPOSITION_LS_RANGE = (0.62, 0.90)

POOLED_ORDER = [model_label(sc, c65, sh)
                for sh in (None, True)
                for sc in ("Large", "SMC")
                for c65 in (False, True)]
WITHIN_ORDER_65 = [model_label(sc, True, sh)
                   for sh in (None, True) for sc in ("Large", "SMC")]
WITHIN_ORDER_ALL8 = [model_label(sc, False, sh)
                     for sh in (None, True)
                     for sc in ("Large", "SMC")] + WITHIN_ORDER_65

# Figure 2 panel order: within each city scale, the All benchmark row
# is placed directly above its shrinking subset, so benchmark-subset
# contrasts read vertically (manuscript Figure 2 caption).
FIG2_LABELS = [model_label(sc, c65, sh)
               for sc in ("Large", "SMC")
               for sh in (None, True)
               for c65 in (False, True)]
FIG2_TAGS = ["(a1)", "(a2)", "(b1)", "(b2)",
             "(c1)", "(c2)", "(d1)", "(d2)"]


# ======================================================================
# Stage 1: read workbook (positional column mapping)
# ======================================================================
def _clean_code(s):
    return (s.astype("string").str.strip()
            .str.replace(r"\.0$", "", regex=True).str.zfill(5))


def read_workbook():
    df = pd.read_excel(XLSX_PATH, sheet_name=SHEET_ANALYSIS, header=0)
    if df.shape[1] != len(ANALYSIS_COLS):
        raise ValueError(f"'{SHEET_ANALYSIS}': {df.shape[1]} cols, "
                         f"expected {len(ANALYSIS_COLS)}")
    df.columns = ANALYSIS_COLS
    df["muni_code"] = _clean_code(df["muni_code"])
    pop = pd.read_excel(XLSX_PATH, sheet_name=SHEET_POP, header=0)
    if pop.shape[1] != len(POP_COLS):
        raise ValueError(f"'{SHEET_POP}': {pop.shape[1]} cols, "
                         f"expected {len(POP_COLS)}")
    pop.columns = POP_COLS
    pop["muni_code"] = _clean_code(pop["muni_code"])
    pop["pop_year"] = pd.to_numeric(pop["pop_year"], errors="coerce")
    print(f"[read] '{SHEET_ANALYSIS}': {len(df)} rows | "
          f"'{SHEET_POP}': {len(pop)} rows")
    return df, pop


# ======================================================================
# Stage 2: calibration + cleaning
# ======================================================================
def shrinkage_table(pop, count_col="population"):
    """SCiRN shrinkage from TOTAL registered population (Eq. 1)."""
    base = (pop[pop["pop_year"] == POP_BASE_YEAR]
            .drop_duplicates("muni_code")
            .set_index("muni_code")[count_col].astype(float))
    end = (pop[pop["pop_year"] == POP_END_YEAR]
           .drop_duplicates("muni_code")
           .set_index("muni_code")[count_col].astype(float))
    t = pd.DataFrame({"p0": base, "p1": end}).dropna()
    t = t[(t["p0"] > 0) & (t["p1"] > 0)]
    t["annual_change_rate"] = ((t["p1"] / t["p0"])
                               ** (1.0 / SHRINK_YEARS) - 1.0)
    t["shrinking"] = (t["annual_change_rate"]
                      <= SHRINK_THRESHOLD).astype(int)
    return t.reset_index()[["muni_code", "annual_change_rate",
                            "shrinking"]]


def identify_large_city_values(muni_type):
    """Identify the CityType value set classified as Large.

    The set must reproduce the manuscript count of large-city
    municipalities EXACTLY; any residual error stops the pipeline
    rather than silently proceeding with a near-match.
    """
    values = sorted(muni_type["city_type"].dropna().unique().tolist())
    best, best_err = None, None
    for r in range(1, len(values)):
        for S in itertools.combinations(values, r):
            err = abs(muni_type["city_type"].isin(S).sum()
                      - TABLE1["Large"]["munis"])
            if best_err is None or err < best_err:
                best, best_err = set(S), err
    if best_err != 0:
        raise ValueError(
            f"CityType calibration failed: best candidate {sorted(best)} "
            f"misses the manuscript large-city count by {best_err}. "
            "Check the workbook's CityType coding.")
    return best


def calibrate_and_clean(df, pop):
    print("[calibrate] value codings")
    for c in SATISFACTION_VARS + ["sex", "age_category", "education",
                                  "household_income", "household_size"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["year"] = pd.to_numeric(df["year"], errors="coerce")

    # satisfaction: raw 1-11 -> 0-10 scale
    for c in SATISFACTION_VARS:
        df[c] = df[c] - 1.0
        out_of_range = int(((df[c] < 0) | (df[c] > 10)).sum())
        if out_of_range:
            print(f"  [warn] {c}: {out_of_range} values outside 0-10 "
                  "after shift -> set NaN")
            df[c] = df[c].where((df[c] >= 0) & (df[c] <= 10))

    # (a) CityType -> Large set (exact match enforced)
    muni_type = df.dropna(subset=["muni_code"]).drop_duplicates(
        "muni_code")[["muni_code", "city_type"]]
    large_values = identify_large_city_values(muni_type)
    print(f"  (a) CityType Large values = {sorted(large_values)}")
    df["city_scale"] = np.where(df["city_type"].isin(large_values),
                                "Large", "SMC")
    df.attrs["large_city_values"] = sorted(map(str, large_values))

    # (b) shrinkage from total registered population (Eq. 1)
    df = df.merge(shrinkage_table(pop, "population"), on="muni_code",
                  how="left")

    # (c) age cohort: 65+ <=> bin code > 10
    df["age65plus"] = (df["age_category"] > 10).astype("Int64")

    # covariate completions
    df["survey_year"] = df["year"]
    df["covid_period"] = df["year"].isin(COVID_YEARS).astype(float)

    # sample definition: complete-case rule (manuscript Section 3.2) --
    # respondent ID plus all analysis variables
    need = (["PID", "year", "age_category", "city_type", "shrinking"]
            + SATISFACTION_VARS
            + ["sex", "education", "household_income", "household_size"])
    frame = df.dropna(subset=need).copy()
    frame.attrs["large_city_values"] = df.attrs["large_city_values"]
    print(f"  complete-case frame: {len(frame)} person-years "
          f"(target {TABLE1['person_years_total']})")

    # jitter: fixed protocol (seed, sequential draws, PRE-sort order);
    # see module docstring, note 1 -- no reported estimate depends on it
    np.random.seed(JITTER_SEED)
    for c in JITTER_ORDER:
        frame[c] = frame[c].astype(float) + np.random.normal(
            0.0, JITTER_SD, size=len(frame))
    frame = frame.rename(columns={"PID": "respondent_id"})
    frame["shrinking"] = frame["shrinking"].astype(int)
    frame["age65plus"] = frame["age65plus"].astype(int)
    frame = frame.sort_values(["respondent_id", "year"]).reset_index(
        drop=True)
    return frame


# ======================================================================
# Stage 3: Gate 1 (Table 1) + table1 export
# ======================================================================
def count_valid(sub) -> int:
    ids = sub["respondent_id"].to_numpy()
    yrs = np.round(sub["year"].to_numpy())
    return int(((ids[1:] == ids[:-1]) & (yrs[1:] == yrs[:-1] + 1)).sum())


def stratum(df, scale=None, cohort65=None, shrink=None):
    m = pd.Series(True, index=df.index)
    if scale is not None:
        m &= df["city_scale"] == scale
    if cohort65 is not None:
        m &= df["age65plus"] == int(cohort65)
    if shrink is not None:
        m &= df["shrinking"] == int(shrink)
    return df[m]


def gate1(df) -> bool:
    ok = True
    rows = []

    def check(name, got, exp, tol=0):
        nonlocal ok
        good = abs(got - exp) <= tol
        ok = ok and good
        rows.append({"item": name, "pipeline": got, "manuscript": exp,
                     "match": "OK" if good else "CHECK"})
        print(f"  {name:52s} got {got:>10} | expected {exp:>10} "
              f"[{'OK' if good else 'CHECK'}]")

    print("[gate 1] manuscript Table 1:")
    check("person-years (total)", len(df), TABLE1["person_years_total"])
    check("unique respondents (total)",
          df["respondent_id"].nunique(), TABLE1["respondents_total"])
    for scale in ("Large", "SMC"):
        for shrink, key in ((None, scale), (True, "Shrink" + scale)):
            sub = stratum(df, scale, None, shrink)
            t = TABLE1[key]
            tag = ("" if shrink is None else "Shrinking ") + scale
            if shrink is None:
                check(f"{tag}: municipalities",
                      sub["muni_code"].nunique(), t["munis"])
                check(f"{tag}: shrinking municipalities",
                      sub.loc[sub["shrinking"] == 1,
                              "muni_code"].nunique(),
                      t["shrinking_munis"])
            check(f"{tag}: person-years", len(sub), t["person_years"])
            check(f"{tag}: unique respondents",
                  sub["respondent_id"].nunique(), t["respondents"])
            med = sub.groupby("muni_code")["annual_change_rate"] \
                     .first().median()
            check(f"{tag}: median annual change (%/yr)",
                  round(100 * med, 2), t["median_r"], tol=0.005)
            check(f"{tag}: valid lag-1 transitions (total)",
                  count_valid(sub), t["trans_total"])
            check(f"{tag}: valid transitions (15-64)",
                  count_valid(stratum(df, scale, False, shrink)),
                  t["trans_1564"])
            check(f"{tag}: valid transitions (65+)",
                  count_valid(stratum(df, scale, True, shrink)),
                  t["trans_65"])
    print(f"  aged 65+ share: {100 * df['age65plus'].mean():.1f}% "
          f"(manuscript {TABLE1['share65']}%)")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TABLE_DIR / "table1_sample_structure.csv",
                              index=False)
    return ok


# ======================================================================
# Stage 4: estimation primitives
# ======================================================================
def transition_mask(sub):
    ids = sub["respondent_id"].to_numpy()
    yrs = np.round(sub["year"].to_numpy())
    valid = np.zeros(len(sub), dtype=bool)
    valid[1:] = (ids[1:] == ids[:-1]) & (yrs[1:] == yrs[:-1] + 1)
    return np.repeat((~valid)[:, None], len(ALL_VARS), axis=1)


def run_model(sub):
    """Returns (val, p, n_valid) with lag-0..1 matrices over ALL_VARS."""
    data = sub[ALL_VARS].to_numpy(dtype=float)
    mask = transition_mask(sub)
    n_valid = int((~mask[1:, 0]).sum())
    if len(sub) < 25 or n_valid < 10:
        shape = (len(ALL_VARS), len(ALL_VARS), TAU_MAX + 1)
        print(f"    [skip] stratum too small (rows={len(sub)}, "
              f"n_valid={n_valid})")
        return np.full(shape, np.nan), np.full(shape, np.nan), n_valid
    data = np.where(np.isnan(data), MISSING_FLAG, data)
    tdf = TigDataFrame(data, mask=mask, var_names=ALL_VARS,
                       missing_flag=MISSING_FLAG)
    pcmci = PCMCI(dataframe=tdf,
                  cond_ind_test=ParCorr(significance="analytic",
                                        mask_type="y"),
                  verbosity=0)
    res = pcmci.run_pcmci(tau_min=TAU_MIN, tau_max=TAU_MAX,
                          pc_alpha=PC_ALPHA)
    return res["val_matrix"], res["p_matrix"], n_valid


def bh_qvalues_lag1(p_matrix):
    """BH q-values over the 182 directed lagged non-self hypotheses."""
    ps, idx = [], []
    for i in range(N_SAT):
        for j in range(N_SAT):
            if i != j:
                ps.append(p_matrix[i, j, 1])
                idx.append((i, j))
    ps = np.asarray(ps)
    m = len(ps)
    order = np.argsort(ps)
    q_sorted = ps[order] * m / np.arange(1, m + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(q_sorted, 1.0)
    out = np.full((N_SAT, N_SAT), np.nan)
    for (i, j), qi in zip(idx, q):
        out[i, j] = qi
    return out


def reported_edges(val, p, q, floor=None):
    floor = MCI_FLOOR if floor is None else floor
    edges = []
    for i in range(N_SAT):
        for j in range(N_SAT):
            if i != j and not np.isnan(q[i, j]) \
                    and q[i, j] < Q_THRESHOLD \
                    and abs(val[i, j, 1]) >= floor:
                edges.append((SATISFACTION_VARS[i], SATISFACTION_VARS[j],
                              float(val[i, j, 1]), float(p[i, j, 1]),
                              float(q[i, j])))
    return sorted(edges, key=lambda e: -abs(e[2]))


def out_strength(edges):
    s = {}
    for (u, _v, mci, *_r) in edges:
        s[u] = s.get(u, 0.0) + abs(mci)
    return s


def all_edge_rows(label, val, p, q, n_valid):
    """All 182 directed lagged hypotheses for one model."""
    rows = []
    for i in range(N_SAT):
        for j in range(N_SAT):
            if i == j:
                continue
            rows.append({"model": label,
                         "source_t_minus_1": SATISFACTION_VARS[i],
                         "target_t": SATISFACTION_VARS[j],
                         "mci": val[i, j, 1], "p": p[i, j, 1],
                         "q": q[i, j], "effective_n": n_valid})
    return rows


def criterion_counts(val, p, q):
    """Edge counts under the four sensitivity criteria."""
    counts = {"p05_mci05": 0, "q05_mci05": 0,
              "p05_mci10": 0, "q05_mci10": 0}
    for i in range(N_SAT):
        for j in range(N_SAT):
            if i == j:
                continue
            a = abs(val[i, j, 1])
            pp, qq = p[i, j, 1], q[i, j]
            if np.isnan(pp) or np.isnan(qq):
                continue
            if pp < .05 and a >= .05:
                counts["p05_mci05"] += 1
            if qq < .05 and a >= .05:
                counts["q05_mci05"] += 1
            if pp < .05 and a >= .10:
                counts["p05_mci10"] += 1
            if qq < .05 and a >= .10:
                counts["q05_mci10"] += 1
    return counts


def contemporaneous_summary(val, p):
    """Verification of the lag-0 layer over the 91 unoriented pairs.

    Pair p-value = max of the two directed lag-0 p-values
    (conservative). Returns counts under uncorrected p < .05 and the
    maximum BH q over the 91 pair hypotheses.
    """
    pair_ps = []
    for i in range(N_SAT):
        for j in range(i + 1, N_SAT):
            pair_ps.append(max(p[i, j, 0], p[j, i, 0]))
    pair_ps = np.asarray(pair_ps)
    m = len(pair_ps)
    order = np.argsort(pair_ps)
    q_sorted = pair_ps[order] * m / np.arange(1, m + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    return {"n_pairs": m,
            "n_p_lt_05": int((pair_ps < 0.05).sum()),
            "max_p": float(pair_ps.max()),
            "max_q_bh": float(min(q_sorted.max(), 1.0))}


def restrict_min_waves(sub):
    counts = sub.groupby("respondent_id")["year"].transform("size")
    return sub[counts >= MIN_WAVES_DEMEANED].copy()


def demean(sub):
    """Person-demean within the estimation subset (docstring, note 3)."""
    out = restrict_min_waves(sub)
    for c in SATISFACTION_VARS:
        out[c] = out[c] - out.groupby("respondent_id")[c].transform("mean")
    return out.sort_values(["respondent_id", "year"]).reset_index(drop=True)


def fisher_z(r1, n1, r2, n2):
    if any(map(math.isnan, (r1, r2))) or min(n1, n2) < 4:
        return float("nan"), float("nan")
    z = (math.atanh(r1) - math.atanh(r2)) / math.sqrt(
        1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    return z, math.erfc(abs(z) / math.sqrt(2.0))


# ======================================================================
# Stage 5: tables
# ======================================================================
def write_table2(frame):
    rows = []
    for scale in ("Large", "SMC"):
        for shrink in (None, True):
            for c65 in (False, True):
                sub = stratum(frame, scale, c65, shrink)
                grp = ("All" if shrink is None else "Shrinking") + \
                      (" Large Cities" if scale == "Large" else " SMCs")
                for c in SATISFACTION_VARS:
                    pm = sub.groupby("respondent_id")[c].mean()
                    within = sub[c] - sub.groupby(
                        "respondent_id")[c].transform("mean")
                    rows.append({
                        "city_group": grp,
                        "age_cohort": "65+" if c65 else "15-64",
                        "domain": c,
                        "mean": round(sub[c].mean(), 2),
                        "sd_overall": round(sub[c].std(ddof=1), 2),
                        "sd_between": round(pm.std(ddof=1), 2),
                        "sd_within": round(within.std(ddof=1), 2)})
    pd.DataFrame(rows).to_csv(TABLE_DIR / "table2_descriptives.csv",
                              index=False)


def write_metadata(input_name, large_city_values):
    meta = {
        "pipeline_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_workbook": input_name,
        "parameters": {
            "pc_alpha": PC_ALPHA, "tau_min": TAU_MIN, "tau_max": TAU_MAX,
            "q_threshold": Q_THRESHOLD, "mci_floor": MCI_FLOOR,
            "sensitivity_floor": SENS_FLOOR,
            "bh_family_size_lag1": N_SAT * (N_SAT - 1),
            "min_waves_demeaned": MIN_WAVES_DEMEANED,
            "jitter_sd": JITTER_SD, "jitter_seed": JITTER_SEED,
            "jitter_order": JITTER_ORDER,
            "shrink_threshold": SHRINK_THRESHOLD,
            "shrink_years": [POP_BASE_YEAR, POP_END_YEAR],
            "shrink_population": "total registered (Japanese+foreign)",
            "mask_type": "y", "cond_ind_test": "ParCorr(analytic)",
            "large_city_citytype_values": large_city_values,
        },
        "software": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "tigramite": getattr(tigramite, "__version__", "unknown"),
        },
    }
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(DATA_DIR / "run_metadata.json", "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print("[meta] data/run_metadata.json written")


# ======================================================================
# Stage 6: figures (hub-centred common layout)
# ======================================================================
HUB = "Life Satisfaction"
WRAPPED = {
    "Life Satisfaction": "Life\nSatisfaction",
    "Finances & Assets": "Finances\n& Assets",
    "Jobs & Wages": "Jobs &\nWages", "Housing": "Housing",
    "Work-Life Balance": "Work-Life\nBalance", "Health": "Health",
    "Education": "Education",
    "Friends & Community": "Friends &\nCommunity",
    "Trust in Government & Courts": "Trust in Gov.\n& Courts",
    "Air & Water (environment)": "Air & Water\n(environment)",
    "Personal Safety": "Personal\nSafety",
    "Parenting Support": "Parenting\nSupport",
    "Caregiving Support": "Caregiving\nSupport",
    "Enjoyment of Life": "Enjoyment\nof Life",
}
SCALE_MAX, WIDTH_PER_MCI = 1.2, 17.0
NODE_R0, NODE_R1 = 0.085, 0.075
CURV_RECIP, CURV_AVOID, AVOID_DIST = 0.18, 0.16, 0.17
CMAP, NORM = cm.Oranges, colors.Normalize(vmin=0.0, vmax=1.2)
COL_POS, COL_NEG, COL_EMPTY = "#4a4a4a", "#c62828", "#f2f2f2"
PANEL_W, PANEL_H, LEGEND_H = 4.3, 4.55, 1.15


def ring_order_2opt(pair_w, initial):
    """Deterministic 2-opt descent minimising |MCI|-weighted angular
    distance; gives the common ring layout shared by all panels."""
    order, n = list(initial), len(initial)

    def cost(o):
        pos = {node: i for i, node in enumerate(o)}
        return sum(w * min(abs(pos[a] - pos[b]),
                           n - abs(pos[a] - pos[b]))
                   for (a, b), w in pair_w.items())

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


def build_positions(pooled_edges, within_edges):
    pair_w = {}
    for edges in list(pooled_edges.values()) + list(within_edges.values()):
        for (u, v, mci, *_r) in edges:
            if HUB in (u, v):
                continue
            k = tuple(sorted((u, v)))
            pair_w[k] = max(pair_w.get(k, 0.0), abs(mci))
    ring = ring_order_2opt(pair_w, [d for d in SATISFACTION_VARS
                                    if d != HUB])
    pos = {HUB: np.array([0.0, 0.0])}
    for i, node in enumerate(ring):
        ang = np.deg2rad(90.0 - i * 360.0 / len(ring))
        pos[node] = np.array([np.cos(ang), np.sin(ang)])
    print("[layout] ring order:", " | ".join(ring))
    return pos


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


def draw_panel(ax, edges, strengths, title, pos, empty_note=None):
    ax.set_xlim(-1.62, 1.62)
    ax.set_ylim(-1.60, 1.55)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=9, fontweight="bold", pad=4)
    patches = {}
    for node, xy in pos.items():
        s = strengths.get(node, 0.0)
        c = Circle(xy, node_radius(s),
                   facecolor=CMAP(NORM(s)) if s > 0 else COL_EMPTY,
                   edgecolor="#888888", linewidth=0.7, zorder=3)
        ax.add_patch(c)
        patches[node] = c
    for node, xy in pos.items():
        if node == HUB:
            ax.text(xy[0], xy[1] - node_radius(strengths.get(node, 0))
                    - 0.055, WRAPPED[node], ha="center", va="top",
                    fontsize=6.4, zorder=5)
        else:
            u = xy / np.linalg.norm(xy)
            lx, ly = xy + u * 0.30
            ha = "left" if lx > 0.45 else "right" if lx < -0.45 \
                else "center"
            va = "bottom" if ly > 0.85 else "top" if ly < -0.85 \
                else "center"
            ax.text(lx, ly, WRAPPED[node], ha=ha, va=va, fontsize=6.4,
                    zorder=5)
    pairset = {(u, v) for (u, v, *_r) in edges}
    for (u, v, mci, *_r) in sorted(edges, key=lambda e: abs(e[2])):
        neg = mci < 0
        ax.add_patch(FancyArrowPatch(
            pos[u], pos[v], patchA=patches[u], patchB=patches[v],
            connectionstyle=f"arc3,rad={edge_rad(u, v, pairset, pos)}",
            arrowstyle="-|>", mutation_scale=9,
            linewidth=abs(mci) * WIDTH_PER_MCI,
            color=COL_NEG if neg else COL_POS,
            linestyle=(0, (4, 2)) if neg else "solid",
            shrinkA=0, shrinkB=0, zorder=2, alpha=0.9))
    if empty_note:
        ax.text(0, -1.50, empty_note, ha="center", va="center",
                fontsize=7, style="italic", color="#666666")


def add_common_legend(fig, fig_h):
    bottom = 0.42 / fig_h
    cax = fig.add_axes([0.13, bottom, 0.34, 0.17 / fig_h])
    cb = fig.colorbar(cm.ScalarMappable(norm=NORM, cmap=CMAP),
                      cax=cax, orientation="horizontal")
    cb.set_ticks([0.0, 0.3, 0.6, 0.9, 1.2])
    cb.ax.tick_params(labelsize=6.5)
    cax.set_title("Out-strength centrality (common scale, Table 5)",
                  fontsize=7, pad=3)
    handles = [Line2D([], [], color=COL_POS, lw=w * WIDTH_PER_MCI,
                      label=f"|MCI| = {w:.2f}") for w in (.10, .15, .20)]
    handles.append(Line2D([], [], color=COL_NEG,
                          lw=0.12 * WIDTH_PER_MCI,
                          linestyle=(0, (4, 2)),
                          label="negative dependency"))
    fig.legend(handles=handles, loc="lower right",
               bbox_to_anchor=(0.95, bottom - 0.10 / fig_h),
               fontsize=6.8, frameon=False, ncol=2,
               title="Lagged dependency strength", title_fontsize=7)


def panel_title(tag, label):
    strat, age = label.rsplit(", aged ", 1)
    cohort = "Older adults aged 65+" if age == "65+" \
        else "Respondents aged 15-64"
    return f"{tag} {strat}\n({cohort})"


def network_figure(edge_dict, labels, tags, pos, fname,
                   empty_note_on_zero=False):
    """Generic n-panel (2-column) network figure."""
    nrow = (len(labels) + 1) // 2
    fig_h = nrow * PANEL_H + LEGEND_H
    fig = plt.figure(figsize=(2 * PANEL_W, fig_h))
    gs = gridspec.GridSpec(nrow, 2, left=0.03, right=0.97,
                           top=1.0 - 0.28 / fig_h,
                           bottom=LEGEND_H / fig_h,
                           wspace=0.05, hspace=0.14)
    for k, lab in enumerate(labels):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        edges = edge_dict.get(lab, [])
        note = ("No edges satisfy the reporting criteria"
                if empty_note_on_zero and not edges else None)
        draw_panel(ax, edges, out_strength(edges),
                   panel_title(tags[k], lab), pos, empty_note=note)
    add_common_legend(fig, fig_h)
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{fname}.{ext}", dpi=300)
    plt.close(fig)
    print(f"[figure] {fname} written ({nrow}x2 panels)")


def forest_figure(contrasts, fname="figure4"):
    fig, axf = plt.subplots(figsize=(8.6, 3.6))
    fig.subplots_adjust(left=0.24, right=0.96, top=0.90, bottom=0.18)
    y = np.arange(len(contrasts))[::-1]
    for i, r in enumerate(contrasts):
        yy = y[i]
        grey = "specificity" in r["stratum"]
        for rr, n, off, filled in [
                (r["shrinking_mci"], r["shrinking_n"], +0.14, True),
                (r["nonshrinking_mci"], r["nonshrinking_n"], -0.14,
                 False)]:
            if math.isnan(rr) or n < 4:
                continue
            se = 1.0 / math.sqrt(n - 3)
            lo = math.tanh(math.atanh(rr) - 1.959964 * se)
            hi = math.tanh(math.atanh(rr) + 1.959964 * se)
            col = "#9e9e9e" if grey else (COL_NEG if filled
                                          else "#1565c0")
            axf.plot([lo, hi], [yy + off, yy + off], color=col, lw=1.4)
            axf.plot(rr, yy + off, marker="o",
                     markerfacecolor=col if filled else "white",
                     markeredgecolor=col, markersize=5.5, zorder=3)
        if math.isnan(r["z"]):
            continue
        ptxt = ("p < 0.001" if r["p"] < 0.001 else f"p = {r['p']:.3f}")
        axf.text(0.30, yy, f"z = {r['z']:.2f}, {ptxt}", fontsize=7,
                 va="center", color="#555555")
    axf.axvline(0, color="#888888", lw=0.8, linestyle=":")
    axf.set_yticks(y)
    axf.set_yticklabels([r.get("label", r["stratum"])
                         for r in contrasts], fontsize=7)
    axf.set_xlim(-0.30, 0.52)
    axf.set_xlabel("Within-person MCI: Caregiving Support (t-1) "
                   "\u2192 Work-Life Balance (t)", fontsize=8)
    axf.tick_params(axis="x", labelsize=7)
    for sp in ("top", "right"):
        axf.spines[sp].set_visible(False)
    axf.legend(handles=[
        Line2D([], [], marker="o", color=COL_NEG, lw=0,
               markerfacecolor=COL_NEG, label="Shrinking"),
        Line2D([], [], marker="o", color="#1565c0", lw=0,
               markerfacecolor="white", label="Non-shrinking")],
        fontsize=7, frameon=False, loc="upper left")
    axf.set_title("Decline contrasts; bars are approximate 95% CIs "
                  "(Fisher z, conditioning-set df not reflected)",
                  fontsize=8, loc="left")
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{fname}.{ext}", dpi=300)
    plt.close(fig)
    print(f"[figure] {fname} written")


# ======================================================================
# Main
# ======================================================================
def parse_args(argv):
    ap = argparse.ArgumentParser(
        description="Reproduction pipeline (see module docstring).")
    ap.add_argument("--xlsx", default=None,
                    help="path to the Cabinet Office workbook")
    ap.add_argument("--outdir", default=".",
                    help="root directory for data/, tables/, figures/")
    ap.add_argument("--force", action="store_true",
                    help="proceed past a failed Gate 1")
    ap.add_argument("--jitter-seed", type=int, default=None,
                    help="alternative jitter seed (jitter-invariance "
                         "check; see docstring, note 1)")
    ap.add_argument("--skip-figures", action="store_true",
                    help="skip figure rendering")
    return ap.parse_args(argv)


def main(argv=None):
    global XLSX_PATH, DATA_DIR, TABLE_DIR, FIG_DIR
    global FORCE_PROCEED, JITTER_SEED, SKIP_FIGURES
    if argv is not None or not hasattr(sys, "ps1"):
        args = parse_args(argv if argv is not None else sys.argv[1:])
        if args.xlsx:
            XLSX_PATH = args.xlsx
        root = Path(args.outdir)
        DATA_DIR, TABLE_DIR, FIG_DIR = (root / "data", root / "tables",
                                        root / "figures")
        FORCE_PROCEED = FORCE_PROCEED or args.force
        SKIP_FIGURES = SKIP_FIGURES or args.skip_figures
        if args.jitter_seed is not None:
            JITTER_SEED = args.jitter_seed
    if "XLSX_PATH_OVERRIDE" in globals():
        XLSX_PATH = globals()["XLSX_PATH_OVERRIDE"]

    t0 = time.time()
    for d in (DATA_DIR, TABLE_DIR, FIG_DIR):
        d.mkdir(parents=True, exist_ok=True)
    df, pop = read_workbook()
    frame = calibrate_and_clean(df, pop)
    write_metadata(Path(XLSX_PATH).name,
                   frame.attrs.get("large_city_values", []))
    if not gate1(frame) and not FORCE_PROCEED:
        print("[stop] Gate 1 FAILED; estimation not run "
              "(use --force to override).")
        return
    write_table2(frame)

    # ---- pooled (8 primary strata) -----------------------------------
    pooled, contemp_rows = {}, []
    all_pooled_rows, sens_rows = [], []
    print("[estimate] pooled models")
    for scale in ("Large", "SMC"):
        for shrink in (None, True):
            for c65 in (False, True):
                lab = model_label(scale, c65, shrink)
                sub = stratum(frame, scale, c65, shrink)
                val, p, n_valid = run_model(sub)
                q = bh_qvalues_lag1(p)
                pooled[lab] = (reported_edges(val, p, q), n_valid)
                contemp_rows.append({"model": lab,
                                     **contemporaneous_summary(val, p)})
                sens_rows.append({"model": lab,
                                  **criterion_counts(val, p, q)})
                if EXPORT_ALL_EDGES:
                    all_pooled_rows += all_edge_rows(lab, val, p, q,
                                                     n_valid)
                key = ("Shrink" if shrink else "") + \
                      ("Large" if scale == "Large" else "SMC")
                exp = TABLE1[key]["trans_65" if c65 else "trans_1564"]
                print(f"  {lab}: n_valid={n_valid} (Table 1: {exp}) "
                      f"[{'OK' if n_valid == exp else 'CHECK'}], "
                      f"{len(pooled[lab][0])} edges "
                      f"({time.time()-t0:.0f}s)")

    print("[gate 2] out-strength sentinels vs Table 5:")
    for (lab, dom), exp in TABLE5_SENTINELS.items():
        s = out_strength(pooled[lab][0]).get(dom, 0.0)
        print(f"  {lab} | {dom}: {s:.3f} vs {exp:.3f} "
              f"[{'OK' if abs(s - exp) < 5e-4 else 'CHECK'}]")

    print("[gate 6] contemporaneous layer (91 unoriented pairs):")
    for r in contemp_rows:
        ok = r["n_p_lt_05"] == N_PAIRS_LAG0
        print(f"  {r['model']}: {r['n_p_lt_05']}/{r['n_pairs']} "
              f"significant (max p = {r['max_p']:.2e}, "
              f"max BH q = {r['max_q_bh']:.2e}) "
              f"[{'OK' if ok else 'CHECK'}]")
    pd.DataFrame(contemp_rows).to_csv(
        DATA_DIR / "contemporaneous_pairs.csv", index=False)
    pd.DataFrame(sens_rows).to_csv(
        DATA_DIR / "sensitivity_edge_counts.csv", index=False)

    # ---- within-person (all 8 primary strata + non-shrinking) --------
    print("[estimate] person-demeaned models")
    within, cache, all_within_rows = {}, {}, []
    for scale in ("Large", "SMC"):
        for shrink in (None, True, False):
            for c65 in (False, True):
                lab = model_label(scale, c65, shrink)
                sub = demean(stratum(frame, scale, c65, shrink))
                val, p, n_valid = run_model(sub)
                q = bh_qvalues_lag1(p)
                cache[(scale, c65, shrink)] = (val, n_valid)
                if EXPORT_ALL_EDGES:
                    all_within_rows += all_edge_rows(lab, val, p, q,
                                                     n_valid)
                if shrink is not None and shrink is False:
                    continue          # non-shrinking: contrasts only
                within[lab] = (reported_edges(val, p, q), n_valid)
                # sensitivity: LS out-strength under the loose floor
                loose = reported_edges(val, p, q, floor=SENS_FLOOR)
                within[lab] += (out_strength(loose).get(HUB, 0.0),)
                exp = WITHIN_PROFILE[lab]
                print(f"  {lab}: n_valid={n_valid}, "
                      f"{len(within[lab][0])} edges (Sec 4.3: {exp}) "
                      f"[{'OK' if len(within[lab][0]) == exp else 'CHECK'}]"
                      f" ({time.time()-t0:.0f}s)")

    print("[gate 3] within-person edge counts vs Sec 4.3 / Table 6:")
    n_zero = sum(1 for lab in WITHIN_ORDER_ALL8
                 if len(within[lab][0]) == 0)
    print(f"  strata with zero edges: {n_zero}/8 (manuscript: 5) "
          f"[{'OK' if n_zero == 5 else 'CHECK'}]")

    print("[gate 4] sensitivity (q<.05, |MCI|>=.05): "
          "Life-Satisfaction out-strength")
    n_ls_zero = sum(1 for lab in WITHIN_ORDER_ALL8
                    if within[lab][2] == 0.0)
    for lab in WITHIN_ORDER_ALL8:
        print(f"  {lab}: LS out-strength (loose) = {within[lab][2]:.3f}")
    print(f"  zero in {n_ls_zero}/8 strata (manuscript: 7) "
          f"[{'OK' if n_ls_zero == 7 else 'CHECK'}]")

    # ---- composition check (pooled spec, 3-wave-plus subsample) ------
    print("[gate 5] composition check (pooled, 3-wave-plus subsample):")
    comp_rows = []
    for scale in ("Large", "SMC"):
        for shrink in (None, True):
            for c65 in (False, True):
                lab = model_label(scale, c65, shrink)
                sub = restrict_min_waves(stratum(frame, scale, c65,
                                                 shrink))
                val, p, n_valid = run_model(sub)
                q = bh_qvalues_lag1(p)
                edges = reported_edges(val, p, q)
                s = out_strength(edges)
                ls = s.get(HUB, 0.0)
                hub_ok = ls > 0 and ls == max(s.values(), default=0.0)
                lo, hi = COMPOSITION_LS_RANGE
                in_range = lo - 0.005 <= ls <= hi + 0.005
                comp_rows.append({"model": lab, "n_valid": n_valid,
                                  "n_edges": len(edges),
                                  "ls_out_strength": f"{ls:.3f}",
                                  "ls_leading_hub": hub_ok})
                print(f"  {lab}: LS out-strength {ls:.3f} "
                      f"(range {lo}-{hi}) hub={hub_ok} "
                      f"[{'OK' if in_range and hub_ok else 'CHECK'}] "
                      f"({time.time()-t0:.0f}s)")
    pd.DataFrame(comp_rows).to_csv(
        DATA_DIR / "composition_check.csv", index=False)

    def coef(scale, c65, shrink, src, tgt):
        val, n = cache[(scale, c65, shrink)]
        i, j = SATISFACTION_VARS.index(src), SATISFACTION_VARS.index(tgt)
        return float(val[i, j, 1]), n

    # ---- contrasts (targeted family m = 3 + specificity) -------------
    slab = {"SMC": "SMCs", "Large": "Large cities"}
    contrasts, m_family = [], 3
    for scale, src, tgt in [
            ("SMC", "Caregiving Support", "Work-Life Balance"),
            ("Large", "Caregiving Support", "Work-Life Balance"),
            ("SMC", "Caregiving Support", "Jobs & Wages")]:
        r_s, n_s = coef(scale, True, True, src, tgt)
        r_n, n_n = coef(scale, True, False, src, tgt)
        z, pv = fisher_z(r_s, n_s, r_n, n_n)
        contrasts.append({
            "stratum": f"{slab[scale]}, aged 65+" + (
                "" if tgt == "Work-Life Balance" else " (CS->JW)"),
            "label": f"{slab[scale]}\n(Respondents aged 65+)",
            "edge": f"{src} -> {tgt}",
            "shrinking_mci": r_s, "shrinking_n": n_s,
            "nonshrinking_mci": r_n, "nonshrinking_n": n_n,
            "z": z, "p": pv,
            "p_bonferroni": min(1.0, pv * m_family)})
        print(f"[contrast] {slab[scale]} 65+ {src}->{tgt}: "
              f"{r_s:.3f} vs {r_n:.3f}; z={z:.2f}, p={pv:.4f}")
    for scale in ("SMC", "Large"):
        r_s, n_s = coef(scale, False, True, "Caregiving Support",
                        "Work-Life Balance")
        r_n, n_n = coef(scale, False, False, "Caregiving Support",
                        "Work-Life Balance")
        z, pv = fisher_z(r_s, n_s, r_n, n_n)
        contrasts.append({
            "stratum": f"{slab[scale]}, aged 15-64 (specificity)",
            "label": f"{slab[scale]}\n(Respondents aged 15-64)\n"
                     "(specificity)",
            "edge": "Caregiving Support -> Work-Life Balance",
            "shrinking_mci": r_s, "shrinking_n": n_s,
            "nonshrinking_mci": r_n, "nonshrinking_n": n_n,
            "z": z, "p": pv, "p_bonferroni": float("nan")})

    # ---- data/ exports (full precision) ------------------------------
    with open(DATA_DIR / "edges_pooled.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "source", "target", "mci", "p", "q",
                    "effective_n"])
        for lab in POOLED_ORDER:
            edges, n_valid = pooled[lab]
            for (u, v, mci, p, q) in edges:
                w.writerow([lab, u, v, repr(mci), repr(p), repr(q),
                            n_valid])
    with open(DATA_DIR / "edges_within.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "source", "target", "mci", "q",
                    "effective_n"])
        for lab in WITHIN_ORDER_ALL8:
            edges, n_valid = within[lab][0], within[lab][1]
            if not edges:
                w.writerow([lab, "", "", "", "", n_valid])
            for (u, v, mci, _p, q) in edges:
                w.writerow([lab, u, v, repr(mci), repr(q), n_valid])
    with open(DATA_DIR / "contrasts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["stratum", "edge", "shrinking_mci", "shrinking_n",
                    "nonshrinking_mci", "nonshrinking_n", "z", "p",
                    "p_bonferroni"])
        for r in contrasts:  # ALL rows, incl. CS->JW and specificity
            w.writerow([r["stratum"], r["edge"],
                        repr(r["shrinking_mci"]), r["shrinking_n"],
                        repr(r["nonshrinking_mci"]),
                        r["nonshrinking_n"], repr(r["z"]),
                        repr(r["p"]), repr(r["p_bonferroni"])])
    if EXPORT_ALL_EDGES:
        pd.DataFrame(all_pooled_rows).to_csv(
            DATA_DIR / "all_edges_pooled.csv", index=False)
        pd.DataFrame(all_within_rows).to_csv(
            DATA_DIR / "all_edges_within.csv", index=False)

    # ---- tables/ ------------------------------------------------------
    pd.DataFrame([{"model": lab,
                   "n_reported_edges": len(pooled[lab][0])}
                  for lab in POOLED_ORDER]).to_csv(
        TABLE_DIR / "table3_edge_counts.csv", index=False)
    rows = []
    for lab in POOLED_ORDER:
        edges, n_valid = pooled[lab]
        for (u, v, mci, p, q) in edges:
            rows.append({"model": lab, "source_t_minus_1": u,
                         "target_t": v, "mci": f"{mci:.3f}",
                         "p": f"{p:.1e}", "q": f"{q:.1e}",
                         "effective_n": n_valid})
    pd.DataFrame(rows).to_csv(TABLE_DIR / "table4_pooled_edges.csv",
                              index=False)
    rows = []
    for lab in POOLED_ORDER:
        for dom, s in sorted(out_strength(pooled[lab][0]).items(),
                             key=lambda kv: -kv[1]):
            rows.append({"model": lab, "domain": dom,
                         "out_strength": f"{s:.3f}"})
    pd.DataFrame(rows).to_csv(TABLE_DIR / "table5_out_strength.csv",
                              index=False)
    rows = []
    for lab in WITHIN_ORDER_65:  # manuscript Table 6 is 65+ only
        edges, n_valid = within[lab][0], within[lab][1]
        if not edges:
            rows.append({"model": lab, "source_t_minus_1":
                         "(no edges satisfy the reporting criteria)",
                         "target_t": "", "mci": "", "q": "",
                         "effective_n": n_valid})
        for (u, v, mci, _p, q) in edges:
            rows.append({"model": lab, "source_t_minus_1": u,
                         "target_t": v, "mci": f"{mci:.3f}",
                         "q": f"{q:.3g}", "effective_n": n_valid})
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / "table6_within_edges.csv", index=False)
    pd.DataFrame([{"stratum": r["stratum"], "edge": r["edge"],
                   "shrinking_mci": f"{r['shrinking_mci']:.3f}",
                   "shrinking_n": r["shrinking_n"],
                   "nonshrinking_mci": f"{r['nonshrinking_mci']:.3f}",
                   "nonshrinking_n": r["nonshrinking_n"],
                   "z": f"{r['z']:.2f}", "p": f"{r['p']:.4f}",
                   "p_bonferroni": (f"{r['p_bonferroni']:.4f}"
                                    if not math.isnan(r["p_bonferroni"])
                                    else "")}
                  for r in contrasts]).to_csv(
        TABLE_DIR / "table7_contrasts.csv", index=False)

    # ---- figures ------------------------------------------------------
    if not SKIP_FIGURES:
        pooled_edges_only = {k: v[0] for k, v in pooled.items()}
        within_edges_only = {k: v[0] for k, v in within.items()}
        pos = build_positions(pooled_edges_only, within_edges_only)
        # Figure 2: all eight pooled strata (benchmark row above its
        # shrinking subset within each scale; 4 rows x 2 columns)
        network_figure(pooled_edges_only, FIG2_LABELS, FIG2_TAGS,
                       pos, "figure2")
        # Figure 3: person-demeaned within-person networks (65+)
        network_figure(within_edges_only, WITHIN_ORDER_65,
                       ["(a)", "(b)", "(c)", "(d)"], pos, "figure3",
                       empty_note_on_zero=True)
        # Figure 4: decline-contrast forest plot (WLB rows only, as in
        # the manuscript; the CS->JW contrast is in Table 7 and
        # data/contrasts.csv)
        forest_figure([r for r in contrasts
                       if "CS->JW" not in r["stratum"]])

    print(f"[done] data/, tables/, figures/ written "
          f"({time.time()-t0:.0f}s). All outputs are aggregate-level "
          "and safe to share; the workbook stays local.")


if __name__ == "__main__":
    main(sys.argv[1:])
