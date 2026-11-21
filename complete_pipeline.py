#!/usr/bin/env python3
"""Complete analysis pipeline (repository release v2.1.1): one command from the survey workbook to
every estimate, table and figure of the article.

  python complete_pipeline.py --xlsx /path/to/workbook.xlsx --outdir results [--workers 4]

The pipeline runs the scripts of this repository in a fixed order and with fixed settings.  The estimation
scripts are the files with which the article's results were produced and are not modified by this driver.

  stage      script(s)                               output directory     article elements
  ---------  --------------------------------------  -------------------  ---------------------------------
  network    validation_analysis.py                  OUT/network          Tables 3-6; Table 7, Block A;
             --window strict3 --design mci                                Figures 2-3; Tables S1, S2, S5, S6
  targeted   validation_analysis.py                  OUT/targeted         Table 7, Block B; Table 8;
             --window lag1 --design lag1only                              Tables S1-S4
  care       care_links.py                           OUT/care             Figure 4; Tables S4(c, d), S7(b), S8
  tables     make_tables.py                          OUT/tables           Tables 1-8 and S1-S8 (tables.xlsx, CSV)
  figures    make_figures.py                         OUT/figures          Figures 2-4

Network run: every lagged value that enters a test (source, history of the target, conditioning variables at
lag 1 and lag 2) belongs to the same respondent in consecutive survey years, because estimation is restricted
to three-consecutive-wave windows.  Targeted run: tests that involve no lag-2 term, on all valid lag-1
transitions.  Each stage resumes from its checkpoints, so an interrupted run can be restarted with the same
command; --stages selects stages (comma-separated) when earlier ones are already complete, and --network-dir /
--targeted-dir point to existing output directories of the two validation_analysis.py runs.

The permutation steps dominate the runtime (several hours sequentially).  --workers N computes their
replicates in N parallel processes per city scale (disjoint --rep-range slices, then one aggregating call);
results are identical to a sequential run because every replicate has its own random stream.

Only aggregate outputs are written.  The individual-level data are not part of the repository.  To check that
the code runs without the data, write a synthetic workbook (random content in the real layout) and run it:

  python complete_pipeline.py --make-synthetic synthetic.xlsx
  python complete_pipeline.py --xlsx synthetic.xlsx --outdir results_synthetic --quick --force
"""
import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

__version__ = "2.1.1"
HERE = Path(__file__).resolve().parent
STAGES = ["network", "targeted", "care", "tables", "figures"]
SCRIPTS = ["validation_analysis.py", "care_links.py", "make_tables.py", "make_figures.py", "complete_pipeline.py"]
RUNS = {  # settings of the two validation_analysis.py runs, as used for the article
    "network": dict(window="strict3", design="mci", steps="gate,networks,diag,common,boot", label_blocks="global", dose=False),
    "targeted": dict(window="lag1", design="lag1only", steps="gate,networks,diag,common,boot,ahiv,window_check,period_check",
                     label_blocks="global,region", dose=True),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def show(cmd):
    return " ".join(str(c) for c in cmd[1:])


def run(cmd, dry):
    print(f"\n[pipeline] {show(cmd)}", flush=True)
    if not dry:
        subprocess.run([str(c) for c in cmd], check=True, cwd=HERE)


def run_parallel(cmds, dry, stagger):
    """Start the commands as parallel processes (staggered, so that their bookkeeping writes do not collide)."""
    procs = []
    for cmd in cmds:
        print(f"\n[pipeline] (worker) {show(cmd)}", flush=True)
        if not dry:
            procs.append(subprocess.Popen([str(c) for c in cmd], cwd=HERE))
            time.sleep(stagger)
    failed = [p.args for p in procs if p.wait() != 0]
    if failed:
        raise SystemExit(f"[pipeline] {len(failed)} worker process(es) failed; rerun the same command to resume")


def slices(n, k):
    k = max(1, min(k, n))
    edges = [round(i * n / k) for i in range(k + 1)]
    return [f"{a}:{b}" for a, b in zip(edges[:-1], edges[1:]) if b > a]


def validation_run(name, a, data, outdir):
    cfg = RUNS[name]
    base = [sys.executable, HERE / "validation_analysis.py", *data, "--outdir", outdir, "--window", cfg["window"],
            "--design", cfg["design"], "--seed", a.seed, "--nboot", a.nboot, "--nrep-time", a.nrep_time,
            "--nrep-label", a.nrep_label, "--nperm-dose", a.nperm_dose] + (["--force"] if a.force else [])
    run(base + ["--steps", cfg["steps"], "--label-blocks", cfg["label_blocks"]], a.dry_run)
    for step, nrep in (("perm_time", a.nrep_time), ("perm_label", a.nrep_label)):
        tail = ["--steps", step, "--label-blocks", cfg["label_blocks"]]
        if a.workers > 1:
            for scale in ("SMC", "Large"):
                run_parallel([base + tail + ["--scales", scale, "--rep-range", s] for s in slices(nrep, a.workers)], a.dry_run, a.stagger)
        run(base + tail, a.dry_run)                     # computes missing replicates (if any) and aggregates
    if cfg["dose"]:
        run(base + ["--steps", "dose", "--label-blocks", cfg["label_blocks"]], a.dry_run)


# ---------------------------------------------------------------------------------------------------------
# synthetic workbook for a dry run of the pipeline (code check only; NOT data, NOT results)
# ---------------------------------------------------------------------------------------------------------
ANALYSIS_HEADERS = [
    "PID", "LifeSatisfaction", "Finances & assets", "Jobs & wages", "Housing", "Work–life balance", "Health", "Education",
    "Friends & community", "Trust in government & courts", "Air & water (environment)", "Personal safety",
    "Parenting support", "Caregiving support", "Enjoyment of life", "CityName", "CityCode-6", "CityCode-5", "CityType",
    "Areakm2", "Maeged_flag", "Metropolice_flag", "Japanese", "Foreign", "Population", "Inflow", "Outflow",
    "Sex", "Age", "Year", "Education.1", "Income", "household",
]
POP_HEADERS = ["調査年 コード", "調査年", "Year", "地域 コード", "CityCode-5", "地域", "Japanese", "Foreign", "Population", "Inflow", "Outflow"]


def synthetic_frames(n_respondents=9000, n_large=60, n_smc=320, seed=1):
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(seed)
    munis = []
    for k in range(n_large + n_smc):
        pref = int(rng.integers(1, 48))
        munis.append({"code": f"{pref:02d}{100 + k:03d}", "type": "1_large" if k < n_large else str(rng.choice(["2_mid", "3_small"])),
                      "rate": float(rng.normal(-0.006, 0.008)), "p0": int(rng.integers(8000, 900000))})
    large = [m for m in munis if m["type"] == "1_large"]
    smc = [m for m in munis if m["type"] != "1_large"]
    rows = []
    for pid in range(1, n_respondents + 1):
        pool = large if rng.random() < 0.5 else smc
        m = pool[int(rng.integers(len(pool)))]
        waves = int(rng.choice([1, 2, 3, 4, 5, 6], p=[0.12, 0.13, 0.25, 0.22, 0.18, 0.10]))
        start = 2019 + int(rng.integers(0, 8 - waves))
        age = int(rng.integers(11, 15)) if rng.random() < 0.38 else int(rng.integers(3, 11))
        sex, edu, inc, hh = int(rng.integers(1, 3)), int(rng.integers(1, 7)), int(rng.integers(1, 11)), int(rng.integers(1, 7))
        trait = rng.normal(6.0, 1.4) + rng.normal(0, 0.9, 14)
        prev = rng.normal(0, 1, 14)
        for t in range(waves):
            cur = 0.25 * prev + 0.35 * rng.normal() + rng.normal(0, 1, 14)
            cur[1:] += 0.30 * prev[0]                    # arbitrary lagged links, so that some edges are reported
            cur[4] -= 0.35 * prev[12]
            codes = np.clip(np.rint(trait + cur), 0, 10).astype(int) + 1          # workbook codes 1..11
            rows.append([pid, *codes.tolist(), "synthetic", m["code"] + "0", m["code"], m["type"], 100.0, 0, 0,
                         m["p0"], 0, m["p0"], 0, 0, sex, age, start + t, edu, inc, hh])
            prev = cur
    data = pd.DataFrame(rows, columns=ANALYSIS_HEADERS)
    pop = []
    for m in munis:
        for year, p in ((2017, m["p0"]), (2023, int(round(m["p0"] * (1.0 + m["rate"]) ** 6)))):
            pop.append([year, f"{year}", year, m["code"], m["code"], "synthetic", p, 0, p, 0, 0])
    return data, pd.DataFrame(pop, columns=POP_HEADERS)


def synthetic_workbook(path, respondents=9000, seed=1):
    """Random respondents, municipalities and scores in the sheet and column layout of data/schema.csv."""
    import pandas as pd
    data, pop = synthetic_frames(respondents, seed=seed)
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        data.to_excel(xw, sheet_name="24_分析データ", index=False)
        pop.to_excel(xw, sheet_name="17_Population", index=False)
    print(f"[synthetic] {len(data):,} person-years, {data.PID.nunique():,} respondents, {pop['CityCode-5'].nunique()} municipalities -> {path}")
    print("[synthetic] code check only: run complete_pipeline.py with --quick --force; nothing from this file may be reported")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", help="survey workbook (sheets 24_分析データ and 17_Population; see data/schema.csv)")
    ap.add_argument("--data-csv"); ap.add_argument("--pop-csv")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--network-dir", default=None, help="existing or new output directory of the network run (default OUT/network)")
    ap.add_argument("--targeted-dir", default=None, help="existing or new output directory of the targeted run (default OUT/targeted)")
    ap.add_argument("--stages", default=",".join(STAGES), help="comma-separated subset of: " + ", ".join(STAGES))
    ap.add_argument("--workers", type=int, default=1, help="parallel processes for the permutation replicates")
    ap.add_argument("--stagger", type=float, default=20.0, help="seconds between the starts of parallel workers")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--nboot", type=int, default=2000, help="municipality block-bootstrap replicates")
    ap.add_argument("--nrep-time", type=int, default=500, help="within-person time-permutation replicates")
    ap.add_argument("--nrep-label", type=int, default=500, help="municipality label-permutation replicates")
    ap.add_argument("--nperm-dose", type=int, default=2000, help="randomisation replicates of the rate interaction")
    ap.add_argument("--nperm-family", type=int, default=2000, help="randomisation replicates of the 13-target screen")
    ap.add_argument("--make-synthetic", metavar="PATH", help="write a synthetic workbook (code check only) to PATH and exit")
    ap.add_argument("--quick", action="store_true", help="tiny replicate counts: checks that the code runs, NOT for results")
    ap.add_argument("--force", action="store_true", help="continue past a failed Gate 1 (Table 1 reproduction)")
    ap.add_argument("--dry-run", action="store_true", help="print the commands without running them")
    a = ap.parse_args(argv)
    if a.make_synthetic:
        return synthetic_workbook(a.make_synthetic)
    if a.quick:
        a.nboot, a.nrep_time, a.nrep_label, a.nperm_dose, a.nperm_family = 40, 4, 4, 400, 40
    if a.xlsx:
        data = ["--xlsx", Path(a.xlsx).resolve()]
    elif a.data_csv and a.pop_csv:
        data = ["--data-csv", Path(a.data_csv).resolve(), "--pop-csv", Path(a.pop_csv).resolve()]
    else:
        ap.error("supply --xlsx, or --data-csv and --pop-csv")
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    if set(stages) - set(STAGES):
        ap.error("unknown stage; choose from " + ", ".join(STAGES))
    out = Path(a.outdir).resolve()
    net = Path(a.network_dir).resolve() if a.network_dir else out / "network"
    strain = Path(a.targeted_dir).resolve() if a.targeted_dir else out / "targeted"
    care, tabs, figs = (out / d for d in ("care", "tables", "figures"))
    t0 = time.time()
    if not a.dry_run:
        out.mkdir(parents=True, exist_ok=True)
        log = out / "pipeline_run.json"
        meta = json.loads(log.read_text(encoding="utf-8")) if log.exists() else {"calls": []}
        meta["calls"].append({
            "utc": datetime.now(timezone.utc).isoformat(), "pipeline_version": __version__, "stages": stages, "workers": a.workers,
            "seed": a.seed, "replicates": {"nboot": a.nboot, "nrep_time": a.nrep_time, "nrep_label": a.nrep_label,
                                           "nperm_dose": a.nperm_dose, "nperm_family": a.nperm_family},
            "quick": a.quick, "force": a.force, "python": sys.version.split()[0], "platform": platform.platform(),
            "packages": {n: importlib.metadata.version(n) for n in ("numpy", "pandas", "scipy", "matplotlib", "openpyxl", "tigramite")},
            "script_sha256": {s: sha256(HERE / s) for s in SCRIPTS}})
        log.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    for stage in [s for s in STAGES if s in stages]:
        print(f"\n[pipeline] ===== stage: {stage} =====", flush=True)
        if stage in RUNS:
            validation_run(stage, a, data, net if stage == "network" else strain)
        elif stage == "care":
            run([sys.executable, HERE / "care_links.py", *data, "--va-dir", HERE, "--strain-dir", strain, "--outdir", care,
                 "--nperm", a.nperm_family, "--seed", a.seed], a.dry_run)
        elif stage == "tables":
            run([sys.executable, HERE / "make_tables.py", *data, "--network-dir", net, "--strain-dir", strain, "--care-dir", care,
                 "--outdir", tabs], a.dry_run)
        elif stage == "figures":
            run([sys.executable, HERE / "make_figures.py", "--network-dir", net, "--strain-dir", strain, "--care-dir", care,
                 "--outdir", figs], a.dry_run)
    if a.quick:
        print("\n[pipeline] --quick was used: these outputs check the code only and must not be reported", flush=True)
    print(f"\n[pipeline] done in {(time.time() - t0) / 60:.1f} min: {out}", flush=True)


if __name__ == "__main__":
    main()
