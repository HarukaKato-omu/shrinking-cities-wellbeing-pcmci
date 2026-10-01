# Population decline and the temporal structure of well-being

Replication code for:

> Kato, H. *Population decline and the temporal structure of well-being:
> Contextual strain amid between-person stability in Japan's shrinking
> cities.* *Cities*, in revision.

This repository contains the Python code used to estimate the temporal
(lag-1) dependency structure among fourteen satisfaction domains in
Japanese municipalities, stratified by city scale (Large cities vs.
Small-to-Medium Cities, SMCs) and age cohort (15–64 vs. 65+), using PCMCI
(`tigramite`) with a within-person transition mask (`mask_type="y"`),
Benjamini–Hochberg multiplicity control with an effect-size floor, a
person-demeaned within-person specification, and decline contrasts.
It has two layers:

1. **`complete_pipeline.py` (v1.0.0, unchanged since release 1.0.0)** goes
   from the survey workbook to every estimate in the article: panel
   construction, six self-verification gates, the 28 PCMCI models, and the
   deterministic generation of **Tables 1–7** and **Figures 2–3**.
2. **`validation_analysis.py` (added in release 1.1.0)** runs on top of the
   pipeline's specification without altering any of its estimates and
   produces the validation results of manuscript Section 3.5: the
   additional columns of Tables 6 and 7, **Table 8**, the inputs of
   **Figure 4** (drawn by `make_figure4.py`), and **Supplementary Tables
   S1–S4**.

(Figure 1, the municipal map, is produced in GIS and is not part of the
code.)

## Repository structure

```
├── complete_pipeline.py         # Full analysis and reporting pipeline (v1.0.0)
├── validation_analysis.py       # Validation layer: Tables 6-8 (added columns), Figure 4 inputs, S1-S4
├── make_figure4.py              # Figure 4 (dose-response) from validation outputs
├── notebooks/
│   └── run_validation_local.ipynb   # Optional: parallel workers on a local machine
├── data/
│   └── schema.csv               # Workbook layout: variable definitions and sources (no raw data)
├── .gitignore                   # Prevents accidental data/output commits
├── requirements.txt
├── CITATION.cff
├── LICENSE
└── README.md
```

Running the code creates output directories (`--outdir`) that are not
tracked in the repository.

## Data availability

The individual-level workbook is **not** redistributed here and cannot
be reconstructed from public sources. The satisfaction items and
respondent covariates come from the Survey on Satisfaction and Quality
of Life (Cabinet Office of Japan, waves 2019–2025), whose individual
data the Cabinet Office provides on application for approved academic
use, subject to conditions that prohibit re-identification and
re-provision to third parties
(<https://www5.cao.go.jp/keizai2/wellbeing/manzoku/index.html>).
Municipal population counts (2017, 2023) come from the Basic Resident
Register (e-Stat of Japan) and are public. `data/schema.csv` documents
both sheets that the code reads — `24_分析データ` (33 columns, read
positionally) and `17_Population` (11 columns) — with each variable's
type, unit, definition, source, access conditions, and whether it is
used, so that researchers who obtain the individual data can assemble an
equivalent workbook and run both scripts without modification.

All outputs are **aggregate-level** (edge lists, centralities,
contrasts, counts, permutation summaries) and safe to share; nothing
person-level is written to disk.

## Requirements

- Python ≥ 3.10
- `numpy`, `pandas`, `scipy`, `openpyxl`, `matplotlib`, `tigramite` ≥ 5.2
  (the article's results were produced with tigramite 5.2.10.1)

```
pip install -r requirements.txt
```

## How to run

### 1. Main pipeline (Tables 1–7, Figures 2–3)

```
python complete_pipeline.py --xlsx /path/to/workbook.xlsx --outdir results
```

Reads the workbook, builds the panel, runs Gate 1 and stops on failure,
then estimates all 28 PCMCI models and writes the outputs. Options:
`--force` (proceed past a failed Gate 1), `--jitter-seed N`
(jitter-invariance check), `--skip-figures`. **Runtime** roughly 30–90
minutes.

The pipeline also writes `figures/figure4.*`, the decline-contrast
forest plot of the first revision. That figure is no longer in the
article; the article's Figure 4 is produced by `make_figure4.py` (below).
The pipeline is kept unchanged so that every number in Tables 1–7 can be
regenerated exactly.

### 2. Validation layer (Tables 6–8, Figure 4, Supplementary S1–S4)

`validation_analysis.py` is standalone (it does not import the pipeline)
but reproduces the pipeline's specification exactly: identical panel
construction and jitter, the same mask, PCMCI parameters, and reporting
criteria. Its first step (Gate 1) checks the 30 sample-structure values of
Table 1 and stops on any mismatch; its target fits reproduce the Table 7
coefficients, and every fixed-conditioning-set partial correlation is
checked against the corresponding Tigramite MCI value.

```
python validation_analysis.py --xlsx /path/to/workbook.xlsx --outdir validation \
    --steps gate,diag,common,boot,ahiv,perm_time,perm_label,window_check,period_check,dose \
    --nboot 2000 --nrep-time 500 --nrep-label 500 --nperm-dose 2000 --seed 2026
python make_figure4.py --indir validation --outdir figures
```

| Step | What it does | Article element |
|---|---|---|
| `gate` | Table 1 reproduction (hard stop unless `--force`) | Table 1 |
| `diag` | wave-count distributions, contemporaneous correlation, implied short-panel bias, own-lag MCI | Table S1 |
| `common` | MCI conditioning sets; re-estimation under the union set, own lag only, and no conditioning | Table 7 note; Table S2 |
| `boot` | municipality block bootstrap of the contrast (conditioning sets fixed) | Table 7, cluster-bootstrap z |
| `perm_time` | within-person time permutation: per-edge null bands for all 182 lagged edges, contrast statistics | Table 6, null band and permutation p; Table S1 |
| `perm_label` | municipality label permutation: single-edge, E3-family (13 contrasts), all-182 and paper-rule p-values; `--label-blocks global,region,prefecture` | Table 7, permutation p |
| `ahiv` | Anderson–Hsiao IV in first differences, binary subsets | Table S4(a) |
| `window_check`, `period_check` | transition-window and period/position sensitivity | Table S3 |
| `dose` | continuous moderation by the municipal population change rate: within OLS with own lag and year fixed effects, cluster SE, randomisation inference, AH-IV with interaction, quintile profile, family test over all 182 edges; then the Figure 4 inputs and the region robustness check | Table 8; Figure 4; Table S4(b, c) |
| `dose_extra` | only the Figure 4 inputs and the region robustness check (a few minutes) | Figure 4; Table 8 note |

Repeating the same command resumes cached fits and permutation
replicates; separate calls per `--scales` merge their tables. Permutation
p-values are (b + 1)/(B + 1). **Runtime** with the counts above: the PCMCI
permutation steps (`perm_time`, `perm_label`) dominate — several hours
sequentially, about one hour with parallel workers (`--rep-range START:END`
slices per step/scale, aggregated by a final call without `--rep-range`;
see `notebooks/run_validation_local.ipynb`). All other steps take minutes.

Options: `--window strict3` (three consecutive waves at all stages, used
only by `window_check`), `--exclude-movers`, `--scales SMC` or `Large`,
`--force`, `--rep-range`.

## Output ↔ article mapping

### `complete_pipeline.py --outdir results`

| Output | Article element |
|---|---|
| `tables/table1_sample_structure.csv` | Table 1 (with the Gate 1 pipeline-vs-manuscript comparison) |
| `tables/table2_descriptives.csv` | Table 2 |
| `tables/table3_edge_counts.csv` | Table 3 |
| `tables/table4_pooled_edges.csv` | Table 4 |
| `tables/table5_out_strength.csv` | Table 5 |
| `tables/table6_within_edges.csv` | Table 6, columns MCI, q-value, effective n |
| `tables/table7_contrasts.csv` | Table 7, columns MCI (n) and Fisher z (p) |
| `figures/figure2.(png\|pdf)` | Figure 2 — pooled networks, all eight strata |
| `figures/figure3.(png\|pdf)` | Figure 3 — person-demeaned within-person networks (65+) |
| `figures/figure4.(png\|pdf)` | decline-contrast forest plot (first revision; not in the article) |
| `data/edges_pooled.csv`, `data/edges_within.csv` | reported edges, full precision |
| `data/all_edges_pooled.csv`, `data/all_edges_within.csv` | all 182 lagged hypotheses per model, incl. non-significant |
| `data/contrasts.csv` | decline contrasts, full precision (targeted family + age-pattern check) |
| `data/sensitivity_edge_counts.csv` | edge counts under four reporting criteria (Section 4.3) |
| `data/composition_check.csv` | pooled models on the 3-wave-plus subsample (Section 4.3) |
| `data/contemporaneous_pairs.csv` | lag-0 verification layer (Section 4.2) |
| `data/run_metadata.json` | all parameters + software versions |

### `validation_analysis.py --outdir validation`

| Output | Article element |
|---|---|
| `G1_table1_check.csv` | Gate 1 (30 Table 1 values) |
| `V3_table6_edges_vs_null_{SMC,Large}.csv` | Table 6, null 95% band and permutation p (all 182 edges per model) |
| `V3_perm_time_{SMC,Large}.json` / `.csv` | Table 6 note (edge counts under the null); Table S1 |
| `V4_cluster_boot.csv` | Table 7, cluster-bootstrap z (p) |
| `V6_perm_label_{SMC,Large}_{global,region,prefecture}.json` / `.csv` | Table 7, permutation p (single edge, E3 family, all 182); region-blocked values in the note |
| `V9_dose_response.csv` | Table 8; Table S4(b) (AH-IV rows) |
| `V9_dose_response_family.csv` | Table S4(c) |
| `V9_dose_response_region.csv` | Table 8 note (region fixed effects; region-blocked randomisation) |
| `V9_figure4_lines.json`, `V9_figure4_quintile_slopes.csv` | inputs of Figure 4 |
| `V1_diagnostics.csv`, `V1_wave_distribution.csv` | Table S1(a) |
| `V1_wave_distribution_all_strata.csv` | Table S1(b) |
| `V1_conditioning_sets.json`, `V2_common_set.csv` | Table 7 note; Table S2 |
| `V7_window_sensitivity.csv`, `V8_period_position.csv`, `V8_period_position_contrasts.csv` | Table S3 |
| `V5_anderson_hsiao.csv` | Table S4(a) |
| `run_metadata.json`, `sample_flow.json` | identity of the run (code and data hashes, seed, settings) and package versions |

### `make_figure4.py --indir validation --outdir figures`

| Output | Article element |
|---|---|
| `figures/figure4.(png\|pdf)` | Figure 4 — dose-response of the care-work dependency by population change rate |

## Reproducibility notes

- The transition mask restricts every lagged test to consecutive-year
  observations of the same respondent; applied to the dependent
  variable (`mask_type="y"`), it governs both MCI estimates and their
  analytic significance (manuscript Section 3.3). Conditioning variables
  at lag 2 (the parents of the source variable shifted by one lag in the
  MCI stage) are read from the preceding row, which for transitions
  without a three-year history belongs to a non-consecutive year or
  another respondent; `validation_analysis.py --steps common` quantifies
  the consequence (Supplementary Table S2, row "own lag only") and
  `--window strict3` re-estimates on three-consecutive-wave windows.
- Reporting criteria: BH q < .05 over the 182 directed lag-1
  hypotheses of each model and |MCI| ≥ .10; the sensitivity criterion
  lowers the floor to .05. `PC_ALPHA = 0.05` is the PCMCI
  condition-selection level.
- A negligible jitter (SD = 1e-5) breaks degenerate ties in covariate
  codes; the seed and draw order are fixed for exact stream
  reproduction, and no reported estimate depends on the seed
  (`--jitter-seed`).
- Person-demeaning is applied within each estimation subset for
  respondents observed in ≥ 3 waves (manuscript Section 3.4).
- The set of `city_type` values classified as Large is identified
  automatically and must reproduce the manuscript count of 122
  large-city municipalities exactly; otherwise the pipeline stops.
- `validation_analysis.py` seeds every permutation replicate from
  (seed, step, scale, replicate), so results are identical whether
  replicates are computed in one process or in parallel slices.
  `run_metadata.json` records the code and data hashes, the seed and
  settings of each run, and the Python and package versions.

## Releases

- **v1.0.0** (2026-09-11): `complete_pipeline.py`; Tables 1–7, Figures 2–4
  of the first revision.
- **v1.1.0** (2026-10-01): adds `validation_analysis.py` and
  `make_figure4.py` (manuscript Section 3.5, Tables 6–8, Figure 4,
  Supplementary Tables S1–S4). `complete_pipeline.py` is unchanged.

## Citation

If you use this code, please cite the paper above and, for the exact
released version of this code, the archived DOI (Zenodo; citation
metadata for this repository is provided in `CITATION.cff`).

## License

MIT — see `LICENSE`. The license applies to the code only and does not
extend to the survey data, which remain subject to the Cabinet Office's
conditions of provision.

## Contact

Haruka Kato, Osaka Metropolitan University
