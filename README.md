# Population decline and the temporal structure of well-being

Replication code for:

> Kato, H. *Population decline and the temporal structure of well-being:
> A candidate contextual strain amid between-person stability in Japan's
> shrinking cities.* *Cities*, in revision.

This repository contains the Python code that estimates the lagged
dependency structure among fourteen satisfaction domains in Japanese
municipalities, stratified by city scale (large cities vs. small and
medium-sized cities, SMCs) and age group (15–64 vs. 65+), with PCMCI
(`tigramite`), Benjamini–Hochberg multiplicity control with an
effect-size floor, a person-demeaned within-person specification, and
the decline contrasts and validation analyses of the article.

**Release v2.1.1** generates every result of the article: Tables 1–8,
Figures 2–4, and Supplementary Tables S1–S8. No estimate from an
earlier release is retained. (Figure 1, the municipal map, is produced
in GIS and is not part of the code.)

## One command

```
pip install -r requirements.txt
python complete_pipeline.py --xlsx /path/to/workbook.xlsx --outdir results --workers 4
```

`complete_pipeline.py` runs the scripts below in a fixed order and with
the settings used for the article, and writes everything under
`results/`. It resumes from checkpoints, so an interrupted run is
restarted with the same command.

| Stage | Script(s) | Output | Article elements |
|---|---|---|---|
| `network` | `validation_analysis.py --window strict3 --design mci` | `results/network/` | Tables 3–6; Table 7, Block A; Figures 2–3; Tables S1, S2, S5, S6 |
| `targeted` | `validation_analysis.py --window lag1 --design lag1only` | `results/targeted/` | Table 7, Block B; Table 8; Tables S1–S4 |
| `care` | `care_links.py` | `results/care/` | Figure 4; Tables S4(c, d), S7(b), S8 |
| `tables` | `make_tables.py` | `results/tables/` | Tables 1–8 and S1–S8 (`tables.xlsx` and CSV files) |
| `figures` | `make_figures.py` | `results/figures/` | Figures 2–4 |

Options: `--stages network,targeted,…` (run a subset), `--workers N`
(parallel processes for the permutation replicates), `--network-dir` /
`--targeted-dir` (existing output directories of the two runs),
`--dry-run` (print the commands), `--force` (continue past a failed
Gate 1), `--quick` (tiny replicate counts, code check only),
`--make-synthetic` (see below), and the replicate counts
`--nboot`, `--nrep-time`, `--nrep-label`, `--nperm-dose`,
`--nperm-family` (defaults 2000, 500, 500, 2000, 2000, as in the
article). Each script can also be run on its own; its docstring gives
the usage.

**Runtime.** The permutation steps of the two runs dominate: several
hours sequentially, and proportionally less with `--workers N`. All
other steps take minutes.

## The two estimation samples

The lag examined is one year in every analysis, and the estimation
sample follows what each model requires (manuscript Section 3.3).

- **Network run (three-consecutive-wave windows).** The MCI test of a
  link from X at t−1 to Y at t conditions on the lagged parents of Y and
  on the lagged parents of X, the latter at t−2. An observation at t is
  therefore retained as a target only if the same respondent was also
  observed in the two preceding survey years. Every lagged value that
  enters a test — as the source, as the history of the target, or as a
  conditioning variable at lag 1 or lag 2 — then belongs to the same
  respondent's consecutive survey history, in the PC stage and in the
  MCI stage alike. All pooled and person-demeaned network models
  (Tables 3–6, Figures 2–3) and Block A of Table 7 are estimated on
  these windows.
- **Targeted run (all valid lag-1 transitions).** Block B of Table 7 is
  a lagged partial correlation conditioned on the lagged parents of the
  target only (`max_conds_px=0`; parent selection with a maximum lag of
  one), and Table 8 / Figure 4 are within-person regressions with the
  target's own lag. Neither contains a lag-2 term, so every value they
  use belongs to a valid consecutive pair. Neither is an MCI statistic.
  Both are re-estimated on the windows (Supplementary Tables S3 and S4).

## Repository structure

```
├── complete_pipeline.py     # one command: runs every stage below
├── validation_analysis.py   # estimation: PCMCI networks, contrasts, bootstrap, permutation tests,
│                            #   Anderson–Hsiao, rate interaction (network run and targeted run)
├── care_links.py            # Caregiving Support links: the interaction model for each of the 13 targets,
│                            #   coefficients within rate quintiles, subset regressions
├── make_tables.py           # Tables 1–8 and Supplementary Tables S1–S8
├── make_figures.py          # Figures 2–4
├── data/
│   └── schema.csv           # workbook layout: variable definitions and sources (no raw data)
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
equivalent workbook and run the pipeline without modification.

All outputs are **aggregate-level** (edge lists, coefficients,
contrasts, counts, permutation summaries) and safe to share; nothing
person-level is written to disk.

### Running the code without the data

```
python complete_pipeline.py --make-synthetic synthetic.xlsx
python complete_pipeline.py --xlsx synthetic.xlsx --outdir results_synthetic --quick --force
```

The synthetic workbook has the layout of the real one and randomly
generated content. It checks that every stage runs; its outputs have no
relation to the survey and must not be reported. `--force` is needed
because Gate 1 compares the sample structure with Table 1.

## Requirements

- Python ≥ 3.10
- `numpy`, `pandas`, `scipy`, `openpyxl`, `matplotlib`, `tigramite` ≥ 5.2

The article's results were produced with Python 3.13.9, numpy 2.3.5,
pandas 2.2.3, scipy 1.17.0 and tigramite 5.2.10.1.

## Output ↔ article mapping

### Tables (`results/tables/`, also as sheets of `tables.xlsx`)

| Output | Article element |
|---|---|
| `T1_sample_structure.csv` | Table 1 |
| `T2_descriptives.csv`, `T2_descriptives_wide.csv` | Table 2 |
| `T3_edge_counts.csv` | Table 3 (main criteria and the two sensitivity criteria of the note) |
| `T4_pooled_edges.csv` | Table 4 |
| `S_depth_life_satisfaction_hub.csv`, `T5_out_strength.csv` | Table 5 (conditioning depths; reported outgoing edges and out-strength) |
| `T6_within_edges.csv`, `T6b_null_summary_and_care_links.csv` | Table 6 and its note (null bands, permutation p, edges under the null) |
| `T7_decline_contrast.csv` | Table 7, Blocks A and B (65+): coefficients, Fisher z, cluster bootstrap, permutation p |
| `T7d_other_contrasts.csv` | Table 7, age-pattern checks (15–64); Supplementary Table S7(a) |
| `T8_dose_response.csv`, `T8c_region_robustness.csv` | Table 8 and its note |
| `T7c_conditioning_conventions.csv` | Supplementary Table S2 |
| `T7e_position_and_period.csv`, `T7f_window_by_design.csv` | Supplementary Table S3 |
| `T7b_anderson_hsiao.csv`, `T8_dose_response.csv` (Anderson–Hsiao rows) | Supplementary Table S4(a, b) |
| `S4c_interaction_placement.csv`, `S4d_rate_region_windows.csv` | Supplementary Table S4(c, d) |
| `S_depth_negative_pooled_edges.csv` | Supplementary Table S5 |
| `T5_out_strength.csv` | Supplementary Table S6 |
| `S7b_jobs_wages_regressions.csv` | Supplementary Table S7(b) |
| `S8a_family_13_targets.csv`, `S8b_quintile_coefficients.csv` | Supplementary Table S8 |

Supplementary Table S1 is read from `V1_diagnostics.csv`,
`V1_wave_distribution*.csv`, `V1_conditioning_sets.json` and
`V3_perm_time_*.json` of the two runs. The remaining sheets of
`tables.xlsx` (`S_two_layers_summary`, `T5b_…`, `S_lag1only_…`) compare
the network run with the lag-1-only network layer of the targeted run;
they are not reported in the article.

### Estimation outputs

| Output | Content |
|---|---|
| `network/`, `targeted/`: `N_all_edges.csv`, `N_network_summary.csv` | all 182 directed lagged links of every model (coefficient, p, q, reported flag, conditioning set) |
| `G1_table1_check.csv` | Gate 1: 30 values of Table 1 |
| `V1_*`, `V2_common_set.csv` | short-panel diagnostics, wave-count distributions, conditioning sets and conventions |
| `V3_perm_time_*`, `V3_table6_edges_vs_null_*` | within-person time permutation (500 replicates per model) |
| `V4_cluster_boot.csv` | municipality block bootstrap of the contrast (2,000 replicates) |
| `V5_anderson_hsiao.csv` | Anderson–Hsiao estimates in first differences |
| `V6_perm_label_*` | municipality label permutation (500 replicates) |
| `V7_window_sensitivity.csv`, `V8_period_position*.csv` | window, run-position and survey-period sensitivity |
| `V9_dose_response*.csv`, `V9_figure4_*` | continuous moderation by the municipal annual population change rate |
| `care/V10_family_<scale>_<sample>.csv` | the interaction model for each of the 13 targets (lag-1 transitions and windows) |
| `care/V11_quintile_slopes.csv`, `care/V11_subset_regressions.csv` | coefficients within rate quintiles; shrinking and non-shrinking subset regressions |
| `run_metadata.json`, `sample_flow.json`, `pipeline_run.json` | identity of each run: code and data hashes, seed, settings, package versions |

### Figures (`results/figures/`)

| Output | Article element |
|---|---|
| `figure2.(png\|pdf)` | Figure 2 — pooled networks, eight strata |
| `figure3.(png\|pdf)` | Figure 3 — person-demeaned networks in SMCs, by age group |
| `figure4.(png\|pdf)` | Figure 4 — Caregiving Support → 13 target domains by the pace of population decline |

## Reproducibility notes

- **Mask.** The window mask (network run) or the transition mask
  (targeted run) is applied to the dependent variable
  (`mask_type="y"`), so the estimates and their analytic significance
  are computed only from valid targets. The stacked person-year data
  are ordered by respondent and survey year.
- **Reporting criteria.** BH q < .05 over the 182 directed lag-1
  hypotheses of each model and |MCI| ≥ .10. `PC_ALPHA = 0.05` is the
  PCMCI condition-selection level; the maximum lag is one.
- **Person-demeaning** is applied within each estimation subset for
  respondents observed in at least three waves.
- **Time permutation.** For each respondent, the complete vector of 14
  satisfaction scores observed in a wave is moved as one unit and the
  vectors are re-ordered across that respondent's waves
  (`permute_within_person`); the 14 variables are never permuted
  independently. Covariates and survey years keep their slots, the data
  are demeaned again, and the full PCMCI procedure is re-run.
- **Random streams.** Every bootstrap and permutation replicate is
  seeded from (seed, step, scale, replicate), so results are identical
  whether replicates are computed in one process or in parallel slices.
  Permutation p-values are (b + 1)/(B + 1).
- **Jitter.** A negligible jitter (SD = 1e-5, fixed seed) breaks
  degenerate ties in covariate codes; no reported estimate depends on
  it.
- **Checks.** Gate 1 stops the run unless the workbook reproduces the
  sample structure of Table 1. Every fixed-conditioning-set partial
  correlation is checked against the corresponding Tigramite value.
  `care_links.py` and `make_tables.py` print a comparison with the
  published values of Table 8, Supplementary Tables S7(b) and S8(b),
  and Table 1.
- **Provenance of the scripts.** `validation_analysis.py` is the file
  with which the article's results were produced; `run_metadata.json`
  records its SHA-256 for each run. `care_links.py` and
  `make_tables.py` consolidate the scripts used for the article (the
  13-target screen and the table assembly) with their estimation code
  and random streams unchanged; on identical inputs they reproduce the
  outputs of those scripts byte for byte. They also compute the
  quintile coefficients and subset regressions of Supplementary
  Tables S7(b) and S8(b) and Tables 1–2. `make_figures.py` redraws
  Figures 2–4 from the aggregate outputs; typographic details of the
  published figures (fonts, label placement) may differ slightly.
  `complete_pipeline.py` adds no estimate: it runs the other scripts
  in order.
- **Script version strings.** Version constants inside individual
  scripts record each file's own history; the repository release is
  v2.1.1.

## Releases

- **v1.0.0** (2026-09-11): `complete_pipeline.py`; Tables 1–7 and
  Figures 2–4 of the first revision.
- **v1.1.0** (2026-10-01): adds `validation_analysis.py` and
  `make_figure4.py` (second revision).
- **v2.1.1** (2026-11-21): third revision. The network models are
  estimated on three-consecutive-wave windows, so that every lagged
  value entering a test belongs to the same respondent's consecutive
  survey history; all tables and figures are regenerated under this
  specification. `complete_pipeline.py` now runs every stage; the
  first-revision pipeline and `make_figure4.py` are superseded and
  remain available under the tags v1.0.0 and v1.1.0.

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
