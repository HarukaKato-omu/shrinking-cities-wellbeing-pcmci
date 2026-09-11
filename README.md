# Population decline and the temporal structure of well-being

Replication code for:

> Kato, H. *Population decline and the temporal structure of well-being:
> Contextual strain amid between-person stability in Japan's shrinking
> cities.* *Cities*, in revision.

[author check: final title — "shrinking cities" vs. "shrinking small and
medium-sized cities"]

This repository contains the Python pipeline used to estimate the
temporal (lag-1) dependency structure among fourteen satisfaction
domains in Japanese municipalities, stratified by city scale (Large
cities vs. Small-to-Medium Cities) and age cohort (15–64 vs. 65+), using
PCMCI (`tigramite`) with a within-person transition mask
(`mask_type="y"`), Benjamini–Hochberg multiplicity control with an
effect-size floor, a person-demeaned within-person specification, and
decline contrasts (Fisher z). The pipeline goes from the survey workbook
to every reported number: panel construction, the six self-verification
gates, estimation, and the deterministic generation of **Tables 1–7**
and **Figures 2–4**. (Figure 1, the municipal map, is produced in GIS
and is not part of the pipeline.)

## Repository structure

```
├── complete_pipeline.py    # Full analysis and reporting pipeline (v1.0.0)
├── data/
│   └── schema.csv          # Workbook layout: variable definitions and sources (no raw data)
├── .gitignore              # Prevents accidental data/output commits
├── requirements.txt
├── CITATION.cff
├── LICENSE
└── README.md
```

Running the pipeline creates `data/`, `tables/`, and `figures/` under
`--outdir`. These are not tracked in the repository.

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
both sheets the pipeline reads — `24_分析データ` (33 columns, read
positionally) and `17_Population` (11 columns) — with each variable's
type, unit, definition, source, access conditions, and whether the
pipeline uses it, so that researchers who obtain the individual data
can assemble an equivalent workbook and run the pipeline without
modification.

All pipeline outputs are **aggregate-level** (edge lists, centralities,
contrasts, counts) and safe to share; nothing person-level is written
to disk.

## Requirements

- Python ≥ 3.10
- Required packages: `numpy`, `pandas`, `openpyxl`, `matplotlib`,
  `tigramite` ≥ 5.2

```
pip install -r requirements.txt
```

## How to run

1. Obtain the individual data, assemble the workbook as documented in
   `data/schema.csv`, and note its path.
2. `python complete_pipeline.py --xlsx /path/to/workbook.xlsx --outdir results`
   — reads the workbook, builds the panel, runs Gate 1 and stops on
   failure, then estimates all 28 PCMCI models and writes the outputs.
3. Options: `--force` (proceed past a failed Gate 1), `--jitter-seed N`
   (jitter-invariance check), `--skip-figures`. From a notebook, set
   `XLSX_PATH_OVERRIDE = "/path/to/workbook.xlsx"` in a prior cell and
   paste the script into one cell.

**Runtime.** Roughly 30–90 minutes (28 PCMCI runs).

## Self-verification gates

The pipeline validates itself against the manuscript at six gates and
stops (Gate 1) or flags `CHECK` (Gates 2–6) on any mismatch:

| Gate | Checks | Against |
|---|---|---|
| 1 | sample structure, transitions, medians | Table 1 |
| 2 | out-strength sentinels | Table 5 |
| 3 | within-person edge counts (0/0/0/0 for ages 15–64; 0/1/1/4 for 65+) | Section 4.3 / Table 6 |
| 4 | Life-Satisfaction out-strength zero in 7/8 demeaned strata under the loose criterion (q < .05, \|MCI\| ≥ .05) | Section 4.3 |
| 5 | composition check: pooled 3-wave-plus models retain the Life-Satisfaction hub (0.62–0.90) | Section 4.3 |
| 6 | all 91 unoriented lag-0 pairs significant in every pooled model | Section 4.2 |

## Output ↔ article mapping

| Output | Article element |
|---|---|
| `tables/table1_sample_structure.csv` | Table 1 (with the Gate 1 pipeline-vs-manuscript comparison) |
| `tables/table2_descriptives.csv` | Table 2 |
| `tables/table3_edge_counts.csv` | Table 3 |
| `tables/table4_pooled_edges.csv` | Table 4 |
| `tables/table5_out_strength.csv` | Table 5 |
| `tables/table6_within_edges.csv` | Table 6 (65+ strata) |
| `tables/table7_contrasts.csv` | Table 7 |
| `figures/figure2.(png\|pdf)` | Figure 2 — pooled networks, all eight strata |
| `figures/figure3.(png\|pdf)` | Figure 3 — person-demeaned within-person networks (65+) |
| `figures/figure4.(png\|pdf)` | Figure 4 — decline-contrast forest plot |
| `data/edges_pooled.csv`, `data/edges_within.csv` | reported edges, full precision |
| `data/all_edges_pooled.csv`, `data/all_edges_within.csv` | all 182 lagged hypotheses per model, incl. non-significant (within file includes non-shrinking subsets) |
| `data/contrasts.csv` | decline contrasts, full precision (targeted family + age-pattern check) |
| `data/sensitivity_edge_counts.csv` | edge counts under four reporting criteria (Section 4.3) |
| `data/composition_check.csv` | pooled models on the 3-wave-plus subsample (Section 4.3) |
| `data/contemporaneous_pairs.csv` | lag-0 verification layer (Section 4.2) |
| `data/run_metadata.json` | all parameters + software versions |

## Reproducibility notes

- The transition mask restricts every lagged test to consecutive-year
  observations of the same respondent; applied to the dependent
  variable (`mask_type="y"`), it governs both MCI estimates and their
  analytic significance (manuscript Section 3.3).
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
- `data/run_metadata.json` records every parameter and the Python and
  package versions of each run.

## Citation

If you use this code, please cite the paper above and, for the exact
released version of this code, the archived DOI (Zenodo; citation
metadata for this repository is provided in `CITATION.cff`).
[author check: add the Zenodo DOI on release]

## License

MIT — see `LICENSE`. The license applies to the code only and does not
extend to the survey data, which remain subject to the Cabinet Office's
conditions of provision.

## Contact

Haruka Kato, Associate Professor, Osaka Metropolitan University
