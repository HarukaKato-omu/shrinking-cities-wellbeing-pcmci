# Population decline and the temporal structure of well-being — reproduction pipeline

Complete, self-verifying analysis and reporting pipeline for:

> Kato, H. *Population decline and the temporal structure of well-being:
> Contextual strain amid between-person stability in Japan's shrinking
>  (in revision).

The pipeline goes from the Cabinet Office survey workbook to every
reported number: panel construction, PCMCI estimation with a
transition mask (`mask_type="y"`), Benjamini–Hochberg multiplicity
control with an effect-size floor, a person-demeaned within-person
specification, decline contrasts (Fisher z), and the deterministic
generation of **Tables 1–7** and **Figures 2–4**. (Figure 1, the
municipal map, is produced in GIS and is not part of the pipeline.)

## Data access

The individual-level workbook (Survey on Satisfaction and Quality of
Life, Cabinet Office of Japan, waves 2019–2025) is provided for
approved academic use and is **not** distributed with this repository.
Researchers can apply for access through the Cabinet Office of Japan's
designated data-use procedures. Municipal population counts (2017,
2023) come from the Basic Resident Register (e-Stat of Japan).

All pipeline outputs are **aggregate-level** (edge lists, centralities,
contrasts, counts) and safe to share; nothing person-level is written
to disk.

## Requirements

Python ≥ 3.10 and:

```
pip install -r requirements.txt
```

(numpy, pandas, openpyxl, matplotlib, tigramite ≥ 5.2)

## Usage

```
python complete_pipeline.py --xlsx /path/to/workbook.xlsx --outdir results
```

Options: `--force` (proceed past a failed Gate 1), `--jitter-seed N`
(jitter-invariance check), `--skip-figures`. Runtime is roughly 30–90
minutes (28 PCMCI runs).

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

## Outputs

| Path | Content |
|---|---|
| `data/edges_pooled.csv`, `data/edges_within.csv` | reported edges, full precision |
| `data/all_edges_pooled.csv`, `data/all_edges_within.csv` | all 182 lagged hypotheses per model, incl. non-significant (within file includes non-shrinking subsets) |
| `data/contrasts.csv` | decline contrasts (targeted family + specificity) |
| `data/sensitivity_edge_counts.csv` | edge counts under four reporting criteria |
| `data/composition_check.csv` | pooled models on the 3-wave-plus subsample |
| `data/contemporaneous_pairs.csv` | lag-0 verification layer |
| `data/run_metadata.json` | all parameters + software versions |
| `tables/table1…table7*.csv` | manuscript Tables 1–7 |
| `figures/figure2–4.(png\|pdf)` | manuscript Figures 2–4 |

## Reproducibility notes

* The transition mask restricts every lagged test to consecutive-year
  observations of the same respondent; applied to the dependent
  variable (`mask_type="y"`), it governs both MCI estimates and their
  analytic significance (manuscript Section 3.3).
* A negligible jitter (SD = 1e-5) breaks degenerate ties in covariate
  codes; the seed and draw order are fixed for exact stream
  reproduction, and no reported estimate depends on the seed
  (`--jitter-seed`).
* Person-demeaning is applied within each estimation subset for
  respondents observed in ≥ 3 waves (manuscript Section 3.4).

## Citing

See `CITATION.cff`. Please cite both the article and, for the exact
released version of this code, the archived DOI (Zenodo).

## License

MIT — see `LICENSE`.
