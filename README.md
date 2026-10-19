# Transcriptomic lift across every clinical baseline

![Example runs](https://github.com/simeonpetrinin/transcriptomic_clin_benchmark/actions/workflows/example.yml/badge.svg) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23483562.svg)](https://doi.org/10.5281/zenodo.23483562)

Code for *The added predictive value of transcriptomic data across every clinical baseline: a pan-cancer survival prediction benchmark* (Petrinin and Baev, in preparation).

For 21 TCGA cancers and 9 external cohorts, every combination of a prespecified set of clinical variables (14,011 clinical subsets in TCGA) is fitted with and without RNA-seq data in repeated nested cross-validation. Models: IPF-Lasso and the Volkmann offset model on the full lattice; CoxBoost, blockForest, IPF-ridge and an RNA-score model on 12 subsets per cancer. Metrics: Uno's C-index, integrated Brier score, D-calibration and ICI.

## Contents

| Path | |
|---|---|
| `master.R` | entry point for the TCGA analysis |
| `external_master.R` | entry point for the external cohorts |
| `paths.R`, `settings.txt` | file paths; number of parallel workers |
| `R/` | data loading, clinical subsets, models, metrics, cross-validation |
| `scripts/` | summaries, rescue run, robustness run, sensitivity analyses, audits |
| `figure_generation/` | Python scripts for Figures 2-7 and Supplementary Figures S1-S4 |
| `example/` | the paper's LAML and GBM results, used by the example runs |
| `data/processed/LAML_modelling.rds`, `GBM_modelling.rds` | modelling data for the example runs |
| `results/robustness_specs.csv` | clinical subsets used in the robustness run |
| `renv.lock`, `requirements.txt` | R and Python package versions |

## Requirements

- R 4.6.1. Package versions are pinned in `renv.lock`; restore them with `renv::restore()`. `master.R` stops if the library does not match the lockfile.
- Python 3.13.9 for the figures: `python -m pip install -r requirements.txt`. Figures 3-7, S3 and S4 use the Inter font (https://rsms.me/inter/); without it, matplotlib falls back to DejaVu Sans and the layout shifts slightly.
- Set the number of parallel workers in `settings.txt`.

## Example runs (TCGA-LAML and TCGA-GBM, a few minutes each)

```r
renv::restore()
source("scripts/run_example.R")
```

Reruns the full LAML analysis (7 clinical subsets plus the RNA-only model, 10 repeats x 5 folds) from the bundled data and compares it with `example/LAML_expected_subset_summary.csv`. For GBM, run `Rscript scripts/run_example.R GBM` from the project root. Each run prints `PASS` when every C-index and IBS value is within 0.001 of the paper's, or stops with an error. Both checks run on GitHub Actions on every push (badge above). Neither cancer needed refitting with `penalized`, so the rescue step is not exercised by these examples.

## Full reproduction

### Data

Place in `data/raw/`:

- `TCGA-CDR.xlsx` (TCGA Pan-Cancer Clinical Data Resource, Liu et al. 2018) and `clinical_PANCAN_patient_with_followup.tsv`, from the [PanCanAtlas publication page](https://gdc.cancer.gov/about-data/publications/pancanatlas)
- `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz`, from the UCSC Xena TCGA Pan-Cancer hub

External cohorts: the GEO series (GSE39582, GSE68465, GSE65858, GSE13507, GSE84437, GSE14520) are downloaded automatically. CGGA (mRNAseq_693, http://www.cgga.org.cn) and Beat AML waves 1-4 (https://biodev.github.io/BeatAML2/) are downloaded by hand into `data/raw/external/CGGA/` and `data/raw/external/BEATAML/`. `R/external_data.R` also declares SCAN-B, METABRIC and E-MTAB-1980, which are not part of the paper.

### TCGA main analysis

```r
source("master.R")
run_audit()        # cancer and clinical-variable selection
run_pipeline()     # prepare, cross-validate and summarise every cancer
```

Single stages: `prepare_cancer("LUAD")`, `run_cv("LUAD")`, `summarise("LUAD")`. The full lattice takes several days on 6 cores.

### Rescue of failed glmnet fits

```r
source("master.R")
source("scripts/run_rescue.R")       # BRCA, CESC, COAD, KIRP, READ, SARC, UCEC
source("scripts/run_rescue_kirc.R")  # KIRC
source("scripts/rebuild_summaries.R"); rebuild_summaries()
```

Writes `<CANCER>_rescue.csv` and rebuilds the summaries with the rescued fits; the `_norescue` files are the sensitivity analysis.

### Secondary analyses

| Analysis | Commands (after `source("master.R")`) |
|---|---|
| robustness (other models) | `full_run_robustness()` |
| strongest and weakest baselines chosen inside the folds | `source("scripts/run_best_worst_subset_clin.R"); full_run_baseline_selection()` |
| partial likelihood and ICI | `source("scripts/run_loglik_ici.R"); full_run_loglik_ici()` |
| published baselines | `source("scripts/build_published_baselines.R"); source("scripts/run_published_baselines.R"); full_run_published()` |
| Volkmann ladder (BRCA) | `source("scripts/run_volkmann_ladder.R"); run_volkmann()` |
| external cohorts | `source("external_master.R"); run_external_panel()` |

## Figures

Run from the project root after the analyses above. Each script reads `results/` and writes to `figures/`.

| Figure | Script |
|---|---|
| 1 | drawn by hand, no script |
| 2 | `python figure_generation/fig2.py` |
| 3 | `python figure_generation/fig3.py` |
| 4 | `python figure_generation/fig4.py` |
| 5 | `python figure_generation/fig5.py` |
| 6 | `python figure_generation/fig6.py` |
| 7 | `python figure_generation/fig7.py` |
| S1 | `python figure_generation/gene_sweep.py` |
| S2 | `python figure_generation/runtime.py` |
| S3, S4 | `python figure_generation/external_ibs_dcal.py` |
| extra outputs, not in the paper: per-cancer Shapley values, specification curves, T/N/M with and without stage | `python figure_generation/variable_importance.py` |

`common.py` holds the data loading shared by Figures 2-4.

## Code archive

This code (release v1.0) is archived on Zenodo: [doi:10.5281/zenodo.23483562](https://doi.org/10.5281/zenodo.23483562).

## Data and results

Processed data, per-fold results and summaries for all cohorts: Zenodo, [doi:10.5281/zenodo.23479392](https://doi.org/10.5281/zenodo.23479392).

## Citation

Petrinin S, Baev M. The added predictive value of transcriptomic data across every clinical baseline: a pan-cancer survival prediction benchmark. In preparation.

## License

MIT (see `LICENSE`). The data are subject to the terms of their original sources (TCGA, GEO, CGGA, Beat AML).
