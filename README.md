# Transcriptomic lift across every clinical baseline

![Example runs](https://github.com/simeonpetrinin/transcriptomic_clin_benchmark/actions/workflows/example.yml/badge.svg) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23786213.svg)](https://doi.org/10.5281/zenodo.23786213)

Code for *The added predictive value of transcriptomic data across every clinical baseline: a pan-cancer survival prediction benchmark* (Petrinin and Baev, in preparation).

For 21 TCGA cancers and 10 external cohorts, every combination of a prespecified set of clinical variables (14,011 clinical subsets in TCGA) is fitted with and without RNA-seq data in repeated nested cross-validation. Models: IPF-Lasso, the Volkmann offset model, IPF-ridge and an RNA-score model on the full lattice; CoxBoost and blockForest (against a clinical random survival forest) on 83 specifications per cancer (the RNA-only model, age alone, the full main set and 80 sampled subsets). Metrics: Uno's C-index, integrated Brier score, D-calibration and ICI.

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
| `results/robustness_specs.csv` | clinical subsets used in the robustness run (CoxBoost, blockForest) |
| `results/robustness_specs_gene_sweep.csv` | the smaller, earlier robustness sample used in the gene-count sweep (Figure S1) |
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

Reruns the full LAML analysis (7 clinical subsets plus the RNA-only model, 10 repeats x 5 folds) from the bundled data twice, first with IPF-Lasso and the offset model, then with IPF-ridge and the RNA-score model, and compares them with `example/LAML_expected_subset_summary.csv` and `example/LAML_expected_rs_subset_summary.csv`. For GBM, run `Rscript scripts/run_example.R GBM` from the project root. Each run prints `PASS` when every C-index and IBS value is within 0.001 of the paper's, or stops with an error. Both checks run on GitHub Actions on every push (badge above). Neither cancer needed refitting with `penalized`, so the rescue step is not exercised by these examples.

## Full reproduction

### Data

Place in `data/raw/`:

- `TCGA-CDR.xlsx` (TCGA Pan-Cancer Clinical Data Resource, Liu et al. 2018) and `clinical_PANCAN_patient_with_followup.tsv`, from the [PanCanAtlas publication page](https://gdc.cancer.gov/about-data/publications/pancanatlas)
- `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz`, from the UCSC Xena TCGA Pan-Cancer hub

External cohorts: the GEO series (GSE39582, GSE68465, GSE65858, GSE13507, GSE84437, GSE14520) are downloaded automatically. CGGA (mRNAseq_693, http://www.cgga.org.cn) and Beat AML waves 1-4 (https://biodev.github.io/BeatAML2/) are downloaded by hand into `data/raw/external/CGGA/` and `data/raw/external/BEATAML/`. METABRIC is the cBioPortal study `brca_metabric` (https://www.cbioportal.org/study/summary?id=brca_metabric); its `data_clinical_patient.txt`, `data_clinical_sample.txt` and `data_mrna_illumina_microarray.txt` go in `data/raw/external/METABRIC/brca_metabric/`.

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
| IPF-ridge and RNA score (full lattice, key `<CANCER>_rs`); CoxBoost and blockForest (robustness subsets, key `<CANCER>_robust`) | `source("scripts/run_weekend.R")`, then the printed steps |
| gene-count sweep (100 and all genes, robustness sample) | `sweep_genes_robust()` |
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
| 8 | `python figure_generation/fig8.py` |
| S1 | `python figure_generation/gene_sweep.py` |
| S2 | `python figure_generation/runtime.py` |
| S3, S4 | `python figure_generation/external_ibs_dcal.py` |
| extra outputs, not in the paper: per-cancer Shapley values, specification curves, T/N/M with and without stage | `python figure_generation/variable_importance.py` |

`common.py` holds the data loading shared by Figures 2-4.

## Code archive

This code is archived on Zenodo. Release v2.1 (the code used in the paper): [doi:10.5281/zenodo.23786213](https://doi.org/10.5281/zenodo.23786213); v2.1 is the code of v2.0.1 with the figures saved at print resolution. Release v2.0.1 (the same code as v2.0): [doi:10.5281/zenodo.23677300](https://doi.org/10.5281/zenodo.23677300). Release v1.0: [doi:10.5281/zenodo.23483562](https://doi.org/10.5281/zenodo.23483562).

## Data and results

Processed data, per-fold results and summaries for all cohorts: Zenodo, [doi:10.5281/zenodo.23651249](https://doi.org/10.5281/zenodo.23651249).

## Citation

Petrinin S, Baev M. The added predictive value of transcriptomic data across every clinical baseline: a pan-cancer survival prediction benchmark. In preparation.

## License

MIT (see `LICENSE`). The data are subject to the terms of their original sources (TCGA, GEO, CGGA, Beat AML).
