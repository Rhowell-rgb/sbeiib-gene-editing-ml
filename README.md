# SBEIIb gene editing + machine learning

Code for a layered machine-learning workflow to select CRISPR/Cas9-edited *SBEIIb* rice lines that balance
nutritional quality (low glycaemic index, high resistant starch) with grain texture.

> **Status: work in progress / preliminary.** The data are unpublished and are not included in this repository.
> Results, thresholds and conclusions may change.

## Study design (summary)
- Three sgRNA target sites in *SBEIIb*: promoter DOF motifs, exon-15 boundary, near the stop codon.
- Three elite backgrounds: IR64, IRRI154 (NSIC Rc222) and Samba Mahsuri; T0 to T4 and backcross generations.
- Early phenotypes (small seed numbers): size-exclusion chromatography (SEC) chain-length profile and texture profile analysis (TPA).
- Later phenotypes: in vitro GI, resistant starch (RS), digestible carbohydrate (DC); RVA, protein and amino acids where available.

## Repository layout
```
python/    build_master.py, make_workbook.py    Data cleaning -> ML-ready master table (+ QC log)
matlab/    S00-S06 scripts, RUN_ALL.m, helpers/ Analysis pipeline
data/      (empty; see data/README.md)
docs/      data_dictionary.csv                  Column definitions
results/   (created by the scripts; git-ignored)
```

## Requirements
- Python 3.9+ with `pip install -r requirements.txt`
- MATLAB R2022a+ with the Statistics and Machine Learning Toolbox (uses `fitrensemble`, `fitrgp`, `fitrnet`, `lasso`, `shapley`, `fitlme`)

## Quick start
```bash
# 1. Build the master table (needs the compiled workbook; see data/README.md)
python python/build_master.py path/to/compiled_SBEIIb_modeling.xlsx data
```
```matlab
% 2. In MATLAB: set the current folder to matlab/ and run the steps one at a time
S01_prepare_data
S02_layer1_genotype_to_phenotype
S03_layer2_nutrition
S04_shap_interpretation
S05_desirability_ranking
S06_mixed_model_baseline
% or:  RUN_ALL
```
Settings (features, targets, CV, desirability thresholds) are in `matlab/S00_config.m`.
For a fast test run set `cfg.nRepeats = 1` and `cfg.models = {'EN','RF'}`.

| Script | Purpose | Main outputs (`results/`) |
|---|---|---|
| `S01_prepare_data` | Load master table, express traits relative to the background control median, CLR-transform SEC fractions, build feature matrices, fix grouped CV folds | `prepared.mat` |
| `S02_layer1_genotype_to_phenotype` | Layer 1: edit features to SEC/TPA; elastic net, random forest, Gaussian process, shallow ANN vs a background+construct baseline | `L1_cv_results.csv`, plots |
| `S03_layer2_nutrition` | Layer 2: GI and RS from genotype, early phenotypes, both, and both + stacked Layer-1 predictions | `L2_cv_results.csv` |
| `S04_shap_interpretation` | Random forest + SHAP + permutation importance | `L2_SHAP_*.csv/png` |
| `S05_desirability_ranking` | Desirability scores (low GI / high RS / texture / protein), Pareto set, line/event/genotype-group rankings | `S05_*.csv/png` |
| `S06_mixed_model_baseline` | Linear mixed model: trait ~ allele doses + background + generation + (1\|T0 event) | `S06_LMM_dose_effects.csv` |

## Methodological rules built into the code
1. **Grouped cross-validation by T0 event (`CVGroup`).** All progeny of one T0 plant stay in the same fold; random K-fold would leak siblings between training and test sets.
2. **Out-of-fold stacking.** Layer-1 predictions passed to Layer 2 come from models that never saw that row.
3. **Preprocessing inside the training fold** (imputation, scaling).
4. **Traits modelled as continuous values relative to the background control median;** scenario groups (e.g. low GI + texture) are derived afterwards with desirability scores, avoiding labels built from the model inputs.
5. **Features excluded on purpose:** GBSS SNPs (identical within each background, hence confounded) and raw residue-loss columns (constant within a construct); dose-weighted domain-loss features are used instead.
6. **Generation is a confounded covariate for some traits** (assay batch/season vs biology). Set `cfg.includeGeneration` to `true`/`false` in `S00_config.m` and report both.
7. **Model performance is judged against a background + construct baseline,** not against zero.

## Reproducibility
Random seed and CV settings are in `S00_config.m` (`cfg.seed`, `cfg.nFolds`, `cfg.nRepeats`).
`build_master.py` writes a QC log listing every data-quality item it detects (mismatched IDs, duplicated records, suspicious controls, unresolved genotype calls).

## License
MIT (see `LICENSE`). Fill in the copyright holder before publishing; check your institution's policy on code licensing.
