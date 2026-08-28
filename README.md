# RWE Programming Portfolio

Auditable observational-research programming portfolio built around a fully synthetic longitudinal gout cohort. The repository demonstrates cohort construction, propensity-score weighting, survival analysis, SQL/Python traceability, OMOP-shaped reconstruction, sensitivity analysis, independent statistical reconciliation, production QC and reproducible reviewer deliverables without using proprietary patient data.

## Provenance

This is a **restored/reconstructed** version of an earlier 2026 portfolio project prepared for observational-research programming roles. Archived application materials preserved detailed specifications and validation claims, but the original source archive was not found. See [`docs/RESTORATION_PROVENANCE.md`](docs/RESTORATION_PROVENANCE.md).

The repository therefore separates **currently executable evidence** from **historical archived evidence**. It does not claim that every historical v0.5.x result has been reproduced byte-for-byte.

## Current executable workflow

The repository now covers four linked layers.

### 1. Cohort and statistical analysis

- deterministic synthetic cohort of **9,184 patients**;
- controlled-vs-uncontrolled gout exposure;
- propensity-score modelling and stabilised IPTW;
- weighted balance and effective-sample-size checks;
- weighted Cox proportional-hazards model for incident CKD;
- an independently implemented Breslow partial-likelihood Cox estimator for coefficient reconciliation.

### 2. RWE validation and sensitivity analysis

- SQLite SQL ↔ pandas patient-level reconciliation;
- source-schema ↔ OMOP-shaped reconstruction;
- 1st/99th-percentile weight trimming;
- induced missingness with complete-case and median-imputed analyses;
- bootstrap uncertainty;
- Schoenfeld-residual proportional-hazards screen;
- synthetic negative-control outcome;
- alternative CKD outcome definition.

### 3. Production programming controls

- parameterised `StudyConfig` with study ID/version, population, exposure/comparator, outcome and thresholds;
- machine-readable cohort attrition ledger;
- **12-check runtime QC registry** with PASS/FAIL status and SHA-256 content fingerprint;
- analysis-dataset specification/data dictionary;
- analysis specification plus Table 1 balance, Table 2 outcomes and Table 3 primary treatment-effect outputs.

### 4. Reviewer-ready delivery

`build_reviewer_bundle()` assembles a traceable package containing:

- study configuration;
- analysis specification;
- cohort attrition table;
- analysis-data dictionary;
- baseline/balance, outcome and effect tables;
- independent Cox reconciliation results;
- runtime QC sign-off manifest;
- self-contained HTML validation report;
- machine-readable bundle inventory.

GitHub Actions runs the validation suite on every push/PR, builds the reviewer bundle and uploads it as the `rwe-reviewer-bundle` workflow artifact.

## Run

```bash
python -m pip install -r requirements.txt
PYTHONPATH=src pytest -q
```

Run the primary workflow:

```bash
PYTHONPATH=src python - <<'PY'
from rwe_programming import run_pipeline
print(run_pipeline())
PY
```

Build the reviewer bundle:

```bash
PYTHONPATH=src python - <<'PY'
from rwe_programming import build_reviewer_bundle
print(build_reviewer_bundle("artifacts/reviewer_bundle", n_boot=50))
PY
```

## Validation policy

The GitHub Actions run for the current commit is the authoritative executable status. Key controls include:

- **0 patient-level discrepancies** between SQLite SQL and pandas cohort reconstruction;
- **0 patient-level discrepancies** after OMOP-shaped decomposition/reconstruction;
- post-IPTW balance checks using absolute standardised mean differences;
- effective-sample-size monitoring before and after weight trimming;
- independent Cox coefficient reconciliation;
- deterministic cohort attrition and analysis-table reconciliation;
- runtime QC manifest generation;
- missing-data, bootstrap, PH-screen, negative-control and alternative-outcome analyses;
- reviewer-bundle generation as a CI artifact.

The fixed restoration seed uses **9,184 synthetic patients**. The earlier restoration baseline produced a maximum absolute weighted SMD of about **0.0034**; exact executable outputs should be read from the current code and workflow artifacts rather than copied from archived application prose.

## Historical archived evidence

Archived application artefacts document that the earlier portfolio version included SQL↔pandas patient-level reconstruction, source-schema↔OMOP-shaped checks, custom Cox-vs-`statsmodels` coefficient reconciliation, missing-data and weight-restriction sensitivity, bootstrap uncertainty, proportional-hazards diagnostics, negative-control analysis and outcome/phenotype sensitivity. They recorded **52 tests and 60 runtime QC checks** in a prior release.

Those historical counts remain provenance evidence only. This restoration deliberately does not present them as current executable counts.

## Why this maps to observational-research programming

The project is organised around tasks common in RWE/RWD programming: translating cohort/exposure/outcome specifications into analysis datasets, documenting attrition, checking time zero and follow-up logic, implementing propensity-score methods and survival models, independently validating key estimates, identifying data/programming discrepancies, producing analysis tables, documenting sensitivity analyses and delivering reviewer-ready QC evidence.

## Data policy

No Amgen data, CPRD data, MarketScan/Optum/Komodo/IQVIA/Symphony claims data, or other proprietary/identifiable healthcare data are included. All patient-level data generated by this repository are synthetic.
