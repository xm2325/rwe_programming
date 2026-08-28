# RWE Programming Portfolio

Auditable observational-research programming portfolio built around a fully synthetic longitudinal gout cohort. The repository demonstrates source-to-analysis cohort construction, propensity-score weighting, survival analysis, SQL/Python traceability, OMOP-shaped reconstruction, sensitivity analysis, independent statistical reconciliation, production QC and reproducible reviewer deliverables without using proprietary patient data.

## Provenance

This is a **restored/reconstructed** version of an earlier 2026 portfolio project prepared for observational-research programming roles. Archived application materials preserved detailed specifications and validation claims, but the original source archive was not found. See [`docs/RESTORATION_PROVENANCE.md`](docs/RESTORATION_PROVENANCE.md).

The repository therefore separates **currently executable evidence** from **historical archived evidence**. It does not claim that every historical v0.5.x result has been reproduced byte-for-byte.

## Current executable workflow — v0.9.0-restored

The repository now covers five linked layers.

### 1. Longitudinal RWD source-to-analysis construction

- six synthetic source domains: `patients`, `enrollment`, `diagnoses`, `labs`, `medications` and `outcomes`;
- explicit index date, **365-day baseline window** and up-to-**5-year follow-up** bounded by enrollment end;
- baseline diabetes, hypertension and gout-flare ascertainment from diagnosis events;
- last observed baseline eGFR and serum-urate extraction from lab events;
- uncontrolled-gout phenotype from baseline serum urate and flare burden;
- incident CKD restricted to events strictly after time zero and before the end of follow-up;
- independent **Python/pandas** and **SQLite SQL** cohort builders;
- patient-level SQL↔Python reconciliation with a maximum numeric-difference check.

### 2. Cohort and statistical analysis

- deterministic final analysis cohort of **9,184 patients**;
- larger synthetic source population with staged age, baseline-eGFR, follow-up and exposure-classification exclusions;
- controlled-vs-uncontrolled gout exposure;
- propensity-score modelling and stabilised IPTW;
- weighted balance and effective-sample-size checks;
- explicit IPTW-weighted Breslow Cox partial likelihood for incident CKD;
- **subject-level sandwich standard error** for the primary Cox interval, with model-based SE retained for comparison;
- IPTW Kaplan–Meier point-estimate curves and 1-, 3- and 5-year survival summaries by exposure group;
- independent unweighted Breslow Cox implementation reconciled against `statsmodels` PHReg on the supported unweighted estimand.

### 3. RWE validation and sensitivity analysis

- flat-cohort SQLite SQL ↔ pandas patient-level reconciliation;
- source-schema ↔ OMOP-shaped reconstruction;
- longitudinal event-table SQL ↔ Python cohort reconciliation;
- 1st/99th-percentile weight trimming;
- induced missingness with complete-case and median-imputed analyses;
- bootstrap uncertainty;
- weighted risk-set proportional-hazards screen;
- synthetic negative-control outcome;
- alternative CKD outcome definition.

### 4. Production programming controls

- parameterised `StudyConfig` with study ID/version, population, exposure/comparator, outcome and thresholds;
- source-to-analysis cohort attrition ledger with non-zero exclusions at every eligibility stage;
- **14-check analysis/runtime QC registry** with PASS/FAIL status and SHA-256 content fingerprint;
- separate **15-check longitudinal source/time-zero QC registry** covering enrollment coverage, foreign keys, baseline laboratory availability, medication timing, post-index outcomes, observed follow-up and reproducible exposure derivation;
- analysis-dataset specification/data dictionary;
- analysis specification plus Table 1 balance, Table 2 outcomes and Table 3 primary treatment-effect outputs;
- Table 3 records robust and model-based uncertainty separately.

### 5. Reviewer-ready delivery

`build_reviewer_bundle()` assembles a traceable package containing:

- study configuration and analysis specification;
- cohort attrition and analysis-data dictionary;
- baseline/balance, outcome and primary-effect tables;
- weighted Kaplan–Meier curve and fixed-time survival summary CSVs;
- all six longitudinal source-domain CSVs plus the derived longitudinal analysis cohort;
- longitudinal source inventory and SQL↔Python builder reconciliation JSON;
- analysis/runtime QC manifest and an independent longitudinal source/time-zero QC sign-off manifest;
- independent Cox reconciliation results;
- self-contained HTML validation report;
- machine-readable bundle inventory that includes its own path.

GitHub Actions runs the validation suite on every push/PR, builds the reviewer bundle and uploads it as the `rwe-reviewer-bundle` workflow artifact.

## Run

```bash
python -m pip install -r requirements.txt
PYTHONPATH=src pytest -q
```

Build and reconcile the longitudinal cohort:

```bash
PYTHONPATH=src python - <<'PY'
from rwe_programming import reconcile_longitudinal_builders
print(reconcile_longitudinal_builders())
PY
```

Run the longitudinal source/time-zero QC registry:

```bash
PYTHONPATH=src python - <<'PY'
from rwe_programming import longitudinal_qc_manifest
print(longitudinal_qc_manifest())
PY
```

Run the primary statistical workflow:

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

- patient-level reconciliation between independently implemented longitudinal SQL and Python cohort builders;
- **0 patient-level discrepancies** between the flat SQLite SQL and pandas cohort reconstruction;
- **0 patient-level discrepancies** after OMOP-shaped decomposition/reconstruction;
- non-zero staged source-to-analysis exclusions with exact reconciliation to 9,184 final patients;
- post-IPTW balance checks using absolute standardised mean differences;
- effective-sample-size monitoring before and after weight trimming;
- sandwich and model-based Cox uncertainty outputs;
- weighted Kaplan–Meier monotonicity and range checks;
- independent Cox coefficient reconciliation;
- two independent QC manifests: analysis/runtime and longitudinal source/time-zero;
- missing-data, bootstrap, PH-screen, negative-control and alternative-outcome analyses;
- reviewer-bundle generation as a CI artifact.

The fixed restoration seed uses **9,184 synthetic analysis patients**. Exact executable outputs should be read from the current code and workflow artifacts rather than copied from archived application prose.

## Historical archived evidence

Archived application artefacts document that the earlier portfolio version included SQL↔pandas patient-level reconstruction, source-schema↔OMOP-shaped checks, custom Cox-vs-`statsmodels` coefficient reconciliation, missing-data and weight-restriction sensitivity, bootstrap uncertainty, proportional-hazards diagnostics, negative-control analysis and outcome/phenotype sensitivity. They recorded **52 tests and 60 runtime QC checks** in a prior release.

Those historical counts remain provenance evidence only. This restoration deliberately does not present them as current executable counts.

## Why this maps to observational-research programming

The project is organised around tasks common in RWE/RWD programming: translating cohort/exposure/outcome specifications into longitudinal source-table logic and analysis datasets, documenting index/baseline/follow-up windows and attrition, checking time zero and source integrity, independently validating SQL and Python cohort construction, implementing propensity-score methods and survival models, producing reviewer-facing tables and weighted survival outputs, documenting sensitivity analyses and delivering traceable QC evidence.

## Data policy

No Amgen data, CPRD data, MarketScan/Optum/Komodo/IQVIA/Symphony claims data, or other proprietary/identifiable healthcare data are included. All patient-level data generated by this repository are synthetic.
