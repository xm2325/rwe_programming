# RWE Programming Portfolio

Auditable observational-research programming portfolio with two deliberately separated evidence layers: a fully reproducible **synthetic longitudinal RWE study package** for end-to-end source-to-outcome programming, and a **public CDC NHANES real-data propensity-score layer** for treatment-model, overlap and balance validation on real participants. No proprietary patient data are used.

## Provenance

This is a **restored/reconstructed** version of an earlier 2026 portfolio project prepared for observational-research programming roles. Archived application materials preserved detailed specifications and validation claims, but the original source archive was not found. See [`docs/RESTORATION_PROVENANCE.md`](docs/RESTORATION_PROVENANCE.md).

The repository separates **currently executable evidence** from **historical archived evidence**. It does not claim that every historical v0.5.x result has been reproduced byte-for-byte.

## Current executable workflow — v0.10.0-restored

### 1. Longitudinal source-to-analysis construction

The synthetic longitudinal study contains six source domains: `patients`, `enrollment`, `diagnoses`, `labs`, `medications` and `outcomes`.

- explicit patient-specific index date, **365-day baseline** and up to **1,825-day follow-up** bounded by enrollment;
- baseline diabetes, hypertension and gout-flare ascertainment from pre-index diagnosis events;
- last qualifying baseline eGFR and serum-urate extraction from laboratory events;
- explicit **baseline CKD exclusion derived from pre-index diagnosis history**;
- uncontrolled-gout phenotype from baseline serum urate and flare burden;
- primary `INCIDENT_CKD` outcome restricted to strictly post-index source events;
- source-derived `INCIDENT_CKD_CONFIRM` event for a stricter confirmed-CKD sensitivity phenotype;
- independently generated source-level `NEGATIVE_CONTROL` post-index outcome with no exposure term;
- longitudinal source population deliberately contains baseline-CKD-ineligible patients while the deterministic eligible analysis cohort remains **9,184**;
- independent **Python/pandas** and **SQLite SQL** cohort builders execute the same eligibility, timing and outcome logic;
- patient-level SQL↔Python reconciliation includes primary, strict-CKD and negative-control endpoints.

### 2. Source-derived propensity weighting and survival analysis

Reviewer-facing primary analyses consume the cohort built from longitudinal source tables. The old flat synthetic cohort remains only as an explicitly labelled regression/reference layer.

- propensity-score modelling and stabilised IPTW;
- weighted covariate-balance and effective-sample-size checks;
- explicit IPTW-weighted Breslow Cox partial likelihood;
- **case-weighted score-residual sandwich standard error conditional on the estimated IPTW**;
- model-based standard error retained separately;
- IPTW Kaplan–Meier point-estimate curves and fixed-time survival summaries;
- optimized cumulative-risk-set Cox implementation reconciled against the slower transparent Python implementation;
- independent unweighted Cox reconciliation against `statsmodels` PHReg.

#### Independent R weighted-Cox validation

A separate GitHub Actions workflow independently fits the same deterministic source-derived weighted analysis in R using:

`survival::coxph(..., weights=stabilized_weight, ties="breslow", robust=TRUE, cluster=patient_id)`

On the fixed 2,000-patient validation cohort:

| Quantity | Python | R `survival` | Absolute difference |
|---|---:|---:|---:|
| Cox coefficient | 0.02017184556 | 0.02017179200 | **5.36 × 10^-8** |
| Robust SE | 0.24591381945 | 0.24591382338 | **3.93 × 10^-9** |

The dedicated workflow now requires **both** coefficient and robust-SE differences to be ≤ `1e-6`.

### 3. Fully source-derived sensitivity suite

The same longitudinal analysis cohort now drives all principal sensitivity analyses:

- 1st/99th-percentile weight trimming;
- induced missingness with complete-case and median-imputed analyses;
- bootstrap HR uncertainty;
- weighted risk-set proportional-hazards screen;
- **source-derived negative-control outcome analysis**;
- **source-derived stricter/confirmed CKD phenotype analysis**.

The reviewer report no longer uses the legacy flat cohort for these sensitivities.

### 4. Production programming controls

- parameterised `StudyConfig` with study ID/version, population, exposure/comparator, outcome and thresholds;
- **18-check longitudinal source/time-zero QC registry** covering keys, enrollment, baseline laboratory availability, medication timing, outcome timing, follow-up, exposure reproducibility, baseline CKD exclusion, strict-CKD consistency and negative-control availability;
- a separate **14-check legacy/reference runtime QC registry**, retained only for regression/provenance checks;
- source-derived analysis-dataset specification and data dictionary;
- Table 1 balance, Table 2 outcomes and Table 3 primary effect outputs;
- robust and model-based Cox uncertainty reported separately;
- end-to-end tests protecting `longitudinal source → eligible cohort → IPTW → Cox/KM → sensitivity → TLF` lineage.

### 5. Reviewer-ready delivery

`build_reviewer_bundle()` assembles a traceable package including:

- study configuration and analysis specification;
- explicitly labelled legacy/reference attrition and runtime-QC files;
- **source-derived** baseline/balance, outcome and primary-effect tables;
- **source-derived** weighted Kaplan–Meier curve and fixed-time survival summary;
- all six longitudinal source-domain CSVs plus the derived eligible analysis cohort;
- longitudinal source inventory and SQL↔Python builder reconciliation;
- 18-check longitudinal QC manifest;
- `source_derived_sensitivity.json`;
- independent Cox reference results;
- self-contained HTML validation report;
- ADRG-like [`docs/REVIEWER_GUIDE.md`](docs/REVIEWER_GUIDE.md), copied into the bundle;
- machine-readable bundle inventory.

The normal `RWE validation` GitHub Actions workflow runs tests, builds the reviewer bundle and uploads it as an artifact.

## Public real-data PS validation — CDC NHANES 2017–2018

A dedicated workflow downloads official CDC NHANES XPT files, verifies the SAS XPORT header, harmonises demographic/questionnaire/laboratory/prescription components, fits treatment propensity models and publishes **aggregate evidence only**. Participant-level analysis rows are removed before artifact upload.

Population: adults reporting doctor-diagnosed gout.

Exposure: current urate-lowering therapy (ULT) identified from prescription generic names: allopurinol, febuxostat, probenecid or pegloticase.

Default PS covariates:

- age;
- sex;
- race/ethnicity;
- BMI;
- diabetes;
- hypertension;
- serum creatinine.

Current serum urate is retained descriptively but deliberately excluded from the treatment PS because current laboratory measurement may be downstream of current ULT in this cross-sectional dataset.

### Executed real-data results

The validated complete-case analysis contains **291** adults with doctor-diagnosed gout: **86** current ULT users and **205** untreated participants.

- fitted ordinary PS range: **0.0245–0.9132**;
- observed treated/untreated common-support interval: **0.0365–0.6433**;
- unweighted max |SMD|: **0.3637**.

| Weighting analysis | Max absolute SMD | ESS |
|---|---:|---:|
| Stabilised IPTW | **0.1799** | **225.5** |
| PS overlap weighting | **0.000096** | **214.7** |
| Survey-weighted PS × IPTW | **0.1392** | **109.9** |
| Survey-weighted PS × overlap | **0.000039** | **83.2** |

The ordinary IPTW result is intentionally retained even though its maximum absolute SMD remains above a conventional 0.10 balance target. **Overlap weighting is a separate sensitivity that targets the common-support population; it is not the same estimand and is not presented as a post-hoc replacement for IPTW.**

The NHANES layer makes no longitudinal or causal treatment-effect claim. Machine-readable diagnostics explicitly retain `causal_effect_claim=false`. See [`docs/REAL_DATA_PS.md`](docs/REAL_DATA_PS.md).

## Run locally

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

Run longitudinal QC:

```bash
PYTHONPATH=src python - <<'PY'
from rwe_programming import longitudinal_qc_manifest
print(longitudinal_qc_manifest())
PY
```

Run the source-derived primary analysis:

```bash
PYTHONPATH=src python - <<'PY'
from rwe_programming import (
    build_analysis_cohort_python,
    make_longitudinal_source_population,
    summarise_analysis,
)
sources = make_longitudinal_source_population()
cohort = build_analysis_cohort_python(sources)
print(summarise_analysis(cohort))
PY
```

Build the reviewer bundle:

```bash
PYTHONPATH=src python - <<'PY'
from rwe_programming import build_reviewer_bundle
print(build_reviewer_bundle("artifacts/reviewer_bundle", n_boot=50))
PY
```

For public real-data validation, use the dedicated `NHANES real-data PS validation` workflow. For independent weighted survival validation, use `Independent R weighted Cox reconciliation`.

## Validation policy

The GitHub Actions runs for the current code are the authoritative executable evidence. Key controls include:

- patient-level reconciliation between independent longitudinal SQL and Python cohort builders;
- source-level baseline CKD exclusions with zero excluded-patient leakage;
- strict-CKD and negative-control source-event consistency;
- source-derived primary TLF and sensitivity lineage;
- optimized weighted-Cox risk-set calculations reconciled against a transparent Python reference;
- independent R reconciliation of the weighted Cox **coefficient and robust SE**;
- weighted Kaplan–Meier monotonicity/range checks;
- balance and ESS monitoring before/after weighting and trimming;
- public real-data PS common-support, weight and balance diagnostics;
- reviewer-bundle generation as a CI artifact.

The fixed restoration seed uses **9,184 synthetic eligible analysis patients**. Exact executable results should be read from current code and workflow artifacts rather than copied from archived application prose.

## Historical archived evidence

Archived application artefacts recorded an earlier release with **52 tests and 60 runtime QC checks**, plus SQL↔pandas reconstruction, OMOP-shaped checks, Cox reconciliation, missingness/weight-restriction sensitivity, bootstrap uncertainty, PH diagnostics, negative controls and outcome/phenotype sensitivity.

Those historical counts remain provenance evidence only. They are not presented as the current executable counts.

## Why this maps to observational-research programming

The repository demonstrates tasks common in RWE/RWD programming: translating cohort/exposure/outcome specifications into longitudinal source-table logic, defining time zero and observation windows, deriving exclusions from source history, independently validating SQL/Python cohort construction, estimating and auditing propensity-score weights, implementing and independently reconciling weighted survival analysis, running source-derived sensitivity analyses, producing reviewer-facing TLFs, and packaging traceable QC evidence. The separate NHANES layer additionally shows that the PS/weighting diagnostics can operate on real public epidemiological data rather than synthetic rows alone.

## Data policy

No Amgen, CPRD, MarketScan, Optum, Komodo, IQVIA, Symphony, or other proprietary/identifiable healthcare data are included. Longitudinal patient-level data generated by this repository are synthetic. Public NHANES participant data are downloaded transiently by the dedicated workflow; raw participant-level rows are excluded from published workflow artifacts and from version control.
