# Reviewer Guide — Auditable Observational RWE Programming Portfolio

## 1. Scope and provenance

This repository is a restored/reconstructed portfolio implementation of an observational-research programming workflow. It does **not** contain Amgen, CPRD, MarketScan, Optum, Komodo, IQVIA, Symphony, or other proprietary patient-level data.

Two evidence layers are deliberately separated:

1. **Synthetic longitudinal study package** — deterministic patient-level source domains are generated so that cohort logic, time zero, eligibility, outcomes, propensity-score weighting, survival analysis, sensitivity analyses, SQL/Python reconciliation, QC and reviewer deliverables can be audited end to end.
2. **Public real-data propensity-score validation** — CDC NHANES 2017–2018 is used only to demonstrate a real participant-level treatment propensity model, overlap, weighting and balance diagnostics. The NHANES workflow makes no longitudinal or causal treatment-effect claim.

Historical application materials recorded larger validation counts in an earlier portfolio version. Those counts remain provenance evidence only; the current executable CI and manifests are authoritative.

## 2. Synthetic longitudinal study question

Among synthetic adults with gout who meet the analysis eligibility criteria, compare uncontrolled gout (UCG) with controlled gout (CG) for time to incident chronic kidney disease (CKD), using stabilised inverse-probability-of-treatment weighting (IPTW).

This is a programming and methods demonstration, not clinical evidence.

## 3. Source domains and data flow

The source layer contains six event/domain tables:

- `patients`
- `enrollment`
- `diagnoses`
- `labs`
- `medications`
- `outcomes`

The reviewer path is:

`source domains → time-zero/baseline/follow-up logic → eligibility → analysis cohort → PS/IPTW → weighted Cox/KM → sensitivity analyses → TLFs/QC/reconciliation`

The default deterministic analysis cohort contains 9,184 eligible synthetic patients. Additional source patients with pre-index baseline CKD are deliberately generated so that the CKD exclusion is derived from diagnosis history rather than asserted after cohort generation.

## 4. Time zero, baseline and follow-up

- Index date: patient-specific synthetic gout study index.
- Baseline lookback: 365 days before index.
- Follow-up: from strictly after index to the earliest of enrollment end or 1,825 days.
- Baseline eGFR and serum urate: last qualifying measurement in the baseline window.
- Baseline diabetes, hypertension and gout flares: derived from pre-index diagnosis events.
- Baseline CKD: derived from pre-index diagnosis history and excluded.
- Incident outcomes: accepted only when strictly post-index and within observed follow-up.

Independent pandas/Python and SQLite SQL cohort builders implement the longitudinal derivation and are reconciled patient by patient.

## 5. Exposure and propensity-score estimand

UCG is derived from baseline serum urate and/or baseline flare burden. The propensity model uses prespecified baseline covariates and produces stabilised IPTW. Balance is assessed with standardised mean differences and effective sample size (ESS).

The target is an IPTW pseudo-population comparison of UCG versus CG within this synthetic study population. Propensity scores are clipped only at the configured numerical guardrails; weight trimming is reported separately as a sensitivity analysis rather than silently changing the primary estimand.

## 6. Outcomes

### Primary outcome

`INCIDENT_CKD` — first qualifying post-index CKD outcome event.

### Stricter alternative phenotype

`INCIDENT_CKD_CONFIRM` — a later source event generated only after a primary CKD event and at least 30 days later. This supports a source-derived stricter/confirmed CKD sensitivity definition.

### Negative-control outcome

`NEGATIVE_CONTROL` — an independently generated post-index outcome process with no exposure term. It is used to test whether the weighting/analysis workflow creates a strong spurious exposure association in an outcome deliberately generated without an exposure effect.

## 7. Primary analysis

Reviewer-facing primary outputs are generated from the longitudinal source-derived cohort:

- propensity-score estimation and stabilised IPTW;
- weighted baseline balance and ESS;
- explicit IPTW-weighted Breslow Cox partial likelihood;
- subject-level sandwich standard error conditional on the estimated IPTW;
- model-based standard error retained separately;
- weighted Kaplan–Meier point estimates and fixed-time survival summaries;
- Table 1 balance, Table 2 outcomes and Table 3 primary effect.

The custom weighted Cox risk-set implementation is reconciled numerically against a slower transparent reference implementation. A separate unweighted custom Cox implementation is reconciled against `statsmodels` PHReg because PHReg does not provide the observation-weight interface needed to serve as the weighted reference.

## 8. Source-derived sensitivity suite

The same longitudinal analysis cohort now drives:

- 1st/99th-percentile weight-trimming sensitivity;
- induced-missingness complete-case versus median-imputation sensitivity;
- bootstrap hazard-ratio uncertainty;
- weighted risk-set proportional-hazards screen;
- source-derived negative-control outcome analysis;
- source-derived stricter/confirmed CKD phenotype analysis.

These are no longer driven by the legacy flat synthetic cohort.

## 9. QC and reconciliation

The longitudinal QC registry checks source keys, enrollment coverage, baseline laboratory availability, medication timing, post-index outcome timing, follow-up bounds, exposure reproducibility, source-derived baseline CKD exclusion, strict CKD phenotype consistency and negative-control outcome availability.

The package also retains clearly labelled legacy/reference checks such as flat-cohort SQL↔pandas and OMOP-shaped reconstruction. These are regression/reference evidence only and should not be interpreted as the primary source-to-analysis path.

## 10. Public NHANES real-data propensity-score layer

A separate workflow downloads CDC NHANES 2017–2018 XPT files from the official public DataFiles endpoint and validates the XPORT file signature before analysis.

Population: adults reporting doctor-diagnosed gout.

Exposure: current urate-lowering therapy (ULT) identified from the prescription-medication questionnaire using allopurinol, febuxostat, probenecid or pegloticase generic drug names.

Default treatment-PS covariates:

- age;
- sex;
- race/ethnicity;
- BMI;
- diabetes;
- hypertension;
- serum creatinine.

Current serum urate is retained as a descriptive variable but intentionally excluded from the default treatment propensity model because medication use and laboratory measurement are cross-sectional and serum urate may be downstream of current ULT.

The real-data workflow produces aggregate evidence only for publication:

- complete-case/treatment counts;
- propensity-score range and common support;
- pre/post-IPTW balance table, including race/ethnicity indicator levels;
- PS overlap bins;
- IPTW and survey×IPTW weight summaries/ESS;
- machine-readable real-data QC manifest.

The CI artifact removes participant-level rows before upload. NHANES is used for treatment-model diagnostics only; the workflow explicitly sets `causal_effect_claim=false`.

## 11. Important limitations

- Synthetic longitudinal results are not clinical findings and must never be presented as Amgen or proprietary RWD results.
- NHANES 2017–2018 is cross-sectional for this use case and does not support the longitudinal CKD treatment-effect analysis implemented in the synthetic study layer.
- The OMOP check is OMOP-shaped reconstruction, not certification of full OMOP CDM compliance.
- The PH diagnostic is a screening diagnostic, not a complete formal weighted proportional-hazards testing framework.
- The primary sandwich variance is conditional on the estimated IPTW and does not currently propagate propensity-model estimation uncertainty.
- Independent external validation of the weighted Cox robust standard error is a planned validation extension; current independent `statsmodels` reconciliation is limited to the supported unweighted estimand.

## 12. What a reviewer should run first

1. Run `pytest -q` with `PYTHONPATH=src`.
2. Build the reviewer bundle with `build_reviewer_bundle()`.
3. Inspect `longitudinal_qc_manifest.json` and `longitudinal_builder_reconciliation.json`.
4. Inspect Table 1/2/3, weighted survival outputs and `source_derived_sensitivity.json`.
5. For public real-data validation, run the dedicated `NHANES real-data PS validation` GitHub Actions workflow and inspect its aggregate artifact.
