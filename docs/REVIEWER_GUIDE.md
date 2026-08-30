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

UCG is derived from baseline serum urate and/or baseline flare burden. This creates an important distinction between **variables that define the exposure phenotype** and **pre-exposure covariates used to model propensity**.

The prespecified synthetic PS adjustment set is:

- age;
- sex indicator;
- baseline diabetes;
- baseline hypertension;
- baseline eGFR.

`baseline_urate` and `prior_flares` remain visible in Table 1 because they are clinically informative baseline characteristics, but they are **not included in the PS model and are not subject to the post-weighting balance QC gate**, because they define UCG itself. Requiring them to balance would amount to conditioning the propensity model on the exposure definition and creates a deterministic-classification/positivity problem rather than a valid confounding adjustment strategy.

The propensity model produces stabilised IPTW and balance QC is applied to the five-variable PS adjustment set using standardised mean differences (SMDs) and effective sample size (ESS). Propensity scores are clipped only at configured numerical guardrails; weight trimming is reported separately as a sensitivity analysis rather than silently changing the primary estimand.

The old flat synthetic dataset is retained only as a regression/reference fixture. Its exposure assignment is stochastic with overlap and is deliberately distinct from the reviewer-facing longitudinal phenotype definition.

## 6. Outcomes

### Primary outcome

`INCIDENT_CKD` — first qualifying post-index CKD outcome event.

### Stricter alternative phenotype

`INCIDENT_CKD_CONFIRM` — a later source event generated only after a primary CKD event and at least 30 days later. This supports a source-derived stricter/confirmed CKD sensitivity definition.

### Negative-control outcome

`NEGATIVE_CONTROL` — an independently generated post-index outcome process with no exposure term. It is used to test whether the weighting/analysis workflow creates a strong spurious exposure association in an outcome deliberately generated without an exposure effect.

## 7. Primary analysis and independent Cox reconciliation

Reviewer-facing primary outputs are generated from the longitudinal source-derived cohort:

- propensity-score estimation and stabilised IPTW;
- weighted balance on the PS adjustment set and ESS;
- reviewer-visible descriptive values for phenotype-defining baseline urate and flare burden;
- explicit IPTW-weighted Breslow Cox partial likelihood;
- case-weighted score-residual sandwich standard error conditional on the estimated IPTW;
- model-based standard error retained separately;
- weighted Kaplan–Meier point estimates and fixed-time survival summaries;
- Table 1 balance/descriptives, Table 2 outcomes and Table 3 primary effect.

The custom weighted Cox risk-set implementation is reconciled numerically against a slower transparent Python reference implementation. A separate GitHub Actions workflow also performs an independent cross-language validation using R `survival::coxph(..., weights=..., ties="breslow", robust=TRUE, cluster=patient_id)` on a deterministic 2,000-patient source-derived analysis cohort.

The validated reconciliation is extremely close:

- Python coefficient: `0.02017184556180873`;
- R coefficient: `0.0201717919974732`;
- absolute coefficient difference: **5.36 × 10^-8**;
- Python robust SE: `0.24591381945097987`;
- R robust SE: `0.245913823383128`;
- absolute robust-SE difference: **3.93 × 10^-9**.

The workflow treats both coefficient and robust-SE reconciliation at `1e-6` tolerance as hard validation gates. A separate unweighted implementation remains reconciled against `statsmodels` PHReg as an additional reference check.

## 8. Source-derived sensitivity suite

The same longitudinal analysis cohort drives:

- 1st/99th-percentile weight-trimming sensitivity;
- induced-missingness complete-case versus median-imputation sensitivity;
- bootstrap hazard-ratio uncertainty;
- weighted risk-set proportional-hazards screen;
- source-derived negative-control outcome analysis;
- source-derived stricter/confirmed CKD phenotype analysis.

These are no longer driven by the legacy flat synthetic cohort. Balance checks after weight trimming use the same prespecified PS adjustment set as the primary propensity model.

## 9. QC and reconciliation

The longitudinal QC registry contains **18 checks** covering source keys, enrollment coverage, baseline laboratory availability, medication timing, post-index outcome timing, follow-up bounds, exposure reproducibility, source-derived baseline CKD exclusion, strict CKD phenotype consistency and negative-control outcome availability.

The package also retains clearly labelled legacy/reference checks such as flat-cohort SQL↔pandas and OMOP-shaped reconstruction. These are regression/reference evidence only and should not be interpreted as the primary source-to-analysis path.

## 10. Public NHANES real-data propensity-score layer

A separate workflow downloads CDC NHANES 2017–2018 XPT files from the official public DataFiles endpoint and validates the SAS XPORT signature before analysis.

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

The currently executed public-data complete-case analysis contains **291** adults with doctor-diagnosed gout: **86** current ULT users and **205** untreated participants. The ordinary stabilised-IPTW analysis reduced maximum absolute SMD from **0.3637** to **0.1799**, which is transparently reported as residual imbalance above a conventional 0.10 target. An explicitly separate propensity-score overlap-weighting sensitivity, which changes the target to the common-support population, produced maximum absolute SMD **0.000096** with ESS **214.7**. The survey-weighted PS analyses produced maximum absolute SMD **0.1392** for survey×IPTW and **0.000039** for survey×overlap, with corresponding ESS **109.9** and **83.2**.

These analyses are reported side by side; overlap weighting is not presented as the same estimand as IPTW or as a post-hoc replacement for an inconvenient IPTW balance result.

The real-data workflow produces aggregate evidence only for publication:

- complete-case/treatment counts;
- propensity-score range and common support;
- unweighted/IPTW/overlap/survey-weighted balance table, including race/ethnicity indicator levels;
- PS overlap bins;
- weight summaries and ESS;
- machine-readable real-data QC manifest.

The CI artifact removes participant-level rows before upload. NHANES is used for treatment-model diagnostics only; the workflow explicitly sets `causal_effect_claim=false`.

## 11. Important limitations

- Synthetic longitudinal results are not clinical findings and must never be presented as Amgen or proprietary RWD results.
- NHANES 2017–2018 is cross-sectional for this use case and does not support the longitudinal CKD treatment-effect analysis implemented in the synthetic study layer.
- The synthetic UCG comparison is a methods/programming estimand; because UCG is a phenotype rather than an intervention, causal-treatment language should be avoided.
- The OMOP check is OMOP-shaped reconstruction, not certification of full OMOP CDM compliance.
- The PH diagnostic is a screening diagnostic, not a complete formal weighted proportional-hazards testing framework.
- The primary sandwich variance is conditional on the estimated IPTW and does not propagate propensity-model estimation uncertainty.
- The independent R reconciliation validates the implemented weighted Breslow coefficient and case-weighted robust sandwich SE for the deterministic validation dataset; it does not by itself validate the scientific exposure model or resolve uncertainty from PS estimation.

## 12. What a reviewer should run first

1. Run `pytest -q` with `PYTHONPATH=src`.
2. Build the reviewer bundle with `build_reviewer_bundle()`.
3. Inspect `longitudinal_qc_manifest.json` and `longitudinal_builder_reconciliation.json`.
4. Inspect Table 1/2/3, weighted survival outputs and `source_derived_sensitivity.json`.
5. Run the dedicated `Independent R weighted Cox reconciliation` workflow and inspect its aggregate reconciliation artifact.
6. Run the dedicated `NHANES real-data PS validation` workflow and inspect its aggregate real-data artifact.
