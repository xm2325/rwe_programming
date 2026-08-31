# Interview Guide — RWE Programming Portfolio

This guide is for explaining the repository in an observational-research programming, RWE/RWD, epidemiology-programming or statistical-programming interview. The goal is not to show every file. The goal is to demonstrate that the study question can be translated into auditable source-to-analysis logic, that the statistical implementation is validated, and that limitations are stated precisely.

## The three pieces of evidence to show first

If there is only time to open three outputs, use these:

1. **`longitudinal_builder_reconciliation.json`** — demonstrates that independent Python/pandas and SQLite SQL cohort builders produce the same patient-level analysis cohort from the longitudinal source domains. This is the strongest programming/QC evidence.
2. **`table1_balance.csv` + `table3_primary_effect.csv`** — demonstrates the source-derived PS/IPTW analysis, makes the PS adjustment set explicit, keeps phenotype-defining variables descriptive, and reports the weighted Cox effect with robust and model-based uncertainty separately.
3. **`source_derived_sensitivity.json`** — demonstrates that trimming, missingness, bootstrap uncertainty, PH screening, negative control and stricter CKD phenotype analyses all run from the same longitudinal source-derived cohort rather than from a disconnected toy dataset.

For an interviewer who wants independent statistical validation rather than another reviewer table, replace item 3 with the dedicated **Independent R weighted Cox reconciliation** artifact. For an interviewer who asks whether anything uses real participant data, show the separate **NHANES real-data PS validation** artifact.

## 30-second opening

> I built this as an auditable observational-research programming portfolio. There are two deliberately separate evidence layers. The first is a synthetic longitudinal study package where I can expose the full source-to-analysis lineage: source tables, eligibility, time zero, propensity weighting, weighted survival analysis, sensitivities, QC and reviewer outputs. The second uses public CDC NHANES participants only to validate the propensity-score, overlap and balance workflow on real data. I keep those layers separate so I do not claim a real longitudinal treatment effect from a cross-sectional public dataset.

## Five-minute version

### 0:00–0:45 — What problem does the repository solve?

The synthetic study asks how to compare uncontrolled versus controlled gout for time to incident CKD in an observational-style longitudinal cohort. It is a programming/methods demonstration, not clinical evidence.

The important design choice is that the analysis does not begin from a ready-made modelling table. It begins from six source domains:

`patients → enrollment → diagnoses → labs → medications → outcomes`

### 0:45–1:45 — Show the source-to-analysis chain

Explain that the program derives:

- patient-specific index dates;
- a 365-day baseline window;
- observed follow-up up to 1,825 days;
- diabetes, hypertension and gout flares from diagnosis history;
- baseline eGFR and serum urate from laboratory records;
- baseline CKD exclusion from pre-index diagnosis history;
- incident CKD only from strictly post-index events within observed enrollment.

Then show **`longitudinal_builder_reconciliation.json`**. The key point is not simply that SQL exists; it is that the SQL and Python implementations are independent cohort builders and are reconciled patient by patient.

### 1:45–2:45 — Explain the propensity-score decision

The UCG phenotype is defined partly by baseline serum urate and flare burden. Therefore those variables are visible in Table 1 but are not placed into the propensity model as if they were ordinary pre-exposure confounders.

The prespecified synthetic PS adjustment set is:

- age;
- sex;
- diabetes;
- hypertension;
- baseline eGFR.

This avoids a deterministic-classification/positivity problem caused by conditioning the treatment model on variables that define the exposure itself. Balance QC applies to the PS adjustment set; urate and flare burden remain reviewer-visible descriptive variables.

### 2:45–3:45 — Explain weighted survival analysis and validation

The primary analysis uses stabilised IPTW and an explicit weighted Breslow Cox partial likelihood. The package reports:

- the coefficient and hazard ratio;
- robust sandwich SE conditional on the estimated IPTW;
- model-based SE separately;
- weighted Kaplan–Meier point estimates;
- ESS and covariate-balance diagnostics.

The optimized Python weighted-Cox implementation is checked against a slower transparent Python reference. More importantly, a dedicated workflow independently fits the deterministic validation cohort in R using `survival::coxph(..., weights=..., ties="breslow", robust=TRUE, cluster=patient_id)`.

Current fixed-cohort reconciliation:

- coefficient absolute difference: **5.36 × 10^-8**;
- robust-SE absolute difference: **3.93 × 10^-9**;
- both are required to be ≤ `1e-6`.

This is the strongest evidence that the implemented weighted survival calculation is numerically reproducible across languages.

### 3:45–4:30 — Show source-derived sensitivities

Open **`source_derived_sensitivity.json`** and explain that the same source-derived longitudinal cohort drives:

- 1st/99th percentile weight trimming;
- induced missingness with complete-case and median-imputed analyses;
- bootstrap HR uncertainty;
- weighted PH screening;
- a source-derived negative-control outcome;
- a stricter confirmed-CKD phenotype.

The important engineering point is lineage: these are not sensitivities run on a separate convenience dataset.

### 4:30–5:00 — Show the real-data boundary

The NHANES workflow uses real public CDC participant data for treatment-PS diagnostics only. In the current complete-case gout sample there are 291 participants, including 86 current ULT users and 205 untreated participants.

Ordinary stabilised IPTW reduces maximum absolute SMD from 0.3637 to 0.1799, which is intentionally reported rather than hidden. A separate overlap-weighting sensitivity gives much tighter balance, but targets the common-support population and is explicitly labelled as a different estimand.

Close with: **the synthetic layer demonstrates end-to-end longitudinal RWE programming; the NHANES layer demonstrates that the PS diagnostics also operate on real public data without pretending NHANES supports the longitudinal CKD effect analysis.**

## Ten-minute version

Use the five-minute structure, then add the following detail.

### A. Study specification and time zero

Describe how the analysis specification is translated into code before modelling:

- eligibility and source-history exclusions;
- index-date construction;
- baseline and follow-up windows;
- exposure/comparator definition;
- outcome timing;
- censoring/observation bounds.

A useful sentence is: **“I treated time zero and cohort construction as programming objects that must be testable, rather than as prose that only appears in a protocol.”**

### B. QC architecture

The longitudinal QC registry contains 18 checks covering keys, enrollment coverage, baseline laboratory availability, medication timing, post-index outcomes, follow-up bounds, exposure reproducibility, baseline CKD exclusion, strict-CKD consistency and negative-control availability.

Explain the difference between:

- **primary source-to-analysis QC**, which protects the reviewer-facing longitudinal study;
- **legacy/reference QC**, retained only for regression/provenance and clearly labelled so it cannot be mistaken for the primary path.

### C. Why independent implementations matter

There are several intentionally different validation mechanisms:

1. Python versus SQLite SQL for cohort construction.
2. Fast versus transparent-reference Python for weighted Cox risk sets.
3. Python versus R `survival::coxph` for weighted Cox coefficient and robust SE.
4. Current executable CI evidence versus archived historical claims.

This separation is useful in regulated or review-heavy work because agreement between genuinely different implementations is stronger evidence than rerunning the same function twice.

### D. Production/reviewer delivery

The reviewer bundle contains analysis specification, data dictionary, source inventory, longitudinal QC, cohort-builder reconciliation, Table 1/2/3, weighted survival outputs, sensitivities, validation report and machine-readable inventory.

Explain that the bundle is built in GitHub Actions after tests pass. This connects code validation to the delivered evidence rather than leaving reviewer files as manually assembled outputs.

### E. Real-data PS layer

For NHANES, explain three deliberate constraints:

1. Current serum urate is descriptive but excluded from the default treatment PS because it may be downstream of current ULT in the cross-sectional timing structure.
2. Ordinary IPTW residual imbalance is reported even though max |SMD| remains above 0.10.
3. Overlap weighting is not described as a rescue of IPTW; it changes the target population/estimand.

This is a good place to demonstrate judgement rather than only coding ability.

## Likely interviewer questions

### “Why use synthetic longitudinal data at all?”

Because a public repository cannot contain proprietary claims/EHR data, while the purpose of this layer is to expose the complete longitudinal programming lineage. Synthetic data allow every source table, exclusion, event window and patient-level reconciliation to be auditable. I then added a separate NHANES layer so the PS diagnostics are also demonstrated on real public participant data.

### “Why not use NHANES for the CKD outcome analysis?”

Because the current NHANES use is cross-sectional for medication and laboratory information. It is suitable for a real-data treatment-model and balance demonstration, not for claiming a longitudinal ULT-to-CKD treatment effect.

### “Why are serum urate and prior flares not in the synthetic PS?”

Because they define the UCG phenotype. Treating exposure-defining variables as ordinary confounders in the propensity model creates deterministic separation/positivity problems. They remain visible descriptively in Table 1, while balance QC is applied to the prespecified pre-exposure adjustment set.

### “Does the robust SE include uncertainty from estimating the propensity score?”

No. The reported sandwich variance is conditional on the estimated IPTW. The repository states that limitation explicitly. The independent R reconciliation validates the implemented weighted Cox and robust sandwich calculation, not the full uncertainty of a two-stage estimated-weight procedure.

### “Why implement Cox yourself instead of only calling a library?”

The goal is not to replace validated libraries. The explicit implementation makes weighted risk-set calculations auditable and testable. It is then independently reconciled against R `survival::coxph`, so the custom implementation is used as transparent programming evidence rather than as an unvalidated substitute for standard software.

### “What would you change with access to production RWD?”

Replace the synthetic source generator with governed claims/EHR ingestion; preserve the same protocol-to-code structure, source-history eligibility logic, time-zero checks, independent cohort validation, PS diagnostics, survival implementation validation, sensitivity framework and reviewer manifests. I would also use the organisation's validated statistical stack and data model rather than treating this public portfolio implementation as production software.

## Claims to make

- The repository demonstrates end-to-end source-to-analysis observational programming.
- SQL and Python longitudinal cohort builders are reconciled patient by patient.
- Weighted Cox coefficient and robust SE are independently reconciled against R on a fixed validation cohort.
- The reviewer-facing sensitivity suite uses the source-derived longitudinal cohort.
- The public NHANES layer uses real CDC participant data for PS/overlap/balance diagnostics.
- Current executable CI is the authoritative validation evidence.

## Claims not to make

- Do not call the synthetic findings clinical evidence.
- Do not imply the repository contains Amgen or proprietary RWD.
- Do not call the NHANES workflow a longitudinal ULT→CKD causal analysis.
- Do not present overlap weighting as the same estimand as ordinary IPTW.
- Do not claim full OMOP CDM compliance; the repository contains OMOP-shaped reference reconstruction only.
- Do not claim the robust SE propagates PS-estimation uncertainty.
- Do not quote historical 52-test/60-QC counts as current executable counts.

## One-sentence close

> The main thing I wanted to demonstrate is not a particular synthetic hazard ratio; it is that I can take an observational study specification, translate it into auditable longitudinal programming, validate the cohort and statistical implementation independently, expose limitations, and deliver reviewer-ready evidence.