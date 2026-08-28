# Real-data propensity-score layer: NHANES

The repository's full source-to-outcome longitudinal workflow remains synthetic so that it is fully public, reproducible and free of restricted patient data. A separate real-data layer uses public NHANES data to demonstrate propensity-score treatment modelling on real participants.

## Why NHANES

NHANES is a public CDC/NCHS survey with linked respondent identifiers across demographic, questionnaire, examination, laboratory and prescription-medication components. The 2017–2018 release includes the variables required for a defensible gout treatment-model demonstration:

- `MCQ160n`: doctor ever told the participant they had gout;
- `RXDDRUG`: standardised generic prescription-drug name;
- `RIDAGEYR`, `RIAGENDR`, `RIDRETH3`: age, sex and race/ethnicity;
- `BMXBMI`: body mass index;
- `LBXSCR`: serum creatinine;
- `LBXSUA`: serum uric acid;
- `DIQ010`: doctor-diagnosed diabetes;
- `BPQ020`: doctor-diagnosed hypertension;
- `WTMEC2YR`: MEC examination survey weight where available.

The downloader in `scripts/download_nhanes_2017_2018.py` retrieves the public XPT component files directly from CDC. Raw NHANES files should remain outside version control.

## Population and exposure

Population: adults reporting doctor-diagnosed gout.

Exposure: current urate-lowering therapy (ULT) use identified from prescription-medication generic names. The initial vocabulary includes allopurinol, febuxostat, probenecid and pegloticase.

The exposure is intended to demonstrate treatment propensity, overlap, weighting and balance on real data. It is not treated as a longitudinal treatment-effect study.

## Default propensity-score covariates

The default model uses:

- age;
- sex;
- race/ethnicity;
- BMI;
- diabetes;
- hypertension;
- serum creatinine.

Current serum urate is deliberately excluded from the default treatment propensity score. NHANES medication use and laboratory measurements are cross-sectional, so serum urate may already have been affected by ULT and cannot safely be assumed to be a pre-treatment confounder. It remains available for descriptive analyses only.

## Outputs

`prepare_nhanes_gout_ps()` harmonises CDC-style component tables into one analysis frame. `fit_nhanes_propensity_score()` estimates treatment propensity and stabilised IPTW. If the NHANES MEC survey weight is available, a normalised survey-weight × IPTW field is also produced for sensitivity work.

`nhanes_ps_diagnostics()` reports:

- treated and untreated counts;
- propensity-score range;
- 99th percentile of stabilised weights;
- effective sample size;
- maximum absolute SMD before and after PS weighting;
- an explicit `causal_effect_claim = false` flag.

## Interpretation boundary

This module provides real-patient evidence that the programming workflow can ingest public epidemiological data, define a clinically meaningful treatment exposure, estimate propensity scores and audit overlap/balance. It should not be described as showing that ULT causes a change in CKD risk, serum urate or another outcome because NHANES does not establish the longitudinal temporal ordering required for that causal claim.

For a real longitudinal treatment-effect extension, a credentialed EHR/claims source such as MIMIC-IV can be considered, subject to its data-use agreement and local-compute restrictions. Restricted patient data must never be committed to this repository or uploaded to third-party services.
