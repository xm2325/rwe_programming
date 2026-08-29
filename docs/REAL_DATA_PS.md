# Real-data propensity-score layer: NHANES

The repository's full source-to-outcome longitudinal workflow remains synthetic so that it is fully public, reproducible and free of restricted patient data. A separate real-data layer uses public CDC NHANES 2017–2018 data to demonstrate propensity-score treatment modelling on real participants.

## Why NHANES

NHANES is a public CDC/NCHS survey with linked respondent identifiers across demographic, questionnaire, examination, laboratory and prescription-medication components. The 2017–2018 release includes the variables required for a defensible gout treatment-model demonstration:

- `MCQ160N`: doctor ever told the participant they had gout;
- `RXDDRUG`: standardised generic prescription-drug name;
- `RIDAGEYR`, `RIAGENDR`, `RIDRETH3`: age, sex and race/ethnicity;
- `BMXBMI`: body mass index;
- `LBXSCR`: serum creatinine;
- `LBXSUA`: serum uric acid;
- `DIQ010`: doctor-diagnosed diabetes;
- `BPQ020`: doctor-diagnosed hypertension;
- `WTMEC2YR`: MEC examination survey weight.

The downloader in `scripts/download_nhanes_2017_2018.py` retrieves public XPT component files from CDC's `Public/2017/DataFiles` endpoint and rejects any response that does not carry a SAS XPORT `HEADER RECORD` signature. Raw participant files remain outside version control.

## Population and exposure

Population: adults reporting doctor-diagnosed gout.

Exposure: current urate-lowering therapy (ULT) identified from prescription-medication generic names. The vocabulary includes allopurinol, febuxostat, probenecid and pegloticase.

The exposure is used to demonstrate treatment propensity, common support, weighting and balance on real data. It is not treated as a longitudinal treatment-effect study.

## Default propensity-score covariates

The prespecified model uses:

- age;
- sex;
- race/ethnicity;
- BMI;
- diabetes;
- hypertension;
- serum creatinine.

Current serum urate is deliberately excluded from the treatment propensity model. NHANES medication use and laboratory measurements are cross-sectional, so serum urate may already have been affected by ULT and cannot safely be assumed to be a pre-treatment confounder. It remains available for descriptive analyses only.

The primary treatment model is a near-unpenalised logistic fit on these covariates. A separate survey-weighted logistic propensity model uses the MEC examination weights as a sensitivity specification.

## Weighting estimands

The module deliberately reports more than one weighting strategy rather than tuning a propensity model until one balance number looks favourable.

### Stabilised IPTW

Stabilised inverse-probability-of-treatment weights target the broader treatment-comparison pseudo-population represented by the fitted PS model. Residual balance is reported even when it exceeds a conventional 0.10 absolute-SMD target.

### Overlap weighting

An explicit overlap-weighting sensitivity assigns `1-PS` to treated participants and `PS` to untreated participants. This changes the target population to participants with the strongest treatment-choice overlap/common support; it must not be presented as the same estimand as IPTW. A survey×overlap sensitivity is constructed from the independently survey-weighted propensity model.

## Executed public-data result

The dedicated GitHub Actions workflow successfully downloaded and signature-verified the official CDC XPT files, ran the analysis, removed participant-level rows and uploaded aggregate evidence.

For the current validated NHANES 2017–2018 complete-case analysis:

- complete-case adults with doctor-diagnosed gout: **291**;
- current ULT: **86**;
- untreated: **205**;
- ULT fraction: **0.2955**;
- fitted PS range: **0.0245–0.9132**;
- observed treated/untreated common-support interval: **0.0365–0.6433**.

Balance and effective-sample-size results:

| Weighting analysis | Max absolute SMD | ESS |
|---|---:|---:|
| Unweighted | 0.3637 | — |
| Stabilised IPTW | **0.1799** | **225.5** |
| PS overlap weights | **0.000096** | **214.7** |
| Survey-weighted PS × IPTW | **0.1392** | **109.9** |
| Survey-weighted PS × overlap | **0.000039** | **83.2** |

These results are intentionally reported side by side. In particular, the ordinary IPTW result does **not** meet a 0.10 maximum absolute-SMD target in this small real-data sample. The overlap result gives near-exact balance on the fitted model covariates because it targets the common-support population; it is not a post hoc replacement for the IPTW estimand.

## Reviewer evidence outputs

`prepare_nhanes_gout_ps()` harmonises the CDC component tables. `fit_nhanes_propensity_score()` estimates the ordinary and survey-weighted treatment propensity models plus IPTW/overlap weights.

The real-data workflow emits aggregate evidence:

- `nhanes_balance_table.csv` — unweighted, IPTW, overlap, survey-IPTW and survey-overlap SMDs, including race/ethnicity indicator levels;
- `nhanes_ps_overlap_bins.csv` — treated/untreated counts across PS bins;
- `nhanes_weight_diagnostics.csv` — weight quantiles and ESS;
- `nhanes_ps_diagnostics.json` — sample, support and headline balance diagnostics;
- `nhanes_real_data_qc_manifest.json` — machine-readable QC sign-off.

The CI workflow removes `nhanes_gout_ps_analysis.csv` before artifact upload, so published workflow evidence is aggregate-only.

## Interpretation boundary

This module provides real-participant evidence that the programming workflow can ingest public epidemiological data, define a clinically meaningful treatment exposure, estimate propensity scores and audit overlap/balance. It should not be described as showing that ULT causes a change in CKD risk, serum urate or another outcome because NHANES does not establish the longitudinal temporal ordering required for that causal claim. The machine-readable diagnostics therefore retain `causal_effect_claim = false`.

For a real longitudinal treatment-effect extension, a credentialed EHR/claims source such as MIMIC-IV can be considered subject to its data-use agreement and local-compute restrictions. Restricted patient data must never be committed to this repository or uploaded to third-party services.
