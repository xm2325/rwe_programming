from __future__ import annotations

import pandas as pd


ANALYSIS_DATA_DICTIONARY = [
    {"variable": "patient_id", "type": "integer", "role": "identifier", "definition": "Synthetic patient identifier", "allowed": "unique positive integer"},
    {"variable": "age", "type": "float", "role": "baseline covariate", "definition": "Age at cohort entry, years", "allowed": ">=18"},
    {"variable": "female", "type": "integer", "role": "baseline covariate", "definition": "Synthetic sex indicator", "allowed": "0/1"},
    {"variable": "diabetes", "type": "integer", "role": "baseline covariate", "definition": "Baseline diabetes indicator", "allowed": "0/1"},
    {"variable": "hypertension", "type": "integer", "role": "baseline covariate", "definition": "Baseline hypertension indicator", "allowed": "0/1"},
    {"variable": "egfr", "type": "float", "role": "baseline covariate", "definition": "Baseline estimated glomerular filtration rate", "allowed": ">=45 in source generator"},
    {"variable": "baseline_urate", "type": "float", "role": "baseline covariate", "definition": "Baseline serum urate", "allowed": "synthetic continuous"},
    {"variable": "prior_flares", "type": "integer", "role": "baseline covariate", "definition": "Prior gout flare count", "allowed": ">=0"},
    {"variable": "ucg", "type": "integer", "role": "exposure", "definition": "Uncontrolled gout exposure indicator", "allowed": "0=controlled, 1=uncontrolled"},
    {"variable": "followup_years", "type": "float", "role": "time", "definition": "Observed time from index to CKD event or censoring", "allowed": ">0 and <=5"},
    {"variable": "ckd_event", "type": "integer", "role": "primary outcome", "definition": "Incident CKD event indicator", "allowed": "0/1"},
    {"variable": "propensity_score", "type": "float", "role": "derived", "definition": "Estimated probability of uncontrolled gout conditional on baseline covariates", "allowed": "[0.01, 0.99]"},
    {"variable": "stabilized_weight", "type": "float", "role": "derived", "definition": "Stabilised inverse-probability-of-treatment weight", "allowed": ">0 finite"},
]


def data_dictionary_frame() -> pd.DataFrame:
    return pd.DataFrame(ANALYSIS_DATA_DICTIONARY)
