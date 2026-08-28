from __future__ import annotations

import pandas as pd

from .pipeline import COVARS, effective_sample_size, fit_weighted_cox, propensity_weights, weighted_smd


def prepare_weighted_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Return an analysis cohort with propensity scores and stabilised IPTW."""
    required = set(COVARS) | {"patient_id", "ucg", "followup_years", "ckd_event"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Analysis cohort is missing required columns: {missing}")
    if {"propensity_score", "stabilized_weight"}.issubset(df.columns):
        return df.copy()
    return propensity_weights(df.copy())


def summarise_analysis(df: pd.DataFrame) -> dict:
    """Run the primary weighting/balance/survival-effect summary on any valid cohort."""
    weighted = prepare_weighted_analysis(df)
    smds = {c: abs(weighted_smd(weighted, c)) for c in COVARS}
    return {
        "n": int(len(weighted)),
        "treated_fraction": float(weighted.ucg.mean()),
        "max_abs_weighted_smd": float(max(smds.values())),
        "effective_sample_size": effective_sample_size(weighted.stabilized_weight),
        "weighted_cox": fit_weighted_cox(weighted),
    }
