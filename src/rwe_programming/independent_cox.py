from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from .pipeline import fit_weighted_cox, make_synthetic_cohort, propensity_weights


def _breslow_loglik(beta: float, df: pd.DataFrame) -> float:
    """Frequency-weighted Breslow partial log-likelihood for one binary exposure."""
    t = df["followup_years"].to_numpy(float)
    event = df["ckd_event"].to_numpy(int)
    x = df["ucg"].to_numpy(float)
    w = df["stabilized_weight"].to_numpy(float)
    expbx = np.exp(beta * x)
    ll = 0.0
    for event_time in np.unique(t[event == 1]):
        deaths = (t == event_time) & (event == 1)
        d_weight = float(w[deaths].sum())
        if d_weight == 0:
            continue
        risk = t >= event_time
        ll += beta * float((w[deaths] * x[deaths]).sum())
        ll -= d_weight * np.log(float((w[risk] * expbx[risk]).sum()))
    return ll


def fit_independent_weighted_cox(df: pd.DataFrame) -> dict[str, float]:
    result = minimize_scalar(lambda b: -_breslow_loglik(float(b), df), bounds=(-3.0, 3.0), method="bounded", options={"xatol": 1e-11})
    if not result.success:
        raise RuntimeError(f"Independent Cox optimisation failed: {result.message}")
    beta = float(result.x)
    return {"coef": beta, "hr": float(np.exp(beta)), "partial_loglik": float(-result.fun)}


def reconcile_cox(n: int = 9184, seed: int = 20260817) -> dict[str, float]:
    df = propensity_weights(make_synthetic_cohort(n=n, seed=seed))
    reference = fit_weighted_cox(df)
    independent = fit_independent_weighted_cox(df)
    return {
        "statsmodels_coef": reference["coef"],
        "independent_coef": independent["coef"],
        "absolute_coef_difference": abs(reference["coef"] - independent["coef"]),
        "statsmodels_hr": reference["hr"],
        "independent_hr": independent["hr"],
    }
