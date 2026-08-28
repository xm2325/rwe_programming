from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .pipeline import fit_weighted_cox, make_synthetic_cohort, propensity_weights


def bootstrap_hr(
    n_boot: int = 100,
    seed: int = 20260828,
    df: pd.DataFrame | None = None,
) -> dict:
    base = make_synthetic_cohort() if df is None else df.copy()
    rng = np.random.default_rng(seed)
    hrs: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(base), len(base))
        sample = base.iloc[idx].copy().reset_index(drop=True)
        sample["patient_id"] = np.arange(1, len(sample) + 1)
        sample = propensity_weights(sample)
        hrs.append(fit_weighted_cox(sample)["hr"])

    values = np.asarray(hrs)
    return {
        "n_boot": n_boot,
        "median_hr": float(np.median(values)),
        "p2_5": float(np.quantile(values, 0.025)),
        "p97_5": float(np.quantile(values, 0.975)),
    }


def proportional_hazards_diagnostic(df: pd.DataFrame | None = None) -> dict:
    base = make_synthetic_cohort() if df is None else df.copy()
    weighted = propensity_weights(base)
    beta = fit_weighted_cox(weighted)["coef"]
    time = weighted.followup_years.to_numpy(float)
    event = weighted.ckd_event.to_numpy(int)
    x = weighted.ucg.to_numpy(float)
    w = weighted.stabilized_weight.to_numpy(float)
    expbx = np.exp(beta * x)

    event_times: list[float] = []
    residuals: list[float] = []
    for i in np.flatnonzero(event == 1):
        risk = time >= time[i]
        denom = float((w[risk] * expbx[risk]).sum())
        risk_mean = float((w[risk] * x[risk] * expbx[risk]).sum() / denom)
        event_times.append(float(time[i]))
        residuals.append(float(x[i] - risk_mean))

    rho, p_value = spearmanr(np.log(np.asarray(event_times)), np.asarray(residuals))
    return {
        "method": "Weighted-risk-set Schoenfeld residual screen versus log event time",
        "n_events": len(event_times),
        "rho": float(rho),
        "p_value": float(p_value),
        "flag_p_lt_0_05": bool(p_value < 0.05),
    }
