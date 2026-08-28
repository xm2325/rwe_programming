from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from statsmodels.duration.hazard_regression import PHReg

from .pipeline import fit_weighted_cox, make_synthetic_cohort, propensity_weights


def _unweighted_breslow_loglik(beta: float, df: pd.DataFrame) -> float:
    time = df["followup_years"].to_numpy(float)
    event = df["ckd_event"].to_numpy(int)
    x = df["ucg"].to_numpy(float)
    expbx = np.exp(beta * x)
    ll = 0.0
    for event_time in np.unique(time[event == 1]):
        deaths = (time == event_time) & (event == 1)
        risk = time >= event_time
        d = int(deaths.sum())
        ll += beta * float(x[deaths].sum())
        ll -= d * np.log(float(expbx[risk].sum()))
    return ll


def fit_independent_unweighted_cox(df: pd.DataFrame) -> dict[str, float]:
    result = minimize_scalar(
        lambda b: -_unweighted_breslow_loglik(float(b), df),
        bounds=(-3.0, 3.0),
        method="bounded",
        options={"xatol": 1e-11},
    )
    if not result.success:
        raise RuntimeError(f"Independent Cox optimisation failed: {result.message}")
    beta = float(result.x)
    return {"coef": beta, "hr": float(np.exp(beta)), "partial_loglik": float(-result.fun)}


def fit_independent_weighted_cox(df: pd.DataFrame) -> dict[str, float]:
    """Compatibility wrapper for the repository's explicit weighted estimator."""
    return fit_weighted_cox(df)


def reconcile_cox(n: int = 9184, seed: int = 20260817) -> dict[str, float | str]:
    raw = make_synthetic_cohort(n=n, seed=seed)
    statsmodels = PHReg(
        endog=raw["followup_years"],
        exog=raw[["ucg"]],
        status=raw["ckd_event"],
        ties="breslow",
    ).fit(disp=False)
    independent = fit_independent_unweighted_cox(raw)
    weighted = fit_weighted_cox(propensity_weights(raw))
    statsmodels_beta = float(statsmodels.params[0])
    return {
        "parity_estimand": "unweighted Cox; Breslow ties",
        "reason": "statsmodels PHReg does not expose observation weights in its documented API",
        "statsmodels_coef": statsmodels_beta,
        "independent_coef": independent["coef"],
        "absolute_coef_difference": abs(statsmodels_beta - independent["coef"]),
        "statsmodels_hr": float(np.exp(statsmodels_beta)),
        "independent_hr": independent["hr"],
        "primary_iptw_hr": weighted["hr"],
    }
