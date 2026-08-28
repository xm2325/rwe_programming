from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr
from statsmodels.duration.hazard_regression import PHReg

from .pipeline import fit_weighted_cox, make_synthetic_cohort, propensity_weights


def bootstrap_hr(n_boot: int = 100, seed: int = 20260828) -> dict:
    base = make_synthetic_cohort()
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


def proportional_hazards_diagnostic() -> dict:
    weighted = propensity_weights(make_synthetic_cohort())
    model = PHReg(
        weighted.followup_years,
        weighted[["ucg"]],
        status=weighted.ckd_event,
        ties="breslow",
        freq_weights=weighted.stabilized_weight,
    )
    result = model.fit(disp=False)
    residuals = np.asarray(result.schoenfeld_residuals)[:, 0]
    mask = (weighted.ckd_event.to_numpy() == 1) & np.isfinite(residuals)
    rho, p_value = spearmanr(
        np.log(weighted.loc[mask, "followup_years"].to_numpy()),
        residuals[mask],
    )
    return {
        "method": "Spearman correlation of Schoenfeld residuals with log event time (screen)",
        "rho": float(rho),
        "p_value": float(p_value),
        "flag_p_lt_0_05": bool(p_value < 0.05),
    }
