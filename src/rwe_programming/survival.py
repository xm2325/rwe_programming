from __future__ import annotations

import numpy as np
import pandas as pd


def weighted_kaplan_meier(
    df: pd.DataFrame,
    time_col: str = "followup_years",
    event_col: str = "ckd_event",
    exposure_col: str = "ucg",
    weight_col: str = "stabilized_weight",
) -> pd.DataFrame:
    """Compute IPTW Kaplan-Meier curves separately by exposure group.

    This function returns point estimates only. It deliberately does not attach
    naive Greenwood confidence intervals because variance estimation for weighted
    KM curves requires additional assumptions.
    """
    rows: list[dict[str, float | int]] = []
    for exposure in sorted(df[exposure_col].dropna().unique()):
        g = df[df[exposure_col] == exposure].copy()
        times = g[time_col].to_numpy(float)
        events = g[event_col].to_numpy(int)
        weights = g[weight_col].to_numpy(float)
        survival = 1.0
        rows.append({
            "ucg": int(exposure),
            "time": 0.0,
            "survival": 1.0,
            "weighted_at_risk": float(weights.sum()),
            "weighted_events": 0.0,
        })
        for t in np.unique(times[events == 1]):
            risk = times >= t
            deaths = (times == t) & (events == 1)
            risk_weight = float(weights[risk].sum())
            death_weight = float(weights[deaths].sum())
            if risk_weight <= 0:
                continue
            survival *= max(0.0, 1.0 - death_weight / risk_weight)
            rows.append({
                "ucg": int(exposure),
                "time": float(t),
                "survival": float(survival),
                "weighted_at_risk": risk_weight,
                "weighted_events": death_weight,
            })
    return pd.DataFrame(rows)


def survival_at_times(curve: pd.DataFrame, times: tuple[float, ...] = (1.0, 3.0, 5.0)) -> pd.DataFrame:
    rows = []
    for exposure, g in curve.groupby("ucg"):
        g = g.sort_values("time")
        for target in times:
            eligible = g[g.time <= target]
            survival = 1.0 if eligible.empty else float(eligible.iloc[-1].survival)
            rows.append({"ucg": int(exposure), "time_years": float(target), "survival": survival})
    return pd.DataFrame(rows)
