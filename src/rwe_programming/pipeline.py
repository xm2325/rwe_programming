from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import expit
from sklearn.linear_model import LogisticRegression

RNG_SEED = 20260817
N_PATIENTS = 9184
COVARS = ["age", "female", "diabetes", "hypertension", "egfr", "baseline_urate", "prior_flares"]


def make_synthetic_cohort(n: int = N_PATIENTS, seed: int = RNG_SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    age = np.clip(rng.normal(58, 12, n), 18, 90)
    female = rng.binomial(1, 0.43, n)
    diabetes = rng.binomial(1, expit(-2.0 + 0.025 * (age - 50)), n)
    hypertension = rng.binomial(1, expit(-1.3 + 0.035 * (age - 50)), n)
    egfr = np.clip(rng.normal(92 - 0.35 * (age - 50) - 7 * diabetes, 14, n), 45, 140)
    baseline_urate = np.clip(rng.normal(7.4 + 0.35 * diabetes + 0.2 * hypertension, 1.1, n), 3.5, 13.0)
    prior_flares = rng.poisson(np.clip(0.7 + 0.18 * (baseline_urate - 6) + 0.25 * diabetes, 0.1, 4), n)
    logit_ucg = -0.6 + 0.28 * (baseline_urate - 7) + 0.18 * prior_flares + 0.25 * diabetes + 0.15 * hypertension
    ucg = rng.binomial(1, expit(logit_ucg), n)
    linpred = 0.025 * (age - 58) + 0.45 * diabetes + 0.35 * hypertension - 0.018 * (egfr - 90) + 0.08 * ucg
    rate = 0.035 * np.exp(linpred)
    event_time = rng.exponential(1 / rate)
    censor_time = rng.uniform(1.0, 5.0, n)
    followup = np.minimum(event_time, censor_time)
    ckd_event = (event_time <= censor_time).astype(int)

    # Negative-control outcome: generated without an exposure effect.
    neg_rate = 0.025 * np.exp(0.02 * (age - 58) + 0.1 * diabetes)
    neg_time = rng.exponential(1 / neg_rate)
    negative_event = (neg_time <= censor_time).astype(int)
    negative_followup = np.minimum(neg_time, censor_time)

    return pd.DataFrame({
        "patient_id": np.arange(1, n + 1),
        "age": age,
        "female": female,
        "diabetes": diabetes,
        "hypertension": hypertension,
        "egfr": egfr,
        "baseline_urate": baseline_urate,
        "prior_flares": prior_flares,
        "ucg": ucg,
        "followup_years": followup,
        "ckd_event": ckd_event,
        "negative_followup_years": negative_followup,
        "negative_event": negative_event,
    })


def propensity_weights(df: pd.DataFrame, covars: list[str] | None = None) -> pd.DataFrame:
    covars = COVARS if covars is None else covars
    model = LogisticRegression(max_iter=2000)
    model.fit(df[covars], df["ucg"])
    ps = np.clip(model.predict_proba(df[covars])[:, 1], 0.01, 0.99)
    p_t = float(df["ucg"].mean())
    sw = np.where(df["ucg"].eq(1), p_t / ps, (1 - p_t) / (1 - ps))
    out = df.copy()
    out["propensity_score"] = ps
    out["stabilized_weight"] = sw
    return out


def weighted_smd(df: pd.DataFrame, col: str, weight_col: str = "stabilized_weight") -> float:
    t = df["ucg"].to_numpy()
    x = df[col].to_numpy(float)
    w = df[weight_col].to_numpy(float)

    def wmean(mask):
        return np.average(x[mask], weights=w[mask])

    def wvar(mask, mean):
        return np.average((x[mask] - mean) ** 2, weights=w[mask])

    m1, m0 = wmean(t == 1), wmean(t == 0)
    v1, v0 = wvar(t == 1, m1), wvar(t == 0, m0)
    denom = np.sqrt((v1 + v0) / 2)
    return 0.0 if denom == 0 else float((m1 - m0) / denom)


def _breslow_components(
    beta: float,
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    weight_col: str,
) -> tuple[float, float, float]:
    time = df[time_col].to_numpy(float)
    event = df[event_col].to_numpy(int)
    x = df["ucg"].to_numpy(float)
    w = df[weight_col].to_numpy(float)
    expbx = np.exp(beta * x)
    loglik = 0.0
    score = 0.0
    information = 0.0

    for event_time in np.unique(time[event == 1]):
        deaths = (time == event_time) & (event == 1)
        risk = time >= event_time
        death_weight = float(w[deaths].sum())
        death_wx = float((w[deaths] * x[deaths]).sum())
        s0 = float((w[risk] * expbx[risk]).sum())
        s1 = float((w[risk] * x[risk] * expbx[risk]).sum())
        s2 = float((w[risk] * x[risk] ** 2 * expbx[risk]).sum())
        mean_x = s1 / s0
        loglik += beta * death_wx - death_weight * np.log(s0)
        score += death_wx - death_weight * mean_x
        information += death_weight * (s2 / s0 - mean_x**2)

    return loglik, score, information


def fit_weighted_cox(
    df: pd.DataFrame,
    time_col: str = "followup_years",
    event_col: str = "ckd_event",
    weight_col: str = "stabilized_weight",
) -> dict[str, float]:
    """Fit a one-exposure IPTW Cox model using a weighted Breslow partial likelihood."""
    result = minimize_scalar(
        lambda b: -_breslow_components(float(b), df, time_col, event_col, weight_col)[0],
        bounds=(-3.0, 3.0),
        method="bounded",
        options={"xatol": 1e-10},
    )
    if not result.success:
        raise RuntimeError(f"Weighted Cox optimisation failed: {result.message}")
    beta = float(result.x)
    _, _, information = _breslow_components(beta, df, time_col, event_col, weight_col)
    if not np.isfinite(information) or information <= 0:
        raise RuntimeError("Weighted Cox information matrix is not positive")
    se = float(np.sqrt(1.0 / information))
    return {
        "coef": beta,
        "hr": float(np.exp(beta)),
        "se": se,
        "ci_low": float(np.exp(beta - 1.96 * se)),
        "ci_high": float(np.exp(beta + 1.96 * se)),
        "method": "IPTW-weighted Cox; Breslow partial likelihood",
    }


def effective_sample_size(weights: pd.Series) -> float:
    return float(weights.sum() ** 2 / weights.pow(2).sum())


def run_pipeline(n: int = N_PATIENTS, seed: int = RNG_SEED) -> dict:
    df = propensity_weights(make_synthetic_cohort(n, seed))
    smds = {c: abs(weighted_smd(df, c)) for c in COVARS}
    return {
        "n": len(df),
        "treated_fraction": float(df.ucg.mean()),
        "max_abs_weighted_smd": max(smds.values()),
        "effective_sample_size": effective_sample_size(df.stabilized_weight),
        "weighted_cox": fit_weighted_cox(df),
    }
