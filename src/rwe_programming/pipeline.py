from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from statsmodels.duration.hazard_regression import PHReg

RNG_SEED = 20260817
N_PATIENTS = 9184


def make_synthetic_cohort(n: int = N_PATIENTS, seed: int = RNG_SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    age = np.clip(rng.normal(58, 12, n), 18, 90)
    female = rng.binomial(1, 0.43, n)
    diabetes = rng.binomial(1, expit(-2.0 + 0.025 * (age - 50)), n)
    hypertension = rng.binomial(1, expit(-1.3 + 0.035 * (age - 50)), n)
    egfr = np.clip(rng.normal(92 - 0.35 * (age - 50) - 7 * diabetes, 14, n), 45, 140)
    baseline_urate = np.clip(
        rng.normal(7.4 + 0.35 * diabetes + 0.2 * hypertension, 1.1, n), 3.5, 13.0
    )
    prior_flares = rng.poisson(
        np.clip(0.7 + 0.18 * (baseline_urate - 6) + 0.25 * diabetes, 0.1, 4), n
    )

    logit_ucg = (
        -0.6
        + 0.28 * (baseline_urate - 7)
        + 0.18 * prior_flares
        + 0.25 * diabetes
        + 0.15 * hypertension
    )
    ucg = rng.binomial(1, expit(logit_ucg), n)

    linpred = (
        0.025 * (age - 58)
        + 0.45 * diabetes
        + 0.35 * hypertension
        - 0.018 * (egfr - 90)
        + 0.08 * ucg
    )
    rate = 0.035 * np.exp(linpred)
    event_time = rng.exponential(1 / rate)
    censor_time = rng.uniform(1.0, 5.0, n)
    followup = np.minimum(event_time, censor_time)
    ckd_event = (event_time <= censor_time).astype(int)

    return pd.DataFrame(
        {
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
        }
    )


def propensity_weights(df: pd.DataFrame) -> pd.DataFrame:
    covars = [
        "age",
        "female",
        "diabetes",
        "hypertension",
        "egfr",
        "baseline_urate",
        "prior_flares",
    ]
    X = df[covars]
    y = df["ucg"]
    model = LogisticRegression(max_iter=2000)
    model.fit(X, y)
    ps = np.clip(model.predict_proba(X)[:, 1], 0.01, 0.99)
    p_t = y.mean()
    sw = np.where(y.eq(1), p_t / ps, (1 - p_t) / (1 - ps))
    out = df.copy()
    out["propensity_score"] = ps
    out["stabilized_weight"] = sw
    return out


def weighted_smd(df: pd.DataFrame, col: str) -> float:
    t = df["ucg"].to_numpy()
    x = df[col].to_numpy(float)
    w = df["stabilized_weight"].to_numpy(float)

    def wmean(mask):
        return np.average(x[mask], weights=w[mask])

    def wvar(mask, mean):
        return np.average((x[mask] - mean) ** 2, weights=w[mask])

    m1, m0 = wmean(t == 1), wmean(t == 0)
    v1, v0 = wvar(t == 1, m1), wvar(t == 0, m0)
    denom = np.sqrt((v1 + v0) / 2)
    return 0.0 if denom == 0 else (m1 - m0) / denom


def fit_weighted_cox(df: pd.DataFrame) -> dict[str, float]:
    model = PHReg(
        endog=df["followup_years"],
        exog=df[["ucg"]],
        status=df["ckd_event"],
        ties="breslow",
        freq_weights=df["stabilized_weight"],
    )
    res = model.fit(disp=False)
    beta = float(res.params[0])
    se = float(res.bse[0])
    return {"coef": beta, "hr": float(np.exp(beta)), "se": se}


def run_pipeline(n: int = N_PATIENTS, seed: int = RNG_SEED) -> dict:
    df = propensity_weights(make_synthetic_cohort(n=n, seed=seed))
    covars = [
        "age",
        "female",
        "diabetes",
        "hypertension",
        "egfr",
        "baseline_urate",
        "prior_flares",
    ]
    smds = {c: abs(weighted_smd(df, c)) for c in covars}
    cox = fit_weighted_cox(df)
    ess = float(
        df["stabilized_weight"].sum() ** 2 / (df["stabilized_weight"] ** 2).sum()
    )
    return {
        "n": int(len(df)),
        "treated_fraction": float(df["ucg"].mean()),
        "max_abs_weighted_smd": float(max(smds.values())),
        "effective_sample_size": ess,
        "weighted_cox": cox,
    }
