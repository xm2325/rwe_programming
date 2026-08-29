from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from sklearn.linear_model import LogisticRegression
from statsmodels.duration.hazard_regression import PHReg

RNG_SEED = 20260817
N_PATIENTS = 9184
COVARS = ["age", "female", "diabetes", "hypertension", "egfr", "baseline_urate", "prior_flares"]


def make_synthetic_cohort(n: int = N_PATIENTS, seed: int = RNG_SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    patient_id = np.arange(1, n + 1)
    age = np.clip(np.rint(rng.normal(58, 12, n)), 18, 90).astype(int)
    female = rng.binomial(1, 0.43, n)
    diabetes = rng.binomial(1, 1 / (1 + np.exp(-(-2.0 + 0.025 * (age - 50)))), n)
    hypertension = rng.binomial(1, 1 / (1 + np.exp(-(-1.3 + 0.035 * (age - 50)))), n)
    egfr = np.clip(rng.normal(92 - 0.35 * (age - 50) - 7 * diabetes, 14, n), 45, 140)
    baseline_urate = np.clip(rng.normal(7.4 + 0.35 * diabetes + 0.2 * hypertension, 1.1, n), 3.5, 13.0)
    prior_flares = rng.poisson(np.clip(0.7 + 0.18 * (baseline_urate - 6) + 0.25 * diabetes, 0.1, 4), n)
    ucg = ((baseline_urate >= 8.0) | (prior_flares >= 2)).astype(int)

    linpred = 0.025 * (age - 58) + 0.45 * diabetes + 0.35 * hypertension - 0.018 * (egfr - 90) + 0.08 * ucg
    event_time = rng.exponential(1 / (0.035 * np.exp(linpred)))
    censor = rng.uniform(2.0, 5.0, n)
    ckd_event = (event_time <= censor).astype(int)
    followup_years = np.minimum(event_time, censor)

    # Kept for legacy/reference validation only. The reviewer-facing negative-control
    # endpoint is generated from the longitudinal outcomes source domain.
    neg_event_time = rng.exponential(1 / 0.025, n)
    negative_event = (neg_event_time <= censor).astype(int)
    negative_followup_years = np.minimum(neg_event_time, censor)

    return pd.DataFrame({
        "patient_id": patient_id,
        "age": age,
        "female": female,
        "diabetes": diabetes,
        "hypertension": hypertension,
        "egfr": egfr,
        "baseline_urate": baseline_urate,
        "prior_flares": prior_flares,
        "ucg": ucg,
        "followup_years": followup_years,
        "ckd_event": ckd_event,
        "negative_followup_years": negative_followup_years,
        "negative_event": negative_event,
    })


def propensity_weights(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    x = out[COVARS]
    model = LogisticRegression(max_iter=3000)
    model.fit(x, out["ucg"])
    ps = np.clip(model.predict_proba(x)[:, 1], 0.01, 0.99)
    p = float(out["ucg"].mean())
    out["propensity_score"] = ps
    out["stabilized_weight"] = np.where(out.ucg.eq(1), p / ps, (1 - p) / (1 - ps))
    return out


def _weighted_mean_var(x: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    mean = float(np.average(x, weights=w))
    var = float(np.average((x - mean) ** 2, weights=w))
    return mean, var


def weighted_smd(df: pd.DataFrame, column: str, weight_col: str = "stabilized_weight") -> float:
    t = df.ucg.eq(1).to_numpy()
    x = df[column].to_numpy(float)
    w = df[weight_col].to_numpy(float)
    m1, v1 = _weighted_mean_var(x[t], w[t])
    m0, v0 = _weighted_mean_var(x[~t], w[~t])
    denom = np.sqrt((v1 + v0) / 2)
    return 0.0 if denom == 0 else float((m1 - m0) / denom)


def _breslow_components_reference(
    beta: float,
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    weight_col: str,
) -> tuple[float, float, float]:
    """Transparent O(events × rows) weighted Breslow reference implementation."""
    time = df[time_col].to_numpy(float)
    event = df[event_col].to_numpy(int)
    x = df["ucg"].to_numpy(float)
    w = df[weight_col].to_numpy(float)
    expbx = np.exp(beta * x)
    loglik = score = information = 0.0
    for t in np.unique(time[event == 1]):
        deaths = (time == t) & (event == 1)
        risk = time >= t
        d_weight = float(w[deaths].sum())
        d_wx = float((w[deaths] * x[deaths]).sum())
        rw = w[risk] * expbx[risk]
        s0 = float(rw.sum())
        s1 = float((rw * x[risk]).sum())
        s2 = float((rw * x[risk] ** 2).sum())
        xbar = s1 / s0
        loglik += beta * d_wx - d_weight * np.log(s0)
        score += d_wx - d_weight * xbar
        information += d_weight * (s2 / s0 - xbar**2)
    return loglik, score, information


def _sorted_risk_arrays(
    beta: float,
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    weight_col: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    time = df[time_col].to_numpy(float)
    event = df[event_col].to_numpy(int)
    x = df["ucg"].to_numpy(float)
    w = df[weight_col].to_numpy(float)
    order = np.argsort(-time, kind="mergesort")
    ts = time[order]
    es = event[order]
    xs = x[order]
    ws = w[order]
    expbx = np.exp(beta * xs)
    group_ends = np.r_[np.flatnonzero(ts[:-1] != ts[1:]), len(ts) - 1]
    return order, ts, es, xs, ws, expbx, group_ends


def _breslow_components(
    beta: float,
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    weight_col: str,
) -> tuple[float, float, float]:
    """O(rows log rows) weighted Breslow components using cumulative risk sums."""
    _, _, es, xs, ws, expbx, group_ends = _sorted_risk_arrays(beta, df, time_col, event_col, weight_col)
    risk0 = np.cumsum(ws * expbx)
    risk1 = np.cumsum(ws * xs * expbx)
    risk2 = np.cumsum(ws * xs**2 * expbx)

    loglik = 0.0
    score = 0.0
    information = 0.0
    start = 0
    for end in group_ends:
        group = slice(start, int(end) + 1)
        death = es[group] == 1
        if death.any():
            wg = ws[group][death]
            xg = xs[group][death]
            death_weight = float(wg.sum())
            death_wx = float((wg * xg).sum())
            s0 = float(risk0[end])
            s1 = float(risk1[end])
            s2 = float(risk2[end])
            mean_x = s1 / s0
            loglik += beta * death_wx - death_weight * np.log(s0)
            score += death_wx - death_weight * mean_x
            information += death_weight * (s2 / s0 - mean_x**2)
        start = int(end) + 1
    return loglik, score, information


def _subject_score_residuals(
    beta: float,
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    weight_col: str,
) -> np.ndarray:
    """Case-weighted Breslow score residuals matching R survival's definition.

    For subject i the unweighted score residual is
      event_i * (x_i - xbar(T_i))
      - exp(beta*x_i) * sum_{t <= T_i} dLambda(t) * (x_i - xbar(t)).
    R's ``residuals.coxph(..., type='dfbeta', weighted=TRUE)`` then multiplies
    this residual by the case weight before applying inverse information.  The
    implementation below follows that ordering explicitly.
    """
    order, _, es, xs, ws, expbx, group_ends = _sorted_risk_arrays(
        beta, df, time_col, event_col, weight_col
    )
    risk0 = np.cumsum(ws * expbx)
    risk1 = np.cumsum(ws * xs * expbx)
    hazard_increment = np.zeros(len(group_ends), dtype=float)
    xbar = np.zeros(len(group_ends), dtype=float)

    start = 0
    group_index = np.empty(len(es), dtype=int)
    for g, end in enumerate(group_ends):
        end = int(end)
        group = slice(start, end + 1)
        group_index[start:end + 1] = g
        death = es[group] == 1
        if death.any():
            death_weight = float(ws[group][death].sum())
            hazard_increment[g] = death_weight / float(risk0[end])
            xbar[g] = float(risk1[end]) / float(risk0[end])
        start = end + 1

    # With descending follow-up times, event times <= a subject's observed time
    # are the current and later groups, hence reverse cumulative sums.
    cum_hazard = np.cumsum(hazard_increment[::-1])[::-1]
    cum_xbar_hazard = np.cumsum((hazard_increment * xbar)[::-1])[::-1]
    g = group_index
    event_contribution = es * (xs - xbar[g])
    compensator = expbx * (xs * cum_hazard[g] - cum_xbar_hazard[g])
    u_sorted = ws * (event_contribution - compensator)

    u = np.empty_like(u_sorted)
    u[order] = u_sorted
    return u


def fit_weighted_cox(
    df: pd.DataFrame,
    time_col: str = "followup_years",
    event_col: str = "ckd_event",
    weight_col: str = "stabilized_weight",
) -> dict[str, float | str]:
    """Fit one-exposure IPTW Cox with model-based and sandwich uncertainty."""
    result = minimize_scalar(
        lambda b: -_breslow_components(float(b), df, time_col, event_col, weight_col)[0],
        bounds=(-3.0, 3.0),
        method="bounded",
        options={"xatol": 1e-10},
    )
    if not result.success:
        raise RuntimeError(f"Weighted Cox optimisation failed: {result.message}")
    beta = float(result.x)
    _, score, information = _breslow_components(beta, df, time_col, event_col, weight_col)
    if not np.isfinite(information) or information <= 0:
        raise RuntimeError("Weighted Cox information matrix is not positive")

    model_se = float(np.sqrt(1.0 / information))
    subject_scores = _subject_score_residuals(beta, df, time_col, event_col, weight_col)
    meat = float(np.sum(subject_scores**2))
    robust_var = meat / (information**2)
    robust_se = float(np.sqrt(max(robust_var, 0.0)))

    return {
        "coef": beta,
        "hr": float(np.exp(beta)),
        "se": robust_se,
        "robust_se": robust_se,
        "model_based_se": model_se,
        "ci_low": float(np.exp(beta - 1.96 * robust_se)),
        "ci_high": float(np.exp(beta + 1.96 * robust_se)),
        "score_at_solution": float(score),
        "method": "IPTW-weighted Cox; Breslow partial likelihood; case-weighted score-residual sandwich variance",
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
