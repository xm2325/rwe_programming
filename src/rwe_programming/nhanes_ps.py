from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression


ULT_DRUG_TERMS = (
    "allopurinol",
    "febuxostat",
    "probenecid",
    "pegloticase",
)

PS_COVARIATES = [
    "age",
    "female",
    "race_ethnicity",
    "bmi",
    "diabetes",
    "hypertension",
    "serum_creatinine",
]


@dataclass(frozen=True)
class NHANESPSDefinition:
    cycle: str = "2017-2018"
    population: str = "Adults with doctor-diagnosed gout"
    exposure: str = "Current urate-lowering therapy use reported in the prescription-medication questionnaire"
    estimand: str = "Propensity-score overlap and covariate-balance demonstration; no causal treatment-effect claim"
    ps_specification: str = "near-unpenalized logistic MLE on prespecified baseline covariates; survey-weighted logistic sensitivity fitted separately"
    gout_variable: str = "MCQ160N"
    age_variable: str = "RIDAGEYR"
    sex_variable: str = "RIAGENDR"
    race_variable: str = "RIDRETH3"
    bmi_variable: str = "BMXBMI"
    creatinine_variable: str = "LBXSCR"
    urate_variable: str = "LBXSUA"
    survey_weight_variable: str = "WTMEC2YR"


def _normalise_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out.columns = [str(c).upper() for c in out.columns]
    return out


def _binary_yes(series: pd.Series) -> pd.Series:
    return series.eq(1).astype(int)


def _ult_patient_ids(rx: pd.DataFrame) -> set[int]:
    rx = _normalise_columns(rx)
    if "RXDDRUG" not in rx or "SEQN" not in rx:
        raise ValueError("RXQ_RX input must contain SEQN and RXDDRUG")
    names = rx["RXDDRUG"].fillna("").astype(str).str.lower()
    mask = np.zeros(len(rx), dtype=bool)
    for term in ULT_DRUG_TERMS:
        mask |= names.str.contains(term, regex=False).to_numpy()
    return set(rx.loc[mask, "SEQN"].astype(int).tolist())


def prepare_nhanes_gout_ps(
    demo: pd.DataFrame,
    mcq: pd.DataFrame,
    rx: pd.DataFrame,
    biopro: pd.DataFrame,
    bmx: pd.DataFrame,
    diabetes: pd.DataFrame | None = None,
    blood_pressure: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the public NHANES gout treatment-PS analysis frame.

    All incoming XPT column names are normalised case-insensitively because SAS/Pandas
    readers may expose questionnaire variable names with different case. Serum urate is
    intentionally descriptive only: current urate can be downstream of current ULT use
    in this cross-sectional cycle and is therefore excluded from the default PS model.
    """
    demo = _normalise_columns(demo)
    mcq = _normalise_columns(mcq)
    rx = _normalise_columns(rx)
    biopro = _normalise_columns(biopro)
    bmx = _normalise_columns(bmx)
    diabetes = None if diabetes is None else _normalise_columns(diabetes)
    blood_pressure = None if blood_pressure is None else _normalise_columns(blood_pressure)

    required_demo = {"SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3"}
    required_mcq = {"SEQN", "MCQ160N"}
    required_biopro = {"SEQN", "LBXSCR", "LBXSUA"}
    required_bmx = {"SEQN", "BMXBMI"}
    for label, frame, required in (
        ("DEMO", demo, required_demo),
        ("MCQ", mcq, required_mcq),
        ("BIOPRO", biopro, required_biopro),
        ("BMX", bmx, required_bmx),
    ):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{label} input missing required columns: {sorted(missing)}")

    cols_demo = ["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3"]
    if "WTMEC2YR" in demo.columns:
        cols_demo.append("WTMEC2YR")
    out = demo[cols_demo].merge(mcq[["SEQN", "MCQ160N"]], on="SEQN", how="inner")
    out = out[out["MCQ160N"].eq(1)].copy()
    out = out.merge(biopro[["SEQN", "LBXSCR", "LBXSUA"]], on="SEQN", how="left")
    out = out.merge(bmx[["SEQN", "BMXBMI"]], on="SEQN", how="left")

    if diabetes is not None and {"SEQN", "DIQ010"}.issubset(diabetes.columns):
        out = out.merge(diabetes[["SEQN", "DIQ010"]], on="SEQN", how="left")
        out["diabetes"] = _binary_yes(out["DIQ010"])
    else:
        out["diabetes"] = 0

    if blood_pressure is not None and {"SEQN", "BPQ020"}.issubset(blood_pressure.columns):
        out = out.merge(blood_pressure[["SEQN", "BPQ020"]], on="SEQN", how="left")
        out["hypertension"] = _binary_yes(out["BPQ020"])
    else:
        out["hypertension"] = 0

    ult_ids = _ult_patient_ids(rx)
    out["ult_use"] = out["SEQN"].astype(int).isin(ult_ids).astype(int)
    out["age"] = pd.to_numeric(out["RIDAGEYR"], errors="coerce")
    out["female"] = out["RIAGENDR"].eq(2).astype(int)
    out["race_ethnicity"] = pd.to_numeric(out["RIDRETH3"], errors="coerce")
    out["bmi"] = pd.to_numeric(out["BMXBMI"], errors="coerce")
    out["serum_creatinine"] = pd.to_numeric(out["LBXSCR"], errors="coerce")
    out["serum_urate"] = pd.to_numeric(out["LBXSUA"], errors="coerce")
    if "WTMEC2YR" in out.columns:
        out["survey_weight"] = pd.to_numeric(out["WTMEC2YR"], errors="coerce")

    keep = ["SEQN", "ult_use", *PS_COVARIATES, "serum_urate"]
    if "survey_weight" in out.columns:
        keep.append("survey_weight")
    return out[keep].sort_values("SEQN").reset_index(drop=True)


def _model_matrix(frame: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    categorical = [c for c in covariates if c == "race_ethnicity"]
    return pd.get_dummies(frame[covariates], columns=categorical, drop_first=False, dtype=float)


def _standardize_matrix(x: pd.DataFrame) -> np.ndarray:
    values = x.to_numpy(float)
    mean = values.mean(axis=0)
    scale = values.std(axis=0)
    scale = np.where(scale > 0, scale, 1.0)
    return (values - mean) / scale


def _fit_ps(x: pd.DataFrame, y: pd.Series, sample_weight: np.ndarray | None = None) -> np.ndarray:
    # C=1e6 approximates an unpenalised logistic MLE while retaining broad sklearn
    # version compatibility. Standardisation avoids numerical scale artefacts.
    model = LogisticRegression(C=1e6, max_iter=5000, solver="lbfgs")
    model.fit(_standardize_matrix(x), y, sample_weight=sample_weight)
    return np.clip(model.predict_proba(_standardize_matrix(x))[:, 1], 0.01, 0.99)


def _stabilized_weights(y: pd.Series, ps: np.ndarray, treated_fraction: float) -> np.ndarray:
    return np.where(
        y.eq(1),
        treated_fraction / ps,
        (1.0 - treated_fraction) / (1.0 - ps),
    )


def fit_nhanes_propensity_score(
    frame: pd.DataFrame,
    covariates: list[str] | None = None,
) -> pd.DataFrame:
    covariates = PS_COVARIATES if covariates is None else covariates
    needed = {"ult_use", *covariates}
    missing = needed.difference(frame.columns)
    if missing:
        raise ValueError(f"PS frame missing required columns: {sorted(missing)}")

    complete = frame.dropna(subset=covariates + ["ult_use"]).copy()
    if complete["ult_use"].nunique() != 2:
        raise ValueError("Both treated and untreated gout patients are required to estimate a propensity score")

    x = _model_matrix(complete, covariates)
    y = complete["ult_use"]
    ps = _fit_ps(x, y)
    treated_fraction = float(y.mean())
    complete["propensity_score"] = ps
    complete["stabilized_weight"] = _stabilized_weights(y, ps, treated_fraction)

    if "survey_weight" in complete.columns:
        survey = complete["survey_weight"].to_numpy(float)
        survey_fit_weight = survey / np.nanmean(survey)
        survey_ps = _fit_ps(x, y, sample_weight=survey_fit_weight)
        survey_treated_fraction = float(np.average(y.to_numpy(float), weights=survey))
        survey_sw = _stabilized_weights(y, survey_ps, survey_treated_fraction)
        combined = survey * survey_sw
        complete["survey_propensity_score"] = survey_ps
        complete["survey_stabilized_weight"] = survey_sw
        complete["survey_iptw_weight"] = combined / np.nanmean(combined)
    return complete.reset_index(drop=True)


def _smd_arrays(x: np.ndarray, z: np.ndarray, w: np.ndarray) -> float:
    def mean_var(mask: np.ndarray) -> tuple[float, float]:
        ww = w[mask]
        xx = x[mask]
        mean = float(np.average(xx, weights=ww))
        var = float(np.average((xx - mean) ** 2, weights=ww))
        return mean, var

    m1, v1 = mean_var(z == 1)
    m0, v0 = mean_var(z == 0)
    denom = np.sqrt((v1 + v0) / 2.0)
    return 0.0 if denom == 0 else float((m1 - m0) / denom)


def nhanes_balance_table(frame: pd.DataFrame, covariates: list[str] | None = None) -> pd.DataFrame:
    covariates = PS_COVARIATES if covariates is None else covariates
    x = _model_matrix(frame, covariates)
    z = frame["ult_use"].to_numpy(int)
    rows = []
    for column in x.columns:
        values = x[column].to_numpy(float)
        before = _smd_arrays(values, z, np.ones(len(frame), dtype=float))
        after = _smd_arrays(values, z, frame["stabilized_weight"].to_numpy(float))
        row = {
            "variable": column,
            "smd_unweighted": before,
            "abs_smd_unweighted": abs(before),
            "smd_iptw": after,
            "abs_smd_iptw": abs(after),
        }
        if "survey_iptw_weight" in frame.columns:
            survey_after = _smd_arrays(values, z, frame["survey_iptw_weight"].to_numpy(float))
            row["smd_survey_iptw"] = survey_after
            row["abs_smd_survey_iptw"] = abs(survey_after)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("abs_smd_iptw", ascending=False).reset_index(drop=True)


def nhanes_overlap_table(frame: pd.DataFrame, n_bins: int = 10) -> pd.DataFrame:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins = pd.cut(frame["propensity_score"], edges, include_lowest=True, right=True)
    temp = pd.DataFrame({"ps_bin": bins, "ult_use": frame["ult_use"].astype(int)})
    table = (
        temp.groupby(["ps_bin", "ult_use"], observed=False)
        .size()
        .unstack(fill_value=0)
        .rename(columns={0: "untreated_n", 1: "treated_n"})
        .reset_index()
    )
    for col in ["untreated_n", "treated_n"]:
        if col not in table:
            table[col] = 0
    table["ps_bin"] = table["ps_bin"].astype(str)
    return table[["ps_bin", "untreated_n", "treated_n"]]


def _ess(weights: pd.Series | np.ndarray) -> float:
    w = np.asarray(weights, dtype=float)
    return float(w.sum() ** 2 / np.square(w).sum())


def nhanes_weight_diagnostics(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name in ["stabilized_weight", "survey_stabilized_weight", "survey_iptw_weight"]:
        if name not in frame.columns:
            continue
        w = frame[name].astype(float)
        rows.append({
            "weight": name,
            "min": float(w.min()),
            "p50": float(w.quantile(0.50)),
            "p95": float(w.quantile(0.95)),
            "p99": float(w.quantile(0.99)),
            "max": float(w.max()),
            "effective_sample_size": _ess(w),
        })
    return pd.DataFrame(rows)


def nhanes_ps_diagnostics(frame: pd.DataFrame, covariates: list[str] | None = None) -> dict:
    balance = nhanes_balance_table(frame, covariates)
    ps = frame["propensity_score"]
    w = frame["stabilized_weight"]
    treated = frame.loc[frame.ult_use.eq(1), "propensity_score"]
    untreated = frame.loc[frame.ult_use.eq(0), "propensity_score"]
    overlap_low = max(float(treated.min()), float(untreated.min()))
    overlap_high = min(float(treated.max()), float(untreated.max()))
    result = {
        "definition": NHANESPSDefinition().__dict__,
        "n_complete": int(len(frame)),
        "treated_n": int(frame["ult_use"].sum()),
        "untreated_n": int((1 - frame["ult_use"]).sum()),
        "treated_fraction": float(frame["ult_use"].mean()),
        "propensity_min": float(ps.min()),
        "propensity_max": float(ps.max()),
        "common_support_low": overlap_low,
        "common_support_high": overlap_high,
        "weight_p99": float(w.quantile(0.99)),
        "effective_sample_size": _ess(w),
        "max_abs_smd_before": float(balance["abs_smd_unweighted"].max()),
        "max_abs_smd_after": float(balance["abs_smd_iptw"].max()),
        "causal_effect_claim": False,
        "interpretation": "Real NHANES treatment-model diagnostics only; cross-sectional timing does not support a longitudinal treatment-effect claim.",
    }
    if "abs_smd_survey_iptw" in balance.columns:
        result["max_abs_smd_survey_iptw"] = float(balance["abs_smd_survey_iptw"].max())
        result["survey_iptw_effective_sample_size"] = _ess(frame["survey_iptw_weight"])
    return result


def nhanes_ps_qc_manifest(frame: pd.DataFrame) -> dict:
    diagnostics = nhanes_ps_diagnostics(frame)
    checks = [
        ("RQ001", "Complete-case analysis contains both treatment groups", diagnostics["treated_n"] > 0 and diagnostics["untreated_n"] > 0),
        ("RQ002", "Propensity scores are finite and strictly inside (0,1)", bool(np.isfinite(frame.propensity_score).all() and frame.propensity_score.between(0, 1, inclusive="neither").all())),
        ("RQ003", "Stabilised IPTW are finite and positive", bool(np.isfinite(frame.stabilized_weight).all() and (frame.stabilized_weight > 0).all())),
        ("RQ004", "Treated and untreated propensity-score ranges have common support", diagnostics["common_support_high"] > diagnostics["common_support_low"]),
        ("RQ005", "Effective sample size is positive and no larger than complete-case N", 0 < diagnostics["effective_sample_size"] <= diagnostics["n_complete"] + 1e-9),
        ("RQ006", "Serum urate is not used in the default treatment propensity model", "serum_urate" not in PS_COVARIATES),
        ("RQ007", "Workflow makes no causal treatment-effect claim", diagnostics["causal_effect_claim"] is False),
    ]
    if "survey_iptw_weight" in frame.columns:
        checks.append(("RQ008", "Survey-weighted PS and combined survey×IPTW weights are finite and positive", bool(np.isfinite(frame.survey_propensity_score).all() and np.isfinite(frame.survey_iptw_weight).all() and (frame.survey_iptw_weight > 0).all())))
    manifest = {
        "status": "PASS" if all(passed for _, _, passed in checks) else "FAIL",
        "checks_total": len(checks),
        "checks_passed": int(sum(passed for _, _, passed in checks)),
        "checks": [
            {"check_id": check_id, "description": description, "passed": bool(passed)}
            for check_id, description, passed in checks
        ],
        "diagnostics": diagnostics,
    }
    canonical = json.dumps(manifest, sort_keys=True).encode()
    manifest["content_sha256"] = hashlib.sha256(canonical).hexdigest()
    return manifest
