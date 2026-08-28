from __future__ import annotations

from dataclasses import dataclass

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
    gout_variable: str = "MCQ160n"
    age_variable: str = "RIDAGEYR"
    sex_variable: str = "RIAGENDR"
    race_variable: str = "RIDRETH3"
    bmi_variable: str = "BMXBMI"
    creatinine_variable: str = "LBXSCR"
    urate_variable: str = "LBXSUA"
    survey_weight_variable: str = "WTMEC2YR"


def _binary_yes(series: pd.Series) -> pd.Series:
    return series.eq(1).astype(int)


def _ult_patient_ids(rx: pd.DataFrame) -> set[int]:
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
    """Build a real-data PS analysis frame from public NHANES component tables.

    The function intentionally keeps serum urate as a descriptive variable only. Because
    NHANES medication use and laboratory measurements are cross-sectional, current serum
    urate may be downstream of urate-lowering treatment and is therefore not included in
    the default treatment propensity-score covariate set.
    """
    required_demo = {"SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3"}
    required_mcq = {"SEQN", "MCQ160n"}
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
    out = demo[cols_demo].merge(mcq[["SEQN", "MCQ160n"]], on="SEQN", how="inner")
    out = out[out["MCQ160n"].eq(1)].copy()
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

    x = pd.get_dummies(complete[covariates], columns=["race_ethnicity"], drop_first=False, dtype=float)
    model = LogisticRegression(max_iter=3000)
    model.fit(x, complete["ult_use"])
    ps = np.clip(model.predict_proba(x)[:, 1], 0.01, 0.99)
    treated_fraction = float(complete["ult_use"].mean())
    sw = np.where(
        complete["ult_use"].eq(1),
        treated_fraction / ps,
        (1.0 - treated_fraction) / (1.0 - ps),
    )
    complete["propensity_score"] = ps
    complete["stabilized_weight"] = sw
    if "survey_weight" in complete.columns:
        combined = complete["survey_weight"].to_numpy(float) * sw
        complete["survey_iptw_weight"] = combined / np.nanmean(combined)
    return complete


def _smd(frame: pd.DataFrame, column: str, weight: str | None = None) -> float:
    x = frame[column].to_numpy(float)
    z = frame["ult_use"].to_numpy(int)
    w = np.ones(len(frame), dtype=float) if weight is None else frame[weight].to_numpy(float)

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


def nhanes_ps_diagnostics(frame: pd.DataFrame, covariates: list[str] | None = None) -> dict:
    covariates = PS_COVARIATES if covariates is None else covariates
    numeric_balance = [c for c in covariates if c != "race_ethnicity"]
    before = {c: abs(_smd(frame, c)) for c in numeric_balance}
    after = {c: abs(_smd(frame, c, "stabilized_weight")) for c in numeric_balance}
    ps = frame["propensity_score"]
    w = frame["stabilized_weight"]
    ess = float(w.sum() ** 2 / np.square(w).sum())
    return {
        "n_complete": int(len(frame)),
        "treated_n": int(frame["ult_use"].sum()),
        "untreated_n": int((1 - frame["ult_use"]).sum()),
        "propensity_min": float(ps.min()),
        "propensity_max": float(ps.max()),
        "weight_p99": float(w.quantile(0.99)),
        "effective_sample_size": ess,
        "max_abs_smd_before": float(max(before.values())) if before else 0.0,
        "max_abs_smd_after": float(max(after.values())) if after else 0.0,
        "causal_effect_claim": False,
        "interpretation": "Real NHANES treatment-model diagnostics only; cross-sectional timing does not support a longitudinal treatment-effect claim.",
    }
