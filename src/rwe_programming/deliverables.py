from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .pipeline import COVARS, fit_weighted_cox, make_synthetic_cohort, propensity_weights, weighted_smd
from .survival import survival_at_times, weighted_kaplan_meier


ANALYSIS_SPEC = {
    "study_design": "retrospective synthetic longitudinal cohort",
    "population": "adult synthetic patients with gout; baseline eGFR >= 45 in the final analysis cohort",
    "index_date": "synthetic cohort entry",
    "exposure": {"0": "controlled gout", "1": "uncontrolled gout"},
    "primary_outcome": "incident CKD event during follow-up",
    "estimand": "hazard ratio for uncontrolled versus controlled gout in the IPTW pseudo-population",
    "confounding_adjustment": "stabilised inverse-probability-of-treatment weighting",
    "propensity_covariates": COVARS,
    "survival_model": "IPTW-weighted Cox proportional hazards using an explicit weighted Breslow partial likelihood",
    "uncertainty": "subject-level sandwich standard error for the primary Cox interval; model-based SE retained for comparison; non-parametric bootstrap sensitivity",
    "weighted_survival": "IPTW Kaplan-Meier point estimates by exposure group; no naive Greenwood CI attached",
    "independent_validation": "custom unweighted Breslow implementation reconciled against statsmodels PHReg on the supported unweighted estimand",
    "qc": [
        "patient-id uniqueness",
        "source-to-analysis attrition reconciliation",
        "propensity bounds",
        "weight positivity",
        "post-weighting SMD",
        "effective sample size",
        "independent Cox reconciliation",
    ],
}


def table1_balance(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in COVARS:
        for exposure in (0, 1):
            g = df[df.ucg == exposure]
            rows.append({
                "covariate": col,
                "ucg": exposure,
                "n": len(g),
                "mean": float(g[col].mean()),
                "sd": float(g[col].std(ddof=1)),
            })
    out = pd.DataFrame(rows)
    smd = {c: abs(weighted_smd(df, c)) for c in COVARS}
    out["abs_weighted_smd"] = out["covariate"].map(smd)
    return out


def table2_outcomes(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for exposure in (0, 1):
        g = df[df.ucg == exposure]
        rows.append({
            "ucg": exposure,
            "n": len(g),
            "events": int(g.ckd_event.sum()),
            "event_percent": float(100 * g.ckd_event.mean()),
            "person_years": float(g.followup_years.sum()),
        })
    return pd.DataFrame(rows)


def table3_primary_effect(df: pd.DataFrame) -> pd.DataFrame:
    cox = fit_weighted_cox(df)
    return pd.DataFrame([{
        "contrast": "UCG vs CG",
        "hazard_ratio": cox["hr"],
        "ci95_low": cox["ci_low"],
        "ci95_high": cox["ci_high"],
        "log_hr": cox["coef"],
        "robust_se": cox["robust_se"],
        "model_based_se": cox["model_based_se"],
    }])


def write_deliverables(output_dir: str | Path, n: int = 9184, seed: int = 20260817) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = propensity_weights(make_synthetic_cohort(n=n, seed=seed))
    curve = weighted_kaplan_meier(df)
    paths = {
        "analysis_spec": out / "analysis_spec.json",
        "table1": out / "table1_balance.csv",
        "table2": out / "table2_outcomes.csv",
        "table3": out / "table3_primary_effect.csv",
        "weighted_km_curve": out / "weighted_km_curve.csv",
        "weighted_survival_summary": out / "weighted_survival_summary.csv",
    }
    paths["analysis_spec"].write_text(json.dumps(ANALYSIS_SPEC, indent=2), encoding="utf-8")
    table1_balance(df).to_csv(paths["table1"], index=False)
    table2_outcomes(df).to_csv(paths["table2"], index=False)
    table3_primary_effect(df).to_csv(paths["table3"], index=False)
    curve.to_csv(paths["weighted_km_curve"], index=False)
    survival_at_times(curve).to_csv(paths["weighted_survival_summary"], index=False)
    return {k: str(v) for k, v in paths.items()}
