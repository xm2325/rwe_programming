from __future__ import annotations

import sqlite3
import numpy as np
import pandas as pd

from .pipeline import (
    COVARS,
    PS_COVARS,
    effective_sample_size,
    fit_weighted_cox,
    make_synthetic_cohort,
    propensity_weights,
    weighted_smd,
)


def sql_pandas_reconciliation(df: pd.DataFrame | None = None) -> dict:
    df = make_synthetic_cohort() if df is None else df.copy()
    conn = sqlite3.connect(":memory:")
    df.to_sql("source_cohort", conn, index=False)
    query = """
        SELECT patient_id, age, female, diabetes, hypertension, egfr,
               baseline_urate, prior_flares, ucg, followup_years, ckd_event
        FROM source_cohort
        WHERE age >= 18 AND egfr >= 45
        ORDER BY patient_id
    """
    sql = pd.read_sql_query(query, conn)
    pandas = (
        df.loc[(df.age >= 18) & (df.egfr >= 45), sql.columns]
        .sort_values("patient_id")
        .reset_index(drop=True)
    )
    numeric = [c for c in sql.columns if c != "patient_id"]
    mismatches = int(
        (~np.isclose(
            sql[numeric].to_numpy(float),
            pandas[numeric].to_numpy(float),
            rtol=0,
            atol=1e-12,
        )).any(axis=1).sum()
    )
    return {
        "rows_sql": len(sql),
        "rows_pandas": len(pandas),
        "patient_level_discrepancies": mismatches,
    }


def omop_shape_reconciliation(df: pd.DataFrame | None = None) -> dict:
    df = make_synthetic_cohort() if df is None else df.copy()
    person = df[["patient_id", "female"]].rename(
        columns={"patient_id": "person_id", "female": "gender_concept_id"}
    )
    measurement = df[["patient_id", "baseline_urate", "egfr"]].rename(
        columns={"patient_id": "person_id"}
    )
    observation = df[["patient_id", "diabetes", "hypertension", "prior_flares", "ucg"]].rename(
        columns={"patient_id": "person_id"}
    )
    reconstructed = (
        person.merge(measurement, on="person_id")
        .merge(observation, on="person_id")
        .rename(columns={"person_id": "patient_id", "gender_concept_id": "female"})
    )
    cols = [
        "patient_id", "female", "baseline_urate", "egfr",
        "diabetes", "hypertension", "prior_flares", "ucg",
    ]
    source = df[cols].sort_values("patient_id").reset_index(drop=True)
    rebuilt = reconstructed[cols].sort_values("patient_id").reset_index(drop=True)
    mismatches = int(
        (~np.isclose(source.to_numpy(float), rebuilt.to_numpy(float), rtol=0, atol=1e-12))
        .any(axis=1)
        .sum()
    )
    return {
        "patient_level_discrepancies": mismatches,
        "n_person": len(person),
        "n_measurement": len(measurement),
        "n_observation": len(observation),
    }


def trim_weights(df: pd.DataFrame, lower: float = 0.01, upper: float = 0.99) -> pd.DataFrame:
    out = df.copy()
    lo, hi = out.stabilized_weight.quantile([lower, upper])
    out["trimmed_weight"] = out.stabilized_weight.clip(lo, hi)
    return out


def missingness_sensitivity(df: pd.DataFrame | None = None, seed: int = 11) -> dict:
    df = make_synthetic_cohort() if df is None else df.copy()
    rng = np.random.default_rng(seed)
    missing = rng.random(len(df)) < (0.08 + 0.05 * df.diabetes)
    analysis = df.copy()
    analysis.loc[missing, "baseline_urate"] = np.nan

    required = [*COVARS, "ucg", "followup_years", "ckd_event"]
    complete_case = propensity_weights(analysis.dropna(subset=required).copy())
    imputed = analysis.copy()
    imputed["baseline_urate"] = imputed["baseline_urate"].fillna(imputed["baseline_urate"].median())
    imputed = propensity_weights(imputed)

    complete_hr = fit_weighted_cox(complete_case)["hr"]
    imputed_hr = fit_weighted_cox(imputed)["hr"]
    return {
        "source": "supplied_analysis_cohort" if df is not None else "legacy_flat_default",
        "missing_fraction": float(missing.mean()),
        "complete_case_n": int(len(complete_case)),
        "complete_case_hr": complete_hr,
        "median_imputed_hr": imputed_hr,
        "hr_abs_difference": abs(complete_hr - imputed_hr),
    }


def weight_trimming_sensitivity(df: pd.DataFrame | None = None) -> dict:
    base = propensity_weights(make_synthetic_cohort() if df is None else df.copy())
    trimmed = trim_weights(base)
    untrimmed = fit_weighted_cox(base)
    trimmed_cox = fit_weighted_cox(trimmed, weight_col="trimmed_weight")
    return {
        "untrimmed_hr": untrimmed["hr"],
        "trimmed_hr": trimmed_cox["hr"],
        "untrimmed_ess": effective_sample_size(base.stabilized_weight),
        "trimmed_ess": effective_sample_size(trimmed.trimmed_weight),
        "max_abs_smd_trimmed": max(
            abs(weighted_smd(trimmed, c, "trimmed_weight")) for c in PS_COVARS
        ),
    }


def negative_control_analysis(df: pd.DataFrame | None = None) -> dict:
    base = make_synthetic_cohort() if df is None else df.copy()
    required = {"negative_followup_years", "negative_event"}
    if not required.issubset(base.columns):
        raise ValueError("negative-control analysis requires source-derived negative outcome columns")
    weighted = propensity_weights(base)
    result = fit_weighted_cox(
        weighted,
        time_col="negative_followup_years",
        event_col="negative_event",
    )
    result["outcome_definition"] = "source-derived NEGATIVE_CONTROL event" if df is not None else "legacy flat negative-control event"
    return result


def outcome_sensitivity(df: pd.DataFrame | None = None) -> dict:
    base = make_synthetic_cohort() if df is None else df.copy()
    weighted = propensity_weights(base)
    primary = fit_weighted_cox(weighted)

    if {"ckd_strict_followup_years", "ckd_strict_event"}.issubset(weighted.columns):
        alt_cox = fit_weighted_cox(
            weighted,
            time_col="ckd_strict_followup_years",
            event_col="ckd_strict_event",
        )
        definition = "source-derived confirmed CKD phenotype"
        alt_events = int(weighted.ckd_strict_event.sum())
    else:
        alternative = weighted.copy()
        alternative["ckd_event_alt"] = (
            (alternative.ckd_event == 1) | (alternative.egfr < 60)
        ).astype(int)
        alternative["followup_alt"] = np.where(
            (alternative.egfr < 60) & (alternative.ckd_event == 0),
            0.5,
            alternative.followup_years,
        )
        alt_cox = fit_weighted_cox(
            alternative,
            time_col="followup_alt",
            event_col="ckd_event_alt",
        )
        definition = "legacy flat alternative CKD definition"
        alt_events = int(alternative.ckd_event_alt.sum())

    return {
        "primary_hr": primary["hr"],
        "alternative_hr": alt_cox["hr"],
        "alternative_events": alt_events,
        "alternative_definition": definition,
        "absolute_difference": abs(primary["hr"] - alt_cox["hr"]),
    }
