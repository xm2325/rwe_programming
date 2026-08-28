from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class LongitudinalStudyWindows:
    baseline_days: int = 365
    followup_days: int = 1825
    min_age: int = 18
    min_baseline_egfr: float = 45.0
    uncontrolled_urate_threshold: float = 8.0
    uncontrolled_flare_threshold: int = 2


def make_longitudinal_sources(n: int = 9184, seed: int = 20260817) -> dict[str, pd.DataFrame]:
    """Generate a deterministic longitudinal cohort whose members are analysis-eligible."""
    rng = np.random.default_rng(seed)
    patient_id = np.arange(1, n + 1)
    index_date = pd.Timestamp("2023-01-01") + pd.to_timedelta(rng.integers(0, 180, n), unit="D")
    age = np.clip(np.rint(rng.normal(58, 12, n)), 18, 90).astype(int)
    female = rng.binomial(1, 0.43, n)
    birth_year = pd.DatetimeIndex(index_date).year - age

    patients = pd.DataFrame({"patient_id": patient_id, "birth_year": birth_year, "female": female})
    enrollment = pd.DataFrame({
        "patient_id": patient_id,
        "enrollment_start": index_date - pd.to_timedelta(rng.integers(400, 900, n), unit="D"),
        "enrollment_end": index_date + pd.to_timedelta(rng.integers(730, 1826, n), unit="D"),
        "index_date": index_date,
    })

    diabetes = rng.binomial(1, 1 / (1 + np.exp(-(-2.0 + 0.025 * (age - 50)))), n)
    hypertension = rng.binomial(1, 1 / (1 + np.exp(-(-1.3 + 0.035 * (age - 50)))), n)
    egfr = np.clip(rng.normal(92 - 0.35 * (age - 50) - 7 * diabetes, 14, n), 45, 140)
    urate = np.clip(rng.normal(7.4 + 0.35 * diabetes + 0.2 * hypertension, 1.1, n), 3.5, 13.0)
    flares = rng.poisson(np.clip(0.7 + 0.18 * (urate - 6) + 0.25 * diabetes, 0.1, 4), n)

    lab_rows: list[dict] = []
    diagnosis_rows: list[dict] = []
    medication_rows: list[dict] = []
    outcome_rows: list[dict] = []
    for i, pid in enumerate(patient_id):
        idx = index_date[i]
        lab_rows.append({"patient_id": pid, "lab_date": idx - pd.Timedelta(days=int(rng.integers(1, 60))), "lab_name": "EGFR", "value": float(egfr[i])})
        lab_rows.append({"patient_id": pid, "lab_date": idx - pd.Timedelta(days=int(rng.integers(1, 60))), "lab_name": "SERUM_URATE", "value": float(urate[i])})
        if diabetes[i]:
            diagnosis_rows.append({"patient_id": pid, "diagnosis_date": idx - pd.Timedelta(days=int(rng.integers(10, 360))), "code": "DIABETES"})
        if hypertension[i]:
            diagnosis_rows.append({"patient_id": pid, "diagnosis_date": idx - pd.Timedelta(days=int(rng.integers(10, 360))), "code": "HYPERTENSION"})
        for _ in range(int(flares[i])):
            diagnosis_rows.append({"patient_id": pid, "diagnosis_date": idx - pd.Timedelta(days=int(rng.integers(1, 365))), "code": "GOUT_FLARE"})
        medication_rows.append({"patient_id": pid, "drug_date": idx - pd.Timedelta(days=int(rng.integers(1, 90))), "drug_name": "ULT", "days_supply": int(rng.integers(30, 91))})

        ucg = int((urate[i] >= 8.0) or (flares[i] >= 2))
        linpred = 0.025 * (age[i] - 58) + 0.45 * diabetes[i] + 0.35 * hypertension[i] - 0.018 * (egfr[i] - 90) + 0.08 * ucg
        event_time_days = int(max(1, rng.exponential(365.25 / (0.035 * np.exp(linpred)))))
        followup_limit = min(1825, int((enrollment.iloc[i].enrollment_end - idx).days))
        if event_time_days <= followup_limit:
            outcome_rows.append({"patient_id": pid, "outcome_date": idx + pd.Timedelta(days=event_time_days), "outcome": "INCIDENT_CKD"})

    return {
        "patients": patients,
        "enrollment": enrollment,
        "diagnoses": pd.DataFrame(diagnosis_rows, columns=["patient_id", "diagnosis_date", "code"]),
        "labs": pd.DataFrame(lab_rows),
        "medications": pd.DataFrame(medication_rows),
        "outcomes": pd.DataFrame(outcome_rows, columns=["patient_id", "outcome_date", "outcome"]),
    }


def make_longitudinal_source_population(
    n_eligible: int = 9184,
    seed: int = 20260817,
    baseline_ckd_exclusion_fraction: float = 0.05,
) -> dict[str, pd.DataFrame]:
    """Generate source tables containing eligible patients plus explicit baseline-CKD exclusions."""
    eligible = make_longitudinal_sources(n=n_eligible, seed=seed)
    n_excluded = max(1, int(np.ceil(n_eligible * baseline_ckd_exclusion_fraction)))
    excluded = make_longitudinal_sources(n=n_excluded, seed=seed + 731)
    offset = n_eligible

    for frame in excluded.values():
        frame["patient_id"] = frame["patient_id"].astype(int) + offset

    index_dates = excluded["enrollment"].set_index("patient_id").index_date
    baseline_ckd = pd.DataFrame({
        "patient_id": index_dates.index.astype(int),
        "diagnosis_date": index_dates.to_numpy() - pd.to_timedelta(100, unit="D"),
        "code": "BASELINE_CKD",
    })
    excluded["diagnoses"] = pd.concat([excluded["diagnoses"], baseline_ckd], ignore_index=True)

    return {
        name: pd.concat([eligible[name], excluded[name]], ignore_index=True)
        for name in eligible
    }


def build_analysis_cohort_python(sources: dict[str, pd.DataFrame], windows: LongitudinalStudyWindows = LongitudinalStudyWindows()) -> pd.DataFrame:
    patients = sources["patients"].copy()
    enrollment = sources["enrollment"].copy()
    diagnoses = sources["diagnoses"].copy()
    labs = sources["labs"].copy()
    outcomes = sources["outcomes"].copy()

    base = patients.merge(enrollment, on="patient_id", how="inner")
    base["age"] = pd.DatetimeIndex(base.index_date).year - base.birth_year
    base["baseline_start"] = base.index_date - pd.to_timedelta(windows.baseline_days, unit="D")
    base["followup_end"] = pd.concat([
        base.index_date + pd.to_timedelta(windows.followup_days, unit="D"),
        base.enrollment_end,
    ], axis=1).min(axis=1)

    lab = labs.merge(base[["patient_id", "baseline_start", "index_date"]], on="patient_id")
    lab = lab[(lab.lab_date >= lab.baseline_start) & (lab.lab_date < lab.index_date)]
    egfr = lab[lab.lab_name.eq("EGFR")].sort_values("lab_date").groupby("patient_id").tail(1).set_index("patient_id").value
    urate = lab[lab.lab_name.eq("SERUM_URATE")].sort_values("lab_date").groupby("patient_id").tail(1).set_index("patient_id").value

    dx = diagnoses.merge(base[["patient_id", "baseline_start", "index_date"]], on="patient_id")
    dx = dx[(dx.diagnosis_date >= dx.baseline_start) & (dx.diagnosis_date < dx.index_date)]
    flare = dx[dx.code.eq("GOUT_FLARE")].groupby("patient_id").size()
    diabetes = dx[dx.code.eq("DIABETES")].groupby("patient_id").size().gt(0)
    hypertension = dx[dx.code.eq("HYPERTENSION")].groupby("patient_id").size().gt(0)
    baseline_ckd = dx[dx.code.eq("BASELINE_CKD")].groupby("patient_id").size().gt(0)

    out = base.copy()
    out["egfr"] = out.patient_id.map(egfr)
    out["baseline_urate"] = out.patient_id.map(urate)
    out["prior_flares"] = out.patient_id.map(flare).fillna(0).astype(int)
    out["diabetes"] = out.patient_id.map(diabetes).fillna(False).astype(int)
    out["hypertension"] = out.patient_id.map(hypertension).fillna(False).astype(int)
    out["baseline_ckd"] = out.patient_id.map(baseline_ckd).fillna(False).astype(int)
    out = out[(out.age >= windows.min_age) & (out.egfr >= windows.min_baseline_egfr) & out.baseline_ckd.eq(0)].copy()
    out["ucg"] = ((out.baseline_urate >= windows.uncontrolled_urate_threshold) | (out.prior_flares >= windows.uncontrolled_flare_threshold)).astype(int)

    event = outcomes[outcomes.outcome.eq("INCIDENT_CKD")].merge(out[["patient_id", "index_date", "followup_end"]], on="patient_id")
    event = event[(event.outcome_date > event.index_date) & (event.outcome_date <= event.followup_end)]
    first_event = event.groupby("patient_id").outcome_date.min()
    out["event_date"] = out.patient_id.map(first_event)
    out["ckd_event"] = out.event_date.notna().astype(int)
    out["analysis_end"] = out.event_date.fillna(out.followup_end)
    out["followup_years"] = (out.analysis_end - out.index_date).dt.days / 365.25
    cols = ["patient_id", "age", "female", "diabetes", "hypertension", "egfr", "baseline_urate", "prior_flares", "ucg", "followup_years", "ckd_event"]
    return out[cols].sort_values("patient_id").reset_index(drop=True)


def load_sources_to_sqlite(sources: dict[str, pd.DataFrame], conn: sqlite3.Connection) -> None:
    for name, frame in sources.items():
        temp = frame.copy()
        for col in temp.columns:
            if "date" in col or col.endswith("_start") or col.endswith("_end"):
                if pd.api.types.is_datetime64_any_dtype(temp[col]):
                    temp[col] = temp[col].dt.strftime("%Y-%m-%d")
        temp.to_sql(name, conn, index=False, if_exists="replace")


def run_sql_cohort_builder(sources: dict[str, pd.DataFrame], sql_path: str | Path) -> pd.DataFrame:
    conn = sqlite3.connect(":memory:")
    load_sources_to_sqlite(sources, conn)
    sql = Path(sql_path).read_text(encoding="utf-8")
    return pd.read_sql_query(sql, conn).sort_values("patient_id").reset_index(drop=True)
