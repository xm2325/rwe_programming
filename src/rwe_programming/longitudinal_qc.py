from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .longitudinal import (
    LongitudinalStudyWindows,
    build_analysis_cohort_python,
    make_longitudinal_source_population,
)


@dataclass(frozen=True)
class LongitudinalQCCheck:
    check_id: str
    description: str
    passed: bool
    observed: str
    expected: str
    severity: str = "ERROR"


def run_longitudinal_qc(
    sources: dict[str, pd.DataFrame] | None = None,
    n: int = 9184,
    seed: int = 20260817,
    windows: LongitudinalStudyWindows = LongitudinalStudyWindows(),
) -> list[LongitudinalQCCheck]:
    sources = make_longitudinal_source_population(n_eligible=n, seed=seed) if sources is None else sources
    patients = sources["patients"]
    enrollment = sources["enrollment"]
    diagnoses = sources["diagnoses"]
    labs = sources["labs"]
    medications = sources["medications"]
    outcomes = sources["outcomes"]
    patient_ids = set(patients.patient_id.tolist())

    e = enrollment.copy()
    baseline_start = e.index_date - pd.to_timedelta(windows.baseline_days, unit="D")
    analysis = build_analysis_cohort_python(sources, windows)

    lab_window = labs.merge(e[["patient_id", "index_date"]], on="patient_id", how="left")
    lab_window["baseline_start"] = lab_window.index_date - pd.to_timedelta(windows.baseline_days, unit="D")
    baseline_labs = lab_window[(lab_window.lab_date >= lab_window.baseline_start) & (lab_window.lab_date < lab_window.index_date)]
    egfr_counts = baseline_labs[baseline_labs.lab_name.eq("EGFR")].groupby("patient_id").size()
    urate_counts = baseline_labs[baseline_labs.lab_name.eq("SERUM_URATE")].groupby("patient_id").size()

    med_window = medications.merge(e[["patient_id", "index_date"]], on="patient_id", how="left")
    med_window["baseline_start"] = med_window.index_date - pd.to_timedelta(windows.baseline_days, unit="D")
    meds_in_baseline = (med_window.drug_date >= med_window.baseline_start) & (med_window.drug_date < med_window.index_date)

    out_window = outcomes.merge(e[["patient_id", "index_date", "enrollment_end"]], on="patient_id", how="left")
    outcomes_after_time_zero = (out_window.outcome_date > out_window.index_date).all() if len(out_window) else True
    outcomes_within_enrollment = (out_window.outcome_date <= out_window.enrollment_end).all() if len(out_window) else True

    source_fk_ok = {
        "diagnoses": set(diagnoses.patient_id.tolist()).issubset(patient_ids),
        "labs": set(labs.patient_id.tolist()).issubset(patient_ids),
        "medications": set(medications.patient_id.tolist()).issubset(patient_ids),
        "outcomes": set(outcomes.patient_id.tolist()).issubset(patient_ids),
    }

    expected_ucg = ((analysis.baseline_urate >= windows.uncontrolled_urate_threshold) | (analysis.prior_flares >= windows.uncontrolled_flare_threshold)).astype(int)
    baseline_ckd_ids = set(diagnoses.loc[diagnoses.code.eq("BASELINE_CKD"), "patient_id"].astype(int).tolist())
    analysis_ids = set(analysis.patient_id.astype(int).tolist())
    excluded_ckd_count = len(baseline_ckd_ids)

    return [
        LongitudinalQCCheck("LQ001", "Patient identifiers are unique", bool(patients.patient_id.is_unique), str(patients.patient_id.nunique()), str(len(patients))),
        LongitudinalQCCheck("LQ002", "Exactly one enrollment row per patient", bool(enrollment.patient_id.is_unique and len(enrollment) == len(patients)), str(len(enrollment)), str(len(patients))),
        LongitudinalQCCheck("LQ003", "Continuous baseline enrollment covers the full lookback", bool((e.enrollment_start <= baseline_start).all()), f"violations={int((e.enrollment_start > baseline_start).sum())}", "0"),
        LongitudinalQCCheck("LQ004", "Enrollment extends beyond time zero", bool((e.enrollment_end > e.index_date).all()), f"violations={int((e.enrollment_end <= e.index_date).sum())}", "0"),
        LongitudinalQCCheck("LQ005", "Diagnosis foreign keys resolve to patients", bool(source_fk_ok["diagnoses"]), str(source_fk_ok["diagnoses"]), "True"),
        LongitudinalQCCheck("LQ006", "Laboratory foreign keys resolve to patients", bool(source_fk_ok["labs"]), str(source_fk_ok["labs"]), "True"),
        LongitudinalQCCheck("LQ007", "Medication foreign keys resolve to patients", bool(source_fk_ok["medications"]), str(source_fk_ok["medications"]), "True"),
        LongitudinalQCCheck("LQ008", "Outcome foreign keys resolve to patients", bool(source_fk_ok["outcomes"]), str(source_fk_ok["outcomes"]), "True"),
        LongitudinalQCCheck("LQ009", "Every patient has a baseline eGFR measurement", bool(len(egfr_counts) == len(patients) and (egfr_counts >= 1).all()), str(len(egfr_counts)), str(len(patients))),
        LongitudinalQCCheck("LQ010", "Every patient has a baseline serum urate measurement", bool(len(urate_counts) == len(patients) and (urate_counts >= 1).all()), str(len(urate_counts)), str(len(patients))),
        LongitudinalQCCheck("LQ011", "Medication records occur strictly before time zero within baseline", bool(meds_in_baseline.all()), f"violations={int((~meds_in_baseline).sum())}", "0"),
        LongitudinalQCCheck("LQ012", "Incident outcomes occur strictly after time zero", bool(outcomes_after_time_zero), f"violations={int((out_window.outcome_date <= out_window.index_date).sum()) if len(out_window) else 0}", "0"),
        LongitudinalQCCheck("LQ013", "Incident outcomes do not extend beyond observed enrollment", bool(outcomes_within_enrollment), f"violations={int((out_window.outcome_date > out_window.enrollment_end).sum()) if len(out_window) else 0}", "0"),
        LongitudinalQCCheck("LQ014", "Analysis follow-up is strictly positive and within horizon", bool((analysis.followup_years > 0).all() and (analysis.followup_years <= windows.followup_days / 365.25 + 1e-12).all()), f"range=[{analysis.followup_years.min():.6f}, {analysis.followup_years.max():.6f}]", f">0 and <={windows.followup_days / 365.25:.6f}"),
        LongitudinalQCCheck("LQ015", "Exposure classification is reproducible from baseline urate/flares", bool(np.array_equal(expected_ucg.to_numpy(), analysis.ucg.to_numpy())), f"mismatches={int((expected_ucg.to_numpy() != analysis.ucg.to_numpy()).sum())}", "0"),
        LongitudinalQCCheck("LQ016", "Baseline CKD patients are present in source history and excluded from analysis", bool(excluded_ckd_count > 0 and baseline_ckd_ids.isdisjoint(analysis_ids) and len(analysis) == n), f"baseline_ckd_source={excluded_ckd_count}; analysis_n={len(analysis)}", f"baseline_ckd_source>0; analysis_n={n}; overlap=0"),
    ]


def longitudinal_qc_manifest(
    sources: dict[str, pd.DataFrame] | None = None,
    n: int = 9184,
    seed: int = 20260817,
    windows: LongitudinalStudyWindows = LongitudinalStudyWindows(),
) -> dict:
    checks = run_longitudinal_qc(sources=sources, n=n, seed=seed, windows=windows)
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(check.passed for check in checks) else "FAIL",
        "checks_total": len(checks),
        "checks_passed": int(sum(check.passed for check in checks)),
        "windows": {
            "baseline_days": windows.baseline_days,
            "followup_days": windows.followup_days,
            "min_age": windows.min_age,
            "min_baseline_egfr": windows.min_baseline_egfr,
            "uncontrolled_urate_threshold": windows.uncontrolled_urate_threshold,
            "uncontrolled_flare_threshold": windows.uncontrolled_flare_threshold,
        },
        "checks": [asdict(check) for check in checks],
    }
    canonical = json.dumps({k: v for k, v in payload.items() if k != "generated_at_utc"}, sort_keys=True).encode()
    payload["content_sha256"] = hashlib.sha256(canonical).hexdigest()
    return payload


def write_longitudinal_qc_manifest(
    path: str | Path,
    sources: dict[str, pd.DataFrame] | None = None,
    n: int = 9184,
    seed: int = 20260817,
    windows: LongitudinalStudyWindows = LongitudinalStudyWindows(),
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(longitudinal_qc_manifest(sources, n, seed, windows), indent=2), encoding="utf-8")
    return path
