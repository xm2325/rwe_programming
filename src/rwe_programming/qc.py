from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .pipeline import COVARS, propensity_weights, weighted_smd
from .source_population import make_synthetic_source_population, select_analysis_cohort
from .study_config import DEFAULT_CONFIG, StudyConfig


@dataclass(frozen=True)
class QCCheck:
    check_id: str
    description: str
    passed: bool
    observed: str
    expected: str
    severity: str = "ERROR"


def cohort_attrition(df: pd.DataFrame, config: StudyConfig = DEFAULT_CONFIG) -> pd.DataFrame:
    steps = []
    current = df.copy()

    def record(order: int, rule: str, frame: pd.DataFrame, previous_n: int | None) -> None:
        n = len(frame)
        excluded = 0 if previous_n is None else previous_n - n
        steps.append({
            "step": order,
            "rule": rule,
            "included_n": int(n),
            "excluded_at_step_n": int(excluded),
            "retained_percent": 100.0 * n / len(df) if len(df) else 0.0,
        })

    record(1, "Synthetic source population", current, None)
    before = len(current)
    current = current[current["age"] >= config.min_age].copy()
    record(2, f"Age >= {config.min_age}", current, before)
    before = len(current)
    current = current[current["egfr"] >= config.min_baseline_egfr].copy()
    record(3, f"Baseline eGFR >= {config.min_baseline_egfr:g}", current, before)
    before = len(current)
    current = current[current["followup_years"] > 0].copy()
    record(4, "Positive follow-up time", current, before)
    before = len(current)
    current = current[current["ucg"].isin([0, 1])].copy()
    record(5, "Valid exposure classification", current, before)
    return pd.DataFrame(steps)


def run_qc_registry(config: StudyConfig = DEFAULT_CONFIG) -> list[QCCheck]:
    source = make_synthetic_source_population(n_eligible=config.n_patients, seed=config.seed)
    raw = select_analysis_cohort(source, config)
    df = propensity_weights(raw)
    smds = {c: abs(weighted_smd(df, c)) for c in COVARS}
    max_smd = max(smds.values())
    ess = float(df["stabilized_weight"].sum() ** 2 / (df["stabilized_weight"] ** 2).sum())
    attrition = cohort_attrition(source, config)
    excluded_total = int(attrition.excluded_at_step_n.sum())

    return [
        QCCheck("QC001", "Expected final analysis patient count", bool(len(raw) == config.n_patients), str(len(raw)), str(config.n_patients)),
        QCCheck("QC002", "Source patient IDs are unique", bool(source.patient_id.is_unique), str(source.patient_id.nunique()), str(len(source))),
        QCCheck("QC003", "No missing source patient IDs", bool(source.patient_id.notna().all()), str(int(source.patient_id.isna().sum())), "0"),
        QCCheck("QC004", "Final exposure is binary", bool(set(raw.ucg.unique()).issubset({0, 1})), str(sorted(raw.ucg.unique().tolist())), "[0, 1]"),
        QCCheck("QC005", "Final follow-up is positive", bool((raw.followup_years > 0).all()), f"min={raw.followup_years.min():.6f}", ">0"),
        QCCheck("QC006", "Final follow-up within configured horizon", bool((raw.followup_years <= config.max_followup_years).all()), f"max={raw.followup_years.max():.6f}", f"<={config.max_followup_years}"),
        QCCheck("QC007", "Propensity scores bounded", bool(df.propensity_score.between(config.propensity_clip_low, config.propensity_clip_high).all()), f"[{df.propensity_score.min():.6f}, {df.propensity_score.max():.6f}]", f"[{config.propensity_clip_low}, {config.propensity_clip_high}]"),
        QCCheck("QC008", "Stabilised weights positive and finite", bool(np.isfinite(df.stabilized_weight).all() and (df.stabilized_weight > 0).all()), f"[{df.stabilized_weight.min():.6f}, {df.stabilized_weight.max():.6f}]", "finite and >0"),
        QCCheck("QC009", "Post-weighting balance threshold", bool(max_smd < config.balance_threshold_abs_smd), f"max abs SMD={max_smd:.6f}", f"<{config.balance_threshold_abs_smd}"),
        QCCheck("QC010", "Effective sample size remains substantial", bool(ess >= 0.70 * len(df)), f"ESS={ess:.1f}", f">={0.70 * len(df):.1f}"),
        QCCheck("QC011", "Attrition ledger reconciles to final cohort", bool(int(attrition.iloc[-1].included_n) == len(raw)), str(int(attrition.iloc[-1].included_n)), str(len(raw))),
        QCCheck("QC012", "Primary outcome is binary", bool(set(raw.ckd_event.unique()).issubset({0, 1})), str(sorted(raw.ckd_event.unique().tolist())), "[0, 1]"),
        QCCheck("QC013", "Source population contains staged exclusions", bool(excluded_total > 0), str(excluded_total), ">0"),
        QCCheck("QC014", "Every eligibility stage excludes records", bool((attrition.iloc[1:].excluded_at_step_n > 0).all()), str(attrition.iloc[1:].excluded_at_step_n.tolist()), "all >0"),
    ]


def qc_manifest(config: StudyConfig = DEFAULT_CONFIG) -> dict:
    checks = run_qc_registry(config)
    payload = {
        "study_id": config.study_id,
        "study_version": config.study_version,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(c.passed for c in checks) else "FAIL",
        "checks_total": len(checks),
        "checks_passed": int(sum(c.passed for c in checks)),
        "checks": [asdict(c) for c in checks],
    }
    canonical = json.dumps(
        {k: v for k, v in payload.items() if k != "generated_at_utc"},
        sort_keys=True,
    ).encode()
    payload["content_sha256"] = hashlib.sha256(canonical).hexdigest()
    return payload


def write_qc_manifest(path: str | Path, config: StudyConfig = DEFAULT_CONFIG) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(qc_manifest(config), indent=2), encoding="utf-8")
    return path
