from __future__ import annotations

from pathlib import Path

import numpy as np

from .longitudinal import (
    build_analysis_cohort_python,
    make_longitudinal_source_population,
    run_sql_cohort_builder,
)


CORE_COLUMNS = [
    "patient_id", "age", "female", "diabetes", "hypertension", "egfr",
    "baseline_urate", "prior_flares", "ucg", "followup_years", "ckd_event",
]


def reconcile_longitudinal_builders(
    n: int = 9184,
    seed: int = 20260817,
    sql_path: str | Path = "sql/longitudinal_analysis_cohort.sql",
) -> dict:
    sources = make_longitudinal_source_population(n_eligible=n, seed=seed)
    py = build_analysis_cohort_python(sources)[CORE_COLUMNS]
    sql = run_sql_cohort_builder(sources, sql_path)[CORE_COLUMNS]

    if len(py) != len(sql):
        return {
            "rows_source_patients": int(len(sources["patients"])),
            "rows_python": len(py),
            "rows_sql": len(sql),
            "patient_level_discrepancies": abs(len(py) - len(sql)),
            "max_absolute_numeric_difference": None,
        }

    id_mismatch = py.patient_id.to_numpy() != sql.patient_id.to_numpy()
    numeric = [c for c in CORE_COLUMNS if c != "patient_id"]
    diff = np.abs(py[numeric].to_numpy(float) - sql[numeric].to_numpy(float))
    row_mismatch = id_mismatch | (~np.isclose(diff, 0.0, rtol=0, atol=1e-10)).any(axis=1)
    return {
        "rows_source_patients": int(len(sources["patients"])),
        "rows_python": len(py),
        "rows_sql": len(sql),
        "baseline_ckd_excluded": int(len(sources["patients"]) - len(py)),
        "patient_level_discrepancies": int(row_mismatch.sum()),
        "max_absolute_numeric_difference": float(diff.max()) if diff.size else 0.0,
        "source_table_rows": {name: int(len(frame)) for name, frame in sources.items()},
    }
