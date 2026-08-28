from pathlib import Path

from rwe_programming.longitudinal import build_analysis_cohort_python, make_longitudinal_sources
from rwe_programming.longitudinal_validation import reconcile_longitudinal_builders


def test_longitudinal_sources_have_expected_domains():
    sources = make_longitudinal_sources(n=500, seed=20260817)
    assert set(sources) == {"patients", "enrollment", "diagnoses", "labs", "medications", "outcomes"}
    assert len(sources["patients"]) == 500
    assert len(sources["enrollment"]) == 500
    assert len(sources["labs"]) == 1000
    assert len(sources["medications"]) == 500
    assert set(sources["labs"].lab_name.unique()) == {"EGFR", "SERUM_URATE"}


def test_python_builder_respects_time_zero_and_windows():
    sources = make_longitudinal_sources(n=500, seed=20260817)
    cohort = build_analysis_cohort_python(sources)
    assert len(cohort) == 500
    assert cohort.patient_id.is_unique
    assert (cohort.age >= 18).all()
    assert (cohort.egfr >= 45).all()
    assert (cohort.followup_years > 0).all()
    assert (cohort.followup_years <= 5.0 + 1e-9).all()
    assert set(cohort.ucg.unique()).issubset({0, 1})
    assert set(cohort.ckd_event.unique()).issubset({0, 1})


def test_sql_and_python_longitudinal_builders_reconcile():
    result = reconcile_longitudinal_builders(
        n=600,
        seed=20260817,
        sql_path=Path("sql/longitudinal_analysis_cohort.sql"),
    )
    assert result["rows_python"] == result["rows_sql"] == 600
    assert result["patient_level_discrepancies"] == 0
    assert result["max_absolute_numeric_difference"] < 1e-10
