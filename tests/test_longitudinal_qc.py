from rwe_programming.longitudinal import make_longitudinal_source_population
from rwe_programming.longitudinal_qc import longitudinal_qc_manifest, run_longitudinal_qc
from rwe_programming.longitudinal_validation import reconcile_longitudinal_builders


def test_longitudinal_qc_manifest_passes():
    sources = make_longitudinal_source_population(n_eligible=300, seed=20260817)
    checks = run_longitudinal_qc(sources=sources, n=300, seed=20260817)
    manifest = longitudinal_qc_manifest(sources=sources, n=300, seed=20260817)
    assert len(checks) == 18
    assert all(check.passed for check in checks)
    assert manifest["status"] == "PASS"
    assert manifest["checks_passed"] == manifest["checks_total"] == 18
    assert len(manifest["content_sha256"]) == 64


def test_longitudinal_sql_python_reconciliation_is_exact_after_ckd_exclusion():
    result = reconcile_longitudinal_builders(n=250, seed=20260817)
    assert result["rows_source_patients"] > 250
    assert result["baseline_ckd_excluded"] > 0
    assert result["rows_python"] == result["rows_sql"] == 250
    assert result["patient_level_discrepancies"] == 0
    assert result["max_absolute_numeric_difference"] < 1e-10
