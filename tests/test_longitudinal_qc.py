from rwe_programming.longitudinal import make_longitudinal_sources
from rwe_programming.longitudinal_qc import longitudinal_qc_manifest, run_longitudinal_qc
from rwe_programming.longitudinal_validation import reconcile_longitudinal_builders


def test_longitudinal_qc_manifest_passes():
    sources = make_longitudinal_sources(n=300, seed=20260817)
    checks = run_longitudinal_qc(sources=sources, n=300, seed=20260817)
    manifest = longitudinal_qc_manifest(sources=sources, n=300, seed=20260817)
    assert len(checks) == 15
    assert all(check.passed for check in checks)
    assert manifest["status"] == "PASS"
    assert manifest["checks_passed"] == manifest["checks_total"] == 15
    assert len(manifest["content_sha256"]) == 64


def test_longitudinal_sql_python_reconciliation_is_exact():
    result = reconcile_longitudinal_builders(n=250, seed=20260817)
    assert result["rows_python"] == result["rows_sql"] == 250
    assert result["patient_level_discrepancies"] == 0
    assert result["max_absolute_numeric_difference"] < 1e-10
