from rwe_programming.report import build_report
from rwe_programming.sensitivity import bootstrap_hr, proportional_hazards_diagnostic
from rwe_programming.validation import (
    missingness_sensitivity,
    negative_control_analysis,
    omop_shape_reconciliation,
    outcome_sensitivity,
    sql_pandas_reconciliation,
    weight_trimming_sensitivity,
)


def test_sql_pandas_zero_discrepancies():
    assert sql_pandas_reconciliation()["patient_level_discrepancies"] == 0


def test_omop_shape_zero_discrepancies():
    assert omop_shape_reconciliation()["patient_level_discrepancies"] == 0


def test_weight_trimming_finite():
    result = weight_trimming_sensitivity()
    assert result["trimmed_hr"] > 0
    assert result["trimmed_ess"] > 0
    assert result["max_abs_smd_trimmed"] < 0.05


def test_missingness_sensitivity_finite():
    result = missingness_sensitivity()
    assert 0 < result["missing_fraction"] < 0.2
    assert result["complete_case_hr"] > 0
    assert result["median_imputed_hr"] > 0


def test_negative_control_near_null():
    result = negative_control_analysis()
    assert 0.7 < result["hr"] < 1.3


def test_outcome_sensitivity_finite():
    result = outcome_sensitivity()
    assert result["primary_hr"] > 0
    assert result["alternative_hr"] > 0


def test_bootstrap_runs():
    result = bootstrap_hr(n_boot=8)
    assert result["n_boot"] == 8
    assert result["p2_5"] <= result["median_hr"] <= result["p97_5"]


def test_ph_diagnostic_runs():
    result = proportional_hazards_diagnostic()
    assert "rho" in result
    assert "p_value" in result


def test_report_builds(tmp_path):
    assert build_report(tmp_path / "report.html", n_boot=5).exists()
