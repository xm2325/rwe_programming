from rwe_programming.metadata import data_dictionary_frame
from rwe_programming.pipeline import make_synthetic_cohort
from rwe_programming.qc import cohort_attrition, qc_manifest, run_qc_registry
from rwe_programming.study_config import DEFAULT_CONFIG, StudyConfig


def test_default_study_config_is_reproducible():
    assert DEFAULT_CONFIG.study_id == "RWE-GOUT-CKD-001"
    assert DEFAULT_CONFIG.n_patients == 9184
    assert DEFAULT_CONFIG.seed == 20260817


def test_attrition_ledger_reconciles():
    raw = make_synthetic_cohort(n=1000, seed=20260817)
    attrition = cohort_attrition(raw, StudyConfig(n_patients=1000))
    assert attrition.iloc[0].included_n == 1000
    assert attrition.iloc[-1].included_n == 1000
    assert attrition.excluded_at_step_n.sum() == 0


def test_runtime_qc_manifest_passes_default_study():
    checks = run_qc_registry()
    manifest = qc_manifest()
    assert len(checks) == 12
    assert all(check.passed for check in checks)
    assert manifest["status"] == "PASS"
    assert manifest["checks_passed"] == manifest["checks_total"] == 12
    assert len(manifest["content_sha256"]) == 64


def test_analysis_data_dictionary_covers_analysis_columns():
    dictionary = data_dictionary_frame()
    expected = {
        "patient_id", "age", "female", "diabetes", "hypertension", "egfr",
        "baseline_urate", "prior_flares", "ucg", "followup_years", "ckd_event",
        "propensity_score", "stabilized_weight",
    }
    assert expected.issubset(set(dictionary.variable))
    assert dictionary.variable.is_unique
