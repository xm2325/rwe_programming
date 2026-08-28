from rwe_programming.metadata import data_dictionary_frame
from rwe_programming.qc import cohort_attrition, qc_manifest, run_qc_registry
from rwe_programming.source_population import make_synthetic_source_population, select_analysis_cohort
from rwe_programming.study_config import DEFAULT_CONFIG, StudyConfig


def test_default_study_config_is_reproducible():
    assert DEFAULT_CONFIG.study_id == "RWE-GOUT-CKD-001"
    assert DEFAULT_CONFIG.study_version == "0.9.0-restored"
    assert DEFAULT_CONFIG.n_patients == 9184
    assert DEFAULT_CONFIG.seed == 20260817


def test_attrition_ledger_has_real_exclusions_and_reconciles():
    config = StudyConfig(n_patients=1000)
    source = make_synthetic_source_population(n_eligible=1000, seed=config.seed)
    analysis = select_analysis_cohort(source, config)
    attrition = cohort_attrition(source, config)
    assert attrition.iloc[0].included_n > 1000
    assert attrition.iloc[-1].included_n == 1000
    assert attrition.excluded_at_step_n.sum() == len(source) - 1000
    assert (attrition.iloc[1:].excluded_at_step_n > 0).all()
    assert len(analysis) == 1000


def test_runtime_qc_manifest_passes_default_study():
    checks = run_qc_registry()
    manifest = qc_manifest()
    assert len(checks) == 14
    assert all(check.passed for check in checks)
    assert manifest["status"] == "PASS"
    assert manifest["checks_passed"] == manifest["checks_total"] == 14
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
