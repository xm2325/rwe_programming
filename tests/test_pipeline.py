from rwe_programming.pipeline import make_synthetic_cohort, propensity_weights, run_pipeline


def test_patient_count_and_unique_ids():
    df = make_synthetic_cohort()
    assert len(df) == 9184
    assert df.patient_id.is_unique


def test_propensity_and_weights_are_finite():
    df = propensity_weights(make_synthetic_cohort())
    assert df.propensity_score.between(0.01, 0.99).all()
    assert (df.stabilized_weight > 0).all()
    assert df.stabilized_weight.notna().all()


def test_balance_improves_to_small_residual_smd():
    result = run_pipeline()
    assert result["max_abs_weighted_smd"] < 0.03


def test_cox_returns_positive_hr():
    result = run_pipeline()
    assert result["weighted_cox"]["hr"] > 0
