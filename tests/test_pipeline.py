import numpy as np

from rwe_programming.pipeline import (
    COVARS,
    PS_COVARS,
    make_synthetic_cohort,
    propensity_weights,
    run_pipeline,
)


def test_patient_count_and_unique_ids():
    df = make_synthetic_cohort()
    assert len(df) == 9184
    assert df.patient_id.is_unique


def test_propensity_and_weights_are_finite():
    df = propensity_weights(make_synthetic_cohort())
    assert df.propensity_score.between(0.01, 0.99).all()
    assert (df.stabilized_weight > 0).all()
    assert df.stabilized_weight.notna().all()


def test_propensity_model_excludes_phenotype_defining_variables():
    assert "baseline_urate" in COVARS
    assert "prior_flares" in COVARS
    assert "baseline_urate" not in PS_COVARS
    assert "prior_flares" not in PS_COVARS

    base = make_synthetic_cohort(n=800)
    changed = base.copy()
    changed["baseline_urate"] = changed["baseline_urate"] + 50.0
    changed["prior_flares"] = changed["prior_flares"] + 100

    ps_base = propensity_weights(base).propensity_score.to_numpy()
    ps_changed = propensity_weights(changed).propensity_score.to_numpy()
    assert np.allclose(ps_base, ps_changed, rtol=0, atol=1e-12)


def test_balance_improves_to_small_residual_smd():
    result = run_pipeline()
    assert result["max_abs_weighted_smd"] < 0.03


def test_cox_returns_positive_hr():
    result = run_pipeline()
    assert result["weighted_cox"]["hr"] > 0
