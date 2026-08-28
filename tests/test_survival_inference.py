import numpy as np

from rwe_programming.pipeline import (
    _breslow_components,
    _breslow_components_reference,
    _subject_score_residuals,
    fit_weighted_cox,
    make_synthetic_cohort,
    propensity_weights,
)
from rwe_programming.survival import survival_at_times, weighted_kaplan_meier


def test_weighted_cox_returns_sandwich_and_model_based_uncertainty():
    df = propensity_weights(make_synthetic_cohort(n=2000, seed=20260817))
    result = fit_weighted_cox(df)
    assert np.isfinite(result["robust_se"])
    assert np.isfinite(result["model_based_se"])
    assert result["robust_se"] > 0
    assert result["model_based_se"] > 0
    assert result["ci_low"] < result["hr"] < result["ci_high"]
    assert abs(result["score_at_solution"]) < 1e-3
    assert "sandwich" in result["method"].lower()


def test_fast_breslow_components_match_reference_and_subject_scores_reconcile():
    df = propensity_weights(make_synthetic_cohort(n=600, seed=20260817))
    for beta in (-0.4, 0.0, 0.35):
        fast = _breslow_components(beta, df, "followup_years", "ckd_event", "stabilized_weight")
        reference = _breslow_components_reference(beta, df, "followup_years", "ckd_event", "stabilized_weight")
        assert np.allclose(fast, reference, rtol=0, atol=1e-10)
        subject_score = _subject_score_residuals(beta, df, "followup_years", "ckd_event", "stabilized_weight").sum()
        assert np.isclose(subject_score, fast[1], rtol=0, atol=1e-10)


def test_weighted_km_is_monotone_and_starts_at_one():
    df = propensity_weights(make_synthetic_cohort(n=1500, seed=20260817))
    curve = weighted_kaplan_meier(df)
    for _, group in curve.groupby("ucg"):
        group = group.sort_values("time")
        assert group.iloc[0].time == 0.0
        assert group.iloc[0].survival == 1.0
        assert ((group.survival.diff().dropna()) <= 1e-12).all()
        assert group.survival.between(0, 1).all()
        assert (group.weighted_at_risk >= 0).all()
        assert (group.weighted_events >= 0).all()


def test_weighted_survival_summary_has_both_groups_and_requested_times():
    df = propensity_weights(make_synthetic_cohort(n=1500, seed=20260817))
    summary = survival_at_times(weighted_kaplan_meier(df), times=(1.0, 3.0, 5.0))
    assert set(summary.ucg) == {0, 1}
    assert set(summary.time_years) == {1.0, 3.0, 5.0}
    assert len(summary) == 6
    assert summary.survival.between(0, 1).all()
