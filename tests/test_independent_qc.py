from rwe_programming.deliverables import ANALYSIS_SPEC, table1_balance, table2_outcomes, table3_primary_effect
from rwe_programming.independent_cox import reconcile_cox
from rwe_programming.pipeline import make_synthetic_cohort, propensity_weights


def test_independent_cox_reconciles_tightly():
    result = reconcile_cox(n=2500)
    assert result["absolute_coef_difference"] < 1e-5


def test_analysis_spec_has_core_estimand_fields():
    for key in ["population", "exposure", "primary_outcome", "estimand", "confounding_adjustment", "survival_model"]:
        assert key in ANALYSIS_SPEC


def test_tlf_tables_are_deterministic_and_complete():
    df = propensity_weights(make_synthetic_cohort(n=1000))
    t1 = table1_balance(df)
    t2 = table2_outcomes(df)
    t3 = table3_primary_effect(df)
    assert len(t1) == 14
    assert t1.abs_weighted_smd.max() < 0.08
    assert t2.n.sum() == 1000
    assert t2.events.sum() == int(df.ckd_event.sum())
    assert len(t3) == 1
    assert t3.hazard_ratio.iloc[0] > 0
