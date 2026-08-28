import json

import numpy as np
import pandas as pd

from rwe_programming.analysis import prepare_weighted_analysis, summarise_analysis
from rwe_programming.deliverables import write_deliverables
from rwe_programming.longitudinal import build_analysis_cohort_python, make_longitudinal_source_population


def test_source_population_builds_exact_requested_analysis_cohort():
    sources = make_longitudinal_source_population(n_eligible=400, seed=20260817)
    cohort = build_analysis_cohort_python(sources)
    assert len(sources["patients"]) > 400
    assert len(cohort) == 400
    assert set(cohort.patient_id).isdisjoint(
        set(sources["diagnoses"].loc[sources["diagnoses"].code.eq("BASELINE_CKD"), "patient_id"])
    )


def test_source_derived_cohort_runs_primary_iptw_analysis():
    sources = make_longitudinal_source_population(n_eligible=500, seed=20260817)
    cohort = build_analysis_cohort_python(sources)
    weighted = prepare_weighted_analysis(cohort)
    result = summarise_analysis(cohort)

    assert len(weighted) == result["n"] == 500
    assert weighted.propensity_score.between(0.01, 0.99).all()
    assert np.isfinite(weighted.stabilized_weight).all()
    assert (weighted.stabilized_weight > 0).all()
    assert np.isfinite(result["max_abs_weighted_smd"])
    assert result["effective_sample_size"] > 0
    assert np.isfinite(result["weighted_cox"]["hr"])
    assert result["weighted_cox"]["ci_low"] < result["weighted_cox"]["hr"] < result["weighted_cox"]["ci_high"]


def test_reviewer_deliverables_use_supplied_source_derived_cohort(tmp_path):
    sources = make_longitudinal_source_population(n_eligible=350, seed=20260817)
    cohort = build_analysis_cohort_python(sources)
    paths = write_deliverables(tmp_path, analysis_df=cohort)

    table1 = pd.read_csv(paths["table1"])
    table2 = pd.read_csv(paths["table2"])
    table3 = pd.read_csv(paths["table3"])
    spec = json.loads((tmp_path / "analysis_spec.json").read_text())

    assert table1.groupby("ucg").n.first().sum() == 350
    assert table2.n.sum() == 350
    assert len(table3) == 1
    assert table3.hazard_ratio.iloc[0] > 0
    assert "longitudinal source tables" in spec["source_to_analysis"]
