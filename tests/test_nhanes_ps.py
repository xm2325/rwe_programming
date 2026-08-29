import numpy as np
import pandas as pd

from rwe_programming.nhanes_ps import (
    PS_COVARIATES,
    fit_nhanes_propensity_score,
    nhanes_balance_table,
    nhanes_overlap_table,
    nhanes_ps_diagnostics,
    nhanes_ps_qc_manifest,
    nhanes_weight_diagnostics,
    prepare_nhanes_gout_ps,
)


def _component_frames(n=80):
    seqn = np.arange(1000, 1000 + n)
    demo = pd.DataFrame({
        "SEQN": seqn,
        "RIDAGEYR": np.linspace(40, 80, n),
        "RIAGENDR": np.where(np.arange(n) % 2 == 0, 1, 2),
        "RIDRETH3": 1 + (np.arange(n) % 5),
        "WTMEC2YR": np.linspace(1000, 3000, n),
    })
    # Deliberately mixed case to protect case-insensitive SAS/XPT handling.
    mcq = pd.DataFrame({"SEQN": seqn, "MCQ160n": 1})
    biopro = pd.DataFrame({
        "SEQN": seqn,
        "LBXSCR": 0.8 + 0.01 * (np.arange(n) % 20),
        "LBXSUA": 5.0 + 0.08 * np.arange(n),
    })
    bmx = pd.DataFrame({"SEQN": seqn, "BMXBMI": 22 + 0.15 * np.arange(n)})
    diq = pd.DataFrame({"SEQN": seqn, "DIQ010": np.where(np.arange(n) % 5 == 0, 1, 2)})
    bpq = pd.DataFrame({"SEQN": seqn, "BPQ020": np.where(np.arange(n) % 3 == 0, 1, 2)})
    rx = pd.DataFrame({
        "SEQN": seqn,
        "RXDDRUG": np.where(np.arange(n) % 2 == 0, "ALLOPURINOL", "LISINOPRIL"),
    })
    return demo, mcq, rx, biopro, bmx, diq, bpq


def test_prepare_nhanes_ps_detects_ult_and_keeps_urate_descriptive_only():
    frame = prepare_nhanes_gout_ps(*_component_frames())
    assert len(frame) == 80
    assert frame.ult_use.sum() == 40
    assert "serum_urate" in frame.columns
    assert "serum_urate" not in PS_COVARIATES
    assert set(frame.ult_use.unique()) == {0, 1}
    assert frame.survey_weight.notna().all()


def test_nhanes_ps_fit_and_aggregate_evidence_are_finite_without_effect_claim():
    frame = prepare_nhanes_gout_ps(*_component_frames())
    weighted = fit_nhanes_propensity_score(frame)
    diagnostic = nhanes_ps_diagnostics(weighted)
    balance = nhanes_balance_table(weighted)
    overlap = nhanes_overlap_table(weighted)
    weights = nhanes_weight_diagnostics(weighted)
    qc = nhanes_ps_qc_manifest(weighted)

    assert weighted.propensity_score.between(0.01, 0.99).all()
    assert (weighted.stabilized_weight > 0).all()
    assert np.isfinite(weighted.survey_iptw_weight).all()
    assert diagnostic["n_complete"] == 80
    assert diagnostic["treated_n"] == 40
    assert diagnostic["untreated_n"] == 40
    assert diagnostic["causal_effect_claim"] is False
    assert diagnostic["effective_sample_size"] > 0
    assert "race_ethnicity_1" in set(balance.variable)
    assert np.isfinite(balance.filter(like="smd").to_numpy(float)).all()
    assert overlap[["treated_n", "untreated_n"]].to_numpy().sum() == 80
    assert set(weights.weight) == {"stabilized_weight", "survey_iptw_weight"}
    assert qc["status"] == "PASS"
    assert qc["checks_passed"] == qc["checks_total"] == 7
    assert len(qc["content_sha256"]) == 64
